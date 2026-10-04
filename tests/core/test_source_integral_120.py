"""Independent mathematical and conservative checks for the new face operator.

Synthetic fixtures are mathematical tests, never replacement site materials.
"""
import copy,json,tempfile,unittest
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import quad
from hmsl.material import Material,read_material,RW
from hmsl.face_integral import LiquidFaceIntegral
from hmsl.solver import ColumnModel,State
from hmsl.project import Project
from hmsl.climate import read_climate
from hmsl.common import InputError,sha256
ROOT=Path(__file__).resolve().parents[1]

class SourceIntegralTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.path=Path(self.tmp.name)
        self.cfg=json.loads((ROOT/'projects/demo/project.json').read_text())
        self.obj=json.loads((ROOT/'projects/demo/materials/SYNTHETIC.json').read_text())
        self.cfg['numerics']['liquid_face_scheme']='source_integral'
    def tearDown(self):self.tmp.cleanup()
    def model(self):return ColumnModel(self.cfg,Material(copy.deepcopy(self.obj)))
    def climate(self,hours=24,wet=False):
        t=pd.date_range('2020-01-01',periods=hours,freq='h')
        d=pd.DataFrame({'date':t.strftime('%Y-%m-%d'),'hour_utc':t.hour,'tas_C':26.498 if wet else 25.,
          'hurs_pct':94.505 if wet else 60.,'pr_mm_h':8.27748 if wet else 0.,'rsds_W_m2':0.,
          'rlds_W_m2':444.718 if wet else 380.,'sfcWind_m_s':1.842 if wet else 2.})
        p=self.path/'forcing.csv';d.to_csv(p,index=False);return read_climate(p)
    def linear_diffusivity_material(self,D=1e-7):
        o=copy.deepcopy(self.obj);o['functions']={}
        lo=.001;hi=float(o['theta_effective']);a=(3.-8.)/(hi-lo)
        const=np.log10(RW*D/(-np.log(10.)*a))
        o['functions']['pC(Theta_l)_de']=[[lo,hi],[8.,3.]]
        o['functions']['lgKl(Theta_l)']=[[lo,hi],[const-8.,const-3.]]
        o['functions']['lgKv(Theta_l)']=[[lo,hi],[-50.,-50.]]
        return Material(o)
    def test_analytic_integral_matches_piecewise_quadrature(self):
        m=self.model().m;op=LiquidFaceIntegral(m);rng=np.random.default_rng(629)
        for low,high in np.sort(rng.uniform(m.lo,m.hi,size=(12,2)),axis=1):
            def f(t):
                i=min(np.searchsorted(op.x,t,side='right')-1,len(op.a)-1)
                return op.pref[i]*np.exp(op.rate[i]*(t-op.x[i]))
            points=op.x[(op.x>low)&(op.x<high)]
            ref=quad(f,low,high,points=points,epsabs=1e-30,epsrel=1e-10,limit=400)[0]
            self.assertAlmostEqual(float(op.integral(low,high))/ref,1.,delta=1e-9)
    def test_constant_diffusivity_primitive_exact(self):
        D=1e-7;m=self.linear_diffusivity_material(D);op=LiquidFaceIntegral(m)
        x=np.linspace(m.lo,m.hi,71)
        np.testing.assert_allclose(op.primitive(x),RW*D*(x-m.lo),rtol=2e-13,atol=1e-18)
    def test_sign_additivity_and_equal_states(self):
        op=self.model()._liquid_integral;a,b,c=np.linspace(op.x[0],op.x[-1],3)
        self.assertEqual(float(op.integral(b,b)),0.)
        self.assertAlmostEqual(float(op.integral(a,b)),-float(op.integral(b,a)),places=15)
        np.testing.assert_allclose(op.integral(a,c),op.integral(a,b)+op.integral(b,c),rtol=1e-13)
    def test_close_values_across_breakpoint(self):
        op=self.model()._liquid_integral;i=len(op.x)//2;x=op.x[i];eps=1e-12
        ref=op._segment(np.array([i-1]),np.array([x-eps]),np.array([x]))[0]
        ref+=op._segment(np.array([i]),np.array([x]),np.array([x+eps]))[0]
        np.testing.assert_allclose(op.integral(x-eps,x+eps),ref,rtol=1e-12)
    def test_nonmonotone_k_not_clipped_to_endpoint_values(self):
        o=copy.deepcopy(self.obj);o['functions']={'pC(Theta_l)_de':[[.001,.17,.34],[8.,5.,2.]],
            'lgKl(Theta_l)':[[.001,.17,.34],[-15.,-8.,-15.]],'lgKv(Theta_l)':[[.001,.34],[-30.,-30.]]}
        o['theta_effective']=.34;o['porosity']=.38
        m=Material(o);op=LiquidFaceIntegral(m);th=np.array([.001,.34]);p,_,k,_,_=m.fields(th,np.full(2,25.))
        self.assertGreater(op.face_k(th,p,k)[0],1e-14)
    def test_zero_pressure_gradient_returns_local_k(self):
        m=self.model().m;op=LiquidFaceIntegral(m);th=np.full(5,.10);p,_,k,_,_=m.fields(th,np.full(5,25.))
        np.testing.assert_allclose(op.face_k(th,p,k),k[:-1],rtol=1e-14)
    def test_constant_k_mean_is_exact(self):
        o=copy.deepcopy(self.obj);o['functions']['lgKl(Theta_l)']=[[.001,.34],[-12.,-12.]]
        m=Material(o);op=LiquidFaceIntegral(m);th=np.array([.01,.08,.12,.3]);p,_,k,_,_=m.fields(th,np.full(4,25.))
        np.testing.assert_allclose(op.face_k(th,p,k),1e-12,rtol=1e-12)
    def test_stationary_nonlinear_diffusion_face_flux(self):
        m=self.model().m;op=LiquidFaceIntegral(m)
        # Invert a linear Kirchhoff potential by bisection, not by re-fitting a law.
        W=np.linspace(float(op.primitive(.02)),float(op.primitive(.29)),22);lo=np.full(22,m.lo);hi=np.full(22,m.hi)
        for _ in range(60):
            mid=(lo+hi)/2;v=op.primitive(mid);lo=np.where(v<W,mid,lo);hi=np.where(v>=W,mid,hi)
        th=(lo+hi)/2;p,_,k,_,_=m.fields(th,np.full(22,25.))
        q=-op.face_k(th,p,k)*np.diff(p)/.03/RW
        np.testing.assert_allclose(q,-np.diff(W)/.03/RW,rtol=1e-11,atol=1e-18)
    def test_hydrostatic_pressure_gradient_cancels(self):
        m=self.linear_diffusivity_material();op=LiquidFaceIntegral(m);dist=.01
        p=-np.linspace(1e4,1e4+RW*9.80665*dist*9,10)
        xy=m.f[m.store];th=np.interp(np.log10(-p),xy[1][::-1],xy[0][::-1]);_,_,k,_,_=m.fields(th,np.full(10,25.))
        q=-op.face_k(th,p,k)*(np.diff(p)/dist+RW*9.80665)/RW
        np.testing.assert_allclose(q,0.,atol=1e-18)
    def test_closed_diffusion_matches_discrete_eigenmode(self):
        cfg=copy.deepcopy(self.cfg);cfg['_closed_test']=True;cfg['_zero_gravity_test']=True;cfg['geometry']['cells']=20
        D=1e-7;m=ColumnModel(cfg,self.linear_diffusivity_material(D));st=State(.17+.01*np.cos(np.pi*m.z/m.H),np.full(m.n,25.),np.zeros(m.n))
        out,h,_,s=m.integrate(self.climate(),st)
        rate=4*D*np.sin(np.pi/(2*m.n))**2/m.dz**2
        exact=.17+.01*np.cos(np.pi*m.z/m.H)*np.exp(-rate*24*3600)
        np.testing.assert_allclose(out.theta,exact,rtol=1e-6,atol=2e-8)
        self.assertLess(s['water_max_closure_error_kg_m2'],1e-7)
    def test_new_scheme_closed_uniform_stationarity(self):
        self.cfg['_closed_test']=True;self.cfg['_zero_gravity_test']=True;m=self.model();st=m.initial();st.sulfate[:]=.03
        out,_,_,_=m.integrate(self.climate(4),st)
        np.testing.assert_allclose(out.matrix(),st.matrix(),rtol=1e-10,atol=1e-10)
    def test_new_scheme_severe_wetting_bounds_and_balances(self):
        m=self.model();st=m.initial();st.theta[:]=m.m.hi-.00005;st.sulfate[:]=.2
        out,h,_,s=m.integrate(self.climate(6,True),st)
        self.assertLessEqual(out.theta.max(),m.m.hi+1e-9)
        self.assertGreaterEqual(out.theta.min(),m.m.lo-1e-9)
        self.assertLess(s['water_max_closure_error_kg_m2'],1e-7)
        self.assertLess(s['sulfate_max_closure_error_mol_m2'],1e-9)
        self.assertGreater(s['seepage_water_kg_m2'],0.)
    def test_new_scheme_jacobian_matches_columns(self):
        self.cfg['geometry']['cells']=10;m=self.model();st=m.initial();st.theta=np.linspace(.04,.2,m.n);st.temperature=np.linspace(20,24,m.n);st.sulfate=np.linspace(.001,.02,m.n)
        y=np.r_[st.matrix().ravel(),np.zeros(m.nl)];forc=m.forcing_array(self.climate(2).data)
        j=m.jacobian(0,y,forc).toarray();base=m.derivatives(0,y,forc)
        for k in range(m.size):
            step=np.sqrt(np.finfo(float).eps)*max(abs(y[k]),[.01,1,.001][k%3]);z=y.copy();z[k]+=step
            ref=(m.derivatives(0,z,forc)-base)/step
            np.testing.assert_allclose(j[:,k],ref,rtol=1e-3,atol=1e-6)
    def test_scheme_included_in_scientific_hash(self):
        p=Project(ROOT/'projects/demo');a=p.physics_hash();p.config['numerics']['liquid_face_scheme']='source_integral'
        self.assertNotEqual(a,p.physics_hash())
    def test_unknown_scheme_rejected(self):
        self.cfg['numerics']['liquid_face_scheme']='invented'
        with self.assertRaises(InputError):self.model()
    def test_no_unsupported_RH_conversion(self):
        o=copy.deepcopy(self.obj);o['functions']={k:v for k,v in o['functions'].items() if 'pC' not in k}
        o['functions']['Theta_l(RH)_de']=[[.01,.5,.99],[.001,.15,.34]]
        with self.assertRaises(InputError):LiquidFaceIntegral(Material(o))
    def test_real_source_unchanged(self):
        path=ROOT/'projects/mohenjo_daro/materials/18_15_IBEET_3.m6'
        if not path.exists():self.skipTest('Original material source not bundled.')
        before=sha256(path);m=read_material(path);op=LiquidFaceIntegral(m)
        op.primitive(np.linspace(m.lo,m.hi,1000));self.assertEqual(before,sha256(path))
        self.assertEqual(before,'b3280e6bfd9fd0baa2db335c18fac10cbb2d34a93e832ebc057efa2dabff0f48')

if __name__=='__main__':unittest.main()
