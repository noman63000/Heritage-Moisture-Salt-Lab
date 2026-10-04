"""Regression tests added after reproducing the user's real 400-cell failure.
These are implementation tests, not field validation or a full climate ensemble.
"""
import copy,json,tempfile,unittest
from pathlib import Path
import numpy as np
import pandas as pd
from hmsl.solver import ColumnModel,State,RW,LEDGERS
from hmsl.material import Material,read_material
from hmsl.climate import read_climate
from hmsl.project import Project
from hmsl.common import InputError,Cancelled
from hmsl.workflow import spinup
ROOT=Path(__file__).resolve().parents[1]

class RepairTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.path=Path(self.tmp.name)
        self.cfg=json.loads((ROOT/'projects/demo/project.json').read_text())
        self.obj=json.loads((ROOT/'projects/demo/materials/SYNTHETIC.json').read_text())
    def tearDown(self):self.tmp.cleanup()
    def climate(self,n=24,**kw):
        n=max(2,n)
        t=pd.date_range('2020-01-01',periods=n,freq='h')
        d=pd.DataFrame(dict(date=t.strftime('%Y-%m-%d'),hour_utc=t.hour,tas_C=26.498,
          hurs_pct=94.505,pr_mm_h=8.27748,rsds_W_m2=0.,rlds_W_m2=444.718,sfcWind_m_s=1.842))
        for k,v in kw.items():d[k]=v
        p=self.path/'weather.csv';d.to_csv(p,index=False);return read_climate(p)
    def model(self,n=16):
        self.cfg['geometry']['cells']=n
        return ColumnModel(self.cfg,Material(copy.deepcopy(self.obj)))
    def test_wet_endpoint_derivative_is_not_outward(self):
        m=self.model(40);st=m.initial();st.theta[:]=m.m.hi;st.sulfate[:]=.2
        y=np.r_[st.matrix().ravel(),np.zeros(m.nl)]
        rhs=m.derivatives(0,y,m.forcing_array(self.climate(1).data))
        self.assertLessEqual(rhs[:m.size].reshape(m.n,3)[:,0].max(),1e-15)
        self.assertGreater(rhs[m.size+2],0.)
        self.assertGreater(rhs[m.size+3],0.)
        self.assertAlmostEqual(np.sum(rhs[:m.size].reshape(m.n,3)[:,0]*m.widths)*RW,rhs[m.size],places=10)
    def test_wetting_from_below_stays_in_range_and_conserves(self):
        m=self.model(16);st=m.initial();st.theta[:]=m.m.hi-.00005;st.sulfate[:]=.2
        final,h,_,s=m.integrate(self.climate(12),st)
        self.assertLessEqual(h.theta_max_m3_m3.max(),m.m.hi+1e-9)
        self.assertGreaterEqual(final.theta.min(),m.m.lo-1e-9)
        self.assertLess(s['water_max_closure_error_kg_m2'],1e-7)
        self.assertLess(s['sulfate_max_closure_error_mol_m2'],1e-9)
        self.assertGreater(s['seepage_water_kg_m2'],0)
        self.assertGreater(s['seepage_salt_mol_m2'],0)
        np.testing.assert_allclose(h.net_water_kg_m2_hour,h.basal_water_kg_m2_hour+h.accepted_rain_kg_m2_hour-h.net_evaporation_kg_m2_hour-h.seepage_water_kg_m2_hour,atol=1e-8)
        np.testing.assert_allclose(h.net_salt_mol_m2_hour,h.basal_salt_mol_m2_hour-h.seepage_salt_mol_m2_hour,atol=1e-10)
    def test_rain_ledger_accounts_for_all_offered_rain(self):
        m=self.model(16);c=self.climate(12)
        _,h,_,_=m.integrate(c)
        offered=(m.H*m.a*m.surface['rain_catch_fraction']+float(m.surface['crown_rain']))*8.27748
        np.testing.assert_allclose(h.accepted_rain_kg_m2_hour+h.unabsorbed_rain_kg_m2_hour,offered,rtol=1e-7,atol=1e-8)
    def test_closed_column_has_no_invented_seepage(self):
        self.cfg['_closed_test']=True;self.cfg['_zero_gravity_test']=True
        m=self.model();st=m.initial();st.theta[:]=m.m.hi;st.sulfate[:]=.2
        f,h,_,s=m.integrate(self.climate(2),st)
        np.testing.assert_allclose(f.matrix(),st.matrix(),atol=1e-10)
        self.assertEqual(s['seepage_water_kg_m2'],0.)
    def test_truncated_material_not_called_saturation(self):
        self.obj['theta_effective']=self.obj['theta_effective']*1.02
        with self.assertRaisesRegex(InputError,'coverage ends'):
            self.model()
    def test_bad_initial_moisture_is_not_clipped(self):
        m=self.model();st=m.initial();st.theta[-1]=m.m.hi+.001
        with self.assertRaisesRegex(InputError,'Initial/checkpoint'):m.integrate(self.climate(1),st)
    def test_water_only_does_not_remove_or_supply_salt(self):
        m=self.model();st=m.initial();st.theta[:]=m.m.hi-.00001;st.sulfate[:]=.2
        f,h,_,_=m.integrate(self.climate(2),st,water_only=True)
        np.testing.assert_allclose(f.sulfate,st.sulfate,atol=1e-10)
        self.assertLess(h.net_salt_mol_m2_hour.abs().max(),1e-10)
    def test_custom_sparse_jacobian_matches_single_column_differences(self):
        m=self.model(10);st=m.initial();st.theta=np.linspace(.06,.2,m.n)
        st.temperature=np.linspace(20,24,m.n);st.sulfate=np.linspace(.001,.01,m.n)
        y=np.r_[st.matrix().ravel(),np.zeros(m.nl)];forc=m.forcing_array(self.climate(1).data)
        j=m.jacobian(0,y,forc).toarray();base=m.derivatives(0,y,forc)
        for k in range(m.size):
            step=np.sqrt(np.finfo(float).eps)*max(abs(y[k]),[.01,1,.001][k%3]);z=y.copy();z[k]+=step
            exact=(m.derivatives(0,z,forc)-base)/step
            np.testing.assert_allclose(j[:,k],exact,rtol=1e-3,atol=1e-6)
    def test_nonuniform_grid_flux_cancellation(self):
        self.cfg['geometry']['mesh_scheme']='boundary_refined';m=self.model(30)
        st=m.initial();st.theta=np.linspace(.04,.2,m.n);st.sulfate=np.linspace(.001,.02,m.n)
        y=np.r_[st.matrix().ravel(),np.zeros(m.nl)];d=m.derivatives(0,y,m.forcing_array(self.climate(1).data))
        dx=d[:m.size].reshape(m.n,3)
        self.assertAlmostEqual(np.dot(dx[:,0],m.widths)*RW,d[m.size],places=9)
        self.assertAlmostEqual(np.dot(dx[:,2],m.widths),d[m.size+1],places=12)
    def test_discrete_water_diffusion_eigenmode(self):
        self.cfg['_closed_test']=True;self.cfg['_zero_gravity_test']=True;m=self.model(20);D=1e-7
        def fields(th,T):
            one=np.ones_like(th);return 1e4*(th-.2),one*.5,one*RW*D/1e4,one*1e-30,one*.6
        m.m.fields=fields
        st=State(.2+.01*np.cos(np.pi*m.z/m.H),np.full(m.n,25.),np.zeros(m.n))
        f,_,_,_=m.integrate(self.climate(24),st)
        decay=4*D*np.sin(np.pi/(2*m.n))**2/m.dz**2
        exact=.2+.01*np.cos(np.pi*m.z/m.H)*np.exp(-decay*24*3600)
        np.testing.assert_allclose(f.theta,exact,rtol=1e-6,atol=1e-8)
    def test_discrete_heat_diffusion_eigenmode(self):
        self.cfg['_closed_test']=True;self.cfg['_zero_gravity_test']=True;m=self.model(20)
        def fields(th,T):
            a=np.ones_like(th);return a*0,a*.5,a*1e-30,a*1e-30,a*.6
        m.m.fields=fields
        st=State(np.full(m.n,.2),25.+np.cos(np.pi*m.z/m.H),np.zeros(m.n))
        f,_,_,_=m.integrate(self.climate(24),st)
        alpha=.6/(m.m.rho*m.m.cp+RW*4180*.2)
        rate=4*alpha*np.sin(np.pi/(2*m.n))**2/m.dz**2
        exact=25.+np.cos(np.pi*m.z/m.H)*np.exp(-rate*24*3600)
        np.testing.assert_allclose(f.temperature,exact,rtol=1e-6,atol=1e-6)
    def test_review_required_for_real_site_boundary_change(self):
        src=ROOT/'projects/mohenjo_daro'
        p=Project(src);p.config['review']['wet_end_closure_approved']=False
        self.assertFalse(p.readiness()['trial_ready'])
    def test_source_arrays_and_bytes_are_not_overwritten(self):
        from hmsl.common import sha256
        src=ROOT/'projects/mohenjo_daro/materials/18_15_IBEET_3.m6'
        if not src.exists():self.skipTest('Real source not present in this installation.')
        before=sha256(src);m=read_material(src)
        m.fields(np.array([m.lo,m.hi]),np.array([25.,25.]))
        self.assertEqual(before,sha256(src));self.assertAlmostEqual(m.hi,.349795)
        self.assertIsInstance(m.obj['functions']['lgKl(Theta_l)'],list)
    def test_spinup_cancellation_resume_preserves_chunks(self):
        import shutil
        src=ROOT/'projects/demo';target=self.path/'spinup_demo'
        shutil.copytree(src,target,ignore=shutil.ignore_patterns('outputs','*.npz','qa.json','spinup.json','spinup_cache','qa_cache','last_job.json'))
        p=Project(target);p.config['_closed_test']=True;p.config['_zero_gravity_test']=True;p.save(p.config)
        flag=[False]
        def progress(fr,msg):
            if '24 /' in msg:flag[0]=True
        with self.assertRaises(Cancelled):spinup(p,'DEMO',progress=progress,cancel=lambda:flag[0])
        cps=list((target/'spinup_cache').glob('*/checkpoint.json'));self.assertEqual(len(cps),1)
        checkpoint=json.loads(cps[0].read_text());self.assertEqual(checkpoint['completed_hours'],24)
        result=spinup(p,'DEMO');self.assertTrue(result['passed']);self.assertEqual(len(result['records']),1)

    def test_material_freezing_is_blocked_even_when_air_is_not_negative(self):
        self.cfg['surface']['sky_view_factor']=1.
        m=self.model();st=m.initial();st.temperature[:]=.1
        with self.assertRaisesRegex(InputError,'material cell fell below'):
            m.integrate(self.climate(2,tas_C=0.,hurs_pct=50.,rlds_W_m2=0.,pr_mm_h=0.),st)
    def test_full_case_preflight_does_not_clip_cold_data(self):
        from types import SimpleNamespace
        from hmsl.preflight import check_cases,require_supported
        c=self.climate(2,tas_C=-1.)
        p=SimpleNamespace(path=self.path,config={'runs':[{'id':'TEST'}]},climate=lambda rid:c)
        r=check_cases(p)
        self.assertFalse(r['all_selected_supported']);self.assertEqual(r['cases'][0]['subzero_air_hours'],2)
        with self.assertRaisesRegex(InputError,'before starting'):require_supported(p,['TEST'])
        self.assertEqual(float(c.data.tas_C.min()),-1.)
    def test_cancel_inside_wet_hour_is_safe(self):
        m=self.model();st=m.initial();st.theta[:]=m.m.hi-.000001
        with self.assertRaises(Cancelled):m.integrate(self.climate(2),st,cancel=lambda:True)

if __name__=='__main__':unittest.main(verbosity=2)
