"""Reduced 1-D vertical finite-volume water/heat/sulfate transport solver.

See docs/METHODS.md for equations, units and limitations. This is not a calibrated
Mohenjo-daro model, DELPHIN replacement, mixed-ion model, or damage predictor.
State per cell: volumetric liquid water, Celsius temperature, TOTAL Na2SO4 formula
moles / m3 of bulk material. Nine extra states independently integrate external
fluxes, including seepage and its dissolved salt. Wet-bound events activate an
explicit exposed-face drainage closure; no endpoint mass clipping is used.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix, csc_matrix, vstack, hstack
from .common import InputError, Cancelled
from .material import RW, RV, Material
LEDGERS=('net_water_kg_m2','net_salt_mol_m2','seepage_water_kg_m2','seepage_salt_mol_m2',
         'unabsorbed_rain_kg_m2','accepted_rain_kg_m2','net_evaporation_kg_m2','basal_water_kg_m2','basal_salt_mol_m2')
G=9.80665; LV=2.45e6; CW=4180.; SIGMA=5.670374419e-8

class PhysicsScopeError(InputError):
    """A physical mechanism is missing; reducing the step is not a remedy."""
    pass

def psat(T):
    return 610.94*np.exp(17.625*np.asarray(T)/(np.asarray(T)+243.04))

def harmonic(a,b): return 2*a*b/np.maximum(a+b,1e-100)

@dataclass
class State:
    theta: np.ndarray
    temperature: np.ndarray
    sulfate: np.ndarray
    @classmethod
    def from_matrix(cls,x): return cls(x[:,0].copy(),x[:,1].copy(),x[:,2].copy())
    def matrix(self): return np.column_stack((self.theta,self.temperature,self.sulfate))

class Solubility:
    """Optional, user-sourced SINGLE SALT equilibrium partition in mol/kg water.
    No metastability, activity model, hydrate-water accounting or damage law.
    """
    def __init__(self,path):
        d=pd.read_csv(path)
        required={'T_C','Na2SO4_mol_kg_water','source'}
        if not required.issubset(d): raise InputError('Solubility file needs T_C, Na2SO4_mol_kg_water, source.')
        x=d.T_C.to_numpy(float); y=d.Na2SO4_mol_kg_water.to_numpy(float)
        if len(x)<2 or not np.isfinite(x).all() or not np.isfinite(y).all() or not np.all(np.diff(x)>0) or np.any(y<=0):
            raise InputError('Solubility must be positive finite values at strictly increasing temperatures.')
        if d.source.isna().any() or d.source.astype(str).str.strip().eq('').any(): raise InputError('Every solubility row needs a source citation.')
        self.x=x; self.y=y
    def concentration(self,T):
        if np.any(T<self.x[0]-1e-6) or np.any(T>self.x[-1]+1e-6):
            raise InputError('Material temperature left the supplied solubility table. Extend it with verified data; no extrapolation is applied.')
        return 1000.*np.interp(T,self.x,self.y)

class ColumnModel:
    def __init__(self, cfg: dict, material: Material, solubility=None):
        self.cfg=cfg; self.m=material; self.solubility=solubility
        g=cfg['geometry']; self.n=int(g['cells']); self.H=float(g['height_m']); self.b=float(g['thickness_m'])
        if self.n<4 or self.n>400 or self.H<=0 or self.b<=0: raise InputError('Use 4-400 cells and positive height/thickness.')
        self.mesh_scheme=g.get('mesh_scheme','uniform')
        eta=np.linspace(0.,1.,self.n+1)
        if self.mesh_scheme=='boundary_refined':
            # Fixed transformation across all QA meshes. Physical endpoints and
            # integration-band boundaries do not move when resolution changes.
            eta=(1.+np.tanh(3.*(2.*eta-1.))/np.tanh(3.))/2.
        elif self.mesh_scheme!='uniform': raise InputError('Unknown mesh scheme.')
        self.edges=self.H*eta; self.widths=np.diff(self.edges)
        self.z=(self.edges[:-1]+self.edges[1:])/2.
        self.dist=np.diff(self.z)
        self.dz=self.H/self.n  # legacy nominal spacing, not used in FV balances
        self.size=3*self.n; self.nl=len(LEDGERS)
        self.faces=float(g['exposed_faces']); self.a=self.faces/self.b
        self.surface=cfg['surface']; self.salt=cfg['salt']
        if self.salt['diffusion_convention'] not in ('pore','bulk'):
            raise InputError('Select and document the sulfate diffusion convention: pore or bulk.')
        self.closed=bool(cfg.get('_closed_test',False))
        self.gravity=0. if cfg.get('_zero_gravity_test',False) else G
        n=self.size+self.nl; j=lil_matrix((n,n),dtype=int)
        for k in range(self.n):
            for kk in range(max(0,k-1),min(self.n,k+2)):
                j[k*3:k*3+3,kk*3:kk*3+3]=1
        j[self.size:,:self.size]=1
        self.jac=j.tocsr()
        self.atol=np.r_[np.tile([1e-10,1e-6,1e-11],self.n),[1e-7,1e-10,1e-7,1e-10,1e-7,1e-7,1e-7,1e-7,1e-10]]
        self._last=None
        self.liquid_face_scheme=cfg.get('numerics',{}).get('liquid_face_scheme','harmonic_legacy')
        if self.liquid_face_scheme not in ('harmonic_legacy','source_integral'):
            raise InputError('Unknown liquid face discretisation.')
        self._liquid_integral=None
        if self.liquid_face_scheme=='source_integral':
            from .face_integral import LiquidFaceIntegral
            self._liquid_integral=LiquidFaceIntegral(self.m)
        self._physical_pattern=self.jac[:self.size,:self.size].tocoo()
        # Physical rows couple only nearest neighbouring cells. The two global
        # flux ledger rows must NOT turn sparse finite differences into N groups.
        # Their Jacobian is derived from the exact conservative row identities;
        # their RHS values remain independently computed from boundary fluxes.
        self._jrow=self._physical_pattern.row
        self._jcol=self._physical_pattern.col
        self._colours=[np.flatnonzero(np.arange(3*self.n)%9==k) for k in range(9)]
        self._capacity_gap=(self.m.sat-self.m.hi)/self.m.sat
        self._events_total=0; self._steps_total=0; self._min_internal_T=float("inf")
        if self._capacity_gap>0.001 and not self.closed:
            raise InputError('Liquid/storage source coverage ends more than 0.1% below effective saturation. Supply complete source functions; no drainage-at-truncated-table assumption is made.')
    def initial(self):
        ini=self.cfg['initial']; T=np.full(self.n,float(ini['temperature_C']))
        th=np.full(self.n,float(self.m.theta_at_rh(float(ini['RH_pct'])/100.,T[0])))
        if np.any(th<=self.m.lo) or np.any(th>=self.m.hi): raise InputError('Initial RH maps outside the shared material-function range. Choose an in-range RH.')
        return State(th,T,np.full(self.n,float(self.salt['initial_sulfate_mol_m3_bulk'])))
    def partition(self,theta,T,total):
        # Clamp only intermediate coefficient evaluation. Accepted states are
        # checked separately and are never mass-clipped.
        th=np.maximum(theta,self.m.lo)
        c=np.maximum(total,0)/th
        if self.solubility is not None: c=np.minimum(c,self.solubility.concentration(T))
        solid=total-th*c
        return c,np.maximum(solid,0)
    def derivatives(self,t,y,forcing,water_only=False):
        X=y[:self.size].reshape(self.n,3); th,T,ns=X.T
        if np.any(T<-50) or np.any(T>100): raise InputError('Solver temperature left -50 to 100 C. Review boundary fluxes and numerical settings.')
        i=min(int(max(t,0)//3600),len(forcing)-1)
        Ta,RHa,rain,sw,lw,wind,supply,csource=forcing[i]
        p,rh,kl,kv,lam=self.m.fields(th,T)
        pv=rh*psat(T); pvair=RHa/100.*psat(Ta)
        ql=np.zeros(self.n+1); qv=np.zeros(self.n+1)
        face_kl=self.face_k(kl) if self._liquid_integral is None else self._liquid_integral.face_k(th,p,kl)
        ql[1:-1]=-face_kl/RW*((p[1:]-p[:-1])/self.dist+RW*self.gravity)
        qv[1:-1]=-self.face_k(kv)*(pv[1:]-pv[:-1])/self.dist/RW
        limit=np.clip(1.-(np.maximum(th,0)/self.m.hi)**4,0,1)
        ql[0]=supply/1000./86400.*limit[0]
        # Supplied flux is PARAMETRIC; no groundwater-head/soil inference.
        exposure=np.full(self.n,self.a)
        exposure[-1]+=1./self.widths[-1]
        hc=5.7+3.8*wind
        beta=hc/(1.2*1005.*RV*(Ta+273.15))  # Lewis number 1, declared reference air density/cp
        conduct=kv+kl*RW*RV*(T+273.15)/np.maximum(pv,1.)
        beta_eff=1./(1./beta+(self.b/4.)/np.maximum(conduct,1e-30))
        evapor=exposure*beta_eff*(pv-pvair)  # kg/m3/s, negative = condensation
        rainfall=np.full(self.n,self.a*self.surface['rain_catch_fraction']*rain/1000./3600.)
        if self.surface['crown_rain']: rainfall[-1]+=rain/1000./3600./self.widths[-1]
        captured=rainfall*limit
        # Do not remove water that is already at the supported dry endpoint.
        # This changes fluxes, NOT stored mass. Internal face fluxes remain paired.
        dry=th<=self.m.lo
        ql[1:-1]=np.where(((ql[1:-1]>0)&dry[:-1])|((ql[1:-1]<0)&dry[1:]),0.,ql[1:-1])
        qv[1:-1]=np.where(((qv[1:-1]>0)&dry[:-1])|((qv[1:-1]<0)&dry[1:]),0.,qv[1:-1])
        if self.closed:
            ql[0]=0; evapor[:]=0; captured[:]=0; exposure[:]=0
        dtheta=(ql[:-1]-ql[1:]+qv[:-1]-qv[1:])/self.widths + captured-evapor/RW
        # Complementarity closure: at the source-supported wet endpoint, positive
        # net influx leaves an EXPOSED face as liquid seepage. No state clipping,
        # material extrapolation, or hidden water deletion is used. This is an
        # explicit reduced free-draining-face assumption, NOT a saturated pressure
        # solve or an intrinsic measured material law. Endpoint gap is reported.
        self._raw_dtheta=dtheta.copy()
        wet=getattr(self,'_active_wet',None)
        if wet is None: wet=(th>=self.m.hi)
        seep=np.where(wet & (exposure>0) & (not self.closed),np.maximum(dtheta,0.),0.)
        dtheta-=seep
        # At the dry endpoint, reduce evaporation (and its latent heat) rather
        # than manufacture water. Remaining unsupported drying will still fail.
        prevented=np.where(dry,np.minimum(np.maximum(-dtheta,0.),np.maximum(evapor/RW,0.)),0.)
        evapor-=RW*prevented
        dtheta+=prevented
        # Reference-temperature sensible heat equation. Heat carried by water
        # advection/hydrate transitions is deliberately outside this reduced model.
        qheat=np.zeros(self.n+1)
        qheat[1:-1]=-self.face_k(lam)*(T[1:]-T[:-1])/self.dist
        sky=self.surface['sky_view_factor']
        effective_lw=sky*lw+(1-sky)*SIGMA*(Ta+273.15)**4
        absorbed=self.surface['absorptivity']*self.surface['solar_projection_factor']*sw
        net=hc*(Ta-T)+absorbed+self.surface['emissivity']*(effective_lw-SIGMA*(T+273.15)**4)
        energy=(qheat[:-1]-qheat[1:])/self.widths + exposure*net-LV*evapor
        dT=energy/(self.m.rho*self.m.cp+RW*CW*np.maximum(th,0))
        c,solid=self.partition(th,T,ns)
        salt_flux=np.zeros(self.n+1)
        face_theta=np.maximum((th[1:]+th[:-1])/2.,0)
        sat_fraction=np.clip(face_theta/self.m.sat,0,1)
        D=float(self.salt['D_sulfate_m2_s'])*sat_fraction**float(self.salt['unsaturated_exponent'])
        if self.salt['diffusion_convention']=='pore': D=D*face_theta
        salt_flux[1:-1]=np.where(ql[1:-1]>=0,c[:-1],c[1:])*ql[1:-1]-D*np.diff(c)/self.dist
        salt_flux[0]=ql[0]*(csource/96.06) # mg SO4/L = g SO4/m3; /96.06 g/mol
        if water_only: salt_flux[:]=0
        salt_loss=np.zeros(self.n) if water_only else seep*c
        dns=(salt_flux[:-1]-salt_flux[1:])/self.widths-salt_loss
        # Each entry below is an independently evaluated external exchange,
        # never inferred from changes in the stored water/salt inventory.
        local=np.zeros((self.n,self.nl))
        local[:,2]=RW*seep; local[:,3]=salt_loss
        local[:,4]=RW*(rainfall-captured) if not self.closed else 0.
        local[:,5]=RW*captured; local[:,6]=evapor
        local[0,7]=RW*ql[0]/self.widths[0]
        local[0,8]=salt_flux[0]/self.widths[0]
        local[:,0]=local[:,7]+local[:,5]-local[:,6]-local[:,2]
        local[:,1]=local[:,8]-local[:,3]
        self._ledger_local=local
        return np.r_[np.column_stack((dtheta,dT,dns)).ravel(),np.sum(local*self.widths[:,None],axis=0)]

    def face_k(self,k):
        """Resistance-weighted harmonic face conductivity on a nonuniform grid."""
        return (self.widths[:-1]+self.widths[1:])/(self.widths[:-1]/k[:-1]+self.widths[1:]/k[1:])

    def jacobian(self,t,y,forcing,water_only=False):
        base=self.derivatives(t,y,forcing,water_only)
        base_local=self._ledger_local.copy()
        n=self.size
        scale=np.tile([0.01,1.,0.001],self.n)
        steps=np.sqrt(np.finfo(float).eps)*np.maximum(np.abs(y[:n]),scale)
        vals=np.zeros(len(self._jrow)); led=np.zeros((self.nl,n))
        for colour,cols in enumerate(self._colours):
            z=y.copy(); z[cols]+=steps[cols]
            diff=self.derivatives(t,z,forcing,water_only)[:n]-base[:n]
            sel=self._jcol%9==colour
            vals[sel]=diff[self._jrow[sel]]/steps[self._jcol[sel]]
            dl=(self._ledger_local-base_local)*self.widths[:,None]
            accum=dl.copy(); accum[1:]+=dl[:-1]; accum[:-1]+=dl[1:]
            led[:,cols]=(accum[cols//3]/steps[cols,None]).T
        physical=csc_matrix((vals,(self._jrow,self._jcol)),shape=(n,n))
        left=vstack([physical,csc_matrix(led)],format='csc')
        return hstack([left,csc_matrix((n+self.nl,self.nl))],format='csc')

    def solve_hour(self,y,forcing_one,water_only,maxstep,rtol,rhs):
        """Event-driven active wet-end constraint, not endpoint state clipping."""
        from types import SimpleNamespace
        z=y.copy(); t0=0.; events_count=0; steps=0
        self._active_wet=np.zeros(self.n,dtype=bool)
        self.derivatives(0.,z,forcing_one,water_only)
        raw=self._raw_dtheta.copy()
        eligible=np.full(self.n,self.a>0) if not self.closed else np.zeros(self.n,bool)
        if not self.closed: eligible[-1]=True
        self._active_wet=eligible & (z[:self.size].reshape(self.n,3)[:,0]>=self.m.hi-1e-12) & (raw>0.)
        while t0<3600.-1e-9:
            active=self._active_wet.copy(); inactive=(~active)&eligible
            def hit(t,x):
                return float(np.min(self.m.hi-x[:self.size].reshape(self.n,3)[inactive,0])) if inactive.any() else 1.
            hit.terminal=True; hit.direction=-1
            def release(t,x):
                if not active.any(): return 1.
                self.derivatives(t,x,forcing_one,water_only)
                return float(np.min(self._raw_dtheta[active]))
            release.terminal=True; release.direction=-1
            out=solve_ivp(rhs,(t0,3600.),z,method='BDF',
                first_step=min(30.,maxstep,3600.-t0),max_step=maxstep,rtol=rtol,
                atol=self.atol,jac=lambda t,x:self.jacobian(t,x,forcing_one,water_only),
                events=[hit,release])
            if not out.success: return out
            self._min_internal_T=min(self._min_internal_T,float(out.y[1:self.size:3].min()))
            if self._min_internal_T < -1e-7:
                raise PhysicsScopeError('A simulated material cell fell below 0 C. Ice/freezing physics is not implemented. The temperature was not clipped and the run is not accepted.')
            steps+=len(out.t)-1
            z=out.y[:,-1]; t0=float(out.t[-1])
            if out.status!=1:break
            events_count+=1
            if events_count>5*self.n+100:
                raise InputError('Wet-end active-set event limit reached. No unconverged hour accepted.')
            if len(out.t_events[0]):
                th=z[:self.size].reshape(self.n,3)[:,0]
                candidates=np.flatnonzero(inactive)
                which=candidates[np.argmin(self.m.hi-th[candidates])]
                self._active_wet[which]=True
            else:
                self.derivatives(t0,z,forcing_one,water_only)
                candidates=np.flatnonzero(active)
                which=candidates[np.argmin(self._raw_dtheta[candidates])]
                self._active_wet[which]=False
        self._active_wet=None
        self._events_total+=events_count; self._steps_total+=steps
        return SimpleNamespace(success=True,y=z[:,None],message='hour complete',
            event_switches=events_count,internal_steps=steps)

    def forcing_array(self,df):
        arr=df[['tas_C','hurs_pct','pr_mm_h','rsds_W_m2','rlds_W_m2','sfcWind_m_s']].to_numpy(float)
        b=self.cfg['boundary']
        s=df['capillary_supply_mm_day'].to_numpy(float) if 'capillary_supply_mm_day' in df else np.full(len(df),b['capillary_supply_mm_day'])
        c=df['sulfate_source_mg_L'].to_numpy(float) if 'sulfate_source_mg_L' in df else np.full(len(df),b['sulfate_source_mg_L'])
        return np.column_stack((arr,s,c))
    def integrate(self, climate, initial=None, cancel=None, progress=None, chunk_callback=None, water_only=False):
        """Integrate with forcing held constant over each climate hour.

        Earlier builds integrated a multi-hour chunk while switching the forcing
        index inside the RHS.  That gives the stiff BDF solver discontinuities at
        every hour boundary and can become fragile on fine meshes.  Here each
        climate hour is a separate IVP; physical state is continuous between
        hours, while the external mass-accounting states accumulate across the
        checkpoint chunk.  This also localises failures to an exact forcing hour.
        """
        if (climate.data.tas_C<0).any():
            raise InputError('Sub-zero climate in this case. Ice/freezing is not modelled; use climate screening or a freeze-thaw capable solver.')
        if initial is None: initial=self.initial()
        state=initial.matrix()
        if state.shape!=(self.n,3): raise InputError('Initial state and current mesh differ; regenerate the spin-up.')
        if (not np.isfinite(state).all() or state[:,0].min()<self.m.lo-1e-12
                or state[:,0].max()>self.m.hi+1e-12 or state[:,2].min()<-1e-12):
            raise InputError('Initial/checkpoint state is outside the supported moisture or nonnegative salt range. No state was clipped.')
        y=np.r_[state.ravel(),np.zeros(self.nl)]
        block=int(self.cfg['numerics']['chunk_hours']); rows=[]; profiles=[]; rejected=0
        maxstep=float(self.cfg['numerics']['max_step_s']); rtol=float(self.cfg['numerics']['rtol'])

        for start in range(0,climate.n,block):
            if cancel and cancel(): raise Cancelled('Run cancelled after the last saved chunk. Use resume to restart at the last checkpoint.')
            end=min(start+block,climate.n)
            frame=climate.data.iloc[start:end]
            forc_all=self.forcing_array(frame)

            # Rebase only the two accounting integrators at each committed chunk.
            y[self.size:]=0.
            current=y[:self.size].reshape(self.n,3)
            init_water=RW*np.sum(current[:,0]*self.widths)
            init_salt=np.sum(current[:,2]*self.widths)

            raw_hours=[]; ext_water_hours=[]; ext_salt_hours=[]; ledger_hours=[]
            for local in range(len(frame)):
                if cancel and cancel(): raise Cancelled('Run cancelled; completed checkpoint chunks are retained.')
                forcing_one=forc_all[local:local+1]
                absolute_row=start+local
                failure_notes=[]; result=None

                def rhs(t,z):
                    if cancel and cancel(): raise Cancelled('Run cancelled; completed checkpoint chunks are retained.')
                    return self.derivatives(t,z,forcing_one,water_only)

                # Several start-step/max-step fallbacks are numerical safeguards,
                # not changes to the physical forcing.  The BDF solver remains
                # adaptive within each one-hour, constant-forcing interval.
                for factor in (1., .25, .0625, .015625):
                    try:
                        trial=self.solve_hour(y,forcing_one,water_only,maxstep*factor,min(rtol,3e-8),rhs)
                    except (Cancelled,PhysicsScopeError):
                        raise
                    except Exception as exc:
                        rejected+=1
                        failure_notes.append(f'{factor:g}x step exception: {type(exc).__name__}: {exc}')
                        continue
                    if not trial.success:
                        rejected+=1
                        failure_notes.append(f'{factor:g}x step solver: {trial.message}')
                        continue
                    raw_one=trial.y[:self.size,-1].reshape(self.n,3)
                    th=raw_one[:,0]; ns=raw_one[:,2]
                    good=(np.isfinite(trial.y).all() and th.min()>=self.m.lo-1e-9 and
                          th.max()<=self.m.hi+1e-9 and ns.min()>=-1e-10)
                    if not good:
                        rejected+=1
                        failure_notes.append(
                            f'{factor:g}x step rejected state: theta [{th.min():.6g},{th.max():.6g}] '
                            f'allowed [{self.m.lo:.6g},{self.m.hi:.6g}], sulfate min {ns.min():.6g}')
                        continue
                    result=trial
                    break

                if result is None:
                    forcing_desc=(f'Tair={forcing_one[0,0]:.3f} C, RH={forcing_one[0,1]:.3f} %, '
                                  f'rain={forcing_one[0,2]:.6g} mm/h, SW={forcing_one[0,3]:.3f} W/m2, '
                                  f'LW={forcing_one[0,4]:.3f} W/m2, wind={forcing_one[0,5]:.3f} m/s, '
                                  f'supply={forcing_one[0,6]:.6g} mm/day, SO4={forcing_one[0,7]:.6g} mg/L')
                    detail=' | '.join(failure_notes[-4:]) if failure_notes else 'no solver diagnostic returned'
                    raise InputError(
                        f'Numerical integration failed at climate row {absolute_row} '
                        f'({frame.timestamp_utc.iloc[local]}), mesh={self.n} cells. {forcing_desc}. '
                        f'Material theta range [{self.m.lo:.6g},{self.m.hi:.6g}]. Solver diagnostics: {detail}')

                y=result.y[:,-1]
                raw_hours.append(y[:self.size].reshape(self.n,3).copy())
                ledger_hours.append(y[self.size:].copy())
                ext_water_hours.append(float(y[self.size])); ext_salt_hours.append(float(y[self.size+1]))
                if progress:
                    done=absolute_row+1
                    # Hourly progress is especially useful for fine-mesh QA and
                    # makes it clear that a long solve is advancing.
                    progress(done/climate.n,f'{done:,} / {climate.n:,} hours; mesh {self.n} cells')

            raw=np.asarray(raw_hours)
            ext_water_series=np.asarray(ext_water_hours)
            ext_salt_series=np.asarray(ext_salt_hours)
            water=RW*(raw[:,:,0]*self.widths).sum(axis=1)
            salt=(raw[:,:,2]*self.widths).sum(axis=1)
            ew=water-init_water-ext_water_series
            es=salt-init_salt-ext_salt_series
            if np.max(np.abs(ew))>max(1e-4,1e-5*max(init_water,1)) or np.max(np.abs(es))>max(1e-7,1e-5*max(init_salt,1)):
                raise InputError('Mass-closure test failed. Outputs are not accepted as a completed run.')

            increments=np.diff(np.vstack([np.zeros((1,self.nl)),np.asarray(ledger_hours)]),axis=0)
            for k in range(len(frame)):
                th,T,ns=raw[k].T; c,solid=self.partition(th,T,ns)
                band=0.10
                edges=np.arange(0.,self.H+band,band)
                if edges[-1]<self.H: edges=np.r_[edges,self.H]
                else: edges[-1]=self.H
                cell_lo=self.edges[:-1]; cell_hi=self.edges[1:]
                band_c=[]
                for lo,hi in zip(edges[:-1],edges[1:]):
                    overlap=np.maximum(0.,np.minimum(cell_hi,hi)-np.maximum(cell_lo,lo))
                    liquid=np.sum(th*overlap)
                    if liquid>1e-15:
                        band_c.append(float(np.sum(c*th*overlap)/liquid))
                c_band_max=max(band_c) if band_c else float(c.max())
                item={'forcing_timestamp_utc':str(frame.timestamp_utc.iloc[k]),'elapsed_hours':start+k+1,'closure_segment_start_hour':start,
                  'theta_min_m3_m3':float(th.min()),'T_min_C':float(T.min()),
                  'theta_mean_m3_m3':float(np.sum(th*self.widths)/self.H),'theta_max_m3_m3':float(th.max()),
                  'T_mean_C':float(np.sum(T*self.widths)/self.H),'T_max_C':float(T.max()),
                  'water_inventory_kg_m2_footprint':float(water[k]),'sulfate_total_mol_m2_footprint':float(salt[k]),
                  'sulfate_cmax_mol_m3_liquid':float(c.max()),
                  'sulfate_cmax_0p10m_bandavg_mol_m3_liquid':c_band_max,
                  'water_closure_error_kg_m2':float(ew[k]),
                  'sulfate_closure_error_mol_m2':float(es[k])}
                if self.solubility is not None: item['equilibrium_solid_mol_m2_footprint']=float(np.sum(solid*self.widths))
                item.update({name+'_hour':float(increments[k,j]) for j,name in enumerate(LEDGERS)})
                rows.append(item)

            profile=pd.DataFrame({'elapsed_hours':end,'height_m':self.z,'theta_m3_m3':raw[-1,:,0],
                'temperature_C':raw[-1,:,1],'sulfate_total_mol_m3_bulk':raw[-1,:,2]})
            profiles.append(profile)
            if chunk_callback: chunk_callback(end,State.from_matrix(raw[-1]),pd.DataFrame(rows[-len(frame):]),profile)
            if progress: progress(end/climate.n,f'{end:,} / {climate.n:,} hours; water closure {np.max(np.abs(ew)):.2e} kg/m2')

        report={'solver':'1-D finite-volume / SciPy BDF; hourly constant-forcing IVPs','cells':self.n,'liquid_face_scheme':self.liquid_face_scheme,'mesh_scheme':self.mesh_scheme,'smallest_cell_m':float(self.widths.min()),'largest_cell_m':float(self.widths.max()),'max_step_s':maxstep,'requested_rtol':rtol,'effective_rtol':min(rtol,3e-8),
          'completed_hours':climate.n,'rejected_hour_attempts':rejected,
          'water_max_closure_error_kg_m2':max(abs(r['water_closure_error_kg_m2']) for r in rows),
          'sulfate_max_closure_error_mol_m2':max(abs(r['sulfate_closure_error_mol_m2']) for r in rows),
          'salt_phase_mode':'optional single-salt equilibrium partition' if self.solubility else 'transport only; all salt treated as dissolved-equivalent; no phase/damage output',
          'material_name':self.m.name,'is_demo':self.m.is_demo,'boundary_closure':'source-supported exposed-face seepage','supported_wet_endpoint':self.m.hi,'effective_saturation':self.m.sat,'coverage_gap_fraction':self._capacity_gap,'wet_end_event_switches':self._events_total,'accepted_internal_steps':self._steps_total,'min_internal_material_temperature_C':self._min_internal_T,
          'seepage_water_kg_m2':sum(r['seepage_water_kg_m2_hour'] for r in rows),
          'seepage_salt_mol_m2':sum(r['seepage_salt_mol_m2_hour'] for r in rows),
          'unabsorbed_rain_kg_m2':sum(r['unabsorbed_rain_kg_m2_hour'] for r in rows)}
        return State.from_matrix(y[:self.size].reshape(self.n,3)),pd.DataFrame(rows),pd.concat(profiles,ignore_index=True),report
