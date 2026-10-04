"""Equilibrium pore-water freeze/thaw extension for the reduced Heritage Lab model.

The checkpoint water state remains three variables per cell.  In this extension
``theta`` is the *total condensable water mass expressed as liquid-water-equivalent
volume per bulk volume*.  At every RHS evaluation it is partitioned algebraically
into mobile liquid and immobile ice-equivalent water using the source retention
curve and the pure-water ice/liquid equilibrium vapour-pressure ratio.

This is an equilibrium phase-change extension, not a frost-damage model.  It adds
latent heat, ice sensible heat, liquid-only transport, and salt exclusion from ice.
It does not model nucleation/supercooling hysteresis, ice-pressure damage, pore
expansion, snow accumulation, or salt freezing-point depression.
"""
from __future__ import annotations
import math
from types import SimpleNamespace
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from hmsl.solver import ColumnModel, State, LEDGERS, PhysicsScopeError, psat, SIGMA, LV, CW
from hmsl.material import RW, RV
from hmsl.common import InputError, Cancelled
from .pore_ice import freezing_limit_rh

VERSION = 'equilibrium-freeze-thaw-1.0'
MIN_PHASE_C = -20.0
LF = 333550.0       # J/kg, water fusion latent heat near 0 C
CI = 2100.0         # J/(kg K), representative ice-Ih heat capacity near 0 C
ICE_EPS = 1e-12

SOURCES = [
    'https://doi.org/10.1016/S0360-1323(00)00066-4',
    'https://doi.org/10.1256/qj.04.94',
    'https://www.iapws.org/relguide/Ice-2009.html',
    'https://www.iapws.org/relguide/Supercooled.html',
]

class EquilibriumFreezeThawModel(ColumnModel):
    """Original 3-state FV model with an algebraic liquid/ice equilibrium split.

    Above the phase boundary the inherited derivative is called verbatim, so the
    warm/liquid calculation follows the original arithmetic path.  When ice is
    present, theta is total-water equivalent; flux coefficients and salt
    concentration use only the algebraically available liquid fraction.
    """
    def __init__(self,*a,**kw):
        super().__init__(*a,**kw)
        # Dense, deterministic equilibrium table.  No fitted material constant is
        # introduced: the only material input is the already frozen retention curve.
        self._phase_T = np.linspace(MIN_PHASE_C, 0.0, 8001)
        rh = freezing_limit_rh(self._phase_T)
        cap = np.asarray(self.m.theta_at_rh(rh, self._phase_T), float)
        cap = np.clip(cap, self.m.lo, self.m.hi)
        # A physical retention curve must give nondecreasing unfrozen capacity with T.
        if np.any(np.diff(cap) < -1e-8):
            raise InputError('Source retention curve gives a non-monotone equilibrium unfrozen-water capacity. Phase extension was not started.')
        cap = np.maximum.accumulate(cap)
        self._phase_cap = cap
        self._phase_dcap = np.maximum(np.gradient(cap, self._phase_T), 0.0)
        self._phase_max_ice = 0.0
        self._phase_max_ice_fraction = 0.0
        self._phase_hours = 0
        self._freeze_thaw_switches = 0
        self._prev_hour_has_ice = False
        self._hour_diag = None

    def unfrozen_capacity(self,T):
        t=np.asarray(T,float)
        if np.any(t < MIN_PHASE_C-1e-9):
            raise PhysicsScopeError(f'Material temperature fell below {MIN_PHASE_C:g} C, outside the declared equilibrium phase-extension range.')
        return np.where(t>=0.0, self.m.hi, np.interp(t,self._phase_T,self._phase_cap))

    def capacity_dT(self,T):
        t=np.asarray(T,float)
        if np.any(t < MIN_PHASE_C-1e-9):
            raise PhysicsScopeError(f'Material temperature fell below {MIN_PHASE_C:g} C, outside the declared equilibrium phase-extension range.')
        return np.where(t>=0.0, 0.0, np.interp(t,self._phase_T,self._phase_dcap))

    def phase_partition(self,theta_total,T):
        total=np.asarray(theta_total,float); temp=np.asarray(T,float)
        cap=self.unfrozen_capacity(temp)
        liquid=np.minimum(total,cap)
        # State validation keeps total in range.  This clip is only a material
        # coefficient-evaluation guard, analogous to the original solver.
        liquid_eval=np.clip(liquid,self.m.lo,self.m.hi)
        ice=np.maximum(total-liquid,0.0)
        return liquid_eval,ice,cap

    def _set_hour_diag_from_solution(self,out):
        if self._hour_diag is None:
            self._hour_diag={'max_ice_theta':0.0,'max_ice_fraction':0.0,'min_liquid_theta':float('inf'),'phase_internal_samples':0}
        for k in range(out.y.shape[1]):
            X=out.y[:self.size,k].reshape(self.n,3); total,T,_=X.T
            liquid,ice,_=self.phase_partition(total,T)
            frac=np.divide(ice,np.maximum(total,1e-30))
            self._hour_diag['max_ice_theta']=max(self._hour_diag['max_ice_theta'],float(ice.max()))
            self._hour_diag['max_ice_fraction']=max(self._hour_diag['max_ice_fraction'],float(frac.max()))
            self._hour_diag['min_liquid_theta']=min(self._hour_diag['min_liquid_theta'],float(liquid.min()))
            self._hour_diag['phase_internal_samples']+=1

    def derivatives(self,t,y,forcing,water_only=False):
        X=y[:self.size].reshape(self.n,3); total,T,ns=X.T
        if np.any(T<-50) or np.any(T>100):
            raise InputError('Solver temperature left -50 to 100 C. Review boundary fluxes and numerical settings.')
        liquid,ice,cap=self.phase_partition(total,T)
        # Exact legacy path whenever no phase is present.  This is important for
        # the already verified 17 non-freezing cases and warm segments.
        if not np.any(ice>ICE_EPS):
            return super().derivatives(t,y,forcing,water_only)

        i=min(int(max(t,0)//3600),len(forcing)-1)
        Ta,RHa,rain,sw,lw,wind,supply,csource=forcing[i]
        p,rh,kl,kv,lam=self.m.fields(liquid,T)
        pv=rh*psat(T); pvair=RHa/100.*psat(Ta)
        ql=np.zeros(self.n+1); qv=np.zeros(self.n+1)
        face_kl=self.face_k(kl) if self._liquid_integral is None else self._liquid_integral.face_k(liquid,p,kl)
        ql[1:-1]=-face_kl/RW*((p[1:]-p[:-1])/self.dist+RW*self.gravity)
        qv[1:-1]=-self.face_k(kv)*(pv[1:]-pv[:-1])/self.dist/RW
        # Pore occupancy uses total condensable water; mobile transport uses liquid.
        limit=np.clip(1.-(np.maximum(total,0)/self.m.hi)**4,0,1)
        ql[0]=supply/1000./86400.*limit[0]
        exposure=np.full(self.n,self.a); exposure[-1]+=1./self.widths[-1]
        hc=5.7+3.8*wind
        beta=hc/(1.2*1005.*RV*(Ta+273.15))
        conduct=kv+kl*RW*RV*(T+273.15)/np.maximum(pv,1.)
        beta_eff=1./(1./beta+(self.b/4.)/np.maximum(conduct,1e-30))
        evapor=exposure*beta_eff*(pv-pvair)
        rainfall=np.full(self.n,self.a*self.surface['rain_catch_fraction']*rain/1000./3600.)
        if self.surface['crown_rain']:
            rainfall[-1]+=rain/1000./3600./self.widths[-1]
        captured=rainfall*limit
        dry=(total<=self.m.lo) & (ice<=ICE_EPS)
        ql[1:-1]=np.where(((ql[1:-1]>0)&dry[:-1])|((ql[1:-1]<0)&dry[1:]),0.,ql[1:-1])
        qv[1:-1]=np.where(((qv[1:-1]>0)&dry[:-1])|((qv[1:-1]<0)&dry[1:]),0.,qv[1:-1])
        if self.closed:
            ql[0]=0; evapor[:]=0; captured[:]=0; exposure[:]=0
        dtotal=(ql[:-1]-ql[1:]+qv[:-1]-qv[1:])/self.widths + captured-evapor/RW
        self._raw_dtheta=dtotal.copy()
        wet=getattr(self,'_active_wet',None)
        if wet is None: wet=(total>=self.m.hi)
        seep=np.where(wet & (exposure>0) & (not self.closed),np.maximum(dtotal,0.),0.)
        dtotal-=seep
        prevented=np.where(dry,np.minimum(np.maximum(-dtotal,0.),np.maximum(evapor/RW,0.)),0.)
        evapor-=RW*prevented; dtotal+=prevented

        qheat=np.zeros(self.n+1)
        qheat[1:-1]=-self.face_k(lam)*(T[1:]-T[:-1])/self.dist
        sky=self.surface['sky_view_factor']
        effective_lw=sky*lw+(1-sky)*SIGMA*(Ta+273.15)**4
        absorbed=self.surface['absorptivity']*self.surface['solar_projection_factor']*sw
        net=hc*(Ta-T)+absorbed+self.surface['emissivity']*(effective_lw-SIGMA*(T+273.15)**4)
        energy=(qheat[:-1]-qheat[1:])/self.widths + exposure*net-LV*evapor
        active=ice>ICE_EPS
        dcap=self.capacity_dT(T)
        sensible=self.m.rho*self.m.cp + RW*(CW*liquid + CI*ice)
        # H = sensible reference enthalpy - rho_w Lf theta_ice.
        # In the equilibrium branch theta_ice = theta_total - theta_u(T), giving
        # C_app*dT/dt = Q + rho_w Lf*dtheta_total/dt.
        capp=sensible + RW*LF*dcap*active
        dT=(energy + RW*LF*dtotal*active)/np.maximum(capp,1e-12)

        c,solid=ColumnModel.partition(self,liquid,T,ns)
        salt_flux=np.zeros(self.n+1)
        face_theta=np.maximum((liquid[1:]+liquid[:-1])/2.,0)
        sat_fraction=np.clip(face_theta/self.m.sat,0,1)
        D=float(self.salt['D_sulfate_m2_s'])*sat_fraction**float(self.salt['unsaturated_exponent'])
        if self.salt['diffusion_convention']=='pore': D=D*face_theta
        salt_flux[1:-1]=np.where(ql[1:-1]>=0,c[:-1],c[1:])*ql[1:-1]-D*np.diff(c)/self.dist
        salt_flux[0]=ql[0]*(csource/96.06)
        if water_only: salt_flux[:]=0
        salt_loss=np.zeros(self.n) if water_only else seep*c
        dns=(salt_flux[:-1]-salt_flux[1:])/self.widths-salt_loss
        local=np.zeros((self.n,self.nl))
        local[:,2]=RW*seep; local[:,3]=salt_loss
        local[:,4]=RW*(rainfall-captured) if not self.closed else 0.
        local[:,5]=RW*captured; local[:,6]=evapor
        local[0,7]=RW*ql[0]/self.widths[0]; local[0,8]=salt_flux[0]/self.widths[0]
        local[:,0]=local[:,7]+local[:,5]-local[:,6]-local[:,2]
        local[:,1]=local[:,8]-local[:,3]
        self._ledger_local=local
        return np.r_[np.column_stack((dtotal,dT,dns)).ravel(),np.sum(local*self.widths[:,None],axis=0)]

    def solve_hour(self,y,forcing_one,water_only,maxstep,rtol,rhs):
        """Original wet-end active set without the old no-ice temperature stop."""
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
            out=solve_ivp(rhs,(t0,3600.),z,method='BDF',first_step=min(30.,maxstep,3600.-t0),
                max_step=maxstep,rtol=rtol,atol=self.atol,
                jac=lambda t,x:self.jacobian(t,x,forcing_one,water_only),events=[hit,release])
            if not out.success:return out
            self._min_internal_T=min(self._min_internal_T,float(out.y[1:self.size:3].min()))
            self._set_hour_diag_from_solution(out)
            steps+=len(out.t)-1; z=out.y[:,-1]; t0=float(out.t[-1])
            if out.status!=1:break
            events_count+=1
            if events_count>5*self.n+100:
                raise InputError('Wet-end active-set event limit reached. No unconverged hour accepted.')
            if len(out.t_events[0]):
                th=z[:self.size].reshape(self.n,3)[:,0]; candidates=np.flatnonzero(inactive)
                which=candidates[np.argmin(self.m.hi-th[candidates])]; self._active_wet[which]=True
            else:
                self.derivatives(t0,z,forcing_one,water_only); candidates=np.flatnonzero(active)
                which=candidates[np.argmin(self._raw_dtheta[candidates])]; self._active_wet[which]=False
        self._active_wet=None; self._events_total+=events_count; self._steps_total+=steps
        return SimpleNamespace(success=True,y=z[:,None],message='hour complete',event_switches=events_count,internal_steps=steps)

    def integrate(self,climate,initial=None,cancel=None,progress=None,chunk_callback=None,water_only=False):
        """Checkpoint-compatible integration with explicit phase diagnostics."""
        if initial is None:initial=self.initial()
        state=initial.matrix()
        if state.shape!=(self.n,3):raise InputError('Initial state and current mesh differ; regenerate the spin-up.')
        if (not np.isfinite(state).all() or state[:,0].min()<self.m.lo-1e-12 or state[:,0].max()>self.m.hi+1e-12 or state[:,2].min()<-1e-12):
            raise InputError('Initial/checkpoint state is outside the supported total-water or nonnegative salt range. No state was clipped.')
        y=np.r_[state.ravel(),np.zeros(self.nl)]
        block=int(self.cfg['numerics']['chunk_hours']);rows=[];profiles=[];rejected=0
        maxstep=float(self.cfg['numerics']['max_step_s']);rtol=float(self.cfg['numerics']['rtol'])
        for start in range(0,climate.n,block):
            if cancel and cancel():raise Cancelled('Run cancelled after the last saved chunk. Use resume to restart at the last checkpoint.')
            end=min(start+block,climate.n);frame=climate.data.iloc[start:end];forc_all=self.forcing_array(frame)
            y[self.size:]=0.;current=y[:self.size].reshape(self.n,3)
            init_water=RW*np.sum(current[:,0]*self.widths);init_salt=np.sum(current[:,2]*self.widths)
            raw_hours=[];ext_water_hours=[];ext_salt_hours=[];ledger_hours=[];phase_hours=[]
            for local in range(len(frame)):
                if cancel and cancel():raise Cancelled('Run cancelled; completed checkpoint chunks are retained.')
                forcing_one=forc_all[local:local+1];absolute_row=start+local;failure_notes=[];result=None
                def rhs(t,z):
                    if cancel and cancel():raise Cancelled('Run cancelled; completed checkpoint chunks are retained.')
                    return self.derivatives(t,z,forcing_one,water_only)
                for factor in (1.,.25,.0625,.015625):
                    self._hour_diag={'max_ice_theta':0.0,'max_ice_fraction':0.0,'min_liquid_theta':float('inf'),'phase_internal_samples':0}
                    try:trial=self.solve_hour(y,forcing_one,water_only,maxstep*factor,min(rtol,3e-8),rhs)
                    except (Cancelled,PhysicsScopeError):raise
                    except Exception as exc:
                        rejected+=1;failure_notes.append(f'{factor:g}x step exception: {type(exc).__name__}: {exc}');continue
                    if not trial.success:
                        rejected+=1;failure_notes.append(f'{factor:g}x step solver: {trial.message}');continue
                    raw_one=trial.y[:self.size,-1].reshape(self.n,3);tot=raw_one[:,0];ns=raw_one[:,2]
                    good=(np.isfinite(trial.y).all() and tot.min()>=self.m.lo-1e-9 and tot.max()<=self.m.hi+1e-9 and ns.min()>=-1e-10)
                    if not good:
                        rejected+=1;failure_notes.append(f'{factor:g}x step rejected state: total theta [{tot.min():.6g},{tot.max():.6g}] allowed [{self.m.lo:.6g},{self.m.hi:.6g}], sulfate min {ns.min():.6g}');continue
                    result=trial;break
                if result is None:
                    detail=' | '.join(failure_notes[-4:]) if failure_notes else 'no solver diagnostic returned'
                    raise InputError(f'Phase-enabled numerical integration failed at climate row {absolute_row} ({frame.timestamp_utc.iloc[local]}), mesh={self.n} cells. Solver diagnostics: {detail}')
                y=result.y[:,-1];X=y[:self.size].reshape(self.n,3);tot,T,ns=X.T;liq,ice,_=self.phase_partition(tot,T)
                has_ice=bool(np.any(ice>ICE_EPS));
                if has_ice:self._phase_hours+=1
                if has_ice!=self._prev_hour_has_ice:self._freeze_thaw_switches+=1
                self._prev_hour_has_ice=has_ice
                self._phase_max_ice=max(self._phase_max_ice,float(ice.max()))
                frac=np.divide(ice,np.maximum(tot,1e-30));self._phase_max_ice_fraction=max(self._phase_max_ice_fraction,float(frac.max()))
                diag=dict(self._hour_diag);diag.update(phase_active_cells=int(np.sum(ice>ICE_EPS)),ice_theta_mean_m3_m3=float(np.sum(ice*self.widths)/self.H),ice_theta_max_m3_m3=float(ice.max()),liquid_theta_mean_m3_m3=float(np.sum(liq*self.widths)/self.H),liquid_theta_min_m3_m3=float(liq.min()),total_theta_mean_m3_m3=float(np.sum(tot*self.widths)/self.H),phase_model=VERSION)
                phase_hours.append(diag);raw_hours.append(X.copy());ledger_hours.append(y[self.size:].copy());ext_water_hours.append(float(y[self.size]));ext_salt_hours.append(float(y[self.size+1]))
                if progress:
                    done=absolute_row+1;progress(done/climate.n,f'{done:,} / {climate.n:,} hours; mesh {self.n} cells; ice max {float(ice.max()):.3g}')
            raw=np.asarray(raw_hours);ext_water_series=np.asarray(ext_water_hours);ext_salt_series=np.asarray(ext_salt_hours)
            water=RW*(raw[:,:,0]*self.widths).sum(axis=1);salt=(raw[:,:,2]*self.widths).sum(axis=1)
            ew=water-init_water-ext_water_series;es=salt-init_salt-ext_salt_series
            if np.max(np.abs(ew))>max(1e-4,1e-5*max(init_water,1)) or np.max(np.abs(es))>max(1e-7,1e-5*max(init_salt,1)):
                raise InputError('Mass-closure test failed in phase-enabled segment. Outputs are not accepted as a completed run.')
            increments=np.diff(np.vstack([np.zeros((1,self.nl)),np.asarray(ledger_hours)]),axis=0)
            for k in range(len(frame)):
                total,T,ns=raw[k].T;liquid,ice,_=self.phase_partition(total,T);c,solid=ColumnModel.partition(self,liquid,T,ns)
                band=.10;edges=np.arange(0.,self.H+band,band)
                if edges[-1]<self.H:edges=np.r_[edges,self.H]
                else:edges[-1]=self.H
                cell_lo=self.edges[:-1];cell_hi=self.edges[1:];band_c=[];band_ice=[]
                for lo,hi in zip(edges[:-1],edges[1:]):
                    overlap=np.maximum(0.,np.minimum(cell_hi,hi)-np.maximum(cell_lo,lo));lv=np.sum(liquid*overlap);span=max(hi-lo,1e-15)
                    if lv>1e-15:band_c.append(float(np.sum(c*liquid*overlap)/lv))
                    band_ice.append(float(np.sum(ice*overlap)/span))
                c_band_max=max(band_c) if band_c else float(c.max()); ice_band_max=max(band_ice) if band_ice else float(ice.max())
                item={'forcing_timestamp_utc':str(frame.timestamp_utc.iloc[k]),'elapsed_hours':start+k+1,'closure_segment_start_hour':start,
                  'theta_min_m3_m3':float(liquid.min()),'theta_mean_m3_m3':float(np.sum(liquid*self.widths)/self.H),'theta_max_m3_m3':float(liquid.max()),
                  'total_theta_min_m3_m3':float(total.min()),'total_theta_mean_m3_m3':float(np.sum(total*self.widths)/self.H),'total_theta_max_m3_m3':float(total.max()),
                  'ice_theta_mean_m3_m3':float(np.sum(ice*self.widths)/self.H),'ice_theta_max_m3_m3':float(ice.max()),'ice_theta_max_0p10m_bandavg_m3_m3':ice_band_max,'phase_active_cells':int(np.sum(ice>ICE_EPS)),
                  'T_min_C':float(T.min()),'T_mean_C':float(np.sum(T*self.widths)/self.H),'T_max_C':float(T.max()),
                  'water_inventory_kg_m2_footprint':float(water[k]),'ice_water_inventory_kg_m2_footprint':float(RW*np.sum(ice*self.widths)),
                  'sulfate_total_mol_m2_footprint':float(salt[k]),'sulfate_cmax_mol_m3_liquid':float(c.max()),'sulfate_cmax_0p10m_bandavg_mol_m3_liquid':c_band_max,
                  'water_closure_error_kg_m2':float(ew[k]),'sulfate_closure_error_mol_m2':float(es[k]),'phase_model':VERSION}
                item.update({k2:v for k2,v in phase_hours[k].items() if k2 not in item})
                if self.solubility is not None:item['equilibrium_solid_mol_m2_footprint']=float(np.sum(solid*self.widths))
                item.update({name+'_hour':float(increments[k,j]) for j,name in enumerate(LEDGERS)});rows.append(item)
            total,T,ns=raw[-1].T;liquid,ice,_=self.phase_partition(total,T)
            profile=pd.DataFrame({'elapsed_hours':end,'height_m':self.z,'theta_m3_m3':liquid,'total_theta_m3_m3':total,'ice_theta_m3_m3':ice,'temperature_C':T,'sulfate_total_mol_m3_bulk':ns})
            profiles.append(profile)
            if chunk_callback:chunk_callback(end,State.from_matrix(raw[-1]),pd.DataFrame(rows[-len(frame):]),profile)
            if progress:progress(end/climate.n,f'{end:,} / {climate.n:,} hours; water closure {np.max(np.abs(ew)):.2e} kg/m2')
        report={'solver':'1-D finite-volume / SciPy BDF + equilibrium pore-water phase partition','phase_model':VERSION,'cells':self.n,'liquid_face_scheme':self.liquid_face_scheme,'mesh_scheme':self.mesh_scheme,'smallest_cell_m':float(self.widths.min()),'largest_cell_m':float(self.widths.max()),'max_step_s':maxstep,'requested_rtol':rtol,'effective_rtol':min(rtol,3e-8),'completed_hours':climate.n,'rejected_hour_attempts':rejected,
          'water_max_closure_error_kg_m2':max(abs(r['water_closure_error_kg_m2']) for r in rows),'sulfate_max_closure_error_mol_m2':max(abs(r['sulfate_closure_error_mol_m2']) for r in rows),
          'phase_active_hours':self._phase_hours,'freeze_thaw_state_switches':self._freeze_thaw_switches,'maximum_ice_water_equivalent_m3_m3':self._phase_max_ice,'maximum_ice_fraction_of_total_water':self._phase_max_ice_fraction,'min_internal_material_temperature_C':self._min_internal_T,
          'phase_assumptions':'Local thermodynamic equilibrium; pure-water freezing relation from retention curve; no salt freezing-point depression, nucleation hysteresis, ice-pressure damage, pore-expansion mechanics, or explicit snow layer.',
          'salt_phase_behavior':'Transported Na2SO4-equivalent remains in liquid phase; freezing concentrates dissolved-equivalent salt. No mixed-electrolyte activity model.',
          'latent_heat_J_kg':LF,'ice_heat_capacity_J_kgK':CI,'phase_temperature_floor_C':MIN_PHASE_C,
          'material_name':self.m.name,'is_demo':self.m.is_demo,'boundary_closure':'source-supported exposed-face seepage','supported_wet_endpoint':self.m.hi,'effective_saturation':self.m.sat,'coverage_gap_fraction':self._capacity_gap,'wet_end_event_switches':self._events_total,'accepted_internal_steps':self._steps_total,
          'seepage_water_kg_m2':sum(r['seepage_water_kg_m2_hour'] for r in rows),'seepage_salt_mol_m2':sum(r['seepage_salt_mol_m2_hour'] for r in rows),'unabsorbed_rain_kg_m2':sum(r['unabsorbed_rain_kg_m2_hour'] for r in rows)}
        return State.from_matrix(y[:self.size].reshape(self.n,3)),pd.DataFrame(rows),pd.concat(profiles,ignore_index=True),report
