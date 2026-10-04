"""Acceptance-only adapter. State equations are byte-unchanged in hmsl.

Original solve_hour wet-end logic is retained verbatim except for dense output
and the temperature-only scope test. Every accepted cold BDF state is screened;
when any accepted temperature is within 0.5 C of zero, dense interpolation is
also checked at intervals of at most 60 s. This is numerical screening, not proof
of an exact continuous minimum or physical validation of the brick's freezing.
"""
from __future__ import annotations
import json
import numpy as np
from scipy.integrate import solve_ivp
from hmsl.solver import ColumnModel, PhysicsScopeError
from hmsl.common import InputError
from .exact_cache import CachedModel
from .scope_model import ColdAirScopeMixin
from .pore_ice import assess, METHOD

class PoreScopeMixin:
    def _begin_cold_hour(self, stamp, forcing):
        self._cold_timestamp=str(stamp)
        self._cold_hour={'cold_material_seen':False,'min_material_T_internal_C':float('inf'),
            'max_pore_ice_saturation':0.0,'max_ambient_frost_saturation':0.0,'sampled_states':0}
        # Classify cold precipitation even if all modeled brick cells are warm.
        if forcing[0]<=0.0 and forcing[2]>0.0:
            raise PhysicsScopeError('COLD_PRECIPITATION_PHASE_UNRESOLVED at '+str(stamp)+
                ': the supplied precipitation has no rain/snow phase label; liquid rain cannot be assumed here. No climate values changed.')

    def _check_cold_solution(self, out, forcing):
        samples=out.y[:self.size]; times=out.t.copy()
        if samples[1::3].min()<=0.5 and out.sol is not None:
            tt=np.unique(np.concatenate([np.linspace(a,b,max(2,int(np.ceil((b-a)/60.))+1))
                for a,b in zip(out.t[:-1],out.t[1:])])) if len(out.t)>1 else out.t
            # Include every accepted solver point plus supplementary interpolation.
            samples=np.column_stack([samples,out.sol(tt)[:self.size]]); times=np.r_[times,tt]
            order=np.argsort(times);samples=samples[:,order];times=times[order]
        X=samples.reshape(self.n,3,-1)
        if not np.isfinite(X).all():raise InputError('Non-finite accepted cold-screen states.')
        h=self._cold_hour
        h['min_material_T_internal_C']=min(h['min_material_T_internal_C'],float(X[:,1,:].min()))
        h['sampled_states']+=int(X.shape[2])
        which=np.flatnonzero(np.any(X[:,1,:]<0,axis=0))
        exposed=np.full(self.n,self.a>0) if not self.closed else np.zeros(self.n,bool)
        if not self.closed:exposed[-1]=True
        for k in which:
            record=assess(self.m,X[:,0,k],X[:,1,k],forcing,exposed)
            h['cold_material_seen']=True
            h['max_pore_ice_saturation']=max(h['max_pore_ice_saturation'],record['maximum_ice_saturation'])
            h['max_ambient_frost_saturation']=max(h['max_ambient_frost_saturation'],record['maximum_ambient_frost_saturation'])
            if not record['allowed']:
                record['forcing_timestamp_utc']=self._cold_timestamp
                record['sample_offset_seconds_in_hour']=float(times[k])
                if record['cell'] is not None:record['height_m']=float(self.z[record['cell']])
                self.last_phase_event=record
                raise PhysicsScopeError('PHASE_REVIEW_REQUIRED: '+json.dumps(record,sort_keys=True)+
                    '. Potential ice/frost is outside this liquid-only model. Not an observed freezing event; no failing hour committed.')

    def _cold_hour_values(self):
        h=self._cold_hour
        return {'cold_scope_min_internal_T_C':h['min_material_T_internal_C'],
                'cold_scope_subzero_material':int(h['cold_material_seen']),
                'cold_scope_max_pore_ice_saturation':h['max_pore_ice_saturation'],
                'cold_scope_max_ambient_frost_saturation':h['max_ambient_frost_saturation'],
                'cold_scope_sampled_state_vectors':h['sampled_states']}

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
                events=[hit,release],dense_output=True)
            if not out.success: return out
            self._min_internal_T=min(self._min_internal_T,float(out.y[1:self.size:3].min()))
            self._check_cold_solution(out, forcing_one[0])
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

class PoreLiquidModel(PoreScopeMixin,ColdAirScopeMixin,ColumnModel):
    pass

class PoreLiquidCachedModel(PoreScopeMixin,ColdAirScopeMixin,CachedModel):
    pass
