"""Hour/chunk integration adapter around the unchanged transport equations.

ColdAirColumnModel retains the original material-temperature-only stop.
PoreLiquidModel overrides solve_hour and attaches the separately reviewed
pore-liquid acceptance screen and per-hour diagnostics. Original scientific
source files are never edited by either adapter.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from hmsl.solver import State, LEDGERS, PhysicsScopeError
from hmsl.common import InputError, Cancelled
from hmsl.material import RW
from hmsl.solver import ColumnModel
from .exact_cache import CachedModel

class ColdAirScopeMixin:
    def integrate(self, climate, initial=None, cancel=None, progress=None, chunk_callback=None, water_only=False):
        """Retained equations and chunking; cold-air and optional scope hooks."""
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
            y[self.size:]=0.
            current=y[:self.size].reshape(self.n,3)
            init_water=RW*np.sum(current[:,0]*self.widths)
            init_salt=np.sum(current[:,2]*self.widths)
            raw_hours=[]; ext_water_hours=[]; ext_salt_hours=[]; ledger_hours=[]; cold_hours=[]
            for local in range(len(frame)):
                if cancel and cancel(): raise Cancelled('Run cancelled; completed checkpoint chunks are retained.')
                forcing_one=forc_all[local:local+1]
                absolute_row=start+local
                failure_notes=[]; result=None
                def rhs(t,z):
                    if cancel and cancel(): raise Cancelled('Run cancelled; completed checkpoint chunks are retained.')
                    return self.derivatives(t,z,forcing_one,water_only)
                for factor in (1., .25, .0625, .015625):
                    if hasattr(self,'_begin_cold_hour'): self._begin_cold_hour(frame.timestamp_utc.iloc[local],forcing_one[0])
                    try:
                        trial=self.solve_hour(y,forcing_one,water_only,maxstep*factor,min(rtol,3e-8),rhs)
                    except (Cancelled,PhysicsScopeError):
                        raise
                    except Exception as exc:
                        rejected+=1; failure_notes.append(f'{factor:g}x step exception: {type(exc).__name__}: {exc}'); continue
                    if not trial.success:
                        rejected+=1; failure_notes.append(f'{factor:g}x step solver: {trial.message}'); continue
                    raw_one=trial.y[:self.size,-1].reshape(self.n,3)
                    th=raw_one[:,0]; ns=raw_one[:,2]
                    good=(np.isfinite(trial.y).all() and th.min()>=self.m.lo-1e-9 and
                          th.max()<=self.m.hi+1e-9 and ns.min()>=-1e-10)
                    if not good:
                        rejected+=1
                        failure_notes.append(f'{factor:g}x step rejected state: theta [{th.min():.6g},{th.max():.6g}] allowed [{self.m.lo:.6g},{self.m.hi:.6g}], sulfate min {ns.min():.6g}')
                        continue
                    result=trial; break
                if result is None:
                    forcing_desc=(f'Tair={forcing_one[0,0]:.3f} C, RH={forcing_one[0,1]:.3f} %, rain={forcing_one[0,2]:.6g} mm/h, SW={forcing_one[0,3]:.3f} W/m2, LW={forcing_one[0,4]:.3f} W/m2, wind={forcing_one[0,5]:.3f} m/s, supply={forcing_one[0,6]:.6g} mm/day, SO4={forcing_one[0,7]:.6g} mg/L')
                    detail=' | '.join(failure_notes[-4:]) if failure_notes else 'no solver diagnostic returned'
                    raise InputError(f'Numerical integration failed at climate row {absolute_row} ({frame.timestamp_utc.iloc[local]}), mesh={self.n} cells. {forcing_desc}. Material theta range [{self.m.lo:.6g},{self.m.hi:.6g}]. Solver diagnostics: {detail}')
                y=result.y[:,-1]
                cold_hours.append(self._cold_hour_values() if hasattr(self,'_cold_hour_values') else {})
                raw_hours.append(y[:self.size].reshape(self.n,3).copy())
                ledger_hours.append(y[self.size:].copy())
                ext_water_hours.append(float(y[self.size])); ext_salt_hours.append(float(y[self.size+1]))
                if progress:
                    done=absolute_row+1; progress(done/climate.n,f'{done:,} / {climate.n:,} hours; mesh {self.n} cells')
            raw=np.asarray(raw_hours); ext_water_series=np.asarray(ext_water_hours); ext_salt_series=np.asarray(ext_salt_hours)
            water=RW*(raw[:,:,0]*self.widths).sum(axis=1); salt=(raw[:,:,2]*self.widths).sum(axis=1)
            ew=water-init_water-ext_water_series; es=salt-init_salt-ext_salt_series
            if np.max(np.abs(ew))>max(1e-4,1e-5*max(init_water,1)) or np.max(np.abs(es))>max(1e-7,1e-5*max(init_salt,1)):
                raise InputError('Mass-closure test failed. Outputs are not accepted as a completed run.')
            increments=np.diff(np.vstack([np.zeros((1,self.nl)),np.asarray(ledger_hours)]),axis=0)
            for k in range(len(frame)):
                th,T,ns=raw[k].T; c,solid=self.partition(th,T,ns)
                band=0.10; edges=np.arange(0.,self.H+band,band)
                if edges[-1]<self.H: edges=np.r_[edges,self.H]
                else: edges[-1]=self.H
                cell_lo=self.edges[:-1]; cell_hi=self.edges[1:]; band_c=[]
                for lo,hi in zip(edges[:-1],edges[1:]):
                    overlap=np.maximum(0.,np.minimum(cell_hi,hi)-np.maximum(cell_lo,lo)); liquid=np.sum(th*overlap)
                    if liquid>1e-15: band_c.append(float(np.sum(c*th*overlap)/liquid))
                c_band_max=max(band_c) if band_c else float(c.max())
                item={'forcing_timestamp_utc':str(frame.timestamp_utc.iloc[k]),'elapsed_hours':start+k+1,'closure_segment_start_hour':start,
                  'theta_min_m3_m3':float(th.min()),'T_min_C':float(T.min()),'theta_mean_m3_m3':float(np.sum(th*self.widths)/self.H),'theta_max_m3_m3':float(th.max()),
                  'T_mean_C':float(np.sum(T*self.widths)/self.H),'T_max_C':float(T.max()),'water_inventory_kg_m2_footprint':float(water[k]),'sulfate_total_mol_m2_footprint':float(salt[k]),
                  'sulfate_cmax_mol_m3_liquid':float(c.max()),'sulfate_cmax_0p10m_bandavg_mol_m3_liquid':c_band_max,
                  'water_closure_error_kg_m2':float(ew[k]),'sulfate_closure_error_mol_m2':float(es[k])}
                item.update(cold_hours[k])
                if self.solubility is not None: item['equilibrium_solid_mol_m2_footprint']=float(np.sum(solid*self.widths))
                item.update({name+'_hour':float(increments[k,j]) for j,name in enumerate(LEDGERS)})
                rows.append(item)
            profile=pd.DataFrame({'elapsed_hours':end,'height_m':self.z,'theta_m3_m3':raw[-1,:,0],'temperature_C':raw[-1,:,1],'sulfate_total_mol_m3_bulk':raw[-1,:,2]})
            profiles.append(profile)
            if chunk_callback: chunk_callback(end,State.from_matrix(raw[-1]),pd.DataFrame(rows[-len(frame):]),profile)
            if progress: progress(end/climate.n,f'{end:,} / {climate.n:,} hours; water closure {np.max(np.abs(ew)):.2e} kg/m2')
        report={'solver':'1-D finite-volume / SciPy BDF; hourly constant-forcing IVPs','cells':self.n,'liquid_face_scheme':self.liquid_face_scheme,'mesh_scheme':self.mesh_scheme,'smallest_cell_m':float(self.widths.min()),'largest_cell_m':float(self.widths.max()),'max_step_s':maxstep,'requested_rtol':rtol,'effective_rtol':min(rtol,3e-8),
          'completed_hours':climate.n,'rejected_hour_attempts':rejected,'water_max_closure_error_kg_m2':max(abs(r['water_closure_error_kg_m2']) for r in rows),'sulfate_max_closure_error_mol_m2':max(abs(r['sulfate_closure_error_mol_m2']) for r in rows),
          'salt_phase_mode':'optional single-salt equilibrium partition' if self.solubility else 'transport only; all salt treated as dissolved-equivalent; no phase/damage output','material_name':self.m.name,'is_demo':self.m.is_demo,
          'boundary_closure':'source-supported exposed-face seepage','supported_wet_endpoint':self.m.hi,'effective_saturation':self.m.sat,'coverage_gap_fraction':self._capacity_gap,'wet_end_event_switches':self._events_total,
          'accepted_internal_steps':self._steps_total,'min_internal_material_temperature_C':self._min_internal_T,'seepage_water_kg_m2':sum(r['seepage_water_kg_m2_hour'] for r in rows),
          'seepage_salt_mol_m2':sum(r['seepage_salt_mol_m2_hour'] for r in rows),'unabsorbed_rain_kg_m2':sum(r['unabsorbed_rain_kg_m2_hour'] for r in rows),
          'cold_air_scope':('Pore-liquid stability screen active; no ice dynamics or frost damage prediction.' if hasattr(self,'_cold_hour_values') else 'Sub-zero outdoor-air forcing permitted; original simulated-material <0 C runtime guard retained.')}
        return State.from_matrix(y[:self.size].reshape(self.n,3)),pd.DataFrame(rows),pd.concat(profiles,ignore_index=True),report

class ColdAirColumnModel(ColdAirScopeMixin, ColumnModel):
    pass

class ColdAirCachedModel(ColdAirScopeMixin, CachedModel):
    pass
