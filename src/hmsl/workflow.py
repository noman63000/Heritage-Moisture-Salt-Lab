"""Auditable work flows used by both the GUI and the command line."""
from __future__ import annotations
import copy, json, os, platform, sys, threading, time
from pathlib import Path
from importlib.metadata import version
import numpy as np
import pandas as pd
from . import __version__
from .common import InputError, Cancelled, atomic_json, digest, load_json, now, sha256, safe_name
from .project import Project
from .climate import climate_screen
from .solver import ColumnModel, State
from .report import write_report



def qa_case_path(project,rid):
    d=project.path/'qa_cases'; d.mkdir(exist_ok=True)
    return d/(safe_name(rid)+'.json')

def save_qa_record(project,rid,record):
    """Keep the latest QA plus a per-case audit copy.

    Per-case copies let later failures on another climate case coexist with a
    previously passing case and make overnight queues resumable/auditable.
    """
    atomic_json(project.path/'qa.json',record)
    atomic_json(qa_case_path(project,rid),record)

def load_case_qa(project,rid):
    p=qa_case_path(project,rid)
    return load_json(p) if p.is_file() else {}

def valid_case_qa(project,rid):
    q=load_case_qa(project,rid)
    if not q.get('passed') or q.get('physics_hash')!=project.physics_hash() or q.get('run_id')!=rid:
        return False
    try:
        run=project.run(rid)
        return q.get('climate_sha256')==sha256(project.file(run['file']))
    except Exception:
        return False

LIMITATIONS=[
 'This is a research implementation, not a calibrated site model or validated DELPHIN substitute.',
 'One vertical dimension, uniform across thickness; lateral exchange is a distributed source. No 2-D wall/soil solution.',
 'Lower water input is a reviewed parametric flux, not a flux derived from the measured groundwater depth.',
 'One transported Na2SO4-equivalent binary salt, not separately charged ions or a NaCl/Na2SO4 mixed electrolyte.',
 'Unsaturated sulfate mobility, diffusion convention, rain capture and radiation projection are explicit modelling assumptions.',
 'No ice, salt-dependent sorption, pore clogging, crystallisation kinetics, hydration-water balance or mechanical damage.',
 'Wet-end exposed-face seepage is a reduced boundary assumption, not a measured drainage law or saturated pressure solution. Source table endpoints are not extrapolated. The small coverage-to-saturation gap is reported.',
 'The heat balance omits water-carried sensible enthalpy and phase-change enthalpy beyond evaporation/condensation.',
 'Optional solubility data provide a single-salt equilibrium partition only; phase cycling and supersaturation require a fuller thermodynamic treatment.',
 'Climate screening cutoffs are not calibrated deterioration thresholds.']

def provenance(project,climate,rid,mode):
    return {'application':'Heritage Moisture & Salt Lab','application_version':__version__,
      'created_utc':now(),'site':project.config['site_name'],'run_id':rid,'mode':mode,
      'is_demo':project.config.get('is_demo',False),'climate_sha256':climate.source_sha256,
      'climate_file':climate.source_name,'calendar':climate.calendar,'physics_hash':project.physics_hash(),
      'config':project.config,'limitations':LIMITATIONS,
      'python':sys.version,'platform':platform.platform(),
      'dependencies':{k:version(k) for k in ('numpy','scipy','pandas','openpyxl')}}

def require_trial(project):
    r=project.readiness()
    if not r['trial_ready']: raise InputError('\n'.join(r['errors']))

def save_state(path,state):
    p=Path(path); temp=p.with_name(p.name+'.tmp.npz')
    np.savez_compressed(temp,theta=state.theta,temperature=state.temperature,sulfate=state.sulfate)
    os.replace(temp,p)

def load_state(path):
    with np.load(path,allow_pickle=False) as z:
        return State(z['theta'].copy(),z['temperature'].copy(),z['sulfate'].copy())

def _folder(p,mode,rid,key):
    dest=p.path/'outputs'/mode/(rid+'_'+key[:12]); dest.mkdir(parents=True,exist_ok=True); return dest

def screen_run(project,rid,progress=None,cancel=None):
    if cancel and cancel():raise Cancelled('Cancelled.')
    if progress:progress(0,'Reading and validating '+rid)
    c=project.climate(rid); key=digest({'climate':c.source_sha256,'calendar':c.calendar,'mode':'screen-v1'})
    folder=_folder(project,'screening',rid,key)
    s,d=climate_screen(c); s['run_id']=rid
    prov=provenance(project,c,rid,'CLIMATE SCREENING ONLY')
    write_report(folder,s,d,prov)
    atomic_json(folder/'completed.json',{'passed':True,'created_utc':now(),'key':key})
    if progress:progress(1,'Screening report saved')
    return str(folder.relative_to(project.path))

def run_transport(project,rid,max_hours=None,require_frozen=False,resume=True,progress=None,cancel=None):
    require_trial(project)
    if require_frozen:
        ready=project.readiness()
        if not ready['frozen_ready']:
            raise InputError('Batch transport needs a current converged water/heat spin-up and reproducibility freeze. Complete those steps first.')
        if not valid_case_qa(project,rid):
            raise InputError(f'{rid} does not have a current passing case-specific Numerical QA record. Run/queue QA for this case before production.')
    c=project.climate(rid)
    if max_hours is not None:
        if int(max_hours)<2: raise InputError('A trial needs at least two hours.')
        c=c.subset(0,min(int(max_hours),c.n))
    m=project.material(); phase=project.phase(); model=ColumnModel(project.config,m,phase)
    mode='SYNTHETIC_DEMO' if project.config.get('is_demo') else ('FROZEN_REDUCED_TRANSPORT' if require_frozen else 'EXPLORATORY_REDUCED_TRANSPORT')
    key=digest({'physics':project.physics_hash(),'climate':c.source_sha256,'calendar':c.calendar,'hours':c.n,'mode':mode})
    folder=_folder(project,mode.lower(),rid,key)
    if resume and (folder/'completed.json').is_file():
        if progress:progress(1,'Already complete; unchanged inputs. Existing result retained.')
        return str(folder.relative_to(project.path))
    prov=provenance(project,c,rid,mode)
    atomic_json(folder/'provenance.json',prov); atomic_json(folder/'resolved_project.json',project.config)
    state=None; completed=0
    if require_frozen:
        state=load_state(project.path/'spinup.npz')
    elif project.readiness()['spinup_valid']:
        state=load_state(project.path/'spinup.npz')
    cp=folder/'checkpoint.json'
    if resume and cp.is_file():
        prev=load_json(cp)
        if prev.get('key')!=key: raise InputError('Checkpoint inputs differ from the current inputs. Use a new output folder.')
        statefile=folder/prev.get('state_file','checkpoint.npz')
        if sha256(statefile)!=prev['state_sha256']: raise InputError('Checkpoint state hash mismatch.')
        completed=int(prev['completed_hours']); state=load_state(statefile)
    if not resume and cp.is_file(): raise InputError('A partial run already exists. Resume it or use another project; partial results are not overwritten.')
    (folder/'chunks').mkdir(exist_ok=True)
    def callback(end,newstate,h,pr):
        end+=completed; h=h.copy(); pr=pr.copy()
        h['elapsed_hours']+=completed; pr['elapsed_hours']+=completed
        h['closure_segment_start_hour']+=completed
        h.to_csv(folder/'chunks'/f'hours_{end:09d}.csv',index=False)
        pr.to_csv(folder/'chunks'/f'profile_{end:09d}.csv',index=False)
        statefile=folder/f'state_{end:09d}.npz'
        save_state(statefile,newstate)
        atomic_json(cp,{'key':key,'completed_hours':end,'state_file':statefile.name,'state_sha256':sha256(statefile),'saved_utc':now()})
        # Remove superseded state files only after the new checkpoint is committed.
        for old in folder.glob('state_*.npz'):
            if old!=statefile: old.unlink(missing_ok=True)
    def prog(f,msg):
        if progress:progress((completed+f*(c.n-completed))/c.n,msg)
    if completed<c.n:
        _,_,_,stats=model.integrate(c.subset(completed),state,cancel,prog,callback)
    else: stats={}
    hourfiles=sorted((folder/'chunks').glob('hours_*.csv'))
    # A checkpoint is authoritative. Ignore an uncommitted chunk left by a crash.
    checkpoint=load_json(cp)
    hourfiles=[f for f in hourfiles if int(f.stem.split('_')[1])<=checkpoint['completed_hours']]
    h=pd.concat([pd.read_csv(f) for f in hourfiles],ignore_index=True)
    if len(h)!=c.n or not np.array_equal(h.elapsed_hours.to_numpy(),np.arange(1,c.n+1)):
        raise InputError('Checkpoint/chunk continuity failed. The run is not marked complete.')
    h.to_csv(folder/'hourly_results.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    pr=pd.concat([pd.read_csv(f) for f in sorted((folder/'chunks').glob('profile_*.csv')) if int(f.stem.split('_')[1])<=c.n],ignore_index=True)
    pr.to_csv(folder/'profiles.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    h['date']=h.forcing_timestamp_utc.str.slice(0,10)
    cols=[s for s in ('theta_mean_m3_m3','theta_max_m3_m3','T_mean_C','T_max_C','sulfate_cmax_0p10m_bandavg_mol_m3_liquid','sulfate_cmax_mol_m3_liquid') if s in h]
    daily=h.groupby('date',sort=False)[cols].mean().reset_index()
    summary={'mode':mode,'run_id':rid,'completed_hours':c.n,'is_demo':project.config.get('is_demo',False),
      'water_inventory_final_kg_m2_footprint':float(h.water_inventory_kg_m2_footprint.iloc[-1]),
      'sulfate_total_final_mol_m2_footprint':float(h.sulfate_total_mol_m2_footprint.iloc[-1]),
      'theta_mean_time_average_m3_m3':float(h.theta_mean_m3_m3.mean()),'temperature_mean_C':float(h.T_mean_C.mean()),
      'maximum_0p10m_bandavg_sulfate_concentration_mol_m3_liquid':float(h.sulfate_cmax_0p10m_bandavg_mol_m3_liquid.max()),
      'diagnostic_single_cell_max_sulfate_concentration_mol_m3_liquid':float(h.sulfate_cmax_mol_m3_liquid.max()),
      'water_max_closure_error_kg_m2':float(h.water_closure_error_kg_m2.abs().max()),
      'sulfate_max_closure_error_mol_m2':float(h.sulfate_closure_error_mol_m2.abs().max()),
      'salt_phase_mode':'user-sourced single-salt equilibrium partition' if phase else 'transport only; no crystal, phase-cycle or damage prediction',
      'seepage_water_total_kg_m2_footprint':float(h.seepage_water_kg_m2_hour.sum()),
      'seepage_sulfate_total_mol_m2_footprint':float(h.seepage_salt_mol_m2_hour.sum()),
      'unabsorbed_rain_total_kg_m2_footprint':float(h.unabsorbed_rain_kg_m2_hour.sum()),
      'wet_end_supported_theta':m.hi,'effective_saturation_theta':m.sat,
      'wet_end_closure':'Source-supported exposed-face drainage; endpoint gap retained; no post-step water clipping.',
      'mass_closure_scope':'independent inventory balance in each committed chunk; physical states persist unchanged across chunks and resumes',
      'material':m.name,'cells':model.n,'max_step_s':project.config['numerics']['max_step_s'],
      'full_requested_climate_case':max_hours is None,'interpretation':'Scenario output, not measured condition or validated damage.'}
    if phase:summary['equilibrium_solid_final_mol_m2_footprint']=float(h.equilibrium_solid_mol_m2_footprint.iloc[-1])
    prov['solver_stats_last_segment']=stats; prov['warnings']=c.warnings+m.warnings
    write_report(folder,summary,daily,prov)
    atomic_json(folder/'completed.json',{'passed':True,'created_utc':now(),'key':key,'hours':c.n})
    if progress:progress(1,'Run complete; report and data saved')
    return str(folder.relative_to(project.path))

def spinup(project,rid,progress=None,cancel=None):
    """Water/heat spin-up with recoverable, input-hashed chunk checkpoints."""
    require_trial(project); c=project.climate(rid); cfg=copy.deepcopy(project.config)
    yr=c.data.date.str.slice(0,4).iloc[0]
    cycle=c.subset(0,int(c.data.date.str.startswith(yr).sum()))
    if not cfg.get('is_demo'):
        expected=360*24 if c.calendar=='360_day' else (365*24 if c.calendar in ('365_day','noleap') else (366 if __import__('calendar').isleap(int(yr)) else 365)*24)
        if cycle.n!=expected or str(cycle.data.date.iloc[0])!=yr+'-01-01':
            raise InputError('Spin-up requires a complete first model-calendar year starting 1 January.')
    physics=project.physics_hash()
    key=digest({'physics':physics,'climate':c.source_sha256,'calendar':c.calendar,'year':str(yr),'method':'spinup-1.1'})
    cache=project.path/'spinup_cache'/key;cache.mkdir(parents=True,exist_ok=True)
    checkpoint=cache/'checkpoint.json'
    model=ColumnModel(cfg,project.material(),None);state=model.initial()
    records=[];passed=False;cycles=int(cfg['numerics']['spinup_max_cycles']);first=0;completed=0
    prev=state.matrix().copy()
    if checkpoint.is_file():
        meta=load_json(checkpoint)
        if meta.get('key')!=key: raise InputError('Spin-up cache identity mismatch.')
        for prefix in ('latest','year_start'):
            if sha256(cache/meta[prefix+'_file'])!=meta[prefix+'_sha256']:
                raise InputError('Spin-up checkpoint hash mismatch; cached state was not used.')
        state=load_state(cache/meta['latest_file']);prev=load_state(cache/meta['year_start_file']).matrix()
        records=meta['records'];first=int(meta['cycle_index']);completed=int(meta['completed_hours'])
    for i in range(first,cycles):
        if cancel and cancel():raise Cancelled('Spin-up cancelled; completed chunks can be resumed.')
        if i!=first or completed==0:prev=state.matrix().copy()
        year_start=cache/f'cycle_{i+1:03d}_start.npz'
        if completed==0:save_state(year_start,State.from_matrix(prev))
        def cb(fr,msg):
            if progress:progress((i+(completed+fr*(cycle.n-completed))/cycle.n)/cycles,f'Spin-up cycle {i+1}/{cycles}: '+msg)
        def commit(end,newstate,h,profile):
            total=completed+end
            dest=cache/f'cycle_{i+1:03d}_hour_{total:06d}.npz';save_state(dest,newstate)
            atomic_json(checkpoint,{'key':key,'cycle_index':i,'completed_hours':total,'records':records,
                'latest_file':dest.name,'latest_sha256':sha256(dest),
                'year_start_file':year_start.name,'year_start_sha256':sha256(year_start),'saved_utc':now()})
            for old in cache.glob('cycle_*_hour_*.npz'):
                if old!=dest:old.unlink(missing_ok=True)
        if completed<cycle.n:
            state,_,_,_=model.integrate(cycle.subset(completed),state,cancel,cb,commit,water_only=True)
        th=float(np.max(np.abs(state.theta-prev[:,0])));te=float(np.max(np.abs(state.temperature-prev[:,1])))
        records.append({'cycle':i+1,'max_theta_change_m3_m3':th,'max_temperature_change_C':te})
        if th<=cfg['numerics']['spinup_theta_tolerance'] and te<=cfg['numerics']['spinup_temperature_tolerance_C']:
            passed=True;break
        completed=0
        next_start=cache/f'cycle_{i+2:03d}_start.npz';save_state(next_start,state)
        atomic_json(checkpoint,{'key':key,'cycle_index':i+1,'completed_hours':0,'records':records,
            'latest_file':next_start.name,'latest_sha256':sha256(next_start),
            'year_start_file':next_start.name,'year_start_sha256':sha256(next_start),'saved_utc':now()})
    state.sulfate[:]=cfg['salt']['initial_sulfate_mol_m3_bulk']
    save_state(project.path/'spinup.npz',state)
    result={'passed':passed,'physics_hash':physics,'baseline_run':rid,
      'baseline_climate_sha256':c.source_sha256,'initial_state_sha256':sha256(project.path/'spinup.npz'),
      'created_utc':now(),'records':records,'cycle_hours':cycle.n,
      'method':'Repeated first full model year, water/heat only, then explicit initial salt. Same resulting state for all climate cases. No salt steady state claimed.',
      'resume':'Input-hashed water/heat checkpoints are retained after committed chunks.'}
    atomic_json(project.path/'spinup.json',result)
    if not passed:raise InputError('Spin-up did not meet the stated tolerances in the allowed cycles. The convergence log is saved. Do not freeze; review the recorded convergence trend.')
    return result

def numerical_qa(project,rid,hours=72,progress=None,cancel=None):
    """Staged, restartable numerical QA.

    Stage 1 checks mesh convergence first because a failed mesh test makes the
    time-step tests unnecessary. Stage 2 checks time-step convergence only after
    all mesh windows pass. Each completed variant is cached against the exact
    physics/climate/window/cell/step identity, so cancellation or a browser/server
    restart does not throw away completed QA work.
    """
    require_trial(project); c=project.climate(rid); cfg=copy.deepcopy(project.config)
    hours=min(max(24,int(hours)),c.n); n=int(cfg['geometry']['cells']); dt=float(cfg['numerics']['max_step_s'])
    if n*2>400: raise InputError('QA doubles the selected mesh; choose at most 200 base cells.')
    tol=float(cfg['numerics']['qa_relative_tolerance'])
    physics=project.physics_hash()
    # Three forcing windows: beginning, wettest and hottest complete window.
    rain=c.data.pr_mm_h.rolling(hours).sum().to_numpy(); hot=c.data.tas_C.rolling(hours).mean().to_numpy()
    starts=sorted(set([0,max(0,int(np.nanargmax(rain))-hours+1),max(0,int(np.nanargmax(hot))-hours+1)]))
    cache_dir=project.path/'qa_cache'; cache_dir.mkdir(exist_ok=True)
    qa_version='1.2.0-source-integral-staged'

    def metrics(h):
        return np.array([h.theta_mean_m3_m3.mean(),h.water_inventory_kg_m2_footprint.iloc[-1],h.T_mean_C.mean(),
                         h.sulfate_total_mol_m2_footprint.iloc[-1],
                         h.sulfate_cmax_0p10m_bandavg_mol_m3_liquid.max()])
    def diagnostic(h):
        return np.array([h.sulfate_cmax_mol_m3_liquid.max()])
    floors=np.array([.01,1.,1.,1e-5,.001])

    # base is reused in both stages; fine tests mesh, half/quarter test time step.
    spec={'base':(n,dt),'mesh_fine':(n*2,dt),'time_half':(n,dt/2.),'time_quarter':(n,dt/4.)}
    potential=len(starts)*len(spec)
    completed_counter=[0]

    def run_variant(start,label):
        nn,step=spec[label]
        key=digest({'qa_version':qa_version,'physics':physics,'climate':c.source_sha256,
                    'calendar':c.calendar,'start':int(start),'hours':int(hours),
                    'label':label,'cells':int(nn),'max_step_s':float(step)})
        cache=cache_dir/(key+'.json')
        if cache.is_file():
            saved=load_json(cache)
            completed_counter[0]+=1
            if progress: progress(min(.99,completed_counter[0]/potential),
                f'QA cached: {label} ({nn} cells) window {str(c.data.date.iloc[start])}; reusing verified completed task')
            return np.array(saved['metrics'],float),np.array(saved['diagnostic'],float),saved.get('runtime_s',0.)
        if cancel and cancel(): raise Cancelled('QA cancelled. Completed QA variants are cached and will be reused on restart.')
        ordinal=completed_counter[0]+1
        if progress: progress(completed_counter[0]/potential,
            f'QA task {ordinal}/{potential}: starting {label} ({nn} cells, max step {step:g} s) for {str(c.data.date.iloc[start])}. Fine meshes can take much longer.')
        v=copy.deepcopy(cfg);v['geometry']['cells']=nn;v['numerics']['max_step_s']=step
        model=ColumnModel(v,project.material(),project.phase())
        began=time.monotonic(); stop=threading.Event()
        # solve_ivp can spend minutes inside one checkpoint block. Emit a status
        # heartbeat so the dashboard does not look frozen while the CPU is working.
        def heartbeat():
            while not stop.wait(15.):
                if progress:
                    elapsed=time.monotonic()-began
                    progress(completed_counter[0]/potential,
                        f'QA task {ordinal}/{potential}: {label} ({nn} cells), worker still active; {elapsed/60:.1f} min elapsed. Awaiting solver progress; a heartbeat alone does not prove advancement.')
        hb=threading.Thread(target=heartbeat,daemon=True);hb.start()
        def cb(fr,msg):
            if progress:
                progress(min(.99,(completed_counter[0]+fr)/potential),
                    f'QA task {ordinal}/{potential} {label} ({nn} cells): '+msg)
        try:
            _,h,_,stats=model.integrate(c.subset(start,start+hours),cancel=cancel,progress=cb)
        except Exception as exc:
            save_qa_record(project,rid,{'passed':False,'stage':'integration_error',
                'physics_hash':physics,'qa_method_version':qa_version,'run_id':rid,
                'climate_sha256':c.source_sha256,'created_utc':now(),'base_cells':n,'fine_cells':2*n,
                'failed_variant':label,'failed_start_row':int(start),'cells':nn,'error':str(exc),
                'relative_tolerance':tol,'windows':results.copy(),
                'scope':'Incomplete QA. No convergence pass is inferred from a partial or failed integration.'})
            raise
        finally:
            stop.set()
        runtime=time.monotonic()-began
        m=metrics(h); d=diagnostic(h)
        atomic_json(cache,{'qa_version':qa_version,'key':key,'label':label,'cells':nn,'max_step_s':step,
            'start_row':int(start),'first_date':str(c.data.date.iloc[start]),'hours':hours,
            'metrics':m.tolist(),'diagnostic':d.tolist(),'runtime_s':runtime,'solver_stats':stats,
            'completed_utc':now()})
        completed_counter[0]+=1
        if progress: progress(min(.99,completed_counter[0]/potential),
            f'QA completed {label} ({nn} cells) for {str(c.data.date.iloc[start])} in {runtime/60:.1f} min; result cached')
        return m,d,runtime

    # Stage 1: mesh first. Do not spend time on time-step refinements if the
    # current mesh is not yet acceptable.
    results=[]
    for start in starts:
        b,bd,bt=run_variant(start,'base'); f,fd,ft=run_variant(start,'mesh_fine')
        mesh_error=float(np.max(np.abs(b-f)/np.maximum(np.abs(f),floors)))
        diagnostic_mesh=float(np.max(np.abs(bd-fd)/np.maximum(np.abs(fd),.001)))
        results.append({'start_row':int(start),'first_date':str(c.data.date.iloc[start]),'hours':hours,
          'passed':mesh_error<=tol,'stage':'mesh','maximum_normalized_errors':{'base_vs_mesh_fine':mesh_error},
          'metric_order':['theta time mean','final water','temperature time mean','final sulfate','peak 0.10 m band-averaged dissolved-equivalent concentration'],
          'metric_vectors':{'base':b.tolist(),'mesh_fine':f.tolist()},
          'runtime_seconds':{'base':bt,'mesh_fine':ft},
          'normalized_errors_by_metric':{'base_vs_mesh_fine':(np.abs(b-f)/np.maximum(np.abs(f),floors)).tolist()},
          'diagnostic_only':{'note':'Single-cell peak concentration is diagnostic only because its physical support changes with mesh refinement.',
             'metric_order':['single-cell peak dissolved-equivalent concentration'],
             'maximum_normalized_errors':{'base_vs_mesh_fine':diagnostic_mesh},
             'metric_vectors':{'base':bd.tolist(),'mesh_fine':fd.tolist()}}})
    mesh_pass=all(x['passed'] for x in results)
    if not mesh_pass:
        record={'passed':False,'stage':'mesh','base_cells':n,'fine_cells':2*n,'qa_method_version':qa_version,'physics_hash':physics,'run_id':rid,
          'climate_sha256':c.source_sha256,'created_utc':now(),'relative_tolerance':tol,'windows':results,
          'scope':'Staged local QA. Mesh convergence is tested first on water, temperature, sulfate inventory and a fixed-support 0.10 m band-averaged dissolved concentration. Time-step tests are intentionally skipped until every mesh window passes. Single-cell concentration peak is diagnostic only. Not site validation or proof of convergence in all 21 future cases.',
          'state_initialisation':'Same physical initial RH, T and total salt on every mesh; no interpolation of a coarse solution.',
          'cache_note':'Completed variants are cached in qa_cache and reused only when their exact physics/climate/window/cell/step key matches.'}
        save_qa_record(project,rid,record)
        raise InputError('Numerical QA mesh stage did not pass. qa.json is saved. Time-step tests were skipped to avoid wasting computation; completed variants are cached for the next attempt.')

    # Stage 2: only after mesh passes, verify time-step independence.
    for item,start in zip(results,starts):
        b=np.array(item['metric_vectors']['base'],float); bd=np.array(item['diagnostic_only']['metric_vectors']['base'],float)
        h,hd,ht=run_variant(start,'time_half'); q,qd,qt=run_variant(start,'time_quarter')
        e1=float(np.max(np.abs(b-h)/np.maximum(np.abs(h),floors)))
        e2=float(np.max(np.abs(h-q)/np.maximum(np.abs(q),floors)))
        d1=float(np.max(np.abs(bd-hd)/np.maximum(np.abs(hd),.001)))
        d2=float(np.max(np.abs(hd-qd)/np.maximum(np.abs(qd),.001)))
        item['normalized_errors_by_metric'].update({'base_vs_time_half':(np.abs(b-h)/np.maximum(np.abs(h),floors)).tolist(),'time_half_vs_time_quarter':(np.abs(h-q)/np.maximum(np.abs(q),floors)).tolist()})
        item['stage']='mesh+time'; item['maximum_normalized_errors'].update({'base_vs_time_half':e1,'time_half_vs_time_quarter':e2})
        item['metric_vectors'].update({'time_half':h.tolist(),'time_quarter':q.tolist()})
        item['runtime_seconds'].update({'time_half':ht,'time_quarter':qt})
        item['diagnostic_only']['maximum_normalized_errors'].update({'base_vs_time_half':d1,'time_half_vs_time_quarter':d2})
        item['diagnostic_only']['metric_vectors'].update({'time_half':hd.tolist(),'time_quarter':qd.tolist()})
        item['passed']=item['passed'] and e1<=tol and e2<=tol
    record={'passed':all(x['passed'] for x in results),'stage':'complete','base_cells':n,'fine_cells':2*n,'qa_method_version':qa_version,
      'physics_hash':physics,'run_id':rid,'climate_sha256':c.source_sha256,'created_utc':now(),
      'relative_tolerance':tol,'windows':results,
      'scope':'Staged local mesh/timestep QA for selected forcing windows only. Pass/fail uses water, temperature, sulfate inventory and a fixed-support 0.10 m band-averaged dissolved concentration. Single-cell concentration peak is diagnostic only. Not site validation or proof of convergence in all 21 future cases; review extreme future cases separately.',
      'state_initialisation':'Same physical initial RH, T and total salt on every mesh; no interpolation of a coarse solution.',
      'cache_note':'Completed variants are cached in qa_cache and reused only when their exact physics/climate/window/cell/step key matches.'}
    save_qa_record(project,rid,record)
    if not record['passed']:raise InputError('Numerical QA time-step stage did not pass. qa.json is saved; completed variants are cached.')
    if progress:progress(1,'Numerical QA passed. Mesh and time-step tests are current; cached task records retained for audit.')
    return record

def prepare_selected_cases(project,ids,hours=72,baseline_run='R00',progress=None,cancel=None):
    """Sequential overnight preparation queue.

    QA is run independently for every selected case and failures are recorded
    without aborting later cases. Water/heat spin-up is deliberately run only
    once on the designated baseline climate so all production cases share one
    conditioned initial state. Inputs are then frozen once, using a passing QA
    case as the numerical reference.
    """
    if not ids: raise InputError('Select at least one climate case for the preparation queue.')
    require_trial(project)
    started=now(); rows=[]; n=len(ids)
    for i,rid in enumerate(ids):
        if cancel and cancel(): raise Cancelled('Preparation queue cancelled. Completed case QA records and caches are retained.')
        def qcb(fr,msg):
            if progress: progress(.88*(i+fr)/max(n,1),f'{rid} QA ({i+1}/{n}): '+msg)
        try:
            numerical_qa(project,rid,hours,qcb,cancel)
            rows.append({'run_id':rid,'qa':'PASS','error':None})
        except Cancelled:
            raise
        except Exception as exc:
            rows.append({'run_id':rid,'qa':'FAIL','error':str(exc)})
            # numerical_qa writes its structured failure record; continue.
            if progress: progress(.88*(i+1)/max(n,1),f'{rid} QA failed; saved record and continuing to the next selected case.')
    # One common baseline spin-up preserves climate-comparison comparability.
    spin_status='REUSED'; spin_error=None
    r=project.readiness()
    if not (r.get('spinup_valid') and (r.get('spinup') or {}).get('baseline_run')==baseline_run):
        try:
            if progress: progress(.90,'QA sweep finished. Starting one shared water/heat spin-up on '+baseline_run+'.')
            spinup(project,baseline_run,lambda fr,msg: progress(.90+.08*fr,'Shared spin-up: '+msg) if progress else None,cancel)
            spin_status='PASS'
        except Cancelled:
            raise
        except Exception as exc:
            spin_status='FAIL'; spin_error=str(exc)
    passing=[x['run_id'] for x in rows if x['qa']=='PASS' and valid_case_qa(project,x['run_id'])]
    freeze_status='SKIPPED'
    freeze_error=None
    ref=baseline_run if baseline_run in passing else (passing[0] if passing else None)
    if ref and spin_status!='FAIL':
        try:
            if progress: progress(.99,'Freezing reviewed inputs once; QA reference '+ref+'.')
            freeze(project,qa_reference_run=ref); freeze_status='PASS'
        except Exception as exc:
            freeze_status='FAIL'; freeze_error=str(exc)
    else:
        freeze_error=spin_error or 'No selected case produced a current passing QA record, so inputs were not frozen.'
    result={'passed':all(x['qa']=='PASS' for x in rows) and freeze_status=='PASS',
      'created_utc':now(),'started_utc':started,'selected_runs':list(ids),'baseline_spinup_run':baseline_run,
      'case_qa':rows,'qa_passed':sum(x['qa']=='PASS' for x in rows),'qa_failed':sum(x['qa']!='PASS' for x in rows),
      'spinup':spin_status,'spinup_error':spin_error,'freeze':freeze_status,'freeze_error':freeze_error,'qa_reference_run':ref,
      'method':'Selected cases are QA-tested sequentially. A failure does not stop later cases. One baseline water/heat spin-up is shared by all cases, then reviewed inputs are frozen once.'}
    atomic_json(project.path/'prep_queue.json',result)
    if progress: progress(1.,f"Preparation queue finished: {result['qa_passed']} QA passed, {result['qa_failed']} failed; spin-up {spin_status}; freeze {freeze_status}.")
    return result

def freeze(project,qa_reference_run=None):
    require_trial(project); r=project.readiness()
    if not r['spinup_valid']: raise InputError('Current inputs need a matching, converged water/heat spin-up before freezing.')
    if project.config.get('is_demo'): raise InputError('Synthetic demonstration projects cannot be certified or frozen as site input sets.')
    ref=qa_reference_run or (r.get('qa') or {}).get('run_id')
    if not ref or not valid_case_qa(project,ref):
        raise InputError('Freeze needs at least one current passing case-specific Numerical QA record. Run QA for R00 or another reference case first.')
    qref=load_case_qa(project,ref)
    runs=[]
    for case in project.config['runs']:
        c=project.climate(case['id']);runs.append({'id':case['id'],'file':case['file'],'sha256':c.source_sha256,'hours':c.n,'calendar':c.calendar})
    if not runs:raise InputError('No climate cases are registered.')
    result={'passed':True,'created_utc':now(),'physics_hash':project.physics_hash(),'run_plan_hash':digest(project.config['runs']),
      'runs':runs,'initial_state_sha256':sha256(project.path/'spinup.npz'),
      'spinup_record_sha256':sha256(project.path/'spinup.json'),'qa_reference_run':ref,
      'qa_reference_sha256':sha256(qa_case_path(project,ref)),
      'scope':'Reproducibility lock of this reduced implementation, not validation or full compliance with every method proposed in the earlier audit workbook.',
      'qa_scope':qref.get('scope'),'phase_mode':'single-salt equilibrium approximation' if project.phase() else 'transport only'}
    atomic_json(project.path/'freeze.json',result)
    return result
