"""One independent case per process; preserves committed work and adds v1.4 phase physics only for unfinished work."""
from __future__ import annotations
import os
for _v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_v]='1'
import json, shutil, sys, time, traceback
from collections import deque
from pathlib import Path
from .safety import SafetyError, FileLock, atomic, now, stamp, read_json, verify_core, inspect_output, mirror_checkpoint, hash_file, checkpoint_pair


def _normalize_phase_exports(folder, legacy_prefix_hours):
    """Add explicit zero-ice fields to retained legacy prefix in final merged exports.

    Old chunk files are never rewritten.  The merged convenience exports are
    normalized after all hours exist, because v1.3 had already certified the
    retained prefix as liquid-stable (no modeled ice state).
    """
    import pandas as pd
    folder=Path(folder)
    hp=folder/'hourly_results.csv.gz'
    out={'legacy_prefix_hours':int(legacy_prefix_hours)}
    if hp.is_file():
        h=pd.read_csv(hp)
        for dst,src in [('total_theta_min_m3_m3','theta_min_m3_m3'),('total_theta_mean_m3_m3','theta_mean_m3_m3'),('total_theta_max_m3_m3','theta_max_m3_m3')]:
            if dst not in h:h[dst]=h[src]
            else:h[dst]=h[dst].fillna(h[src])
        for col in ('ice_theta_mean_m3_m3','ice_theta_max_m3_m3','ice_water_inventory_kg_m2_footprint'):
            if col not in h:h[col]=0.0
            else:h[col]=h[col].fillna(0.0)
        if 'phase_active_cells' not in h:h['phase_active_cells']=0
        else:h['phase_active_cells']=h['phase_active_cells'].fillna(0).astype(int)
        if 'phase_model' not in h:h['phase_model']='legacy-liquid-prefix'
        else:h['phase_model']=h['phase_model'].fillna('legacy-liquid-prefix')
        h.to_csv(hp,index=False,compression={'method':'gzip','mtime':0})
        out.update(phase_active_hours=int((h.phase_active_cells>0).sum()),
                   maximum_ice_theta_m3_m3=float(h.ice_theta_max_m3_m3.max()),
                   maximum_ice_water_inventory_kg_m2_footprint=float(h.ice_water_inventory_kg_m2_footprint.max()),
                   minimum_material_temperature_C=float(h.T_min_C.min()))
        diag=[c for c in ['forcing_timestamp_utc','elapsed_hours','T_min_C','theta_mean_m3_m3','total_theta_mean_m3_m3','ice_theta_mean_m3_m3','ice_theta_max_m3_m3','ice_water_inventory_kg_m2_footprint','phase_active_cells','sulfate_cmax_mol_m3_liquid','phase_model'] if c in h]
        h[diag].to_csv(folder/'phase_diagnostics.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    pp=folder/'profiles.csv.gz'
    if pp.is_file():
        p=pd.read_csv(pp)
        if 'total_theta_m3_m3' not in p:p['total_theta_m3_m3']=p['theta_m3_m3']
        else:p['total_theta_m3_m3']=p['total_theta_m3_m3'].fillna(p['theta_m3_m3'])
        if 'ice_theta_m3_m3' not in p:p['ice_theta_m3_m3']=0.0
        else:p['ice_theta_m3_m3']=p['ice_theta_m3_m3'].fillna(0.0)
        p.to_csv(pp,index=False,compression={'method':'gzip','mtime':0})
    return out


def run_case(root,site,rid,events,pause,session_dir,auto_cache=True):
    """Windows-spawn-safe worker.  Completed cases are never recalculated."""
    root=Path(root).resolve();session_dir=Path(session_dir);sys.path.insert(0,str(root))
    from hmsl.project import Project
    from hmsl.common import Cancelled,digest
    from hmsl.workflow import load_state,run_transport
    import hmsl.workflow as workflow
    original_model=workflow.ColumnModel; original_write_report=workflow.write_report
    started=time.perf_counter();folder=None;receipt=None;ready_to_pause=False
    last_sent=0.;initial_done=0;durable_done=0;total=None;recent=deque(maxlen=121)
    def send(**kw):events.put(dict(run_id=rid,updated_utc=now(),**kw))
    def cancelled():return bool(ready_to_pause)
    try:
        with FileLock(root/'performance_runtime'/'locks'/(site+'_'+rid+'.lock')):
            verify_core(root)
            if pause.is_set():raise Cancelled('Paused before starting this case.')
            send(status='VERIFYING',message='Checking frozen inputs, source identity, committed checkpoint and phase-extension policy.')
            p=Project(root/'projects'/site)
            if not p.readiness()['frozen_ready']:raise SafetyError('Existing spin-up/freeze is not valid for this configuration. No settings or QA records were changed.')
            c=p.climate(rid);total=c.n
            original_climate=p.climate;p.climate=lambda case:c if case==rid else original_climate(case)
            physics=p.physics_hash()
            from .scope_review import require_review
            scope_review=require_review(root,site,physics)
            key=digest({'physics':physics,'climate':c.source_sha256,'calendar':c.calendar,'hours':c.n,'mode':'FROZEN_REDUCED_TRANSPORT'})
            folder=p.path/'outputs'/'frozen_reduced_transport'/(rid+'_'+key[:12]);parent=folder.parent
            other=[d for d in parent.glob(rid+'_*') if d!=folder and (d/'checkpoint.json').is_file()] if parent.is_dir() else []
            if other and not (folder/'checkpoint.json').is_file():raise SafetyError('A differently keyed previous output exists for '+rid+'. No restart at hour zero was allowed.')
            if (folder/'provenance.json').is_file():
                audit=inspect_output(folder,deep=True)
                if audit['physics_hash']!=physics or audit['climate_sha256']!=c.source_sha256:raise SafetyError('Existing provenance does not match the current run identity.')
                if audit.get('key') and audit['key']!=key:raise SafetyError('Checkpoint key does not match this exact case/duration/mode.')
                initial_done=durable_done=audit['saved_hours']
                if initial_done>total:raise SafetyError('Checkpoint extends beyond this climate case.')
                if audit['status']=='COMPLETED':
                    send(status='COMPLETED',saved_hours=initial_done,done_hours=initial_done,total_hours=total,output=str(folder.relative_to(p.path)),message='Already completed; verified output retained and not recalculated.')
                    return
            elif folder.exists() and any(folder.iterdir()):raise SafetyError('Output folder has files but no provenance. No overwrite or new run was allowed.')
            if not workflow.valid_case_qa(p,rid):raise SafetyError(rid+' has no current passing case-specific QA record. No production run was started.')
            if pause.is_set():raise Cancelled('Paused after validation; no simulation hour was advanced.')
            if initial_done:
                cp,_=checkpoint_pair(folder);state=load_state(folder/cp['state_file']);mirror_checkpoint(folder,root/'performance_recovery'/site/rid)
            else:state=load_state(p.path/'spinup.npz')
            # Once v1.4 starts an unfinished case, the checkpoint theta variable is
            # total water-equivalent.  At migration it is safe because v1.3 only
            # committed liquid-stable prefixes before a PHASE_REVIEW_REQUIRED stop.
            from .phase_change import EquilibriumFreezeThawModel, VERSION as PHASE_VERSION, SOURCES as PHASE_SOURCES
            workflow.ColumnModel=EquilibriumFreezeThawModel
            folder.mkdir(parents=True,exist_ok=True)
            marker=folder/'phase_extension_v140.json'
            phase_meta={'phase_extension':'Heritage Lab 1.4.0','phase_model':PHASE_VERSION,'started_utc':now(),'run_id':rid,
                'legacy_prefix_hours':initial_done,'checkpoint_migration':'At the first v1.4 resume, legacy theta equals total water because the retained v1.3 prefix was accepted only in the liquid-stable domain. Subsequent v1.4 checkpoint theta stores total condensable water as liquid-water-equivalent volume.',
                'warm_path':'Calls the original derivative verbatim whenever no equilibrium ice is present.','sources':PHASE_SOURCES,
                'limitations':['Local thermodynamic equilibrium; no nucleation/supercooling hysteresis.','Pure-water phase relation; Na2SO4 freezing-point depression/activity is not credited.','No ice-pressure damage, pore-volume expansion mechanics, or explicit snow accumulation.','Thermal conductivity remains the source material function evaluated at liquid water content; no separate ice-conductivity mixing law is introduced.'],
                'policy_review':scope_review}
            atomic(marker,phase_meta)
            receipt=folder/('performance_execution_'+stamp()+'.json')
            data={'extension':'Heritage performance/resume 1.4.0 equilibrium phase extension','started_utc':now(),'run_id':rid,'original_physics_hash':physics,'original_checkpoint_key':key,'resume_hour':initial_done,'status':'RUNNING','original_scientific_files':verify_core(root),'extension_sources':{f.name:hash_file(f) for f in Path(__file__).parent.glob('*.py')},'threads_per_process':1,'engine':'equilibrium_phase_reference','legacy_warm_prefix_hours':initial_done,'scope_review':scope_review,'note':'Original hmsl core files are byte-identical. Effective equations are extended only when equilibrium ice is present; this is not an unchanged-equations claim.'}
            atomic(receipt,data);atomic(session_dir/(rid+'_phase.json'),phase_meta)
            send(status='RUNNING',done_hours=initial_done,saved_hours=initial_done,total_hours=total,engine='equilibrium_phase_reference',message=(f'Resuming from saved hour {initial_done:,}; committed prefix retained. ' if initial_done else 'Starting unfinished case. ')+'Equilibrium liquid/ice partition and fusion latent heat are active only when needed.')
            def progress(fr,msg):
                nonlocal ready_to_pause,last_sent,durable_done
                current=time.perf_counter();done=max(initial_done,min(total,int(round(fr*total))))
                is_commit='water closure' in msg or msg.startswith('Run complete')
                if is_commit and (folder/'checkpoint.json').exists():
                    durable_done=read_json(folder/'checkpoint.json')['completed_hours'];mirror_checkpoint(folder,root/'performance_recovery'/site/rid)
                    if pause.is_set():ready_to_pause=True
                if current-last_sent<1 and not is_commit:return
                last_sent=current;recent.append((current,done));ds=recent[-1][1]-recent[0][1];dt=recent[-1][0]-recent[0][0]
                rate=dt/ds if ds>0 and dt>0 else None;remaining=rate*(total-done) if rate else None
                if shutil.disk_usage(root).free<512*1024**2:pause.set()
                send(status='PAUSING' if pause.is_set() else ('FINALIZING' if done>=total else 'RUNNING'),done_hours=done,saved_hours=durable_done,total_hours=total,seconds_per_model_hour=rate,remaining_seconds_estimate=remaining,elapsed_seconds=current-started,engine='equilibrium_phase_reference',output=str(folder.relative_to(p.path)),message=('All requested hours are saved; finalizing exports. ' if done>=total else ('Pause requested: finishing this checkpoint chunk safely. ' if pause.is_set() else ''))+msg)
            def report_with_phase(folder_arg,summary,daily,provenance):
                phase_summary=_normalize_phase_exports(folder_arg,initial_done)
                summary.update({'water_phase_model':'Local-equilibrium liquid/ice partition derived from source retention curve; fusion latent heat included.','phase_extension_version':'1.4.0 / '+PHASE_VERSION,'legacy_liquid_prefix_hours':initial_done,**phase_summary})
                provenance['phase_extension']=phase_meta
                provenance['solver_scope_note']='The original core files remain byte-identical, but cold phase-active states use the v1.4 equilibrium phase equations. Warm/no-ice derivative evaluations follow the original arithmetic path.'
                limitations=[x for x in provenance.get('limitations',[]) if not x.startswith('No ice,') and 'phase-change enthalpy beyond evaporation/condensation' not in x]
                limitations += phase_meta['limitations'] + ['Water/ice phase change is modeled for transport/energy only; no freeze-thaw mechanical damage index is produced.']
                provenance['limitations']=limitations
                original_write_report(folder_arg,summary,daily,provenance)
                atomic(Path(folder_arg)/'phase_extension_summary.json',phase_summary|{'phase_model':PHASE_VERSION,'legacy_prefix_hours':initial_done,'sources':PHASE_SOURCES})
            workflow.write_report=report_with_phase
            result=run_transport(p,rid,require_frozen=True,resume=True,progress=progress,cancel=cancelled)
            final=inspect_output(folder,deep=True);data.update(status='COMPLETED',finished_utc=now(),saved_hours=final['saved_hours']);atomic(receipt,data)
            phase_meta.update(status='COMPLETED',finished_utc=now(),saved_hours=final['saved_hours']);atomic(marker,phase_meta)
            send(status='COMPLETED',saved_hours=final['saved_hours'],done_hours=final['saved_hours'],total_hours=total,output=result,elapsed_seconds=time.perf_counter()-started,message='Completed with the equilibrium pore-water phase extension; retained prefix and final files verified.')
    except Cancelled as exc:
        saved=read_json(folder/'checkpoint.json')['completed_hours'] if folder and (folder/'checkpoint.json').is_file() else initial_done
        if receipt and receipt.is_file():
            data=read_json(receipt);data.update(status='PAUSED',finished_utc=now(),saved_hours=saved);atomic(receipt,data)
        send(status='PAUSED',saved_hours=saved,done_hours=saved,total_hours=total,message=str(exc))
    except Exception as exc:
        log=session_dir/(rid+'_error.log');log.write_text(traceback.format_exc(),encoding='utf-8')
        if receipt and receipt.is_file():
            data=read_json(receipt);data.update(status='FAILED',finished_utc=now(),error=str(exc));atomic(receipt,data)
        saved=durable_done
        try:
            if folder and (folder/'checkpoint.json').is_file():saved=read_json(folder/'checkpoint.json')['completed_hours']
        except Exception:pass
        send(status='FAILED',saved_hours=saved,done_hours=saved,total_hours=total,message=str(exc),error=str(exc),log=log.name,elapsed_seconds=time.perf_counter()-started)
    finally:
        workflow.ColumnModel=original_model;workflow.write_report=original_write_report
        events.put({'run_id':rid,'worker_exit':True,'updated_utc':now()})
