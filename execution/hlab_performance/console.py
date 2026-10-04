"""Independent case scheduler with resumable original outputs (execution only)."""
from __future__ import annotations
import os
for _v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import argparse, copy, csv, http.server, io, json, multiprocessing as mp, platform
import queue, secrets, shutil, sys, threading, time, traceback, urllib.parse, webbrowser, zipfile
from collections import deque
from pathlib import Path
from . import __version__
from .scope_review import reviewed,require_review,TEXT,POLICY
from .resources import MAX_WORKERS, CasePause, resolve_workers, resources
from .safety import FileLock, SafetyError, atomic, hash_file, inspect_output, now, read_json, safe_path, stamp, verify_core, checkpoint_pair

ACTIVE = {'STARTING','VERIFYING','BENCHMARKING','RUNNING','PAUSING','FINALIZING'}
TERMINAL = {'COMPLETED','FAILED','BLOCKED','PAUSED','INTERRUPTED','CHECK_NEEDED','OTHER_VERSION'}

class Manager:
    def __init__(self, root, site):
        self.root=Path(root).resolve(); self.site=site; sys.path.insert(0,str(self.root))
        from hmsl.project import Project
        self.p=Project(safe_path(self.root/'projects',site)); verify_core(self.root)
        self.config=self.p.config; self.physics=self.p.physics_hash()
        self.ids=[r['id'] for r in self.config['runs']]
        if len(set(self.ids)) != len(self.ids): raise SafetyError('Project has duplicate run IDs.')
        self.meta={r['id']:r for r in self.config['runs']}
        self.lock=threading.RLock(); self.active=False; self.ctx=mp.get_context('spawn')
        self.pause=self.ctx.Event(); self.workers={}; self.case_pauses={}
        self.target_workers=4; self.worker_mode=4
        self.resource_state=resources(self.root)
        self.queue_path=self.p.path/'performance_queue.json'
        self.state=read_json(self.queue_path) if self.queue_path.is_file() else {
            'status':'IDLE','cases':[],'message':'Ready to verify and resume saved work.'}
        if self.state.get('status') in ('RUNNING','PAUSING','RESOURCE_WAIT'):
            self.state['status']='INTERRUPTED'
            self.state['message']='Previous console stopped. Saved output records, not old live counters, determine progress.'
        self.history_path=safe_path(self.root/'performance_runtime'/'case_status',site)
        self._receipt_cache={}
        self.inventory=[]; self.refresh()
        inv={x['run_id']:x for x in self.inventory}
        for x in self.state.get('cases',[]):
            actual=inv.get(x.get('run_id'),{})
            if actual.get('status') in ('CHECK_NEEDED','OTHER_VERSION'):
                x.update(status=actual['status'],message=actual.get('message','Saved output needs inspection.'))
            elif actual.get('status')=='COMPLETED':
                x.update(status='COMPLETED',saved_hours=actual['saved_hours'],done_hours=actual['saved_hours'])
            else:
                if x.get('status') in ACTIVE:
                    x.update(status='INTERRUPTED',message='Previous worker stopped; saved checkpoint is the restart point.')
                elif x.get('status')=='COMPLETED':
                    x.update(status='CHECK_NEEDED',message='Old queue says complete but current completion evidence is missing; inspect before running.')
                x['saved_hours']=actual.get('saved_hours',0)
                x['done_hours']=x['saved_hours']
            x.pop('remaining_seconds_estimate',None); x.pop('seconds_per_model_hour',None)
        if any(x.get('status')=='CHECK_NEEDED' for x in self.state.get('cases',[])):
            self.state['status']='FINISHED_WITH_ISSUES'

    def refresh(self):
        ready=self.p.readiness(); rows=[]
        freeze_path=self.p.path/'freeze.json'
        totals={x['id']:x.get('hours') for x in read_json(freeze_path).get('runs',[])} if freeze_path.is_file() else {}
        for r in self.config['runs']:
            entry={'run_id':r['id'],'label':r.get('model','')+' / '+r.get('scenario','')+' / '+r.get('period',''),
                   'status':'NOT_STARTED','saved_hours':0,'total_hours':totals.get(r['id']),
                   'qa_status':next((x.get('qa_status','UNKNOWN') for x in ready['runs'] if x['id']==r['id']),'UNKNOWN'), 'message':''}
            mode_root=self.p.path/'outputs'/'frozen_reduced_transport'
            for legacy in (self.p.path/'outputs',mode_root):
                if (legacy/'provenance.json').is_file():
                    pv=read_json(legacy/'provenance.json')
                    if pv.get('run_id')==r['id'] and pv.get('mode')=='FROZEN_REDUCED_TRANSPORT':
                        entry.update(status='CHECK_NEEDED',message='Nonstandard flat output folder retained. Automatic relocation or hour-zero replacement is blocked.')
            candidates=[]
            for d in mode_root.glob(r['id']+'_*'):
                if not (d/'provenance.json').is_file(): continue
                try:
                    x=inspect_output(d,deep=False)
                    if x['physics_hash']==self.physics and x['climate_sha256']==r.get('expected_sha256'):
                        x['folder']=str(d.relative_to(self.p.path)); candidates.append(x)
                    elif x['status'] in ('COMPLETED','PARTIAL'):
                        entry.update(status='OTHER_VERSION',message='Differently keyed results exist; no overwrite or new hour-zero run is allowed.')
                except Exception as exc: entry.update(status='CHECK_NEEDED',message=str(exc))
            if len(candidates)>1:
                entry.update(status='CHECK_NEEDED',message='Multiple current-identity output folders require inspection.')
            elif candidates and entry['status']!='CHECK_NEEDED': entry.update(candidates[0])
            # Retain failures/blocks across later queues, keyed to the same science
            # and climate. Actual completion/checkpoint evidence takes precedence.
            receipt_path=safe_path(self.history_path,r['id']+'.json')
            if receipt_path.is_file() and entry['status'] not in ('COMPLETED','CHECK_NEEDED','OTHER_VERSION'):
                try:
                    old=read_json(receipt_path)
                    if old.get('physics_hash')==self.physics and old.get('climate_sha256')==r.get('expected_sha256'):
                        if old.get('status') in ('FAILED','BLOCKED','PAUSED','INTERRUPTED'):
                            entry.update(status=old['status'],message=old.get('message',''),last_attempt_utc=old.get('recorded_utc'))
                        elif old.get('status')=='COMPLETED':
                            entry.update(status='CHECK_NEEDED',message='Previous receipt says complete, but current final completion evidence is missing.')
                except (OSError,ValueError): pass
            rows.append(entry)
        self.inventory=rows; self.ready=ready['frozen_ready']
        return rows

    def snapshot(self):
        with self.lock:
            self.resource_state=resources(self.root)
            s=copy.deepcopy(self.state)
            inv=copy.deepcopy(self.inventory)
            jobs={r['run_id']:r for r in s.get('cases',[])}
            combined=[dict(r,**{k:v for k,v in jobs.get(r['run_id'],{}).items() if k!='run_id'}) for r in inv]
            counts={}
            for r in combined: counts[r['status']]=counts.get(r['status'],0)+1
            running=sum(r['status'] in ACTIVE for r in combined) if self.active else 0
            rates=[1/r['seconds_per_model_hour'] for r in s.get('cases',[])
                   if self.active and r.get('status')=='RUNNING' and (r.get('seconds_per_model_hour') or 0)>0]
            return {'title':'Heritage Lab | All-cases parallel & resume','version':__version__,
                    'scope_reviewed':reviewed(self.root,self.site,self.physics),'scope_policy':POLICY,'scope_text':TEXT,
                    'reference_science_version':'1.2.0 / 1.2.1 warm/liquid core + 1.4.0 equilibrium phase extension','site':self.config['site_name'],
                    'site_id':self.site,'ready':self.ready,'active':self.active,'inventory':inv,'queue':s,
                    'cores_reported':os.cpu_count(),'resources':dict(self.resource_state),
                    'free_disk_GB':self.resource_state['free_disk_GB'],'max_workers':MAX_WORKERS,
                    'active_workers':running,'target_workers':self.target_workers,
                    'case_counts':counts,'recent_combined_model_hours_per_second':sum(rates) or None,
                    'backup_receipt':read_json(self.root/'performance_installation.json') if (self.root/'performance_installation.json').is_file() else None}

    def persist(self):
        atomic(self.queue_path,self.state)
        for r in self.state.get('cases',[]):
            if r.get('status') not in TERMINAL: continue
            rid=r['run_id']
            record={'run_id':rid,'physics_hash':self.physics,
                    'climate_sha256':self.meta.get(rid,{}).get('expected_sha256'),
                    'status':r['status'],'saved_hours':r.get('saved_hours',0),
                    'message':r.get('message',''),'session':self.state.get('session')}
            signature=json.dumps(record,sort_keys=True)
            if self._receipt_cache.get(rid)!=signature:
                record['recorded_utc']=now()
                atomic(safe_path(self.history_path,rid+'.json'),record)
                self._receipt_cache[rid]=signature

    def start(self,ids,workers=4,auto_cache=True,pore_scope_approved=False):
        with self.lock:
            if self.active: raise SafetyError('A queue is active. Use its worker limit controls, or pause safely first.')
            if not isinstance(ids,list) or any(not isinstance(i,str) for i in ids): raise SafetyError('Select registered cases.')
            ids=list(dict.fromkeys(ids))
            if not ids or any(i not in self.ids for i in ids): raise SafetyError('Select registered case IDs.')
            nworkers=resolve_workers(workers,len(ids))
            verify_core(self.root); self.refresh()
            if not self.ready: raise SafetyError('Matching frozen preparation is not valid; no approval records were changed.')
            require_review(self.root,self.site,self.physics,approve=pore_scope_approved)
            inventory={x['run_id']:x for x in self.inventory}
            ids.sort(key=lambda i:(0 if inventory[i].get('saved_hours',0)>0 and inventory[i]['status']!='COMPLETED' else 1,self.ids.index(i)))
            session_dir=self.root/'performance_runtime'/'sessions'/stamp(); session_dir.mkdir(parents=True)
            if self.queue_path.is_file():
                atomic(session_dir/'previous_queue.json',read_json(self.queue_path))
            rows=[]
            for i in ids:
                r=dict(inventory[i]); r['done_hours']=r['saved_hours']
                if r['status']=='COMPLETED': r['message']='Completed output retained; not recalculated.'
                elif r['status'] not in ('CHECK_NEEDED','OTHER_VERSION'):
                    r.update(status='PENDING',message='Queued for independent verification and resume.')
                rows.append(r)
            self.target_workers=nworkers; self.worker_mode=workers
            self.state={'status':'RUNNING','started_utc':now(),'updated_utc':now(),
                        'session':str(session_dir.relative_to(self.root)),'workers_requested':nworkers,
                        'worker_mode':workers,'auto_cache':bool(auto_cache),'selected_runs':ids,'cases':rows,
                        'message':'Each eligible case runs independently; failures and blocks do not stop other cases.',
                        'execution_extension_version':__version__}
            self.pause.clear(); self.case_pauses={}; self.workers={}; self.active=True; self.persist()
            self.thread=threading.Thread(target=self._supervise,args=(session_dir,bool(auto_cache)),daemon=True)
            self.thread.start()

    def set_workers(self,value):
        with self.lock:
            if not self.active: raise SafetyError('Start a queue first, or choose the limit before starting.')
            self.target_workers=resolve_workers(value,len(self.state['selected_runs']))
            self.worker_mode=value
            self.state.update(workers_requested=self.target_workers,worker_mode=value,updated_utc=now())
            self.state['message']='Worker limit updated. Raising it fills available slots; lowering it does not kill active workers. Pause individual cases to release slots safely.'
            self.persist()

    def pause_case(self,rid):
        with self.lock:
            if not self.active: raise SafetyError('There is no active queue.')
            rows={x['run_id']:x for x in self.state['cases']}
            if rid not in rows: raise SafetyError('This case is not in the queue.')
            r=rows[rid]
            if rid in self.case_pauses and r['status'] in ACTIVE:
                self.case_pauses[rid].set()
                r.update(status='PAUSING',message='This case will stop after its next committed checkpoint; other cases continue.')
            elif r['status']=='PENDING':
                r.update(status='PAUSED',message='Held before launch; no new calculation was started.')
            else: raise SafetyError('This case is not running or waiting.')
            self.persist()

    def request_pause(self):
        with self.lock:
            if self.active:
                self.pause.set()
                self.state.update(status='PAUSING',message='Waiting for active workers to commit their next chunks. Pending cases will not start.')
                self.persist()

    def _supervise(self,session_dir,auto_cache):
        from .worker import run_case
        events=self.ctx.Queue(); pending=deque(x['run_id'] for x in self.state['cases'] if x['status']=='PENDING')
        rows={x['run_id']:x for x in self.state['cases']}
        last_save=0.; last_resource=0.; next_launch=0.; resource=resources(self.root)
        def accept(msg):
            rid=msg.get('run_id')
            if rid in rows and not msg.get('worker_exit'):
                with self.lock: rows[rid].update(msg)
        try:
            while pending or self.workers:
                now_mono=time.monotonic()
                if now_mono-last_resource>=1:
                    resource=resources(self.root); last_resource=now_mono
                    with self.lock: self.resource_state=resource
                    if resource['critical']:
                        self.pause.set()
                        with self.lock: self.state['message']=resource['critical_reason']+' Finishing committed chunks, not killing workers.'
                # Stagger spawn slightly to avoid a simultaneous import/I/O burst.
                if pending and not self.pause.is_set() and resource['can_launch'] and now_mono>=next_launch:
                    with self.lock:
                        if len(self.workers)<self.target_workers:
                            rid=pending.popleft()
                            if rows[rid]['status']=='PENDING':
                                event=self.ctx.Event(); self.case_pauses[rid]=event
                                proc=self.ctx.Process(target=run_case,args=(str(self.root),self.site,rid,events,
                                      CasePause(self.pause,event),str(session_dir),auto_cache),daemon=True)
                                try:
                                    proc.start(); self.workers[rid]=(proc,time.monotonic())
                                    rows[rid].update(status='STARTING',worker_pid=proc.pid,message='Starting isolated worker and verifying existing files.')
                                except Exception as exc:
                                    self.case_pauses.pop(rid,None)
                                    rows[rid].update(status='FAILED',message='Worker could not start: '+str(exc))
                                next_launch=now_mono+.15
                try: accept(events.get(timeout=.10))
                except queue.Empty: pass
                # Drain a bounded burst so 21 workers do not build a progress backlog.
                for _ in range(256):
                    try: accept(events.get_nowait())
                    except queue.Empty: break
                for rid,(proc,born) in list(self.workers.items()):
                    if not proc.is_alive():
                        proc.join()
                        # Child Queue feeder is flushed before normal exit.
                        while True:
                            try: accept(events.get_nowait())
                            except queue.Empty: break
                        with self.lock:
                            del self.workers[rid]; self.case_pauses.pop(rid,None)
                            if rows[rid]['status'] not in TERMINAL:
                                rows[rid].update(status='INTERRUPTED',message=f'Worker exited with code {proc.exitcode}; saved checkpoint remains authoritative.')
                        proc.close()
                if self.pause.is_set() and not self.workers: break
                if now_mono-last_save>=1:
                    with self.lock:
                        self.state['updated_utc']=now()
                        self.state['status']='PAUSING' if self.pause.is_set() else ('RESOURCE_WAIT' if pending and not resource['can_launch'] else 'RUNNING')
                        self.state['resource_hold']=resource['launch_hold_reason'] if pending and not resource['can_launch'] else ''
                        self.persist()
                    last_save=now_mono
            with self.lock:
                failed=[r['run_id'] for r in rows.values() if r['status'] in ('FAILED','BLOCKED','INTERRUPTED','CHECK_NEEDED','OTHER_VERSION')]
                unfinished=any(r['status'] in ('PAUSED','PENDING') for r in rows.values())
                self.state.update(status='PAUSED' if self.pause.is_set() or unfinished else ('FINISHED_WITH_ISSUES' if failed else 'COMPLETED'),
                                  finished_utc=now(),failed_or_blocked=failed,resource_hold='',
                                  message='Queue paused; keep all output/checkpoint files and resume remaining cases.' if self.pause.is_set() or unfinished else 'Queue finished. Completed, blocked and failed cases are reported separately.')
                self.persist()
        except Exception as exc:
            self.pause.set()
            with self.lock: self.state.update(status='FAILED',message=str(exc)); self.persist()
            (session_dir/'supervisor_error.log').write_text(traceback.format_exc(),encoding='utf-8')
        finally:
            for proc,_ in list(self.workers.values()): proc.join()
            events.close(); events.join_thread()
            try:
                with self.lock: self.workers={}; self.case_pauses={}; self.refresh()
            except Exception: pass
            finally:
                with self.lock: self.active=False

    def export_csv(self):
        s=self.snapshot(); jobs={r['run_id']:r for r in s['queue'].get('cases',[])}
        fields=['run_id','label','status','done_hours','total_hours','saved_hours','saved_utc','seconds_per_model_hour','remaining_seconds_estimate','message']
        out=io.StringIO(newline=''); w=csv.DictWriter(out,fieldnames=fields,extrasaction='ignore'); w.writeheader()
        for base in s['inventory']:
            r=dict(base,**{k:v for k,v in jobs.get(base['run_id'],{}).items() if k!='run_id'})
            # Defend spreadsheets from formula injection in labels/errors.
            for key,v in r.items():
                if isinstance(v,str) and v.startswith(('=','+','-','@')): r[key]="'"+v
            w.writerow(r)
        return out.getvalue()

    def diagnostic(self,rid):
        if rid not in self.ids: raise SafetyError('Unknown run ID.')
        folders=[]
        for d in (self.p.path/'outputs'/'frozen_reduced_transport').glob(rid+'_*'):
            if (d/'provenance.json').is_file() and read_json(d/'provenance.json').get('physics_hash')==self.physics: folders.append(d)
        if len(folders)!=1: raise SafetyError('This case needs exactly one current-version output folder.')
        d=folders[0]; cp,state=checkpoint_pair(d); buf=io.BytesIO()
        with zipfile.ZipFile(buf,'w',zipfile.ZIP_DEFLATED) as z:
            z.writestr('checkpoint.json',json.dumps(cp,indent=2)); z.writestr(cp['state_file'],state)
            for n in ['phase_extension_v140.json','phase_extension_summary.json','phase_diagnostics.csv.gz','cold_scope_policy.json','cold_scope_summary.json','cold_scope_report.html','provenance.json','resolved_project.json','summary.json']:
                if (d/n).is_file(): z.write(d/n,n)
            for n in ['project.json','freeze.json','spinup.json','performance_queue.json']:
                if (self.p.path/n).is_file(): z.write(self.p.path/n,n)
            for pre in ('hours','profile'):
                path=d/'chunks'/f'{pre}_{cp["completed_hours"]:09d}.csv'
                if path.is_file(): z.write(path,'last_chunk/'+path.name)
            z.writestr('machine.json',json.dumps({'python':sys.version,'platform':platform.platform(),
                  'resources':resources(self.root),'diagnostic_saved_utc':now(),
                  'note':'Latest state/chunk snapshot, not a full project backup.'},indent=2))
        return buf.getvalue()


def serve(root,site,port=8765,open_browser=True):
    root=Path(root).resolve(); token=secrets.token_urlsafe(32)
    manager=Manager(root,site)
    class Handler(http.server.BaseHTTPRequestHandler):
        def log_message(self,*args): pass
        def respond(self,code,obj,mime='application/json',name=None):
            data=json.dumps(obj,allow_nan=False).encode() if isinstance(obj,(list,dict)) else (obj.encode() if isinstance(obj,str) else obj)
            self.send_response(code); self.send_header('Content-Type',mime); self.send_header('Content-Length',str(len(data)))
            self.send_header('Cache-Control','no-store'); self.send_header('X-Content-Type-Options','nosniff')
            if name: self.send_header('Content-Disposition','attachment; filename="'+name+'"')
            self.end_headers()
            try: self.wfile.write(data)
            except (BrokenPipeError,ConnectionResetError): pass
        def valid(self): return self.headers.get('Host','').split(':')[0] in ('127.0.0.1','localhost')
        def do_GET(self):
            if not self.valid(): return self.respond(403,{'error':'Localhost only.'})
            u=urllib.parse.urlparse(self.path); q=urllib.parse.parse_qs(u.query)
            try:
                if u.path=='/': return self.respond(200,Path(__file__).with_name('dashboard.html').read_text(encoding='utf-8').replace('__TOKEN__',token),'text/html; charset=utf-8')
                if u.path=='/api/state': return self.respond(200,manager.snapshot())
                if u.path=='/api/queue.json': return self.respond(200,manager.snapshot()['queue'],name='performance_queue.json')
                if u.path=='/api/progress.csv': return self.respond(200,manager.export_csv(),'text/csv; charset=utf-8','case_progress.csv')
                if u.path=='/api/diagnostic':
                    rid=q.get('run_id',[''])[0]
                    return self.respond(200,manager.diagnostic(rid),'application/zip',rid+'_progress_diagnostic.zip')
                return self.respond(404,{'error':'Not found'})
            except Exception as exc: return self.respond(400,{'error':str(exc)})
        def do_POST(self):
            if not self.valid() or not secrets.compare_digest(self.headers.get('X-Token',''),token):
                return self.respond(403,{'error':'Reload this dashboard to refresh local authorization.'})
            try:
                n=int(self.headers.get('Content-Length','0'))
                if n<0 or n>65536: raise SafetyError('Request too large.')
                p=json.loads(self.rfile.read(n) or b'{}')
                if not isinstance(p,dict): raise SafetyError('Invalid request.')
                if self.path=='/api/start': manager.start(p.get('ids',[]),p.get('workers',4),p.get('auto_cache',True),p.get('pore_scope_approved',False))
                elif self.path=='/api/workers': manager.set_workers(p.get('workers'))
                elif self.path=='/api/pause': manager.request_pause()
                elif self.path=='/api/pause_case': manager.pause_case(p.get('run_id'))
                elif self.path=='/api/refresh':
                    if manager.active: raise SafetyError('Full refresh is available when paused. Live progress already updates automatically.')
                    manager.refresh()
                else: raise SafetyError('Unknown action.')
                self.respond(200,manager.snapshot())
            except Exception as exc: self.respond(400,{'error':str(exc)})
    with FileLock(root/'performance_runtime'/'console.lock'):
        server=http.server.ThreadingHTTPServer(('127.0.0.1',int(port)),Handler)
        address=f'http://127.0.0.1:{server.server_port}'
        print('Heritage Lab All-cases Parallel '+__version__+': '+address,flush=True)
        print('Keep this terminal open. Use Pause ALL safely before closing. Do not run the old console at the same time.',flush=True)
        if open_browser: webbrowser.open(address)
        if os.name=='nt':
            try: ctypes_execution_state(0x80000001)
            except Exception: pass
        try: server.serve_forever(poll_interval=.3)
        except KeyboardInterrupt:
            print('Pausing after the next saved checkpoint. Please wait...',flush=True); manager.request_pause()
            while manager.active: time.sleep(.5)
        finally:
            server.server_close()
            if os.name=='nt':
                try: ctypes_execution_state(0x80000000)
                except Exception: pass

def ctypes_execution_state(value):
    import ctypes
    return ctypes.windll.kernel32.SetThreadExecutionState(value)

if __name__=='__main__':
    mp.freeze_support()
    ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); ap.add_argument('--site',default='mohenjo_daro')
    ap.add_argument('--port',type=int,default=8765); ap.add_argument('--no-browser',action='store_true'); args=ap.parse_args()
    serve(args.root,args.site,args.port,not args.no_browser)
