#!/usr/bin/env python3
"""Local-only dashboard. No cloud account, database or Streamlit installation.
Bind to 127.0.0.1 only. Mutating requests require a per-process CSRF token.
"""
from __future__ import annotations
import argparse, hashlib, http.server, json, mimetypes, os, secrets, shutil, sys, threading, time, traceback, urllib.parse, webbrowser, zipfile
from pathlib import Path
from hmsl import __version__
# Small sparse solves should not spawn a large BLAS thread pool on a laptop.
for _name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[_name]=os.environ.get('HMSL_NUM_THREADS','1')
from hmsl.preflight import check_cases,require_supported
from hmsl.common import InputError, Cancelled, atomic_json, safe_name, safe_path, now, load_json, sha256
from hmsl.project import Project, new_project, download_selected_material, SOURCE_MEMBER
from hmsl.climate import read_climate
from hmsl.material import read_material
from hmsl.solver import Solubility
from hmsl.workflow import screen_run, run_transport, numerical_qa, spinup, freeze, prepare_selected_cases

ROOT=Path(__file__).resolve().parent
TOKEN=secrets.token_urlsafe(32)
JOBS_LOCK=threading.Lock(); CANCEL=threading.Event()
JOB={'status':'idle','progress':0.,'message':'No active task','outputs':[]}

def project(s): return Project(safe_path(ROOT/'projects',safe_name(s)))

def update(**kw):
    with JOBS_LOCK: JOB.update(kw)

def launch_job(payload):
    global JOB
    p=project(payload['project']); action=payload['action']; ids=payload.get('ids',[])
    with JOBS_LOCK:
        if JOB['status']=='running': raise InputError('A task is already running. Wait or cancel it before starting another.')
        JOB={'status':'running','progress':0.,'message':'Starting...','outputs':[],'project':p.config['site_id'],'action':action,'started':now()}
    CANCEL.clear()
    def work():
        outputs=[]
        def cb(fr,msg): update(progress=float(fr),message=msg)
        try:
            if action=='source': result=download_selected_material(p);update(message=result['message'])
            elif action=='preflight':
                result=check_cases(p,progress=cb,cancel=CANCEL.is_set);n=sum(not x['full_case_supported'] for x in result['cases']);update(message=f'Compatibility check complete: {n} full climate cases require unsupported physics or input repair. Open the compatibility report.' if n else 'Selected files are compatible with the basic temperature scope; QA and spin-up are still required.')
            elif action=='validate':
                result=p.readiness(thorough=True);atomic_json(p.path/'data_validation.json',result)
                errors=[r for r in result['runs'] if r['status'] in ('MISSING','ERROR')]
                if errors: raise InputError(f'{len(errors)} forcing files are missing or invalid. See data_validation.json / Data tab.')
                update(message='All registered climate files passed hourly/calendar/checksum validation. Material and review status remain separate.')
            elif action=='freeze':
                result=freeze(p);update(message='Reduced-model inputs frozen. This is a reproducibility lock, not site validation.')
            elif action=='prep_queue':
                if not ids: raise InputError('Select at least one climate case for the overnight preparation queue.')
                result=prepare_selected_cases(p,ids,int(payload.get('hours',72)),payload.get('baseline_id','R00'),cb,CANCEL.is_set)
                update(message=f"Overnight preparation finished all automatic stages: {result['qa_passed']} QA passed, {result['qa_failed']} failed; water/heat spin-up {result['spinup']}; freeze {result['freeze']}. Failed QA cases were recorded and did not stop later cases.")
            elif action in ('qa','spinup'):
                if len(ids)!=1: raise InputError('Select exactly one climate case for QA or spin-up.')
                if action=='qa': result=numerical_qa(p,ids[0],int(payload.get('hours',72)),cb,CANCEL.is_set)
                else: result=spinup(p,ids[0],cb,CANCEL.is_set)
                update(message=action+' passed; record saved.')
            elif action in ('screen','trial','batch'):
                if not ids: raise InputError('Select at least one climate case.')
                if action=='batch' and not payload.get('confirm_full'): raise InputError('Confirm that you want full-length, frozen reduced-model runs.')
                queue_rows=[]
                for idx,rid in enumerate(ids):
                    if CANCEL.is_set(): raise Cancelled('Cancelled; completed cases and checkpoints are retained.')
                    def one_cb(fr,msg): cb((idx+fr)/len(ids),f'{rid} ({idx+1}/{len(ids)}): '+msg)
                    try:
                        if action=='batch': require_supported(p,[rid],progress=None,cancel=CANCEL.is_set)
                        if action=='screen': rel=screen_run(p,rid,one_cb,CANCEL.is_set)
                        else: rel=run_transport(p,rid,max_hours=int(payload.get('hours',72)) if action=='trial' else None,
                           require_frozen=action=='batch',progress=one_cb,cancel=CANCEL.is_set)
                        outputs.append(rel);queue_rows.append({'run_id':rid,'status':'PASS','output':rel,'error':None});update(outputs=outputs.copy())
                    except Cancelled:
                        raise
                    except Exception as exc:
                        if action!='batch': raise
                        queue_rows.append({'run_id':rid,'status':'FAIL','output':None,'error':str(exc)})
                        errdir=p.path/'queue_errors';errdir.mkdir(exist_ok=True)
                        (errdir/(safe_name(rid)+'_production.log')).write_text(traceback.format_exc(),encoding='utf-8')
                        one_cb(1.,'FAILED: '+str(exc)+' Continuing to next selected case.')
                if action=='batch':
                    atomic_json(p.path/'production_queue.json',{'created_utc':now(),'selected_runs':list(ids),'cases':queue_rows,
                      'passed':sum(x['status']=='PASS' for x in queue_rows),'failed':sum(x['status']!='PASS' for x in queue_rows),
                      'method':'Selected full cases run sequentially. A failed case is recorded and does not stop later selected cases.'})
                    update(message=f"Production queue finished: {sum(x['status']=='PASS' for x in queue_rows)} passed, {sum(x['status']!='PASS' for x in queue_rows)} failed.")
            else: raise InputError('Unknown task.')
            update(status='complete',progress=1.,finished=now(),outputs=outputs)
        except Cancelled as e: update(status='cancelled',message=str(e),finished=now(),outputs=outputs)
        except Exception as e:
            log=p.path/'last_error.log';log.write_text(traceback.format_exc(),encoding='utf-8')
            update(status='error',message=str(e),finished=now(),outputs=outputs)
        with JOBS_LOCK: atomic_json(p.path/'last_job.json',dict(JOB))
    threading.Thread(target=work,daemon=True).start()

def import_file(p,kind,name,temp,query):
    if p.file('freeze.json').exists(): raise InputError('This project is frozen. Clone it before importing different files.')
    name=safe_name(name)
    if name.lower().endswith('.zip'):
        count=0
        with zipfile.ZipFile(temp) as z:
            entries=z.infolist()
            if len(entries)>5000 or sum(x.file_size for x in entries)>3*1024**3: raise InputError('Archive is too large or contains too many files.')
            # Read into a staging directory. No archive member path is used directly.
            for item in entries:
                if item.is_dir():continue
                base=Path(item.filename.replace('\\','/')).name
                if kind=='material' and base!=SOURCE_MEMBER:continue
                if kind=='climate' and not (base.endswith(('.csv','.csv.gz')) and ('hourly_' in base or 'hourly_model_ready' in base) and 'disaggregation' not in base):continue
                if kind not in ('material','climate'):raise InputError('ZIP import is available only for climate or the selected brick material archive.')
                if item.file_size>600*1024**2:raise InputError('One archive member is too large.')
                dst=p.path/('materials' if kind=='material' else 'climate')/safe_name(base)
                with z.open(item) as f: data=f.read()
                if dst.exists():
                    if hashlib.sha256(data).hexdigest()==sha256(dst):count+=1;continue
                    raise InputError(f'{base} already exists with different bytes. No overwrite was performed.')
                dst.write_bytes(data); count+=1
        if not count: raise InputError('No matching files were found. Import an individual climate file or the selected material archive.')
        if kind=='material':p.material()
        return f'{count} archive files imported or already present. Registered Mohenjo filenames are linked automatically. For a new site, import individual climate files to register cases.'
    target=p.path/('climate' if kind=='climate' else 'materials')/name
    target.parent.mkdir(exist_ok=True)
    if target.exists():raise InputError('A file with this name already exists. Keep frozen source bytes or import a differently named file into a new project.')
    allowed={'climate':('.csv','.csv.gz','.xlsx'),'material':('.m6','.json'),'solubility':('.csv',)}
    if kind not in allowed or not name.lower().endswith(allowed[kind]):raise InputError('File extension is not supported for the selected import type.')
    shutil.move(temp,target);cfg=p.config
    try:
        if kind=='climate':
            cal=query.get('calendar',['gregorian'])[0];c=read_climate(target,cal)
            if not any(Path(r['file']).name==name for r in cfg['runs']):
                existing={r['id'] for r in cfg['runs']}; num=1
                while f'U{num:02d}' in existing:num+=1
                cfg['runs'].append({'id':f'U{num:02d}','model':query.get('model',['user'])[0],
                  'scenario':query.get('scenario',['custom'])[0],'period':query.get('period',['user supplied'])[0],
                  'calendar':cal,'file':'climate/'+name,'expected_sha256':c.source_sha256})
        elif kind=='material':
            m=read_material(target,cfg.get('retention_branch','de'))
            if not cfg.get('is_demo') and m.is_demo:raise InputError('Synthetic material is only allowed in the separate demo project.')
            cfg['material_file']='materials/'+name
        else:
            Solubility(target);cfg['salt']['solubility_file']='materials/'+name
        p.save(cfg)
    except Exception:
        target.unlink(missing_ok=True)
        raise
    return 'Imported and validated. Source bytes are retained; review settings before transport.'

class Handler(http.server.BaseHTTPRequestHandler):
    def log_message(self,fmt,*args):
        if 'GET /api/job' not in (fmt%args): super().log_message(fmt,*args)
    def valid_host(self):
        host=self.headers.get('Host','').split(':')[0]
        return host in ('127.0.0.1','localhost')
    def send_data(self,code,data,kind='application/json',extra=None):
        if isinstance(data,(dict,list)):data=json.dumps(data,allow_nan=False).encode()
        elif isinstance(data,str):data=data.encode('utf-8')
        self.send_response(code);self.send_header('Content-Type',kind+'; charset=utf-8' if kind.startswith(('text/','application/json')) else kind)
        self.send_header('Content-Length',str(len(data)));self.send_header('X-Content-Type-Options','nosniff');self.send_header('Cache-Control','no-store')
        if extra:
            for k,v in extra.items():self.send_header(k,v)
        self.end_headers();self.wfile.write(data)
    def json_body(self):
        n=int(self.headers.get('Content-Length','0'))
        if n>1024*1024:raise InputError('JSON request is too large.')
        return json.loads(self.rfile.read(n))
    def do_GET(self):
        if not self.valid_host():return self.send_data(403,{'error':'Localhost access only.'})
        u=urllib.parse.urlparse(self.path);q=urllib.parse.parse_qs(u.query)
        try:
            if u.path=='/':
                page=(ROOT/'web/index.html').read_text(encoding='utf-8').replace('__TOKEN__',TOKEN)
                return self.send_data(200,page,'text/html')
            if u.path=='/api/session': return self.send_data(200,{'token':TOKEN,'version':__version__})
            if u.path=='/api/projects':
                return self.send_data(200,[{'id':d.name,'name':load_json(d/'project.json')['site_name']} for d in sorted((ROOT/'projects').iterdir()) if (d/'project.json').exists()])
            if u.path=='/api/job':
                with JOBS_LOCK:return self.send_data(200,dict(JOB))
            if u.path=='/api/state':
                p=project(q['project'][0]);return self.send_data(200,{'config':p.config,'readiness':p.readiness(),
                  'compatibility':load_json(p.path/'preflight.json') if (p.path/'preflight.json').exists() else None,
                  'prep_queue':load_json(p.path/'prep_queue.json') if (p.path/'prep_queue.json').exists() else None,
                  'production_queue':load_json(p.path/'production_queue.json') if (p.path/'production_queue.json').exists() else None})
            if u.path=='/api/results':
                p=project(q['project'][0]);rows=[]
                for f in sorted((p.path/'outputs').glob('*/*/summary.json')):
                    s=load_json(f);rows.append({'folder':str(f.parent.relative_to(p.path)).replace('\\','/'),'summary':s})
                return self.send_data(200,rows)
            if u.path=='/api/result':
                p=project(q['project'][0]);folder=p.file(q['folder'][0]);import pandas as pd
                d=pd.read_csv(folder/'daily_summary.csv');stride=max(1,len(d)//900);d=d.iloc[::stride]
                return self.send_data(200,{'summary':load_json(folder/'summary.json'),'daily':json.loads(d.to_json(orient='records'))})
            if u.path=='/file':
                p=project(q['project'][0]);path=p.file(q['path'][0])
                if not path.is_file():raise InputError('Output file not found.')
                extra={'Content-Disposition':'attachment; filename="'+path.name+'"'} if q.get('download') else None
                return self.send_data(200,path.read_bytes(),mimetypes.guess_type(path.name)[0] or 'application/octet-stream',extra)
            for prefix,base in [('/docs/',ROOT/'docs'),('/templates/',ROOT/'templates'),('/source_records/',ROOT/'source_records')]:
                if u.path.startswith(prefix):
                    f=safe_path(base,urllib.parse.unquote(u.path[len(prefix):]));return self.send_data(200,f.read_bytes(),mimetypes.guess_type(f.name)[0] or 'application/octet-stream')
            self.send_data(404,{'error':'Not found.'})
        except Exception as e:self.send_data(400,{'error':str(e)})
    def do_POST(self):
        if not self.valid_host() or self.headers.get('X-HMSL-Token')!=TOKEN:return self.send_data(403,{'error':'Local authorization failed. Reload the dashboard.'})
        origin=self.headers.get('Origin')
        if origin and urllib.parse.urlparse(origin).hostname not in ('127.0.0.1','localhost'):
            return self.send_data(403,{'error':'Cross-origin requests are not allowed.'})
        u=urllib.parse.urlparse(self.path);q=urllib.parse.parse_qs(u.query)
        try:
            if u.path=='/api/cancel':CANCEL.set();return self.send_data(200,{'ok':True})
            if u.path=='/api/upload':
                with JOBS_LOCK:
                    if JOB['status']=='running':raise InputError('Wait for the active job before changing inputs.')
                p=project(q['project'][0]);n=int(self.headers.get('Content-Length','0'))
                if n<=0 or n>600*1024**2:raise InputError('Upload must be 1 byte to 600 MB.')
                tmp=p.path/'.incoming_upload';remaining=n
                with open(tmp,'wb') as f:
                    while remaining:
                        b=self.rfile.read(min(1024*1024,remaining))
                        if not b:raise InputError('Incomplete file upload.')
                        f.write(b);remaining-=len(b)
                try:message=import_file(p,q['kind'][0],q['filename'][0],tmp,q)
                finally:tmp.unlink(missing_ok=True)
                return self.send_data(200,{'message':message})
            b=self.json_body()
            if u.path=='/api/job':launch_job(b);return self.send_data(200,{'ok':True})
            with JOBS_LOCK:
                if JOB['status']=='running':raise InputError('Wait for the active job before changing project inputs.')
            if u.path=='/api/create':
                clone=project(b['clone']).path if b.get('clone') else None
                return self.send_data(200,{'id':new_project(ROOT/'projects',b['name'],clone)})
            if u.path=='/api/settings':
                p=project(b['project']);cfg=p.config
                for group in ('geometry','surface','boundary','salt','initial','numerics','review'):
                    for key,val in b.get('values',{}).get(group,{}).items():
                        if key not in cfg[group]:raise InputError('Unknown configuration field: '+key)
                        if isinstance(val,str) and key not in ('provenance','diffusion_convention','solubility_file','liquid_face_scheme'):raise InputError('Expected a number/boolean for '+key)
                        cfg[group][key]=val
                p.config=cfg;errs=p.validate_settings()
                if errs:raise InputError('\n'.join(errs))
                p.save(cfg);return self.send_data(200,{'message':'Settings saved. Any earlier spin-up or QA with different inputs is invalidated automatically.'})
            self.send_data(404,{'error':'Unknown endpoint.'})
        except Exception as e:self.send_data(400,{'error':str(e)})

def main():
    a=argparse.ArgumentParser();a.add_argument('--port',type=int,default=8765);a.add_argument('--no-browser',action='store_true');args=a.parse_args()
    try:server=http.server.ThreadingHTTPServer(('127.0.0.1',args.port),Handler)
    except OSError:
        server=http.server.ThreadingHTTPServer(('127.0.0.1',0),Handler)
    url=f'http://127.0.0.1:{server.server_port}'
    print('\nHeritage Moisture & Salt Lab '+__version__+'\n'+url+'\nKeep this window open. Ctrl+C stops the application.\n',flush=True)
    if not args.no_browser:threading.Timer(.8,lambda:webbrowser.open(url)).start()
    try:server.serve_forever()
    except KeyboardInterrupt:CANCEL.set();server.server_close()

if __name__=='__main__':main()
