"""Durable sidecar records, verified backups and read-only checkpoint auditing."""
from __future__ import annotations
import csv, hashlib, json, os, shutil, socket, time, uuid
from datetime import datetime, timezone
from pathlib import Path

class SafetyError(RuntimeError): pass

def now(): return datetime.now(timezone.utc).isoformat()
def stamp(): return datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')+'_'+uuid.uuid4().hex[:6]
def read_json(path):
    with open(path,encoding='utf-8') as f:return json.load(f)
def hash_file(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def atomic(path,obj):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    temp=p.with_name(p.name+'.'+uuid.uuid4().hex+'.tmp')
    try:
        with open(temp,'w',encoding='utf-8',newline='\n') as f:
            json.dump(obj,f,indent=2,allow_nan=False);f.flush();os.fsync(f.fileno())
        os.replace(temp,p)
    finally:temp.unlink(missing_ok=True)
def safe_path(root,relative):
    root=Path(root).resolve();p=(root/str(relative)).resolve()
    if not p.is_relative_to(root):raise SafetyError('A path would leave the selected folder.')
    return p

def verify_core(root,manifest=None):
    manifest=manifest or read_json(Path(__file__).with_name('expected_core.json'))
    bad=[n for n,h in manifest.items() if not (Path(root)/n).is_file() or hash_file(Path(root)/n)!=h]
    if bad:raise SafetyError('This add-on supports the verified 1.2.0/1.2.1 scientific core only. Different/missing files: '+', '.join(bad))
    return dict(manifest)

def port_open(port=8765):
    try:
        with socket.create_connection(('127.0.0.1',int(port)),timeout=.4):return True
    except OSError:return False

class FileLock:
    """Kernel-held cross-platform lock. Stale files are not stale locks."""
    def __init__(self,path):self.path=Path(path);self.f=None
    def __enter__(self):
        self.path.parent.mkdir(parents=True,exist_ok=True)
        self.f=open(self.path,'a+b')
        self.f.seek(0,2)
        if self.f.tell()==0:self.f.write(b'0');self.f.flush()
        self.f.seek(0)
        try:
            if os.name=='nt':
                import msvcrt
                msvcrt.locking(self.f.fileno(),msvcrt.LK_NBLCK,1)
            else:
                import fcntl
                fcntl.flock(self.f.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        except OSError as exc:
            self.f.close();self.f=None
            raise SafetyError('Another process is already using this project/case. No second writer was started.') from exc
        return self
    def __exit__(self,*args):
        if self.f:
            if os.name=='nt':
                import msvcrt
                self.f.seek(0);msvcrt.locking(self.f.fileno(),msvcrt.LK_UNLCK,1)
            else:
                import fcntl
                fcntl.flock(self.f.fileno(),fcntl.LOCK_UN)
            self.f.close();self.f=None

def _identity(p):
    s=p.stat();return (s.st_size,s.st_mtime_ns)

def verified_copy(src,dst):
    src=Path(src);dst=Path(dst)
    if src.is_symlink():raise SafetyError('Symbolic links are not followed in a safety backup: '+str(src))
    before=_identity(src);dst.parent.mkdir(parents=True,exist_ok=True)
    h=hashlib.sha256()
    with open(src,'rb') as a,open(dst,'xb') as b:
        for block in iter(lambda:a.read(1024*1024),b''):
            h.update(block);b.write(block)
        b.flush();os.fsync(b.fileno())
    if _identity(src)!=before or hash_file(src)!=h.hexdigest():raise SafetyError('Source changed during backup. Stop the running model and retry. Incomplete backup retained: '+str(dst.parent))
    if hash_file(dst)!=h.hexdigest():raise SafetyError('Backup checksum verification failed: '+str(dst))
    shutil.copystat(src,dst)
    return {'bytes':before[0],'sha256':h.hexdigest()}

def backup_projects(root,progress=print):
    root=Path(root).resolve();source=root/'projects'
    if not source.is_dir():raise SafetyError('The selected application has no projects folder.')
    files=sorted(p for p in source.rglob('*') if p.is_file())
    if any(p.is_symlink() for p in source.rglob('*')):raise SafetyError('A project contains symbolic links. Create an ordinary independent project copy first.')
    size=sum(p.stat().st_size for p in files)
    reserve=512*1024**2
    if shutil.disk_usage(root).free<size+reserve:raise SafetyError(f'Not enough free disk space for a verified full project backup ({size/1024**3:.2f} GB plus 0.5 GB reserve). No application files changed.')
    dest=root/'performance_backups'/stamp();dest.mkdir(parents=True)
    records={};progress(f'Backing up {len(files):,} files ({size/1024**3:.2f} GB). Please wait.')
    for i,p in enumerate(files):
        rel=p.relative_to(root).as_posix();records[rel]=verified_copy(p,dest/rel)
        if i%250==0:progress(f'Verified backup {i+1:,}/{len(files):,} files')
    # Recheck inventory and metadata after copying to detect concurrent writers.
    current=sorted(p.relative_to(root).as_posix() for p in source.rglob('*') if p.is_file())
    if current!=list(records):raise SafetyError('Project file list changed during backup. Original files are untouched; stop the model and retry.')
    for rel,record in records.items():
        if (root/rel).stat().st_size!=record['bytes'] or hash_file(root/rel)!=record['sha256']:
            raise SafetyError('A project file changed before backup verification completed: '+rel)
    atomic(dest/'BACKUP_COMPLETE.json',{'complete':True,'created_utc':now(),'file_count':len(records),'total_bytes':size,'files':records})
    return dest

def checkpoint_pair(folder):
    """Read a matching snapshot; retries handle a live writer advancing its state."""
    folder=Path(folder)
    for _ in range(3):
        try:
            cp=read_json(folder/'checkpoint.json')
            if Path(cp['state_file']).name!=cp['state_file']:raise SafetyError('Invalid state filename.')
            state=(folder/cp['state_file']).read_bytes()
            if hashlib.sha256(state).hexdigest()!=cp['state_sha256']:raise SafetyError('Checkpoint state checksum mismatch.')
            if read_json(folder/'checkpoint.json')!=cp:time.sleep(.03);continue
            return cp,state
        except FileNotFoundError:
            time.sleep(.03)
    raise SafetyError('No stable matching checkpoint/state pair could be read.')

def audit_chunks(folder,completed):
    folder=Path(folder);expected=1;count=0
    for p in sorted((folder/'chunks').glob('hours_*.csv')):
        endpoint=int(p.stem.split('_')[-1])
        if endpoint>completed:continue
        last=None
        with open(p,newline='',encoding='utf-8') as f:
            reader=csv.DictReader(f)
            if 'elapsed_hours' not in (reader.fieldnames or []):raise SafetyError('Chunk lacks elapsed_hours: '+p.name)
            for row in reader:
                hour=int(row['elapsed_hours'])
                if hour!=expected:raise SafetyError(f'Chunk continuity failed at {p.name}: expected hour {expected}, found {hour}. Nothing was overwritten.')
                expected+=1;count+=1;last=hour
        if last!=endpoint:raise SafetyError('Chunk endpoint/content mismatch: '+p.name)
        profile=folder/'chunks'/('profile_'+f'{endpoint:09d}'+'.csv')
        if not profile.is_file():raise SafetyError('Missing saved profile chunk: '+profile.name)
    if count!=int(completed):raise SafetyError(f'Only {count:,} committed output rows exist for checkpoint {completed:,}. Restore matching chunks before resuming.')
    return count

def inspect_output(folder,deep=False):
    folder=Path(folder);prov=read_json(folder/'provenance.json')
    result={'run_id':prov.get('run_id'),'folder':str(folder),'physics_hash':prov.get('physics_hash'),'climate_sha256':prov.get('climate_sha256'),'mode':prov.get('mode'),'saved_hours':0,'status':'NOT_STARTED'}
    if (folder/'checkpoint.json').is_file():
        cp,state=checkpoint_pair(folder);result.update(saved_hours=int(cp['completed_hours']),key=cp['key'],saved_utc=cp.get('saved_utc'),status='PARTIAL')
        if deep:audit_chunks(folder,result['saved_hours'])
    if (folder/'completed.json').is_file():
        done=read_json(folder/'completed.json')
        if not done.get('passed') or done.get('key')!=result.get('key') or done.get('hours')!=result['saved_hours']:
            raise SafetyError('Completion marker and checkpoint do not match: '+str(folder))
        for name in ['summary.json','hourly_results.csv.gz','profiles.csv.gz','daily_summary.csv','report.html']:
            if not (folder/name).is_file():raise SafetyError('Completion marker exists but final output is missing: '+name)
        if (folder/'phase_extension_v140.json').is_file():
            phasefile=folder/'phase_extension_summary.json'
            if not phasefile.is_file():raise SafetyError('Phase-extension completion lacks its phase summary. No successful status assigned.')
            ph=read_json(phasefile)
            if ph.get('phase_model')!='equilibrium-freeze-thaw-1.0':raise SafetyError('Phase-extension summary has an unexpected model identity.')
        elif (folder/'cold_scope_policy.json').is_file():
            scopefile=folder/'cold_scope_summary.json'
            if not scopefile.is_file():raise SafetyError('New-scope completion lacks its cold-scope audit. No successful status assigned.')
            sc=read_json(scopefile)
            if sc.get('passed_under_screen') is not True or sc.get('hours_total')!=done['hours']:
                raise SafetyError('Cold-scope audit does not match the completed run.')
        result.update(status='COMPLETED',total_hours=done['hours'])
    return result

def mirror_checkpoint(folder,recovery):
    folder=Path(folder);recovery=Path(recovery);cp,state=checkpoint_pair(folder)
    dest=recovery/f"hour_{int(cp['completed_hours']):09d}"
    if (dest/'RECOVERY_COMPLETE.json').is_file():return str(dest)
    dest.mkdir(parents=True,exist_ok=True)
    statepath=dest/cp['state_file']
    with open(statepath,'wb') as f:f.write(state);f.flush();os.fsync(f.fileno())
    if hash_file(statepath)!=cp['state_sha256']:raise SafetyError('Recovery state verification failed.')
    atomic(dest/'checkpoint.json',cp)
    for prefix in ('hours','profile'):
        name=f'{prefix}_{cp["completed_hours"]:09d}.csv';src=folder/'chunks'/name
        if src.is_file():
            target=dest/'chunks'/name
            if not target.exists():verified_copy(src,target)
    for name in ('provenance.json','resolved_project.json'):
        if (folder/name).is_file() and not (dest/name).exists():verified_copy(folder/name,dest/name)
    atomic(dest/'RECOVERY_COMPLETE.json',{'complete':True,'saved_hours':cp['completed_hours'],'state_sha256':cp['state_sha256'],'note':'Latest state and last chunk only. Earlier chunks remain in the original outputs folder/full safety backup.','created_utc':now()})
    complete=sorted(p for p in recovery.glob('hour_*') if (p/'RECOVERY_COMPLETE.json').is_file())
    for old in complete[:-2]:shutil.rmtree(old)
    return str(dest)
