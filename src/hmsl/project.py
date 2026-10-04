from __future__ import annotations
import copy, io, json, shutil, zipfile, hashlib, urllib.request
from pathlib import Path
import numpy as np
from .common import InputError, atomic_json, digest, load_json, safe_name, safe_path, sha256, now
from .material import read_material
from .climate import read_climate
from .solver import Solubility

ARCHIVE_URL='https://zenodo.org/records/5656966/files/Material_files_for_DELPHIN_software_HAM.zip?download=1'
ARCHIVE_MD5='609b7de53b838c0868814e00ff288ecd'
SOURCE_MEMBER='18_15_IBEET_3.m6'

class Project:
    def __init__(self,path):
        self.path=Path(path).resolve(); self.config=load_json(self.path/'project.json')
        self.config.setdefault('numerics',{}).setdefault('liquid_face_scheme','harmonic_legacy')
        # In-memory schema migration only. User approves and saves explicitly.
        self.config.setdefault('review',{}).setdefault('wet_end_closure_approved',bool(self.config.get('is_demo',False)))
    def file(self,rel): return safe_path(self.path,rel)
    def run(self,rid):
        matches=[r for r in self.config['runs'] if r['id']==rid]
        if len(matches)!=1: raise InputError('Run ID not found or not unique.')
        return matches[0]
    def climate(self,rid):
        r=self.run(rid)
        return read_climate(self.file(r['file']),r['calendar'],r.get('column_map'),r.get('expected_sha256'))
    def material(self):
        m=read_material(self.file(self.config['material_file']),self.config.get('retention_branch','de'))
        identity=self.config.get('require_material_identity')
        if identity and identity not in m.name: raise InputError(f'This project requires material {identity}; clone the project to change the material identity.')
        if not self.config.get('is_demo') and m.is_demo: raise InputError('Synthetic material cannot run a real-site project. Use the separate Demo project.')
        if identity=='18_15':
            targets={'density_kg_m3':1831.3,'specific_heat_J_kgK':873.9,'porosity':.38,'lambda_W_mK':.626}
            for k,val in targets.items():
                if not np.isclose(float(m.obj[k]),val,rtol=.005,atol=.0001):
                    raise InputError(f'The downloaded {k} does not match the selected workbook material within rounding tolerance. Review source versions; no scalar was overwritten.')
        return m
    def phase(self):
        s=self.config['salt'].get('solubility_file')
        return Solubility(self.file(s)) if s else None
    def physics_hash(self):
        c=copy.deepcopy(self.config)
        c.pop('runs',None)
        c.pop('site_name',None)
        # Actual implementation bytes, including the numerical face operator.
        # Scientific revisions must not inherit a previous numerical certificate.
        c['implementation_files']={n:sha256(Path(__file__).parent/n)
            for n in ('solver.py','material.py','face_integral.py','climate.py','workflow.py','project.py','common.py')}
        for key in ('material_file',):
            p=self.file(c[key]); c[key+'_sha256']=sha256(p) if p.is_file() else 'MISSING'
        s=c['salt'].get('solubility_file')
        if s: c['salt']['solubility_sha256']=sha256(self.file(s)) if self.file(s).is_file() else 'MISSING'
        return digest(c)
    def save(self,c):
        if self.file('freeze.json').is_file(): raise InputError('This project is frozen. Clone it to change inputs; do not edit the production base.')
        atomic_json(self.path/'project.json',c); self.config=c
    def validate_settings(self):
        c=self.config; errors=[]
        fields=[('geometry','height_m',.01,100),('geometry','thickness_m',.01,10),('geometry','cells',4,400),('geometry','exposed_faces',0,2),
          ('surface','emissivity',.01,1),('surface','absorptivity',0,1),('surface','solar_projection_factor',0,2),
          ('surface','sky_view_factor',0,1),('surface','rain_catch_fraction',0,1),
          ('boundary','capillary_supply_mm_day',0,1000),('boundary','sulfate_source_mg_L',0,100000),
          ('salt','D_sulfate_m2_s',1e-15,1e-7),('salt','unsaturated_exponent',0,10),('salt','initial_sulfate_mol_m3_bulk',0,1e5),
          ('initial','RH_pct',.01,99.999),('initial','temperature_C',0,80),
          ('numerics','max_step_s',1,3600),('numerics','rtol',1e-10,1e-3),('numerics','chunk_hours',1,168),
          ('numerics','spinup_max_cycles',1,100),('numerics','spinup_theta_tolerance',1e-8,.01),
          ('numerics','spinup_temperature_tolerance_C',1e-5,1),('numerics','qa_relative_tolerance',1e-5,.1)]
        for a,b,low,high in fields:
            try:
                val=float(c[a][b])
                if not np.isfinite(val) or not low<=val<=high: raise ValueError()
            except (TypeError,ValueError,KeyError): errors.append(f'{a}.{b} must be between {low} and {high}.')
        for group,key in [('geometry','cells'),('geometry','exposed_faces'),('numerics','chunk_hours'),('numerics','spinup_max_cycles')]:
            try:
                if float(c[group][key])!=int(c[group][key]): errors.append(group+'.'+key+' must be an integer.')
            except (ValueError,TypeError,KeyError,OverflowError): pass
        if c.get('numerics',{}).get('liquid_face_scheme','harmonic_legacy') not in ('harmonic_legacy','source_integral'): errors.append('Choose a supported liquid face discretisation.')
        if c['salt']['diffusion_convention'] not in ('pore','bulk'): errors.append('Choose pore or bulk diffusion convention.')
        if c.get('require_material_identity')=='18_15' and c['salt']['D_sulfate_m2_s']!=2.91e-10:
            errors.append('The locked base sulfate coefficient is 2.91e-10 m2/s. Clone a sensitivity project to change it.')
        return errors
    def readiness(self, thorough=False):
        out={'site':self.config['site_name'],'is_demo':self.config.get('is_demo',False),'runs':[],
             'errors':self.validate_settings(),'warnings':[], 'material_ok':False,'review_ok':False}
        ph=self.physics_hash()
        for r in self.config['runs']:
            f=self.file(r['file']); item=dict(r); item['present']=f.is_file(); item['status']='Present' if f.is_file() else 'MISSING'
            if thorough and f.is_file():
                try:
                    climate=self.climate(r['id']); item['hours']=climate.n; item['status']='Validated'; item['warnings']=climate.warnings
                except Exception as e: item['status']='ERROR'; item['error']=str(e)
            qpath=self.path/'qa_cases'/(safe_name(r['id'])+'.json')
            q=load_json(qpath) if qpath.is_file() else {}
            qvalid=False
            if q:
                try:
                    qvalid=bool(q.get('passed')) and q.get('physics_hash')==ph and q.get('run_id')==r['id'] and f.is_file() and q.get('climate_sha256')==sha256(f)
                except Exception:qvalid=False
            item['qa_valid']=qvalid
            item['qa_status']='PASS' if qvalid else ('STALE' if q.get('passed') else ('FAIL' if q else 'NOT RUN'))
            item['qa_record']='qa_cases/'+safe_name(r['id'])+'.json' if q else None
            out['runs'].append(item)
        try:
            m=self.material()
            if self.config['numerics']['liquid_face_scheme']=='source_integral':
                from .face_integral import LiquidFaceIntegral
                LiquidFaceIntegral(m)
            out['material_ok']=True; out['material_name']=m.name
            out['warnings']+=m.warnings
            gap=(m.sat-m.hi)/m.sat
            if gap>0.001:
                out['errors'].append('Source functions end more than 0.1% below effective saturation. Complete source coverage is required; a truncated material table is not treated as saturation.')
            out['wet_end']={'supported_theta':m.hi,'effective_saturation':m.sat,'relative_gap':gap,
                'closure':'Positive net water influx at the source-supported wet endpoint leaves exposed cells as liquid seepage, carrying the computed dissolved-equivalent salt concentration. No saturated pressure or 2-D runoff model.'}
            if self.config['salt'].get('solubility_file'): self.phase()
        except Exception as e: out['errors'].append(str(e))
        out['review_ok']=all(self.config.get('review',{}).get(k,False) for k in ('geometry_boundary_approved','diffusion_convention_approved','reduced_scope_approved'))
        if not self.config.get('is_demo') and not self.config['review'].get('wet_end_closure_approved',False):
            out['review_ok']=False
            out['errors'].append('Update 1.1: review the exposed-face wet-end drainage assumption in Model settings, then Save reviewed settings. Original source tables are unchanged.')
        if not out['review_ok']: out['errors'].append('Review geometry/exposure, parametric lower supply, sulfate convention and the reduced model scope in Settings.')
        if not self.config['salt'].get('solubility_file'):
            out['warnings'].append('Transport-only mode: no sulfate equilibrium solubility table has been supplied. No crystal/damage results will be generated.')
        for name in ('qa','spinup','freeze'):
            f=self.file(name+'.json'); val=load_json(f) if f.is_file() else {}
            out[name]=val
            out[name+'_valid']=val.get('physics_hash')==ph and bool(val.get('passed',False))
        if out['qa_valid']:
            try:
                run=self.run(out['qa']['run_id'])
                out['qa_valid']=sha256(self.file(run['file']))==out['qa'].get('climate_sha256')
            except Exception: out['qa_valid']=False
        if out['spinup_valid']:
            state=self.file('spinup.npz')
            out['spinup_valid']=state.is_file() and sha256(state)==out['spinup'].get('initial_state_sha256')
            try:
                run=self.run(out['spinup']['baseline_run'])
                out['spinup_valid']=out['spinup_valid'] and sha256(self.file(run['file']))==out['spinup'].get('baseline_climate_sha256')
            except Exception: out['spinup_valid']=False
        if out['freeze_valid']:
            out['freeze_valid']=out['freeze'].get('run_plan_hash')==digest(self.config['runs'])
            if out['freeze_valid'] and out['freeze'].get('initial_state_sha256'):
                state=self.file('spinup.npz')
                out['freeze_valid']=state.is_file() and sha256(state)==out['freeze'].get('initial_state_sha256')
            ref=out['freeze'].get('qa_reference_run')
            refhash=out['freeze'].get('qa_reference_sha256')
            if out['freeze_valid'] and ref and refhash:
                qpath=self.path/'qa_cases'/(safe_name(ref)+'.json')
                out['freeze_valid']=qpath.is_file() and sha256(qpath)==refhash
        prep=self.path/'prep_queue.json'; out['prep_queue']=load_json(prep) if prep.is_file() else {}
        out['trial_ready']=out['material_ok'] and out['review_ok'] and not out['errors']
        # A later failed QA on another climate case must not invalidate a valid
        # reproducibility freeze. Production eligibility is checked per case.
        out['frozen_ready']=out['trial_ready'] and out['spinup_valid'] and out['freeze_valid']
        return out

def new_project(root,name,clone_from=None):
    sid=safe_name(name).lower(); dest=Path(root)/sid
    if dest.exists(): raise InputError('A project with this name already exists.')
    if clone_from:
        old=Project(clone_from); c=copy.deepcopy(old.config)
    else:
        old=Project(Path(root)/'mohenjo_daro'); c=copy.deepcopy(old.config)
        c['runs']=[]; c['material_file']='materials/YOUR_MATERIAL.m6'; c['is_demo']=False
        c['require_material_identity']=None; c['material_source']='Supply a citation for your selected material'
        c['site_latitude']=None; c['site_longitude']=None
        c['salt']['D_sulfate_m2_s']=None
        c['salt']['provenance']='Supply material/salt-specific evidence and explain the concentration and flux convention.'
        c['salt']['solubility_file']=None
        c['boundary']['capillary_supply_mm_day']=None
        c['boundary']['sulfate_source_mg_L']=None
        c['boundary']['provenance']='Enter local evidence or explicitly justified scenario assumptions.'
    c['site_id']=sid; c['site_name']=str(name)
    c['input_selection']='User site project; Mohenjo proxy properties are not automatically valid for another site'
    c['review']={k:False for k in c['review']}
    dest.mkdir(parents=True); (dest/'climate').mkdir(); (dest/'materials').mkdir()
    if clone_from:
        for sub in ('climate','materials'):
            for f in (old.path/sub).glob('*'):
                if f.is_file(): shutil.copy2(f,dest/sub/f.name)
        c['require_material_identity']=None # explicit sensitivity clone, no longer R00-R20 base
        c['runs']=[dict(r,id='S_'+r['id']) for r in c['runs']]
    atomic_json(dest/'project.json',c)
    return sid

def download_selected_material(project):
    if project.file('freeze.json').exists(): raise InputError('Frozen project; source files cannot be replaced.')
    try:
        req=urllib.request.Request(ARCHIVE_URL,headers={'User-Agent':'HeritageMoistureSaltLab/1.0 (research source download)'})
        with urllib.request.urlopen(req,timeout=45) as r: data=r.read(10*1024*1024+1)
    except Exception as e:
        raise InputError('Source download could not be completed. Open docs/SOURCE_DATA_LINKS.html, download Material_files_for_DELPHIN_software_HAM.zip, and import that ZIP or its 18_15_IBEET_3.m6 file. Network error: '+str(e)) from e
    if len(data)>10*1024*1024 or hashlib.md5(data).hexdigest()!=ARCHIVE_MD5:
        raise InputError('Archive checksum/size differs from Zenodo version 5656966. No material was installed; review source version.')
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        matches=[n for n in z.namelist() if Path(n).name==SOURCE_MEMBER]
        if len(matches)!=1: raise InputError('Selected material is absent or duplicated in archive.')
        raw=z.read(matches[0])
    p=project.path/'materials'/SOURCE_MEMBER
    p.parent.mkdir(exist_ok=True); p.write_bytes(raw)
    atomic_json(project.path/'materials/source_receipt.json',{'url':ARCHIVE_URL,'archive_md5':ARCHIVE_MD5,'sha256':sha256(p),'retrieved_utc':now(),'license':'CC BY 4.0; Freimanis et al. historic Latvian brick dataset'})
    m=project.material()
    return {'material':m.name,'sha256':sha256(p),'message':'Source bytes verified and constitutive tables parsed. Review assumptions and run numerical QA before freezing.'}
