"""Read source tables without fitting or mixing material identities.

DELPHIN ASCII .m6 is a container, not an executable model. This adapter supports
isotropic/U-direction tabulated functions only. Unsupported constitutive data
raise an error; they are not silently discarded. Retention hysteresis is not
modelled: the branch selection is explicit in the project.
"""
from __future__ import annotations
import re, json
from pathlib import Path
import numpy as np
from .common import InputError, sha256
RW=1000.; RV=461.5

def interp(x,xy): return np.interp(x,np.asarray(xy[0],float),np.asarray(xy[1],float))

def table(xy,name):
    a=np.array(xy,dtype=float)
    if a.ndim!=2 or a.shape[0]!=2 or a.shape[1]<2 or not np.isfinite(a).all():
        raise InputError(f'Invalid {name} table: need two finite equal-length arrays.')
    if not np.all(np.diff(a[0])>0): raise InputError(f'{name} x-values must be strictly increasing.')
    return a

def clean_id(s):
    s=re.sub(r'\s+','',s).replace('Ol','Theta_l').replace('My(Phi)','mew(RH)')
    return s

def parse_m6(path, branch='de'):
    raw=Path(path).read_bytes()
    if b'\x00' in raw: raise InputError('Binary .m6 is not supported. Export a DELPHIN ASCII material file or use the tabulated JSON template.')
    try: text=raw.decode('utf-8-sig')
    except UnicodeDecodeError: text=raw.decode('cp1252')
    lines=[ln.split('#',1)[0].strip() for ln in text.splitlines()]
    functions={}; scalars={}; section=''; i=0; anisotropic=False
    while i<len(lines):
        ln=lines[i]; i+=1
        if not ln: continue
        if ln.startswith('['):
            section=ln.strip('[] ').replace(' ','')
            if section.endswith(('_V','_W')): anisotropic=True
            continue
        if ln.startswith('FUNCTION') and '=' in ln:
            key=clean_id(ln.split('=',1)[1]); vectors=[]
            while len(vectors)<2 and i<len(lines):
                v=lines[i]; i+=1
                if not v: continue
                try: vals=[float(t.replace('D','E').replace('d','e')) for t in v.split()]
                except ValueError as e: raise InputError(f'Malformed {key} table. Each vector must occupy one line, as in the ASCII specification.') from e
                vectors.append(vals)
            if key in functions: raise InputError(f'Duplicate function {key}; choose one explicit direction/representation.')
            functions[key]=table(vectors,key).tolist()
        elif '=' in ln:
            k,v=ln.split('=',1); k=k.strip()
            if k.endswith(('_V','_W')): anisotropic=True
            k=k.removesuffix('_U')
            match=re.match(r'^\s*([-+0-9.eEdD]+)(?:\s|$)',v)
            if match:
                try: scalars[k]=float(match.group(1).replace('D','E').replace('d','e'))
                except ValueError: pass
    if anisotropic: raise InputError('Anisotropic material file detected. This reduced solver requires an explicitly selected isotropic/U-only export.')
    aliases={'RHO':'density_kg_m3','CE':'specific_heat_J_kgK','THETA_POR':'porosity','OPOR':'porosity',
             'THETA_EFF':'theta_effective','OEFF':'theta_effective','THETA_CAP':'theta_capillary',
             'OCAP':'theta_capillary','LAMBDA':'lambda_W_mK','MEW':'mu','AW':'Aw_kg_m2_s05'}
    out={aliases[k]:v for k,v in scalars.items() if k in aliases}
    for k in ['density_kg_m3','specific_heat_J_kgK','porosity','theta_effective','lambda_W_mK']:
        if k not in out: raise InputError(f'.m6 lacks required scalar {k}. The workbook cannot supply a missing saturation by inference.')
    out.update(name=Path(path).stem,source_sha256=sha256(path),is_demo=False,functions=functions,
               retention_branch=branch,source_file=Path(path).name)
    return out

class Material:
    def __init__(self, obj):
        self.obj=obj; self.f=obj['functions']; self.name=obj['name']; self.is_demo=bool(obj.get('is_demo',False))
        for key in ('density_kg_m3','specific_heat_J_kgK','porosity','theta_effective','lambda_W_mK'):
            val=obj.get(key)
            if val is None or not np.isfinite(float(val)) or float(val)<=0: raise InputError(f'Material needs a positive {key}.')
        self.rho=float(obj['density_kg_m3']); self.cp=float(obj['specific_heat_J_kgK'])
        self.por=float(obj['porosity']); self.sat=float(obj['theta_effective']); self.lam=float(obj['lambda_W_mK'])
        if not (0<self.sat<=self.por<1): raise InputError('Require 0 < effective saturation <= porosity < 1.')
        for k,xy in self.f.items(): self.f[k]=table(xy,k).tolist()
        branch=obj.get('retention_branch','de')
        candidates=[f'Theta_l(pC)_{branch}',f'pC(Theta_l)_{branch}',f'Theta_l(RH)_{branch}']
        if branch=='de': candidates += ['Ol(pC)','Ol(RH)','Theta_l(pC)','Theta_l(RH)','pC(Theta_l)']
        candidates=[k for k in candidates if k in self.f]
        if not candidates:
            raise InputError(f'No retention representation for branch {branch} was found.')
        # DELPHIN material files may legitimately contain the same retention law in
        # both forward and inverse tabulations, e.g. Theta_l(pC)_de and
        # pC(Theta_l)_de.  These are representations of one constitutive law, not
        # two competing laws.  Accept the pair only after checking that they are
        # mutually consistent, then use Theta_l(pC) as the canonical form.
        direct=f'Theta_l(pC)_{branch}'; inverse=f'pC(Theta_l)_{branch}'
        if direct in candidates and inverse in candidates:
            fd=np.asarray(self.f[direct],dtype=float); fi=np.asarray(self.f[inverse],dtype=float)
            if not np.all(np.diff(fd[0])>0) or not np.all(np.diff(fd[1])<0):
                raise InputError(f'{direct} must have increasing pC and decreasing moisture.')
            if not np.all(np.diff(fi[0])>0) or not np.all(np.diff(fi[1])<0):
                raise InputError(f'{inverse} must have increasing moisture and decreasing pC.')
            lo=max(float(fd[1].min()),float(fi[0].min())); hi=min(float(fd[1].max()),float(fi[0].max()))
            if hi<=lo:
                raise InputError(f'{direct} and {inverse} have no overlapping moisture range.')
            sample=np.linspace(lo,hi,25)
            pc_from_direct=np.interp(sample,fd[1][::-1],fd[0][::-1])
            pc_from_inverse=np.interp(sample,fi[0],fi[1])
            span=max(float(np.ptp(np.r_[fd[0],fi[1]])),1.0)
            tol=max(0.10,0.03*span)
            if np.max(np.abs(pc_from_direct-pc_from_inverse))>tol:
                raise InputError(f'{direct} and {inverse} are both present but are not mutually consistent. Refusing to choose silently.')
            candidates=[k for k in candidates if k!=inverse]
        if len(candidates)!=1:
            raise InputError(f'Multiple independent retention representations for branch {branch} are present. Found {candidates}. Select one source branch explicitly.')
        self.store=candidates[0]; xy=np.asarray(self.f[self.store])
        th=xy[0] if self.store.startswith('pC(') else xy[1]
        if self.store.startswith('Theta_l(RH)'):
            if xy[0,0]<0 or xy[0,-1]>1: raise InputError('Sorption RH abscissa must use 0-1, not percent.')
            direction=1
        elif self.store.startswith('pC('):
            direction=1
            if not np.all(np.diff(xy[1])<0): raise InputError('Capillary suction pC must decrease strictly with increasing moisture.')
        else: direction=-1
        if not np.all(np.diff(th)*direction>0):
            raise InputError('Retention moisture values must be strictly monotonic for reversible interpolation; export an invertible source table rather than silently dropping plateaus.')
        self.lo=max(float(min(th)),1e-10); self.hi=min(float(max(th)),self.sat)
        if self.hi-self.lo<0.001: raise InputError('Moisture storage range is too narrow for this transport solver.')
        lk=[k for k in ('lgKl(Theta_l)','lgDl(Theta_l)') if k in self.f]
        if len(lk)!=1: raise InputError('Exactly one full liquid function lgKl(Theta_l) or lgDl(Theta_l) is required; KLEFF/Aw alone is not a curve.')
        self.liquid=lk[0]
        if 'lgKv(Theta_l)' not in self.f and 'mew(RH)' not in self.f and not obj.get('mu'):
            raise InputError('Material needs lgKv(Theta_l), mew(RH) or a scalar mu.')
        allowed={self.store,'lgKl(Theta_l)','lgDl(Theta_l)','lgKv(Theta_l)','mew(RH)','lambda(Theta_l)','lambda(T)'}
        known_branches={f'{name}_{b}' for name in ('Theta_l(pC)','pC(Theta_l)','Theta_l(RH)') for b in ('ad','de')}
        unsupported=set(self.f)-allowed-known_branches
        if unsupported: raise InputError('Unsupported material functions: '+', '.join(sorted(unsupported))+'. No functions were silently discarded.')
        for k in [self.liquid]+[x for x in ('lgKv(Theta_l)','lambda(Theta_l)') if x in self.f]:
            x=np.asarray(self.f[k][0])
            self.lo=max(self.lo,float(x[0])); self.hi=min(self.hi,float(x[-1]))
        if self.lo>=self.hi: raise InputError('Material storage/transport tables have no common moisture range.')
        self.warnings=['Single retention branch; no hysteresis or salt-dependent sorption.',
          'Moisture/thermal functions are source functions. This solver is not DELPHIN and has not been cross-validated against it.']
        # Cache numeric arrays separately; original source object and bytes stay intact.
        self.f={k:np.asarray(v,dtype=float) for k,v in self.f.items()}
        if self.hi<self.sat:
            self.warnings.append(f'Common wet-end table coverage is {self.hi:g}; source effective saturation is {self.sat:g}. This distinction is retained and reported.')
        if 'lgDl(Theta_l)' in self.f: self.warnings.append('Kl is derived from Dl and the source retention derivative for gravity; pressure-driven comparisons require QA.')
    def theta_at_rh(self,rh,T=20.):
        x=np.clip(np.asarray(rh,float),1e-12,1-1e-12)
        if 'RH)' in self.store: return interp(x,self.f[self.store])
        pc=np.log10(-RW*RV*(T+273.15)*np.log(x))
        if self.store.startswith('pC('):
            xy=np.asarray(self.f[self.store]); return np.interp(pc,xy[1][::-1],xy[0][::-1])
        return interp(pc,self.f[self.store])
    def pressure(self,theta,T):
        th=np.clip(theta,self.lo,self.hi)
        xy=np.asarray(self.f[self.store])
        if self.store.startswith('pC('): return -10.**interp(th,self.f[self.store])
        if 'pC)' in self.store: return -10.**np.interp(th,xy[1][::-1],xy[0][::-1])
        rh=np.clip(np.interp(th,xy[1],xy[0]),1e-12,1-1e-12)
        return RW*RV*(T+273.15)*np.log(rh)
    def fields(self,theta,T):
        th=np.clip(theta,self.lo,self.hi); p=self.pressure(th,T)
        rh=np.exp(np.maximum(p/(RW*RV*(T+273.15)),-700))
        if self.liquid.startswith('lgKl'): kl=10.**interp(th,self.f[self.liquid])
        else:
            eps=1e-7
            a=np.clip(th-eps,self.lo,self.hi); b=np.clip(th+eps,self.lo,self.hi)
            dp=(self.pressure(b,T)-self.pressure(a,T))/np.maximum(b-a,1e-15)
            kl=RW*10.**interp(th,self.f[self.liquid])/np.maximum(dp,1e-12)
        if 'lgKv(Theta_l)' in self.f: kv=10.**interp(th,self.f['lgKv(Theta_l)'])
        else:
            mu=interp(rh,self.f['mew(RH)']) if 'mew(RH)' in self.f else float(self.obj['mu'])
            # Declared reference air permeability; not a fitted material parameter.
            kv=np.full_like(th,2.0e-10)/mu
        lt='lambda(T)' in self.f; lw='lambda(Theta_l)' in self.f
        if lt:
            tx=self.f['lambda(T)'][0]
            if np.any(T<tx[0]) or np.any(T>tx[-1]): raise InputError('Temperature is outside the source thermal-conductivity table. No extrapolation was applied.')
        lam=np.full_like(th,self.lam)
        if lt and lw: lam=self.lam*interp(th,self.f['lambda(Theta_l)'])*interp(T,self.f['lambda(T)'])
        elif lt: lam=interp(T,self.f['lambda(T)'])
        elif lw: lam=interp(th,self.f['lambda(Theta_l)'])
        if np.any(kv<=0) or np.any(kl<=0) or np.any(lam<=0): raise InputError('Non-positive transport function.')
        return p,rh,kl,kv,lam

def read_material(path, branch='de'):
    p=Path(path)
    if not p.is_file():
        raise InputError(f'Missing material curve file {p.name}. Use Download selected brick source or import the .m6 file. No replacement curves will be generated.')
    if p.suffix.lower()=='.m6': obj=parse_m6(p,branch)
    elif p.suffix.lower()=='.json':
        obj=json.loads(p.read_text()); obj['source_sha256']=sha256(p); obj['retention_branch']=branch
    else: raise InputError('Material must be tabulated .json or ASCII .m6.')
    return Material(obj)
