import copy, json, sys, time
from pathlib import Path
import numpy as np

# This verification script is run against a full verified Heritage Lab app by
# setting HERITAGE_APP_ROOT. It intentionally does not bundle climate archives.
ROOT=Path(__import__('os').environ.get('HERITAGE_APP_ROOT','')).resolve()
if not ROOT.is_dir():
    raise SystemExit('Set HERITAGE_APP_ROOT to a verified Heritage Lab application folder.')
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'payload'))
from hmsl.project import Project
from hmsl.workflow import load_state
from hmsl.solver import ColumnModel, PhysicsScopeError
from hlab_performance.phase_change import EquilibriumFreezeThawModel

p=Project(ROOT/'projects/mohenjo_daro'); seed=load_state(p.path/'spinup.npz')
records=[]
def check(name,ok,detail=''):
    records.append({'name':name,'passed':bool(ok),'detail':str(detail)})
    if not ok: raise AssertionError(name+': '+str(detail))

# 1: no-ice arithmetic identity.
c=p.climate('R00').subset(0,24)
a=ColumnModel(copy.deepcopy(p.config),p.material(),p.phase()).integrate(c,seed)[0]
b=EquilibriumFreezeThawModel(copy.deepcopy(p.config),p.material(),p.phase()).integrate(c,seed)[0]
check('Warm 24 h final state is bit-for-bit identical',np.array_equal(a.matrix(),b.matrix()))

# 2: actual R14 cold window that v1.3 stopped; v1.4 forms and thaws ice.
c14=p.climate('R14').subset(175011,175011+288)
t=time.time();m=EquilibriumFreezeThawModel(copy.deepcopy(p.config),p.material(),p.phase());st,h,pr,stats=m.integrate(c14,seed)
check('R14 288 h cold window completes',len(h)==288)
check('R14 window activates phase model',stats['phase_active_hours']>0,stats['phase_active_hours'])
check('R14 water closure',float(h.water_closure_error_kg_m2.abs().max())<1e-8)
check('R14 sulfate closure',float(h.sulfate_closure_error_mol_m2.abs().max())<1e-10)
check('R14 actually forms ice-equivalent water',stats['maximum_ice_water_equivalent_m3_m3']>0,stats['maximum_ice_water_equivalent_m3_m3'])
records[-1]['seconds']=time.time()-t

# 3: checkpoint/resume equivalence around same phase event.
m1=EquilibriumFreezeThawModel(copy.deepcopy(p.config),p.material(),p.phase()); mid=m1.integrate(c14.subset(0,144),seed)[0]
m2=EquilibriumFreezeThawModel(copy.deepcopy(p.config),p.material(),p.phase()); split=m2.integrate(c14.subset(144,288),mid)[0]
check('R14 split resume final state is bit-for-bit identical',np.array_equal(st.matrix(),split.matrix()))

# 4: three other representative cold windows complete with the phase-capable engine.
for rid,start in [('R17',78407),('R18',113183),('R20',122714)]:
    cc=p.climate(rid).subset(start,start+288); mm=EquilibriumFreezeThawModel(copy.deepcopy(p.config),p.material(),p.phase()); ss,hh,pp,rr=mm.integrate(cc,seed)
    check(rid+' representative cold window completes',len(hh)==288)
    check(rid+' water closure',float(hh.water_closure_error_kg_m2.abs().max())<1e-8)

# 5: phase equilibrium increases liquid salt concentration because ice excludes salt.
mm=EquilibriumFreezeThawModel(copy.deepcopy(p.config),p.material(),p.phase())
tot=np.full(mm.n,0.025); temp=np.full(mm.n,-4.0); liq,ice,_=mm.phase_partition(tot,temp)
ns=np.full(mm.n,0.1); c,_=ColumnModel.partition(mm,liq,temp,ns)
check('Cold equilibrium has positive ice at theta=0.025,T=-4C',float(ice.max())>0)
check('Salt concentration rises when part of pore water is ice',float(c.mean())>float((ns/tot).mean()))

# 6: declared phase range remains guarded.
try:
    mm.phase_partition(tot,np.full(mm.n,-20.1)); ok=False
except PhysicsScopeError: ok=True
check('Below -20 C remains outside declared phase range',ok)

out=Path(__file__).resolve().parents[1]/'verification'/'phase_extension_checks.json'
out.write_text(json.dumps({'passed':all(x['passed'] for x in records),'checks':records},indent=2))
print(json.dumps({'passed':True,'checks':len(records),'r14_phase_hours':stats['phase_active_hours'],'r14_max_ice':stats['maximum_ice_water_equivalent_m3_m3']},indent=2))
