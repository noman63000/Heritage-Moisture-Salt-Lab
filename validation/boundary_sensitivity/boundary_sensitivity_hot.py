import copy, sys, time
from pathlib import Path
import numpy as np, pandas as pd
ROOT=Path('model/Heritage_Lab_v1.1.0_REPAIR/payload').resolve(); sys.path.insert(0,str(ROOT))
from hmsl.project import Project
from hmsl.solver import ColumnModel
p=Project(ROOT/'projects/q1_sensitivity'); m=p.material(); rows=[]
for rid in ['R00','R08']:
    c=p.climate(rid); win=72
    # hottest complete 72-h window, excluding any with subzero irrelevant
    arr=c.data.tas_C.rolling(win).mean().to_numpy(); start=max(0,int(np.nanargmax(arr))-win+1); sub=c.subset(start,start+win)
    for supply in [0.1,0.2,0.4]:
        cfg=copy.deepcopy(p.config); cfg['boundary']['capillary_supply_mm_day']=supply
        model=ColumnModel(cfg,m,None); st=model.initial(); t=time.time(); final,h,profiles,rep=model.integrate(sub,st)
        rows.append({'run_id':rid,'window_type':'hottest_72h','start_timestamp':sub.data.timestamp_utc.iloc[0], 'end_timestamp':sub.data.timestamp_utc.iloc[-1],
        'supply_mm_day':supply,'mean_theta':h.theta_mean_m3_m3.mean(),'final_theta':h.theta_mean_m3_m3.iloc[-1], 'max_theta':h.theta_max_m3_m3.max(),
        'mean_T_C':h.T_mean_C.mean(),'max_T_C':h.T_max_C.max(),'final_water_kg_m2':h.water_inventory_kg_m2_footprint.iloc[-1],
        'final_salt_mol_m2':h.sulfate_total_mol_m2_footprint.iloc[-1], 'mean_cmax_band':h.sulfate_cmax_0p10m_bandavg_mol_m3_liquid.mean(), 'final_cmax_band':h.sulfate_cmax_0p10m_bandavg_mol_m3_liquid.iloc[-1],
        'basal_water_kg_m2':h.basal_water_kg_m2_hour.sum(),'basal_salt_mol_m2':h.basal_salt_mol_m2_hour.sum(), 'accepted_rain_kg_m2':h.accepted_rain_kg_m2_hour.sum(),'net_evap_kg_m2':h.net_evaporation_kg_m2_hour.sum(),
        'runtime_s':time.time()-t,'water_closure_max':h.water_closure_error_kg_m2.abs().max(),'salt_closure_max':h.sulfate_closure_error_mol_m2.abs().max()})
        print(rows[-1],flush=True)
out=pd.DataFrame(rows); out.to_csv('boundary_sensitivity_hot72.csv',index=False)
rels=[]
for rid,g in out.groupby('run_id'):
    b=g[g.supply_mm_day==0.2].iloc[0]
    for _,r in g.iterrows():
      d=r.to_dict()
      for k in ['mean_theta','final_theta','final_water_kg_m2','final_salt_mol_m2','mean_cmax_band','final_cmax_band']:
        d[k+'_pct_vs_nominal']=100*(r[k]/b[k]-1)
      rels.append(d)
pd.DataFrame(rels).to_csv('boundary_sensitivity_hot72_relative.csv',index=False)
