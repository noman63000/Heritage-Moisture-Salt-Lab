from pathlib import Path
import pandas as pd, numpy as np, hashlib, json, glob, subprocess
BASE=Path('/mnt/data'); OUT=BASE/'hourly_future_v6'; NEX=BASE/'cmip6_processed/Mohenjo_daro_NEX_CMIP6_daily_all_models.csv.gz'; CAP=1200.0

def sha256(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()

nex=pd.read_csv(NEX,dtype={'date':str})
# direct lookup by filename key
groups={(m,s,p):g.reset_index(drop=True) for (m,s,p),g in nex.groupby(['model','scenario','period'],sort=False)}
files=sorted([Path(x) for x in glob.glob(str(OUT/'Mohenjo_daro_hourly_*.csv.gz')) if 'mapping' not in x])
qa=[]; man=[]; sums=[]; adj=[]
for p in files:
 key=None
 for k in groups:
  if p.name==f'Mohenjo_daro_hourly_{k[0]}_{k[1]}_{k[2]}.csv.gz': key=k; break
 if key is None: continue
 g=groups[key]; m,s,per=key; df=pd.read_csv(p); n=len(g)
 if len(df)!=n*24: raise RuntimeError(f'row mismatch {p.name}: {len(df)} vs {n*24}')
 # gzip integrity
 subprocess.run(['gzip','-t',str(p)],check=True)
 dates=df.date.astype(str).to_numpy().reshape(n,24); hrs=df.hour_utc.to_numpy().reshape(n,24)
 if not np.all(dates[:,0]==g.date.astype(str).to_numpy()): raise RuntimeError('date mismatch '+p.name)
 if not np.all(hrs==np.arange(24)[None,:]): raise RuntimeError('hour mismatch '+p.name)
 A={c:df[c].to_numpy(float).reshape(n,24) for c in ['tas_C','hurs_pct','pr_mm_h','rsds_W_m2','rlds_W_m2','sfcWind_m_s','dewpoint_C','vapor_pressure_Pa']}
 errs={'tmean':np.abs(A['tas_C'].mean(1)-g.tas_C.to_numpy()),'tmin':np.abs(A['tas_C'].min(1)-g.tasmin_C.to_numpy()),'tmax':np.abs(A['tas_C'].max(1)-g.tasmax_C.to_numpy()),'rh':np.abs(A['hurs_pct'].mean(1)-g.hurs_pct.to_numpy()),'pr':np.abs(A['pr_mm_h'].sum(1)-g.pr_mm_day.to_numpy()),'sw':np.abs(A['rsds_W_m2'].mean(1)-g.rsds_W_m2.to_numpy()),'lw':np.abs(A['rlds_W_m2'].mean(1)-g.rlds_W_m2.to_numpy()),'w':np.abs(A['sfcWind_m_s'].mean(1)-g.sfcWind_m_s.to_numpy())}
 capped=(A['rsds_W_m2'].max(1)>=CAP-1e-5)
 qa.append({'model':m,'scenario':s,'period':per,'days':n,'hourly_rows':len(df),'max_abs_tmean_error_C':errs['tmean'].max(),'max_abs_tmin_error_C':errs['tmin'].max(),'max_abs_tmax_error_C':errs['tmax'].max(),'max_abs_rhmean_error_pct':errs['rh'].max(),'max_abs_precip_sum_error_mm':errs['pr'].max(),'max_abs_rsds_mean_error_W_m2':errs['sw'].max(),'max_abs_rlds_mean_error_W_m2':errs['lw'].max(),'max_abs_wind_mean_error_m_s':errs['w'].max(),'hourly_tas_min_C':A['tas_C'].min(),'hourly_tas_max_C':A['tas_C'].max(),'hourly_rh_min_pct':A['hurs_pct'].min(),'hourly_rh_max_pct':A['hurs_pct'].max(),'hourly_pr_max_mm_h':A['pr_mm_h'].max(),'hourly_rsds_max_W_m2':A['rsds_W_m2'].max(),'hourly_rlds_min_W_m2':A['rlds_W_m2'].min(),'hourly_rlds_max_W_m2':A['rlds_W_m2'].max(),'hourly_wind_max_m_s':A['sfcWind_m_s'].max(),'negative_values_count':int((A['pr_mm_h']<0).sum()+(A['rsds_W_m2']<0).sum()+(A['rlds_W_m2']<0).sum()+(A['sfcWind_m_s']<0).sum()),'rh_out_of_bounds_count':int(((A['hurs_pct']<0)|(A['hurs_pct']>100)).sum()),'days_shortwave_capped':int(capped.sum())})
 sums.append({'model':m,'scenario':s,'period':per,'days':n,'mean_tas_C':A['tas_C'].mean(),'absolute_max_tas_C':A['tas_C'].max(),'absolute_min_tas_C':A['tas_C'].min(),'mean_hurs_pct':A['hurs_pct'].mean(),'total_precip_mm':A['pr_mm_h'].sum(),'mean_rsds_W_m2':A['rsds_W_m2'].mean(),'mean_rlds_W_m2':A['rlds_W_m2'].mean(),'mean_wind_m_s':A['sfcWind_m_s'].mean(),'max_hourly_precip_mm_h':A['pr_mm_h'].max(),'max_hourly_rsds_W_m2':A['rsds_W_m2'].max()})
 man.append({'model':m,'scenario':s,'period':per,'file':p.name,'rows':len(df),'size_bytes':p.stat().st_size,'sha256':sha256(p),'calendar':g.calendar.iloc[0],'grid_lat':g.grid_lat.iloc[0],'grid_lon':g.grid_lon.iloc[0]})
 adj.append({'model':m,'scenario':s,'period':per,'days_shortwave_capped':int(capped.sum()),'cap_W_m2':CAP,'max_after_W_m2':float(A['rsds_W_m2'].max())})
 print('checked',p.name,'rows',len(df),'capped days',int(capped.sum()),flush=True)

q=pd.DataFrame(qa); ma=pd.DataFrame(man); su=pd.DataFrame(sums); ad=pd.DataFrame(adj)
q.to_csv(OUT/'Mohenjo_daro_hourly_forcing_QA.csv',index=False,float_format='%.10g')
ma.to_csv(OUT/'Mohenjo_daro_hourly_forcing_manifest.csv',index=False)
su.to_csv(OUT/'Mohenjo_daro_hourly_forcing_period_summary.csv',index=False,float_format='%.6f')
ad.to_csv(OUT/'Mohenjo_daro_shortwave_cap_adjustments.csv',index=False,float_format='%.6f')
mp=OUT/'Mohenjo_daro_hourly_disaggregation_metadata.json'; meta=json.load(open(mp))
meta['shortwave_physical_cap_W_m2']=CAP
meta['shortwave_cap_method']='Scaled ERA5 shortwave values above 1200 W/m2 were capped and removed energy redistributed across that day\'s positive-radiation hours in proportion to remaining capacity. This preserves the exact NEX daily shortwave mean and night zeros.'
meta['shortwave_days_adjusted_total']=int(ad.days_shortwave_capped.sum())
meta['qa_after_final_write']={'hourly_files':len(ma),'hourly_rows_total':int(ma.rows.sum()),'max_daily_constraint_error':{c:float(q[c].max()) for c in q.columns if c.startswith('max_abs_')},'max_hourly_shortwave_W_m2':float(q.hourly_rsds_max_W_m2.max()),'negative_values_count_total':int(q.negative_values_count.sum()),'rh_out_of_bounds_count_total':int(q.rh_out_of_bounds_count.sum())}
with open(mp,'w') as f: json.dump(meta,f,indent=2)
print('\nFINAL')
print('files',len(ma),'rows',int(ma.rows.sum()),'capped days',int(ad.days_shortwave_capped.sum()))
print(q.filter(regex='max_abs').max().to_string())
print('max hourly: T',q.hourly_tas_max_C.max(),'RH',q.hourly_rh_max_pct.max(),'rain',q.hourly_pr_max_mm_h.max(),'SW',q.hourly_rsds_max_W_m2.max(),'LW',q.hourly_rlds_max_W_m2.max(),'wind',q.hourly_wind_max_m_s.max())
print('negative',q.negative_values_count.sum(),'rh_oob',q.rh_out_of_bounds_count.sum())
