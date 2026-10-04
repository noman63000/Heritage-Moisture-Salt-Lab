from pathlib import Path
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
import zlib, hashlib, json, shutil, subprocess

BASE=Path('/mnt/data')
NEX=BASE/'cmip6_processed/Mohenjo_daro_NEX_CMIP6_daily_all_models.csv.gz'
ERA_RAW=BASE/'Mohenjo_daro_ERA5Land_2010_2023_raw_merged.csv.gz'
OUT=BASE/'hourly_future_v6'
OUT.mkdir(exist_ok=True)
FEATURE_COLS=['tas_C','dtr_C','hurs_pct','pr_mm_day','rsds_W_m2','rlds_W_m2','sfcWind_m_s']
WEIGHTS=np.array([1,1,1,1,1,1,1,0.35],float); K=5; WET=0.01

def sha256(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()

def sigmoid(x): return 1/(1+np.exp(-np.clip(x,-60,60)))
def logit(x):
 x=np.clip(x,1e-8,1-1e-8); return np.log(x/(1-x))
def rh_magnus(T,Td): return np.clip(100*np.exp(17.625*Td/(243.04+Td)-17.625*T/(243.04+T)),0,100)
def dewpoint(T,rh):
 r=np.clip(rh/100,1e-8,1); g=np.log(r)+17.625*T/(243.04+T); return 243.04*g/(17.625-g)
def vp(T,rh): return rh/100*610.94*np.exp(17.625*T/(243.04+T))

def mean_logit(base,target,it=32):
 tgt=np.clip(target,0,1); out=np.empty_like(base,float)
 z=tgt<=1e-12; o=tgt>=1-1e-12; m=~(z|o); out[z]=0; out[o]=1
 if m.any():
  L=logit(base[m]); lo=np.full(L.shape[0],-30.); hi=np.full(L.shape[0],30.); tm=tgt[m]
  for _ in range(it):
   c=(lo+hi)/2; mm=sigmoid(L+c[:,None]).mean(1); q=mm<tm; lo[q]=c[q]; hi[~q]=c[~q]
  out[m]=sigmoid(L+((lo+hi)/2)[:,None])
 return out

def temp_shape(a):
 mn=a.min(); mx=a.max()
 if mx-mn<1e-9:return np.full(24,.5)
 s=(a-mn)/(mx-mn); i=(s>1e-10)&(s<1-1e-10)
 if not i.any():return s
 L=logit(s[i]); lo=-30.; hi=30.
 for _ in range(40):
  c=(lo+hi)/2; f=s.copy(); f[i]=sigmoid(L+c)
  if f.mean()<.5:lo=c
  else:hi=c
 f=s.copy(); f[i]=sigmoid(L+(lo+hi)/2); return f

# ERA5
e=pd.read_csv(ERA_RAW,parse_dates=['valid_time']); T=e.t2m.to_numpy()-273.15; Td=e.d2m.to_numpy()-273.15
eh=pd.DataFrame({'date':e.valid_time.dt.strftime('%Y-%m-%d'),'month':e.valid_time.dt.month,'day':e.valid_time.dt.day,'hour':e.valid_time.dt.hour,
 'tas_C':T,'hurs_pct':rh_magnus(T,Td),'pr_mm_h':e.tp.to_numpy()*1000,'rsds_W_m2':e.ssrd.to_numpy()/3600,
 'rlds_W_m2':e.strd.to_numpy()/3600,'sfcWind_m_s':np.hypot(e.u10.to_numpy(),e.v10.to_numpy())})
assert (eh.groupby('date').size()==24).all()
ed=eh.groupby('date').agg(month=('month','first'),day=('day','first'),tas_C=('tas_C','mean'),tasmax_C=('tas_C','max'),tasmin_C=('tas_C','min'),hurs_pct=('hurs_pct','mean'),pr_mm_day=('pr_mm_h','sum'),rsds_W_m2=('rsds_W_m2','mean'),rlds_W_m2=('rlds_W_m2','mean'),sfcWind_m_s=('sfcWind_m_s','mean')).reset_index()
ed['dtr_C']=ed.tasmax_C-ed.tasmin_C; ed['era_idx']=np.arange(len(ed)); era_dates=ed.date.to_numpy()
def pm(c): return eh.pivot(index='date',columns='hour',values=c).reindex(ed.date).to_numpy(float)
Tm=pm('tas_C'); RHm=pm('hurs_pct')/100; Pm=pm('pr_mm_h'); SWm=pm('rsds_W_m2'); LWm=pm('rlds_W_m2'); Wm=pm('sfcWind_m_s')
TS=np.vstack([temp_shape(x) for x in Tm])
SWr=SWm/np.where(SWm.mean(1)[:,None]>1e-12,SWm.mean(1)[:,None],1)
LWr=LWm/np.where(LWm.mean(1)[:,None]>1e-12,LWm.mean(1)[:,None],1)
Wr=Wm/np.where(Wm.mean(1)[:,None]>1e-12,Wm.mean(1)[:,None],1)
Ps=Pm.sum(1); Pf=np.divide(Pm,Ps[:,None],out=np.zeros_like(Pm),where=Ps[:,None]>1e-12)
ES={}
for m,g in ed.groupby('month'):
 gg=g.reset_index(drop=True); pct=[gg[c].rank(pct=True,method='average').to_numpy() for c in FEATURE_COLS]; df=(gg.day.to_numpy()-.5)/gg.day.max(); X=np.column_stack(pct+[df])*WEIGHTS
 wet=np.where(gg.pr_mm_day.to_numpy()>WET)[0]; Xw=np.column_stack([pct[3],pct[2],pct[4],df])
 ES[m]=(gg,cKDTree(X),wet,cKDTree(Xw[wet]) if len(wet) else None)

f=pd.read_csv(NEX,dtype={'date':str}); f['dtr_C']=f.tasmax_C-f.tasmin_C
maps=[]; qas=[]; mans=[]; sums=[]
for (model,scenario,period),g0 in f.groupby(['model','scenario','period'],sort=False):
 g=g0.reset_index(drop=True); n=len(g); ai=np.empty(n,int); pi=np.empty(n,int); ad=np.empty(n,float); ar=np.empty(n,int)
 for m,inds0 in g.groupby('month').groups.items():
  inds=np.array(list(inds0),int); gm=g.loc[inds]; pct=[gm[c].rank(pct=True,method='average').to_numpy() for c in FEATURE_COLS]; df=(gm.day.to_numpy()-.5)/gm.day.max(); Xq=np.column_stack(pct+[df])*WEIGHTS
  gg,tree,wet,wtree=ES[m]; kk=min(K,len(gg)); dist,nn=tree.query(Xq,k=kk); dist=np.atleast_2d(dist).reshape(len(inds),kk); nn=np.atleast_2d(nn).reshape(len(inds),kk)
  keys=[f'{model}|{scenario}|{period}|{d}' for d in gm.date.to_numpy()]
  choices=np.fromiter((zlib.crc32(k.encode())%kk for k in keys),dtype=int,count=len(keys))
  rows=np.arange(len(inds)); locs=nn[rows,choices].astype(int)
  era_ids=gg.era_idx.to_numpy(dtype=int)
  ai[inds]=era_ids[locs]; ad[inds]=dist[rows,choices]; ar[inds]=choices+1; pi[inds]=ai[inds]
  need=(gm.pr_mm_day.to_numpy()>1e-12) & (gg.pr_mm_day.to_numpy()[locs]<=WET)
  if need.any() and wtree is not None:
   qmat=np.column_stack([pct[3],pct[2],pct[4],df])[need]
   kw=min(K,len(wet)); _,nw=wtree.query(qmat,k=kw); nw=np.asarray(nw)
   if kw==1: nw=nw.reshape(-1,1)
   keys_need=np.asarray(keys,dtype=object)[need]
   c2=np.fromiter((zlib.crc32((k+'|pr').encode())%kw for k in keys_need),dtype=int,count=len(keys_need))
   rr=np.arange(len(keys_need)); wetloc=wet[nw[rr,c2].astype(int)]
   pi[inds[need]]=era_ids[wetloc]
 # hourly matrices
 tmin=g.tasmin_C.to_numpy()[:,None]; tmax=g.tasmax_C.to_numpy()[:,None]; Th=tmin+(tmax-tmin)*TS[ai]
 resid=g.tas_C.to_numpy()-Th.mean(1)
 # vectorized tiny mean correction using all non-extreme hours; NEX residuals are <<0.001 C
 masks=(TS[ai]>1e-10)&(TS[ai]<1-1e-10); cnt=masks.sum(1); adj=np.divide(resid*24,cnt,out=np.zeros(n),where=cnt>0); Th+=masks*adj[:,None]
 RH=mean_logit(RHm[ai],g.hurs_pct.to_numpy()/100)*100
 Ph=Pf[pi]*g.pr_mm_day.to_numpy()[:,None]
 bad=(g.pr_mm_day.to_numpy()>1e-12)&(Ph.sum(1)<=1e-12)
 if bad.any(): Ph[bad,12:15]=g.loc[bad,'pr_mm_day'].to_numpy()[:,None]*np.array([.25,.5,.25])
 SW=SWr[ai]*g.rsds_W_m2.to_numpy()[:,None]; LW=LWr[ai]*g.rlds_W_m2.to_numpy()[:,None]; WW=Wr[ai]*g.sfcWind_m_s.to_numpy()[:,None]
 TD=dewpoint(Th,RH); VP=vp(Th,RH)
 # write compact file; metadata is in manifest
 out=pd.DataFrame({'date':np.repeat(g.date.to_numpy(),24),'hour_utc':np.tile(np.arange(24),n),'tas_C':Th.ravel(),'hurs_pct':RH.ravel(),'dewpoint_C':TD.ravel(),'vapor_pressure_Pa':VP.ravel(),'pr_mm_h':Ph.ravel(),'rsds_W_m2':SW.ravel(),'rlds_W_m2':LW.ravel(),'sfcWind_m_s':WW.ravel()})
 p=OUT/f'Mohenjo_daro_hourly_{model}_{scenario}_{period}.csv.gz'
 # Reuse complete outputs from an interrupted earlier run; rewrite only incomplete/missing files.
 expected_rows=len(out)
 complete=False
 if p.exists() and p.stat().st_size>1000000:
  try:
   rr=subprocess.run(['bash','-lc',f"gzip -dc '{p}' | wc -l"],capture_output=True,text=True,check=True)
   complete=(int(rr.stdout.strip())==expected_rows+1)
  except Exception:
   complete=False
 if not complete:
  tmp=p.with_suffix('')  # .csv
  out.to_csv(tmp,index=False,float_format='%.6f')
  subprocess.run(['gzip','-1','-f',str(tmp)],check=True)
  made=Path(str(tmp)+'.gz')
  if made!=p:
   if p.exists(): p.unlink()
   made.rename(p)
 # mapping vectorized
 maps.append(pd.DataFrame({'date':g.date,'model':model,'scenario':scenario,'period':period,'calendar':g.calendar,'analog_date_ERA5':era_dates[ai],'precip_analog_date_ERA5':era_dates[pi],'analog_neighbor_rank':ar,'analog_distance_quantile_space':ad}))
 errs={'tmean':np.abs(Th.mean(1)-g.tas_C.to_numpy()),'tmin':np.abs(Th.min(1)-g.tasmin_C.to_numpy()),'tmax':np.abs(Th.max(1)-g.tasmax_C.to_numpy()),'rh':np.abs(RH.mean(1)-g.hurs_pct.to_numpy()),'pr':np.abs(Ph.sum(1)-g.pr_mm_day.to_numpy()),'sw':np.abs(SW.mean(1)-g.rsds_W_m2.to_numpy()),'lw':np.abs(LW.mean(1)-g.rlds_W_m2.to_numpy()),'w':np.abs(WW.mean(1)-g.sfcWind_m_s.to_numpy())}
 qas.append({'model':model,'scenario':scenario,'period':period,'days':n,'hourly_rows':n*24,'max_abs_tmean_error_C':errs['tmean'].max(),'max_abs_tmin_error_C':errs['tmin'].max(),'max_abs_tmax_error_C':errs['tmax'].max(),'max_abs_rhmean_error_pct':errs['rh'].max(),'max_abs_precip_sum_error_mm':errs['pr'].max(),'max_abs_rsds_mean_error_W_m2':errs['sw'].max(),'max_abs_rlds_mean_error_W_m2':errs['lw'].max(),'max_abs_wind_mean_error_m_s':errs['w'].max(),'hourly_tas_min_C':Th.min(),'hourly_tas_max_C':Th.max(),'hourly_rh_min_pct':RH.min(),'hourly_rh_max_pct':RH.max(),'hourly_pr_max_mm_h':Ph.max(),'hourly_rsds_max_W_m2':SW.max(),'hourly_rlds_min_W_m2':LW.min(),'hourly_rlds_max_W_m2':LW.max(),'hourly_wind_max_m_s':WW.max(),'unique_general_analogs':len(np.unique(ai)),'unique_precip_analogs':len(np.unique(pi)),'analog_distance_median':np.median(ad),'analog_distance_p95':np.quantile(ad,.95)})
 sums.append({'model':model,'scenario':scenario,'period':period,'days':n,'mean_tas_C':Th.mean(),'absolute_max_tas_C':Th.max(),'absolute_min_tas_C':Th.min(),'mean_hurs_pct':RH.mean(),'total_precip_mm':Ph.sum(),'mean_rsds_W_m2':SW.mean(),'mean_rlds_W_m2':LW.mean(),'mean_wind_m_s':WW.mean(),'max_hourly_precip_mm_h':Ph.max()})
 mans.append({'model':model,'scenario':scenario,'period':period,'file':p.name,'rows':len(out),'size_bytes':p.stat().st_size,'sha256':sha256(p),'calendar':g.calendar.iloc[0],'grid_lat':g.grid_lat.iloc[0],'grid_lon':g.grid_lon.iloc[0]})
 print('done',model,scenario,period,n,'days',round(p.stat().st_size/1e6,2),'MB',flush=True)

mapping=pd.concat(maps,ignore_index=True); mapping.to_csv(OUT/'Mohenjo_daro_hourly_disaggregation_daily_mapping.csv.gz',index=False,compression={'method':'gzip','compresslevel':1},float_format='%.6f')
qa=pd.DataFrame(qas); qa.to_csv(OUT/'Mohenjo_daro_hourly_forcing_QA.csv',index=False,float_format='%.10g')
man=pd.DataFrame(mans); man.to_csv(OUT/'Mohenjo_daro_hourly_forcing_manifest.csv',index=False)
su=pd.DataFrame(sums); su.to_csv(OUT/'Mohenjo_daro_hourly_forcing_period_summary.csv',index=False,float_format='%.6f')
meta={'method':'ERA5-Land same-month quantile-analog temporal disaggregation of daily NEX-GDDP-CMIP6 v2 forcing','era5_baseline':'2010-2023 hourly, ERA5-Land grid 27.3N 68.1E','nex_grid':'27.375N 68.125E','hour_utc_note':'Hourly profiles are retained on the NEX model day in UTC labels; Pakistan Standard Time is UTC+5. The ERA5 radiation profile preserves observed solar timing.','analog_features':FEATURE_COLS+['within_month_day_fraction'],'analog_feature_weights':WEIGHTS.tolist(),'nearest_neighbors':K,'neighbor_selection':'deterministic CRC32 selection among 5 nearest ERA5 analogs','temperature':'ERA5 hourly temperature shape normalized to exact NEX daily Tmin/Tmax and corrected to exact NEX daily mean','relative_humidity':'ERA5 hourly RH pattern shifted in logit space to exact NEX daily mean while bounded 0-100%','precipitation':'NEX daily total distributed by same-month ERA5 wet-day hourly fractions; exact daily total preserved','shortwave':'ERA5 hourly shortwave shape multiplicatively scaled to exact NEX daily mean; night zeros preserved','longwave':'ERA5 hourly longwave shape multiplicatively scaled to exact NEX daily mean','wind':'ERA5 hourly wind-speed shape multiplicatively scaled to exact NEX daily mean','dewpoint_and_vapor_pressure':'derived from hourly T and RH using Magnus relation','caution':'Statistical temporal disaggregation: future subdaily storm sequencing/intensity is not independently projected. Run the nonlinear deterioration model separately for each GCM/scenario/period, then summarize the deterioration outputs across the ensemble.','hourly_output_files':len(man),'hourly_rows_total':int(man.rows.sum()),'daily_mapping_rows':len(mapping)}
with open(OUT/'Mohenjo_daro_hourly_disaggregation_metadata.json','w') as fh: json.dump(meta,fh,indent=2)
print('TOTAL',int(man.rows.sum()),'hourly rows',flush=True)
print(qa.filter(regex='max_abs').max().to_string(),flush=True)
print('MAX PHYSICAL')
print(qa[['hourly_tas_min_C','hourly_tas_max_C','hourly_rh_min_pct','hourly_rh_max_pct','hourly_pr_max_mm_h','hourly_rsds_max_W_m2','hourly_rlds_min_W_m2','hourly_rlds_max_W_m2','hourly_wind_max_m_s']].agg(['min','max']).to_string(),flush=True)
