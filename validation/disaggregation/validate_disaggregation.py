import pandas as pd, numpy as np, zlib, json
from scipy.spatial import cKDTree
from pathlib import Path

P=Path('/mnt/data/q1_revision/Mohenjo_daro_ERA5Land_2010_2023_hourly_model_ready.csv.gz')
e=pd.read_csv(P, parse_dates=['timestamp_utc'])
e['year']=e.timestamp_utc.dt.year;e['month']=e.timestamp_utc.dt.month;e['day']=e.timestamp_utc.dt.day;e['date_s']=e.timestamp_utc.dt.strftime('%Y-%m-%d');e['hour']=e.timestamp_utc.dt.hour
# arrays / daily targets
ed=e.groupby('date_s').agg(year=('year','first'),month=('month','first'),day=('day','first'),tas_C=('tas_C','mean'),tasmax_C=('tas_C','max'),tasmin_C=('tas_C','min'),hurs_pct=('hurs_pct','mean'),pr_mm_day=('pr_mm_h','sum'),rsds_W_m2=('rsds_W_m2','mean'),rlds_W_m2=('rlds_W_m2','mean'),sfcWind_m_s=('sfcWind_m_s','mean')).reset_index()
ed['dtr_C']=ed.tasmax_C-ed.tasmin_C
FEATURE_COLS=['tas_C','dtr_C','hurs_pct','pr_mm_day','rsds_W_m2','rlds_W_m2','sfcWind_m_s']
WEIGHTS=np.array([1,1,1,1,1,1,1,0.35],float);K=5;WET=0.01

def pm(col): return e.pivot(index='date_s',columns='hour',values=col).reindex(ed.date_s).to_numpy(float)
Tm=pm('tas_C'); RHm=pm('hurs_pct')/100; Pm=pm('pr_mm_h'); SWm=pm('rsds_W_m2'); LWm=pm('rlds_W_m2'); Wm=pm('sfcWind_m_s')

def sigmoid(x): return 1/(1+np.exp(-np.clip(x,-60,60)))
def logit(x): x=np.clip(x,1e-8,1-1e-8); return np.log(x/(1-x))
def mean_logit(base,target,it=32):
    tgt=np.clip(target,0,1); out=np.empty_like(base,float); z=tgt<=1e-12;o=tgt>=1-1e-12;m=~(z|o);out[z]=0;out[o]=1
    if m.any():
        L=logit(base[m]);lo=np.full(L.shape[0],-30.);hi=np.full(L.shape[0],30.);tm=tgt[m]
        for _ in range(it):
            c=(lo+hi)/2;mm=sigmoid(L+c[:,None]).mean(1);q=mm<tm;lo[q]=c[q];hi[~q]=c[~q]
        out[m]=sigmoid(L+((lo+hi)/2)[:,None])
    return out

def temp_shape(a):
    mn=a.min();mx=a.max()
    if mx-mn<1e-9:return np.full(24,.5)
    s=(a-mn)/(mx-mn);i=(s>1e-10)&(s<1-1e-10)
    if not i.any():return s
    L=logit(s[i]);lo=-30.;hi=30.
    for _ in range(40):
        c=(lo+hi)/2;f=s.copy();f[i]=sigmoid(L+c)
        if f.mean()<.5:lo=c
        else:hi=c
    f=s.copy();f[i]=sigmoid(L+(lo+hi)/2);return f

TS=np.vstack([temp_shape(x) for x in Tm])
SWr=SWm/np.where(SWm.mean(1)[:,None]>1e-12,SWm.mean(1)[:,None],1)
LWr=LWm/np.where(LWm.mean(1)[:,None]>1e-12,LWm.mean(1)[:,None],1)
Wr=Wm/np.where(Wm.mean(1)[:,None]>1e-12,Wm.mean(1)[:,None],1)
Ps=Pm.sum(1);Pf=np.divide(Pm,Ps[:,None],out=np.zeros_like(Pm),where=Ps[:,None]>1e-12)

pred={k:np.full_like(v,np.nan,dtype=float) for k,v in [('tas_C',Tm),('hurs_pct',RHm*100),('pr_mm_h',Pm),('rsds_W_m2',SWm),('rlds_W_m2',LWm),('sfcWind_m_s',Wm)]}
analog_dates=np.empty(len(ed),dtype=object)

for yr in sorted(ed.year.unique()):
    target_idxs=np.flatnonzero(ed.year.to_numpy()==yr)
    lib_idxs=np.flatnonzero(ed.year.to_numpy()!=yr)
    for m in range(1,13):
        tind=target_idxs[ed.month.to_numpy()[target_idxs]==m]
        lind=lib_idxs[ed.month.to_numpy()[lib_idxs]==m]
        if len(tind)==0: continue
        gtar=ed.iloc[tind].copy(); glib=ed.iloc[lind].copy()
        # independent percentile feature scaling within target and candidate sets, matching production method
        pt=[gtar[c].rank(pct=True,method='average').to_numpy() for c in FEATURE_COLS]
        pl=[glib[c].rank(pct=True,method='average').to_numpy() for c in FEATURE_COLS]
        dft=(gtar.day.to_numpy()-.5)/gtar.day.max(); dfl=(glib.day.to_numpy()-.5)/glib.day.max()
        Xt=np.column_stack(pt+[dft])*WEIGHTS; Xl=np.column_stack(pl+[dfl])*WEIGHTS
        tree=cKDTree(Xl);kk=min(K,len(lind));dist,nn=tree.query(Xt,k=kk)
        if kk==1: nn=nn[:,None]
        keys=[f'LOYO|{yr}|{d}' for d in gtar.date_s]
        choices=np.fromiter((zlib.crc32(k.encode())%kk for k in keys),dtype=int,count=len(keys))
        locs=nn[np.arange(len(tind)),choices].astype(int); ai=lind[locs]; pi=ai.copy()
        # precipitation wet analog fallback, matching production dimensions
        wetloc=np.flatnonzero(glib.pr_mm_day.to_numpy()>WET)
        need=(gtar.pr_mm_day.to_numpy()>1e-12)&(glib.pr_mm_day.to_numpy()[locs]<=WET)
        if need.any() and len(wetloc):
            Xlw=np.column_stack([pl[3],pl[2],pl[4],dfl])[wetloc]
            treew=cKDTree(Xlw); Xtw=np.column_stack([pt[3],pt[2],pt[4],dft])[need]
            kw=min(K,len(wetloc));_,nw=treew.query(Xtw,k=kw)
            if kw==1:nw=nw[:,None]
            kn=np.asarray(keys,dtype=object)[need]
            c2=np.fromiter((zlib.crc32((k+'|pr').encode())%kw for k in kn),dtype=int,count=len(kn))
            pick=wetloc[nw[np.arange(len(kn)),c2].astype(int)]
            pi[need]=lind[pick]
        # reconstruction
        tmin=gtar.tasmin_C.to_numpy()[:,None];tmax=gtar.tasmax_C.to_numpy()[:,None]
        Th=tmin+(tmax-tmin)*TS[ai];resid=gtar.tas_C.to_numpy()-Th.mean(1);masks=(TS[ai]>1e-10)&(TS[ai]<1-1e-10);cnt=masks.sum(1);adj=np.divide(resid*24,cnt,out=np.zeros(len(tind)),where=cnt>0);Th+=masks*adj[:,None]
        RH=mean_logit(RHm[ai],gtar.hurs_pct.to_numpy()/100)*100
        Ph=Pf[pi]*gtar.pr_mm_day.to_numpy()[:,None]
        SW=SWr[ai]*gtar.rsds_W_m2.to_numpy()[:,None];LW=LWr[ai]*gtar.rlds_W_m2.to_numpy()[:,None];WW=Wr[ai]*gtar.sfcWind_m_s.to_numpy()[:,None]
        for k,a in [('tas_C',Th),('hurs_pct',RH),('pr_mm_h',Ph),('rsds_W_m2',SW),('rlds_W_m2',LW),('sfcWind_m_s',WW)]: pred[k][tind]=a
        analog_dates[tind]=ed.date_s.to_numpy()[ai]

# metrics
obs={'tas_C':Tm,'hurs_pct':RHm*100,'pr_mm_h':Pm,'rsds_W_m2':SWm,'rlds_W_m2':LWm,'sfcWind_m_s':Wm}
rows=[]
for k in obs:
    o=obs[k].ravel();p=pred[k].ravel();rmse=float(np.sqrt(np.mean((p-o)**2)));mae=float(np.mean(np.abs(p-o)));bias=float(np.mean(p-o));r=float(np.corrcoef(o,p)[0,1])
    rows.append(dict(variable=k,rmse=rmse,mae=mae,bias=bias,pearson_r=r,obs_mean=float(o.mean()),pred_mean=float(p.mean()),obs_p95=float(np.quantile(o,.95)),pred_p95=float(np.quantile(p,.95)),obs_p99=float(np.quantile(o,.99)),pred_p99=float(np.quantile(p,.99))))
metrics=pd.DataFrame(rows)
# mean diurnal cycle errors and correlation for each variable
cycles=[]
for k in obs:
    om=obs[k].mean(0);pm=pred[k].mean(0); cycles.append(dict(variable=k,diurnal_rmse=float(np.sqrt(np.mean((pm-om)**2))),diurnal_r=float(np.corrcoef(om,pm)[0,1]),diurnal_max_abs=float(np.max(np.abs(pm-om)))))
diurnal=pd.DataFrame(cycles)
# precip distribution-specific validation preserving daily totals
wet_obs=(Pm>0.1);wet_pred=(pred['pr_mm_h']>0.1)
wet_hours_obs=wet_obs.sum(1);wet_hours_pred=wet_pred.sum(1)
wetdays=ed.pr_mm_day.to_numpy()>0.1
precip_summary={
 'daily_total_rmse_mm':float(np.sqrt(np.mean((pred['pr_mm_h'].sum(1)-Pm.sum(1))**2))),
 'wet_day_count':int(wetdays.sum()),
 'wet_hours_per_wet_day_obs_mean':float(wet_hours_obs[wetdays].mean()),
 'wet_hours_per_wet_day_pred_mean':float(wet_hours_pred[wetdays].mean()),
 'wet_hours_per_wet_day_mae':float(np.mean(np.abs(wet_hours_pred[wetdays]-wet_hours_obs[wetdays]))),
 'hourly_wet_frequency_obs_pct':float(wet_obs.mean()*100),
 'hourly_wet_frequency_pred_pct':float(wet_pred.mean()*100),
 'hourly_p95_wet_obs':float(np.quantile(Pm[Pm>0.1],.95)),
 'hourly_p95_wet_pred':float(np.quantile(pred['pr_mm_h'][pred['pr_mm_h']>0.1],.95)),
 'hourly_p99_wet_obs':float(np.quantile(Pm[Pm>0.1],.99)),
 'hourly_p99_wet_pred':float(np.quantile(pred['pr_mm_h'][pred['pr_mm_h']>0.1],.99)),
 'hourly_max_obs':float(Pm.max()),'hourly_max_pred':float(pred['pr_mm_h'].max())}
# extremes / daily shape retention
extreme={
 'tas_hours_gt40_obs':int((Tm>40).sum()),'tas_hours_gt40_pred':int((pred['tas_C']>40).sum()),
 'tas_hours_gt45_obs':int((Tm>45).sum()),'tas_hours_gt45_pred':int((pred['tas_C']>45).sum()),
 'rh_hours_gt80_obs':int((obs['hurs_pct']>80).sum()),'rh_hours_gt80_pred':int((pred['hurs_pct']>80).sum()),
 'rsds_hours_gt800_obs':int((SWm>800).sum()),'rsds_hours_gt800_pred':int((pred['rsds_W_m2']>800).sum())}
# month-level distribution stats
monthly=[]
for m in range(1,13):
    idx=ed.month.to_numpy()==m
    for k in obs:
        o=obs[k][idx].ravel();p=pred[k][idx].ravel(); monthly.append(dict(month=m,variable=k,rmse=float(np.sqrt(np.mean((p-o)**2))),r=float(np.corrcoef(o,p)[0,1]),obs_mean=float(o.mean()),pred_mean=float(p.mean())))
monthly=pd.DataFrame(monthly)

out=Path('/mnt/data/q1_revision/disaggregation_validation');out.mkdir(exist_ok=True)
metrics.to_csv(out/'loyo_hourly_metrics.csv',index=False);diurnal.to_csv(out/'loyo_diurnal_metrics.csv',index=False);monthly.to_csv(out/'loyo_monthly_metrics.csv',index=False)
pd.DataFrame({'date':ed.date_s,'year':ed.year,'analog_date':analog_dates}).to_csv(out/'loyo_analog_mapping.csv',index=False)
with open(out/'summary.json','w') as f: json.dump({'precipitation':precip_summary,'extremes':extreme,'days':len(ed),'hours':len(ed)*24,'method':'leave-one-year-out cross-validation; daily ERA5 targets reconstructed from analog days in other years only'},f,indent=2)
print(metrics.to_string(index=False))
print('\nDiurnal:\n',diurnal.to_string(index=False))
print('\nPrecip:\n',json.dumps(precip_summary,indent=2))
print('\nExtremes:\n',json.dumps(extreme,indent=2))
