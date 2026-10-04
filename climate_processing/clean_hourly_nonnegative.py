from pathlib import Path
import pandas as pd, numpy as np, glob, subprocess, json
OUT=Path('/mnt/data/hourly_future_v6'); CAP=1200.0
stats=[]
for pstr in sorted(glob.glob(str(OUT/'Mohenjo_daro_hourly_*.csv.gz'))):
    if 'mapping' in pstr: continue
    p=Path(pstr); df=pd.read_csv(p); n=len(df)//24
    P=df.pr_mm_h.to_numpy(float).reshape(n,24)
    SW=df.rsds_W_m2.to_numpy(float).reshape(n,24)
    pneg=int((P<0).sum()); swneg=int((SW<0).sum())
    # precipitation: clamp deaccumulation artifacts to zero and rescale to exact previous daily total
    ptarget=P.sum(1)
    P=np.maximum(P,0)
    ps=P.sum(1)
    positive=ptarget>1e-12
    P[positive]*=(ptarget[positive]/ps[positive])[:,None]
    P[~positive]=0
    # shortwave: clamp tiny negative deaccumulation artifacts to zero and restore daily total
    swtarget=SW.sum(1)
    SW=np.maximum(SW,0)
    ss=SW.sum(1)
    good=swtarget>1e-12
    SW[good]*=(swtarget[good]/ss[good])[:,None]
    SW[~good]=0
    # enforce 1200 cap while preserving total and daylight zeros
    mask=SW.max(1)>CAP+1e-9
    for i in np.where(mask)[0]:
        row=SW[i].copy(); daylight=row>1e-12; target=swtarget[i]
        x=np.minimum(row,CAP); deficit=target-x.sum(); capacity=np.where(daylight,CAP-x,0); cs=capacity.sum()
        if deficit>cs+1e-6: raise RuntimeError(f'cap infeasible {p.name} {i}')
        if deficit>0: x += deficit*capacity/cs
        err=target-x.sum(); avail=np.where(daylight & (x<CAP-1e-9))[0]
        if abs(err)>1e-9 and len(avail): x[avail]+=err/len(avail)
        SW[i]=x
    df['pr_mm_h']=P.ravel(); df['rsds_W_m2']=SW.ravel()
    tmp=p.with_suffix(''); df.to_csv(tmp,index=False,float_format='%.6f'); subprocess.run(['gzip','-1','-f',str(tmp)],check=True)
    stats.append({'file':p.name,'negative_precip_values_fixed':pneg,'negative_shortwave_values_fixed':swneg,'days_shortwave_capped_after_cleaning':int(mask.sum()),'max_pr_mm_h':float(P.max()),'max_rsds_W_m2':float(SW.max())})
    print(p.name,'pneg',pneg,'swneg',swneg,'capdays',int(mask.sum()),flush=True)
pd.DataFrame(stats).to_csv(OUT/'Mohenjo_daro_physical_cleaning_log.csv',index=False)
mp=OUT/'Mohenjo_daro_hourly_disaggregation_metadata.json'; meta=json.load(open(mp))
meta['physical_cleaning']='Tiny negative ERA5-Land de-accumulation artifacts in precipitation and shortwave radiation were clamped to zero. Each affected future day was then renormalized to preserve the exact NEX daily precipitation total or shortwave mean. The 1200 W/m2 shortwave cap was then re-applied with energy redistribution.'
meta['baseline_artifact_note']='ERA5-Land source contained small negative de-accumulated values (baseline raw minima: precipitation about -3.55e-8 m per hour; shortwave about -4 J/m2 per hour). These are numerical/de-accumulation artifacts rather than physical negative precipitation/solar flux.'
with open(mp,'w') as f: json.dump(meta,f,indent=2)
