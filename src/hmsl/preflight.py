"""Full-case compatibility checks. Valid weather is not necessarily model-compatible."""
from .common import InputError,Cancelled,atomic_json,now

def check_cases(project,ids=None,save=True,progress=None,cancel=None):
    ids=ids or [r['id'] for r in project.config['runs']]
    rows=[]
    for index,rid in enumerate(ids):
        if cancel and cancel():raise Cancelled('Compatibility check cancelled. No transport case was started.')
        if progress:progress(index/len(ids),f'Compatibility check {index+1}/{len(ids)}: reading {rid}')
        try:
            c=project.climate(rid);d=c.data;mask=d.tas_C<0
            row={'run_id':rid,'climate_sha256':c.source_sha256,'hours':c.n,
                 'min_air_temperature_C':float(d.tas_C.min()),'subzero_air_hours':int(mask.sum()),
                 'full_case_supported':not bool(mask.any()),
                 'first_subzero_timestamp':str(d.loc[mask,'timestamp_utc'].iloc[0]) if mask.any() else None,
                 'reason':'Ice/freezing physics is outside this reduced model; temperatures were not altered.' if mask.any() else 'No sub-zero forcing detected. Runtime checks and QA are still required.'}
        except Exception as e:row={'run_id':rid,'full_case_supported':False,'reason':str(e)}
        rows.append(row)
    result={'created_utc':now(),'all_selected_supported':all(x['full_case_supported'] for x in rows),'cases':rows,
      'scope':'Preflight compatibility, NOT scientific validation. Negative air temperature is the existing conservative ice-free model scope guard. A supported case still needs QA, spin-up and checks during integration.',
      'instruction':'Do not delete sub-zero hours, replace temperatures with zero, or silently omit GCMs from the planned ensemble. Unsupported cases need a documented change of model scope or a freeze-thaw-capable model.'}
    if save:atomic_json(project.path/'preflight.json',result)
    return result

def require_supported(project,ids,progress=None,cancel=None):
    result=check_cases(project,ids,progress=progress,cancel=cancel)
    bad=[x['run_id'] for x in result['cases'] if not x['full_case_supported']]
    if bad:
        raise InputError('Full-case preflight stopped before starting any selected simulation. Unsupported cases: '+', '.join(bad)+'. Open the compatibility report in Run & check. Valid sub-zero climate is outside this ice-free solver; no temperatures were clipped.')
    return result
