"""Policy-specific audit without editing prior chunks or declaring ice absent as fact."""
from pathlib import Path
import pandas as pd
import numpy as np
from html import escape
from .safety import atomic,read_json
from .pore_ice import METHOD,MAX_ICE_SATURATION,MIN_SCREEN_C,SOURCES

def summarize_scope(folder,resume_hour=0):
    folder=Path(folder);h=pd.read_csv(folder/'hourly_results.csv.gz')
    col='cold_scope_subzero_material'
    available=h[col].notna() if col in h else pd.Series(False,index=h.index)
    checked=h.loc[available]
    # Legacy records may only precede newly screened records and must already
    # belong to the resumed prefix. A missing fresh audit is not a success.
    first=int(np.flatnonzero(available.to_numpy())[0]) if available.any() else len(h)
    if first > int(resume_hour) or not bool(available.iloc[first:].all()):
        raise ValueError('Missing or noncontiguous new scope columns; completion refused.')
    required=['cold_scope_min_internal_T_C','cold_scope_subzero_material',
              'cold_scope_max_pore_ice_saturation','cold_scope_max_ambient_frost_saturation',
              'cold_scope_sampled_state_vectors']
    if len(checked):
        if any(c not in checked for c in required) or not np.isfinite(checked[required].to_numpy(dtype=float)).all():
            raise ValueError('Incomplete/nonfinite cold-scope audit; completion refused.')
        if not checked[col].isin([0,1]).all() or (checked.cold_scope_sampled_state_vectors < 1).any():
            raise ValueError('Invalid cold-scope state counts; completion refused.')
        if (checked.cold_scope_min_internal_T_C < MIN_SCREEN_C).any():
            raise ValueError('A saved state is outside the cold-scope domain; completion refused.')
    record={'policy':METHOD,'passed_under_screen':True,'interpretation':'Numerical no-ice admissibility under a retained single-branch moisture model; not measured freezing or field validation.',
       'hours_total':len(h),'hours_with_new_screen_columns':int(available.sum()),
       'legacy_hours_without_new_screen_columns':int((~available).sum()),
       'legacy_note':'Old chunks were retained, not rechecked in this run. Only checkpoints from the verified original temperature-guard core are reused. Old completion does not become cold-physics validation.',
       'subzero_material_hours_in_checked_part':int(checked[col].sum()) if len(checked) else 0,
       'minimum_internal_material_temperature_C_in_checked_part':float(checked.cold_scope_min_internal_T_C.min()) if len(checked) else None,
       'maximum_pore_ice_saturation_in_checked_part':float(checked.cold_scope_max_pore_ice_saturation.max()) if len(checked) else None,
       'maximum_ambient_frost_saturation_in_checked_part':float(checked.cold_scope_max_ambient_frost_saturation.max()) if len(checked) else None,
       'relative_ice_saturation_limit':MAX_ICE_SATURATION,'minimum_allowed_screen_temperature_C':MIN_SCREEN_C,
       'salt_freezing_depression_used':False,'resume_hour_this_execution':resume_hour,'sources':SOURCES}
    if len(checked) and (checked.cold_scope_max_pore_ice_saturation.max()>=MAX_ICE_SATURATION or checked.cold_scope_max_ambient_frost_saturation.max()>=MAX_ICE_SATURATION):
        raise ValueError('The scope audit found a non-admissible saved cold state. Do not use completion for research.')
    atomic(folder/'cold_scope_summary.json',record)
    labels=[('Total hours',record['hours_total']),('Newly screened hours',record['hours_with_new_screen_columns']),
      ('Retained legacy hours (not retrospectively rescreened)',record['legacy_hours_without_new_screen_columns']),
      ('Subzero-material hours in screened part',record['subzero_material_hours_in_checked_part']),
      ('Minimum sampled material temperature (C)',record['minimum_internal_material_temperature_C_in_checked_part']),
      ('Maximum pore ice-saturation ratio',record['maximum_pore_ice_saturation_in_checked_part']),
      ('Maximum ambient frost-saturation ratio',record['maximum_ambient_frost_saturation_in_checked_part'])]
    rows=''.join('<tr><td>'+escape(str(k))+'</td><td>'+escape(str(v))+'</td></tr>' for k,v in labels)
    links=''.join('<p><a href="'+escape(u,quote=True)+'">'+escape(u)+'</a></p>' for u in SOURCES)
    text=('<!doctype html><meta charset="utf-8"><title>Cold-scope audit</title>'
      '<style>body{max-width:950px;margin:35px auto;padding:20px;font:16px/1.55 system-ui;color:#183b42}'
      'td{padding:9px;border-bottom:1px solid #ccc}.note{padding:18px;background:#fff1d9}</style>'
      '<h1>Revised pore-liquid scope: completed under the acceptance screen</h1>'
      '<p class="note">This is not a freeze-thaw simulation or measured proof that the real brick stayed ice-free. '
      'The transport equations remain the original reduced model. The previous temperature-only restriction was replaced '
      'by a separately reviewed temperature/moisture stability screen. Legacy warm prefixes were retained, not rechecked.</p>'
      '<table>'+rows+'</table><h2>Screening assumptions</h2>'
      '<p>For each subzero material state, pore RH from the source retention curve is divided by the ice/liquid saturation '
      'vapour-pressure ratio (Murphy and Koop, equations 7 and 10). Both pore ice saturation and ambient frost saturation '
      'at a cold exposed cell must remain below 0.99. The 1% buffer is a conservative screening choice, not a confidence '
      'level or measured brick parameter. Cold-screen domain: -20 C or warmer. No salt freezing-point credit is used.</p>'
      '<p>Accepted BDF states are checked, with dense-output samples at no more than 60-second spacing near zero. '
      'This is sampled numerical screening, not an exact continuous-extremum proof. Unknown precipitation phase in '
      'subzero air and rain on a cold exposed surface are not accepted. No ice fraction, latent freezing heat, frost '
      'deposition, salt activity correction, or mechanical damage is calculated. Source retention transferred to cold '
      'conditions is an assumption requiring review.</p><p><a href="cold_scope_summary.json">Machine-readable audit</a> | '
      '<a href="cold_scope_policy.json">Approval and method identity</a></p><h2>Sources</h2>'+links)
    (folder/'cold_scope_report.html').write_text(text,encoding='utf-8')
    main=folder/'report.html'
    if main.is_file():
        text=main.read_text(encoding='utf-8')
        banner=('<aside style="padding:18px;border:2px solid #987125;background:#fff1d9">'
          '<b>Revised cold-scope acceptance policy.</b> This run uses temperature plus pore humidity, not a full freezing model. '
          '<a href="cold_scope_report.html">Read the scope audit and limitations</a>.</aside>')
        end=text.lower().rfind('</body>')
        text=text[:end]+banner+text[end:] if end>=0 else text+banner
        main.write_text(text,encoding='utf-8')
    return record
