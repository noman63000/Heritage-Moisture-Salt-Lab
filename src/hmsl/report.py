from __future__ import annotations
import html,json
from pathlib import Path
import pandas as pd
from .common import atomic_json

def sparkline(values):
    v=[float(x) for x in values if pd.notna(x)]
    if not v:return ''
    step=max(1,len(v)//900); v=v[::step]
    lo,hi=min(v),max(v); scale=hi-lo if hi>lo else 1
    pts=' '.join(f'{30+i*840/max(len(v)-1,1):.1f},{160-(x-lo)/scale*130:.1f}' for i,x in enumerate(v))
    return f'<svg role="img" aria-label="Series from first to last record" viewBox="0 0 900 190"><path d="M30 15 V160 H870" fill="none" stroke="currentColor" opacity=".25"/><polyline points="{pts}" fill="none" stroke="#06776c" stroke-width="2"/><text x="32" y="180" font-size="12">{lo:.3g} to {hi:.3g}; chronological records</text></svg>'

def write_report(folder,summary,daily,provenance):
    p=Path(folder); atomic_json(p/'summary.json',summary); atomic_json(p/'provenance.json',provenance)
    daily.to_csv(p/'daily_summary.csv',index=False)
    display={k:v for k,v in summary.items() if not isinstance(v,(dict,list))}
    table=''.join('<tr><th>'+html.escape(k.replace('_',' '))+'</th><td>'+html.escape(str(v))+'</td></tr>' for k,v in display.items())
    curves=[]
    for k in ['temperature_mean_C','rain_mm','RH_mean_pct','theta_mean_m3_m3','T_mean_C','sulfate_cmax_0p10m_bandavg_mol_m3_liquid','sulfate_cmax_mol_m3_liquid']:
        if k in daily:curves.append('<h2>'+html.escape(k.replace('_',' '))+'</h2>'+sparkline(daily[k]))
    doc='''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>Heritage Moisture &amp; Salt Lab - report</title><style>body{font:16px/1.55 system-ui,sans-serif;max-width:1050px;margin:40px auto;padding:0 24px;color:#16323b;background:#f5f8f8}h1{font-size:32px}h2{margin-top:32px}table{width:100%;border-collapse:collapse;background:white}th,td{text-align:left;padding:10px 14px;border-bottom:1px solid #dce5e5}th{width:50%}.notice{padding:16px;background:#fff3da;border-left:4px solid #c58315}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:white;padding:20px}svg{width:100%;background:white;border-radius:12px}</style></head><body><p>HERITAGE MOISTURE &amp; SALT LAB / 1.0.2</p><h1>'''+html.escape(provenance.get('site',''))+'''</h1><div class="notice">Research scenario output. Climate screening, transport-only output, optional equilibrium partition and synthetic demonstrations are distinct. None is a calibrated damage prediction. Check the mode, input provenance and QA scope below.</div><h2>Summary</h2><table>'''+table+'</table>'+''.join(curves)+'<h2>Provenance and limitations</h2><pre>'+html.escape(json.dumps(provenance,indent=2))+'</pre></body></html>'
    (p/'report.html').write_text(doc,encoding='utf-8')
