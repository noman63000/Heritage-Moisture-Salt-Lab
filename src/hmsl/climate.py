"""Strict hourly input handling, including CF 365_day and 360_day calendars.

Dates are labels. Simulation time is model-calendar elapsed hours, never a
Gregorian reindex of a no-leap time series. Each row describes one hour [t,t+1).
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import calendar, re
import numpy as np
import pandas as pd
from .common import InputError, sha256

REQUIRED = ['tas_C','hurs_pct','pr_mm_h','rsds_W_m2','rlds_W_m2','sfcWind_m_s']
ALIASES = {
 'tas_C':['temperature_C','T_air_C','temp_C'],
 'hurs_pct':['RH_pct','relative_humidity_pct'],
 'pr_mm_h':['precipitation_mm_h','rain_mm_h'],
 'rsds_W_m2':['shortwave_W_m2','solar_W_m2'],
 'rlds_W_m2':['longwave_W_m2'],
 'sfcWind_m_s':['wind_m_s','wind_speed_m_s'],
 'timestamp_utc':['timestamp','datetime','time']}
SUPPORTED = {'standard','gregorian','proleptic_gregorian','365_day','noleap','360_day'}

@dataclass
class Climate:
    data: pd.DataFrame
    calendar: str
    source_sha256: str
    source_name: str
    warnings: list[str]
    @property
    def n(self): return len(self.data)
    def subset(self, start=0, stop=None):
        return Climate(self.data.iloc[start:stop].reset_index(drop=True), self.calendar,
                       self.source_sha256, self.source_name, list(self.warnings))

def ordinal_hours(dates, hours, cal):
    parts = dates.str.extract(r'^(\d{4})-(\d{2})-(\d{2})$')
    if parts.isna().any().any():
        raise InputError('Dates must use YYYY-MM-DD, with hour_utc in a separate column or an hourly timestamp.')
    y,m,d = [parts[i].astype(int).to_numpy() for i in range(3)]
    h = np.asarray(hours, dtype=float)
    if not np.isfinite(h).all() or np.any((h<0)|(h>23)|(h!=np.floor(h))):
        raise InputError('hour_utc must contain integers from 0 to 23.')
    if np.any((y<1)|(m<1)|(m>12)|(d<1)):
        raise InputError('Invalid calendar date.')
    if cal in ('365_day','noleap','360_day'):
        lengths = np.array([30]*12 if cal=='360_day' else [31,28,31,30,31,30,31,31,30,31,30,31])
        if np.any(d > lengths[m-1]):
            raise InputError(f'Invalid date for {cal}; do not insert leap days or coerce a 360-day date.')
        starts = np.r_[0,np.cumsum(lengths)]
        return ((y*int(lengths.sum())+starts[m-1]+d-1)*24+h).astype(np.int64)
    try:
        t = pd.to_datetime(dates, format='%Y-%m-%d', errors='raise')
        return t.astype('int64').to_numpy()//3_600_000_000_000+h.astype(np.int64)
    except Exception as e:
        raise InputError(f'Invalid Gregorian date: {e}') from e

def read_climate(path, calendar_name='gregorian', column_map=None, expected_sha256=None):
    p=Path(path)
    if not p.is_file(): raise InputError(f'Missing climate file: {p.name}. Import it in the Data tab.')
    cal=calendar_name.lower()
    if cal not in SUPPORTED: raise InputError('Calendar must be Gregorian, standard, proleptic_gregorian, 365_day, noleap or 360_day.')
    actual=sha256(p)
    if expected_sha256 and actual!=expected_sha256:
        raise InputError(f'Checksum mismatch: {p.name}. This is not the forcing file recorded in the frozen manifest.')
    if p.suffix.lower()=='.xlsx':
        book=pd.ExcelFile(p, engine='openpyxl')
        if 'Climate' not in book.sheet_names:
            raise InputError("An input workbook needs a sheet named Climate, with row-1 column headings. The large audit workbook is not a climate series; use the supplied CSV.gz files.")
        df=pd.read_excel(book,sheet_name='Climate')
    else:
        try: df=pd.read_csv(p,compression='infer')
        except Exception as e: raise InputError(f'Could not read CSV/CSV.gz: {e}') from e
    df.columns=df.columns.astype(str).str.strip()
    if column_map: df=df.rename(columns=column_map)
    for target, alternatives in ALIASES.items():
        matches=[s for s in alternatives if s in df and target not in df]
        if len(matches)>1: raise InputError(f'Ambiguous aliases for {target}; rename the correct column.')
        if matches: df=df.rename(columns={matches[0]:target})
    missing=[s for s in REQUIRED if s not in df]
    if missing: raise InputError('Missing climate columns: '+', '.join(missing)+'. Units are explicit; Kelvin, fractional RH and daily rain are not silently converted.')
    if 'date' in df and 'hour_utc' in df:
        dates=df['date'].astype(str).str.slice(0,10)
        hours=pd.to_numeric(df['hour_utc'],errors='coerce')
    elif 'timestamp_utc' in df:
        stamps=df['timestamp_utc'].astype(str)
        parts=stamps.str.extract(r'^(\d{4}-\d{2}-\d{2})[ T](\d{2}):00(?::00)?(?:Z|\+00:00)?$')
        if parts.isna().any().any():
            raise InputError('Timestamps must be exact hourly labels, such as 2010-01-01 00:00:00. Convert local time to UTC explicitly first.')
        dates=parts[0]; hours=parts[1].astype(int)
    else: raise InputError('Provide timestamp_utc, or both date and hour_utc.')
    if len(df)<2: raise InputError('At least two hourly rows are required.')
    times=ordinal_hours(dates,hours,cal)
    if not np.all(np.diff(times)==1):
        bad=int(np.flatnonzero(np.diff(times)!=1)[0])+2
        raise InputError(f'Non-consecutive, duplicate or unordered hour at data row {bad}. No gap filling has been applied. Check the selected calendar.')
    for k in REQUIRED+['capillary_supply_mm_day','sulfate_source_mg_L']:
        if k not in df: continue
        df[k]=pd.to_numeric(df[k],errors='coerce')
        if not np.isfinite(df[k]).all(): raise InputError(f'{k} contains missing or non-numeric values.')
    if ((df.hurs_pct<0)|(df.hurs_pct>100)).any(): raise InputError('RH must be in percent, within 0 to 100.')
    if ((df.tas_C < -70)|(df.tas_C > 80)).any(): raise InputError('Temperature outside -70 to 80 C. Check units and source.')
    for k in REQUIRED[2:]+['capillary_supply_mm_day','sulfate_source_mg_L']:
        if k in df and (df[k]<0).any(): raise InputError(f'{k} contains negative values. Correct the source explicitly; values are not clipped.')
    df['date']=dates.to_numpy(); df['hour_utc']=hours.astype(int).to_numpy()
    df['timestamp_utc']=df['date']+' '+df['hour_utc'].map(lambda x:f'{x:02d}:00:00')
    warnings=[]
    if df.hurs_pct.max()<=1: warnings.append('All RH values are <=1%. Verify that fractional RH was not supplied.')
    if 'capillary_supply_mm_day' in df or 'sulfate_source_mg_L' in df:
        warnings.append('Time-varying boundary columns are present; these override the reviewed constant boundary values.')
    if (df.tas_C<0).any(): warnings.append('Air temperature is below freezing. The transport solver has no ice model and will refuse this case; climate screening remains available.')
    return Climate(df.reset_index(drop=True),cal,actual,p.name,warnings)

def climate_screen(c: Climate):
    d=c.data
    daily=d.groupby('date',sort=False).agg(temperature_mean_C=('tas_C','mean'),temperature_max_C=('tas_C','max'),
        temperature_min_C=('tas_C','min'),RH_mean_pct=('hurs_pct','mean'),rain_mm=('pr_mm_h','sum'),
        shortwave_mean_W_m2=('rsds_W_m2','mean'),longwave_mean_W_m2=('rlds_W_m2','mean'),wind_mean_m_s=('sfcWind_m_s','mean'),hours=('tas_C','size')).reset_index()
    # This is an explicitly user-independent screening convention, NOT a salt equilibrium curve.
    dry=d.hurs_pct.to_numpy()<75.0
    crossings=int(np.sum(dry[1:]!=dry[:-1]))
    summary={'mode':'CLIMATE SCREENING ONLY - not a masonry simulation','hours':c.n,'calendar':c.calendar,
      'start':str(d.timestamp_utc.iloc[0]),'end':str(d.timestamp_utc.iloc[-1]),
      'temperature_mean_C':float(d.tas_C.mean()),'temperature_min_C':float(d.tas_C.min()),'temperature_max_C':float(d.tas_C.max()),
      'RH_mean_pct':float(d.hurs_pct.mean()),'precipitation_total_mm':float(d.pr_mm_h.sum()),
      'precipitation_mean_per_model_day_mm':float(d.pr_mm_h.sum()/(c.n/24)),
      'RH75_crossings_screening_only':crossings,'hours_RH_ge_80_pct':int((d.hurs_pct>=80).sum()),
      'hours_air_T_ge_40_C':int((d.tas_C>=40).sum()),'hours_air_T_below_0_C':int((d.tas_C<0).sum()),
      'source_sha256':c.source_sha256,'warnings':c.warnings,
      'interpretation':'75% and 80% RH and 40 C are transparent screening cutoffs, not calibrated deterioration thresholds. Crossings are not crystallisation events.'}
    return summary,daily
