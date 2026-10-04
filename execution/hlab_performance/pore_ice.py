"""No-ice admissibility screen, not an ice/freeze-thaw transport solver.

A source retention/Kelvin pore relative humidity is compared with the ratio of
saturation vapour pressures over planar ice and supercooled pure liquid water.
Murphy & Koop (2005), QJRMS 131:1539-1565, DOI 10.1256/qj.04.94, eqs. 7 and 10.
The freezing-limit-potential interpretation is described by Kuenzel (1994) and
WUFI's technical explanation (2025/2026). No salt freezing-point credit is used.

All three safeguards are mandatory: dry/unfrozen-pore criterion, no ambient
frost deposition potential on cold exposed cells, and no cold precipitation of
unresolved phase. A 1% relative-saturation buffer is a conservative *screening*
choice, not a measured material property or damage threshold.
"""
from __future__ import annotations
import numpy as np
from hmsl.material import RW, RV
from hmsl.solver import psat

METHOD = 'pore-liquid-scope-1.0'
MAX_ICE_SATURATION = 0.99
MIN_SCREEN_C = -20.0  # Declared mild-cold applicability limit; no extrapolation beyond it.
SOURCES = [
 'https://doi.org/10.1256/qj.04.94',
 'https://www.wufi-forum.com/viewtopic.php?t=2278',
 'https://doi.org/10.1029/2018WR023221',
 'https://doi.org/10.1016/j.buildenv.2024.111399',
]

def log_e_ice(T_K):
    T = np.asarray(T_K, dtype=float)
    if not np.isfinite(T).all() or np.any(T <= 110.0):
        raise ValueError('Ice vapour-pressure temperature is outside the source range.')
    return 9.550426 - 5723.265/T + 3.53068*np.log(T) - 0.00728332*T

def log_e_liquid(T_K):
    T = np.asarray(T_K, dtype=float)
    if not np.isfinite(T).all() or np.any(T <= 123.0) or np.any(T >= 332.0):
        raise ValueError('Liquid vapour-pressure temperature is outside 123 < T < 332 K.')
    return (54.842763 - 6763.22/T - 4.210*np.log(T) + 0.000367*T
        + np.tanh(0.0415*(T-218.8))*(53.878-1331.22/T-9.44523*np.log(T)+0.014025*T))

def freezing_limit_rh(T_C):
    """Equilibrium pore RH threshold over liquid, not outside-air humidity."""
    T = np.asarray(T_C, dtype=float) + 273.15
    return np.exp(log_e_ice(T)-log_e_liquid(T))

def assess(material, theta, temperature_C, forcing, exposed):
    """Return a conservative local admissibility record; never mutate states.

    theta and temperature: one same-sized vector of cell states.
    forcing: [air C, air RH %, rain mm/h, SW, LW, wind, supply, source SO4].
    For screening a cold exposed cell, its temperature is the *existing reduced
    model* surface-temperature approximation; no unmodelled face depth is claimed.
    """
    th=np.asarray(theta,dtype=float);T=np.asarray(temperature_C,dtype=float)
    ex=np.broadcast_to(np.asarray(exposed,dtype=bool),T.shape)
    if th.shape!=T.shape or not np.isfinite(th).all() or not np.isfinite(T).all():
        raise ValueError('Non-finite or mismatched state in cold-scope check.')
    if np.any(th<material.lo-1e-9) or np.any(th>material.hi+1e-9):
        raise ValueError('A cold-scope check cannot certify an out-of-range water state.')
    Ta, RH, rain = map(float,forcing[:3])
    if not all(np.isfinite([Ta,RH,rain])) or not 0<=RH<=100 or rain<0:
        raise ValueError('Invalid external forcing in cold-scope check.')
    cold=T<0.0
    out={'method':METHOD,'allowed':True,'cold_cells':int(cold.sum()),
         'minimum_material_temperature_C':float(T.min()),'maximum_ice_saturation':0.0,
         'maximum_ambient_frost_saturation':0.0,'maximum_pore_RH_cold':None,
         'minimum_freezing_limit_RH':None,'reason':None,'cell':None,
         'ice_saturation_acceptance_limit':MAX_ICE_SATURATION}
    if Ta<=0.0 and rain>0.0:
        out.update(allowed=False,reason='COLD_PRECIPITATION_PHASE_UNRESOLVED')
        return out
    if not cold.any():return out
    if T[cold].min()<MIN_SCREEN_C:
        out.update(allowed=False,reason='COLD_SCREEN_TEMPERATURE_RANGE',cell=int(np.argmin(T)))
        return out
    idx=np.flatnonzero(cold);tc=T[cold];thc=th[cold]
    p=material.pressure(thc,tc)
    rhp=np.exp(np.maximum(p/(RW*RV*(tc+273.15)),-700))
    rhl=freezing_limit_rh(tc)
    ice_sat=rhp/rhl
    # The original equation's ambient vapour pressure is retained. Also check a
    # Murphy-Koop ambient reference where covered, taking the larger (safer) one.
    pvair=RH/100.0*float(psat(Ta))
    if 123.0<Ta+273.15<332.0:
        pvair=max(pvair,RH/100.0*float(np.exp(log_e_liquid(Ta+273.15))))
    frost=pvair/np.exp(log_e_ice(tc+273.15))
    out.update(maximum_ice_saturation=float(ice_sat.max()),
               maximum_ambient_frost_saturation=float(frost[ex[cold]].max()) if ex[cold].any() else 0.0,
               maximum_pore_RH_cold=float(rhp.max()),minimum_freezing_limit_RH=float(rhl.min()))
    bad=np.flatnonzero(ice_sat>=MAX_ICE_SATURATION)
    if len(bad):out.update(allowed=False,reason=('PORE_ICE_POSSIBLE' if ice_sat[bad[0]]>=1.0 else 'PORE_ICE_MARGIN_UNCERTAIN'),cell=int(idx[bad[0]]))
    elif np.any(ex[cold] & (frost>=MAX_ICE_SATURATION)):
        bad=np.flatnonzero(ex[cold] & (frost>=MAX_ICE_SATURATION))
        out.update(allowed=False,reason=('SURFACE_FROST_POSSIBLE' if frost[bad[0]]>=1.0 else 'SURFACE_FROST_MARGIN_UNCERTAIN'),cell=int(idx[bad[0]]))
    elif rain>0 and ex[cold].any():
        out.update(allowed=False,reason='RAIN_ON_COLD_SURFACE',cell=int(idx[np.flatnonzero(ex[cold])[0]]))
    if out['cell'] is not None:
        j=out['cell'];out.update(material_temperature_C=float(T[j]),theta_m3_m3=float(th[j]))
    return out
