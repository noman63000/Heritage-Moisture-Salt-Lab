"""Exact arithmetic cache for values independent of the evolving state."""
import inspect,textwrap
import numpy as np
from hmsl.solver import ColumnModel as Reference, RV, SIGMA, psat

src=textwrap.dedent(inspect.getsource(Reference.derivatives))
replacements=[
 ("Ta,RHa,rain,sw,lw,wind,supply,csource=forcing[i]", "Ta,RHa,rain,sw,lw,wind,supply,csource=forcing[i]\n    con=self._constant_forcing(forcing[i])"),
 ("pv=rh*psat(T); pvair=RHa/100.*psat(Ta)","pv=rh*psat(T); pvair=con['pvair']"),
 ("exposure=np.full(self.n,self.a)\n    exposure[-1]+=1./self.widths[-1]", "exposure=con['exposure'].copy()"),
 ("hc=5.7+3.8*wind\n    beta=hc/(1.2*1005.*RV*(Ta+273.15))", "hc=con['hc']\n    beta=con['beta']"),
 ("rainfall=np.full(self.n,self.a*self.surface['rain_catch_fraction']*rain/1000./3600.)\n    if self.surface['crown_rain']: rainfall[-1]+=rain/1000./3600./self.widths[-1]", "rainfall=con['rainfall']"),
 ("sky=self.surface['sky_view_factor']\n    effective_lw=sky*lw+(1-sky)*SIGMA*(Ta+273.15)**4\n    absorbed=self.surface['absorptivity']*self.surface['solar_projection_factor']*sw", "effective_lw=con['effective_lw']\n    absorbed=con['absorbed']"),
]
for before,after in replacements:
    if src.count(before)!=1: raise RuntimeError('Reference derivatives no longer match tested performance adapter.')
    src=src.replace(before,after)
namespace=Reference.derivatives.__globals__.copy()
exec(compile(src,'<exact_constant_cache>','exec'),namespace)
class CachedModel(Reference):
    derivatives=namespace['derivatives']
    def _constant_forcing(self,row):
        key=row.tobytes()
        if getattr(self,'_forcing_key',None)==key:return self._forcing_constants
        Ta,RHa,rain,sw,lw,wind,supply,csource=row
        pvair=RHa/100.*psat(Ta)
        exposure=np.full(self.n,self.a);exposure[-1]+=1./self.widths[-1]
        hc=5.7+3.8*wind
        beta=hc/(1.2*1005.*RV*(Ta+273.15))
        rainfall=np.full(self.n,self.a*self.surface['rain_catch_fraction']*rain/1000./3600.)
        if self.surface['crown_rain']:rainfall[-1]+=rain/1000./3600./self.widths[-1]
        sky=self.surface['sky_view_factor']
        effective_lw=sky*lw+(1-sky)*SIGMA*(Ta+273.15)**4
        absorbed=self.surface['absorptivity']*self.surface['solar_projection_factor']*sw
        self._forcing_key=key
        self._forcing_constants={'pvair':pvair,'exposure':exposure,'hc':hc,'beta':beta,'rainfall':rainfall,'effective_lw':effective_lw,'absorbed':absorbed}
        return self._forcing_constants
