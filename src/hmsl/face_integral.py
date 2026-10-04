"""Source-integral conductivity mean for homogeneous, pressure-based material.

K_face = integral_{p_L}^{p_R} K(p) dp / (p_R-p_L).
It multiplies the SAME pressure-plus-gravity difference as the legacy scheme.
Source pC(theta) and log10 K(theta) are piecewise linear, so K dp is
integrated analytically on their union of breakpoints, with no re-fitting.
This is a numerical discretisation, not a new measured constitutive law.
RH-based retention and a diffusivity-only source are deliberately not adapted.
"""
from __future__ import annotations
import numpy as np
from .common import InputError

class LiquidFaceIntegral:
    def __init__(self, material):
        m=material
        if m.liquid!='lgKl(Theta_l)' or 'RH' in m.store:
            raise InputError('Source-integral faces need pressure-based pC retention and lgKl(Theta_l). Select the documented harmonic legacy scheme for other source representations; none are silently converted.')
        xy=m.f[m.store]
        if m.store.startswith('pC('): tx,pc=xy
        else: tx,pc=xy[1][::-1],xy[0][::-1]
        kx,lk=m.f[m.liquid]
        self.x=np.unique(np.r_[m.lo,m.hi,tx[(tx>m.lo)&(tx<m.hi)],kx[(kx>m.lo)&(kx<m.hi)]])
        self.pc=np.interp(self.x,tx,pc);self.lk=np.interp(self.x,kx,lk)
        self.a=np.diff(self.pc)/np.diff(self.x)
        self.b=np.diff(self.lk)/np.diff(self.x)
        self.rate=np.log(10.)*(self.a+self.b)
        self.pref=-np.log(10.)*self.a*10.**(self.pc[:-1]+self.lk[:-1])
        total=self._segment(np.arange(len(self.a)),self.x[:-1],self.x[1:])
        if np.any(total<=0) or not np.isfinite(total).all():
            raise InputError('Non-positive/non-finite source conductivity integral. Check source table coverage.')
        self.cum=np.r_[0.,np.cumsum(total)]

    def _segment(self,i,lo,hi):
        dx=hi-lo;z=self.rate[i]*dx
        rel=np.ones_like(z,dtype=float);nz=np.abs(z)>1e-7
        rel[nz]=np.expm1(z[nz])/z[nz]
        rel[~nz]=1.+z[~nz]/2.+z[~nz]**2/6.
        return self.pref[i]*np.exp(self.rate[i]*(lo-self.x[i]))*dx*rel

    def integral(self,left,right):
        """Signed integral K dp; stable even for almost equal cell states."""
        left,right=np.broadcast_arrays(np.asarray(left,float),np.asarray(right,float))
        shape=left.shape;left=left.ravel();right=right.ravel()
        # Same coefficient-evaluation bounds as material.fields; accepted solver
        # states are checked independently and are NOT clipped by this operator.
        left=np.clip(left,self.x[0],self.x[-1]);right=np.clip(right,self.x[0],self.x[-1])
        lo=np.minimum(left,right);hi=np.maximum(left,right)
        il=np.clip(np.searchsorted(self.x,lo,side='right')-1,0,len(self.a)-1)
        ih=np.clip(np.searchsorted(self.x,hi,side='right')-1,0,len(self.a)-1)
        val=np.zeros_like(lo);same=il==ih
        val[same]=self._segment(il[same],lo[same],hi[same])
        cross=~same
        if np.any(cross):
            i=il[cross];j=ih[cross]
            val[cross]=((self.cum[j]-self.cum[i+1])+
                        self._segment(i,lo[cross],self.x[i+1])+
                        self._segment(j,self.x[j],hi[cross]))
        return (np.sign(right-left)*val).reshape(shape)

    def primitive(self,theta):
        return self.integral(np.full_like(np.asarray(theta,float),self.x[0]),theta)

    def face_k(self,theta,p,k):
        theta=np.asarray(theta,float);p=np.asarray(p,float);k=np.asarray(k,float)
        dp=np.diff(p)
        out=(k[:-1]+k[1:])/2.
        use=np.abs(dp)>1e-8*np.maximum(np.maximum(np.abs(p[:-1]),np.abs(p[1:])),1.)
        out[use]=self.integral(theta[:-1][use],theta[1:][use])/dp[use]
        if not np.isfinite(out).all() or np.any(out<=0):
            raise InputError('Invalid source-integral face conductivity; no substituted coefficient was used.')
        return out
