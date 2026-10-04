import numpy as np, pandas as pd, json, math
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

L=0.20
D=6.0e-10
N=200  # 1 mm cells, matching published DELPHIN setup
DX=L/N
centers=(np.arange(N)+0.5)*DX

def sorption(phi):
    phi=np.asarray(phi,float)
    return 116.0/(1.0-(1.0/0.118)*np.log(phi))**0.869

wi=float(sorption(0.95))
w_left=float(sorption(0.45))  # exterior at x=0
w_right=float(sorption(0.65)) # interior at x=L

# conservative cell-centred FVM with Dirichlet boundary states at the physical faces.
def rhs(t,w):
    flux=np.empty(N+1)
    flux[0] = -D*(w[0]-w_left)/(DX/2.0)
    flux[1:N] = -D*(w[1:]-w[:-1])/DX
    flux[N] = -D*(w_right-w[-1])/(DX/2.0)
    return -(flux[1:]-flux[:-1])/DX

hours=[100,300,1000]
t_eval=np.array(hours,dtype=float)*3600.0
out=solve_ivp(rhs,(0,t_eval[-1]),np.full(N,wi),method='BDF',t_eval=t_eval,rtol=3e-8,atol=1e-9,max_step=1800.0)
assert out.success, out.message

# Analytical Fourier-series solution to the same benchmark diffusion problem.
# w_ss is the steady linear profile; u=w-w_ss has homogeneous Dirichlet boundaries.
def analytical(x,t,nterms=4000):
    x=np.asarray(x,float)
    wss=w_left+(w_right-w_left)*x/L
    # Initial deviation is a + b*x, integrate analytically for sine coefficients.
    # Bn = 2/L int_0^L [wi-w_left-(w_right-w_left)x/L] sin(n*pi*x/L) dx
    n=np.arange(1,nterms+1,dtype=float)
    k=n*np.pi/L
    a=wi-w_left
    b=-(w_right-w_left)/L
    # integrals: I1=int sin(kx) dx=(1-(-1)^n)/k
    # I2=int x sin(kx) dx= -L*(-1)^n/k
    sign=(-1.0)**n
    I1=(1.0-sign)/k
    I2=-L*sign/k
    B=2.0/L*(a*I1+b*I2)
    decay=np.exp(-D*k*k*t)
    return wss + np.sin(np.outer(x,k))@(B*decay)

rows=[]
profiles=[]
for j,h in enumerate(hours):
    num=out.y[:,j]
    ana=analytical(centers,h*3600.0)
    diff=num-ana
    rmse=float(np.sqrt(np.mean(diff**2)))
    mae=float(np.mean(np.abs(diff)))
    maxabs=float(np.max(np.abs(diff)))
    rng=max(float(np.max(ana)-np.min(ana)),1e-12)
    nrmse=rmse/rng*100
    maxrange=maxabs/rng*100
    # relative error where analytical moisture > 1 kg/m3
    mare=float(np.mean(np.abs(diff)/np.maximum(np.abs(ana),1e-12))*100)
    rows.append(dict(hours=h,rmse_kg_m3=rmse,mae_kg_m3=mae,max_abs_kg_m3=maxabs,nrmse_range_pct=nrmse,max_abs_range_pct=maxrange,mare_pct=mare))
    for x,nv,av in zip(centers,num,ana): profiles.append(dict(hours=h,x_m=x,numerical_w_kg_m3=nv,analytical_w_kg_m3=av,error=nv-av))

summary=pd.DataFrame(rows)
prof=pd.DataFrame(profiles)
summary.to_csv('/mnt/data/q1_revision/benchmark_work/hamstad2_summary.csv',index=False)
prof.to_csv('/mnt/data/q1_revision/benchmark_work/hamstad2_profiles.csv',index=False)
meta={
 'benchmark':'HAMSTAD Benchmark 2 adapter',
 'geometry_m':L,'cells':N,'dx_m':DX,'D_w_m2_s':D,
 'initial_RH':0.95,'exterior_RH':0.45,'interior_RH':0.65,
 'initial_w_kg_m3':wi,'exterior_w_kg_m3':w_left,'interior_w_kg_m3':w_right,
 'solver':'SciPy BDF','rtol':3e-8,'atol':1e-9,'max_step_s':1800,
 'scope':'Benchmark-specific slab/Dirichlet adapter uses the production core finite-volume divergence and BDF integration pattern; it does not test the Mohenjo-daro thickness-averaged distributed exposure boundary.'
}
open('/mnt/data/q1_revision/benchmark_work/hamstad2_meta.json','w').write(json.dumps(meta,indent=2))

fig,ax=plt.subplots(figsize=(7.2,4.6))
for j,h in enumerate(hours):
    num=out.y[:,j]; ana=analytical(centers,h*3600.0)
    ax.plot(centers,num,label=f'{h} h numerical')
    ax.plot(centers,ana,'--',linewidth=1.2,label=f'{h} h analytical')
ax.set_xlabel('Depth from exterior (m)')
ax.set_ylabel('Moisture content (kg m$^{-3}$)')
ax.set_title('HAMSTAD Benchmark 2: numerical core vs analytical solution')
ax.legend(ncol=2,fontsize=8)
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig('/mnt/data/q1_revision/new_figures/fig_hamstad2_benchmark.png',dpi=300)
print(summary.to_string(index=False))
print(json.dumps(meta,indent=2))
