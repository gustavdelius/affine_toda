import numpy as np, time, sys
from sol_boot import boot_coef
from rhonp import rho_plus, blk
om=float(sys.argv[2]) if len(sys.argv)>2 else 0.37
n=int(sys.argv[1]); T=(n+1)*om/2
sig=lambda m: blk(T/2,m,T)**2*blk(T/2-om/2,m,T)*blk(T/2+om/2,m,T)
for th in (0.31+0.07j,-0.2+0.15j):
    t=time.time()
    cd,rd=boot_coef(n,om,th,0,order=-1); cp,rp=boot_coef(n,om,th,1,order=-1)
    P=cp/cd; mu=-1j*T*th/np.pi; mus=[mu+T-k*om for k in range(1,n+1)]
    Y=lambda m: rho_plus(m,n,om)*sig(m)
    print(n,om,th,'res',f'{max(rd,rp):.1e}','  rho sigma: Y/(P prod Y) =',Y(mu)/(P*np.prod([Y(m) for m in mus])),
          ' sigma/prod sigma =',sig(mu)/np.prod([sig(m) for m in mus]), f'{time.time()-t:.0f}s',flush=True)
