import numpy as np, sys
from rhom import rho_minus, cproj_minus
from rhonp import blk
from sol_boot import boot_coef
om=0.37
for n in (1,2,3,4,5):
    T=(n+1)*om/2
    R=lambda mu: rho_minus(mu-T/2+0.5,n,om)*rho_minus(mu+T/2-0.5,n,om)*cproj_minus(mu,om,T)
    print(n,'R_-/((om/2)(T-om/2)) =',[np.round(R(mu)/(blk(om/2,mu,T)*blk(T-om/2,mu,T)),8) for mu in (0.13+0.21j,-0.3+0.1j)])
for n in (2,3,4):
    T=(n+1)*om/2; th=0.31+0.07j; mu=-1j*T*th/np.pi
    P=boot_coef(n,om,th,-1,order=-1)[0]/boot_coef(n,om,th,0,order=-1)[0]
    Y=lambda m: rho_minus(m,n,om)
    mus=[mu+T-k*om for k in range(1,n+1)]
    print(n,'s=-1 bootstrap Y/(P prod Y) with bare rho_-:',Y(mu)/(P*np.prod([Y(m) for m in mus])))
