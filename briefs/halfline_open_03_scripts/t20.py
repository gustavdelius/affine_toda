import numpy as np
from rhom import rho_minus, cproj_minus
from fitR import fitfun
om=0.37
for n in (1,2,3,4):
    T=(n+1)*om/2
    R=lambda mu: rho_minus(mu-T/2+0.5,n,om)*rho_minus(mu+T/2-0.5,n,om)*cproj_minus(mu,om,T)
    d,z,p,sv=fitfun(R,T)
    print(n,'T',round(T,4),'zeros',z,'poles',p, sv)
    print('   const check: R/prod',)
