import numpy as np
from scipy.special import loggamma as lg
from rhonp import rho_plus, blk, cproj_plus
from gz import A1, rich
def F11(mu,h,lam,K=4000):
    def part(KK):
        j=np.arange(1,KK+1); a=h*j*lam
        return np.sum(lg(mu+a-(h-1)*lam+1)+lg(mu+a-lam)+lg(-mu+a-h*lam+1)+lg(-mu+a)-lg(-mu+a-(h-1)*lam+1)-lg(-mu+a-lam)-lg(mu+a-h*lam+1)-lg(mu+a))
    return -np.exp(rich(part,K))
def Aplus(mu,n,om):
    T=(n+1)*om/2
    return A1(mu,n,om)*rho_plus(mu,n,om)*blk(T/2,mu,T)**2*blk(T/2-om/2,mu,T)*blk(T/2+om/2,mu,T)
def KB1(mu,n,om):
    h=n+1; T=h*om/2
    return Aplus(mu-T/2+0.5,n,om)*F11(2*mu,h,om)*Aplus(mu+T/2-0.5,n,om)*cproj_plus(mu,om,T)
def target(mu,n,om):
    T=(n+1)*om/2
    return blk(-n*om/2,mu,T)*blk(n*om/2-0.5,mu,T)*blk(0.5,mu,T)
worst=0
for om in (0.37,0.83):
  for n in (1,2,3,4,5):
    for mu in (0.13+0.21j,-0.31+0.1j):
        r=KB1(mu,n,om)/target(mu,n,om); worst=max(worst,abs(r-1)); print(om,n,mu,r)
print('worst',worst)
