import numpy as np
from scipy.special import loggamma as lg
from rhonp import rho_plus, blk
from rhom import rho_minus
from gz import rich
def a0(mu,n,om,s=1,K=4000):
    T=(n+1)*om/2; sh=0 if s==1 else 0.5
    G=lambda j,a: lg(mu+2*j*T+a)-lg(-mu+2*j*T+a)
    def part(KK):
        j=np.arange(1,KK+1)
        return np.sum(G(j,om/2-1.5*T+sh)+G(j,1-om/2-T/2+sh)-G(j,om/2-T/2+sh)-G(j,1-om/2-1.5*T+sh))
    pre=np.sin(np.pi*om/2)/(np.sin(np.pi*(mu-T/2)) if s==1 else np.cos(np.pi*(mu-T/2)))
    return pre*np.exp(rich(part,K))
for om in (0.37,0.9):
  for n in (2,3,4,5):
    T=(n+1)*om/2
    for mu in (0.13+0.21j,-0.3+0.1j):
        rp=rho_plus(mu,n,om)*blk(T/2,mu,T)**2*blk(-T/2,mu,T)   # A^(+)/(sigma_G a1): rho (T/2)^2 (-T/2) = rho (T/2)
        rm=-rho_minus(mu,n,om)*blk(-T/2,mu,T)
        print(om,n,mu,' s=+1 ratio',np.round(rp/a0(mu,n,om,1),8),'  s=-1 ratio',np.round(rm/a0(mu,n,om,-1),8))
