import numpy as np
from scipy.special import loggamma as lg
from rhonp import blk
from rhom import rho_minus
from gz import rich
def a0m(mu,n,om,K=4000):
    T=(n+1)*om/2
    G=lambda j,a: lg(mu+2*j*T+a)-lg(-mu+2*j*T+a)
    def part(KK):
        j=np.arange(1,KK+1)
        return np.sum(G(j,om/2-1.5*T+0.5)+G(j,0.5-om/2-T/2)-G(j,om/2-T/2+0.5)-G(j,0.5-om/2-1.5*T))
    return np.sin(np.pi*om/2)/np.cos(np.pi*(mu-T/2))*np.exp(rich(part,K))
for om in (0.37,0.9):
  for n in (2,3,4,5):
    T=(n+1)*om/2
    out=[]
    for mu in (0.13+0.21j,-0.3+0.1j):
        r=-blk(-T/2,mu,T)*rho_minus(mu,n,om)/a0m(mu,n,om)
        out.append((np.round(r,6), np.round(r/blk(T/2,mu,T),6), np.round(r*blk(T/2,mu,T),6)))
    print(om,n,out)
