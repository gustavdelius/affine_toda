import numpy as np
from scipy.special import loggamma as lg
from rhonp import rho_plus, blk
from gz import A1, rich
lam=0.37; n=2; T=1.5*lam
def G(j,m,a): return lg(m+3*j*lam+a)-lg(-m+3*j*lam+a)
def A_G(mu,K=4000):
    def part(KK):
        j=np.arange(1,KK+1)
        s1=np.sum(G(j,mu,-7*lam/4)+G(j,mu,-5*lam/4+1)-G(j,mu,-lam/4)-G(j,mu,-11*lam/4+1))
        s2=np.sum(G(2*j,2*mu,-1.5*lam)+G(2*j,2*mu,-4.5*lam+1)-G(2*j,2*mu,-2.5*lam)-G(2*j,2*mu,-3.5*lam+1))
        return s1+s2
    return np.sin(np.pi*lam/2)/np.sin(np.pi*(mu-0.75*lam))*np.exp(rich(part,K))
def mine(mu):
    sig=blk(T/2,mu,T)**2*blk(T/2-lam/2,mu,T)*blk(T/2+lam/2,mu,T)
    return A1(mu,n,lam)*rho_plus(mu,n,lam)*sig
for mu in (0.13+0.21j,-0.31+0.1j,0.05-0.3j):
    sp=blk(lam/4,mu,T)*blk(5*lam/4,mu,T)
    print('A_G^+ sigma^+ / A^(+) =',A_G(mu)*sp/mine(mu), '   with mu->-mu:',A_G(-mu)*blk(lam/4,-mu,T)*blk(5*lam/4,-mu,T)/mine(mu))
