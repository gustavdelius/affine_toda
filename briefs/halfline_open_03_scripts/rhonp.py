import numpy as np
from scipy.special import loggamma as lg
def lpsi(c,mu): return lg(c+mu)-lg(c-mu)
def log_rho_o(mu,a,b,T,K):
    k=np.arange(2*K); c=(2*k+1)*a; sg=(-1.0)**k
    return np.sum(sg*(lpsi(b+c,mu)+lpsi(1-T+c,mu)-lpsi(1-b+c,mu)-lpsi(T+c,mu)))
def rho_plus(mu,n,om,K=4000):
    h=n+1; T=h*om/2; a=T/2; b=om/2
    lo=2*log_rho_o(mu,a,b,T,2*K)-log_rho_o(mu,a,b,T,K)
    return np.sin(np.pi*b)/np.sin(np.pi*(a+mu))*np.exp(lo)
def blk(x,mu,T): return np.sin(np.pi*(mu+x)/(2*T))/np.sin(np.pi*(mu-x)/(2*T))
def cproj_plus(mu,om,T):
    return -np.cos(np.pi*mu)*np.cos(np.pi*(mu+om/2))*np.cos(np.pi*(mu-T))/(np.sin(np.pi*om/2)**2*np.cos(np.pi*(mu-om/2)))
def Rplus(mu,n,om,rho=rho_plus):
    T=(n+1)*om/2
    return rho(mu-T/2+0.5,n,om)*rho(mu+T/2-0.5,n,om)*cproj_plus(mu,om,T)
def target_plus(mu,n,om):
    h=n+1; T=h*om/2
    return blk(n*om/2-0.5,mu,T)*blk(0.5,mu,T)*blk(om/2+0.5,mu,T)/blk(0.5+h*om/2,mu,T)
if __name__=='__main__':
    import time; t=time.time()
    om=0.37
    for n in (1,2):
        mu=0.13+0.21j; print(n, Rplus(mu,n,om)/target_plus(mu,n,om))
    print(time.time()-t)
