import numpy as np
from rhonp import Rplus, target_plus
def fitfun(F,T,dmax=12,seed=2):
    rng=np.random.default_rng(seed); m=6*dmax
    mus=rng.uniform(-T,T,size=m)+1j*rng.uniform(-0.3,0.3,size=m)*T
    taus=np.exp(1j*np.pi*mus/T); vals=np.array([F(m_) for m_ in mus])
    for d in range(1,dmax):
        M=np.hstack([np.vander(taus,d+1,increasing=True), -vals[:,None]*np.vander(taus,d+1,increasing=True)])
        _,sv,Vh=np.linalg.svd(M)
        if sv[-1]/sv[0]<1e-9:
            v=Vh[-1].conj(); p=v[:d+1]; q=v[d+1:]
            z=np.roots(p[::-1]); pl=np.roots(q[::-1])
            # mu-values: tau = e^{i pi mu/T} -> mu = T*angle/pi (mod 2T); report |tau| too
            conv=lambda r: sorted([(round(T*np.angle(x)/np.pi % (2*T),5), round(abs(x),4)) for x in r if 1e-6<abs(x)<1e6])
            return d, conv(z), conv(pl), sv[-1]/sv[0]
    return None
if __name__=='__main__':
 for om in (0.37,):
  for n in (1,2,3):
    T=(n+1)*om/2
    print('n=',n,'om=',om,'T=',T, '2T=',2*T)
    print('  R     :',fitfun(lambda m: Rplus(m,n,om),T))
    print('  target:',fitfun(lambda m: target_plus(m,n,om),T))
