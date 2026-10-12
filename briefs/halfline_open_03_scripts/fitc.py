import numpy as np
from proj import cproj
om=0.37
def fit(n,s,dmax=14):
    rng=np.random.default_rng(3)
    m=4*dmax+10
    ths=rng.normal(size=m)*0.6+1j*rng.normal(size=m)*0.5
    us=np.exp(om*ths); cs=np.array([cproj(n,s,om,t)[0] for t in ths])
    for d in range(1,dmax):
        for e in range(0,d+1):   # allow overall power u^e shift: P(u)/Q(u)
            pass
        M=np.hstack([np.vander(us,d+1,increasing=True), -cs[:,None]*np.vander(us,d+1,increasing=True)])
        _,sv,Vh=np.linalg.svd(M)
        if sv[-1]/sv[0]<1e-11:
            v=Vh[-1].conj(); p=v[:d+1]; qq=v[d+1:]
            return d, np.roots(p[::-1]), np.roots(qq[::-1]), sv[-2:]/sv[0]
for n in (2,3):
  for s in (1,-1):
    d,z,p,sv=fit(n,s)
    print(n,s,'deg',d,sv)
    f=lambda r: [ (round(abs(x),5), round(np.angle(x)/np.pi,5)) for x in sorted(r,key=lambda x:(abs(x),np.angle(x)))]
    print(' zeros (|u|, arg/pi):',f(z)); print(' poles:',f(p))
