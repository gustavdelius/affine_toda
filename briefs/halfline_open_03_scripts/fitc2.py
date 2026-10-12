import numpy as np
from proj import cproj
def fit(n,s,om,dmax=8):
    T=(n+1)*om/2
    rng=np.random.default_rng(5)
    m=4*dmax+10
    ths=rng.normal(size=m)*0.5/T+1j*rng.uniform(-1,1,size=m)*np.pi/(2*T)
    ts=np.exp(2*T*ths); cs=np.array([cproj(n,s,om,t)[0] for t in ths])
    for d in range(1,dmax):
        M=np.hstack([np.vander(ts,d+1,increasing=True), -cs[:,None]*np.vander(ts,d+1,increasing=True)])
        _,sv,Vh=np.linalg.svd(M)
        if sv[-1]/sv[0]<1e-11:
            v=Vh[-1].conj(); p=v[:d+1]; qq=v[d+1:]
            return d, np.roots(p[::-1]), np.roots(qq[::-1]), p[-1]/qq[-1]
om=0.37
for n in (2,3,4,5):
  T=(n+1)*om/2
  for s in (1,-1):
    d,z,p,lead=fit(n,s,om)
    # t = e^{2T th}; t = e^{2 pi i x} <-> 2T th = 2 pi i x ; report log(t)/(i pi) = 2x 
    f=lambda r: sorted([ (round(abs(x),6), round((np.angle(x)/np.pi)%2,6)) for x in r])
    print(n,s,'deg',d,'lead',np.round(lead,6),' zeros',f(z),' poles',f(p))
    print('    om, T, 1/2-T, T-om ...', om%2, T%2, (1-T)%2, (T-om)%2, (2*T)%2,(-T+om)%2)
