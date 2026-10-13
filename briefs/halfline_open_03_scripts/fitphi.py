import numpy as np
from phi import Phi
om=0.37
from boundary_common import q_phys
q=q_phys(om); w=np.sqrt(-q+0j)
def U(n,s,th):
    T=(n+1)*om/2
    return (s*np.cosh(2*T*th)-np.cos(np.pi*T))/(2*np.sin(np.pi*om/2)**2)
from phi import U as Unum
for n in (1,2,3): 
  for s in (1,-1):
    th=0.31+0.2j; print('U check',n,s,abs(U(n,s,th)-Unum(n,s,om,th)[0]))
# rational fit of Phi in u = e^{om th}, as P(u)/Q(u), deg d
def fit(n,s,d):
    rng=np.random.default_rng(1)
    ths=rng.normal(size=4*d+10)*0.8+1j*rng.normal(size=4*d+10)*1.0
    us=np.exp(om*ths); ph=np.array([Phi(n,s,om,t)[0] for t in ths])
    # P(u) - ph Q(u) = 0, unknowns p0..pd, q0..qd
    M=np.hstack([np.vander(us,d+1,increasing=True), -ph[:,None]*np.vander(us,d+1,increasing=True)])
    _,sv,Vh=np.linalg.svd(M); v=Vh[-1].conj()
    return sv[-3:]/sv[0], v[:d+1], v[d+1:]
for n in (1,2,3):
  for s in (1,-1):
    for d in range(1,2*(n+1)+3):
        sv,p,qq=fit(n,s,d)
        if sv[-1]<1e-10:
            # strip common factors: roots
            rp=np.roots(p[::-1]); rq=np.roots(qq[::-1])
            print(n,s,'deg',d,'sv',sv)
            print('  zeros', np.round(np.sort_complex(rp),6))
            print('  poles', np.round(np.sort_complex(rq),6))
            break
