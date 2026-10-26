import numpy as np
from itertools import product
R=np.array([[-1.,0,1],[1,-1,0],[0,1,-1]])
s=np.sqrt(3)/2; ma=np.sqrt(3)
def data(x,chi):
    z=np.exp(ma*x); v=chi+2*np.pi*np.arange(3)/3
    t=1+2*np.cos(v)*z/(1+s)+(1-s)*z*z/(1+s)
    tx=ma*(2*np.cos(v)*z/(1+s)+2*(1-s)*z*z/(1+s))
    phi=-np.log(t)@R; phix=-(tx/t)@R
    return phi,phix
for chi in [.25,2.,-2.]:
    phi,phix=data(0,chi); e=np.exp(R@phi/2)
    vals=sorted((np.linalg.norm(phix+(np.array(c)*e)@R),c) for c in product([-1,1],repeat=3))
    printed={.25:(-1,1,-1),2.:(-1,-1,1),-2.:(1,-1,-1)}[chi]
    print('chi',chi,'actual',vals[0],'printed residual',np.linalg.norm(phix+(np.array(printed)*e)@R))
    assert vals[0][0]<1e-12

print('PASS: direct boundary equation identifies the sign patterns independently')
