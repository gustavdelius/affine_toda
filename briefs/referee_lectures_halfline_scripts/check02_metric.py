import numpy as np
from scipy.integrate import quad
R=np.array([[-1.,0,1],[1,-1,0],[0,1,-1]])
s=np.sqrt(3)/2; ma=np.sqrt(3)
def norm(z,chi):
    v=chi+2*np.pi*np.arange(3)/3
    tau=((1-z)**2+s*(1-z*z)+4*z*np.cos(v/2)**2)/(1+s)
    dt=-2*np.sin(v)*z/(1+s)
    vec=-(dt/tau)@R
    return vec@vec/(ma*z) if z else 0.
for chi in [0.,.25,.6,1.,np.pi/3-.001]:
    I,err=quad(norm,0,1,args=(chi,),epsabs=1e-10,epsrel=1e-10,points=[.9,.99,.9999,.999999])
    print(f'chi={chi:.12g} I={I:.12g} I*sqrt(3)={I*ma:.12g} quad_error={err:.2g}')
    assert abs(I-6/ma)<1e-7
print('PASS: I=6/(beta^2*m1), for beta=m=1, across each sampled interval')
