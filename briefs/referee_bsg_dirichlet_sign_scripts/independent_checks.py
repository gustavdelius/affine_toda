"""Independent checks from SG PDE and the displayed Gamma product; no book imports."""
import numpy as np
from scipy.integrate import quad
from scipy.special import digamma, loggamma
from itertools import product
from math import pi

def tau_derivatives(x,t,theta,a,charge):
    g=np.cosh(theta);h=np.sinh(theta);v=np.tanh(theta);w=np.tanh(theta/2)**2;r=np.tan(a/4)
    D=np.sqrt((w+r*r)/(v*v*w*(1+r*r*w)))
    z=np.array([charge*D*np.exp(g*x-h*t),charge*D*np.exp(g*x+h*t),r/w*np.exp(x)])
    gx=np.array([g,g,1.]);gt=np.array([-h,h,0.])
    A={(0,1):v*v,(0,2):w,(1,2):w}
    tau=tx=tt=txx=ttt=0j
    for bits in product((0,1),repeat=3):
        term=(1j)**sum(bits)
        for j in range(3):
            if bits[j]:term*=z[j]
        for (j,k),val in A.items():
            if bits[j] and bits[k]:term*=val
        kx=np.dot(bits,gx);kt=np.dot(bits,gt)
        tau+=term;tx+=kx*term;tt+=kt*term;txx+=kx*kx*term;ttt+=kt*kt*term
    return tau,tx,tt,txx,ttt

def ux(t,th,a,s):
    z,zx,*_=tau_derivatives(0,t,th,a,s)
    return 4*(zx/z).imag

def phase_derivative(eta,theta,lam,N):
    # d_eta log sigma, directly differentiating the Gamma product in ch.30.
    l=np.arange(1,N+1,dtype=float);t=-1j*lam*theta/pi;x=eta/pi
    def logslope(t):
        return (digamma(.5+x+(2*l-1)*lam-t)-digamma(.5-x+(2*l-1)*lam-t)
               +digamma(.5-x+(2*l-2)*lam-t)-digamma(.5+x+2*l*lam-t)).sum()/pi
    return (-np.tan(eta)+np.tan(eta-1j*lam*theta)+logslope(t)-logslope(-t)).imag

if __name__=='__main__':
    print('Exact tau: PDE and boundary residuals')
    errp=errb=0
    for theta in [.4,.8,1.3]:
      for a in [.6,1.4,-.9]:
       for s in [1,-1]:
        for x,t in [(-2.,-1.),(-.3,.7),(-3.,2.)]:
         z,zx,zt,zxx,ztt=tau_derivatives(x,t,theta,a,s)
         residual=4*((ztt/z-(zt/z)**2)-(zxx/z-(zx/z)**2)).imag+np.sin(4*np.angle(z))
         errp=max(errp,abs(residual))
        for t in [-8.,-.2,4.]:
         z,*_=tau_derivatives(0,t,theta,a,s)
         errb=max(errb,abs(np.exp(1j*(4*np.angle(z)-a))-1))
    print(f'PASS PDE max residual {errp:.3e}; boundary max residual {errb:.3e}')
    for theta in [.4,.8,1.3]:
      for a in [.6,1.4,-.9]:
        ground=2*np.sin(a/2);T=32/np.sinh(theta)
        ints=[quad(lambda t:ux(t,theta,a,s)-ground,-T,T,epsabs=1e-10)[0] for s in [1,-1]]
        print(f'theta={theta:.1f} velocity={np.tanh(theta):.6f} a={a:+.1f} I+={ints[0]:.10f} I-={ints[1]:.10f} diff-8pi={ints[0]-ints[1]-8*pi:.2e}')
        if a>0:
         target=-sum(ints)/(8*pi)
         vals=[phase_derivative((lam+1)*a/2,theta,lam,4000) for lam in [30+np.sqrt(2),100+np.sqrt(2),300+np.sqrt(2)]]
         print(f'  sigma derivative target={target:.8f} lambda=31,101,301: {vals}')
    theta=.8;a=1.4;lam=300+np.sqrt(2)
    print('Gamma product cutoff:',[(N,phase_derivative((lam+1)*a/2,theta,lam,N)) for N in [1000,4000,16000]])
