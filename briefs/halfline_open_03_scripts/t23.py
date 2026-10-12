import numpy as np
from rhonp import rho_plus
from gz import A1
from rhonp import blk
# Neumann: Delta t = ln(1-s^2)/(m sinh th), s = sin(pi/h)/cosh th; C=+1: 2 ln(1+s)/(m sinh th) <-> (beta^2/4h) d(-i ln A_d)/dth -> ln(1+s)
for n in (2,3):
    h=n+1; th=0.6
    s=np.sin(np.pi/h)/np.cosh(th)
    tgt=0.5*np.log(1-s**2)
    for L in (40.,80.,160.):
        T=h*L/2
        def lnA(t):
            mu=-1j*T*t/np.pi+0j
            sig=blk(T/2,mu,T)**2*blk(T/2-L/2,mu,T)*blk(T/2+L/2,mu,T)
            return np.log(A1(mu,n,L,K=3000)*rho_plus(mu,n,L,K=3000)*sig)
        e=1e-5; d=(lnA(th+e)-lnA(th-e))/(2*e)
        beta2=4*np.pi/(L+1)
        print(n,L,'(b^2/4h) d(-i lnA+)/dth =',np.real(-1j*d)*beta2/(4*h),' target (1/2)ln(1-s^2) =',tgt)
