import numpy as np
from gz import *
lam=0.37; n=1; T=lam
for th in (0.3+0.2j,-0.5+0.4j):
    u=-1j*th; mu=-1j*T*th/np.pi
    eta=np.pi/2-1e-7
    off=np.exp(logR0(u,lam)+logsigma(eta,u,lam)+logsigma(0,u,lam))*(-np.sin(lam*u)*np.cos(lam*u))/np.cos(eta)
    print(th,'GZ diag-point offdiag / A_d =',off/A1(mu,n,lam))
