import numpy as np
from gz import logR0, logsigma, A1
from rhom import rho_minus
lam=0.37; n=1; T=lam
for th in (0.3+0.2j,-0.5+0.4j,0.8-0.1j):
    u=-1j*th; mu=-1j*T*th/np.pi
    eta=np.pi/2-1e-7; vth=1j*np.pi*lam/2
    AG=np.exp(logR0(u,lam)+logsigma(eta,u,lam)+logsigma(1j*vth,u,lam))/(np.cos(eta)*np.cosh(vth))*1j*np.sin(lam*u)*np.sin(np.pi*lam/2)
    print(th,'A_GZ(s=-1)/(A_d rho_-) =',AG/(A1(mu,n,lam)*rho_minus(mu,n,lam)))
