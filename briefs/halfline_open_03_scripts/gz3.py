import numpy as np
from gz import *
from rhonp import target_plus, Rplus
lam=0.37; n=1; T=lam; eta=np.pi*(lam+1)/2
def X(mu):
    th=1j*np.pi*mu/T; u=-1j*th
    return A_gz(u,lam,eta)*np.cos(lam*u)/(A1(mu,n,lam)*rho_plus(mu,n,lam))
mu=0.13+0.21j
print('X unitarity', X(mu)*X(-mu), ' crossing', X(T/2-mu)/X(T/2+mu))
print('period 2T', X(mu+2*T)/X(mu), ' period 1', X(mu+1)/X(mu), 'period 2',X(mu+2)/X(mu))
# Ghoshal at Neumann vs DG311*target
def gh(mu,eta,vth=0):
    u=np.pi*mu/T   # u=-i theta, mu=-iT theta/pi -> u = pi mu/T
    l=lam
    pre=(np.cos(u/2+np.pi/(4*l))*np.cos(u/2-np.pi/4-np.pi/(4*l))*np.sin(u/2+np.pi/4))/(np.cos(u/2-np.pi/(4*l))*np.cos(u/2+np.pi/4+np.pi/(4*l))*np.sin(u/2-np.pi/4))
    return pre*(np.cos(eta/l)-np.sin(u))/(np.cos(eta/l)+np.sin(u))*(np.cosh(vth/l)-np.sin(u))/(np.cosh(vth/l)+np.sin(u))
for mu in (0.13+0.21j,-0.3+0.1j):
    print('Ghoshal_N/Ghoshal_diag =',gh(mu,eta)/gh(mu,np.pi/2),' target =',target_plus(mu,n,lam), ' R=',Rplus(mu,n,lam))
