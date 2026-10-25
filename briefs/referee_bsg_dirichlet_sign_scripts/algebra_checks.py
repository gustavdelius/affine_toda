import numpy as np
from scipy.special import loggamma
from scipy.integrate import quad
from pathlib import Path
import re,subprocess
pi=np.pi

def sigma(x,u,L,N=20000):
 l=np.arange(1,N+1);t=L*u/pi;z=x/pi
 def f(t):
  return (loggamma((.5+z+(2*l-1)*L-t).astype(complex))+loggamma((.5-z+(2*l-1)*L-t).astype(complex))-loggamma((.5-z+(2*l-2)*L-t).astype(complex))-loggamma((.5+z+2*l*L-t).astype(complex))).sum()
 return np.cos(x)/np.cos(x+L*u)*np.exp(f(t)-f(-t))
L=2+np.sqrt(2);eta=2.2;nu=(eta-pi/2)/L
for ep in [1e-3,1e-4,1e-5]:
 u=nu+ep;z=sigma(eta,u,L)
 print('pole',ep,'P+',abs(np.cos(L*u+eta)*z),'P-',abs(np.cos(L*u-eta)*z))
print('sigma even ratio',sigma(-eta,-.8j,L)/sigma(eta,-.8j,L))
for a in [.6,1.4,pi]:
 E=[]
 for excited in [False,True]:
  r=np.tan(a/4) if not excited else 1/np.tan(a/4)
  def energy(x):
   u=4*np.arctan(r*np.exp(x));ux=4*r*np.exp(x)/(1+r*r*np.exp(2*x))
   return .5*ux*ux+1-np.cos(u)
  E.append(quad(energy,-50,0,epsabs=1e-11)[0])
 print('static',a,E,'gap',E[1]-E[0],'predicted',8*np.cos(a/2))
for L in [np.sqrt(2),np.e,2+np.sqrt(2)]:
 f=lambda z:-np.cos((L+1)*np.arccos(z+0j))/np.sin(pi*L)
 h=1e-6
 print('node Taylor',L,abs(f(0)-1/(2*np.cos(pi*L/2))),abs((f(h)-f(-h))/(2*h)+(L+1)/(2*np.sin(pi*L/2))))
for td,a in [(.6,.8),(-.35,2.3)]:
 C0=np.cosh(td+1j*a/2);C1=np.cosh(td-1j*a/2)
 print('fusion residual',abs(np.cos(-a/2-1j*td)-C1),abs(np.cos(-a/2+1j*td)-C0))
B=.73;b=np.sqrt(8*pi*B/(2-B));phi=.37;L=-2/B
eta=-(L+1)*1j*b*phi/2
print('E_D continuation residual',abs(B*eta/pi-1j*(1-B/2)*b*phi/pi))
# Minimal two-component double-row example: sign flip conjugates K by C.
C=np.array([[0,1],[1,0]])
def K(e,v,th,L):
 u=-1j*th
 return np.array([[np.cos(L*u)*np.cos(e)*np.cosh(v)-np.sin(L*u)*np.sin(e)*np.sinh(v),-np.sin(L*u)*np.cos(L*u)],[-np.sin(L*u)*np.cos(L*u),np.cos(L*u)*np.cos(e)*np.cosh(v)+np.sin(L*u)*np.sin(e)*np.sinh(v)]])
A=K(.9,.7,.5,np.sqrt(2));Am=K(-.9,.7,.5,np.sqrt(2));N=K(pi*(np.sqrt(2)+1)/2,0,.5,np.sqrt(2))
print('C conjugation / Neumann commutation',np.max(abs(Am-C@A@C)),np.max(abs(N-C@N@C)))
print('Neumann double-row trace,det sign invariance',abs(np.trace(A@N)-np.trace(Am@N)),abs(np.linalg.det(A@N)-np.linalg.det(Am@N)))
# A strip zero-mode computation fixes the sign of its winding label.
beta=1.3;dphi=.4;R=np.sqrt(4*pi)/beta
for Q in [-1,0,1]:
 direct=(dphi+2*pi*Q/beta)**2/(2*pi)
 nlie=.5*(dphi/np.sqrt(pi)+Q*R)**2
 print('strip winding',Q,'direct',direct,'NLIE formula',nlie)
print('NLIE H+>0,H->0: m=-S; one hole S=1 has m=-1 (not +1)')
