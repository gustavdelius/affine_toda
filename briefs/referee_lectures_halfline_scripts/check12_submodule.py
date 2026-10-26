import numpy as np
from math import pi,sqrt
lam=sqrt(2);eta=2.1;v=.43;q=-np.exp(-1j*pi*lam)
e1=-np.cos(eta-1j*v)/np.sin(pi*lam);e0=-np.cos(eta+1j*v)/np.sin(pi*lam)
def matrices(y):return (np.array([[e1*q,y],[1/y,e1/q]]),np.array([[e0/q,1/y],[y,e0*q]]))
y=1j*np.exp(-1j*eta);vec=np.array([1,np.exp(-v)])
betas=[-np.cos(pi*(lam+1)-eta+1j*v)/np.sin(pi*lam),-np.cos(pi*(lam+1)-eta-1j*v)/np.sin(pi*lam)]
for i,(A,b) in enumerate(zip(matrices(y),betas)):print('12 eigenvalue residual node',1-i,np.linalg.norm(A@vec-b*vec))
for y2 in [-np.exp(2j*eta),-np.exp(-2j*eta),-np.exp(2*v),-np.exp(-2*v)]:
    y=np.sqrt(y2+0j);A,D=matrices(y);vals,vs=np.linalg.eig(A)
    residual=min(np.linalg.norm(D@vs[:,j]-np.vdot(vs[:,j],D@vs[:,j])*vs[:,j]) for j in range(2))
    print('12 branch y2',y2,'common eigenvector residual',residual);assert residual<1e-12
print('PASS 12: physical-pole eigenvalues correct; four generic reducibility branches')
