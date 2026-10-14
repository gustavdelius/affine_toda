import numpy as np
from b5_boundary_krein import transfer, sector_eigs, N3
Z=np.zeros((3,3))
for i in range(3): Z[(i+1)%3,i]=1
Zi=Z.T
om=0.6; ths=[0.83,0.31]
for lw,rw in (('U','N'),('N','N')):
    Ts=[transfer(j,('V','W'),ths,om,lw,rw) for j in range(2)]
    for name,ZZ in (('Z x Z',np.kron(Z,Z)),('Z x Z^-1',np.kron(Z,Zi))):
        print(lw+rw,name,'commutator',max(np.abs(T@ZZ-ZZ@T).max() for T in Ts))
# eigenvalue moduli in each sector definition, near the book's threshold
from b5_boundary_krein import Zop
def maxbreak_with(ZZ, om, lw, rw):
    m=0
    for t1 in (10.0,12.0):
        for t2 in (1e-4,1e-3,1e-2):
            Ts=[transfer(j,('V','W'),[t1,t2],om,lw,rw) for j in range(2)]
            lam=sector_eigs(Ts,ZZ)
            m=max(m,np.max(np.abs(np.abs(lam)-1)))
    return m
for lw,rw,oms in (('U','N',(0.6,0.66,0.67,0.69,0.70)),('N','N',(0.45,0.50,0.51,0.53,0.54))):
    for om in oms:
        print(lw+rw, om, 'ZxZ sector:', f"{maxbreak_with(np.kron(Z,Z),om,lw,rw):.1e}", ' ZxZ^-1 sector:', f"{maxbreak_with(np.kron(Z,Zi),om,lw,rw):.1e}")
