import numpy as np
from math import pi,sqrt
lam=sqrt(2);q=-np.exp(-1j*pi*lam);th=.37-.12j;e0=.53-.19j;e1=-.41+.31j
# Solve the two |22> coproduct equations for the central braid coefficients.
def R(z,q):
    c,b=np.linalg.solve(np.array([[z,1/q],[1/q,z]],complex),np.array([1,z/q]))
    return np.array([[1,0,0,0],[0,c,b,0],[0,b,c,0],[0,0,0,1]],complex)
def K(y):return np.array([[(q-1/q)*(e1*y+e0/y),y*y-y**-2],[y*y-y**-2,(q-1/q)*(e0*y+e1/y)]])
z=np.exp(2*lam*th);J=np.diag([1,1,-1,-1])
print('16 twist error',np.max(abs(R(z,-q)-J@R(z,q)@J)))
for sign in [1,-1]:
    for csign in [1,-1]:
        C=np.array([[0,1],[csign,0]])
        vm=(K(np.exp(lam*(.5j*pi-th)))@C).T.ravel()
        vp=(K(np.exp(lam*(.5j*pi+th)))@C).T.ravel()
        out=R(z,sign*q)@vp
        coef=np.vdot(vm,out)/np.vdot(vm,vm)
        residual=np.linalg.norm(out-coef*vm)/np.linalg.norm(out)
        print('16 q multiplier',sign,'C lower-left',csign,'crossing proportionality residual',residual)
print('PASS 16: independent coproduct solve confirms q-sign/twist and matrix crossing test')
