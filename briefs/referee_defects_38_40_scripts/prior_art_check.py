"""Closed-path formula derived independently, compared with Rossi--Weston (5.5).
No book imports. quantum_checks is this referee's independent implementation.
"""
import numpy as np
from quantum_checks import d_operator,relative

def closed(q,ys,yd,offset):
    N=len(ys);out=np.zeros((2**N,2**N),complex)
    spins=[np.array([1-2*((c>>(N-1-j))&1) for j in range(N)]) for c in range(2**N)]
    for col,a in enumerate(spins):
        m=sum(a)
        for row,b in enumerate(spins):
            if sum(b)!=m:continue
            omega=sum(a[i]*b[j]-b[i]*a[j] for i in range(N) for j in range(i+1,N))/4
            value=q**(offset*m/2+omega)
            for j,z in enumerate(np.asarray(ys)/yd):
                flip=(1-a[j]*b[j])//2
                value*=(1j*z*q**(-.5))**flip/(1+z*z)
            out[row,col]=value
    return out

if __name__=='__main__':
    q=np.exp(.71j);yd=.79;s=.3
    for N in range(2,7):
        for homogeneous in (True,False):
            ys=np.full(N,1.13) if homogeneous else np.exp(np.linspace(-.43,.51,N))
            err=relative(closed(q,ys,yd,s),d_operator(q,ys,yd,s))
            print('N=',N,'homogeneous=',homogeneous,'closed formula residual=',err)
            assert err<1e-12
    print('Rossi--Weston (5.5): q_RW=q^-1, s0=q^-s, z_RW/w_RW=-(y/y_d)^2,')
    print('rho_tilde=(1+(y/y_d)^2)^-1 and sqrt(z_RW*q_RW/w_RW)=i*(y/y_d)*q^-1/2.')
    print('With delta(q_RW^m) removed, (5.5) is the homogeneous closed formula above.')
