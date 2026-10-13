"""Stable linearized reflection matrix: integrate from x=-X (asymptotic modes) to x=0, impose y'(0) = -H y(0)."""
import numpy as np
from scipy.integrate import solve_ivp
from halfline_min import roots, basis
def refl_matrix2(n,taufun,C,Om,X=16.0,rtol=1e-12,atol=1e-14,x0=0.0,Xcheck=None):
    h=n+1; al=roots(n); Q=basis(n); A=al@Q
    ma=np.array([2*np.sin(np.pi*a/h) for a in range(1,h)])
    def psi(x):
        t=taufun(x); return np.log(np.roll(t,1)*np.roll(t,-1))-2*np.log(t)
    def M2(x):
        e=np.exp(psi(x)); return (A.T*e)@A
    H=0.5*(A.T*(np.asarray(C)*np.exp(psi(x0)/2)))@A
    w=np.exp(2j*np.pi/h); Fm=np.array([[w**(k*a) for k in range(h)] for a in range(1,h)])/np.sqrt(h)
    U=(Q.T@Fm.T)                       # n x n_ch : y-vector of channel a
    open_=[a for a in range(n) if Om>ma[a]+1e-12]; closed=[a for a in range(n) if a not in open_]
    modes=[]  # (kind, a)
    cols=[]
    xs=-X
    def init(a,kind):
        if kind=='in': p=np.sqrt(Om**2-ma[a]**2); f=np.exp(1j*p*xs); fp=1j*p*f
        elif kind=='out': p=np.sqrt(Om**2-ma[a]**2); f=np.exp(-1j*p*xs); fp=-1j*p*f
        else: k=np.sqrt(ma[a]**2-Om**2); f=np.exp(k*xs); fp=k*f
        return U[:,a]*f, U[:,a]*fp
    sols=[]
    for kind,lst in (('in',open_),('out',open_),('cl',closed)):
        for a in lst:
            y0,yp0=init(a,kind)
            def rhs(x,Y):
                y=Y[:n]; yp=Y[n:]
                return np.concatenate([yp,(M2(x)-Om**2*np.eye(n))@y])
            s=solve_ivp(rhs,[xs,x0],np.concatenate([y0,yp0]).astype(complex),method='DOP853',rtol=rtol,atol=atol)
            Y=s.y[:,-1]; sols.append((kind,a,Y[n:]+H@Y[:n]))
    nin=len(open_)
    Rin=np.array([r for k,a,r in sols if k=='in']).T            # n x nin
    Rrest=np.array([r for k,a,r in sols if k!='in']).T          # n x n
    sol=np.linalg.solve(Rrest,-Rin)                              # n x nin ; first nin rows = B
    K=sol[:nin,:]
    return K,[a+1 for a in open_]
