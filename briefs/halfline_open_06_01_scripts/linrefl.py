"""Linearized reflection matrix about a static real-coupling vacuum (tau-form), a_n^(1), m=beta=1.
   y'' = (M^2(x) - Om^2) y,  y'(0) = -H y(0),  H = (1/2) sum_j C_j e^{psi_j(0)/2} A_j A_j^T,  M^2(x) = sum_j e^{psi_j(x)} A_j A_j^T."""
import numpy as np
from scipy.integrate import solve_ivp
from halfline_min import roots, basis
from tauvac import Tau, cs_from_p

def refl_matrix(n,taufun,C,Om,X=22.0,rtol=1e-11,atol=1e-13,x0=0.0):
    """taufun(x)->tau array (h,) at point x (real positive).  returns (K, labels) K in Fourier basis of open channels."""
    h=n+1; al=roots(n); Q=basis(n); A=al@Q
    ma=np.array([2*np.sin(np.pi*a/h) for a in range(1,h)])
    def psi(x):
        t=taufun(x); return np.log(np.roll(t,1)*np.roll(t,-1))-2*np.log(t)
    def M2(x):
        e=np.exp(psi(x)); return (A.T*e)@A
    ps0=psi(x0); H=0.5*(A.T*(np.asarray(C)*np.exp(ps0/2)))@A
    def rhs(x,Y):
        y=Y[:n*n].reshape(n,n); yp=Y[n*n:].reshape(n,n)
        return np.concatenate([yp.ravel(),((M2(x)-Om**2*np.eye(n))@y).ravel()])
    Y0=np.concatenate([np.eye(n).ravel(),(-H).ravel()])
    sol=solve_ivp(rhs,[x0,-X],Y0,method='DOP853',rtol=rtol,atol=atol)
    Y=sol.y[:,-1]; y=Y[:n*n].reshape(n,n); yp=Y[n*n:].reshape(n,n)
    w=np.exp(2j*np.pi/h); F=np.array([[w**(k*a) for k in range(h)] for a in range(1,h)])/np.sqrt(h)   # (n x h) f_a
    P=np.conj(F)@Q                                              # c_a = f_a^dagger Q y
    c=P@y; cp=P@yp                                              # n_ch x n (columns = basis solutions)
    xs=-X
    open_=[a for a in range(n) if Om>ma[a]]
    closed=[a for a in range(n) if Om<=ma[a]]
    pa={a:np.sqrt(Om**2-ma[a]**2) for a in open_}; ka={a:np.sqrt(ma[a]**2-Om**2) for a in closed}
    # constraints: G_a = (c - c'/kappa)/2 = 0
    if closed:
        G=np.array([(c[a]-cp[a]/ka[a])/2 for a in closed])
        u,s,vh=np.linalg.svd(G); nn=vh.conj().T[:,len(closed):]   # null space (n x #open)
    else: nn=np.eye(n)
    # amplitudes: c = A e^{ipx}+B e^{-ipx}; c' = ip(A e^{ipx}-B e^{-ipx})
    Aa=[];Bb=[]
    for a in open_:
        p=pa[a]; ca=c[a]@nn; cpa=cp[a]@nn
        Aa.append((ca+cpa/(1j*p))/2*np.exp(-1j*p*xs)); Bb.append((ca-cpa/(1j*p))/2*np.exp(1j*p*xs))
    Aa=np.array(Aa); Bb=np.array(Bb)
    K=Bb@np.linalg.inv(Aa)
    return K,[a+1 for a in open_]

if __name__=='__main__':
    from halfline_min import roots
    def blk(x,th,h): return np.sinh(th/2+1j*np.pi*x/(2*h))/np.sinh(th/2-1j*np.pi*x/(2*h))
    for n in (2,3,4):
        h=n+1
        for Cc in (1,-1):
            Om=2*np.sin(np.pi/h)*np.cosh(0.6); th=0.6
            K,lab=refl_matrix(n,lambda x:np.ones(h),[Cc]*h,Om)
            print(n,Cc,lab,np.round(K,6).tolist(), 'formula',-(blk(1,th,h)*blk(h-1,th,h))**(-Cc))
