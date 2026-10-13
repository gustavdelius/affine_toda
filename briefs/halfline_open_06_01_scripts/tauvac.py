"""Static all-species Hirota solutions of real-coupling a_n^(1) (m=beta=1) and the reduced (Bogomolny) energy
   F_P(c) = W + B_P,  W = -2 sum_j tau_j'(0)/tau_j(0)   (energy of the static BPS solution, Bowcock-Perkins (3.1)).
"""
import numpy as np, itertools
from halfline_min import roots, basis

class Tau:
    def __init__(self,n):
        self.n=n; self.h=h=n+1
        self.ma=np.array([2*np.sin(np.pi*a/h) for a in range(1,h)])
        a=np.arange(1,h)
        self.Aab=np.array([[ (np.sin(np.pi*(p-q)/(2*h))**2/np.sin(np.pi*(p+q)/(2*h))**2) for q in a] for p in a])
        self.mus=np.array(list(itertools.product([0,1],repeat=n)))             # 2^n x n
        # interaction log
        self.lnA=np.zeros(len(self.mus))
        for k,mu in enumerate(self.mus):
            s=0.0
            for p in range(n):
                for q in range(p+1,n):
                    if mu[p] and mu[q]:
                        s+=np.log(self.Aab[p,q]) if self.Aab[p,q]>0 else -np.inf
            self.lnA[k]=s
        self.AA=np.exp(self.lnA)
        self.al=roots(n)
    def coeffs(self,cs):
        """cs: complex array of length n, cs[a-1]=c_a (c_{h-a}=conj c_a required for reality)."""
        h=self.h; jj=np.arange(h)
        w=np.exp(2j*np.pi/h)
        f=np.array([cs[a]*w**(jj*(a+1)) for a in range(self.n)])              # n x h
        G=np.where(self.mus[None,:,:]==1, f.T[:,None,:], 1.0)                  # h x K x n
        return G.prod(axis=2)*self.AA[None,:]
    def mass_exp(self):
        return self.mus@self.ma                                              # exponent rate per term
    def tau(self,cs,x):
        """returns tau (h x len(x)), dtau, ddtau (real parts; imag part should vanish)"""
        T=self.coeffs(cs); r=self.mass_exp()
        E=np.exp(np.outer(r,x))                                                # terms x len(x)
        t=T@E; d=(T*r[None,:])@E; dd=(T*(r**2)[None,:])@E
        return t,d,dd
def conj_cs(n,ds,chis,mid=1.0):
    """ds, chis for a=1..floor(n/2) ; mid = real coefficient for a=h/2 (h even)."""
    h=n+1; cs=np.zeros(n,complex)
    for a in range(1,(h+1)//2):
        if 2*a==h: continue
        cs[a-1]=ds[a-1]*np.exp(1j*chis[a-1]); cs[h-a-1]=np.conj(cs[a-1])
    if h%2==0: cs[h//2-1]=mid
    return cs
def nparams(n):
    h=n+1; return 2*((h-1)//2)+(1 if h%2==0 else 0)
def cs_from_p(n,p):
    """p = [ln d_1, chi_1, ln d_2, chi_2, ..., (mid_real)] """
    h=n+1; k=(h-1)//2; cs=np.zeros(n,complex)
    for a in range(1,k+1):
        d=np.exp(p[2*a-2]); cs[a-1]=d*np.exp(1j*p[2*a-1]); cs[h-a-1]=np.conj(cs[a-1])
    if h%2==0: cs[h//2-1]=p[-1]
    return cs

def field_quantities(T,cs,x):
    t,d,dd=T.tau(cs,x); return t.real,d.real,dd.real,np.abs(t.imag).max()

def F_of(T,cs,C,x):
    """returns (F, W, psi0, feasible, tau0)  using grid x (increasing, last = 0)"""
    t,d,dd=T.tau(cs,x)
    t=t.real; d=d.real
    feas=(t>0).all()
    t0=t[:,-1]; d0=d[:,-1]
    h=T.h
    W=-2*np.sum(d0/t0)
    eh=np.sqrt(np.abs(np.roll(t0,1)*np.roll(t0,-1)))/np.abs(t0)
    B=2*np.sum(np.array(C)*(eh-1))
    psi0=2*np.log(np.maximum(eh,1e-300))
    return W+B,W,psi0,feas,t0
