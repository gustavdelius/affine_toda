"""Shared static energy minimization for real-coupling a_n^(1) Toda on x<=0 (book conventions, m=beta=1).
 E[phi] = int_{-L}^0 (phi'^2/2 + sum_j (e^{alpha_j.phi}-1)) dx + 2 sum_j C_j (e^{alpha_j.phi(0)/2}-1),  phi(-L)=0.
 alpha_j = e_j - e_{j+1 mod h}, j=0..n.  Energies in units m/beta^2.  Book -> Bowcock-Perkins: A_j=-C_j, E_book = -2(n+1) C_BP - 2 sum C_j.
"""
import numpy as np, itertools
from scipy.optimize import minimize

def roots(n):
    h=n+1; e=np.eye(h); return np.array([e[j]-e[(j+1)%h] for j in range(h)])

def basis(n):
    h=n+1
    # orthonormal basis of sum-zero subspace
    M=np.eye(h)-np.ones((h,h))/h
    u,s,vt=np.linalg.svd(M); return u[:,:n]          # h x n

def grid(L=25.0, n_wall=60, dx0=1e-3, dxmax=0.15, q=1.12):
    """nonuniform grid on [-L,0]: geometric refinement towards 0, then growth up to dxmax."""
    pts=[0.0]; dx=dx0
    while pts[-1]>-L:
        pts.append(pts[-1]-dx); dx=min(dx*q,dxmax)
    x=np.array(pts[::-1]); x[0]=-L if x[0]<-L else x[0]
    return x

class Prob:
    def __init__(self,n,C,L=25.0,**gk):
        self.n=n; self.h=n+1; self.C=np.array(C,float)
        self.al=roots(n); self.Q=basis(n); self.A=self.al@self.Q      # (h x n): alpha_j.phi = A y
        self.x=grid(L,**gk); self.hh=np.diff(self.x); self.N=len(self.x)
    def unpack(self,v): return v.reshape(self.N-1,self.n)            # y_1..y_N ; y_0 = 0
    def E_grad(self,v):
        Y=np.vstack([np.zeros((1,self.n)),self.unpack(v)])
        d=np.diff(Y,axis=0)/self.hh[:,None]
        psi=Y@self.A.T; ep=np.exp(np.clip(psi,-200,200)); V=(ep-1).sum(1)
        w=np.zeros(self.N); w[:-1]+=self.hh/2; w[1:]+=self.hh/2
        E=0.5*(self.hh*(d**2).sum(1)).sum()+(w*V).sum()
        g=np.zeros_like(Y)
        g[:-1]-=d; g[1:]+=d
        g+= (w[:,None]*ep)@self.A
        eh=np.exp(np.clip(psi[-1]/2,-200,200))
        E+=2*(self.C*(eh-1)).sum()
        g[-1]+= (self.C*eh)@self.A
        return E,g[1:].ravel()
    def phi0(self,v): return self.unpack(v)[-1]
    def psi0(self,v): return self.A@self.phi0(v)
    def run(self,v0,maxiter=20000):
        r=minimize(self.E_grad,v0,jac=True,method='L-BFGS-B',options=dict(maxiter=maxiter,maxfun=maxiter*2,ftol=1e-15,gtol=1e-9,maxcor=40))
        return r
    def random_start(self,rng,amp=2.0,decay=None):
        x=self.x[1:]; n=self.n
        k=rng.choice([0.5,1,2]) if decay is None else decay
        prof=np.exp(k*x)[:,None]
        a=rng.normal(size=n)*amp*rng.choice([0.2,1,3])
        return (prof*a[None,:]).ravel()

def E_book_from_BP(n,Cbp,C): return -2*(n+1)*Cbp-2*np.sum(C)

if __name__=='__main__':
    import sys
    n=2; P=Prob(n,(1,-1,-1)); rng=np.random.default_rng(0)
    v=P.random_start(rng)
    # gradient check
    e0,g=P.E_grad(v); i=7; dv=np.zeros_like(v); dv[i]=1e-6
    print('gradcheck',(P.E_grad(v+dv)[0]-P.E_grad(v-dv)[0])/2e-6,g[i], P.N)
