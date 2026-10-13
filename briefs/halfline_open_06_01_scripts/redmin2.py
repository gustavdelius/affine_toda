"""Reduced (Bogomolny) minimization over the static all-species Hirota family: F_P(c) = W(phi0)+B_P(phi0),
   feasibility = all tau_j > 0 on x<=0 (dense grid + parabola refinement of local minima).  Infimum may sit at the boundary of the feasible region."""
from tauvac import *
from scipy.optimize import minimize
import sys, time, multiprocessing as mp, pickle
XG=np.unique(np.concatenate([np.linspace(-45,-14,60),np.linspace(-14,-0.3,1800)[:-1],-np.geomspace(0.3,1e-6,300),[0.0]]))
def tau_min(T,cs):
    t=T.tau(cs,XG)[0].real
    m=t.min(axis=1).copy()
    y0,y1,y2=t[:,:-2],t[:,1:-1],t[:,2:]
    x0,x1,x2=XG[:-2],XG[1:-1],XG[2:]
    den=(x0-x1)*(x0-x2)*(x1-x2)
    A=(x2*(y1-y0)+x1*(y0-y2)+x0*(y2-y1))/den
    B=(x2**2*(y0-y1)+x1**2*(y2-y0)+x0**2*(y1-y2))/den
    Cc=(x1*x2*(x1-x2)*y0+x2*x0*(x2-x0)*y1+x0*x1*(x0-x1)*y2)/den
    with np.errstate(all='ignore'):
        xm=-B/(2*A); ym=A*xm**2+B*xm+Cc
    ok=(A>0)&(xm>x0)&(xm<x2)&(y1<=y0)&(y1<=y2)
    ym=np.where(ok,ym,np.inf)
    return np.minimum(m,ym.min(axis=1))
def Fval(T,cs,C):
    t,d,dd=T.tau(cs,np.array([0.0])); t0=t[:,0].real; d0=d[:,0].real
    eh=np.sqrt(np.abs(np.roll(t0,1)*np.roll(t0,-1)))/np.abs(t0)
    return -2*np.sum(d0/t0)+2*np.sum(np.array(C)*(eh-1)), t0, eh
BIG=1e3
def obj_factory(T,C):
    n=T.n
    def obj(p):
        if np.any(np.abs(p)>40): return BIG*10
        with np.errstate(all='ignore'):
            cs=cs_from_p(n,p)
            m=tau_min(T,cs)
            if not np.all(np.isfinite(m)): return BIG*10
            if m.min()<=1e-6: return BIG+BIG*(1e-6-m.min())
            F,t0,eh=Fval(T,cs,C)
        return F if np.isfinite(F) else BIG*10
    return obj
def one_start(args):
    n,C,seed=args
    T=Tau(n); h=n+1; rng=np.random.default_rng(seed); obj=obj_factory(T,C)
    for _ in range(400):
        p=np.zeros(nparams(n))
        for a in range(1,(h-1)//2+1):
            p[2*a-2]=rng.uniform(-5,1.5); p[2*a-1]=rng.uniform(0,2*np.pi)
        if h%2==0: p[-1]=rng.choice([-1,1])*np.exp(rng.uniform(-5,1.0))
        if obj(p)<BIG: break
    else: return None
    r=minimize(obj,p,method='Nelder-Mead',options=dict(xatol=1e-9,fatol=1e-11,maxiter=2500,maxfev=3000,adaptive=True))
    return (r.fun,r.x)
def job(args):
    n,C,ns=args
    res=[one_start((n,C,s)) for s in range(ns)]
    res=[r for r in res if r is not None and r[0]<BIG]; res.sort(key=lambda z:z[0])
    return (n,C,res)
if __name__=='__main__':
    n=int(sys.argv[1]); ns=int(sys.argv[2])
    seen=set(); jobs=[]
    for C in itertools.product([1,-1],repeat=n+1):
        orb=min(tuple(np.roll(C,k)) for k in range(n+1))
        if orb in seen: continue
        seen.add(orb); jobs.append((n,C,ns))
    out=[]
    with mp.Pool(6) as pool:
        for n_,C,res in pool.imap_unordered(job,jobs):
            T=Tau(n_)
            if not res: print(C,'none',flush=True); continue
            F,p=res[0]; cs=cs_from_p(n_,p); Fv,t0,eh=Fval(T,cs,C)
            print(C,'prod',int(np.prod(C)),'F=%.5f'%F,'tau0',np.round(t0,4),'n_near',sum(abs(r[0]-F)<1e-3 for r in res),'of',len(res),'p',np.round(p,3),flush=True)
            out.append((C,res[:5]))
    pickle.dump(out,open('red2_%d.pkl'%n,'wb'))
