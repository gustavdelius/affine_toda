from tauvac import *
from scipy.optimize import minimize
import sys
XG=np.concatenate([np.linspace(-30,-3,120),-np.geomspace(3,1e-4,160)[1:],[0.0]])
def make_obj(T,C,xg=XG,big=1e3):
    def obj(p):
        if np.any(np.abs(p[0::2][:len(p)//2*1])>25) : pass
        cs=cs_from_p(T.n,p)
        F,W,psi0,feas,t0=F_of(T,cs,C,xg)
        if not feas or not np.isfinite(F): return big
        return F
    return obj
def search(n,C,nstart=60,seed=0,xg=XG):
    T=Tau(n); np_=nparams(n); h=n+1; rng=np.random.default_rng(seed); obj=make_obj(T,C,xg)
    best=[]
    for s in range(nstart):
        p=np.zeros(np_)
        for a in range(1,(h-1)//2+1):
            p[2*a-2]=rng.uniform(-4,1.5); p[2*a-1]=rng.uniform(0,2*np.pi)
        if h%2==0: p[-1]=rng.choice([-1,1])*np.exp(rng.uniform(-4,1.0))
        if obj(p)>=1e3: continue
        r=minimize(obj,p,method='Nelder-Mead',options=dict(xatol=1e-10,fatol=1e-12,maxiter=4000,maxfev=6000,adaptive=True))
        best.append((r.fun,r.x))
    best.sort(key=lambda z:z[0]); return T,best
if __name__=='__main__':
    n=int(sys.argv[1])
    for C in itertools.product([1,-1],repeat=n+1):
        T,best=search(n,C,nstart=30)
        if not best: print(C,'none feasible'); continue
        F,p=best[0]; cs=cs_from_p(n,p); Fv,W,psi0,feas,t0=F_of(T,cs,C,XG)
        print(C,'prod',np.prod(C),'F=%.5f'%F,'psi0',np.round(psi0,2),'tau0',np.round(t0,4),'n_good',sum(abs(b[0]-F)<1e-4 for b in best),'of',len(best))
