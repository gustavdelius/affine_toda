from tauvac import *
from scipy.optimize import minimize
import sys, time
XG=np.concatenate([np.linspace(-40,-12,60),np.linspace(-12,-0.05,2400)[:-1],-np.geomspace(0.05,1e-5,100),[0.0]])
def make(T,C,eta):
    def F(p):
        cs=cs_from_p(T.n,p); t,d,dd=T.tau(cs,np.array([0.0])); t0=t[:,0].real; d0=d[:,0].real
        eh=np.sqrt(np.abs(np.roll(t0,1)*np.roll(t0,-1)))/np.abs(t0)
        return -2*np.sum(d0/t0)+2*np.sum(np.array(C)*(eh-1))
    def g(p):
        cs=cs_from_p(T.n,p); t=T.tau(cs,XG)[0].real
        return t.min(axis=1)-eta                      # per node min over x
    return F,g
def random_feasible(T,rng,g,tries=2000):
    n=T.n; h=n+1
    for _ in range(tries):
        p=np.zeros(nparams(n))
        for a in range(1,(h-1)//2+1):
            p[2*a-2]=rng.uniform(-5,1.0); p[2*a-1]=rng.uniform(0,2*np.pi)
        if h%2==0: p[-1]=rng.choice([-1,1])*np.exp(rng.uniform(-5,0.7))
        if (g(p)>0).all(): return p
    return None
def bnds(n):
    h=n+1; b=[]
    for a in range(1,(h-1)//2+1): b+= [(-14,3.0),(None,None)]
    if h%2==0: b.append((-8,8))
    return b
def solve(n,C,nstart=12,seed=0,etas=(1e-1,1e-2,1e-3,1e-4)):
    T=Tau(n); rng=np.random.default_rng(seed); out=[]
    for s in range(nstart):
        F,g=make(T,C,etas[0]); p=random_feasible(T,rng,g)
        if p is None: continue
        ok=True
        for eta in etas:
            F,g=make(T,C,eta)
            if not (g(p)>0).all(): 
                # shrink toward feasibility: just continue; SLSQP handles
                pass
            r=minimize(F,p,method='SLSQP',constraints=[{'type':'ineq','fun':g}],bounds=bnds(n),options=dict(maxiter=300,ftol=1e-12))
            if (g(r.x)>-1e-9).all(): p=r.x
        F,g=make(T,C,etas[-1]); cs=cs_from_p(n,p); t=T.tau(cs,np.array([0.0]))[0][:,0].real
        out.append((F(p),p,t,(g(p)>-1e-7).all()))
    out=[o for o in out if o[3]]
    out.sort(key=lambda z:z[0]); return T,out
if __name__=='__main__':
    n=int(sys.argv[1]); C=tuple(int(c) for c in sys.argv[2].split(','))
    t0=time.time(); T,out=solve(n,C,int(sys.argv[3]) if len(sys.argv)>3 else 12)
    for F,p,t,_ in out[:5]: print('%.5f'%F,np.round(t,5),np.round(p,3))
    print('best of',len(out),'time',time.time()-t0)
XV=np.concatenate([-np.geomspace(50,1e-7,200000),[0.0]])
def validate(T,p):
    cs=cs_from_p(T.n,p); t=T.tau(cs,XV)[0].real; return t.min()
