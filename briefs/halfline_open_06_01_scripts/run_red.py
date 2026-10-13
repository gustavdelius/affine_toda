import sys
from redmin import *
n=int(sys.argv[1]); ns=int(sys.argv[2]); allp=len(sys.argv)>3
seen=set()
for C in itertools.product([1,-1],repeat=n+1):
    orb=min(tuple(np.roll(C,k)) for k in range(n+1))
    if not allp and orb in seen: continue
    seen.add(orb)
    T,best=search(n,C,nstart=ns)
    if not best: print(C,'none feasible',flush=True); continue
    F,p=best[0]; cs=cs_from_p(n,p); Fv,W,psi0,feas,t0=F_of(T,cs,C,XG)
    print(C,'prod',np.prod(C),'F=%.5f'%F,'psi0',np.round(psi0,1),'tau0',np.round(t0,4),'|c|',np.round(np.abs(cs[:(n+1)//2]),3),'arg',np.round(np.angle(cs[:(n+1)//2]),3),'mid',np.round(cs[(n+1)//2-1].real,3) if (n+1)%2==0 else '', 'n_good',sum(abs(b[0]-F)<1e-4 for b in best),'of',len(best),flush=True)
