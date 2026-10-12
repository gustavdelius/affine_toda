from halfline_min import *
import sys,time
n=int(sys.argv[1]); pats=[tuple(int(c) for c in s.split(',')) for s in sys.argv[2:]]
for C in pats:
    P=Prob(n,C,dx0=1e-3); res=[]
    for s in range(8):
        r=P.run(P.random_start(np.random.default_rng(s)),maxiter=6000); res.append((r.fun,np.round(P.psi0(r.x),1)))
    res.sort(key=lambda z:z[0]); print(C,'best %.4f'%res[0][0],res[0][1],'worst %.4f'%res[-1][0],flush=True)
