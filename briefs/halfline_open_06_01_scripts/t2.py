from halfline_min import *
C=(1,1,-1)
for dx0 in (1e-2,1e-3,1e-4,1e-5,1e-6,1e-7):
    P=Prob(2,C,dx0=dx0); best=None
    for s in range(6):
        r=P.run(P.random_start(np.random.default_rng(s)))
        if best is None or r.fun<best.fun: best=r
    print(dx0,P.N,best.fun,np.round(P.psi0(best.x),2), best.nit)
print(-6*(3**.5-1)-2)
