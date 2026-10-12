from halfline_min import *
import time
n=2; P=Prob(n,(1,-1,-1)); rng=np.random.default_rng(0)
v=P.random_start(rng)
e0,g=P.E_grad(v)
for i in (len(v)-1,len(v)-2,len(v)-30):
    dv=np.zeros_like(v); dv[i]=1e-6
    print((P.E_grad(v+dv)[0]-P.E_grad(v-dv)[0])/2e-6,g[i])
# a_2 all patterns
for C in itertools.product([1,-1],repeat=3):
    P=Prob(2,C); res=[]
    for s in range(12):
        rng=np.random.default_rng(s); t=time.time()
        r=P.run(P.random_start(rng)); res.append((r.fun,np.round(P.psi0(r.x),2)))
    res.sort(key=lambda z:z[0])
    print(C, 'prod',np.prod(C), 'best',res[0][0], res[0][1],'worst',res[-1][0], 'time',time.time()-t)
