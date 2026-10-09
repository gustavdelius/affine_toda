import numpy as np, itertools
from krein_bethe import R, embed
rng=np.random.default_rng(11)
for n in [3,4]:
  for N in [3,4] if n==3 else [3]:
    for mu in [0.5,1.5,2.5]:
        q=np.exp(1j*mu); worst={}
        for t in range(25):
            x=np.exp(rng.normal(size=N)*2)
            Tm=np.eye(n**N,dtype=complex)
            for j in range(1,N): Tm=embed(R(n,x[0]/x[j],q),0,j,n,N)@Tm
            # sectors labelled by sorted colour content (weight conservation)
            sectors={}
            for idx in itertools.product(range(n),repeat=N):
                sectors.setdefault(tuple(sorted(idx)),[]).append(np.ravel_multi_index(idx,[n]*N))
            for key,ids in sectors.items():
                ev=np.linalg.eigvals(Tm[np.ix_(ids,ids)])
                ncol=len(set(key)); worst[ncol]=max(worst.get(ncol,0),np.abs(np.abs(ev)-1).max())
        print(f"n={n} N={N} mu={mu}: max||s|-1| by number of distinct colours:", {k:float(f'{v:.2e}') for k,v in sorted(worst.items())})
