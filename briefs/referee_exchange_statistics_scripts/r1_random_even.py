# Independent check of prp-exchange-twist with random grading-preserving R's.
# Graded embedding built from graded adjacent swaps (super tensor product), not from the book's sign formula.
import numpy as np, itertools, sys
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
def basis(d,N): return list(itertools.product(range(d),repeat=N))
def op_from_fn(d,N,f):
    B=basis(d,N); idx={b:i for i,b in enumerate(B)}; M=np.zeros((d**N,d**N),complex)
    for i,b in enumerate(B):
        for (c,amp) in f(b): M[idx[c],i]+=amp
    return M
def swap(d,N,k,g,graded):
    def f(b):
        c=list(b); c[k],c[k+1]=c[k+1],c[k]
        s=(-1)**(g[b[k]]*g[b[k+1]]) if graded else 1
        return [(tuple(c),s)]
    return op_from_fn(d,N,f)
def R2_embed12(d,N,R):  # R on factors 0,1 (no sign issue: adjacent, leftmost)
    return np.kron(R,np.eye(d**(N-2)))
def embed(d,N,R,j,g,graded):
    # move factor j to position 1 by adjacent swaps, act, move back
    S=np.eye(d**N)
    for k in range(j-1,0,-1): S=swap(d,N,k,g,graded)@S   # bring j -> 1
    return np.linalg.inv(S)@R2_embed12(d,N,R)@S
def rand_even(d,g):
    R=np.zeros((d*d,d*d),complex)
    for a1,a2,b1,b2 in itertools.product(range(d),repeat=4):
        if (g[a1]+g[a2]-g[b1]-g[b2])%2==0: R[b1*d+b2,a1*d+a2]=rng.normal()+1j*rng.normal()
    return R
def gg(d,g):
    return np.diag([(-1)**(g[b1]*g[b2]) for b1 in range(d) for b2 in range(d)])
worst=0
for d,g in [(2,[0,1]),(3,[0,1,0]),(3,[1,0,1]),(4,[0,1,1,0])]:
    for N in [2,3,4]:
        Rs=[rand_even(d,g) for _ in range(N)]
        T=np.eye(d**N); Tg=np.eye(d**N)
        for j in range(1,N):
            T=embed(d,N,Rs[j],j,g,False)@T
            Tg=embed(d,N,gg(d,g)@Rs[j],j,g,True)@Tg
        D=np.diag([(-1)**(g[b[0]]*(sum(g[x] for x in b)-1)) for b in basis(d,N)])
        r=np.abs(Tg-D@T).max(); rR=np.abs(Tg-T@D).max()
        worst=max(worst,r)
        print(d,g,N,'|Tgr - D T| =%.1e'%r,' |Tgr - T D| =%.1e'%rR)
print('worst',worst)
# second sentence: weights with (-1)^g = eps*exp(i pi gamma.mu); check Tgr = eps^(G-1) Omega_alpha T within G sector
