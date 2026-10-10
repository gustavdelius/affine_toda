# a_n^{(1)} worked example: fusion point, crossing point of V_1, consistency with
# the fused antisoliton multiplet, in ch.7 conventions (Delta(e)=e(x)k+1(x)e).
import numpy as np, itertools
from functools import reduce
rng=np.random.default_rng(1)
def kron(*a): return reduce(np.kron,a)
def vec_rep(N,q,x):
    E=lambda i,j: np.eye(N)[:,[i]]@np.eye(N)[[j],:]
    e=[x*E(N-1,0)]+[E(i,i+1) for i in range(N-1)]
    f=[E(0,N-1)/x]+[E(i+1,i) for i in range(N-1)]
    k=[np.diag([q**(-1 if a==0 else (1 if a==N-1 else 0)) for a in range(N)])]
    k+= [np.diag([q**(1 if a==i else (-1 if a==i+1 else 0)) for a in range(N)]) for i in range(N-1)]
    return e,f,k
def check_rel(e,f,k,q):
    N=len(e)
    for i in range(N):
        c=e[i]@f[i]-f[i]@e[i]-(k[i]-np.linalg.inv(k[i]))/(q-1/q)
        assert np.allclose(c,0)
        for j in range(N):
            if i!=j: assert np.allclose(e[i]@f[j]-f[j]@e[i],0)
def cop(reps,gen):  # iterated coproduct of generator index g on tensor product of reps
    n=len(reps)
    def D(which,i):
        tot=0
        for j in range(n):
            if which=='e':
                fac=[np.eye(len(r[0][0])) for r in reps[:j]]+[reps[j][0][i]]+[reps[m][2][i] for m in range(j+1,n)]
            else:
                fac=[np.linalg.inv(reps[m][2][i]) for m in range(j)]+[reps[j][1][i]]+[np.eye(len(r[0][0])) for r in reps[j+1:]]
            tot=tot+kron(*fac)
        return tot
    return D
def fused(N,q,xs):
    reps=[vec_rep(N,q,x) for x in xs]; n=len(xs)
    D=cop(reps,None)
    Es=[D('e',i) for i in range(N)]; Fs=[D('f',i) for i in range(N)]
    Ks=[kron(*[r[2][i] for r in reps]) for i in range(N)]
    # highest weight vector: weight eps_1+..+eps_n, killed by e_1..e_{N-1}
    idx=[t for t in itertools.product(range(N),repeat=n) if sorted(t)==list(range(n))]
    B=np.zeros((N**n,len(idx)),complex)
    for c,t in enumerate(idx):
        m=0
        for a in t: m=m*N+a
        B[m,c]=1
    M=np.vstack([Es[i]@B for i in range(1,N)])
    u,s,vh=np.linalg.svd(M); v=vh.conj().T[:,-1]
    assert s[-1]<1e-9, s
    hw=B@v
    # generate module by f_1..f_{N-1}
    vecs=[hw]
    while True:
        new=[]
        for w in vecs:
            for i in range(1,N):
                y=Fs[i]@w
                if np.linalg.norm(y)>1e-10: new.append(y)
        U,S,_=np.linalg.svd(np.array(vecs+new).T)
        r=int((S>1e-9).sum())
        nv=list(U[:,:r].T)
        if len(nv)==len(vecs): vecs=nv; break
        vecs=nv
    Wb=np.array(vecs).T
    P=Wb@np.linalg.pinv(Wb)
    inv=max(np.linalg.norm((np.eye(N**n)-P)@G@Wb) for G in Es+Fs)
    return Wb,Es,Fs,Ks,inv
def restrict(Wb,G): return np.linalg.pinv(Wb)@G@Wb
def cyc_trace_dualtype(e):  # tr(e_0 e_{N-1} ... e_1) for V^*-type reps
    N=len(e); M=e[0]
    for i in range(N-1,0,-1): M=M@e[i]
    return np.trace(M)
for N in (3,4,5):
    n=N-1; h=N
    q=np.exp(1j*rng.uniform(0.3,1.2)); x=np.exp(1j*rng.uniform(0,6))*1.3
    check_rel(*vec_rep(N,q,x),q)
    # fusion point: q-antisym subspace of V1(x1)xV1(x2) invariant?
    for zc,name in [(q**2,'q^2'),(q**-2,'q^-2'),(-q**2,'-q^2')]:
        _,_,_,_,inv=fused(N,q,[x*np.sqrt(zc),x/np.sqrt(zc)])
        print(f"N={N}: Lambda^2 invariant at x1/x2={name}: {inv<1e-8}")
    # dual of V1(x)
    e,f,k=vec_rep(N,q,x)
    ed=[(-e[i]@np.linalg.inv(k[i])).T for i in range(N)]
    t_dual=cyc_trace_dualtype(ed)
    # fused V_n inside V1(y a^{n-1}) x ... x V1(y a^{1-n}), a = +1/q or -1/q; spectral centre y
    for sgn in (+1,-1):
        a=sgn/q
        y=1.0+0.37j
        xs=[y*a**(n+1-2*j) for j in range(1,n+1)]
        Wb,Es,Fs,Ks,inv=fused(N,q,xs)
        er=[restrict(Wb,G) for G in Es]
        t_f=cyc_trace_dualtype(er)   # proportional to y (degree 1 in e_0)
        ratio=t_dual/(t_f/y)          # y such that V1(x)^* = Vn^fused(y)
        xc=ratio/x
        cands={'q^h':q**h,'q^-h':q**-h,'(-q)^h':(-q)**h,'(-q)^-h':(-q)**-h}
        m=[k_ for k_,v in cands.items() if abs(v-xc)<1e-8]
        print(f"N={N} a={"+1/q" if sgn>0 else "-1/q"}: fused-invariant {inv<1e-8}; V1(x)^* = Vn_fused(x*xc), xc matches {m}")
