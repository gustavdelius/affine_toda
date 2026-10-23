"""Krein-unitarity of Bethe-Yang transfer matrices built from Jimbo's U_q(sl_n^) R-matrix, |q|=1, real x.
Normalized Rt(x) = R(x)/(xq - 1/(xq)).  Claims:
 (1) Rt(y)^dag = Rt(1/conj y) for all complex y;
 (2) (C x C) Rt_12 (C x C) = Rt_21, C = colour reversal;
 (3) unitarity Rt_12(x) Rt_21(1/x) = 1;
 => T_1 = Rt_{1N}(x_1/x_N) ... Rt_{12}(x_1/x_2) satisfies T^dag = C^{(N)} T^{-1} C^{(N)}."""
import numpy as np, itertools
rng=np.random.default_rng(0)
def R(n,x,q):
    d=n*n; M=np.zeros((d,d),complex)
    for i in range(n):
        for j in range(n):
            if i==j: M[i*n+i,i*n+i]=x*q-1/(x*q)
            else:
                M[i*n+j,i*n+j]=x-1/x
                M[i*n+j,j*n+i]=(q-1/q)*(x if i<j else 1/x)
    return M/(x*q-1/(x*q))
def swap(n):
    P=np.zeros((n*n,n*n))
    for i in range(n):
        for j in range(n): P[i*n+j,j*n+i]=1
    return P
def embed(M,i,j,n,N):
    # act with two-site matrix M on sites i,j of V^{(x)N} (M indexed (a_i,a_j))
    D=n**N; out=np.zeros((D,D),complex)
    for idx in itertools.product(range(n),repeat=N):
        col=np.ravel_multi_index(idx,[n]*N)
        a,b=idx[i],idx[j]
        for a2 in range(n):
            for b2 in range(n):
                v=M[a2*n+b2,a*n+b]
                if v!=0:
                    l=list(idx); l[i]=a2; l[j]=b2
                    out[np.ravel_multi_index(l,[n]*N),col]+=v
    return out
for n in [2,3,4]:
    mu=rng.uniform(0.2,3.0); q=np.exp(1j*mu); P=swap(n)
    C=np.eye(n)[::-1]; CC=np.kron(C,C)
    e1=max(np.abs(R(n,y,q).conj().T-R(n,1/np.conj(y),q)).max() for y in rng.normal(size=5)+1j*rng.normal(size=5))
    x=np.exp(rng.normal()); e2=np.abs(CC@R(n,x,q)@CC-P@R(n,x,q)@P).max()
    e3=np.abs(R(n,x,q)@(P@R(n,1/x,q)@P)-np.eye(n*n)).max()
    print(f"n={n} mu={mu:.3f}: hermiticity {e1:.1e}, colour reversal {e2:.1e}, unitarity {e3:.1e}")
    for N in [3,4] if n<4 else [3]:
        th=rng.normal(size=N)*1.5; T=2.0; xs=np.exp(2*T*th)
        Tm=np.eye(n**N,dtype=complex)
        for j in range(1,N): Tm=embed(R(n,xs[0]/xs[j],q),0,j,n,N)@Tm
        CN=C
        for _ in range(N-1): CN=np.kron(CN,C)
        ek=np.abs(Tm.conj().T-CN@np.linalg.inv(Tm)@CN).max()
        ev=np.linalg.eigvals(Tm); mods=np.sort(np.abs(ev))
        # pairing check: multiset of |s| invariant under |s|->1/|s|, and s-> 1/conj s
        pair=max(min(abs(e-1/np.conj(f)) for f in ev) for e in ev)
        print(f"   N={N}: Krein T^dag=C T^-1 C err {ek:.1e}; max||s|-1| {np.abs(mods-1).max():.3f}; pairing err {pair:.1e}")
