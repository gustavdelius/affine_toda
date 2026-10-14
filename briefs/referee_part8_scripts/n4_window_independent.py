"""Independent check of the N = 4 window of chapter 48 (referee_part8_report.md, Finding 1): T_j built by tensordot
for all four solitons, a_2^(1) fundamental solitons (3 colours), q = exp(i mu).  Prints the largest ||s|-1| in the
three-colour sectors and in the others, over all T_j and 40 rapidity samples, and max |[T_1, T_j]|."""
import numpy as np, itertools
def Rt(n,x,q):  # Jimbo R, normalized, as a 4-index tensor R[a',b',a,b]
    R=np.zeros((n,n,n,n),complex)
    for i in range(n):
        for j in range(n):
            if i==j: R[i,i,i,i]=1
            else:
                R[i,j,i,j]=(x-1/x)/(x*q-1/(x*q))
                R[j,i,i,j]=(q-1/q)*(x if i<j else 1/x)/(x*q-1/(x*q))
    return R
def apply2(T,R,i,j,N):
    # T has shape (n,)*N x (n,)*N ; left-multiply by R acting on factors i,j
    n=R.shape[0]; T=T.reshape((n,)*N+(-1,))
    T=np.moveaxis(T,[i,j],[0,1])
    T=np.tensordot(R,T,axes=([2,3],[0,1]))
    T=np.moveaxis(T,[0,1],[i,j])
    return T.reshape(n**N,-1)
def Tj(n,N,xs,q,j):
    T=np.eye(n**N,dtype=complex)
    order=[(j+k)%N for k in range(1,N)]
    for k in order: T=apply2(T,Rt(n,xs[j]/xs[k],q),j,k,N)
    return T
rng=np.random.default_rng(42)
n,N=3,4
sec={}
for idx in itertools.product(range(n),repeat=N):
    sec.setdefault(tuple(sorted(idx)),[]).append(np.ravel_multi_index(idx,[n]*N))
for mu in [1.0,2.0,2.09,2.0944,2.0954,2.1,2.3,2.6,2.9,3.1]:
    q=np.exp(1j*mu); w3=0; w2=0; comm=0
    for t in range(40):
        xs=np.exp(rng.normal(size=N)*1.5)
        Ts=[Tj(n,N,xs,q,j) for j in range(N)]
        comm=max(comm,max(np.abs(Ts[0]@Ts[j]-Ts[j]@Ts[0]).max() for j in range(1,N)))
        for key,ids in sec.items():
            for T in Ts:
                d=np.abs(np.abs(np.linalg.eigvals(T[np.ix_(ids,ids)]))-1).max()
                if len(set(key))==3: w3=max(w3,d)
                else: w2=max(w2,d)
    print(f"mu={mu:.4f} Re q={np.cos(mu):+.3f}: max dev 3-colour {w3:.1e}, <=2 colours {w2:.1e}, max|[T1,Tj]| {comm:.1e}")
