# Jimbo R-matrix of U_q(sl_3^) vector rep, solved numerically from the intertwining equations,
# for the ch.7 coproduct and for the opposite coproduct (Part VI).  Eigenvalue on Lambda^2.
import numpy as np
exec(open('an_cross.py').read().split('def fused')[0])
N=3; q=np.exp(0.83j); rng=np.random.default_rng(5)
def Rcheck(z,op=False):
    A=vec_rep(N,q,z); B=vec_rep(N,q,1.0)
    def D(r1,r2,kind,i):
        e1,f1,k1=r1; e2,f2,k2=r2; I=np.eye(N)
        if not op:
            return np.kron(e1[i],k2[i])+np.kron(I,e2[i]) if kind=='e' else np.kron(f1[i],I)+np.kron(np.linalg.inv(k1[i]),f2[i])
        else:  # opposite coproduct: Delta^op(e)=k(x)e + e(x)1
            return np.kron(k1[i],e2[i])+np.kron(e1[i],I) if kind=='e' else np.kron(I,f2[i])+np.kron(f1[i],np.linalg.inv(k2[i]))
    rows=[]
    for kind in 'ef':
        for i in range(N):
            L=D(A,B,kind,i); Rr=D(B,A,kind,i)   # Rc L = Rr Rc
            rows.append(np.kron(np.eye(N*N),L.T)-np.kron(Rr,np.eye(N*N)))
    M=np.vstack(rows); u,s,vh=np.linalg.svd(M)
    Rc=vh.conj()[-1].reshape(N*N,N*N)
    return Rc/Rc[0,0], s[-2:]
for op in (False,True):
    z=1.7*np.exp(0.4j)
    Rc,s=Rcheck(z,op)
    ev=np.linalg.eigvals(Rc)
    pred=(1-q**2*z)/(z-q**2) if not op else (1-q**-2*z)/(z-q**-2)
    print("opposite coproduct" if op else "ch.7 coproduct", "| null-space sing. vals",np.round(s,12),
          "| eigenvalues:",np.round(sorted(ev,key=lambda w:abs(w-1)),6)," predicted rho:",np.round(pred,6))
