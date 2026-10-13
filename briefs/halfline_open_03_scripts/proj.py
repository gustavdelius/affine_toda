import numpy as np
from nk_common import k_book, fast_rcheck
from boundary_common import q_phys, an_vector, cop, nodes_of, kron
def singlet(A,B):
    d=A['dim']*B['dim']
    M=np.vstack([cop(A,B,kd,j)-(np.eye(d) if kd=='k' else 0) for j in nodes_of(A) for kd in 'efk'])
    _,sv,Vh=np.linalg.svd(M); w=Vh[-1].conj()
    return sv[-1]/sv[0], w/w[int(np.argmax(np.abs(w)>1e-8))]
def cproj(n,s,om,th,ca=False,diag=False):
    """matrix part of the B1 reflection: (1(x)K(th_a)) R (1(x)K(th_b)) on the singlet; returns (coef, residual)."""
    q=q_phys(om); h=n+1; u=np.pi*(1-2/(h*om)); N=n+1
    ta,tb=th-0.5j*u, th+0.5j*u
    Va=an_vector(n,q,np.exp(om*ta),ca); Vb=an_vector(n,q,np.exp(om*tb),not ca)
    r_in,w=singlet(Va,Vb)
    r_out,wo=singlet(an_vector(n,q,np.exp(-om*tb),ca),an_vector(n,q,np.exp(-om*ta),not ca))
    KV=lambda t: np.eye(N) if diag else k_book(n,om,t,s)
    KW=lambda t: np.eye(N) if diag else k_book(n,om,t,s).T
    Kb = KW(tb) if not ca else KV(tb)     # second constituent reflects first
    Ka = KV(ta) if not ca else KW(ta)
    v1=kron(np.eye(N),Kb)@w
    R=fast_rcheck(Va,an_vector(n,q,np.exp(-om*tb),ca))
    v=kron(np.eye(N),Ka)@(R@v1)
    c=np.vdot(wo,v)/np.vdot(wo,wo)
    return c, np.linalg.norm(v-c*wo)/np.linalg.norm(v), max(r_in,r_out)
if __name__=='__main__':
    for n in (2,3,4,5):
        for s in (1,-1):
            for om in (0.37,1.83):
                for th in (0.31+0.07j,):
                    c1,r1,rs=cproj(n,s,om,th); c2,r2,_=cproj(n,s,om,th,ca=True)
                    print(n,s,om,'c=',np.round(c1,8),'res',f'{r1:.1e}','other order',np.round(c2,8),f'{r2:.1e}',f'{rs:.0e}')
