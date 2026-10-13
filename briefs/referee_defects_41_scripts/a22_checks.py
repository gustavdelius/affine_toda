"""Tzitzeica construction directly from CZ2011 (4.5),(4.6),(4.9).

No textbook imports and no finite truncation of the defect module.
"""
import itertools
import numpy as np
from independent_ar import tensor,finite_hom

W=np.array([1,0,-1])


def vec(q,x):
    K1=np.diag([q,1,1/q]); K0=np.diag([q**-2,1,q**2])
    E1=np.diag(np.ones(2),1).astype(complex); F1=(q+1/q)*E1.T
    E0=np.zeros((3,3),complex); E0[2,0]=x; F0=E0.T/x**2
    return {('k',0):K0@K0,('k',1):K1@K1,
            ('e',0):K0@E0,('e',1):K1@E1,
            ('f',0):F0@np.linalg.inv(K0),('f',1):F1@np.linalg.inv(K1)}


def F(n,q,a,b,control=None):
    if control=='one': return 1
    if control=='sign': return (-1.)**n
    return a*q**(-2*n-1)+b*q**(2*n+1)


def we(n,j,q,a,b,control=None):
    if j==1: return n+1,q**(n+1)
    return n-2,q**(-2*(n-2))*F(n,q,a,b,control)*F(n-1,q,a,b,control)


def qb(n,k,q):
    def fact(n): return np.prod([(q**i-q**-i)/(q-1/q) for i in range(1,n+1)])
    return fact(n)/(fact(k)*fact(n-k))


def serre(q,a,b,control=None):
    errs=[]
    for j,k,m in [(0,1,2),(1,0,5)]:
      for n in range(-3,4):
        terms=[]
        for l in range(m+1):
            P=n; c=1
            for g in reversed([j]*(m-l)+[k]+[j]*l):
                P,t=we(P,g,q,a,b,control); c*=t
            terms.append((-1)**l*qb(m,l,q**(4 if j==0 else 1))*c)
        errs.append(abs(sum(terms))/max(1,sum(map(abs,terms))))
    return max(errs)


def gen(Q,j,q,x,a,b,target=False):
    v=vec(q,x); A=np.zeros((3,3),complex)
    for i in range(3):
        n=Q-W[i]; _,c=we(n,j,q,a,b)
        kw=q**((-4 if j==0 else 2)*n)
        A[:,i]+=v['e',j][:,i]*(kw if target else 1)
        A[i,i]+=c*(1 if target else v['k',j][i,i])
    return A


def seed(q,x,a,b):
    L=[[],[]]
    for target in [0,1]:
      for order in [(0,1,1),(1,0,1),(1,1,0)]:
        P=0; M=np.eye(3,dtype=complex)
        for j in order:
            M=gen(P,j,q,x,a,b,target)@M; P+=(-2 if j==0 else 1)
        L[target].append(M)
    basis=np.eye(9).reshape(9,3,3)
    eq=[]
    for A,B in zip(*L):
        C=np.einsum('ij,njk->nik',B,basis)-np.einsum('nij,jk->nik',basis,A)
        eq.append(C.reshape(9,-1).T)
    _,sv,vh=np.linalg.svd(np.concatenate(eq),full_matrices=False)
    return vh[-1].conj().reshape(3,3),sv


class T:
    def __init__(self,q,x,a,b):
        self.q,self.x,self.a,self.b=q,x,a,b
        M,self.sv=seed(q,x,a,b); self.cache={0:M}
    def __call__(self,Q):
        Q=int(Q)
        if Q not in self.cache:
            p=Q-1 if Q>0 else Q+1
            M=self(p)
            if Q>0:
                A=gen(p,1,self.q,self.x,self.a,self.b)
                B=gen(p,1,self.q,self.x,self.a,self.b,True)
                self.cache[Q]=B@M@np.linalg.inv(A)
            else:
                A=gen(Q,1,self.q,self.x,self.a,self.b)
                B=gen(Q,1,self.q,self.x,self.a,self.b,True)
                self.cache[Q]=np.linalg.solve(B,M@A)
        return self.cache[Q]
    def trans(self,i,n):
        M=self(n+W[i])
        return [(j,n+W[i]-W[j],M[j,i]) for j in range(3)]


def triangle(q,x1,x2,a,b):
    v1=vec(q,x1); v2=vec(q,x2)
    ww=np.array([[i+j] for i in W for j in W])
    R,sv=finite_hom(tensor(v1,v2),tensor(v2,v1),ww,ww)
    R=R.reshape(3,3,3,3)
    A=T(q,x1,a,b); B=T(q,x2,a,b); error=0
    for n in [-2,0,2]:
      for i,j in itertools.product(range(3),repeat=2):
        lhs={};rhs={}
        def add(d,k,c):d[k]=d.get(k,0)+c
        for l,n1,t2 in B.trans(j,n):
          for k,n2,t1 in A.trans(i,n1):
            for c,d in itertools.product(range(3),repeat=2):
              add(lhs,(c,d,n2),R[c,d,k,l]*t1*t2)
        for l,k in itertools.product(range(3),repeat=2):
          for d,n1,t1 in A.trans(k,n):
            for c,n2,t2 in B.trans(l,n1):
              add(rhs,(c,d,n2),t2*t1*R[l,k,i,j])
        scale=max(max(map(abs,lhs.values())),1e-12)
        error=max(error,max(abs(lhs.get(k,0)-rhs.get(k,0)) for k in lhs.keys()|rhs.keys())/scale)
    return error,sv[-2:]


def crossing(q,x,a,b,C):
    xc=-q**-6; A=T(q,x,a,b); B=T(q,x/xc,a,b)
    lhs=[];rhs=[]
    for Q in range(-2,3):
      M=A(Q)
      for j,i in itertools.product(range(3),repeat=2):
        n=Q-W[i]; ai=2-i; bj=2-j
        inv=np.linalg.inv(B(n+W[bj]))
        # The crossed representation's LEFT dual is the original one:
        # C maps V(x/xc)^* to V(x), hence C U^{-t_V} C^{-1}.
        lhs.append(M[j,i]); rhs.append(inv[ai,bj]*C[j,bj]/C[i,ai])
    lhs=np.array(lhs);rhs=np.array(rhs)
    c=np.vdot(rhs,lhs)/np.vdot(rhs,rhs)
    return np.linalg.norm(lhs-c*rhs)/np.linalg.norm(lhs),c


def run():
    q=np.exp(-.77j);x=1.31
    C=np.array([[0,0,1],[0,-1/q,0],[q**-2,0,0]],complex)
    v=vec(q,x); vdual={}
    for j in range(2):
        vdual['k',j]=np.linalg.inv(v['k',j]).T
        vdual['e',j]=(-np.linalg.inv(v['k',j])@v['e',j]).T
        vdual['f',j]=(-v['f',j]@v['k',j]).T
    vc=vec(q,x*(-q**-6))
    print('dual C residual',max(np.linalg.norm(C@vdual[key]-vc[key]@C) for key in v))
    cartan=np.array([[2,-1],[-4,2]]);err=0
    for i,j in itertools.product(range(2),repeat=2):
        qi=q**(4 if i==0 else 1)
        for g,s in [('e',1),('f',-1)]:
            err=max(err,np.linalg.norm(v['k',i]@v[g,j]@np.linalg.inv(v['k',i])-qi**(s*cartan[i,j])*v[g,j]))
        comm=v['e',i]@v['f',j]-v['f',j]@v['e',i]
        rhs=(v['k',i]-np.linalg.inv(v['k',i]))/(qi-1/qi) if i==j else 0
        err=max(err,np.linalg.norm(comm-rhs))
    print('finite representation relations',err)
    for a,b in [(1,.7),(1,0),(0,1),(.4+.3j,-1.1)]:
        M=T(q,x,a,b);er=0
        for Q in range(-5,6):
          for j in range(2):
            A=gen(Q,j,q,x,a,b);B=gen(Q,j,q,x,a,b,True)
            lhs=M(Q+(-2 if j==0 else 1))@A;rhs=B@M(Q)
            er=max(er,np.linalg.norm(lhs-rhs)/max(np.linalg.norm(lhs),1e-12))
        print('F params, Serre,loop nullity/gap,intertwining',a,b,serre(q,a,b),sum(M.sv<1e-10),M.sv[-2],er)
        print('triangle',triangle(q,1.13,.74,a,b))
        print('crossing',crossing(q,x,a,b,C))
        qm=-q;Cm=np.array([[0,0,1],[0,-1/qm,0],[qm**-2,0,0]],complex)
        print('minus q own C / wrong C',crossing(qm,x,a,b,Cm)[0],crossing(qm,x,a,b,C)[0])
    print('Serre failing controls',serre(q,1,.7,'one'),serre(q,1,.7,'sign'))


if __name__=='__main__':run()
