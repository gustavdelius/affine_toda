"""I+ vector solution; source convention checks; separately labelled R comparison."""
import itertools
import sys
from pathlib import Path
import numpy as np
from independent_ar import weights,roots,sector_generator,solve_sector,rcheck,triangle


class PlusVector:
    def __init__(self,h,q,y,mu,s):
        self.h,self.q,self.y,self.mu=h,q,y,mu
        self.ww=weights(h);self.base=s+self.ww[0]
        self.zero,sv=solve_sector(self.base,q,y,mu,plus=True)
        self.sv=sv;self.cache={tuple([0]*h):self.zero}
    def __call__(self,Q):
        delta=np.round(Q-self.base).astype(int);key=tuple(delta)
        assert np.max(abs(Q-self.base-delta))<1e-7
        if key in self.cache:return self.cache[key]
        curr=self.base.copy();M=self.zero.copy()
        for a,n in enumerate(np.cumsum(delta)[:-1]):
          for _ in range(abs(n)):
            nxt=curr+np.sign(n)*roots(self.h)[a]
            P=curr if n>0 else nxt
            A=sector_generator(P,a,self.q,self.y,self.mu,plus=True)
            B=sector_generator(P,a,self.q,self.y,self.mu,plus=True,target=True)
            M=B@M@np.linalg.inv(A) if n>0 else np.linalg.solve(B,M@A)
            curr=nxt
        self.cache[key]=M
        return M
    def trans(self,i,N):
        M=self(N+self.ww[i])
        return [(j,N+self.ww[i]-self.ww[j],M[j,i]) for j in range(self.h)]


def plus_triangle(h,q,y1,y2,mu,s):
    A=PlusVector(h,q,y1,mu,s);B=PlusVector(h,q,y2,mu,s)
    R=rcheck(h,q,y1,y2).reshape(h,h,h,h);err=0
    for i,j in itertools.product(range(h),repeat=2):
        lhs={};rhs={}
        def add(d,k,v):d[k]=d.get(k,0)+v
        for b,N1,t2 in B.trans(j,s):
          for a,N2,t1 in A.trans(i,N1):
            for c,d in itertools.product(range(h),repeat=2):
                add(lhs,(c,d,*np.round(N2,8)),R[c,d,a,b]*t1*t2)
        for b,a in itertools.product(range(h),repeat=2):
          for d,N1,t1 in A.trans(a,s):
            for c,N2,t2 in B.trans(b,N1):
                add(rhs,(c,d,*np.round(N2,8)),t2*t1*R[b,a,i,j])
        scale=max(map(abs,lhs.values()))
        err=max(err,max(abs(lhs.get(k,0)-rhs.get(k,0)) for k in lhs.keys()|rhs.keys())/scale)
    return err


def run():
    w=1.37;q=-np.exp(-1j*np.pi*w);mu=-np.exp(-1j*np.pi*w/2)/(q-1/q)
    for h in [3,4]:
        s=np.arange(h)*.123;s-=s.mean();A=PlusVector(h,q,1.2,mu,s);err=0
        for delta in [np.zeros(h),roots(h)[0],-roots(h)[-1]]:
          Q=A.base+delta
          for a in range(h):
            G=sector_generator(Q,a,q,1.2,mu,plus=True)
            H=sector_generator(Q,a,q,1.2,mu,plus=True,target=True)
            err=max(err,np.linalg.norm(A(Q+roots(h)[a])@G-H@A(Q)))
        print('I+ V h,nullity,entries/column,intertwining,triangle',h,sum(A.sv<1e-10),
              (abs(A.zero)>1e-9).sum(axis=0),err,plus_triangle(h,q,1.43,.87,mu,s))
    # Conventional weights give CZ2007 D as the transpose of the book's D.
    for h in [3,4,5]:
        lambdas=[np.zeros(h)]+[sum(np.eye(h)[:a])-a*np.ones(h)/h for a in range(1,h)]+[np.zeros(h)]
        D=2*sum(np.outer(lambdas[a]-lambdas[a+1],lambdas[a]) for a in range(1,h))
        DCZ=2*sum(np.outer(lambdas[a],lambdas[a]-lambdas[a+1]) for a in range(1,h))
        x=np.exp(-.38);zerr=0
        for a in range(1,h):
            cz=(x+1j*np.exp(-1j*np.pi*a/h))/(x+1j*np.exp(1j*np.pi*a/h))
            book=(x-1j*np.exp(1j*np.pi*(h-a)/h))/(x-1j*np.exp(-1j*np.pi*(h-a)/h))
            zerr=max(zerr,abs(cz-book))
        print('CZ2007 h,D transpose and delay relabelling residual',h,np.linalg.norm(DCZ-D.T),zerr)
    # Book import ONLY for a comparison with an independently constructed R.
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'textbook_code'))
    import boundary_common as bc
    for h in [3,4]:
      for b1,b2 in [(False,False),(True,True),(False,True)]:
        A=bc.an_vector(h-1,q,1.43,b1);B=bc.an_vector(h-1,q,.87,b2)
        ref,nul=bc.rcheck(A,B);own=rcheck(h,q,1.43,.87,b1,b2)
        factor=np.vdot(own,ref)/np.vdot(own,own)
        print('COMPARISON ONLY PartVI R h,b1,b2,error',h,b1,b2,np.linalg.norm(ref-factor*own)/np.linalg.norm(ref))

if __name__=='__main__':run()
