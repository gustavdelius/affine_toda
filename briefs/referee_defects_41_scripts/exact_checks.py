"""Exact finite-field rank certificates and convention/exercise checks.

A nonzero minor mod 1009 is a nonzero rational minor. The identity gives the
matching upper bound on the commutator rank. No defect-lattice truncation.
"""
import itertools
import numpy as np
from independent_ar import solve_sector,weights

P=1009
def pw(a,b):return pow(int(a),int(b),P)
def rank(A):
    A=np.array(A,dtype=np.int64)%P;r=0
    for c in range(A.shape[1]):
        piv=next((i for i in range(r,len(A)) if A[i,c]),None)
        if piv is None:continue
        A[[r,piv]]=A[[piv,r]];A[r]=A[r]*pw(A[r,c],-1)%P
        for i in range(r+1,len(A)): A[i]=(A[i]-A[i,c]*A[r])%P
        r+=1
        if r==len(A):break
    return r
def inv(A):
    n=len(A);A=np.concatenate([A%P,np.eye(n,dtype=np.int64)],axis=1)
    for c in range(n):
        j=next(j for j in range(c,n) if A[j,c])
        A[[c,j]]=A[[j,c]];A[c]=A[c]*pw(A[c,c],-1)%P
        for i in range(n):
            if i!=c:A[i]=(A[i]-A[i,c]*A[c])%P
    return A[:,n:]
def rep(h,y,bar=False):
    I=np.eye(h,dtype=np.int64);R={};sgn=-1 if bar else 1
    for a in range(h):
        b=(a+1)%h;E=np.zeros((h,h),np.int64);E[a,b]=1
        R['k',a]=np.diag([pw(2,sgn*k) for k in I[a]-I[b]])
        R['e',a]=y*(E.T if bar else E)%P
        R['f',a]=pw(y,-1)*(E if bar else E.T)%P
    return R,sgn*I
def tensor(A,B):
    d1=len(A['k',0]);d2=len(B['k',0]);R={}
    for _,j in [k for k in A if k[0]=='k']:
        R['k',j]=np.kron(A['k',j],B['k',j])%P
        R['e',j]=(np.kron(A['e',j],np.eye(d2,dtype=np.int64))+np.kron(A['k',j],B['e',j]))%P
        R['f',j]=(np.kron(A['f',j],inv(B['k',j]))+np.kron(np.eye(d1,dtype=np.int64),B['f',j]))%P
    return R
def certificate(h,bar=False,plus=False,two=False):
    V,wt=rep(h,5*pw(4,-1)%P,bar)
    if two:
        B,wb=rep(h,9*pw(5,-1)%P);V=tensor(V,B)
        wt=np.array([u+v for u in wt for v in wb])
    n=len(wt);Q=wt[0].copy();mu=3*pw(7,-1)%P
    al=np.eye(h,dtype=np.int64)-np.roll(np.eye(h,dtype=np.int64),-1,axis=0)
    def A(Q,a):
        out=np.zeros((n,n),np.int64)
        for i in range(n):
            N=Q-wt[i];kw=pw(2,np.dot(al[a],N));f=mu*pw(2,N[(a+1)%h])%P
            if plus:
                out[:,i]=V['e',a][:,i];out[i,i]+=f*V['k',a][i,i]
            else:
                out[:,i]=V['f',a][:,i]*pw(kw,-1);out[i,i]+=f
        return out%P
    equations=[];sgn=1 if plus else -1;I=np.eye(n,dtype=np.int64)
    for a,b in itertools.combinations(range(h),2):
        ab=A(Q+sgn*al[b],a)@A(Q,b)%P
        ba=A(Q+sgn*al[a],b)@A(Q,a)%P
        L=inv(ba)@ab%P
        equations.append((np.kron(L.T,I)-np.kron(I,L))%P)
    rr=rank(np.concatenate(equations))
    assert rr==n*n-1
    return n,rr,n*n-rr

def run():
    for h in [3,4]:
      for bar,plus in itertools.product([False,True],repeat=2):
        print('exact transported-loop commutant h,bar,plus: dim/rank/nullity',h,bar,plus,certificate(h,bar,plus))
    print('exact hypothesis(ii) h3 V tensor V: dim/rank/nullity',certificate(3,two=True))
    # Diagonal gauge for orientation -1: g(N)=q^(-N.N/2) makes the
    # off-diagonal entry q X. Thus the same CZ xhat requires mu_minus=mu/q.
    w=1.37;q=-np.exp(-1j*np.pi*w);X=1.23*np.exp(-1j*np.pi*w/2)
    mu=-np.exp(-1j*np.pi*w/2)/(q-1/q);Q=np.array([.2,-.05,-.15]);ww=weights(3,True)
    M,sv=solve_sector(Q,q,1.23,mu,bar=True,orient=-1)
    M*=q**(Q[0]-ww[0,0])/M[0,0]
    expected=np.zeros_like(M)
    for i in range(3):
        N=Q-ww[i];j=(i+1)%3
        expected[i,i]=q**N[i];expected[j,i]=X*q**(N[i]-N[j])
    print('orientation-minus closed form residual',np.linalg.norm(M-expected))
    # Congruence of the shifted defect rapidity does not select a quantum width.
    h=6;p=1;a=1;td=.7;w=1.37
    for k in [0,1]:
        shift=2j*np.pi*p/h+2j*np.pi*k/(h*w)
        pole=td+shift-1j*(np.pi/2-np.pi*a/h)-1j*np.pi/(h*w)
        print('classical-equivalent representatives k,pole,E,matrix-cycle',k,pole,np.cosh(pole),np.exp(-h*w*shift))
    # At an exact sector boundary the Hirota solution is singular. A distance
    # epsilon inside the boundary gives the following finite threshold.
    for eps in [.03,.2,np.pi/3]:
        gamma=np.pi/3
        threshold=-.5*np.log(np.sin(2*gamma)/np.tan(eps)-np.cos(2*gamma))
        print('sector boundary distance, theta_star - theta_d',eps,threshold)
    # CZ2011 Gamma arguments fix xi/2, while the printed pole uses xi.
    xi=.8
    for m in [1,2]:
        print('CZ2011 m,corrected Gamma argument,printed Gamma argument',m,
              .5-(xi/2+m*np.pi/3)/xi+m*np.pi/(3*xi),
              .5-(xi+m*np.pi/3)/xi+m*np.pi/(3*xi))

if __name__=='__main__':run()
