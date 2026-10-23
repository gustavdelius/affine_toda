"""Independent charge-sector computations; no imports from code.

All matrices use outgoing rows, incoming columns. Indices here run 0..h-1;
node a joins a to a+1, so the affine node is h-1.
"""
import itertools
from functools import lru_cache
import numpy as np
from scipy.linalg import null_space


def weights(h, bar=False):
    return (-1 if bar else 1)*(np.eye(h)-np.ones((h,h))/h)


def roots(h):
    return np.eye(h)-np.roll(np.eye(h), -1, axis=0)


def rep(h,q,y,bar=False):
    ans={}
    for a in range(h):
        b=(a+1)%h
        E=np.zeros((h,h),complex); E[a,b]=1
        H=np.eye(h)[a]-np.eye(h)[b]
        ans['k',a]=np.diag(q**((-1 if bar else 1)*H))
        ans['e',a]=y*(E.T if bar else E)
        ans['f',a]=(E if bar else E.T)/y
    return ans


def tensor(A,B):
    da=A['k',0].shape[0]; db=B['k',0].shape[0]
    C={}
    for g,j in A:
        if g=='k': C[g,j]=np.kron(A[g,j],B[g,j])
        if g=='e': C[g,j]=np.kron(A[g,j],np.eye(db))+np.kron(A['k',j],B[g,j])
        if g=='f': C[g,j]=np.kron(A[g,j],np.linalg.inv(B['k',j]))+np.kron(np.eye(da),B[g,j])
    return C


def finite_hom(A,B,wa,wb):
    """Solve B T = T A, retaining all equal-weight matrix entries."""
    na=len(wa); nb=len(wb)
    pairs=[(b,a) for b in range(nb) for a in range(na) if np.max(abs(wb[b]-wa[a]))<1e-8]
    basis=np.zeros((len(pairs),nb,na),complex)
    for n,(b,a) in enumerate(pairs): basis[n,b,a]=1
    equations=[]
    for key in A:
        eq=np.einsum('ij,njk->nik',B[key],basis)-np.einsum('nij,jk->nik',basis,A[key])
        equations.append(eq.reshape(len(pairs),-1).T)
    eq=np.concatenate(equations)
    u,s,vh=np.linalg.svd(eq,full_matrices=False)
    z=vh[-1].conj(); T=np.einsum('n,nij->ij',z,basis)
    return T,s


@lru_cache(None)
def rcheck(h,q,y1,y2,bar1=False,bar2=False):
    A=rep(h,q,y1,bar1); B=rep(h,q,y2,bar2)
    w1=weights(h,bar1); w2=weights(h,bar2)
    T,s=finite_hom(tensor(A,B),tensor(B,A),
                   np.array([a+b for a in w1 for b in w2]),
                   np.array([b+a for b in w2 for a in w1]))
    assert s[-1]<1e-10 and s[-2]>1e-5, s[-2:]
    # Overall normalization cancels from all the checks below.
    T/=np.max(abs(T))
    return T


def module_action(N,a,q,mu,plus=False,orient=1):
    h=len(N); b=(a+1)%h
    return N+(1 if plus else -1)*roots(h)[a], mu*q**(orient*N[b])


def sector_generator(Q,a,q,y,mu,bar=False,plus=False,orient=1,target=False):
    """Delta generator between Q and Q +/- alpha; both blocks indexed by V."""
    h=len(Q); w=weights(h,bar); v=rep(h,q,y,bar)
    g='e' if plus else 'f'
    A=np.zeros((h,h),complex)
    for i in range(h):
        N=Q-w[i]
        k=q**np.dot(roots(h)[a],N)
        _,t=module_action(N,a,q,mu,plus,orient)
        for j in range(h):
            if plus:
                A[j,i]+=v[g,a][j,i]*(k if target else 1)
            else:
                A[j,i]+=v[g,a][j,i]*(1 if target else 1/k)
        A[i,i]+=t*((1 if plus else 1/v['k',a][i,i]) if target
                    else (v['k',a][i,i] if plus else 1))
    return A


def sector_T(Q,q,X,bar=False,plus=False):
    h=len(Q); w=weights(h,bar); M=np.zeros((h,h),complex)
    for i in range(h):
        N=Q-w[i]
        if plus:
            assert bar
            M[i,i]=q**(-N[i]); j=(i-1)%h
            M[j,i]=X*q**(N[i]-N[j]-1)
        elif not bar:
            M[i,i]=q**(-N[i]); M[(i-1)%h,i]=X
        else:
            for j in range(h):
                k=(i+j)%h
                M[k,i]=X**j*q**sum(N[(i+t)%h] for t in range(j+1))
    return M


def solve_sector(Q,q,y,mu,bar=False,plus=False,orient=1):
    """No defect cutoff: solve intertwining of closed transport loops at Q."""
    h=len(Q); al=roots(h); sign=1 if plus else -1
    loops=[[],[]]
    for target in [False,True]:
        def A(a,P): return sector_generator(P,a,q,y,mu,bar,plus,orient,target)
        for a,b in itertools.combinations(range(h),2):
            ab=A(a,Q+sign*al[b])@A(b,Q)
            ba=A(b,Q+sign*al[a])@A(a,Q)
            loops[target].append(np.linalg.solve(ba,ab))
    basis=np.eye(h*h).reshape(h*h,h,h)
    rows=[]
    for A,B in zip(*loops):
        eq=np.einsum('ij,njk->nik',B,basis)-np.einsum('nij,jk->nik',basis,A)
        rows.append(eq.reshape(h*h,-1).T)
    _,sv,vh=np.linalg.svd(np.concatenate(rows),full_matrices=False)
    return vh[-1].conj().reshape(h,h),sv


def transitions(i,N,q,X,bar=False,plus=False):
    h=len(N); w=weights(h,bar); M=sector_T(N+w[i],q,X,bar,plus)
    return [(j,N+w[i]-w[j],M[j,i]) for j in range(h) if abs(M[j,i])>1e-13]


def triangle(h,q,y1,y2,mu,bar1=False,bar2=False,plus=False):
    R=rcheck(h,q,y1,y2,bar1,bar2).reshape(h,h,h,h)
    X1=(q-1/q)*mu/y1 if plus else -(q-1/q)*mu*y1
    X2=(q-1/q)*mu/y2 if plus else -(q-1/q)*mu*y2
    err=0
    for N in [np.zeros(h), np.arange(h)-(h-1)/2]:
      for i,j in itertools.product(range(h),repeat=2):
        lhs={}; rhs={}
        def add(d,k,c): d[k]=d.get(k,0)+c
        for b,N1,t2 in transitions(j,N,q,X2,bar2,plus):
          for a,N2,t1 in transitions(i,N1,q,X1,bar1,plus):
            for c,d in itertools.product(range(h),repeat=2):
              add(lhs,(c,d,*np.round(N2,8)),R[c,d,a,b]*t1*t2)
        for b,a in itertools.product(range(h),repeat=2):
          for d,N1,t1 in transitions(a,N,q,X1,bar1,plus):
            for c,N2,t2 in transitions(b,N1,q,X2,bar2,plus):
              add(rhs,(c,d,*np.round(N2,8)),t2*t1*R[b,a,i,j])
        scale=max(max(map(abs,lhs.values())),1)
        err=max(err,max(abs(lhs.get(k,0)-rhs.get(k,0)) for k in lhs.keys()|rhs.keys())/scale)
    return err


def crossing(Q,q,X,physical_shift,bar=False,C=None):
    h=len(Q); w=weights(h,bar); wo=-w
    if C is None: C=np.eye(h)
    # Here C is diagonal. A separate general contraction is used for a22.
    M=sector_T(Q,q,X,bar); cross=np.zeros_like(M)
    for j,i in itertools.product(range(h),repeat=2):
        Nin=Q-w[i]; Nout=Q-w[j]
        charge=Nin+wo[j]
        inv=np.linalg.inv(sector_T(charge,q,X*physical_shift,not bar))
        cross[j,i]=C[i,i]/C[j,j]*inv[i,j]
    factor=np.vdot(cross,M)/np.vdot(cross,cross)
    return factor,np.linalg.norm(M-factor*cross)/np.linalg.norm(M)


def serde(h,q,orient=1,plus=False):
    N=np.arange(h)-(h-1)/2; mu=0.71
    def word(seq):
        P=N.copy(); out=1
        for a in reversed(seq):
            P,c=module_action(P,a,q,mu,plus,orient); out*=c
        return out
    errs=[]
    for a in range(h):
      for b in range(h):
        if a==b: continue
        adj=(a-b)%h in [1,h-1]
        if adj:
            errs.append(abs(word([a,a,b])-(q+1/q)*word([a,b,a])+word([b,a,a])))
        else: errs.append(abs(word([a,b])-word([b,a])))
    return max(errs)


def run():
    omega=1.37; q=-np.exp(-1j*np.pi*omega); mu=-np.exp(-1j*np.pi*omega/2)/(q-1/q)
    y=1.19; X=-(q-1/q)*mu*y
    for h in [3,4,5]:
        print('Serre h',h,'both orientations/Borels',max(serde(h,q,o,p) for o in [-1,1] for p in [False,True]),'plain control',serde(h,q,0))
        Q=np.arange(h)*0.173; Q-=Q.mean()
        for bar,plus in [(False,False),(True,False),(True,True)]:
            xp=(q-1/q)*mu/y if plus else X
            maxerr=0
            for a in range(h):
                A=sector_generator(Q,a,q,y,mu,bar,plus)
                B=sector_generator(Q,a,q,y,mu,bar,plus,target=True)
                M=sector_T(Q,q,xp,bar,plus)
                M2=sector_T(Q+(1 if plus else -1)*roots(h)[a],q,xp,bar,plus)
                maxerr=max(maxerr,np.linalg.norm(M2@A-B@M))
            print('intertwines h,bar,plus',h,bar,plus,maxerr)
            if h<=4:
                T,sv=solve_sector(Q,q,y,mu,bar,plus)
                c=np.vdot(T,M)/np.vdot(T,T)
                print('loop solve nullity/gap/closed-form residual',sum(sv<1e-10),sv[-2],np.linalg.norm(c*T-M))
        if h<=4:
          for bar1,bar2 in [(False,False),(True,True),(False,True)]:
            print('triangle',h,bar1,bar2,triangle(h,q,1.43,.87,mu,bar1,bar2))
          print('triangle plus barbar',h,triangle(h,q,1.43,.87,mu,True,True,True))
          for bar in [False,True]:
            T,sv=solve_sector(Q,q,y,mu,bar,orient=-1)
            print('orientation -1 h,bar,nonzeros per column',h,bar,(abs(T)>1e-9).sum(axis=0).tolist(), 'nullity',sum(sv<1e-10))
        for sign in [1,-1]:
          for cflag in [0,1]:
            C=np.diag((-1.)**np.arange(h)) if cflag else np.eye(h)
            es=[]; fs=[]
            for bar in [False,True]:
              for s in [np.zeros(h),weights(h)[0],Q]:
                c,e=crossing(weights(h,bar)[0]+s,sign*q,X,-q,bar,C)
                es.append(e)
                expected=1-q*X**h if bar else 1+(-q)**(h-1)*X**h
                if sign==1 and cflag==0: fs.append(abs(c-expected))
            print('crossing h,qsign,Cprime,maxerr,scalarerr',h,sign,cflag,max(es),max(fs,default=0))
        M=sector_T(Q,q,X)
        H=M.conj().T@M; H-=np.diag(np.diag(H))
        print('nonunitary offdiag',h,np.max(abs(H)))
        pole=1j*np.pi/2-1j*(np.pi*(h-1)/h+np.pi/(h*omega))
        XP=np.exp(omega*(pole-1j*np.pi/2))
        MP=sector_T(Q,q,XP)
        det=q**(h-1)*(1+(-1/q)**(h-1)*X**h)
        print('det/rank',h,abs(np.linalg.det(M)-det),np.linalg.svd(MP,compute_uv=False).tolist())


if __name__=='__main__': run()
