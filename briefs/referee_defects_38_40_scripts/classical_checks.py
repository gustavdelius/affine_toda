"""Independent roots, classification, r-matrix, delays and reality checks."""
import itertools
from collections import Counter
import numpy as np
import sympy as s

def basis(h):
    U=np.zeros((h,h-1))
    for j in range(1,h):
        U[:j,j-1]=1/np.sqrt(j*(j+1)); U[j,j-1]=-j/np.sqrt(j*(j+1))
    return U

def an(r):
    h=r+1; U=basis(h); I=np.eye(h)
    roots=np.array([I[j]-I[(j+1)%h] for j in range(h)])@U
    P=U.T@np.roll(I,-1,axis=0)@U
    D=2*np.linalg.inv(np.eye(r)-P)
    return U,roots,D,np.eye(r)-D

def diagrams():
    data={}
    for r in range(1,5):data[f'a{r}(1)']=an(r)[1]
    for r in (3,4):
        I=np.eye(r)
        data[f'b{r}(1)']=np.array([-I[0]-I[1]]+[I[j]-I[j+1] for j in range(r-1)]+[I[-1]])
    for r in (2,3,4):
        I=np.eye(r)
        data[f'c{r}(1)']=np.array([-2*I[0]]+[I[j]-I[j+1] for j in range(r-1)]+[2*I[-1]])/np.sqrt(2)
    I=np.eye(4)
    data['d4(1)']=np.array([-I[0]-I[1],I[0]-I[1],I[1]-I[2],I[2]-I[3],I[2]+I[3]])
    data['f4(1)']=np.array([-I[0]-I[1],I[1]-I[2],I[2]-I[3],I[3],(I[0]-I[1]-I[2]-I[3])/2])
    short=np.array([np.sqrt(2/3),0]); long=np.array([-np.sqrt(1.5),1/np.sqrt(2)])
    data['g2(1)']=np.array([-3*short-2*long,short,long])
    for src,dst in [('b3(1)','a5(2)'),('b4(1)','a7(2)'),('c2(1)','d3(2)'),('c3(1)','d4(2)'),('c4(1)','d5(2)'),('f4(1)','e6(2)'),('g2(1)','d4(3)')]:
        A=data[src]; dual=2*A/np.sum(A*A,axis=1)[:,None]
        data[dst]=dual*np.sqrt(2/np.max(np.sum(dual*dual,axis=1)))
    for r in range(1,5):
        I=np.eye(r)
        data[f'a{2*r}(2)']=np.array([-2*I[0]]+[I[j]-I[j+1] for j in range(r-1)]+[I[-1]])/np.sqrt(2)
    return data

def scan():
    outcomes=Counter(); passed=[]
    for name,A in diagrams().items():
        G=A@A.T; N=len(A); r=A.shape[1]
        for p in itertools.permutations(range(N)):
            if not np.allclose(G[np.ix_(p,p)],G,atol=1e-10):continue
            if any(j==p[j] for j in range(N)):outcomes['fixed node']+=1; continue
            P=(np.linalg.pinv(A)@A[list(p)]).T
            if np.linalg.matrix_rank(np.eye(r)-P,tol=1e-8)<r:outcomes['singular']+=1; continue
            D=2*np.linalg.inv(np.eye(r)-P); M=A@D@A.T
            forbidden=[(i,j) for i in range(N) for j in range(N) if j not in (i,p[i]) and abs(M[i,j])>1e-8]
            if forbidden:
                outcomes['cross terms']+=1
                print(name,p,'forbidden count',len(forbidden),'A D A^T=',np.round(M,5).tolist())
            else:outcomes['pass']+=1; passed.append((name,p))
    print('scan:',len(diagrams()),'diagrams,',sum(outcomes.values()),'automorphisms,',dict(outcomes),'passes=',passed,flush=True)
    assert outcomes==Counter({'fixed node':64,'singular':5,'cross terms':2,'pass':7})
    # The abstract proof: transitivity forces equal lengths. A finite tree
    # with >=3 vertices is not vertex transitive. On a cycle the equation
    # A D A^T has nearest-neighbour support only for either basic rotation.
    for r in range(2,10):
        _,A,D,E=an(r); v=np.arange(1,r+1)/7
        h=r+1; U=basis(h)
        cubic=lambda u:np.sum((U@u)**3)
        assert abs(cubic(D@v)+cubic(D.T@v))<1e-9
    print('cubic cancellation independently tested ranks 2..9',flush=True)

def rmat(U,z,b):
    h,r=U.shape; out=np.zeros((h*h,h*h),complex)
    for j in range(r):out+=(z**h+1)*np.kron(np.diag(U[:,j]),np.diag(U[:,j]))
    for i in range(h):
        for j in range(i+1,h):
            e=np.zeros((h,h));e[i,j]=1
            out+=2*(z**(j-i)*np.kron(e,e.T)+z**(h-j+i)*np.kron(e.T,e))
    return b*b*out/(4*(z**h-1))

def bracket_check():
    rng=np.random.default_rng(4038)
    for r in range(1,5):
        U,A,D,E=an(r); h=r+1
        ph=rng.normal(size=r)/4; ps=rng.normal(size=r)/4
        worst=0
        for typ,b in itertools.product((1,-1),(.73,.73j)):
            Et=typ*E; Dt=typ*D; K=np.block([[Et,Dt],[-Dt.T,Et]]); J=np.linalg.inv(K)
            left=-b/2*(U@np.c_[Et,Dt]);right=b/2*(U@np.c_[Dt.T,-Et])
            v=np.r_[ph,ps]; S=np.roll(np.eye(h),1,axis=1)
            def L(lam):
                M=np.eye(h)+.71/lam*S.T if typ==1 else sum((lam/.71)**k*np.linalg.matrix_power(S,k) for k in range(h))
                return np.exp(left@v)[:,None]*M*np.exp(right@v)[None,:]
            lam1=.83+.13j;lam2=1.49-.17j;L1=L(lam1);L2=L(lam2)
            pb=np.zeros((h*h,h*h),complex)
            for i,j,k,l in itertools.product(range(h),repeat=4):
                pb[i*h+k,j*h+l]=L1[i,j]*L2[k,l]*(left[i]+right[j])@J@(left[k]+right[l])
            R=rmat(U,lam1/lam2,b); LL=np.kron(L1,L2); rhs=R@LL-LL@R
            err=np.linalg.norm(pb-rhs)/(1+np.linalg.norm(rhs));worst=max(worst,err)
        print('rank',r,'defect Sklyanin residual',worst,flush=True)
        assert worst<1e-11
        # Verify the local bulk bracket and hence the sign for transport.
        b=.73;mom=rng.normal(size=r)
        def ax(lam):
            X=b/2*np.diag(U@mom); derivatives=[]
            for a in A:derivatives.append(None)
            for j,a in enumerate(A):
                e=np.zeros((h,h));e[j,(j+1)%h]=1
                X=X+.5*np.exp(b*a@ph/2)*(lam*e+e.T/lam)
            ds=[]
            for k in range(r):
                m=np.zeros((h,h),complex)
                for j,a in enumerate(A):
                    e=np.zeros((h,h));e[j,(j+1)%h]=1
                    m+=b/4*a[k]*np.exp(b*a@ph/2)*(lam*e+e.T/lam)
                ds.append(m)
            return X,ds,[b/2*np.diag(U[:,k]) for k in range(r)]
        X,dx,dp=ax(lam1);Y,dy,dq=ax(lam2)
        pb=sum(np.kron(dx[k],dq[k])-np.kron(dp[k],dy[k]) for k in range(r))
        R=rmat(U,lam1/lam2,b);XY=np.kron(X,np.eye(h))+np.kron(np.eye(h),Y)
        err=np.linalg.norm(pb-(R@XY-XY@R))
        print('rank',r,'bulk + sign residual',err,'; transport generator -a_x has - sign',flush=True)
        assert err<1e-11

def delay_reality():
    maxerr=0.; maxidentity=0.; wrong=0.
    for r in range(1,5):
        U,A,D,E=an(r);h=r+1
        for a in range(1,h):
            theta=.83;sig=.67;td=-np.log(sig)
            za=lambda aa:(np.exp(-theta)-1j*sig*np.exp(1j*np.pi*aa/h))/(np.exp(-theta)-1j*sig*np.exp(-1j*np.pi*aa/h))
            z=za(a);eps=np.exp(1j*np.pi*(2*a-h)/(2*h))
            maxidentity=max(maxidentity,abs(z-(-np.exp(-theta)-eps*sig)/(-np.exp(-theta)+sig/eps)))
            maxidentity=max(maxidentity,abs(z*za(h-a)-(np.cosh(theta-td)+np.sin(np.pi*a/h))/(np.cosh(theta-td)-np.sin(np.pi*a/h))))
            for t in (-.8,.2,1.1):
                wave=np.exp(2j*np.pi*a*np.arange(h)/h)*np.exp(-2*np.sin(np.pi*a/h)*np.sinh(theta)*t+.17+.39j)
                def sol(z):
                    w=z*wave; f=-A.T@np.log(1+w)
                    p=-A.T@(w/(1+w)); return f,-2*np.sin(np.pi*a/h)*np.sinh(theta)*p,2*np.sin(np.pi*a/h)*np.cosh(theta)*p
                ph,pt,px=sol(1);ps,st,sx=sol(z)
                def sew(Db):
                    Eb=np.eye(r)-Db
                    f=np.exp(A@(Db.T@ph+Db@ps)/2);g=np.exp(A@Db@(ph-ps)/2)
                    gp=(Db@A.T@(sig*f)+Db.T@A.T@(g/sig))/2
                    gs=(Db.T@A.T@(sig*f)-Db.T@A.T@(g/sig))/2
                    return np.r_[Eb@pt+Db@st-gp,Db.T@pt-Eb@st+gs]
                maxerr=max(maxerr,np.max(np.abs(sew(D)-np.r_[px,sx])))
                if r==2:wrong=max(wrong,np.max(np.abs(sew(np.eye(r)-.6*E)-np.r_[px,sx])))
        # Imaginary-coupling potential and conjugation, including complex sigma.
        ph=np.arange(r)*.1+.15j;ps=np.arange(r)*.13-.21j;b=.8j;sig=.7+.23j
        def potential(ph,ps,sig):
            return sum(sig*(np.exp(b*A@(D.T@ph+D@ps)/2)-1)+(np.exp(b*A@D@(ph-ps)/2)-1)/sig)/(b*b)
        err=abs(np.conj(potential(ph,ps,sig))-potential(-ph.conj(),-ps.conj(),sig.conjugate()))
        assert err<1e-12
    print('a_r soliton sewing ranks 1..4=',maxerr,'delay identities=',maxidentity,'nonintegrable-D soliton sewing control=',wrong,flush=True)
    assert maxerr<1e-10 and maxidentity<1e-10 and wrong>.01
    print('Theta potential identity ranks 1..4: <1e-12 (complex sigma)',flush=True)

if __name__=='__main__':
    scan();bracket_check();delay_reality()
