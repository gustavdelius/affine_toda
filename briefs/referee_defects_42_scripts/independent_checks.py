"""Independent chapter-42 checks; no textbook_code imports. Run with BLAS threads=1."""
import numpy as np
import scipy.linalg as la
import sympy as sp
from scipy.sparse import csc_matrix, hstack, vstack
from scipy.sparse.linalg import lsqr

LAM = 1.37
q = -np.exp(-1j*np.pi*LAM)

def transmission(y, d, N, sign=-1, s=0.0):
    ns = np.arange(-N, N+1); w = len(ns)
    A = np.diag(q**(-ns if sign == -1 else ns))
    D = np.diag(q**(ns+s if sign == -1 else -ns-s))
    B = np.zeros((w,w), complex); C = B.copy()
    for j,n in enumerate(ns):
        if j:
            B[j-1,j] = -1j*y/d*q**(n+s) if sign == -1 else 1j*d/y*q**(n-1)
        if j+1<w:
            C[j+1,j] = -1j*y/d*q**(-n) if sign == -1 else 1j*d/y*q**(-n-1-s)
    return np.block([[A,B],[C,D]])

def rep(y,j=1):
    ep = np.array([[0,1],[0,0]],complex)
    e = y*(ep if j==1 else ep.T)
    f = (ep.T if j==1 else ep)/y
    k = np.diag([q,q**-1] if j==1 else [q**-1,q])
    return e,f,k

def doubled(y,d,N,sign=-1):
    L = transmission(y,d,N,sign); Lm = transmission(1/y,d,N,sign)
    Z = np.zeros_like(L)
    return np.block([[Z,la.inv(Lm)],[L,Z]])

def ls_test(N,sign,kind,tw,td=.15,full=False,margin=2):
    w=2*N+1; ns=np.arange(-N,N+1); I=np.eye(w); z=np.zeros((w,w))
    # All matrix entries are free when full=True. Otherwise use precisely the forced charge grade.
    shifts = (1,0 if tw==0 else 2) if kind=='e' else ((-1 if tw==0 else 1),0 if tw==0 else 2)
    entries=[]
    for name,shift in enumerate(shifts):
        for c in range(w):
            for r in (range(w) if full else [c+shift]):
                if 0<=r<w: entries.append((name,r,c))
    masks=np.tile(abs(ns)<=N-margin,4)
    def X(y,Y,kap,constant):
        e,f,k=rep(y); et,ft,kt=rep(1/y,1-tw)
        if kind=='e':
            a=np.kron(e,I)*constant+np.kron(k,Y)
            b=np.kron(np.eye(2),Y)+np.kron(et,kap)
        else:
            a=np.kron(f,kap)+np.kron(np.eye(2),Y)
            b=np.kron(la.inv(kt),Y)+np.kron(ft,I)*constant
        return la.block_diag(a,b)
    mats=[]; targets=[]
    for th in (-.35,.21,.64):
        y=np.exp(LAM*th); K=doubled(y,np.exp(LAM*td),N,sign)
        def res(Y,kap,c):
            M=K@X(y,Y,kap,c)-X(1/y,Y,kap,c)@K
            return M[np.ix_(masks,masks)].ravel()
        targets.append(-res(z,z,1))
        cols=[]
        for name,r,c in entries:
            E=z.copy(); E[r,c]=1
            col=res(E if name==0 else z,E if name==1 else z,0)
            col[abs(col)<1e-14]=0
            cols.append(csc_matrix(col[:,None]))
        mats.append(hstack(cols,format='csc'))
    A=vstack(mats,format='csc'); b=np.concatenate(targets)
    sol=lsqr(A,b,atol=1e-12,btol=1e-12,iter_lim=6000)
    return la.norm(A@sol[0]-b)/la.norm(b)

def borels():
    N=8; w=2*N+1; ns=np.arange(-N,N+1); I=np.eye(w)
    inner=np.tile(abs(ns)<=N-4,2); sl=np.ix_(inner,inner)
    shift=lambda h,fn: np.array([[fn(n) if m==n+h else 0 for n in ns] for m in ns],complex)
    d=np.exp(LAM*.23); s=.3; mu=1j/(q-q**-1)
    k=np.diag(q**(2*ns+s)); f1=mu/d*shift(-1,lambda n:1); f0=f1.T
    a=1j*d*q**(-s)/(q*q-1); b=1j*d*q**s/(q*q-1)
    t1=shift(2,lambda n:q**(-2*n-2-s)); t0=shift(-2,lambda n:q**(2*n-2+s))
    e1=shift(1,lambda n:a*q**(-2*n)); e0=shift(-1,lambda n:b*q**(2*n))
    mid=np.ix_(abs(ns)<=N-4,abs(ns)<=N-4)
    rels=[t1@t0-I,t1@e1-q**2*e1@t1,t1@e0-q**-2*e0@t1,
          e1@e0-a*b*q*q*I,k@t1-q**4*t1@k]
    serre=e1@e1@e1@e0-(q*q+1+q**-2)*e1@e1@e0@e1+(q*q+1+q**-2)*e1@e0@e1@e1-e0@e1@e1@e1
    print('second-module relations/Serre:',max(abs(R[mid]).max() for R in rels+[serre]))
    errs=[]
    for y in (np.exp(.41*LAM),np.exp((-.3+.17j)*LAM)):
        L=transmission(y,d,N,-1,s)
        for j,F,K,Y,kap in ((1,f1,k,e1,t1),(0,f0,la.inv(k),e0,t0)):
            ep,fp,kp=rep(y,j); et,_,kt=rep(y,1-j)
            r=L@(np.kron(fp,la.inv(K))+np.kron(np.eye(2),F))-(np.kron(la.inv(kp),F)+np.kron(fp,I))@L
            errs.append(abs(r[sl]).max())
            r=L@(np.kron(ep,I)+np.kron(kp,Y))-(np.kron(np.eye(2),Y)+np.kron(et,kap))@L
            errs.append(abs(r[sl]).max())
            r=L@np.kron(kp,kap)-np.kron(kt,kap)@L
            errs.append(abs(r[sl]).max())
    print('explicit diagonal b-minus and twisted b-plus intertwining:',max(errs))

def fusion():
    y,d,r=sp.symbols('y d r',nonzero=True)
    H=lambda z:sp.Matrix([[1,-sp.I*z],[-sp.I*z,1]])
    K=sp.diag(y*r+1/y/r,y/r+r/y)
    B=sp.simplify(H(1/y/d).inv()*K*H(y/d))
    expected=sp.Matrix([[K[0,0]+K[1,1]/d**2,-sp.I*r/d*(y**2-y**-2)],[-sp.I/r/d*(y**2-y**-2),K[1,1]+K[0,0]/d**2]])/(1+1/y**2/d**2)
    print('fusion symbolic block:',sp.simplify(B-expected)==sp.zeros(2))
    eta=.37; td=.21; d0=np.exp(LAM*td)
    def gz(th,v):
        y=np.exp(LAM*th); c1=np.cos(eta-1j*v);c0=np.cos(eta+1j*v)
        return np.array([[-2j*(c1*y+c0/y),y*y-y**-2],[y*y-y**-2,-2j*(c0*y+c1/y)]])
    inv=lambda M:np.array([M[0,0]/M[1,1],M[0,1]*M[1,0]/(M[0,0]*M[1,1])])
    errs=[];bad1=[];bad2=[]
    for N in (5,7,9):
        w=2*N+1
        for th in (-.35,.21+.1j,.64):
            y0=np.exp(LAM*th); R=np.exp(1j*eta)
            D=np.diag([y0*R+1/y0/R,y0/R+R/y0])
            F=la.solve(transmission(1/y0,d0,N),np.kron(D,np.eye(w))@transmission(y0,d0,N))
            for n in range(-N+2,N-2):
                ix=[n+N,w+n+N+1]; f=F[np.ix_(ix,ix)]
                errs.append(abs(inv(f)-inv(gz(th,LAM*td))).max())
                bad1.append(abs(inv(f)-inv(gz(th,-LAM*td))).max())
                bad2.append(abs(inv(f)-inv(gz(th,td))).max())
    print('GZ invariants: correct/max, negative-sign/min, no-lambda/min:',max(errs),min(bad1),min(bad2))
    for th,td in ((.41,.32),(-.27,.61),(.21+.15j,-.38)):
        T=lambda th:-1j*np.sinh((th-td)/2-1j*np.pi/4)/np.sinh((th-td)/2+1j*np.pi/4)
        F=2j*td/np.pi
        block=lambda x:np.sinh(th/2+1j*np.pi*x/4)/np.sinh(th/2-1j*np.pi*x/4)
        val=T(th)/T(-th)
        print('sinh-Gordon fusion:',abs(val-1/(block(1-F)*block(1+F))))

def classical():
    rng=np.random.default_rng(424242); beta=.73; sigma=.81; lam=.94+.28j
    for rank,C in ((1,[.61,-1.17]),(1,[0,0]),(2,[1,-1,1]),(2,[0,0,0])):
        h=rank+1; Id=np.eye(h)
        roots=[Id[-1]-Id[0]]+[Id[i-1]-Id[i] for i in range(1,h)]
        weights=[np.zeros(h)]+[sum(Id[:i])-i/h*np.ones(h) for i in range(1,h)]+[np.zeros(h)]
        D=2*sum(np.outer(weights[i]-weights[i+1],weights[i]) for i in range(1,h))
        # All fields lie in the sum-zero root subspace; extend D as identity on its orthogonal complement.
        D=D+np.ones((h,h))/h; E=Id-D
        Es=[np.outer(Id[-1],Id[0])]+[np.outer(Id[i-1],Id[i]) for i in range(1,h)]
        def L(z,ph,ps):
            A=np.diag(np.exp(-beta/2*(E@ph+D@ps)))
            B=np.diag(np.exp(beta/2*(D.T@ph-E@ps)))
            return A@(Id+sigma/z*sum(e.T for e in Es))@B
        def at(z,ph,px):
            return beta/2*np.diag(px)+sum(np.exp(beta*a@ph/2)*(z*e-e.T/z)/2 for a,e in zip(roots,Es))
        # Solve the field-independent boundary equation at each independent exponential coefficient.
        eq=[]
        for a,e,c in zip(roots,Es,C):
            A=-c*np.diag(a)+lam*e-e.T/lam
            B=-c*np.diag(a)+e.T/lam-lam*e
            eq.append(np.kron(A.T,Id)+np.kron(Id,B))
        basis=la.null_space(np.vstack(eq)); assert basis.shape[1]==1
        K=basis[:,0].reshape((h,h),order='F')
        ph,ps,psit=rng.normal(size=(3,h))*.3
        ph-=ph.mean();ps-=ps.mean();psit-=psit.mean()
        gp=np.zeros(h);gs=np.zeros(h)
        for a in roots:
            f=np.exp(beta/2*a@(D.T@ph+D@ps));g=np.exp(beta/2*a@D@(ph-ps))
            gp+=(sigma*f*(D@a)+g/sigma*(D.T@a))/(2*beta)
            gs+=(sigma*f*(D.T@a)-g/sigma*(D.T@a))/(2*beta)
        gb=sum(c*a*np.exp(beta/2*a@ps)/beta for a,c in zip(roots,C))
        phit=la.solve(D.T,-gb+E@psit-gs); phix=E@phit+D@psit-gp
        def dressed(ph,ps,mode):
            A=L(1/lam,ph,ps);B=L(lam,ph,ps)
            if mode=='bare':return K
            if mode=='inverse':A=la.inv(A);B=la.inv(B)
            return A.T@K@B
        vals=[]
        for mode in ('inverse','direct','bare'):
            for step in ((1e-4,1e-5,1e-6) if mode=='inverse' else (1e-5,)):
                dot=(dressed(ph+step*phit,ps+step*psit,mode)-dressed(ph-step*phit,ps-step*psit,mode))/(2*step)
                K0=dressed(ph,ps,mode);rhs=K0@at(lam,ph,phix)+at(1/lam,ph,phix).T@K0
                vals.append(la.norm(dot-rhs)/max(la.norm(rhs),1e-12))
        print('classical rank/C/inv(h=1e-4,1e-5,1e-6),direct,bare:',rank,C,['%.3g'%v for v in vals])
    phi,psi,b,t=sp.symbols('phi psi b t',real=True)
    potential=sp.exp(-t)*sp.cos(b*(phi+psi)/2)+sp.exp(t)*sp.cos(b*(phi-psi)/2)
    osc=sp.exp(-sp.I*b*phi/2)*sp.cos(b*psi/2-sp.I*t)+sp.exp(sp.I*b*phi/2)*sp.cos(b*psi/2+sp.I*t)
    print('BDG potential exact:',sp.simplify((potential-osc).rewrite(sp.exp))==0)
    flow=-.5*(np.exp(1j*.31/2)*np.sin(.72/2+1j*.43)+np.exp(-1j*.31/2)*np.sin(.72/2-1j*.43))
    print('BDG flow of p-q at beta=m0=1:',flow)

def crossing():
    a,b,c,d,theta=1.3,2.1,1.7,.9,.37+.11j
    F=lambda outer,z:(a*c/outer**2)**(1j*z/(2*np.pi))
    left=F(b,1j*np.pi/2-theta)*np.sqrt(d/b)
    right=(b*d/(a*c))**(1j*theta/np.pi)*F(d,1j*np.pi/2+theta)
    wrong=(a*c/(b*d))**(1j*theta/np.pi)*F(d,1j*np.pi/2+theta)
    print('boundary crossing gauge identity/control:',abs(left-right),abs(left-wrong))

if __name__=='__main__':
    for sign in (-1,1):
        for kind,tw in (('e',0),('e',1),('f',0),('f',1)):
            vals=[ls_test(N,sign,kind,tw) for N in (5,7,9)]
            print('doubled LS sign/kind/twist N=5,7,9:',sign,kind,tw,['%.8g'%v for v in vals],flush=True)
    for sign,kind,tw in ((-1,'e',0),(-1,'f',1),(1,'f',0),(1,'e',1)):
        print('unrestricted ALL entries N=3:',sign,kind,tw,ls_test(3,sign,kind,tw,full=True,margin=1),flush=True)
    for sign,kind,tw,direction in ((-1,'e',0,1),(-1,'f',1,-1),(1,'e',1,-1),(1,'f',0,1)):
        print('limits sign/kind/twist/direction:',sign,kind,tw,direction,[ls_test(5,sign,kind,tw,td=direction*t) for t in (2,4,6)],flush=True)
    borels();fusion();classical();crossing()
