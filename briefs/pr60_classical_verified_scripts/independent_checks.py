"""PR60 independent classical checks. No repository/referee imports."""
import itertools
import numpy as np
import sympy as s
from scipy.integrate import quad

def check(name, ok, detail=''):
    print(('PASS ' if ok else 'FAIL ')+name+(' '+str(detail) if detail else ''), flush=True)
    if not ok: raise AssertionError(name)

def zmat(M): return all(s.simplify(s.expand(x))==0 for x in M)

def rank_one():
    u,v,l,k,c=s.symbols('u v l k c', nonzero=True)
    # u=exp(b(phi-psi)/(2sqrt2)), v=exp(-b(phi+psi)/(2sqrt2)).
    def L(t): return s.Matrix([[u,v/t],[1/(v*t),1/u]])
    H=s.diag(1,-1)/s.sqrt(2)
    e=s.Matrix([[0,1],[0,0]]); f=e.T
    Pi=s.kronecker_product(e,f)+s.kronecker_product(f,e)
    P=s.kronecker_product(L(l),L(k))
    r=c/(2*((l/k)**2-1))*(((l/k)**2+1)*s.kronecker_product(H,H)+2*l/k*Pi)
    # {u,v}=c*u*v/2 where c=-b^2 epsilon/2.
    PB=s.zeros(4)
    for i,j,a,b in itertools.product(range(2),repeat=4):
        PB[2*i+a,2*j+b]=c*u*v/2*(s.diff(L(l)[i,j],u)*s.diff(L(k)[a,b],v)-s.diff(L(l)[i,j],v)*s.diff(L(k)[a,b],u))
    check('A1 quadratic bracket, all 16 entries',zmat(PB-r*P+P*r))
    R=s.Matrix(4,4,s.symbols('r:16'))
    for aa,bb in [(2,3),(2,2),(2,-2)]:
        equations=[]
        for expr in R*P.subs({l:aa,k:bb})-P.subs({l:aa,k:bb})*R:
            equations+=s.Poly(s.expand(expr*u**2*v**2),u,v).coeffs()
        M,_=s.linear_eq_to_matrix(equations,list(R)); dim=16-M.rank()
        check('A1 commutant dimension '+str((aa,bb)),dim==(1 if aa*aa!=bb*bb else 2),dim)
        if bb==aa:
            C=s.diag(0,1,1,0)
            check('Reject reported exceptional generator diag(0,1,1,0)',not zmat(C*P.subs({l:aa,k:bb})-P.subs({l:aa,k:bb})*C))
        if bb==-aa:
            C=s.diag(0,1,-1,0)
            check('Reject reported exceptional generator diag(0,1,-1,0)',not zmat(C*P.subs({l:aa,k:bb})-P.subs({l:aa,k:bb})*C))

def spin3():
    # Dimensionless Phi=beta phi, Psi=beta psi; m retained symbolically.
    p,q,pt,qt,ptt,qtt,m,sg=s.symbols('p q pt qt ptt qtt m sg', real=True, nonzero=True)
    # use p=2(a+b),q=2(a-b) at the final trigonometric simplification.
    px=qt-m*(sg*s.sin((p+q)/2)+s.sin((p-q)/2)/sg)
    qx=pt+m*(sg*s.sin((p+q)/2)-s.sin((p-q)/2)/sg)
    def dt(x): return s.diff(x,p)*pt+s.diff(x,q)*qt+s.diff(x,pt)*ptt+s.diff(x,qt)*qtt
    a,b=s.symbols('a b',real=True)
    for sign in [1,-1]:
        dp=pt+sign*px; dq=qt+sign*qx
        ddp=2*ptt+2*sign*dt(px)+m*m*s.sin(p)
        ddq=2*qtt+2*sign*dt(qx)+m*m*s.sin(q)
        jp=sign*(m*m*s.cos(p)*dp**2-ddp**2+dp**4/4)/(8*m**3)
        jq=sign*(m*m*s.cos(q)*dq**2-ddq**2+dq**4/4)/(8*m**3)
        C=s.cos((p+sign*q)/2)
        D=sign*sg**sign/(2*m*m)*C*dp*dq+sg**(3*sign)*(C-s.Rational(2,3)*C**3)
        expr=s.expand((jq-jp+dt(D)).subs({p:a+b,q:a-b}))
        expr=s.trigsimp(s.expand_trig(expr),method='fu')
        check('Spin '+str(-3*sign)+' full defect jet identity',expr==0)
        # The bulk continuity relation reduces to d_- T4=d_+ X2.
        y,y1,y2=s.symbols('y y1 y2')
        dminusT=2*y2*(-m*m*s.cos(y)*y1)-y1**3*(-m*m*s.sin(y))
        dplusX=s.diff(-m*m*s.cos(y)*y1*y1,y)*y1+s.diff(-m*m*s.cos(y)*y1*y1,y1)*y2
        check('Spin '+str(-3*sign)+' bulk continuity',s.expand(dminusT-dplusX)==0)
    # Static soliton profile after y=m(x cosh theta-t sinh theta).
    for theta in [-0.7,0.0,1.1]:
        for sign in [1,-1]:
            def density(y):
                sech=1/np.cosh(y); der1=2*sech*np.exp(-sign*theta)
                der2=-2*sech*np.tanh(y)*np.exp(-2*sign*theta)
                cosine=1-2*sech**2
                return (der2**2-der1**4/4+cosine*der1**2)/(8*np.cosh(theta))
            val=quad(density,-35,35,epsabs=1e-12)[0]
            target=-2/3*np.exp(-3*sign*theta)
            check('Soliton spin charge '+str((theta,sign)),abs(val-target)<1e-10,f'{val:.12g}')

def matrix_general():
    rng=np.random.default_rng(6001)
    worst=0.
    for h in range(2,6):
        S=np.roll(np.eye(h),-1,axis=0) # S e_j=e_{j-1}
        Q=np.eye(h)-np.ones((h,h))/h
        D=2*np.linalg.pinv(np.eye(h)-S); E=Q-D
        rho=np.zeros((h,h))
        for j in range(h): rho[j,(1-j)%h]=-1
        check('Cayley and reflection h='+str(h),np.max(abs(D+D.T-2*Q))<1e-12 and np.max(abs(rho@D@rho@S.T+D))<1e-12)
        for coupling in [0.87,0.87j]:
            b=coupling; m=1.13; sg=0.71; lam=1.17+0.32j
            ph=Q@rng.normal(size=h); ps=Q@rng.normal(size=h)
            vt=Q@rng.normal(size=h); wt=Q@rng.normal(size=h)
            fj=np.exp(b*(ps-np.roll(ph,-1))); gj=np.exp(b*(ph-ps))
            for typ in [-1,1]:
                DD=D if typ==-1 else -D; EE=E if typ==-1 else -E
                a=sg if typ==-1 else 1/sg
                Bp=m/b*(-a*np.roll(fj,1)+gj/a)
                Bq=m/b*(a*fj-gj/a)
                Bp=Q@Bp; Bq=Q@Bq
                px=EE@vt+DD@wt-Bp; qx=DD.T@vt-EE@wt+Bq
                left=np.exp(-b*(EE@ph+DD@ps)/2)
                right=np.exp(b*(DD.T@ph-EE@ps)/2)
                tilde=np.eye(h,dtype=complex)
                if typ==-1: tilde+=sg/lam*S.T
                else:
                    for j in range(1,h): tilde+=(lam/sg)**j*np.linalg.matrix_power(S,j)
                L=left[:,None]*tilde*right[None,:]
                for pm in [1,-1]:
                    pd=vt+pm*px; qd=wt+pm*qx
                    dL=(-b*(EE@pd+DD@qd)/2)[:,None]*L+L*(b*(DD.T@pd-EE@qd)/2)[None,:]
                    Ap=pm*b/2*np.diag(pd)+pm*m*lam**pm*(np.diag(np.exp(b*(ph-np.roll(ph,-1))/2))@S if pm==1 else S.T@np.diag(np.exp(b*(ph-np.roll(ph,-1))/2)))
                    Aq=pm*b/2*np.diag(qd)+pm*m*lam**pm*(np.diag(np.exp(b*(ps-np.roll(ps,-1))/2))@S if pm==1 else S.T@np.diag(np.exp(b*(ps-np.roll(ps,-1))/2)))
                    worst=max(worst,np.max(abs(dL-L@Aq+Ap@L)))
                if typ==1:
                    M=(1-(lam/sg)**h)*np.linalg.inv(L)
                    wanted=np.diag(np.sqrt(gj))-lam/sg*np.diag(np.sqrt(fj))@S
                    # exp roots, not branch-dependent sqrt, in this small test.
                    check('Type I+ sparse inverse h='+str(h)+' b='+str(b),np.max(abs(M-wanted))<1e-10)
                # Verify action involution f -> g and g -> f.
                newph=rho@ph; newps=rho@S.T@ps
                ff=np.exp(b*(newps-np.roll(newph,-1))); gg=np.exp(b*(newph-newps))
                check('Action involution h='+str(h)+' type='+str(typ)+' b='+str(b),np.max(abs(ff-gj[(-np.arange(h))%h]))<1e-10 and np.max(abs(gg-fj[(-np.arange(h))%h]))<1e-10)
    check('All-rank both-type both-coupling gauge checks',worst<1e-10,f'max {worst:.3e}')

def diagrams():
    data=[]
    for n in range(1,5):
        eye=np.eye(n+1); roots=eye-np.roll(eye,-1,axis=0)
        data.append((f'a{n}1',roots))
    def B(n):
        e=np.eye(n); return np.array([-e[0]-e[1]]+[e[j]-e[j+1] for j in range(n-1)]+[e[-1]])
    def C(n):
        e=np.eye(n); return np.array([-2*e[0]]+[e[j]-e[j+1] for j in range(n-1)]+[2*e[-1]])
    for n in [3,4]: data.append((f'b{n}1',B(n)))
    for n in [2,3,4]: data.append((f'c{n}1',C(n)))
    e=np.eye(4)
    data.append(('d41',np.array([-e[0]-e[1],e[0]-e[1],e[1]-e[2],e[2]-e[3],e[2]+e[3]])))
    F=np.array([-e[0]-e[1],e[1]-e[2],e[2]-e[3],e[3],(e[0]-e[1]-e[2]-e[3])/2])
    a=np.array([1.,-1,0]); b=np.array([-2.,1,1]); G=np.array([-3*a-2*b,a,b])
    data.extend([('f41',F),('g21',G)])
    for n in range(1,5):
        e=np.eye(n); roots=np.array([-2*e[0]]+[e[j]-e[j+1] for j in range(n-1)]+[e[-1]])
        data.append((f'a{2*n}2',roots))
    def dual(R): return R/(R*R).sum(axis=1)[:,None]
    for n in [2,3,4]: data.append((f'a{2*n-1}2',dual(B(n))))
    for n in [3,4]: data.append((f'd{n+1}2',dual(C(n))))
    data.extend([('e62',dual(F)),('d43',dual(G))])
    counts=[0,0,0,0]; bad=[]; total=0
    for name,R in data:
        gram=R@R.T; h=len(R); rank=np.linalg.matrix_rank(R)
        for p in itertools.permutations(range(h)):
            if np.max(abs(gram-gram[np.ix_(p,p)]))>1e-10: continue
            total+=1
            if any(p[i]==i for i in range(h)): counts[0]+=1; continue
            # Work on root span using SVD; solve R P^T=R[p].
            _,_,vh=np.linalg.svd(R,full_matrices=False); basis=vh[:rank].T
            roots=R@basis; P=np.linalg.lstsq(roots,roots[list(p)],rcond=None)[0].T
            if np.linalg.matrix_rank(np.eye(rank)-P,tol=1e-9)<rank: counts[1]+=1; continue
            D=2*np.linalg.inv(np.eye(rank)-P); bil=roots@D@roots.T
            forbidden=sum(abs(bil[i,j])>1e-9 for i in range(h) for j in range(h) if j not in [i,p[i]])
            if forbidden: counts[2]+=1; bad.append((name,p,forbidden))
            else: counts[3]+=1
    check('Affine diagram scan',len(data)==23 and total==78 and counts==[64,5,2,7],f'23 diagrams; {total} symmetries; {counts}; failures={bad}')

def scalar_and_scattering():
    P,Q,a0,a1,b0,b1=s.symbols('P Q a0 a1 b0 b1')
    product=2*(a1*(P+Q)+a0)*(b1*(P-Q)+b0)
    vl=2*a1*b1*P**2+2*(a1*b0+a0*b1)*P+2*a0*b0
    vr=2*a1*b1*Q**2+2*(-a1*b0+a0*b1)*Q
    check('General polynomial scalar branch',s.expand(product-vl+vr)==0)
    check('Opposite-slope counterexample to old replacement',s.expand(product.subs({a1:1,a0:0,b1:0,b0:1})-2*P-2*Q)==0)
    check('Nonzero constant product counterexample',product.subs({a1:0,b1:0,a0:1,b0:1})==2)
    t,sg=s.symbols('t sg',positive=True)
    z=(1+sg*t)/(1-sg*t)
    check('Delay two equations',s.simplify((1-z)/t+sg*(1+z))==0 and s.simplify(-t*(1+z)-(1-z)/sg)==0)
    T=(1+s.I*sg*t)/(1-s.I*sg*t)
    check('Wave unitarity and delay continuation',s.simplify(T*s.conjugate(T)-1)==0 and s.simplify(z.subs(t,s.I*t)-T)==0)
    check('Reject reported z=+coth',float(z.subs({t:2,sg:1}))==-3)

if __name__=='__main__':
    rank_one(); spin3(); matrix_general(); diagrams(); scalar_and_scattering()
