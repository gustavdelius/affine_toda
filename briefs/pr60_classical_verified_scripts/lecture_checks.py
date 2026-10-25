"""Fresh lecture and higher-rank checks, definitions only."""
import itertools
import numpy as np
import sympy as s
from independent_checks import check, zmat

def bianchi():
    u,v,w,a,b,d=s.symbols('u v w a b d',nonzero=True)
    C=(a+b)/(a-b); y=C*(v-w)/(1+v*w); z=s.cancel((u+y)/(1-u*y))
    def sinus(t,r,sgn):
        return s.factor(2*(t+sgn*r)*(1-sgn*t*r)/((1+t*t)*(1+r*r)))
    for pm in [1,-1]:
        dv=d-2*a*sinus(u,v,1) if pm==1 else -d-2/a*sinus(u,v,-1)
        dw=d-2*b*sinus(u,w,1) if pm==1 else -d-2/b*sinus(u,w,-1)
        dz=sum(s.factor(s.diff(z,x)*(1+x*x)/(1+z*z))*der for x,der in [(u,d),(v,dv),(w,dw)])
        for f,der,par in [(v,dv,b),(w,dw,a)]:
            residual=der-dz-2*par*sinus(f,z,1) if pm==1 else der+dz+2/par*sinus(f,z,-1)
            check('Bianchi full arbitrary-background edge '+str((pm,par)),s.cancel(residual)==0)

def geometry():
    x,t=s.symbols('x t',real=True)
    r=s.Matrix([s.sech(x)*s.cos(t),s.sech(x)*s.sin(t),x-s.tanh(x)])
    rx=r.diff(x); rt=r.diff(t)
    check('Pseudosphere metric ss',s.trigsimp(rx.dot(rx)-s.tanh(x)**2)==0)
    check('Pseudosphere metric st',s.trigsimp(rx.dot(rt))==0)
    check('Pseudosphere metric tt',s.trigsimp(rt.dot(rt)-s.sech(x)**2)==0)
    # Determinants of first and second fundamental forms; avoid normalization sqrt.
    cross=rx.cross(rt); EG=s.simplify(rx.dot(rx)*rt.dot(rt)-rx.dot(rt)**2)
    K=s.simplify((cross.dot(r.diff(x,2))*cross.dot(r.diff(t,2))-cross.dot(r.diff(x,t))**2)/EG**2)
    check('Pseudosphere curvature away from rim',s.trigsimp(K+1)==0)
    # At x=0 rx=0: explicitly exclude this singular locus.
    check('Pseudosphere rim singularity',all(v==0 for v in rx.subs(x,0)))
    alpha=s.symbols('alpha',positive=True)
    # SG for tangent-quarter variable Y: (1+Y²)(Ytt-Yxx)-2Y(Yt²-Yx²)+Y(1-Y²)=0.
    Y=s.cot(alpha)*s.sin(t*s.sin(alpha))/s.cosh(x*s.cos(alpha))
    resid=(1+Y**2)*(s.diff(Y,t,2)-s.diff(Y,x,2))-2*Y*(s.diff(Y,t)**2-s.diff(Y,x)**2)+Y*(1-Y**2)
    resid=s.trigsimp(s.expand_trig(resid))
    check('Breather PDE',resid==0)

def lecture_lax():
    p,q,dp,dq,l,m,sg,c=s.symbols('p q dp dq l m sg c',nonzero=True,real=True)
    # Algebra in exponential variables e^{i p/4},e^{i q/4}.
    u,v=s.symbols('u v',nonzero=True)
    H=s.diag(1,-1); X=s.Matrix([[0,1],[1,0]])
    L=s.Matrix([[u/v,sg/(l*u*v)],[sg*u*v/l,v/u]])
    def deriv(M,a,b): return M.diff(u)*s.I*u*a/4+M.diff(v)*s.I*v*b/4
    sinplus=(u*u*v*v-1/(u*u*v*v))/(2*s.I)
    sinminus=(u*u/(v*v)-v*v/(u*u))/(2*s.I)
    for pm in [1,-1]:
        Ap=pm*s.I/4*dp*H+pm*m*l**pm/2*s.Matrix([[0,u**(2*pm)],[u**(-2*pm),0]])
        Aq=pm*s.I/4*dq*H+pm*m*l**pm/2*s.Matrix([[0,v**(2*pm)],[v**(-2*pm),0]])
        residual=deriv(L,dp,dq)-L*Aq+Ap*L
        subst={dq:dp+2*sg*m*sinplus} if pm==1 else {dq:-dp+2*m/sg*sinminus}
        check('Lecture6 Lax lightcone '+str(pm),zmat(residual.subs(subst)))
    # Bulk r-matrix bracket (lecture dimensionless field, beta²=1).
    l1,l2,vel=s.symbols('l1 l2 vel',nonzero=True)
    def ax(lam): return s.I*vel*H/4+m/4*(lam*s.Matrix([[0,u*u],[u**-2,0]])+lam**-1*s.Matrix([[0,u**-2],[u*u,0]]))
    e=s.Matrix([[0,1],[0,0]]); f=e.T
    z=l1/l2; rr=-((z*z+1)*s.kronecker_product(H,H)+4*z*(s.kronecker_product(e,f)+s.kronecker_product(f,e)))/(16*(z*z-1))
    A=ax(l1); B=ax(l2); PB=s.zeros(4)
    for i,j,k,n in itertools.product(range(2),repeat=4):
        PB[2*i+k,2*j+n]=s.I*u/4*(s.diff(A[i,j],u)*s.diff(B[k,n],vel)-s.diff(A[i,j],vel)*s.diff(B[k,n],u))
    total=s.kronecker_product(A,s.eye(2))+s.kronecker_product(s.eye(2),B)
    check('Lecture6 bulk r sign all16',zmat(PB-rr*total+total*rr))
    check('Lecture6 determinant',s.simplify(L.det()-(1-sg**2/l**2))==0)

def higher_rank_pb():
    rng=np.random.default_rng(602)
    worst=0.
    for h in range(2,6):
        # Orthonormal basis of traceless hyperplane.
        raw=np.eye(h)[:,:h-1]-np.eye(h)[:,[-1]]
        basis=np.linalg.qr(raw)[0]; rnk=h-1
        S=np.roll(np.eye(h),-1,axis=0)
        D=basis.T@(2*np.linalg.pinv(np.eye(h)-S))@basis
        E=np.eye(rnk)-D
        for typ in [-1,1]:
            DD=D if typ==-1 else -D; EE=E if typ==-1 else -E
            sympl=np.block([[EE,DD],[-DD.T,EE]]); bracket=np.linalg.inv(sympl)
            coeffL=-np.hstack([basis@EE,basis@DD])/2
            coeffR=np.hstack([basis@DD.T,-basis@EE])/2
            coord=rng.normal(size=2*rnk); sg=0.71; l1=1.13+0.27j;l2=0.83-0.11j
            for b in [0.92,0.92j]:
                def make(lam):
                    if typ==-1: tt=np.eye(h,dtype=complex)+sg/lam*S.T
                    else: tt=sum((lam/sg)**j*np.linalg.matrix_power(S,j) for j in range(h))
                    coeff=b*(coeffL[:,None,:]+coeffR[None,:,:])
                    mat=np.exp(coeff@coord)*tt
                    return mat,coeff*mat[:,:,None]
                L,dL=make(l1); M,dM=make(l2)
                pb=np.einsum('ija,ab,klb->ikjl',dL,bracket,dM).reshape(h*h,h*h)
                z=l1/l2
                cartan=np.diag((np.eye(h)-np.ones((h,h))/h).reshape(h*h))
                rr=(z**h+1)*cartan.astype(complex)
                for i in range(h):
                    for j in range(i+1,h):
                        rr[i*h+j,j*h+i]+=2*z**(j-i)
                        rr[j*h+i,i*h+j]+=2*z**(h-j+i)
                rr*=b*b/(4*(z**h-1))
                product=np.kron(L,M)
                residual=np.max(abs(pb-rr@product+product@rr))/(1+np.max(abs(pb)))
                worst=max(worst,residual)
    check('Higher-rank Poisson bracket h2..5 both types/couplings',worst<1e-11,f'max relative {worst:.3e}')

if __name__=='__main__':
    bianchi();geometry();lecture_lax();higher_rank_pb()
