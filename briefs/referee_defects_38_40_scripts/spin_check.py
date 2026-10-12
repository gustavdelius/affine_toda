"""Independent spin-two derivation. No imports from the book's code.

Use d_e = d_t + e d_x (without sqrt(2)).  T=C(d_e phi)^3+
(d_e^2 phi)^T B d_e phi, Theta=-(d_e phi)^T B grad V.
The conserved quantity is integral(T-Theta) - e F, where d_t F is
the left minus right flux T+Theta.  Thus book spin is -2e.
"""
from pathlib import Path
import itertools
import numpy as np
import sympy as s
from scipy.integrate import quad_vec

HERE=Path(__file__).resolve().parent
U=s.Matrix([[1/s.sqrt(2),1/s.sqrt(6)],[-1/s.sqrt(2),1/s.sqrt(6)],[0,-2/s.sqrt(6)]])
roots=[U.T*(s.eye(3)[:,j]-s.eye(3)[:,(j+1)%3]) for j in range(3)]
P=U.T*s.Matrix([[0,1,0],[0,0,1],[1,0,0]])*U
D=s.simplify(2*(s.eye(2)-P).inv()); B=s.eye(2)-D
coords=s.symbols('p r v w'); vel=s.symbols('p1 r1 v1 w1'); acc=s.symbols('p2 r2 v2 w2')
phi=s.Matrix(coords[:2]); psi=s.Matrix(coords[2:]); vp=s.Matrix(vel[:2]); vs=s.Matrix(vel[2:])
sig=s.symbols('sigma',nonzero=True)
Dt=lambda expr: sum(s.diff(expr,x)*y for xs,ys in ((coords,vel),(vel,acc)) for x,y in zip(xs,ys))
grad=lambda expr,x:s.Matrix([s.diff(expr,t) for t in x])
V=lambda f:sum(s.exp((a.T*f)[0])-1 for a in roots)
u=s.Matrix(s.symbols('u0:2'))
c=s.symbols('c0:4')
cubic=c[0]*u[0]**3+3*c[1]*u[0]**2*u[1]+3*c[2]*u[0]*u[1]**2+c[3]*u[1]**3
eq=[]
for a in roots:
    eq+=s.Poly(s.expand((grad(cubic,u).T*a)[0]-2*(u.T*B*a)[0]*(a.T*u)[0]),*u).coeffs()
sol=s.solve(eq,c); cubic=s.factor(cubic.subs(sol))
C=lambda v:cubic.subs(dict(zip(u,v)),simultaneous=True)
print('D =',D,'B =',B,flush=True)
print('C(u) =',cubic,'; BCZ tensor contraction A(u,u,u)=3 C(u)',flush=True)
print('C - (2/3)sum(weights.u)^3 =',s.simplify(C(u)-s.Rational(2,3)*sum(t**3 for t in U*u)),flush=True)
print('C(Du)+C(D^T u) =',s.simplify(C(D*u)+C(D.T*u)),flush=True)
print('C(Du) =',s.factor(C(D*u)),flush=True)

def build(e,typ=1,scale=1):
    Db=s.eye(2)-scale*B
    Eb=s.eye(2)-Db
    sd=sig if typ==1 else 1/sig
    pot=sum(sd*(s.exp((a.T*(Db.T*phi+Db*psi))[0]/2)-1)
            +(s.exp((a.T*Db*(phi-psi))[0]/2)-1)/sd for a in roots)
    xp=typ*(Eb*vp+Db*vs)-grad(pot,phi)
    xs=typ*(Db.T*vp-Eb*vs)+grad(pot,psi)
    G=0
    for sign,f,v,a,x in ((1,phi,vp,s.Matrix(acc[:2]),xp),(-1,psi,vs,s.Matrix(acc[2:]),xs)):
        light=v+e*x
        second=2*a+2*e*s.Matrix([Dt(t) for t in x])+grad(V(f),f)
        G+=sign*(C(light)+(second.T*B*light)[0]-(light.T*B*grad(V(f),f))[0])
    return s.expand(G),xp,xs

def correction(G):
    # Integrate the acceleration coefficient first, then a closed one-form
    # along the coordinate axes.  All cancellations are exact SymPy zeros.
    def tidy(expr):
        table={}
        for term in s.Add.make_args(s.expand(expr)):
            ex=0; coeff=1
            for factor in s.Mul.make_args(term):
                if factor.func==s.exp:ex+=factor.args[0]
                else:coeff*=factor
            ex=s.expand(ex)
            table[ex]=table.get(ex,0)+coeff
        return s.Add(*(s.factor(co)*s.exp(ex) for ex,co in table.items() if s.expand(co)!=0))
    print('  collecting acceleration coefficient',flush=True)
    az=[tidy(s.diff(G,a)) for a in acc]
    assert all(s.diff(x,v)==0 for x in az for v in vel)
    Fv=sum(x*v for x,v in zip(az,vel))
    rem=tidy(G-Dt(Fv))
    bz=[tidy(s.diff(rem,v)) for v in vel]
    assert tidy(rem-sum(x*v for x,v in zip(bz,vel)))==0
    assert all(tidy(s.diff(bz[i],coords[j])-s.diff(bz[j],coords[i]))==0 for i in range(4) for j in range(i))
    print('  integrating closed field one-form',flush=True)
    F0=0
    for j,q in enumerate(coords):
        integrand=tidy(bz[j]-s.diff(F0,q))
        for term in s.Add.make_args(s.expand(integrand)):
            if term==0:continue
            derivative=tidy(s.diff(term,q))
            ratio=s.cancel(derivative/term)
            F0+=term/ratio if ratio!=0 else term*q
    F0=tidy(F0)
    F=Fv+F0
    residual=tidy(G-Dt(F))
    assert residual==0
    return F,az,F0

def soliton(x,t,theta=.9,xi=.2+.45j,z=1):
    # b=1 with complex fields; analytic continuation gives b=i beta.
    R=np.asarray(roots,dtype=float).reshape(3,2)
    om=np.exp(2j*np.pi*np.arange(3)/3)
    a=np.sqrt(3)*np.cosh(theta); d=-np.sqrt(3)*np.sinh(theta)
    k=z*np.exp(a*x+d*t+xi)*om
    f=-R.T@np.log(1+k)
    first=-R.T@(k/(1+k)); second=-R.T@(k/(1+k)**2)
    return f,d*first,d*d*second,a*first,first,second

def exact_solution_checks(e,F,xp,xs,typ):
    # Time reversal turns an I_- solution with sigma'=1/sigma into I_+.
    sigma=.71; sigma0=sigma if typ==1 else 1/sigma
    theta=.9; a=np.pi/3; z=(np.exp(-theta)-1j*sigma0*np.exp(1j*a))/(np.exp(-theta)-1j*sigma0*np.exp(-1j*a))
    fn=s.lambdify((*coords,*vel,sig),F,'numpy',cse=True)
    sew=s.lambdify((*coords,*vel,sig),s.Matrix.vstack(xp,xs),'numpy',cse=True)
    cb=s.lambdify(tuple(u),cubic,'numpy')
    Br=np.asarray(B,dtype=float); R=np.asarray(roots,dtype=float).reshape(3,2)
    bulk=[]; corrected=[]; defect=[]; sewerr=0
    for t in (-1.3,-.5,0.,.8,1.6):
        l=soliton(0,typ*t,theta,z=1); r=soliton(0,typ*t,theta,z=z)
        args=(*l[0],*r[0],*(typ*l[1]),*(typ*r[1]),sigma)
        sewerr=max(sewerr,np.max(np.abs(np.asarray(sew(*args)).ravel()-np.r_[l[3],r[3]])))
        def density(x,side):
            f,fd,fdd,fx,first,second=soliton(x,typ*t,theta,z=1 if side==0 else z)
            coeff=typ*(-np.sqrt(3)*np.sinh(theta))+e*np.sqrt(3)*np.cosh(theta)
            light=coeff*first; light2=coeff**2*second
            gv=R.T@np.exp(R@f)
            return cb(*light)+light2@Br@light+light@Br@gv
        qleft=quad_vec(lambda x:density(x,0),-24,0,epsabs=1e-10,epsrel=1e-10)[0]
        qright=quad_vec(lambda x:density(x,1),0,24,epsabs=1e-10,epsrel=1e-10)[0]
        Q=qleft+qright; bd=-e*fn(*args)
        bulk.append(Q); defect.append(bd); corrected.append(Q+bd)
    spread=lambda seq:float(np.max(np.abs(np.asarray(seq)-seq[0])))
    print(f'  soliton: sewing={sewerr:.2e}, bulk spread={spread(bulk):.5g}, corrected spread={spread(corrected):.2e}',flush=True)
    assert sewerr<1e-10 and spread(corrected)<1e-8

if __name__=='__main__':
    for typ,e in itertools.product((1,-1),(-1,1)):
        G,xp,xs=build(e,typ)
        F,az,F0=correction(G)
        name=f'type_{typ:+d}_spin_{-2*e:+d}'
        (HERE/(name+'.txt')).write_text('B_defect = '+s.sstr(-e*F)+'\n\nLaTeX:\n'+s.latex(-e*F)+'\n')
        print(name+': exact G-dF/dt = 0; local correction written',flush=True)
        print('  acceleration coefficient =',az,flush=True)
        print('  field-only primitive =',F0,flush=True)
        exact_solution_checks(e,F,xp,xs,typ)
    # Control has the same antisymmetry conditions but a non-integrable D.
    # Evaluate the Euler derivative at constant nonzero velocities, zero
    # accelerations. A nonzero value rules out any local primitive.
    for e in (-1,1):
        G,_,_=build(e,1,s.Rational(3,5))
        third=s.symbols('p3 r3 v3 w3'); fourth=s.symbols('p4 r4 v4 w4')
        jets=[coords,vel,acc,third,fourth]
        d=lambda ex:sum(s.diff(ex,x)*y for xs,ys in zip(jets,jets[1:]) for x,y in zip(xs,ys))
        EL=[s.diff(G,q)-d(s.diff(G,v))+d(d(s.diff(G,a))) for q,v,a in zip(coords,vel,acc)]
        f=s.lambdify((*sum((list(x) for x in jets),[]),sig),EL,'numpy',cse=True)
        vals=[.1,-.2,.3,.15,.4,-.1,.2,.5]+[0.]*12+[.71]
        err=np.max(np.abs(f(*vals)))
        print(f'nonintegrable D, spin {-2*e:+d}: Euler residual={err:.6g}',flush=True)
        assert err>1e-4
