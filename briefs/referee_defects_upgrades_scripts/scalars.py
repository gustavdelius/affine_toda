"""Independent scalar algebra: arbitrary Fourier eigenvalue and fusing phases.
No book imports. Gamma calculations use explicit factors read from the PDFs.
"""
import sympy as S
import numpy as np
from scipy.special import loggamma
import cmath

c,p,d,m=S.symbols('c p d m',nonzero=True)
R,T=S.symbols('R T')
I=S.I
# c=e^(i gamma), p=e^theta, d=e^theta_d. Hessian from f,g covectors.
rot=c*c
D=2/(1-rot); DT=2/(1-1/rot); E=1-D
mass=m*(c-1/c)/I
energy=mass*(p+1/p)/2; momentum=mass*(p-1/p)/2
for typ in (1,-1):
    sigma=1/d if typ==1 else d
    expected=(d/p+c)/(d/p+1/c) if typ==1 else (p/d-c)/(p/d-1/c)
    for side in ('left','right'):
        phi,psi=(1+R,T) if side=='left' else (T,1+R)
        phix,psix=(I*momentum*(1-R),I*momentum*T) if side=='left' else (-I*momentum*T,-I*momentum*(1-R))
        gradphi=m*(sigma*(phi-psi/rot)+(phi-psi)/sigma)
        gradpsi=m*(sigma*(psi-rot*phi)-(phi-psi)/sigma)
        eq1=phix+I*typ*energy*(E*phi+D*psi)+gradphi
        eq2=psix+I*typ*energy*(DT*phi-E*psi)-gradpsi
        mat,vec=S.linear_eq_to_matrix([S.cancel(eq1),S.cancel(eq2)],[R,T])
        values=mat.inv()*vec
        sol={R:S.cancel(values[0]),T:S.cancel(values[1])}
        target=expected if side=='left' else 1/expected.subs(p,1/p)
        assert S.cancel(sol[R])==0 and S.cancel(sol[T]-target)==0
        print(f'linear type={typ:+d}, {side}: R=0; transmission difference=0')

x,u,v=S.symbols('x u v',nonzero=True)
Tm=lambda x,u:(x+u)/(x+1/u)
Tp=lambda x,u:(1-x*u)/(1-x/u)
for name,F in [('minus',Tm),('plus',Tp)]:
    identities=[F(x,u)*F(-x,-1/u)-1,
                F(x/v,u)*F(x*u,v)-F(x,u*v),
                F(-x*v,u)*F(-x/u,v)-F(x,-u*v)]
    assert [S.cancel(a) for a in identities]==[0]*3
    print(f'particle type={name}: crossing and both fusing identities = [0,0,0]')

F=(-I*x-1/u)/(-I*x-u)
Fbar=(-I*x+u)/(-I*x+1/u)
assert S.cancel(Fbar-F.subs({I:-I,u:1/u},simultaneous=True))==0
norm=S.cancel(F*F.subs({I:-I,u:1/u},simultaneous=True))
assert S.cancel(norm-(x*x+1+I*x*(u-1/u))/(x*x+1-I*x*(u-1/u)))==0
print('TBA |F|^2 = (x^2+1-2x sin(gamma))/(x^2+1+2x sin(gamma)); F_half=(x-1)/(x+1)')

def classical(D,h):
    u=cmath.exp(1j*np.pi/h); x=cmath.exp(-D)
    return (x+u)/(x+1/u)
def simple(D):return (cmath.sinh(D)+.3j)/(cmath.sinh(D)-.3j)
for h in (3,4,5,6,8):
    a=np.prod([classical(.41-1j*np.pi+2j*np.pi*k/h,h) for k in range(h)])
    b=np.prod([simple(.41-1j*np.pi+2j*np.pi*k/h) for k in range(h)])
    assert abs(a-1)<1e-12
    print(f'h={h} CDD product: f=Tcl_1 error {abs(a-1):.2e}; simple sinh factor {b:.10f}')
for omega in (np.sqrt(2),np.sqrt(5)):
    h=4; a=2; ub=np.pi/2-np.pi/(h*omega)
    rapidities=[.41-1j*ub+1j*np.pi*(2*k-h+a+1)/h for k in range(h-a)]
    rapidities += [.41+1j*ub+1j*np.pi*(2*k-a+1)/h for k in range(a)]
    print(f'B1^(2) multiplier, omega={omega:.10f}: {np.prod([classical(D,h) for D in rapidities]):.10f}')

def lg(z):return loggamma(np.asarray(z,dtype=complex))
def gammahat(z,a,K):
    k=np.arange(K,dtype=float)
    ans=0j
    for b,sgn in ((1,1),(2,1),(4,-1),(5,-1)):
        base=.5+(6*k+b)*a
        ans+=sgn*np.sum(lg(base+z)+lg(base-z)-2*lg(base))
    return ans
def extrap(z,a,K):
    # Cancel 1/K and 1/K^2 tails with three truncations.
    L=[gammahat(z,a,K*2**j) for j in range(3)]
    return (L[0]-6*L[1]+8*L[2])/3
for xi in (.93,1.73):
    a=np.pi/(3*xi); q=1j*cmath.exp(1j*np.pi*a)
    for K in (256,1024):
        C=.5*cmath.exp(-2*extrap(a,a,K))
        def rho(z):return C*cmath.exp(1j*np.pi*z+extrap(z,a,K))
        bs=[]; cr=[]
        for z0 in (.11+.07j,-.23+.31j,.05-.4j):
            xx=cmath.exp(-2j*np.pi*z0)
            bs.append(abs((1+xx)*rho(z0-a)*rho(z0+a)/rho(z0)-1))
            cr.append(abs(rho(z0)*rho(z0+3*a)*(1+q**-4*xx)*(1-q**-2*xx)-1))
        print(f'CZ2011 xi={xi}, K={K}, common C={C:.9g}: full bootstrap error={max(bs):.2e}, crossing error={max(cr):.2e}')
        assert max(bs+cr)<3e-6

# Exact finite-product telescoping as a multiset of formal Gamma arguments.
from collections import Counter
def gamma_terms(terms):
    out=Counter()
    for base,signz,power in terms:out[(S.expand(base),signz)]+=power
    return {k:v for k,v in out.items() if v}
w,zsym=S.symbols('omega z')
for h in (3,4,7):
    for K in (1,4,9):
        A=lambda k:S.Rational(1,2)+(h*k+S.Rational(h-1,2))*w
        printed=[(A(0),-1,1)]
        equal=[]
        for k in range(1,K+1):printed += [(A(k),-1,1),(A(k-1),1,1),(A(k-1)+w,-1,-1),(A(k-1)+w,1,-1)]
        for k in range(K):equal += [(A(k),-1,1),(A(k),1,1),(A(k)+w,-1,-1),(A(k)+w,1,-1)]
        equal.append((A(K),-1,1))
        assert gamma_terms(printed)==gamma_terms(equal)
print('CZ2009a: exact finite-product identity including Gamma(A_K-z), h=3,4,7; K=1,4,9')
print('PASS. Algebraic derivations for all h and the Gamma endpoint are in the report.')
