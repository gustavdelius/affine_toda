"""Independent referee calculations; no imports from the book's check scripts.

Exterior algebra: sorted generator tuples, left derivatives, coefficients on left.
Formal superjets u_k = D_+^k Phi have parity k mod 2; D_+ is odd.
"""
import sympy as s
import numpy as np
import mpmath as mp

I=s.I
def clean(a):
    return {k:s.expand(v) for k,v in a.items() if s.expand(v)!=0}
def add(*args):
    out={}
    for a in args:
        for k,v in a.items(): out[k]=out.get(k,0)+v
    return clean(out)
def scale(a,c): return clean({k:c*v for k,v in a.items()})
def mul(a,b):
    out={}
    for k,v in a.items():
        for l,w in b.items():
            if set(k)&set(l): continue
            sign=(-1)**sum(x>y for x in k for y in l)
            t=tuple(sorted(k+l));out[t]=out.get(t,0)+sign*v*w
    return clean(out)
def gen(k):return {(k,):s.Integer(1)}
def scalar(x):return {():s.sympify(x)}
def left(a,k):
    out={}
    for t,c in a.items():
        if k in t:
            j=t.index(k);out[t[:j]+t[j+1:]]=(-1)**j*c
    return clean(out)
def zero(a):return all(s.simplify(v.rewrite(s.exp))==0 for v in a.values())
def check(label,ok,detail=''):
    print(('PASS ' if ok else 'FAIL ')+label+(' '+str(detail) if str(detail) else ''),flush=True)
    if not ok:raise AssertionError(label)

# Superspace coefficient calculation, independent exterior generators.
t1,t2,p,q,pminus,qplus=map(gen,range(6))
phip,phim,F,mu,beta,a,phi=s.symbols('phip phim F mu beta a phi',real=True)
Dp=add(scale(p,I),scale(t1,-I*phip),scale(t2,I*F),mul(mul(t1,t2),pminus))
# pminus here denotes d_-p; its occurrence in D+ is unused by the top kinetic term.
Dm=add(scale(q,I),scale(t1,-I*F),scale(t2,-I*phim),scale(mul(mul(t1,t2),qplus),-1))
# Correct jets in top components: D+Phi has d_+q, D-Phi has -d_-p.
Dp=add(scale(p,I),scale(t1,-I*phip),scale(t2,I*F),mul(mul(t1,t2),qplus))
Dm=add(scale(q,I),scale(t1,-I*F),scale(t2,-I*phim),scale(mul(mul(t1,t2),pminus),-1))
def integrate(expr):
    return clean({k[2:]:I*v for k,v in expr.items() if k[:2]==(0,1)})
kin=integrate(scale(mul(Dp,Dm),I/2))
expected=add(scalar((phip*phim+F**2)/2),scale(mul(p,pminus),I/2),scale(mul(q,qplus),I/2))
check('Berezin kinetic coefficients including +F^2/2',zero(add(kin,scale(expected,-1))),kin)
nil=scale(add(mul(t1,p),mul(t2,q),scale(mul(t1,t2),F)),I*a*beta)
potential=integrate(scale(add(scalar(1),nil,scale(mul(nil,nil),s.Rational(1,2))),mu*s.exp(a*beta*phi)/beta**2))
pot_expected=add(scalar(-mu*a*F*s.exp(a*beta*phi)/beta),scale(mul(p,q),I*mu*a*a*s.exp(a*beta*phi)))
check('odd exponential: auxiliary minus, Yukawa plus',zero(add(potential,scale(pot_expected,-1))))
check('even spurion gives negative potential',integrate(scale(mul(t1,t2),I*mu/beta**2))==scalar(-mu/beta**2))
m=s.symbols('m',positive=True)
W=4*m*s.cosh(beta*phi/s.sqrt(2))/beta**2
V=s.diff(W,phi)**2/2
check('SSG mass squared and Yukawa mass',s.simplify(s.diff(V,phi,2).subs(phi,0)-4*m*m)==0 and s.diff(W,phi,2).subs(phi,0)==2*m)
check('SSG bosonic truncation',s.simplify(s.expand_trig(V)-2*m*m*(s.cosh(s.sqrt(2)*beta*phi)-1)/beta**2)==0)

# Matrix Chevalley relations.
def E(i,j):
    a=s.zeros(3);a[i-1,j-1]=1;return a
h=s.diag(1,-1,0)/2
e=[-(E(2,3)+E(3,1))/2,E(1,3)+E(3,2)]
f=[E(1,3)-E(3,2),(E(3,1)-E(2,3))/2]
check('all SSG Chevalley relations',all(e[i]*f[j]+f[j]*e[i]==((-1 if i==0 else 1)*h if i==j else s.zeros(3)) for i in range(2) for j in range(2)) and all(h*e[i]-e[i]*h==(-1 if i==0 else 1)*e[i]/2 for i in range(2)))
check('osp ansatz obstruction', (E(3,1)-E(2,3))**2!=s.zeros(3))

# On-shell superjet calculus. f=Phi, u2,u4,... even; u1,u3,... exterior.
z=s.symbols('z');M=s.symbols('M');ev={0:z,**{k:s.symbols('u'+str(k)) for k in range(2,10,2)}}
def jet(k):return scalar(ev[k]) if k%2==0 else gen(k)
def deriv(expr,images,parity):
    out={}
    for odd,c in expr.items():
        for k,x in ev.items():
            dc=s.diff(c,x)
            if dc!=0:out=add(out,mul(scale(images[k],dc),{odd:1}))
        for j,k in enumerate(odd):
            v=mul(mul({odd[:j]:1},images[k]),{odd[j+1:]:1})
            out=add(out,scale(v,c*(-1)**(parity*j)))
    return clean(out)
plus={k:jet(k+1) for k in range(8)}
def P(expr):return deriv(expr,plus,1)
# i D+D-Phi=-M sinh(Phi), hence D-D+Phi=-i M sinh(Phi).
minus={1:scalar(-I*M*s.sinh(z))}
for k in range(2,8):minus[k]=scale(P(minus[k-1]),-1)
# Currents do not depend on z except through partners, so image 0 unused in D-.
minus[0]={}
T=add(mul(jet(1),jet(2)),scale(jet(3),-1))
Tb=scale(mul(jet(1),scalar(s.exp(-z))),I*M)
J=mul(T,P(T))
Jb=scale(mul(add(mul(jet(1),mul(jet(2),jet(2))),scale(mul(jet(2),jet(3)),-1)),scalar(s.exp(-z))),I*M)
for label,current,partner in [('superspin 3/2',T,Tb),('superspin 7/2',J,Jb)]:
    residual=add(deriv(current,minus,1),scale(P(partner),-1))
    check(label+' current on shell',zero(residual))
check('printed partner -1/2 of correct first partner fails',not zero(add(deriv(T,minus,1),scale(P(Tb),s.Rational(1,2)))))

# Root-derived masses, cubic and quartic couplings, tree amplitude.
mp.mp.dps=35
for n in (2,3):
    basis=np.eye(n)
    for family in ('A','C'):
        if family=='A':
            roots=[-basis[0]-basis[1]]+[basis[j]-basis[j+1] for j in range(n-1)]+[2*basis[-1]]
            coeff=[1]+[1]+[2]*(n-2)+[1]
        else:
            roots=[-2*basis[0]]+[basis[j]-basis[j+1] for j in range(n-1)]+[2*basis[-1]]
            coeff=[1]+[2]*(n-1)+[1]
        H=sum(c*np.outer(r,r) for c,r in zip(coeff,roots))
        vals,vecs=np.linalg.eigh(H)
        expected=sorted(([2]+[8*np.sin(np.pi*j/(2*n-1))**2 for j in range(1,n)]) if family=='A' else [8*np.sin(np.pi*j/(2*n))**2 for j in range(1,n+1)])
        check(f'{family} n={n} Hessian masses',np.max(abs(vals-expected))<1e-12,vals)
        if family=='C':continue
        v=vecs[:,np.argmin(abs(vals-2))]
        c4=sum(c*(r@v)**4 for c,r in zip(coeff,roots))
        c3=np.array([sum(c*(r@v)**2*(r@vecs[:,j]) for c,r in zip(coeff,roots)) for j in range(n)])
        for theta in (0.71,1.19):
            mandel_s=4+4*np.cosh(theta);mandel_u=8-mandel_s
            amp=-c4-sum(g*g*(1/(mandel_s-mass)-1/mass+1/(mandel_u-mass)) for g,mass in zip(c3,vals))
            tree=1j*amp/(8*np.sinh(theta)) # coefficient of beta_DGPZ^2
            hh=2*n-1
            def logPder(th):
                eta=-1j*hh*th/(2*mp.pi)
                return (mp.digamma(eta/hh)-mp.digamma(eta/hh+mp.mpf('.5')))/hh
            bder=(logPder(theta)+logPder(1j*mp.pi-theta)-2j*mp.pi/mp.sinh(hh*theta))/(4*mp.pi)
            check(f'A n={n} tree fixes gamma negative theta={theta}',abs(tree-bder)<1e-12,f'error={float(abs(tree-bder)):.2g}; C4={c4:.8g}')

def ratios(th,n,omega):
    hh=2*n-1;g=-np.sin(np.pi*omega)
    R=1j*g/np.sinh(hh*th/2);D=(-1)**(n+1)*1j*g/np.cosh(hh*th/2)
    B=1+2j*g/np.sinh(hh*th);F=-1+2j*g/np.sinh(hh*th)
    return np.array([[B,0,0,-D],[0,R,1,0],[0,1,R,0],[D,0,0,F]],complex)
for n in (2,3):
    omega=np.sqrt(2)/5;x=.61;y=.83
    A=ratios(x,n,omega);B=ratios(y,n,omega);C=ratios(x+y,n,omega);eye=np.eye(2)
    lhs=np.kron(B,eye)@np.kron(eye,C)@np.kron(A,eye)
    rhs=np.kron(eye,A)@np.kron(C,eye)@np.kron(eye,B)
    check(f'n={n} braid YBE',np.max(abs(lhs-rhs))<1e-12,np.max(abs(lhs-rhs)))
    factor=1+(np.sin(np.pi*omega)/np.sinh((2*n-1)*x/2))**2
    check(f'n={n} algebraic unitarity',np.max(abs(A@ratios(-x,n,omega)-factor*np.eye(4)))<1e-12)
    Bv,Fv,Rv,Dv=A[0,0],A[3,3],A[1,1],A[3,0]
    check(f'n={n} free fermion correct sign',abs(Bv*Fv+1+Dv**2-Rv**2)<1e-12)
    check(f'n={n} printed B.1 fails',abs(Bv*Fv+1-Dv**2-Rv**2)>1e-5)
    check(f'n={n} Majorana amplitude crossing',abs(ratios(1j*np.pi-x,n,omega)[3,0]-1j*Rv)<1e-12)

def block(x,th,h,B):
    def elementary(q):return np.sinh(th/2+1j*np.pi*q/(2*h))/np.sinh(th/2-1j*np.pi*q/(2*h))
    return elementary(x-1)*elementary(x+1)/(elementary(x-1+B)*elementary(x+1-B))
for n in (2,3,4,5):
    def channels(th):
        starts=(1,3) if n%2==0 else (3,1)
        return [-np.prod([block(k,th,2*n,np.sqrt(2)/2) for k in range(start,2*n,4)]) for start in starts]
    a=np.array(channels(.91));b=np.array(channels(-.91));c=np.array(channels(1j*np.pi-.91))
    check(f'C n={n} fermion unitarity/crossing',max(abs(a*b-1))<1e-12 and max(abs(c-(a[::-1] if n%2==0 else a)))<1e-12)

# Exceptional roots, operator masses and canonical cubic Yukawa sign.
alpha=s.symbols('alpha',positive=True);G=s.diag(-(1+alpha)/2,s.Rational(1,2),alpha/2)
roots=[s.Matrix(r) for r in [(-1,-1,-1),(-1,1,1),(1,-1,1),(1,1,-1)]]
MF=sum((r*r.T*G for r in roots),s.zeros(3))
check('D(2,1;alpha) null roots and fermion eigenvalues',all((r.T*G*r)[0]==0 for r in roots) and MF==s.diag(-2*(1+alpha),2,2*alpha))
# For twisted theory source 3.32 cubic cross terms +2m barpsi_i psi_j phi_k.
# Reverse overall negative signature and normalize using (sqrt(1+a),1,sqrt(a)).
normal=s.sqrt(1+alpha)*s.sqrt(alpha)
check('D-family canonical cubic coefficient is negative',s.simplify(-2*(2*s.Symbol('mu'))*normal+4*s.Symbol('mu')*s.sqrt(alpha*(1+alpha)))==0)
print('Independent algebraic checks complete.')
