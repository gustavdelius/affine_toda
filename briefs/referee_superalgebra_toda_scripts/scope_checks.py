"""Small independent algebra checks for source-scope questions."""
import sympy as s
def check(label,ok):
    print(('PASS ' if ok else 'FAIL ')+label,flush=True)
    assert ok
x,ell=s.symbols('x ell');y=s.Function('y')(x)
# Ito-Zhu 3.7 lower two rows: (d+l/x)y=z, (d-l/x)z=p^2 y.
z=s.diff(y,x)+ell*y/x
eliminated=s.diff(z,x)-ell*z/x
check('ODE elimination gives centrifugal term l(l+1)/x^2',s.simplify(eliminated-s.diff(y,x,2)+ell*(ell+1)*y/x**2)==0)
# Nepomechie 19 consistency, prior to fixing the remaining spin-three constant.
phi,alpha,phi0,C=s.symbols('phi alpha phi0 C',real=True)
N=alpha*s.sin((phi-phi0)/2)-4*s.sin(phi/2)
den=alpha*s.cos((phi-phi0)/2)-4*s.cos(phi/2)+C
g=s.sqrt(8*den);f=N/g
check('boundary supercharge equations g*f=N and dg/dphi=-2f',s.simplify(g*f-N)==0 and s.simplify(s.diff(g,phi)+2*f)==0)
# D(2,1;a) kinetic signature and degenerations, source field coordinates.
a=s.symbols('a',real=True)
for name,value,signs in [('positive',2,(1,-1,-1)),('middle',-s.Rational(1,2),(-1,1,-1)),('negative',-2,(-1,-1,1))]:
    kinetic=[-1/(1+a),s.Integer(1),1/a]
    result=[(c*t).subs(a,value) for c,t in zip(kinetic,signs)]
    check('exceptional twisted kinetic signature '+name,all(c>0 for c in result) or all(c<0 for c in result))
# Free-fermion threshold: C-family nth boson has M=2 sqrt(2).
check('C nth boson is at two-fermion threshold, u=0',s.simplify(8-4*2*s.cos(0/2)**2)==0)
# Gandenberger N=1 factor-by-factor identity after h=2n+1, lambda=1/h.
n,j,u=s.symbols('n j u');h=2*n+1;lam=1/h
source=[u+j-lam,u+j-2*n*lam+1,u+j-(n+s.Rational(1,2))*lam,u+j-(n+s.Rational(1,2))*lam+1,
        -u+j-(2*n+1)*lam,-u+j+1,-u+j-(n+s.Rational(3,2))*lam,-u+j-(n-s.Rational(1,2))*lam+1]
target=[u+j-1/h,u+j+1/h,u+j-s.Rational(1,2),u+j+s.Rational(1,2),-u+j-1,-u+j+1,-u+j-1/h-s.Rational(1,2),-u+j+1/h+s.Rational(1,2)]
check('N=1 Gamma numerator arguments coincide; denominator follows u->-u',all(s.simplify(v-w)==0 for v,w in zip(source,target)))
# Scalar free field need not decouple: cross kinetic action du dv+ exp(v) W(q).
q,v=s.symbols('q v');W=s.exp(v)*s.cosh(q)
check('gl extension can source interacting fields through a free null coordinate',s.diff(W,v,q)!=0)
