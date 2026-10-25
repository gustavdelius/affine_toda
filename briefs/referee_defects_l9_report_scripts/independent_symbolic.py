"""Independent checks from printed definitions; no book code imported."""
import sympy as S
I=S.I
p,d,b,m,s,v,a=S.symbols('p d beta m sigma v a', nonzero=True)
# p=exp(i beta phi/2), d=exp(i beta phi_D/2).
cp=(p*d+1/(p*d))/2; cm=(p/d+d/p)/2
sp=(p*d-1/(p*d))/(2*I); sm=(p/d-d/p)/(2*I)
sinphi=(p**2-p**-2)/(2*I); cosphi=(p**2+p**-2)/2
sinD=(d**2-d**-2)/(2*I)
def dxfield(f): return S.diff(f,p)*I*b*p/2
def dt(f): return dxfield(f)*v+S.diff(f,v)*a
F=-m/b*(s*sp+sm/s) # phi_x
W=v+m/b*(s*sp-sm/s) # psi_x
Dm=s*b*b/(2*m*m)*cp*(v+F)*W+s**3*(cp-S.Rational(2,3)*cp**3)
Dp=b*b/(2*m*m*s)*cm*(v-F)*W+s**-3*(cm-S.Rational(2,3)*cm**3)
C=-b/m*sinD*W
def j(e):
    first=v+e*F
    second=2*a+2*e*dxfield(F)*v+m*m/b*sinphi
    return e*b*b/(8*m**3)*(m*m*cosphi*first**2-second**2+b*b/4*first**4)
def check(name,z):
    z=S.simplify(S.factor(z))
    print(('PASS ' if z==0 else 'FAIL ')+name+': '+str(z),flush=True)
check('fused spin3, symbolic beta,m,sigma,phi_D',dt(Dm+Dp+C)-j(1)-j(-1))
check('omitting Dirichlet term has predicted nonzero residual',dt(Dm+Dp)-j(1)-j(-1)-b/m*sinD*dt(W))
print('omitted-term residual = beta/m sin(beta phi_D) [a + (m/2)(sigma c_plus-c_minus/sigma) v]',flush=True)
check('Dirichlet endpoint current', b*b/(8*m**3)*((m*m/b*sinD-2*a)**2-(m*m/b*sinD+2*a)**2)+b/m*sinD*a)
# Constant normalization of sigma's Gamma pairs, using exponentials as symbols.
L,X,H=S.symbols('L X H',nonzero=True)
check('four Gamma exponentials',H*(X**-1+X-X*L-X**-1/L)+H*(L**S.Rational(1,2)-L**-S.Rational(1,2))*(X*L**S.Rational(1,2)-X**-1*L**-S.Rational(1,2)))
y,z,E=S.symbols('y yd exp_i_eta',nonzero=True)
A=lambda t:S.Matrix([[1,-I*t],[-I*t,1]])
kp=y*E+1/(y*E); km=y/E+E/y
M=A(1/(y*z)).inv()*S.diag(kp,km)*A(y/z)
cu=(y+1/y)/2; su=(y-1/y)/(2*I)
ce=(E+1/E)/2; se=(E-1/E)/(2*I)
ch=(z+1/z)/2; sh=(z-1/z)/2
G=S.Matrix([[cu*ce*ch-su*se*sh,-su*cu],[-su*cu,cu*ce*ch+su*se*sh]])
Lam=S.diag(1,-E)
check('sector matrix and conjugation',sum(list((M-4/z/(1+1/(y*z)**2)*Lam.inv()*G*Lam).applyfunc(S.factor))))
x,r=S.symbols('x r',nonzero=True)
T=(1+x*r)/(1+x/r); U=(1-x/r)/(1-x*r)
check('folding ratio',T/U-(1-x*x*r*r)/(1-x*x/r**2))
check('opposite parameters transparent',T*U.subs(x,-x)-1)
check('equal parameters block',T*U-((x-1/x)-(r-1/r))/((x-1/x)+(r-1/r)))
check('continued folded soliton delay',(T*U).subs(x,-I*x)-((x+1/x)+(r-1/r)/I)/((x+1/x)-(r-1/r)/I))
bb,theta=S.symbols('b_ly theta_b')
check('Lee-Yang exponent',(-I*(bb+3)*S.pi/5).subs(bb,3+6*I*theta/S.pi)-S.Rational(6,5)*(theta-I*S.pi))
lam,B=S.symbols('lambda B')
check('boost continuation',(1/(1-1/(lam+1))).subs(lam,-2/B)-(1-B/2))
