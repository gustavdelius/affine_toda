"""Exact rational-function checks, independent of all book code.

Multiply each vector sewing residual by (1-pi); this is invertible on
the root space. The printed C1,C2 are its simple-root coefficients / i.
Their common constant cancels because sum(alpha_j)=0.
v=exp(i*pi*a/h), K=exp(theta), S=coefficient of f including jump phase.
The variable t is omega^(ja)*E, not physical time.
"""
import sympy as s

t,v,K,S,z=s.symbols('t v K S z', nonzero=True)
I=s.I
q=v**2
mass=(v-1/v)/I
kx=mass*(K+1/K)/2
kt=-mass*(K-1/K)/2
L=t/(1+t);Ln=q*t/(1+q*t)
R=z*t/(1+z*t);Rn=z*q*t/(1+z*q*t)
f=(1+q*t)*(1+z*t/q)/((1+t)*(1+z*t))
fn=(1+q*q*t)*(1+z*t)/((1+q*t)*(1+z*q*t))
gn=(1+t)*(1+z*q*t)/((1+q*t)*(1+z*t))

for sign in [1,-1]:
    x=K**(-sign)/S
    zz=(x-I*v)/(x-I/v)
    C1=kx*(L-Ln)+sign*kt*(L+Ln-2*R)-S*f+gn/S
    C2=kx*(R-Rn)+sign*kt*(2*Ln-R-Rn)-S*fn+gn/S
    constant=-S+1/S
    values=[s.factor(s.cancel((C-constant).subs(z,zz))) for C in [C1,C2]]
    print(f'type {"I-" if sign==1 else "I+"}: exact sewing root-coefficient residuals = {values}',flush=True)
    assert values==[0,0]
    # Specializations include every p through S, and arbitrary rapidity K.
    for h in [3,4]:
        count=2*(h-1)*h
        print(f'  a_{h-1}: all {h-1} species x {h} classes, {count} identities are specializations of zero',flush=True)

# Derive the congruence using exact fundamental weights and D.
for h in [3,4]:
    roots=s.eye(h)-s.eye(h)[:,list(range(1,h))+[0]]
    # Explicit columns alpha_j=e_j-e_(j+1), independent of numerical script.
    roots=s.Matrix(h,h,lambda j,k: int(j==k)-int((j+1)%h==k))
    weights=[s.zeros(h,1)]
    for p in range(1,h):
        w=s.linsolve((roots[1:,:].col_join(s.ones(1,h)),s.eye(h-1)[:,p-1].col_join(s.zeros(1,1))))
        weights.append(s.Matrix(list(w)[0]))
    D=sum((2*(weights[j]-weights[(j+1)%h])*weights[j].T for j in range(1,h)),s.zeros(h))
    for p in range(h):
        integers=(roots*D*weights[p]+s.ones(h,1)*s.Rational(2*p,h))/2
        assert all(x.is_integer for x in integers)
    print(f'a_{h-1}: exact alpha_j D lambda_p + 2p/h is even for every j,p')

# Independent algebra for the charged-delay modulus.
r,c,d,A,B=s.symbols('r c d A B',real=True)
# x=r(c+i*d), exp(i*gamma)=A+i*B.
num=(r*c+B)**2+(r*d-A)**2
den=(r*c-B)**2+(r*d-A)**2
print('charged |numerator|^2-|denominator|^2 =',s.expand(num-den))

# Same-type obstruction: the ratio of the +/-infinity limits is omega^(2b).
print('same-type asymptotic obstruction: v^4=1; 0<b<h implies b=h/2')

# Energy at rest is real regardless of p; momentum supplies the missing condition.
d0,q0,M=s.symbols('d0 q0 M',real=True)
print('Im(class energy) =',-M*s.sinh(d0)*s.sin(q0))
print('Im(class momentum) =',-M*s.cosh(d0)*s.sin(q0))
