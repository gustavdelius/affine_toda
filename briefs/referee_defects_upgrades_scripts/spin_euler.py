"""Independent a3 calculation in a rational Hadamard basis; no book imports.

Uses Laurent monomials z_i=exp(field_i/2), not exponent simplification or
the book's homotopy primitive. Euler derivatives test exactness directly.
"""
import sympy as S
from itertools import product
from time import monotonic

start=monotonic()
H=S.Matrix([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/2
roots=[(H.row(j)-H.row((j+1)%4)).T for j in range(4)]
P=S.zeros(4)
for j in range(4): P[(j-1)%4,j]=1
D=2*(S.eye(3)-H.T*P*H).inv()
E=S.eye(3)-D
z=S.symbols('z0:6', nonzero=True)
vel=S.symbols('v0:6'); acc=S.symbols('a0:6'); jerk=S.symbols('j0:6')
sig=S.symbols('sigma', nonzero=True)
def expand(x): return S.expand(x)
def dq(x,i): return z[i]*S.diff(x,z[i])/2
def dt(x):
    return expand(sum(dq(x,i)*vel[i]+S.diff(x,vel[i])*acc[i]+S.diff(x,acc[i])*jerk[i] for i in range(6)))
def mono(c):
    assert all((2*v).is_Integer for v in c), c
    return S.prod(z[i]**(2*c[i]) for i in range(6))
def bulk_exp(a,side): return mono(list(a if side==0 else S.zeros(3,1))+list(a if side==1 else S.zeros(3,1)))
def cubic(u): return 2*sum((H.row(k)*u)[0]**3 for k in range(4))

# Coefficient form of d_-T-d_+Theta; tests normalization of the current.
u=S.Matrix(S.symbols('u0:3'))
bulk=[expand(3*2*sum((H.row(k)*u)[0]**2*(H.row(k)*a)[0] for k in range(4))
             +2*(a.T*(3*E)*u)[0]*(a.T*u)[0]) for a in roots]
assert bulk==[0]*4
bad_bulk=[expand(3*2*sum((H.row(k)*u)[0]**2*(H.row(k)*a)[0] for k in range(4))
                 +2*(a.T*E*u)[0]*(a.T*u)[0]) for a in roots]
assert any(bad_bulk)
print('bulk coefficient identities: 4 exact zeros; B=E control fails',flush=True)

def flux(typ,eps,Dpot):
    Ep=S.eye(3)-Dpot
    sigma=sig if typ==1 else 1/sig
    potential=0
    for a in roots:
        cf=list(Dpot*a/2)+list(Dpot.T*a/2)
        cg=list(Dpot.T*a/2)+list(-Dpot.T*a/2)
        potential+=sigma*mono(cf)+mono(cg)/sigma
    gp=S.Matrix([dq(potential,i) for i in range(3)])
    gq=S.Matrix([dq(potential,i+3) for i in range(3)])
    vp=S.Matrix(vel[:3]); vq=S.Matrix(vel[3:])
    Fx=[typ*(Ep*vp+Dpot*vq)-gp,typ*(Dpot.T*vp-Ep*vq)+gq]
    G=0
    for side,sign in [(0,1),(1,-1)]:
        v=S.Matrix(vel[3*side:3*side+3]); a2=S.Matrix(acc[3*side:3*side+3])
        grad=sum((a*bulk_exp(a,side) for a in roots),S.zeros(3,1))
        light=v+eps*Fx[side]
        light2=2*a2+2*eps*Fx[side].applyfunc(dt)+grad
        theta=sum((a.T*(3*E)*light)[0]*bulk_exp(a,side) for a in roots)
        G+=sign*(cubic(light)+(light2.T*(3*E)*light)[0]+theta)
    return expand(G)

for typ,eps in product([1,-1],[-1,1]):
    G=flux(typ,eps,D)
    euler=[expand(dq(G,i)-dt(S.diff(G,vel[i]))+dt(dt(S.diff(G,acc[i])))) for i in range(6)]
    assert euler==[0]*6
    # A nonzero spin-three cubic invariant ensures this is a nontrivial bulk current.
    print(f'type={typ:+d}, epsilon={eps:+d}, spin={-2*eps:+d}: six Euler derivatives = 0; '
          f'flux has {len(S.Add.make_args(G))} Laurent monomials; {monotonic()-start:.1f}s',flush=True)

# Preserves D+D^T=2 but violates the Toda compatibility, so this is a meaningful control.
Dbad=D+S.Matrix([[0,1,0],[-1,0,0],[0,0,0]])
Gbad=flux(1,1,Dbad)
control=expand(dq(Gbad,0)-dt(S.diff(Gbad,vel[0]))+dt(dt(S.diff(Gbad,acc[0]))))
assert control!=0
point={**dict(zip(z,[1,2,3,4,5,6])),**dict(zip(vel,[2,-1,3,1,-2,4])),
       **dict(zip(acc,[1,2,-1,0,3,-2])),**dict(zip(jerk,[0,1,0,2,0,3])),sig:S.Rational(7,5)}
print('nonintegrable D: Euler derivative 0 is nonzero; value at rational point =',control.subs(point),flush=True)
print('PASS: independent exact a3 flux test, both types/chiralities, symbolic nonzero sigma.',flush=True)
