"""Six-vertex R-matrix of U_q(sl2 hat): chapters 7 and 8.

Solves Jimbo's equations (sec-jimbo-equations) for two spin-1/2 evaluation
representations (eq-evaluation-rep) and checks every printed property of
eq-six-vertex: the matrix itself, unitarity (eq-r-unitarity), the braid
Yang-Baxter equation in the printed index order (eq-ybe-braid), the spin-0
eigenvalue (eq-six-vertex-eigenvalue), the crossing data of sec-qg-crossing
(x_c = q^2, C, c_u) and the sine-Gordon dictionary of sec-six-vertex.
"""
import sympy as sp

q, x, y = sp.symbols('q x y', nonzero=True)
E = sp.Matrix([[0, 1], [0, 0]]); F = sp.Matrix([[0, 0], [1, 0]])
K = sp.diag(q, 1/q); I2 = sp.eye(2); kp = sp.kronecker_product
P = sp.Matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


def gens(z):
    return {'e1': E, 'f1': F, 'k1': K, 'e0': z*F, 'f0': E/z, 'k0': K.inv()}


def cop(gx, gy, i):
    e, f, k = 'e'+i, 'f'+i, 'k'+i
    return {e: kp(gx[e], gx[k]) + kp(I2, gy[e]),
            f: kp(gx[f], I2) + kp(gx[k].inv(), gy[f]),
            k: kp(gx[k], gy[k])}


# 1. Jimbo's equations
syms = sp.symbols('r0:16'); R = sp.Matrix(4, 4, syms); eqs = []
for i in '01':
    A, B = cop(gens(x), gens(y), i), cop(gens(y), gens(x), i)
    for key in A:
        eqs += list(R*A[key] - B[key]*R)
sol = sp.solve(eqs, syms, dict=True)[0]
Rj = sp.simplify((R.subs(sol)/R.subs(sol)[0, 0]).subs(y, 1))


def Rc(u):  # eq-six-vertex
    return sp.Matrix([[1, 0, 0, 0],
                      [0, u*(q**2-1)/(q**2-u), q*(u-1)/(u-q**2), 0],
                      [0, q*(u-1)/(u-q**2), (q**2-1)/(q**2-u), 0],
                      [0, 0, 0, 1]])


report('eq-six-vertex solves Jimbo equations', sp.simplify(Rj - Rc(x)) == sp.zeros(4))
report('unitarity Rc(z)Rc(1/z)=1', sp.simplify(Rc(x)*Rc(1/x)) == sp.eye(4))
x1, x2, x3 = sp.symbols('x1 x2 x3', nonzero=True)
A = lambda M: kp(M, I2); B = lambda M: kp(I2, M)
lhs = A(Rc(x2/x3))*B(Rc(x1/x3))*A(Rc(x1/x2)); rhs = B(Rc(x1/x2))*A(Rc(x1/x3))*B(Rc(x2/x3))
report('braid YBE in printed order', sp.simplify(lhs - rhs) == sp.zeros(8))
ev = {sp.factor(k): v for k, v in Rc(x).eigenvals().items()}
report('eigenvalues 1 (x3) and (1-q^2 z)/(z-q^2)',
       ev.get(1) == 3 and any(sp.simplify(k - (1-q**2*x)/(x-q**2)) == 0 for k in ev))

# 2. crossing (eq-r-crossing)
Rr = lambda u: P*Rc(u)


def t1(M):
    N = sp.zeros(4)
    for a in range(2):
        for b in range(2):
            for c in range(2):
                for d in range(2):
                    N[2*c+b, 2*a+d] = M[2*a+b, 2*c+d]
    return N


C = sp.Matrix([[0, 1], [-q, 0]]); cu = q*(x-1)/(q**2*x-1)
report('crossing: x_c=q^2, C=[[0,1],[-q,0]], c_u=q(z-1)/(q^2 z-1)',
       sp.simplify(t1(Rr(x).inv()) - cu*kp(C.inv(), I2)*Rr(x*q**2)*kp(C, I2)) == sp.zeros(4))

# 3. sine-Gordon dictionary: z=e^{2 lam th}, q=-e^{i pi lam}
lam, th = sp.symbols('lambda theta', positive=True)
sub = {x: sp.exp(2*lam*th), q: -sp.exp(sp.I*sp.pi*lam)}
S = Rr(x).subs(sub)
ST = sp.sinh(lam*th)/sp.sinh(lam*(sp.I*sp.pi - th))
SR = sp.sinh(sp.I*sp.pi*lam)/sp.sinh(lam*(sp.I*sp.pi - th))
num = lambda e: complex(sp.N(e.subs({lam: sp.Rational(37, 23), th: sp.Rational(3, 7)})))
report('b = S_T/S_0', abs(num(S[1, 1] - ST)) < 1e-12)
report('off-diagonal = e^{-+lam th} S_R/S_0',
       abs(num(S[1, 2] - sp.exp(-lam*th)*SR)) < 1e-12 and abs(num(S[2, 1] - sp.exp(lam*th)*SR)) < 1e-12)
# gauge D(theta_i)=diag(e^{-lam th_i/2}, e^{lam th_i/2}) on each particle removes e^{-+lam theta}
t1_, t2_ = sp.symbols('t1 t2', real=True)
D = lambda t: sp.diag(sp.exp(-lam*t/2), sp.exp(lam*t/2))
Sg = (kp(D(t1_), D(t2_))**-1)*S.subs(th, t1_-t2_)*kp(D(t1_), D(t2_))
vals = {lam: sp.Rational(37, 23), t1_: sp.Rational(5, 7), t2_: sp.Rational(2, 7)}
SRn = SR.subs(th, t1_-t2_)
report('gauge to principal gradation gives S_R/S_0 in both off-diagonal entries',
       abs(complex(sp.N((Sg[1, 2] - SRn).subs(vals)))) < 1e-12 and abs(complex(sp.N((Sg[2, 1] - SRn).subs(vals)))) < 1e-12)
