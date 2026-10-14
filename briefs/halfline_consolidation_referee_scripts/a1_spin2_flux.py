"""Referee check of the spin +-2 remark after prp-no-dirichlet (ch. 30).

a_2 Toda, m=beta=1, V = sum_j exp(alpha_j.phi), x^pm=(x^0 pm x^1)/sqrt2,
equation d_+d_- phi = -V'/2.  Solve d_- T3 = d_+ Theta1 for
T3 = (1/3)A_abc d phi_a d phi_b d phi_c + B_ab d^2 phi_a d phi_b,
Theta1 = d_+ phi_a G_a(phi), then evaluate the boundary flux of
Q2 + c Q_{-2} at a Dirichlet wall (d_0^k phi(0)=0, d_1^2 phi(0)=V'(phi_0)).
"""
import sympy as sp

p1, p2 = sp.symbols('p1 p2')            # phi components
phi = sp.Matrix([p1, p2])
al = [sp.Matrix([sp.sqrt(2), 0]), sp.Matrix([-1 / sp.sqrt(2), sp.sqrt(sp.Rational(3, 2))])]
al = [-(al[0] + al[1])] + al             # alpha_0, alpha_1, alpha_2
ex = [sp.exp((a.T * phi)[0]) for a in al]
V = sum(ex)
Vp = sp.Matrix([sp.diff(V, x) for x in (p1, p2)])

# jets: d_+ phi = u, d_+^2 phi = w  (independent symbols for the + chirality)
u = sp.Matrix(sp.symbols('u1 u2'))
w = sp.Matrix(sp.symbols('w1 w2'))
# unknowns
A = {}
for i in range(2):
    for j in range(i, 2):
        for k in range(j, 2):
            A[(i, j, k)] = sp.Symbol(f'A{i}{j}{k}')


def Aget(i, j, k):
    return A[tuple(sorted((i, j, k)))]


b = sp.Symbol('b')
B = sp.Matrix([[0, b], [-b, 0]])
g = sp.symbols('g0:6')
G = sp.Matrix([g[0] * ex[0] + g[1] * ex[1] + g[2] * ex[2], g[3] * ex[0] + g[4] * ex[1] + g[5] * ex[2]])

T3 = sum(sp.Rational(1, 3) * Aget(i, j, k) * u[i] * u[j] * u[k] for i in range(2) for j in range(2) for k in range(2)) \
    + sum(B[i, j] * w[i] * u[j] for i in range(2) for j in range(2))

# d_- acting: d_- phi = v (other chirality, appears only via d_- of functions of phi)
v = sp.Matrix(sp.symbols('v1 v2'))
# d_- u = -V'/2 ; d_- w = d_+(-V'/2) = -(1/2) V'' u
Vpp = sp.hessian(V, (p1, p2))
dm_u = -Vp / 2
dm_w = -(Vpp * u) / 2


def dminus(expr):
    r = 0
    for i in range(2):
        r += sp.diff(expr, u[i]) * dm_u[i] + sp.diff(expr, w[i]) * dm_w[i]
        r += sp.diff(expr, [p1, p2][i]) * v[i]
    return r


def dplus(expr):
    # d_+ of a function of phi and u: phi -> u, u -> w
    r = 0
    for i in range(2):
        r += sp.diff(expr, [p1, p2][i]) * u[i] + sp.diff(expr, u[i]) * w[i]
    return r


Theta1 = (u.T * G)[0]
eq = sp.expand(dminus(T3) - dplus(Theta1))
# collect coefficients of monomials in u, w, v and exponentials
E = sp.symbols('E0:3')
eqE = eq
for k in range(3):
    pass
# substitute exponentials by monomials in independent symbols
x1, x2 = sp.symbols('x1 x2', positive=True)
subs = {p1: sp.log(x1), p2: sp.log(x2)}
eqs = sp.Poly(sp.expand(sp.simplify(eq.subs(subs))), *u, *w, *v).coeffs()
# each coeff is a combination of x1^a x2^b with irrational exponents -> collect by expanding powsimp
conds = []
for c in eqs:
    c = sp.expand(sp.powsimp(sp.expand(c), force=True))
    terms = sp.Add.make_args(c)
    byexp = {}
    for t in terms:
        coeff, rest = t.as_independent(x1, x2)
        byexp[rest] = byexp.get(rest, 0) + coeff
    conds += list(byexp.values())
unk = list(A.values()) + [b] + list(g)
sol = sp.solve(conds, unk, dict=True)
print("solutions:", sol)
s = sol[0]
free = [x for x in unk if x not in s]
print("free parameters:", free)
# pick normalization: set free parameters generic
vals = {x: 1 for x in free}
T3s = T3.subs(s).subs(vals)
Th1 = Theta1.subs(s).subs(vals)
print("T3 nonzero cubic:", sp.simplify(T3s.subs({w[0]: 0, w[1]: 0})) != 0, " B:", sp.simplify(B.subs(s).subs(vals)))

# Dirichlet flux.  At x=0: d_0 phi=0, d_0^2 phi=0, d_1 phi = U, d_0 d_1 phi = Ud, d_1^2 phi = V'.
U = sp.Matrix(sp.symbols('U1 U2'))
Ud = sp.Matrix(sp.symbols('Ud1 Ud2'))
r2 = sp.sqrt(2)
plus = {u[0]: U[0] / r2, u[1]: U[1] / r2, w[0]: (2 * Ud[0] + Vp[0]) / 2, w[1]: (2 * Ud[1] + Vp[1]) / 2}
minus = {u[0]: -U[0] / r2, u[1]: -U[1] / r2, w[0]: (-2 * Ud[0] + Vp[0]) / 2, w[1]: (-2 * Ud[1] + Vp[1]) / 2}
c = sp.Symbol('c')
Tm3 = T3s       # parity image: same coefficients with + -> -
flux = (T3s + Th1).subs(plus, simultaneous=True) - c * (Tm3 + Th1).subs(minus, simultaneous=True)
flux = sp.expand(flux)
# For flux = dSigma/dt with Sigma = Sigma(U): the Ud-part must be exact and the rest must vanish
Udpart = sp.Matrix([sp.diff(flux, Ud[i]) for i in range(2)])
curl = sp.simplify(sp.diff(Udpart[1], U[0]) - sp.diff(Udpart[0], U[1]))
rest = sp.simplify(flux.subs({Ud[0]: 0, Ud[1]: 0}))
print("curl of Ud-coefficient (must vanish):", sp.factor(curl))
print("Ud-free part (must vanish if c makes curl 0):", sp.factor(rest.subs(c, 1)))
print("c=-1 curl:", sp.simplify(curl.subs(c, -1)))
ok = sp.simplify(curl.subs(c, -1)) != 0 and sp.simplify(rest.subs(c, 1)) != 0 and sp.solve(curl, c) == [1]
print("PASS" if ok else "FAIL", "no spin-2 combination survives a Dirichlet wall for a_2")
