"""Checks for Lecture 1 of the defect lecture notes, "Mirrors and windows"
(lecture-notes/defects/01-mirrors-and-windows.qmd).

Written from the definitions of the lecture, independently of code/defect_axioms.py.
One PASS or FAIL line per check; exit status 1 on any FAIL.

  (a) The delta-function defect of a free boson: matching conditions give R and T; unitarity
      relations, crossing, |R|^2+|T|^2 = 1, bound state with the quantum-mechanical energy,
      the regimes of the poles.
  (b) The exchange algebra of a defect: applying the defect relation twice gives the unitarity
      relations R(t)R(-t)+T(t)Tb(-t)=1 and R(t)T(-t)+T(t)Rb(-t)=0 (symbolic rewriting).
  (c) The no-go argument: reduce A_a(t1) A_b(t2) D in two orders with a symbolic rewriting engine;
      the differences of the coefficients are R_b T_a [S_ab(t1+t2)-S_ab(t1-t2)] and
      R_a T_b [1 - S_ab(t1-t2) S_ba(t1+t2)], the both-transmitted coefficients agree identically,
      and the both-reflected ones differ by a multiple of S_ab - S_ba.
  (d) The triangle relation: index form and braid form agree for random T; T = Rcheck(theta-theta_d)
      solves it (a particle is a defect); a perturbation does not.
  (e) Defects: T(theta) T(theta + i pi) = 1 from unitarity and crossing; a bulk amplitude at a
      shifted rapidity, S(theta - theta_0), is a defect (unitarity, crossing, bootstrap).
  (f) The Lee-Yang defect T_b = f_{2/3}(theta - theta_b), theta_b = i pi (3-b)/6: axioms, the
      bound-state pole of Tbar_b, the bootstrap T_b(theta) S(theta+iu) = T_{b+2}(theta), the energies
      m sin(b pi/6) = m cosh theta_b, Hermitian analyticity exactly for Re b = 3 mod 6.
  (g) The Verlinde lines of the Lee-Yang model M(2,5): eigenvalues (1 -+ sqrt 5)/2, fusion Phi x Phi = 1 + Phi.
  (h) The Ising defect of Delfino, Mussardo and Simonetti (exercise): unitarity, crossing, g -> 4/g,
      decoupling at g = +-2.
  (i) The twist defect exp(i alpha mu) solves the triangle relation for a charge-conserving Rcheck.
  (j) The delta defect on a circle (exercise): the Bethe-Yang matrix has eigenvalues T+R and T-R = 1, and the
      even states of the direct solution satisfy e^{ipL}(T+R) = 1.

Runtime: a few seconds.
"""
import cmath
import itertools
import math

import numpy as np
import sympy as sp

OK = True


def report(name, ok):
    global OK
    OK = OK and bool(ok)
    print(("PASS " if ok else "FAIL ") + name)


rng = np.random.default_rng(1)

# ---------------------------------------------------------------------------------------------
# (a) delta-function defect of a free boson
# ---------------------------------------------------------------------------------------------
# L = 1/2 (d phi)^2 - 1/2 m^2 phi^2 - 1/2 g delta(x) phi^2  =>  phi continuous, phi'(0+) - phi'(0-) = g phi(0).
m, g, th = sp.symbols('m g theta')
p = m * sp.sinh(th)
R, T = sp.symbols('R T')
# left: e^{ipx} + R e^{-ipx}, right: T e^{ipx} (common factor e^{-iEt})
eqs = [sp.Eq(1 + R, T),                                  # continuity
       sp.Eq(sp.I * p * T - (sp.I * p - sp.I * p * R), g * T)]  # jump of the derivative
sol = sp.solve(eqs, [R, T], dict=True)[0]
Tdel = sp.simplify(sol[T])
Rdel = sp.simplify(sol[R])
report("(a) delta defect: T = 2ip/(2ip-g), R = g/(2ip-g)",
       sp.simplify(Tdel - 2 * sp.I * p / (2 * sp.I * p - g)) == 0
       and sp.simplify(Rdel - g / (2 * sp.I * p - g)) == 0)
Tm = Tdel.subs(th, -th)
Rm = Rdel.subs(th, -th)
report("(a) unitarity relations R(t)R(-t)+T(t)T(-t)=1, R(t)T(-t)+T(t)R(-t)=0 (parity: Tbar=T, Rbar=R)",
       sp.simplify(Rdel * Rm + Tdel * Tm - 1) == 0 and sp.simplify(Rdel * Tm + Tdel * Rm) == 0)
report("(a) crossing T(i pi - theta) = T(theta)",
       sp.simplify((Tdel.subs(th, sp.I * sp.pi - th) - Tdel).rewrite(sp.exp)) == 0)
vals = []
for _ in range(5):
    mm, gg, tt = 1.0, rng.normal() * 3, rng.normal() * 2
    Tn = complex(Tdel.subs({m: mm, g: gg, th: tt}))
    Rn = complex(Rdel.subs({m: mm, g: gg, th: tt}))
    vals.append(abs(abs(Tn) ** 2 + abs(Rn) ** 2 - 1))
report("(a) |R|^2 + |T|^2 = 1 for real theta and g", max(vals) < 1e-12)
# pole of T: 2ip = g, i.e. sinh theta = -i g/(2m).  For -2m<g<0: theta = iu, sin u = -g/(2m) in the
# physical strip of the defect, energy m cos u = sqrt(m^2 - g^2/4), the bound state of -phi''+g delta phi = -kappa^2 phi.
ok = True
for gg in (-0.3, -1.0, -1.7):
    u = math.asin(-gg / 2)
    ok &= abs(complex(sp.denom(sp.together(Tdel)).subs({m: 1, g: gg, th: sp.I * u}))) < 1e-12
    kappa = -gg / 2                       # psi ~ e^{-kappa|x|}, jump -2 kappa = g
    ok &= abs(math.cos(u) - math.sqrt(1 - kappa ** 2)) < 1e-14 and 0 < u < math.pi / 2
report("(a) bound state for -2m<g<0 at theta = iu, sin u = -g/2m, energy m cos u = sqrt(m^2-g^2/4)", ok)
ok = True
for gg in (0.5, 1.5):        # virtual states on the negative imaginary axis
    u = math.asin(gg / 2)
    s = cmath.sinh(-1j * u)
    ok &= abs(2j * s - gg) < 1e-12
for gg in (2.5, 4.0):        # g > 2m: theta = -i pi/2 +- arccosh(g/2)
    for sgn in (1, -1):
        z = -1j * math.pi / 2 + sgn * math.acosh(gg / 2)
        ok &= abs(2j * cmath.sinh(z) - gg) < 1e-12
for gg in (-2.5, -4.0):      # g < -2m: theta = +i pi/2 +- arccosh(|g|/2): tachyonic, E = m cosh theta imaginary
    for sgn in (1, -1):
        z = 1j * math.pi / 2 + sgn * math.acosh(-gg / 2)
        ok &= abs(2j * cmath.sinh(z) - gg) < 1e-12 and abs(cmath.cosh(z).real) < 1e-12
report("(a) pole regimes: virtual states for 0<g<2m, Im theta = -pi/2 for g>2m, Im theta = pi/2 (imaginary E) for g<-2m", ok)

# ---------------------------------------------------------------------------------------------
# (b), (c) symbolic rewriting in the exchange algebra of a defect
# ---------------------------------------------------------------------------------------------
# Tokens: ('A', label, rapidity_expr) and ('D',).  Numeric values decide the ordering only.
t1, t2 = sp.symbols('theta1 theta2', positive=True)
NUM = {t1: 0.9, t2: 0.4}

Sf = {}
def S(a, b, u):
    key = (a, b)
    if key not in Sf:
        Sf[key] = sp.Function(f'S_{a}{b}')
    return Sf[key](u)

def F(name, a, u):
    return sp.Function(f'{name}_{a}')(u)

def num(x):
    return float(sp.N(x.subs(NUM)))

def moves(word):
    """All applicable rewriting moves: list of (kind, position)."""
    out = []
    for i in range(len(word) - 1):
        x, y = word[i], word[i + 1]
        if x[0] == 'A' and y[0] == 'A' and num(x[2]) > num(y[2]):
            out.append(('bulk', i))
        if x[0] == 'A' and y[0] == 'D' and num(x[2]) > 0:
            out.append(('defL', i))
        if x[0] == 'D' and y[0] == 'A' and num(y[2]) < 0:
            out.append(('defR', i))
    return out

def apply(word, mv):
    kind, i = mv
    w = list(word)
    if kind == 'bulk':
        (_, a, u), (_, b, v) = w[i], w[i + 1]
        return [(S(a, b, sp.simplify(u - v)), tuple(w[:i] + [('A', b, v), ('A', a, u)] + w[i + 2:]))]
    if kind == 'defL':
        (_, a, u) = w[i]
        return [(F('R', a, u), tuple(w[:i] + [('A', a, -u), ('D',)] + w[i + 2:])),
                (F('T', a, u), tuple(w[:i] + [('D',), ('A', a, u)] + w[i + 2:]))]
    if kind == 'defR':
        (_, a, u) = w[i + 1]
        return [(F('Rb', a, -u), tuple(w[:i] + [('D',), ('A', a, -u)] + w[i + 2:])),
                (F('Tb', a, -u), tuple(w[:i] + [('A', a, u), ('D',)] + w[i + 2:]))]
    raise ValueError

def reduce_word(word, prefer):
    """Reduce to out-states, at each step taking the first move of the preferred kind if any."""
    terms = {word: sp.Integer(1)}
    done = {}
    while terms:
        new = {}
        for w, c in terms.items():
            mv = moves(w)
            if not mv:
                done[w] = done.get(w, 0) + c
                continue
            pick = [x for x in mv if x[0] in prefer] or mv
            for c2, w2 in apply(w, pick[0]):
                new[w2] = new.get(w2, 0) + c * c2
        terms = new
    return done

# (b) unitarity: A(t)D = R(t) A(-t) D + T(t) D A(t); expand A(-t)D and D A(t) again with the relations at -t.
tt = sp.Symbol('t')
lhs_R = F('R', 'a', tt) * F('R', 'a', -tt) + F('T', 'a', tt) * F('Tb', 'a', -tt)
lhs_X = F('R', 'a', tt) * F('T', 'a', -tt) + F('T', 'a', tt) * F('Rb', 'a', -tt)
# second application, by hand from the two relations evaluated at -t:
#  A(-t)D = R(-t) A(t) D + T(-t) D A(-t);   D A(t) = D A(-(-t)) = Rb(-t) D A(-t) + Tb(-t) A(t) D
coef_AD = F('R', 'a', tt) * F('R', 'a', -tt) + F('T', 'a', tt) * F('Tb', 'a', -tt)
coef_DA = F('R', 'a', tt) * F('T', 'a', -tt) + F('T', 'a', tt) * F('Rb', 'a', -tt)
report("(b) A(t)D re-expanded: coefficient of A(t)D is RR(-)+T Tb(-), of D A(-t) is R T(-)+T Rb(-)",
       sp.simplify(coef_AD - lhs_R) == 0 and sp.simplify(coef_DA - lhs_X) == 0)

# (c) the no-go reductions: A_a(t1) A_b(t2) D, theta1 > theta2 > 0
start = (('A', 'a', t1), ('A', 'b', t2), ('D',))
first_b = reduce_word(start, prefer=('defL', 'defR'))     # particle b meets the defect first
first_S = reduce_word(start, prefer=('bulk',))            # the particles scatter first
words = set(first_b) | set(first_S)
diff = {w: sp.simplify(first_b.get(w, 0) - first_S.get(w, 0)) for w in words}
w_bRaT = (('A', 'b', -t2), ('D',), ('A', 'a', t1))
w_aRbT = (('A', 'a', -t1), ('D',), ('A', 'b', t2))
w_TT = (('D',), ('A', 'b', t2), ('A', 'a', t1))
w_RR = (('A', 'a', -t1), ('A', 'b', -t2), ('D',))
report("(c) exactly four out-states", words == {w_bRaT, w_aRbT, w_TT, w_RR})
e1 = F('R', 'b', t2) * F('T', 'a', t1) * (S('a', 'b', t1 + t2) - S('a', 'b', t1 - t2))
report("(c) b reflected, a transmitted: difference = R_b T_a [S_ab(t1+t2) - S_ab(t1-t2)]",
       sp.simplify(diff[w_bRaT] - e1) == 0)
e2 = F('R', 'a', t1) * F('T', 'b', t2) * (1 - S('a', 'b', t1 - t2) * S('b', 'a', t1 + t2))
report("(c) a reflected, b transmitted: difference = R_a T_b [1 - S_ab(t1-t2) S_ba(t1+t2)]",
       sp.simplify(diff[w_aRbT] - e2) == 0)
report("(c) both transmitted: the two orders agree identically", diff[w_TT] == 0)
rr = sp.simplify(diff[w_RR].replace(Sf[('b', 'a')], Sf[('a', 'b')]))
report("(c) both reflected: difference vanishes when S_ba = S_ab (parity)", rr == 0 and diff[w_RR] != 0)
# print the individual coefficients for the lecture text
print("     b first : coefficient of A_b(-t2) D A_a(t1) =", first_b[w_bRaT])
print("     S first : coefficient of A_b(-t2) D A_a(t1) =", first_S[w_bRaT])
print("     b first : coefficient of A_a(-t1) D A_b(t2) =", first_b[w_aRbT])
print("     S first : coefficient of A_a(-t1) D A_b(t2) =", first_S[w_aRbT])
print("     b first : coefficient of A_a(-t1) A_b(-t2) D =", first_b[w_RR])
print("     S first : coefficient of A_a(-t1) A_b(-t2) D =", first_S[w_RR])
# numeric illustration: delta defect with the Lee-Yang S-matrix is inconsistent; with S = 1 it is fine
def fx(x, z):
    return (cmath.sinh(z) + 1j * math.sin(math.pi * x)) / (cmath.sinh(z) - 1j * math.sin(math.pi * x))
Tn = lambda z: complex(Tdel.subs({m: 1, g: 0.7, th: z}))
Rn = lambda z: complex(Rdel.subs({m: 1, g: 0.7, th: z}))
a1, a2 = 0.9, 0.4
bad = Rn(a2) * Tn(a1) * (fx(2 / 3, a1 + a2) - fx(2 / 3, a1 - a2))
report("(c) delta-defect amplitudes with S = f_{2/3}: the mixed coefficients disagree (|diff| > 0.1)", abs(bad) > 0.1)

# ---------------------------------------------------------------------------------------------
# (d) the triangle relation: index form = braid form; a particle is a defect
# ---------------------------------------------------------------------------------------------
eta = 0.37 + 0.21j   # generic, not a root of unity

def Rmat(u):
    """Six-vertex R(u) on C^2 (x) C^2, R12(u-v)R13(u)R23(v) = R23(v)R13(u)R12(u-v)."""
    a, b, c = cmath.sinh(u + eta), cmath.sinh(u), cmath.sinh(eta)
    return np.array([[a, 0, 0, 0], [0, b, c, 0], [0, c, b, 0], [0, 0, 0, a]], dtype=complex)

P = np.zeros((4, 4))
for i in range(2):
    for j in range(2):
        P[2 * j + i, 2 * i + j] = 1

def Rcheck(u):
    return P @ Rmat(u)

I2 = np.eye(2)
def on12(X):
    return np.kron(X, I2)
def on23(X):
    return np.kron(I2, X)

# R-check of the physical index convention: Rcheck(v_a x v_b) = sum S_ab^cd v_d x v_c
def S_index(u):
    Rc = Rcheck(u)
    Sx = np.zeros((2, 2, 2, 2), dtype=complex)
    for a, b, c, d in itertools.product(range(2), repeat=4):
        Sx[a, b, c, d] = Rc[2 * d + c, 2 * a + b]
    return Sx

def triangle_braid(Tfun, u1, u2):
    lhs = on23(Rcheck(u1 - u2)) @ on12(Tfun(u1)) @ on23(Tfun(u2))
    rhs = on12(Tfun(u2)) @ on23(Tfun(u1)) @ on12(Rcheck(u1 - u2))
    return np.max(np.abs(lhs - rhs)), lhs - rhs

# write the index form explicitly with distinct letters to avoid label clashes
def triangle_index_residual(Tfun, u1, u2):
    Sx = S_index(u1 - u2)
    def Tidx(u):
        M = Tfun(u)
        X = np.zeros((2, 2, 2, 2), dtype=complex)   # X[a, al, be, b] = T_{a al}^{be b}
        for a, al, be, b in itertools.product(range(2), repeat=4):
            X[a, al, be, b] = M[2 * be + b, 2 * a + al]
        return X
    T1, T2 = Tidx(u1), Tidx(u2)
    # LHS[a,b,al; ga,e,f] = sum_{c,d,be} S[a,b,c,d] T1[c,al,be,f] T2[d,be,ga,e]
    lhs = np.einsum('ABCD,CLMF,DMGE->ABLGEF', Sx, T1, T2)
    # RHS[a,b,al; ga,e,f] = sum_{c,d,be} T2[b,al,be,d] T1[a,be,ga,c] S[c,d,f,e]
    rhs = np.einsum('BLMD,AMGC,CDFE->ABLGEF', T2, T1, Sx)
    return lhs - rhs

def braid_as_index(D):
    """Read the braid residual (map V1 x V2 x W -> W x V2 x V1) in the index labels [a,b,al; ga,e,f]."""
    X = np.zeros((2,) * 6, dtype=complex)
    for a, b, al, ga, e, f in itertools.product(range(2), repeat=6):
        X[a, b, al, ga, e, f] = D[4 * ga + 2 * e + f, 4 * a + 2 * b + al]
    return X

u1, u2 = 0.83 + 0.1j, -0.41 + 0.05j
ok = True
for _ in range(3):
    Ms = [rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)) for _ in range(2)]
    Tr = lambda u, Ms=Ms: Ms[0] + u * Ms[1]
    _, D = triangle_braid(Tr, u1, u2)
    # scatter-first is the left side of the index form and the right side of the braid form
    ok &= np.max(np.abs(braid_as_index(D) + triangle_index_residual(Tr, u1, u2))) < 1e-12
    ok &= np.max(np.abs(braid_as_index(D))) > 1e-3
report("(d) index form and braid form have opposite residuals, entry by entry, for random T", ok)
thd = 0.29 - 0.17j
Tpart = lambda u: Rcheck(u - thd)
r1, _ = triangle_braid(Tpart, u1, u2)
report("(d) a particle is a defect: T = Rcheck(theta - theta_d) solves the triangle relation", r1 < 1e-13)
pert = np.zeros((4, 4), dtype=complex); pert[1, 1] = 0.05
r2, _ = triangle_braid(lambda u: Rcheck(u - thd) + pert, u1, u2)
report("(d) control: a perturbed T violates it", r2 > 1e-3)

# ---------------------------------------------------------------------------------------------
# (e) defects: T(theta)T(theta+i pi) = 1; shifted bulk amplitudes are defects
# ---------------------------------------------------------------------------------------------
z = sp.symbols('z')
def f_sym(x, w):
    return (sp.sinh(w) + sp.I * sp.sin(sp.pi * x)) / (sp.sinh(w) - sp.I * sp.sin(sp.pi * x))
x0, th0 = sp.Rational(2, 3), sp.Symbol('theta0')
Tw = f_sym(x0, z - th0)
Tbw = lambda w: f_sym(x0, w + th0)
def zero(expr):
    for _ in range(4):
        vals = {z: complex(rng.normal(), rng.normal()), th0: complex(rng.normal(), rng.normal())}
        if abs(complex(sp.N(expr.subs(vals)))) > 1e-10:
            return False
    return True
report("(e) standing particle S(theta-theta0): unitarity T(theta) Tbar(-theta) = 1", zero(Tw * Tbw(-z) - 1))
report("(e) crossing T(theta) = Tbar(i pi - theta)", zero(Tw - Tbw(sp.I * sp.pi - z)))
report("(e) hence T(theta) T(theta + i pi) = 1", zero(Tw * Tw.subs(z, z + sp.I * sp.pi) - 1))
report("(e) bootstrap T(theta + i pi/3) T(theta - i pi/3) = T(theta) for S = f_{2/3}",
       zero(Tw.subs(z, z + sp.I * sp.pi / 3) * Tw.subs(z, z - sp.I * sp.pi / 3) - Tw))

# ---------------------------------------------------------------------------------------------
# (f) the Lee-Yang defect
# ---------------------------------------------------------------------------------------------
def f23(w):
    return fx(2 / 3, w)
def thb(b):
    return 1j * math.pi * (3 - b) / 6
def Tb(b, w):
    return f23(w - thb(b))
def Tbb(b, w):
    return f23(w + thb(b))
ok = True
for _ in range(5):
    b = complex(rng.normal() * 3, rng.normal())
    w = complex(rng.normal(), rng.normal())
    ok &= abs(Tb(b, w) * Tbb(b, -w) - 1) < 1e-12
    ok &= abs(Tb(b, w) - Tbb(b, 1j * math.pi - w)) < 1e-12
    ok &= abs(Tb(b, w + 1j * math.pi / 3) * Tb(b, w - 1j * math.pi / 3) - Tb(b, w)) < 1e-12
report("(f) Lee-Yang defect: unitarity, crossing, bootstrap for complex b", ok)
# bound-state pole of Tbar_b at theta = i pi (b+1)/6
ok = True
for b in (-0.5, 0.3, 1.4):
    u = math.pi * (b + 1) / 6
    eps = 1e-7
    ok &= abs(Tbb(b, 1j * u + eps)) > 1e5 and 0 < u < math.pi / 2
    # bootstrap: excited defect T_b(theta) S(theta + i u) = T_{b+2}(theta)
    for w in (0.3 + 0.1j, -0.7 + 0.4j):
        ok &= abs(Tb(b, w) * f23(w + 1j * u) - Tb(b + 2, w)) < 1e-12
        ok &= abs(Tbb(b, w) * f23(w - 1j * u) - Tbb(b + 2, w)) < 1e-12
    # energies: eps(b) = m sin(b pi/6) = m cosh(theta_b); eps(b+2)-eps(b) = m cos u
    ok &= abs(math.sin(b * math.pi / 6) - cmath.cosh(thb(b))) < 1e-14
    ok &= abs(math.sin((b + 2) * math.pi / 6) - math.sin(b * math.pi / 6) - math.cos(u)) < 1e-14
report("(f) pole of Tbar_b at i pi(b+1)/6; excited defect = family at b+2; energies m sin(b pi/6) = m cosh theta_b", ok)
ok = True
for br in (3.0, -3.0, 9.0, 0.5, 2.0, 1.0):
    for bi in (0.0, 0.7):
        b = complex(br, bi)
        herm = all(abs(Tbb(b, w) - Tb(b, -w.conjugate()).conjugate()) < 1e-10
                   for w in (0.4 + 0.2j, -1.1 + 0.3j, 0.8))
        expect = abs(((br - 3) / 6) - round((br - 3) / 6)) < 1e-12
        ok &= (herm == expect)
        if expect:
            ok &= abs(abs(Tb(b, 0.77)) - 1) < 1e-12
report("(f) Hermitian analyticity (and |T|=1 at real theta) exactly when Re b = 3 mod 6", ok)

# ---------------------------------------------------------------------------------------------
# (g) Verlinde lines of M(2,5)
# ---------------------------------------------------------------------------------------------
# Minimal model M(p,p') with p=2, p'=5: Kac labels (1,s), s = 1..4 with (1,s) ~ (1,5-s); take (1,1) and (1,2).
pp, pq = 2, 5
def Smm(r, s, rho, sig):
    return (2 * math.sqrt(2 / (pp * pq)) * (-1) ** (1 + s * rho + r * sig)
            * math.sin(math.pi * pp / pq * r * rho) * math.sin(math.pi * pq / pp * s * sig))
# the standard formula is written for labels (r,s) with 1<=r<=p'-1, 1<=s<=p-1; use (r,s) = (1,1), (2,1)
labels = [(1, 1), (2, 1)]
Smat = np.array([[Smm(*a, *b) for b in labels] for a in labels])
hw = lambda r, s: ((pp * r - pq * s) ** 2 - (pp - pq) ** 2) / (4 * pp * pq)
weights = [hw(*a) for a in labels]
ok = abs(weights[0]) < 1e-14 and abs(weights[1] + 0.2) < 1e-14
ok &= np.max(np.abs(Smat @ Smat - np.eye(2))) < 1e-12
d = [Smat[1, i] / Smat[0, i] for i in range(2)]
ok &= abs(d[0] - (1 - math.sqrt(5)) / 2) < 1e-12 and abs(d[1] - (1 + math.sqrt(5)) / 2) < 1e-12
N = lambda i, j, k: sum(Smat[i, l] * Smat[j, l] * Smat[k, l] / Smat[0, l] for l in range(2))
ok &= abs(N(1, 1, 0) - 1) < 1e-12 and abs(N(1, 1, 1) - 1) < 1e-12
report("(g) M(2,5): S^2 = 1, D_Phi = (1-sqrt5)/2 on h=0 and (1+sqrt5)/2 on h=-1/5, Phi x Phi = 1 + Phi", ok)

# ---------------------------------------------------------------------------------------------
# (h) the Ising defect (exercise)
# ---------------------------------------------------------------------------------------------
def ising(gg, w):
    chi = -2 * math.atan(gg / 2)
    s, c = math.sin(chi), math.cos(chi)
    Tn = c * cmath.sinh(w) / (cmath.sinh(w) - 1j * s)
    Rn = 1j * s * cmath.cosh(w) / (cmath.sinh(w) - 1j * s)
    return Tn, Rn
ok = True
for gg in (0.3, 1.2, -0.8, 3.7):
    chi = -2 * math.atan(gg / 2)
    ok &= abs(math.sin(chi) + gg / (1 + gg ** 2 / 4)) < 1e-14 and abs(math.cos(chi) - (1 - gg ** 2 / 4) / (1 + gg ** 2 / 4)) < 1e-14
    for w in (0.4 + 0.3j, -0.9 + 0.1j):
        T1, R1 = ising(gg, w)
        T2, R2 = ising(gg, -w)
        T3, R3 = ising(gg, 1j * math.pi - w)
        ok &= abs(R1 * R2 + T1 * T2 - 1) < 1e-12 and abs(R1 * T2 + T1 * R2) < 1e-12
        ok &= abs(T1 - T3) < 1e-12 and abs(R1 + R3) < 1e-12
        T4, R4 = ising(4 / gg, w)
        ok &= abs(T4 + T1) < 1e-12 and abs(R4 - R1) < 1e-12
for gg in (2.0, -2.0):
    T1, R1 = ising(gg, 0.6)
    ok &= abs(T1) < 1e-14 and abs(abs(R1) - 1) < 1e-14
report("(h) Ising defect: unitarity, T(i pi - t) = T(t), R(i pi - t) = -R(t), g -> 4/g, T = 0 at g = +-2", ok)

# ---------------------------------------------------------------------------------------------
# (i) the twist defect
# ---------------------------------------------------------------------------------------------
alpha = 0.73
tw = np.diag([np.exp(1j * alpha), np.exp(-1j * alpha)])   # charges +-1 of the two states
# W = C: T = tw on V, and the braid relation becomes [Rcheck, tw x tw] = 0.
Rc = Rcheck(0.52 + 0.3j)
report("(i) twist defect: Rcheck commutes with e^{i alpha mu} x e^{i alpha mu}",
       np.max(np.abs(Rc @ np.kron(tw, tw) - np.kron(tw, tw) @ Rc)) < 1e-14)

# ---------------------------------------------------------------------------------------------
# (j) the delta defect on a circle (exercise): Bethe-Yang with reflection
# ---------------------------------------------------------------------------------------------
# outgoing (right-mover on the right, left-mover on the left) = M incoming, incoming = e^{ipL} outgoing,
# M = [[T, R], [R, T]] (parity).  Eigenvalues T+R = (2ip+g)/(2ip-g) and T-R = 1.  The even family must agree
# with the direct solution phi = cos(p(|x| - L/2)) on [-L/2, L/2]: 2p sin(pL/2) = g cos(pL/2).
from scipy.optimize import brentq
ok = True
gg, Lr = 0.8, 7.3
Tq = lambda pp: 2j * pp / (2j * pp - gg)
Rq = lambda pp: gg / (2j * pp - gg)
for pp in (0.37, 1.21):
    M = np.array([[Tq(pp), Rq(pp)], [Rq(pp), Tq(pp)]])
    ev = sorted(np.linalg.eigvals(M), key=lambda z: z.real)
    want = sorted([(2j * pp + gg) / (2j * pp - gg), 1.0], key=lambda z: complex(z).real)
    ok &= all(abs(a1 - a2) < 1e-12 for a1, a2 in zip(ev, want))
direct = lambda pp: 2 * pp * math.sin(pp * Lr / 2) - gg * math.cos(pp * Lr / 2)
roots = []
grid = np.linspace(0.05, 3.0, 3000)
for x0, x1 in zip(grid[:-1], grid[1:]):
    if direct(x0) * direct(x1) < 0:
        roots.append(brentq(direct, x0, x1))
for r in roots:
    ok &= abs(cmath.exp(1j * r * Lr) * (Tq(r) + Rq(r)) - 1) < 1e-9
    ok &= abs(math.tan(r * Lr / 2) - gg / (2 * r)) < 1e-6 * (1 + abs(math.tan(r * Lr / 2)))
ok &= len(roots) >= 3
report(f"(j) delta defect on a circle: eigenvalues T+R, T-R=1; the even states ({len(roots)} found) solve e^(ipL)(T+R)=1", ok)

raise SystemExit(0 if OK else 1)
