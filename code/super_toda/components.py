"""Component Lagrangians of Toda theories based on Lie superalgebras
(sec-super-toda): auxiliary-field elimination, bosonic sectors, masses.

Every check starts from the superspace Lagrangian

    L = (i/2) D_+Phi.D_-Phi + sum_{j odd} mu_j exp(alpha_j.Phi)
        + i th+ th- sum_{j even} mu_j exp(alpha_j.Phi)          (Evans-Madsen (3.2))

with the conventions of superspace.py, integrates over theta with
int d^2theta th+ th- = i, eliminates the auxiliary field sigma and compares
with the formulas displayed in the chapter.  The coupling enters as
Phi -> beta Phi, L -> L/beta^2, which only rescales mu.

Checks:
 1. Generic root data (rank 2, two odd roots and one even root, symbolic Gram
    matrix): the component Lagrangian is eq-super-components, i.e.
    Evans-Madsen (3.4).
 2. Super sinh-Gordon, C(2)^(2) with alpha_1 = -alpha_0, alpha_1^2 = 1/2:
    eq-ssg-components; the bosonic truncation is the a_1^(1) Lagrangian of
    eq-imaginary-potential with m = mu/2 and the same beta; boson and fermion
    masses are both mu = 2m.
 3. For every purely odd root system used in the chapter, the bosonic mass
    matrix is the square of the fermionic one (exr-super-masses).
 4. B(1,1)^(1): two sinh-Gordon sectors of masses 2m and m with kinetic terms
    of opposite sign (Evans-Madsen (5.4)); Penati-Zanon mass ratio m_2^2/m_1^2 = 4.
 5. D(2,1;alpha)^(1) (Evans-Madsen 1999): masses m(1+alpha), m, m alpha; kinetic
    coefficients -1/(1+alpha), 1, 1/alpha; the twisted reality condition (1)
    makes the kinetic term definite for alpha > 0; the Yukawa couplings of the
    canonically normalized fields are 4 mu sqrt(alpha(1+alpha)) psibar_i psi_j phi_k.
 6. A(1,1)^(1) (sl(2|2)^(1)): with twisted reality and gamma = 1 the bosonic
    sector is a sinh-Gordon plus a sine-Gordon theory of equal mass, with a
    positive kinetic term (Evans-Madsen (7.10)-(7.12)).
 7. The Delius-Grisaru-Penati-Zanon families A^(2)(0,2n-1) and C^(2)(n+1) for
    n = 2, 3, 4: bosonic masses and fermion masses in the ratios of their
    eqs. (2.2) and (3.2).
 8. The second supersymmetry condition (lambda_i + lambda_j) a_ij = 0,
    lambda_j^2 = -1 (Evans-Hollowood 1992 (2.5)): solutions for A(1,1)^(1),
    none for C(2)^(2), B(1,1)^(1) and D(2,1;alpha)^(1).
"""
import itertools
import os
import sys
import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from superspace import (G, Jets, sexp, berezin, I, TP, TM, dot, superspace_lagrangian,  # noqa: E402
                        eliminate_sigma, psi_lin, expected_components, components,
                        potential, mass_matrices)

FAILS = []


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok:
        FAILS.append(name)


# =========================================================== 1. generic data
g11, g12, g22 = sp.symbols('g11 g12 g22')
a1, a2, b1, b2, c1, c2 = sp.symbols('a1 a2 b1 b2 c1 c2')
m1, m2, m3 = sp.symbols('mu1 mu2 mu3')
Gram = sp.Matrix([[g11, g12], [g12, g22]])
roots = [(a1, a2), (b1, b2), (c1, c2)]
odd = [True, True, False]
mus = [m1, m2, m3]
J, Lc, sol = components(Gram, roots, odd, mus)
exp_L = expected_components(J, Gram, roots, odd, mus)
diff = (Lc - exp_L).map(lambda c: sp.simplify(sp.expand(c)))
report('generic rank-2 data: int d^2theta and sigma elimination give eq-super-components '
       '(Evans-Madsen (3.4)), including the fermion kinetic and Yukawa signs', diff.is_zero())
phis = [J.phi(0), J.phi(1)]
sig_expect = [sum(mm * rr[a] * sp.exp(dot(Gram, rr, phis)) for rr, o, mm in zip(roots, odd, mus) if o)
              for a in range(2)]
report('the auxiliary field is sigma = sum_odd mu_j alpha_j exp(alpha_j.phi) (as a vector)',
       all(sp.simplify(sol[J.sig(a)] - sig_expect[a]) == 0 for a in range(2)))

# ================================================ 2. super sinh-Gordon C(2)^(2)
beta, mu, m = sp.symbols('beta mu m', positive=True)
Gram1 = sp.Matrix([[1]])
alpha = sp.Rational(1, 2) ** sp.Rational(1, 2)          # alpha^2 = 1/2
# coupling: Phi -> beta Phi, L -> L/beta^2  <=>  roots beta*alpha, mu -> mu/beta^2
rts = [(beta * alpha,), (-beta * alpha,)]
J, Lc, _ = components(Gram1, rts, [True, True], [mu / beta**2, mu / beta**2], name='s')
Lc = Lc.map(lambda c: sp.expand(c * beta**2))          # undo the overall 1/beta^2
# field normalization: Phi = beta Phi_can, so components scale by beta; compare
# after the substitution phi -> beta phi, psi -> beta psi in the expected form
phi, dpp, dmp = J.phi(0), J.phi(0, 1, 0), J.phi(0, 0, 1)
pp, pm = J.psi('p', 0), J.psi('m', 0)
ppd, pmd = J.psi('p', 0, 0, 1), J.psi('m', 0, 1, 0)
expect = (G.scalar(beta**2 * dpp * dmp / 2
                   - mu**2 / 4 * (sp.exp(sp.sqrt(2) * beta * phi) + sp.exp(-sp.sqrt(2) * beta * phi) - 2))
          + (pp * ppd).map(lambda c: I * beta**2 / 2 * c) + (pm * pmd).map(lambda c: I * beta**2 / 2 * c)
          + (pp * pm).map(lambda c: I * mu * beta**2 * sp.cosh(beta * phi / sp.sqrt(2)) * c))
# Lc is beta^2 x (canonical Lagrangian in the rescaled fields): the fields of Lc are the
# superfield components of Phi_can, the exponent carries beta.  Compare directly:
expect_can = (G.scalar(dpp * dmp / 2
                       - mu**2 / (4 * beta**2) * (sp.exp(sp.sqrt(2) * beta * phi)
                                                  + sp.exp(-sp.sqrt(2) * beta * phi) - 2))
              + (pp * ppd).map(lambda c: I / 2 * c) + (pm * pmd).map(lambda c: I / 2 * c)
              + (pp * pm).map(lambda c: I * mu * sp.cosh(beta * phi / sp.sqrt(2)) * c))
Lcan = Lc.map(lambda c: sp.expand(c / beta**2))
ok = (Lcan - expect_can).map(lambda c: sp.simplify(sp.expand(sp.expand(c).rewrite(sp.exp)))).is_zero()
report('super sinh-Gordon: L = 1/2 d+phi d-phi + (i/2)(psi+ d-psi+ + psi- d+psi-) '
       '- (mu^2/4beta^2)(e^{sqrt2 beta phi} + e^{-sqrt2 beta phi} - 2) + i mu psi+ psi- cosh(beta phi/sqrt2) '
       '(eq-ssg-components)', ok)
U = potential(J, Lcan, 1)
V_book = (m**2 / beta**2) * (sp.exp(sp.sqrt(2) * beta * phi) + sp.exp(-sp.sqrt(2) * beta * phi) - 2)
report('bosonic truncation psi = 0 is the a_1^(1) theory of eq-imaginary-potential with m = mu/2 '
       'and the same beta', sp.simplify(U - V_book.subs(m, mu / 2)) == 0)
mB2 = sp.diff(U, phi, 2).subs(phi, 0)
report('super sinh-Gordon: boson mass^2 = mu^2 = (2m)^2, the a_1^(1) particle mass of sec-tree-level',
       sp.simplify(mB2 - mu**2) == 0)
# fermion: L_ferm = (i/2)(psi+ d-psi+ + psi- d+psi-) + i M psi+ psi-  gives d+d- psi = -M^2 psi
report('super sinh-Gordon: fermion mass = mu (Yukawa coefficient at phi = 0), equal to the boson mass',
       sp.simplify(sp.Abs(mu * sp.cosh(0)) - sp.sqrt(mB2)) == 0)


# ======================================= 3. boson mass matrix = (fermion)^2
def check_square(label, Gram_, roots_, mu_):
    grad, MB2, MF, *_ = mass_matrices(Gram_, roots_, [True] * len(roots_), mu_)
    ok = all(sp.simplify(x) == 0 for x in grad) and sp.simplify(MB2 - MF * MF) == sp.zeros(*MB2.shape)
    report(f'{label}: phi = 0 is stationary and M_B^2 = M_F^2 (exr-super-masses)', ok)
    return MB2, MF


# C(2)^(2)
check_square('C(2)^(2)', sp.Matrix([[1]]), [(alpha,), (-alpha,)], [1, 1])
# B(1,1)^(1) in the Evans-Madsen Gram matrix [[1,-1,-1],[-1,0,2],[-1,2,0]]: realize the
# roots in R^2 with an indefinite metric; alpha_1 + alpha_2 + ... relation 2 a1 + a2 + a3 = 0
GB = sp.Matrix([[1, 0], [0, -1]])
# find vectors: a1 = (1, 0) gives a1^2 = 1; a2 = (x, y), a3 = (u, v)
x, y, u, v = sp.symbols('x y u v')
eqs = [x**2 - y**2, u**2 - v**2, x - (-1), u - (-1), x * u - y * v - 2]
solB = sp.solve(eqs, [x, y, u, v], dict=True)
s0 = solB[0]
rB = [(1, 0), (s0[x], s0[y]), (s0[u], s0[v])]
cart = sp.Matrix(3, 3, lambda i, j: dot(GB, rB[i], rB[j]))
report('B(1,1)^(1): the realization reproduces the Evans-Madsen Cartan matrix and 2a1+a2+a3 = 0',
       cart == sp.Matrix([[1, -1, -1], [-1, 0, 2], [-1, 2, 0]])
       and all(sp.simplify(2 * sp.Integer(rB[0][k]) + rB[1][k] + rB[2][k]) == 0 for k in range(2)))
MB2, MF = check_square('B(1,1)^(1)', GB, rB, [2, 1, 1])
ev = sorted([sp.simplify(e) for e in MB2.eigenvals()], key=lambda q: float(q))
report('B(1,1)^(1): boson masses^2 in the ratio 1 : 4 (masses m and 2m; Penati-Zanon m_1^2 = 2, m_2^2 = 8)',
       sp.simplify(ev[1] / ev[0] - 4) == 0)
# opposite signs: the two mass eigenvectors have kinetic norms of opposite sign
vecs = MB2.eigenvects()
norms = {sp.simplify(val): sp.sign(sp.simplify((vec[0].T * GB * vec[0])[0, 0])) for val, mult, vec in vecs}
report('B(1,1)^(1): the heavy and light sinh-Gordon sectors have kinetic terms of opposite sign '
       '(Evans-Madsen (5.4))', norms[ev[0]] * norms[ev[1]] == -1)

# ================================================== 5. D(2,1;alpha)^(1)
al = sp.symbols('alpha', positive=True)
Ge = sp.diag(-(1 + al) / 2, sp.Rational(1, 2), al / 2)          # eps_i^2, Evans-Madsen (2.8)
eD = [(-1, 1, 1), (1, -1, 1), (1, 1, -1), (-1, -1, -1)]        # alpha_1, alpha_2, alpha_3, alpha_0
cartD = sp.Matrix(4, 4, lambda i, j: sp.simplify(dot(Ge, eD[i], eD[j])))
report('D(2,1;alpha): alpha_i^2 = 0, a_12 = alpha, a_13 = 1, a_23 = -1-alpha (Evans-Madsen (2.7)) '
       'and the affine Cartan matrix (2.12)',
       all(cartD[i, i] == 0 for i in range(4)) and cartD[0, 1] == al and cartD[0, 2] == 1
       and sp.simplify(cartD[1, 2] + 1 + al) == 0)
mu_ = sp.symbols('mu', positive=True)
MB2, MF = check_square('D(2,1;alpha)^(1)', Ge, eD, [mu_] * 4)
evs = sp.simplify(sp.Matrix(MB2).diagonalize()[1])
masses2 = sorted([sp.factor(evs[i, i]) for i in range(3)], key=lambda q: str(q))
target = sorted([sp.factor((2 * mu_)**2 * (1 + al)**2), sp.factor((2 * mu_)**2), sp.factor((2 * mu_)**2 * al**2)],
                key=lambda q: str(q))
report('D(2,1;alpha)^(1): boson masses m(1+alpha), m, m alpha with m = 2 mu (Evans-Madsen (3.23))',
       all(sp.simplify(p - q) == 0 for p, q in zip(masses2, target)))
fm = sorted([sp.simplify(MF[i, i]) for i in range(3)], key=lambda q: str(q))
report('D(2,1;alpha)^(1): fermion mass matrix is diagonal in the eps basis with entries '
       '-m(1+alpha)... up to sign: |M_F| = m(1+alpha), m, m alpha',
       MF.is_diagonal() and sorted([sp.simplify(q**2) for q in fm], key=str)
       == sorted([sp.simplify(q) for q in target], key=str))
# kinetic signs: standard reality -> Gram eps_i^2 = (-(1+a)/2, 1/2, a/2): indefinite.
# twisted reality (1): Phi_2 -> i Phi_2, Phi_3 -> i Phi_3 multiplies eps_2^2, eps_3^2 by -1.
tw = {1: sp.diag(1, -1, -1), 2: sp.diag(-1, 1, -1), 3: sp.diag(-1, -1, 1)}
sig = {}
for k, T in tw.items():
    Gt = sp.diag(*[Ge[i, i] * (1 if T[i, i] == 1 else -1) for i in range(3)])
    sig[k] = [sp.sign(Gt[i, i].subs(al, sp.Rational(7, 10))) for i in range(3)]
report('D(2,1;alpha)^(1), alpha > 0: only the twisted reality condition (1) gives a definite '
       '(negative) kinetic term (Evans-Madsen table after (3.28))',
       sig[1] == [-1, -1, -1] and len(set(sig[2])) > 1 and len(set(sig[3])) > 1)
# Yukawa cubic couplings.  The quadratic fermion terms are
# i mu sum_j (alpha_j.psi+)(alpha_j.psi-) exp(alpha_j.phi); expand to first order in phi.
J, Lc, _ = components(Ge, eD, [True] * 4, [mu_] * 4, name='d')
p0, m1_, ph2 = J.psi('p', 0), J.psi('m', 1), J.phi(2)
mono = list((p0 * m1_).t)[0]
coefL = sp.simplify(sp.diff(Lc.t[mono] * sp.sign(1), ph2).subs({J.phi(b): 0 for b in range(3)}))
# Evans-Madsen normalize kinetic terms as d+phi d-phi and i psi+ d-psi+, twice ours.  Our field
# with kinetic term (1/2) Gram_aa (d phi^a)^2 is phi^a = sqrt(2/|Gram_aa|) phi_EM.
gEM = sp.simplify(coefL * sp.sqrt(8) / sp.sqrt(sp.Abs(Ge[0, 0] * Ge[1, 1] * Ge[2, 2])))
# psibar_1 psi_2 = i(psi1+ psi2- + psi2+ psi1-), so g psibar_1 psi_2 phi_3 has coefficient i g
report('D(2,1;alpha)^(1): in Evans-Madsen field normalization the Yukawa coupling of '
       'psibar_1 psi_2 phi_3 has modulus 4 mu sqrt(alpha(1+alpha)) (Evans-Madsen after (3.32))',
       sp.simplify(sp.Abs(gEM) - 4 * mu_ * sp.sqrt(al * (1 + al))) == 0)

# ======================================================= 6. A(1,1)^(1), N = 2
# simple roots: four null roots on a circle; Evans-Madsen Cartan matrix (5.? Section 5 (4)):
CA = sp.Matrix([[0, 1, 0, -1], [1, 0, -1, 0], [0, -1, 0, 1], [-1, 0, 1, 0]])
GA = sp.diag(sp.Rational(1, 2), -sp.Rational(1, 2))
# phi_1 + phi_3 = phi_2 + phi_4 = 0: alpha_3 = -alpha_1, alpha_4 = -alpha_2
rA = [(1, 1), (1, -1), (-1, -1), (-1, 1)]
cartA = sp.Matrix(4, 4, lambda i, j: dot(GA, rA[i], rA[j]))
report('A(1,1)^(1): realization of the Evans-Madsen Cartan matrix with two null eigenvectors',
       cartA == CA and CA.rank() == 2)
gam = sp.symbols('gamma', positive=True)
for gval, lab in ((1, 'gamma = 1'), (gam, 'generic gamma')):
    MB2, MF = check_square(f'A(1,1)^(1), {lab}', GA, rA, [gval, 1 / gval, gval, 1 / gval])
    evA = [sp.simplify(e) for e in MB2.eigenvals()]
    report(f'A(1,1)^(1), {lab}: the two bosons have equal mass, independent of gamma '
           '(Evans-Hollowood 1992 sec. 2, Evans-Madsen (5.7))', len(set(evA)) == 1)
# twisted reality: sigma exchanges alpha_1 <-> alpha_2? Evans-Madsen: v -> i w with v along the
# negative-norm direction; with GA = diag(1/2, -1/2) the twist phi^2 -> i phi^2 makes the Gram
# matrix positive definite and turns the sinh-Gordon of the second field into sine-Gordon.
J, Lc, _ = components(GA, rA, [True] * 4, [1, 1, 1, 1], name='a')
Ub = potential(J, Lc, 2)
w = sp.symbols('w', real=True)
Utw = sp.simplify(Ub.subs(J.phi(1), I * w).rewrite(sp.cos))
p0 = J.phi(0)
report('A(1,1)^(1) with twisted reality: U = sinh-Gordon(phi^1) + sine-Gordon(w) and the kinetic '
       'term becomes positive (Evans-Madsen (7.10))',
       sp.simplify(Utw - (Utw.subs(w, 0) + Utw.subs(p0, 0) - Utw.subs({w: 0, p0: 0}))) == 0
       and sp.simplify(sp.diff(Utw, w, 2).subs({w: 0, p0: 0}) - sp.diff(Utw, p0, 2).subs({w: 0, p0: 0}) * 1) == 0
       and sp.simplify(Utw.subs(p0, 0) - Utw.subs({p0: 0, w: w + 2 * sp.pi})) == 0
       and (-GA[1, 1]) > 0)


# ===================================== 7. Delius-Grisaru-Penati-Zanon families
def dgpz_family(kind, n):
    """Root data in the normalization of Delius-Grisaru-Penati-Zanon: orthonormal
    eps_i, bosonic roots of squared length 2, fermionic roots of squared length 1."""
    E = sp.eye(n)
    e = [E[:, i] for i in range(n)]
    if kind == 'C':     # C^(2)(n+1): odd -eps_1, even eps_i - eps_{i+1}, odd eps_n
        roots_ = [-e[0]] + [e[i] - e[i + 1] for i in range(n - 1)] + [e[n - 1]]
        odd_ = [True] + [False] * (n - 1) + [True]
        mu_odd = sp.sqrt(2)          # stationarity: mu_odd^2 alpha^2 = Kac label 2 of 2 alpha
        mus_ = [mu_odd] + [2] * (n - 1) + [mu_odd]
    else:               # A^(2)(0,2n-1): even -eps_1-eps_2, even eps_i - eps_{i+1}, odd eps_n
        roots_ = [-e[0] - e[1]] + [e[i] - e[i + 1] for i in range(n - 1)] + [e[n - 1]]
        odd_ = [False] * n + [True]
        lab = [1, 1] + [2] * (n - 2)          # a_{2n-1}^(2) labels of alpha_0..alpha_{n-1}
        mus_ = lab + [sp.sqrt(2)]
    return sp.eye(n), [tuple(rr) for rr in roots_], odd_, mus_


for n in (2, 3, 4):
    for kind in ('A', 'C'):
        Gd, rd, od, md = dgpz_family(kind, n)
        grad, MB2, MF, U, J, Lc = mass_matrices(Gd, rd, od, md)
        evB = sorted(np.linalg.eigvals(np.array(sp.N(MB2), dtype=float)).real)
        evF = sorted(abs(x) ** 2 for x in np.linalg.eigvals(np.array(sp.N(MF), dtype=float)) if abs(x) > 1e-12)
        if kind == 'A':
            h = 2 * n - 1
            target_B = sorted([8 * np.sin(a * np.pi / h) ** 2 for a in range(1, n)] + [2.0])
            target_F = [2.0]
        else:
            h = 2 * n
            target_B = sorted([8 * np.sin(a * np.pi / h) ** 2 for a in range(1, n)] + [8.0])
            target_F = [2.0, 2.0]
        okB = np.allclose(evB, target_B, atol=1e-10)
        okF = np.allclose(evF, target_F, atol=1e-10)
        stat = all(abs(float(sp.N(gg))) < 1e-12 for gg in grad)
        name = 'A^(2)(0,2n-1)' if kind == 'A' else 'C^(2)(n+1)'
        report(f'{name}, n = {n}: phi = 0 stationary; boson and fermion masses^2 are those of '
               f'Delius-Grisaru-Penati-Zanon (2.2)/(3.2) with h = {h}', okB and okF and stat)


# ================================= 8. the second supersymmetry (N = 2) condition
def n2_solutions(C):
    sols = []
    for lams in itertools.product([1j, -1j], repeat=C.shape[0]):
        if all(abs((lams[i] + lams[j]) * complex(C[i, j])) < 1e-12
               for i in range(C.shape[0]) for j in range(C.shape[0])):
            sols.append(lams)
    return sols


solsA = n2_solutions(CA)
report('A(1,1)^(1): exactly two solutions lambda_j = +-i(-1)^j of (lambda_i+lambda_j)a_ij = 0 '
       '(Evans-Hollowood 1992 (2.9))',
       len(solsA) == 2 and all(all(abs(l[j] - l[0] * (-1) ** j) < 1e-12 for j in range(4)) for l in solsA))
report('no second supersymmetry for C(2)^(2), B(1,1)^(1) or D(2,1;alpha)^(1) (alpha = 0.7)',
       not n2_solutions(sp.Matrix([[sp.Rational(1, 2), -sp.Rational(1, 2)], [-sp.Rational(1, 2), sp.Rational(1, 2)]]))
       and not n2_solutions(sp.Matrix([[1, -1, -1], [-1, 0, 2], [-1, 2, 0]]))
       and not n2_solutions(cartD.subs(al, sp.Rational(7, 10))))

sys.exit(1 if FAILS else 0)
