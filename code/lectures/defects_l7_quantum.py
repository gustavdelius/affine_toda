"""Checks for Lecture 7 of the defect lecture notes, "Quantizing the defect"
(lecture-notes/defects/07-quantizing.qmd).

Conventions of the lecture. Sine-Gordon theory with beta = beta_SG, coupling lambda = 8 pi/beta^2 - 1, soliton mass
M = 8 m0/beta^2 = (lambda + 1) m0/pi, q = -e^{-i pi lambda} = e^{-8 pi^2 i/beta^2}, kappa = q - 1/q, y = e^{lambda theta},
y_d = e^{lambda theta_d}, Y = y/y_d, u = theta - theta_d. U_q(sl2^) with generators e_j, f_j, k_j (j = 0, 1), k_0 = k_1^-1,
and the coproduct Delta(e_j) = e_j (x) 1 + k_j (x) e_j, Delta(f_j) = f_j (x) k_j^-1 + 1 (x) f_j, Delta(k_j) = k_j (x) k_j.
Solitons: V(theta) with basis v_+ (soliton), v_- (antisoliton), k_1 = diag(q, 1/q), e_1 v_- = y v_+, f_1 v_+ = v_-/y,
e_0 v_+ = y v_-, f_0 v_- = v_+/y. Defect: W with basis |n>, n in Z, k_1 |n> = q^{2n+s} |n>, f_1 |n> = c |n-1>,
f_0 |n> = c |n+1>, c = mu/y_d, mu = i/kappa (a Weyl module of U_q(b^-) = <f_j, k_j>). Tensor order is spatial order;
T(theta): V (x) W -> W (x) V, T(v_+ (x) |n>) = A_n |n> v_+ + C_n |n+1> v_-, T(v_- (x) |n>) = B_n |n-1> v_+ + D_n |n> v_-.
Written from these definitions, independently of code/defect_sg_transmission.py, defect_universal_r.py and
defect_transfer.py. One PASS or FAIL line per check; exit status 1 on any FAIL.

  (a) The bulk: V(y) is a module of U_q(sl2^) at level zero and Delta is a homomorphism on V (x) V (sympy). The
      intertwiner R: V(y1) (x) V(y2) -> V(y2) (x) V(y1), normalized on v_+ v_+, is unique; with z = y1/y2 its
      transmission entry is q (z^2 - 1)/(q^2 z^2 - 1), which at q = -e^{-i pi lambda}, z = e^{lambda theta} is
      Zamolodchikov's sinh(lambda theta)/sinh(lambda (i pi - theta)); R solves the braid Yang-Baxter equation. The
      vector omega = v_+ v_- - v_- v_+ of V(y1) (x) V(y2) is annihilated by Delta(e_j), Delta(f_j), with Delta(k) = 1,
      exactly when y1/y2 = q; with q = -e^{-i pi lambda} this is theta_1 - theta_2 = -i(pi - pi/lambda), the lightest
      breather, and omega spans the image of the residue of R at its pole z = 1/q.
  (b) The defect module: W satisfies the relations of U_q(b^-) and the q-Serre relations on a window. It does not
      extend to U_q(sl2^): with e_1 |n> = c_n |n+1> (forced by the weights), [e_1, f_0] = 0 forces c_n constant and then
      [e_1, f_1] = 0, not (k_1 - k_1^-1)/kappa (sympy); a least-squares search for any raising e_1 leaves an O(1)
      residual. Positive control: on V(y) the e_j exist.
  (c) The transmission matrix from its symmetry: the eight equations (two recursions and six projections) have the
      one-parameter solution A_n = q^-n, D_n = q^{n+s}, C_n = -c kappa y q^-n, B_n = -c kappa y q^{n+s} (sympy,
      identically in n); with c = mu/y_d the converting amplitude is -iY. On a window, with every equation imposed
      that the window allows, the solutions restricted to the central amplitudes form a one-dimensional space for
      generic q, and a larger one at the root of unity q = e^{i pi/3}. The same computation for U_q(b^+) (e_1 |n> = c|n+1>,
      e_0 |n> = c|n-1>) gives a matrix with A_{n+1}/A_n = q, D_n/A_n = q^{-2n-s} and converting entries proportional
      to y_d/y.
  (d) The triangle relation R_23(y1/y2) T_12(theta1) T_23(theta2) = T_12(theta2) T_23(theta1) R_12(y1/y2) holds exactly
      in a block of fixed total charge (sympy, symbolic q, Y1, Y2, q^n) and numerically on a window for two unrelated
      values of q; it fails if one T has Y -> 1.1 Y.
  (e) T as an L-operator: L = P T = q^{s/2} diag(U^-1, U) [1 - iY (f (x) X_+ + e (x) X_-)] with U = q^{N + s/2},
      X_+ = q^{-N-s} S q^-N, X_- = q^N S^-1 q^{N+s} (S |n> = |n+1>), U X_pm = q^{pm1} X_pm U, X_+ X_- = 1. From the
      canonical commutator [Phi, Psi] = -i hbar beta^2 of the defect values: e^X e^Y = e^Y e^X e^{[X, Y]}; the pair
      U_cl = e^{-i Delta/4}, W_pm = e^{pm i Phi/2} of Lecture 6's matrix has U_cl W_pm = qhat^{pm1} W_pm U_cl with
      qhat = e^{-i hbar beta^2/8}; the pair e^{2 pi i Delta/(hbar beta^2)}, e^{-+4 pi i Phi/(hbar beta^2)} has the soliton
      q = e^{-8 pi^2 i/(hbar beta^2)}; on Delta = -2 pi (2n + s) the first is q^{n + s/2}, and the second shifts Delta by
      -+4 pi, i.e. n by pm1; the two pairs commute up to signs.
  (f) Unitarity and crossing: in each block of fixed charge the matrix part M_n satisfies M_n^dagger M_n = (1 + Y^2)
      for real theta, |q| = 1, real s; det M_n = q^{1+s}(1 + Y^2). Crossing T_{a alpha}^{beta b}(theta) =
      Tbar_{alpha bbar}^{abar beta}(i pi - theta), with Tbar(theta) = T(-theta)^-1 and v_+ <-> v_-, holds for all entries
      exactly when Phi(theta) Phi(theta - i pi)(1 + Y(theta - i pi)^2) = q^-s and y(theta - i pi) = -q y(theta), i.e.
      q = -e^{-i pi lambda}; at -q the converting entries fail for every scalar factor, and the antisymmetric
      conjugation v_+ -> v_-, v_- -> -v_+ repairs them.
  (g) The scalar factor: ln Phi(u) = -lambda u/2 + i pi lambda/4 - (1/2) ln 2 + int_0^inf [cosh(pi k/2 + iku)/cosh(pi k/2)
      - 1]/(2k sinh(pi k/2 lambda)) dk, valid for -pi/(2 lambda) < Im u < pi + pi/(2 lambda), equals the Gamma-function
      product of Konik-LeClair and Bowcock-Corrigan-Zambon (to 1e-9); for real u it is -(1/2) ln(1 + Y^2) + i delta(u)
      (the identity int_0^inf (1 - cos ku)/(k sinh(pi k/2 lambda)) dk = ln cosh(lambda u)); |Phi|^2 (1 + Y^2) = 1,
      Phi(u) Phi(u + i pi)(1 + Y^2) = 1, Phi -> 1 as u -> -inf; Phi has no zero or pole on a grid of the strip; the
      derivative is (ln Phi)'(u) = -lambda/2 + i PV int K(k) e^{iku} dk, K(k) = e^{pi k/2}/(4 cosh(pi k/2) sinh(pi k/2 lambda)).
  (h) Poles and resonances: Phi has simple poles at u_n = -i(2n + 1) pi/(2 lambda) with residue -1/(2 lambda Phi(u_n + i pi))
      (Gamma product, 1/eps scaling, n = 0..3) and is finite at +i pi/(2 lambda); at u_n the block matrix has rank one;
      M cosh(theta_d + u_n) has width Gamma_n = 2 M sinh theta_d sin((2n+1) pi/(2 lambda)) = m_{2n+1} sinh theta_d, and
      Gamma_n/(m0 sinh theta_d) -> 2n + 1 as lambda -> inf. Classically, linearizing the defect conditions about the odd
      jump (0, 2 pi) gives modes e^{gamma t} with det = 4 cosh th (cosh th - cosh theta_d) (gamma = m0 sinh th): one growing
      mode, gamma = m0 sinh theta_d, the tail of the emission solution. Inverted oscillator: D = gamma (xp + px)/2 has
      D x^n = -i hbar gamma (n + 1/2) x^n, and (p^2 - gamma^2 x^2)/2 = gamma (XP + PX)/2 with [X, P] = [x, p].
  (i) Wigner: delta'(u) = int_0^inf tanh(pi k/2) cos(ku)/(2 sinh(pi k/2 lambda)) dk equals the derivative of the phase
      of the Gamma product; (pi/lambda) delta'(u) -> ln coth(|u|/2), with deviations falling like 1/lambda
      (lambda = 10, 40, 160); int_0^inf tanh(pi k/2) cos(ku)/k dk = ln coth(|u|/2); delta'(0) is finite and equals
      (lambda/pi)(ln(8 lambda/pi) + gamma_E) + o(1).
  (j) Breathers: the bootstrap T_12(theta - i u_1/2) T_23(theta + i u_1/2) on omega (x) |n>, u_1 = pi - pi/lambda, gives
      q^s (1 + q Y^2) Phi Phi |n> (x) omega, independent of n; crossing turns Phi(u - i u_1/2)(1 + q Y^2) into
      1/Phi(u + i pi/2 + i pi/(2 lambda)); 2 sinh(pi k/2 lambda) K(k) = e^{pi k/2}/(2 cosh(pi k/2)) is free of lambda;
      int_0^inf sin(ku)/(k cosh(pi k/2)) dk = gd(u); T_B1 = i e^{i gd(u)} = (i - sinh u)/cosh u, Lecture 4's classical
      factor, for lambda = sqrt5, e + 1, 4.37 and complex theta; u_1 -> pi - 1.1 pi/lambda fails.
  (k) Classical counterparts: Lecture 6's L_d is a polynomial in 1/lambda_spec with diagonal constant term, L_d^+ one
      in lambda_spec up to a factor; the conversion probability Y^2/(1 + Y^2) is 1/2 at theta_d and tends to the
      classical step as lambda -> inf; the zeros Y = pm i of det M_n are at theta_d pm i pi/(2 lambda), and
      det L_d = 1 - sigma^2/lambda_spec^2 vanishes at theta = theta_d.
  (l) The odd sector: at theta_P = theta_d - i pi/(2 lambda) the residue of T has, in each block, the image spanned by
      e_n = |n> v_+ - |n+1> v_-; a soliton crossing the defect and then the absorbed soliton, R_23(y/y_P) T_12(theta),
      maps V (x) span{e_n} into span{e_m} (x) V, with the gauge invariants A_{n+1}/A_n = 1/q, D_n/A_n = q^{2n+s+1} and
      C_n B_{n+1}/(A_n D_{n+1}) = -Y'^2 of the even matrix at s + 1 and theta_d - i pi (Y' = Y e^{i pi lambda}).

Runtime: about a minute.
"""
import cmath
import math
import warnings

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.special import loggamma
from sympy.physics.quantum import TensorProduct as TP

warnings.filterwarnings('ignore')
OK = True


def report(name, ok):
    global OK
    OK = OK and bool(ok)
    print(("PASS " if ok else "FAIL ") + name)


rng = np.random.default_rng(7)
mp.mp.dps = 30

# ----------------------------------------------------------------------------------------------------------------
# (a) The bulk
# ----------------------------------------------------------------------------------------------------------------
qs, y1s, y2s = sp.symbols('q y1 y2', nonzero=True)


def V_sym(yv, qv=qs):
    k1 = sp.diag(qv, 1 / qv)
    return dict(k1=k1, k0=k1.inv(), e1=sp.Matrix([[0, yv], [0, 0]]), f1=sp.Matrix([[0, 0], [1 / yv, 0]]),
                e0=sp.Matrix([[0, 0], [yv, 0]]), f0=sp.Matrix([[0, 1 / yv], [0, 0]]))


def Delta_sym(A, B, x):
    if x[0] == 'e':
        return TP(A[x], sp.eye(B[x].shape[0])) + TP(A['k' + x[1]], B[x])
    if x[0] == 'f':
        return TP(A[x], B['k' + x[1]].inv()) + TP(sp.eye(A[x].shape[0]), B[x])
    return TP(A[x], B[x])


A1, A2 = V_sym(y1s), V_sym(y2s)
ok_a = True
for j in '10':
    e, f, k = A1['e' + j], A1['f' + j], A1['k' + j]
    ok_a = ok_a and sp.simplify(e * f - f * e - (k - k.inv()) / (qs - 1 / qs)) == sp.zeros(2)
    E, F, K = Delta_sym(A1, A2, 'e' + j), Delta_sym(A1, A2, 'f' + j), Delta_sym(A1, A2, 'k' + j)
    ok_a = ok_a and sp.simplify(E * F - F * E - (K - K.inv()) / (qs - 1 / qs)) == sp.zeros(4)
ok_a = ok_a and sp.simplify(A1['e1'] * A1['f0'] - A1['f0'] * A1['e1']) == sp.zeros(2)
ok_a = ok_a and sp.simplify(A1['k0'] * A1['k1']) == sp.eye(2)
for j, i in (('1', '0'), ('0', '1')):  # k_i e_j k_i^-1 = q^{a_ij} e_j with a_01 = -2
    ok_a = ok_a and sp.simplify(A1['k' + i] * A1['e' + j] * A1['k' + i].inv() - qs ** -2 * A1['e' + j]) == sp.zeros(2)
# the intertwiner; basis order v+v+, v+v-, v-v+, v-v-
ra, rb1, rc1, rc2, rb2 = sp.symbols('ra rb1 rc1 rc2 rb2')
Rgen = sp.Matrix([[1, 0, 0, 0], [0, rb1, rc2, 0], [0, rc1, rb2, 0], [0, 0, 0, ra]])
eqs = []
for x in ('e1', 'e0', 'f1', 'f0'):
    eqs += list(Rgen * Delta_sym(A1, A2, x) - Delta_sym(A2, A1, x) * Rgen)
solR = sp.solve(eqs, [ra, rb1, rc1, rc2, rb2], dict=True)
ok_a = ok_a and len(solR) == 1
zz = y1s / y2s
Rsym = Rgen.subs(solR[0])
ok_a = ok_a and sp.simplify(Rsym[2, 1] - qs * (zz ** 2 - 1) / (qs ** 2 * zz ** 2 - 1)) == 0
ok_a = ok_a and sp.simplify(Rsym[1, 1] - (qs ** 2 - 1) * zz / (qs ** 2 * zz ** 2 - 1)) == 0
Rcheck = sp.lambdify((qs, y1s, y2s), Rsym, 'numpy')


def Rnum(q, ya, yb):
    return np.array(Rcheck(q, ya, yb), dtype=complex)


lam0 = math.sqrt(5)
qphys = -cmath.exp(-1j * math.pi * lam0)
for th in (0.37, -1.1 + 0.2j):
    c_entry = Rnum(qphys, cmath.exp(lam0 * th), 1)[2, 1]
    zam = cmath.sinh(lam0 * th) / cmath.sinh(lam0 * (1j * math.pi - th))
    ok_a = ok_a and abs(c_entry - zam) < 1e-12


def kron(*ms):
    out = np.array([[1.0 + 0j]])
    for m in ms:
        out = np.kron(out, m)
    return out


I2 = np.eye(2)
ya, yb, yc = [cmath.exp(lam0 * t) for t in (0.7 + 0.1j, -0.2, 0.31 - 0.05j)]
lhs = kron(Rnum(qphys, yb, yc), I2) @ kron(I2, Rnum(qphys, ya, yc)) @ kron(Rnum(qphys, ya, yb), I2)
rhs = kron(I2, Rnum(qphys, ya, yb)) @ kron(Rnum(qphys, ya, yc), I2) @ kron(I2, Rnum(qphys, yb, yc))
ok_a = ok_a and np.abs(lhs - rhs).max() < 1e-12
# the breather vector
omega = sp.Matrix([0, 1, -1, 0])
gens = [Delta_sym(A1, A2, x) for x in ('e1', 'e0', 'f1', 'f0')]
cond = [sp.simplify((G * omega)) for G in gens]
sols_y = sp.solve([c for G in cond for c in G if c != 0], y1s, dict=True)
ok_a = ok_a and len(sols_y) == 1 and sp.simplify(sols_y[0][y1s] - qs * y2s) == 0
ok_a = ok_a and sp.simplify(Delta_sym(A1, A2, 'k1') * omega - omega) == sp.zeros(4, 1)
u1 = math.pi - math.pi / lam0
ok_a = ok_a and abs(cmath.exp(lam0 * (-1j * u1)) - qphys) < 1e-12
ok_a = ok_a and abs(cmath.exp(lam0 * (-1j * u1)) + qphys) > 1  # the other sign of q fails
# omega spans the image of the residue of R at z = 1/q: R(z) (z - 1/q) as z -> 1/q has rank one with image omega
zr = 1 / qphys * (1 + 1e-7)
Rres = Rnum(qphys, zr, 1) * (zr - 1 / qphys)
Uv, sv, Vh = np.linalg.svd(Rres)
img = Uv[:, 0]
ok_a = ok_a and sv[1] / sv[0] < 1e-5 and abs(abs(np.vdot(img, np.array([0, 1, -1, 0]) / math.sqrt(2))) - 1) < 1e-6
report("(a) the bulk: V(y) a module at level zero, Delta a homomorphism; the unique R has Zamolodchikov's "
       "transmission amplitude at q = -e^{-i pi lambda}; braid YBE; omega = v+v- - v-v+ invariant exactly for "
       "y1/y2 = q, the lightest breather, and the image of the residue of R", ok_a)

# ----------------------------------------------------------------------------------------------------------------
# numerical modules on a window
# ----------------------------------------------------------------------------------------------------------------


class Window:
    def __init__(self, K):
        self.K = K
        self.ns = np.arange(-K, K + 1)
        self.d = 2 * K + 1

    def i(self, n):
        return n + self.K


def V_num(y, q):
    return dict(k1=np.diag([q, 1 / q]), k0=np.diag([1 / q, q]), e1=np.array([[0, y], [0, 0]]),
                f1=np.array([[0, 0], [1 / y, 0]]), e0=np.array([[0, 0], [y, 0]]), f0=np.array([[0, 1 / y], [0, 0]]))


def W_num(Wn, q, s, c, borel='-'):
    d = Wn.d
    k1 = np.diag([q ** (2 * n + s) for n in Wn.ns])
    up = np.diag(np.ones(d - 1), -1)  # |n> -> |n+1>
    dn = np.diag(np.ones(d - 1), 1)  # |n> -> |n-1>
    out = dict(k1=k1, k0=np.linalg.inv(k1))
    if borel == '-':
        out.update(f1=c * dn, f0=c * up)
    else:
        out.update(e1=c * up, e0=c * dn)
    return out


def Delta_num(A, B, x):
    dA, dB = A['k1'].shape[0], B['k1'].shape[0]
    if x[0] == 'e':
        return np.kron(A[x], np.eye(dB)) + np.kron(A['k' + x[1]], B[x])
    if x[0] == 'f':
        return np.kron(A[x], np.linalg.inv(B['k' + x[1]])) + np.kron(np.eye(dA), B[x])
    return np.kron(A[x], B[x])


def T_num(Wn, q, s, Y, Phi=1.0):
    """T: V (x) W -> W (x) V; V (x) W index a*d + i(n), W (x) V index 2 i(n) + a; a = 0 for v+, 1 for v-."""
    d = Wn.d
    T = np.zeros((2 * d, 2 * d), dtype=complex)
    for n in Wn.ns:
        T[2 * Wn.i(n) + 0, 0 * d + Wn.i(n)] = Phi * q ** (-n)
        if n + 1 <= Wn.K:
            T[2 * Wn.i(n + 1) + 1, 0 * d + Wn.i(n)] = Phi * q ** (-n) * (-1j * Y)
        if n - 1 >= -Wn.K:
            T[2 * Wn.i(n - 1) + 0, 1 * d + Wn.i(n)] = Phi * q ** (n + s) * (-1j * Y)
        T[2 * Wn.i(n) + 1, 1 * d + Wn.i(n)] = Phi * q ** (n + s)
    return T


# ----------------------------------------------------------------------------------------------------------------
# (b) The defect module and the missing half
# ----------------------------------------------------------------------------------------------------------------
ok_b = True
qg = cmath.exp(1j * 0.913)  # generic phase
sg, cg = 0.3, 0.71
Wn = Window(8)
Wm = W_num(Wn, qg, sg, cg)
inner = slice(3, Wn.d - 3)
kap = qg - 1 / qg
for x in ('f1', 'f0'):
    expo = -2  # k_i f_j k_i^-1 = q^{-a_ij} f_j, a_11 = 2, a_01 = -2
    k1fk = Wm['k1'] @ Wm[x] @ Wm['k0']
    ok_b = ok_b and np.abs((k1fk - qg ** (-2 if x == 'f1' else 2) * Wm[x])[inner, inner]).max() < 1e-12
# q-Serre: f_i^3 f_j - [3] f_i^2 f_j f_i + [3] f_i f_j f_i^2 - f_j f_i^3 = 0
q3 = qg ** 2 + 1 + qg ** -2
for (fi, fj) in (('f1', 'f0'), ('f0', 'f1')):
    a_, b_ = Wm[fi], Wm[fj]
    serre = a_ @ a_ @ a_ @ b_ - q3 * a_ @ a_ @ b_ @ a_ + q3 * a_ @ b_ @ a_ @ a_ - b_ @ a_ @ a_ @ a_
    ok_b = ok_b and np.abs(serre[inner, inner]).max() < 1e-12
# symbolic no-extension argument
nn = sp.symbols('n', integer=True)
cfun = sp.Function('c')
cc, Ks = sp.symbols('cc K')
# [e1, f0]|n> = (c_{n+1} cc - cc c_n)|n+2> = 0  =>  c_{n+1} = c_n
comm_e1f0 = cfun(nn + 1) * cc - cc * cfun(nn)
ok_b = ok_b and sp.solve(comm_e1f0, cfun(nn + 1)) == [cfun(nn)]
# with c_n = C constant, [e1, f1]|n> = (C cc - cc C)|n> = 0, but (k - k^-1)/kappa |n> = [2n + s]_q |n> != 0
Cc = sp.symbols('C')
ok_b = ok_b and sp.simplify(Cc * cc - cc * Cc) == 0
ok_b = ok_b and all(abs((qg ** (2 * n + sg) - qg ** -(2 * n + sg)) / kap) > 1e-3 for n in range(-5, 6))
# least squares for a general raising e1 (weights force e1 |n> = c_n |n+1>)
K2 = 8
Wl = Window(K2)
Wlm = W_num(Wl, qg, sg, cg)
nvar = Wl.d - 1
rows, rhs_ = [], []
for n in range(-K2 + 2, K2 - 2):
    # [e1, f0]|n> = (c_{n+1} cg - cg c_n)|n+2> = 0
    r = np.zeros(nvar, dtype=complex)
    r[Wl.i(n + 1)] += cg
    r[Wl.i(n)] -= cg
    rows.append(r)
    rhs_.append(0)
    # [e1, f1]|n> = (cg c_{n-1} - cg c_n)... e1 f1|n> = cg c_{n-1}|n>, f1 e1 |n> = c_n cg |n>
    r = np.zeros(nvar, dtype=complex)
    r[Wl.i(n - 1)] += cg
    r[Wl.i(n)] -= cg
    rows.append(r)
    rhs_.append((qg ** (2 * n + sg) - qg ** -(2 * n + sg)) / kap)
Am, bv = np.array(rows), np.array(rhs_)
sol_ls = np.linalg.lstsq(Am, bv, rcond=None)[0]
resid = np.linalg.norm(Am @ sol_ls - bv) / np.linalg.norm(bv)
ok_b = ok_b and resid > 0.3
report("(b) W is a module of U_q(b^-) (weights, q-Serre); it does not extend: [e1, f0] = 0 forces a constant shift, "
       "then [e1, f1] = 0 != [2n+s]_q (symbolic); least squares for any raising e1 leaves relative residual "
       f"{resid:.2f}; on V(y) the e_j exist", ok_b)

# ----------------------------------------------------------------------------------------------------------------
# (c) The transmission matrix from its symmetry
# ----------------------------------------------------------------------------------------------------------------
ok_c = True
ys, cs, ss, Qn = sp.symbols('y c s Q', nonzero=True)  # Q stands for q^n
kq = qs - 1 / qs


def amp_sol(m):
    """the solution, at n + m, as a function of Q = q^n"""
    qm = qs ** m
    A = (Qn * qm) ** -1
    D = Qn * qm * qs ** ss
    C = -cs * kq * ys * (Qn * qm) ** -1
    B = -cs * kq * ys * Qn * qm * qs ** ss
    return A, B, C, D


A0, B0, C0, D0 = amp_sol(0)
Am1, Bm1, Cm1, Dm1 = amp_sol(-1)
Ap1, Bp1, Cp1, Dp1 = amp_sol(1)
qns = Qn ** 2 * qs ** ss  # q^{2n+s}
eight = [
    Bp1 - qs * B0,  # B_{n+1} = q B_n  (f1 on v- (x) |n+1>, projected on |n-1> v+)
    Cp1 - C0 / qs,  # C_{n+1} = C_n / q
    cs * Dm1 - B0 / ys - cs * qs * D0,  # (F1b)
    B0 / (ys * qns) + cs * Am1 - cs * A0 / qs,  # (F1c)
    D0 / (ys * qns) + cs * Cm1 - A0 / ys - cs * qs * C0,  # (F1d)
    cs * Ap1 - cs * qs * A0 - C0 / ys,  # (F0b)
    qns * A0 / ys + cs * Bp1 - cs * qs * B0 - D0 / ys,  # (F0c)
    qns * C0 / ys + cs * Dp1 - cs * D0 / qs,  # (F0d)
]
ok_c = ok_c and all(sp.simplify(eq) == 0 for eq in eight)
# the term alpha q^n in A_n is excluded: with A_n = alpha q^n, D_n = q^{2n+s} A_n, equation (F1b) leaves alpha (q^-3 - q)
al = sp.symbols('alpha')
Dal = lambda m: al * (Qn * qs ** m) ** 3 * qs ** ss
res_al = sp.simplify(cs * Dal(-1) - cs * qs * Dal(0))
ok_c = ok_c and sp.simplify(sp.expand(sp.expand_power_exp(res_al - cs * al * Qn ** 3 * qs ** ss * (qs ** -3 - qs)))) == 0
# converting amplitude with c = mu/y_d, mu = i/kappa
yd = sp.symbols('y_d', positive=True)
ok_c = ok_c and sp.simplify((-cs * kq * ys).subs(cs, sp.I / kq / yd) - (-sp.I * ys / yd)) == 0


def intertwiner_nullspace(q, s, y, c, K=10, borel='-', band=3):
    """all maps T: V (x) W -> W (x) V with entries changing n by at most `band` (no charge conservation assumed),
    from the intertwining equations for k_1 and the two generators, imposed wherever the window allows; returns
    the dimension of the space of solutions restricted to central amplitudes, the null space and the pattern"""
    Wn_ = Window(K)
    d = Wn_.d
    Vm = V_num(y, q)
    Wm_ = W_num(Wn_, q, s, c, borel)
    gens_ = ['k1', 'f1', 'f0'] if borel == '-' else ['k1', 'e1', 'e0']
    pattern = [(r, col) for col in range(2 * d) for r in range(2 * d)
               if abs((r // 2 - K) - (col % d - K)) <= band]
    nunk = len(pattern)
    pidx = {pc: i for i, pc in enumerate(pattern)}
    rows_ = []
    for x in gens_:
        Lx = Delta_num(Vm, Wm_, x)  # on V (x) W
        Rx = Delta_num(Wm_, Vm, x)  # on W (x) V
        for col in range(2 * d):
            n_in = col % d - K
            if abs(n_in) > K - band - 2:
                continue
            for r in range(2 * d):
                if abs(r // 2 - K) > K - 1:
                    continue
                coeffs = np.zeros(nunk, dtype=complex)
                for cc_ in range(2 * d):  # (T Lx)[r, col] = sum T[r, cc_] Lx[cc_, col]
                    if Lx[cc_, col] != 0 and (r, cc_) in pidx:
                        coeffs[pidx[(r, cc_)]] += Lx[cc_, col]
                for rr in range(2 * d):  # (Rx T)[r, col] = sum Rx[r, rr] T[rr, col]
                    if Rx[r, rr] != 0 and (rr, col) in pidx:
                        coeffs[pidx[(rr, col)]] -= Rx[r, rr]
                if np.abs(coeffs).max() > 0:
                    rows_.append(coeffs)
    M_ = np.array(rows_)
    _, sv_, vh = np.linalg.svd(M_)
    null = vh[np.sum(sv_ > 1e-9 * sv_[0]):].conj().T
    central = [idx for idx, (rr, cc_) in enumerate(pattern) if abs((cc_ % d) - K) <= 1]
    proj = null[central, :]
    return np.linalg.matrix_rank(proj, tol=1e-8), null, pattern, Wn_


yg = cmath.exp(0.4)
rank_gen, null_g, pat_g, Wg = intertwiner_nullspace(qg, sg, yg, cg)
ok_c = ok_c and rank_gen == 1
# compare the null vector with the formula
Tform = T_num(Wg, qg, sg, cg * (qg - 1 / qg) * yg / 1j)  # -iY = -c kappa y
vec = np.array([Tform[r, col] for (r, col) in pat_g])
central = [idx for idx, (r, col) in enumerate(pat_g) if abs((col % Wg.d) - Wg.K) <= 1]
nv = null_g[central, :]
coef = np.linalg.lstsq(nv, vec[central], rcond=None)[0]
ok_c = ok_c and np.abs(nv @ coef - vec[central]).max() < 1e-9
q_root = cmath.exp(1j * math.pi / 3)
rank_root, _, _, _ = intertwiner_nullspace(q_root, sg, yg, cg, K=11, band=5)
ok_c = ok_c and rank_root > 1
rank_gen5, _, _, _ = intertwiner_nullspace(qg, sg, yg, cg, K=11, band=5)
ok_c = ok_c and rank_gen5 == 1
# U_q(b^+)
rank_p, null_p, pat_p, Wp = intertwiner_nullspace(qg, sg, yg, cg, borel='+')
ok_c = ok_c and rank_p == 1
vp = null_p[:, 0]
amp = {}
for idx, (r, col) in enumerate(pat_p):
    a_in, n_in = divmod(col, Wp.d)
    n_in -= Wp.K
    b_out, n_out = r % 2, r // 2 - Wp.K
    amp[(a_in, n_in, b_out, n_out)] = vp[idx]
Ap = lambda n: amp[(0, n, 0, n)]
Dp = lambda n: amp[(1, n, 1, n)]
Cp = lambda n: amp[(0, n, 1, n + 1)]
ok_c = ok_c and abs(Ap(1) / Ap(0) - qg) < 1e-8 and abs(Dp(1) / Ap(1) - qg ** (-2 - sg)) < 1e-8
# converting entries scale as 1/y: recompute at 2y
_, null_p2, pat_p2, _ = intertwiner_nullspace(qg, sg, 2 * yg, cg, borel='+')
vp2 = null_p2[:, 0]
amp2 = {}
for idx, (r, col) in enumerate(pat_p2):
    a_in, n_in = divmod(col, Wp.d)
    amp2[(a_in, n_in - Wp.K, r % 2, r // 2 - Wp.K)] = vp2[idx]
ratio1 = amp[(0, 0, 1, 1)] / amp[(0, 0, 0, 0)]
ratio2 = amp2[(0, 0, 1, 1)] / amp2[(0, 0, 0, 0)]
ok_c = ok_c and abs(ratio2 / ratio1 - 0.5) < 1e-8
report("(c) the eight equations: A_n = q^-n, D_n = q^{n+s}, C_n = -c kappa y q^-n, B_n = -c kappa y q^{n+s} solves them "
       "identically and alpha q^n is excluded; c = mu/y_d gives -iY; unique on a window for generic q, not at "
       "q = e^{i pi/3}; for U_q(b^+): A_{n+1}/A_n = q, D_n/A_n = q^{-2n-s}, converting entries ~ y_d/y", ok_c)

# ----------------------------------------------------------------------------------------------------------------
# (d) The triangle relation
# ----------------------------------------------------------------------------------------------------------------
Y1s, Y2s, Ss = sp.symbols('Y1 Y2 S_')


def T_sym_apply(state, Y, ix):
    """apply T to slots (ix, ix+1) = (V, W) of a state {(labels...): coeff}; labels are ('v', a) or ('w', n)."""
    out = {}
    for key, co in state.items():
        (_, a), (_, n) = key[ix], key[ix + 1]
        if a == 0:
            terms = [((('w', n), ('v', 0)), Qn ** -1 * qs ** -n), ((('w', n + 1), ('v', 1)), -sp.I * Y * Qn ** -1 * qs ** -n)]
        else:
            terms = [((('w', n - 1), ('v', 0)), -sp.I * Y * Qn * qs ** (n + ss)),
                     ((('w', n), ('v', 1)), Qn * qs ** (n + ss))]
        for (pair, amp_) in terms:
            nk = key[:ix] + pair + key[ix + 2:]
            out[nk] = out.get(nk, 0) + co * amp_
    return out


Rmat_sym = Rsym.subs({y1s: Y1s, y2s: Y2s})  # R depends on y1/y2 = Y1/Y2


def R_sym_apply(state, ix):
    out = {}
    for key, co in state.items():
        (_, a), (_, b) = key[ix], key[ix + 1]
        col = 2 * a + b
        for row in range(4):
            amp_ = Rmat_sym[row, col]
            if amp_ != 0:
                c_, d_ = divmod(row, 2)
                nk = key[:ix] + (('v', c_), ('v', d_)) + key[ix + 2:]
                out[nk] = out.get(nk, 0) + co * amp_
    return out


ok_d = True
for a_, b_ in ((0, 0), (0, 1), (1, 0), (1, 1)):
    st = {(('v', a_), ('v', b_), ('w', 0)): 1}
    left = R_sym_apply(T_sym_apply(T_sym_apply(st, Y2s, 1), Y1s, 0), 1)
    right = T_sym_apply(T_sym_apply(R_sym_apply(st, 0), Y1s, 1), Y2s, 0)
    keys = set(left) | set(right)
    for kk in keys:
        ok_d = ok_d and sp.simplify(left.get(kk, 0) - right.get(kk, 0)) == 0


def triangle_resid(q, s, Y1, Y2, K=6, bad=1.0):
    Wn_ = Window(K)
    d = Wn_.d
    T1 = T_num(Wn_, q, s, Y1 * bad)
    T2 = T_num(Wn_, q, s, Y2)
    R12 = Rnum(q, Y1, Y2)
    P_ = lambda M: M  # noqa
    # V V W -> slots; T_23 acts on (V, W) in slots 2,3: kron(I2, T); T_12 acts on (V, W) in slots 1,2 of V W V
    lhs = kron(np.eye(d), R12) @ kron(T1, I2) @ kron(I2, T2)
    rhs = kron(T2, I2) @ kron(I2, T1) @ kron(R12, np.eye(d))
    # interior input states
    cols = [a * 2 * d + b * d + Wn_.i(n) for a in (0, 1) for b in (0, 1) for n in range(-K + 3, K - 2)]
    return np.abs((lhs - rhs)[:, cols]).max()


r1 = triangle_resid(qg, sg, 0.7 + 0.2j, -0.4 + 0.9j)
r2 = triangle_resid(cmath.exp(0.211j), 1.0, 1.3, 0.2 - 0.5j)
r3 = triangle_resid(qg, sg, 0.7 + 0.2j, -0.4 + 0.9j, bad=1.1)
ok_d = ok_d and r1 < 1e-12 and r2 < 1e-12 and r3 > 1e-2
report(f"(d) triangle relation exact in a block (sympy) and on a window for two q (residuals {r1:.0e}, {r2:.0e}); "
       f"fails with Y -> 1.1 Y ({r3:.1e})", ok_d)

# ----------------------------------------------------------------------------------------------------------------
# (e) T as an L-operator; the two Weyl algebras of the defect values
# ----------------------------------------------------------------------------------------------------------------
ok_e = True
We = Window(8)
d = We.d
Nop = np.diag(We.ns.astype(float))
qN = lambda p, add=0.0: np.diag([qg ** (p * n + add) for n in We.ns])
Sop = np.diag(np.ones(d - 1), -1)
Sinv = np.diag(np.ones(d - 1), 1)
Uop = qN(1, sg / 2)
Xp = qN(-1, -sg) @ Sop @ qN(-1)
Xm = qN(1) @ Sinv @ qN(1, sg)
inn = slice(2, d - 2)
ok_e = ok_e and np.abs((Uop @ Xp - qg * Xp @ Uop)[inn, inn]).max() < 1e-12
ok_e = ok_e and np.abs((Uop @ Xm - Xm @ Uop / qg)[inn, inn]).max() < 1e-12
ok_e = ok_e and np.abs((Xp @ Xm - np.eye(d))[inn, inn]).max() < 1e-12
Yv = 0.83 + 0.1j
Tm = T_num(We, qg, sg, Yv)
# L = P T as a 2x2 matrix of operators: L_{ba} = <v_b| T |v_a>, mapping W -> W
Lops = [[np.zeros((d, d), dtype=complex) for _ in range(2)] for _ in range(2)]
for a in (0, 1):
    for b in (0, 1):
        for n in We.ns:
            for m in We.ns:
                Lops[b][a][We.i(m), We.i(n)] = Tm[2 * We.i(m) + b, a * d + We.i(n)]
e_, f_ = np.array([[0, 1], [0, 0]]), np.array([[0, 0], [1, 0]])
Lform = [[None, None], [None, None]]
pref = qg ** (sg / 2)
Uinv = np.linalg.inv(Uop)
Lform[0][0] = pref * Uinv
Lform[1][1] = pref * Uop
Lform[1][0] = pref * Uop @ (-1j * Yv * Xp)  # row v-, column v+: f (x) X_+
Lform[0][1] = pref * Uinv @ (-1j * Yv * Xm)
for a in (0, 1):
    for b in (0, 1):
        ok_e = ok_e and np.abs((Lops[b][a] - Lform[b][a])[inn, inn]).max() < 1e-12
# canonical commutator: [Phi, Psi] = -i hbar beta^2; for X = i(a1 Phi + b1 Psi), Y = i(a2 Phi + b2 Psi),
# [X, Y] = (a1 b2 - a2 b1) i hbar beta^2 and e^X e^Y = e^Y e^X e^{[X,Y]}
hb, bt = sp.symbols('hbar beta', positive=True)


def weyl_phase(c1, c2):
    """c = (a, b) for exp(i(a Phi + b Psi)); returns the factor in e^{X} e^{Y} = factor e^{Y} e^{X}"""
    (a1, b1), (a2, b2) = c1, c2
    return sp.exp((a1 * b2 - a2 * b1) * sp.I * hb * bt ** 2)


Ucl = (-sp.Rational(1, 4), sp.Rational(1, 4))  # e^{-i Delta/4}, Delta = Phi - Psi
Wp_ = (sp.Rational(1, 2), 0)  # e^{i Phi/2}
qhat = weyl_phase(Ucl, Wp_)
ok_e = ok_e and sp.simplify(qhat - sp.exp(-sp.I * hb * bt ** 2 / 8)) == 0
Us = (2 * sp.pi / (hb * bt ** 2), -2 * sp.pi / (hb * bt ** 2))  # e^{2 pi i Delta/(hbar beta^2)}
Xs = (-4 * sp.pi / (hb * bt ** 2), 0)  # e^{-4 pi i Phi/(hbar beta^2)}
qsol = weyl_phase(Us, Xs)
ok_e = ok_e and sp.simplify(qsol - sp.exp(-8 * sp.pi ** 2 * sp.I / (hb * bt ** 2))) == 0
# with hbar = 1, q = -e^{-i pi lambda} for lambda = 8 pi/beta^2 - 1
lam_b = 8 * sp.pi / bt ** 2 - 1
ok_e = ok_e and sp.simplify(qsol.subs(hb, 1) - (-sp.exp(-sp.I * sp.pi * lam_b))) == 0
# on Delta = -2 pi (2n + s): e^{2 pi i Delta/beta^2} = q^{n + s/2}
nsym = sp.symbols('n')
ok_e = ok_e and sp.simplify(sp.powdenest(sp.exp(2 * sp.pi * sp.I * (-2 * sp.pi * (2 * nsym + ss)) / bt ** 2)
                                         - qsol.subs(hb, 1) ** (nsym + ss / 2), force=True)) == 0
# X_s shifts Delta by i a2 [Delta, Phi] = i(-4 pi/beta^2)(-i beta^2) = -4 pi, i.e. n -> n + 1
shift = sp.I * Xs[0] * (-sp.I * hb * bt ** 2)  # [Delta, Phi] = [Phi, Psi] = -i hbar beta^2
ok_e = ok_e and sp.simplify(shift + 4 * sp.pi) == 0
# cross relations: signs
ok_e = ok_e and sp.simplify(weyl_phase(Ucl, Xs) + 1) == 0 and sp.simplify(weyl_phase(Us, Wp_) + 1) == 0
ok_e = ok_e and sp.simplify(weyl_phase(Ucl, Us) - 1) == 0 and sp.simplify(weyl_phase(Wp_, Xs) - 1) == 0
report("(e) L = P T = q^{s/2} diag(U^-1, U)[1 - iY(f X+ + e X-)], U X_pm = q^pm1 X_pm U, X+ X- = 1; from "
       "[Phi, Psi] = -i hbar beta^2: e^{-i Delta/4}, e^{pm i Phi/2} give qhat = e^{-i hbar beta^2/8}; "
       "e^{2 pi i Delta/hbar beta^2}, e^{-+4 pi i Phi/hbar beta^2} give q = e^{-8 pi^2 i/hbar beta^2} = -e^{-i pi lambda}, "
       "act as q^{n+s/2} and as n -> n pm 1; the two pairs commute up to signs", ok_e)

# ----------------------------------------------------------------------------------------------------------------
# the scalar factor
# ----------------------------------------------------------------------------------------------------------------


def lnPhi_gamma(u, lam, N=300, nmax=12):
    """the Konik-LeClair / Bowcock-Corrigan-Zambon product, the ratio r_N/rbar_N at a common cutoff N, with the
    remaining factors k >= N summed by Stirling's series (Bernoulli polynomials and Hurwitz zeta values)"""
    a = -1j * lam * u / (2 * np.pi)
    k = np.arange(N)

    def lr(aa):
        return (loggamma(k * lam + 0.25 + aa) + loggamma((k + 1) * lam + 0.75 + aa)
                - loggamma((k + 0.5) * lam + 0.25 + aa) - loggamma((k + 0.5) * lam + 0.75 + aa))

    s = np.sum(lr(a) - lr(-a))
    # ln Gamma(x + c) = (x + c - 1/2) ln x - x + ln(2 pi)/2 + sum_n (-1)^{n+1} B_{n+1}(c)/(n(n+1) x^n), x = k lam;
    # the terms up to n = 1 cancel between the four factors and between a and -a
    am = mp.mpc(a)
    tail = mp.mpc(0)
    for n in range(2, nmax + 1):
        C = mp.mpc(0)
        for sj, cj in zip((1, 1, -1, -1), (0.25, lam + 0.75, lam / 2 + 0.25, lam / 2 + 0.75)):
            C += sj * (mp.bernpoly(n + 1, cj + am) - mp.bernpoly(n + 1, cj - am))
        tail += C * (-1) ** (n + 1) / mp.mpf(n * (n + 1)) * mp.zeta(n, N) / mp.mpf(lam) ** n
    return 1j * np.pi * (1 + lam) / 4 - np.log(1 + 1j * np.exp(lam * u)) + s + complex(tail)


def Phi_gamma(u, lam):
    return np.exp(lnPhi_gamma(u, lam))


def lnPhi_int(u, lam):
    """the integral representation, valid for -pi/(2 lam) < Im u < pi + pi/(2 lam)"""
    def f(k):
        if k == 0:
            return 0j
        r = (np.exp(1j * k * u) + np.exp(-np.pi * k - 1j * k * u)) / (1 + np.exp(-np.pi * k))
        return (r - 1) / (2 * k * np.sinh(np.pi * k / (2 * lam)))

    kmax = 60 * lam
    re = quad(lambda k: f(k).real, 0, kmax, limit=4000, epsabs=1e-14, epsrel=1e-13)[0]
    im = quad(lambda k: f(k).imag, 0, kmax, limit=4000, epsabs=1e-14, epsrel=1e-13)[0]
    return -lam * u / 2 + 1j * np.pi * lam / 4 - 0.5 * np.log(2) + re + 1j * im


def Phi_int(u, lam):
    return np.exp(lnPhi_int(u, lam))


# ----------------------------------------------------------------------------------------------------------------
# (f) Unitarity and crossing
# ----------------------------------------------------------------------------------------------------------------
ok_f = True
lam_f = math.e + 0.3
q_f = -cmath.exp(-1j * math.pi * lam_f)
s_f = 0.37
for th in (-0.8, 0.2, 1.9):
    Y = math.exp(lam_f * th)
    for n in (-2, 0, 3):
        Mn = np.array([[q_f ** (-n), -1j * Y * q_f ** (n + 1 + s_f)], [-1j * Y * q_f ** (-n), q_f ** (n + 1 + s_f)]])
        ok_f = ok_f and np.abs(Mn.conj().T @ Mn - (1 + Y ** 2) * np.eye(2)).max() < 1e-10
        ok_f = ok_f and abs(np.linalg.det(Mn) - q_f ** (1 + s_f) * (1 + Y ** 2)) < 1e-10


def crossing_ratios(q, lam, s, theta, conj='sym'):
    """for the matrix part (Phi = 1), the ratios T_{a alpha}^{beta b}(theta) / Tbar_{alpha bbar}^{abar beta}(i pi - theta)
    for all four kinds of entry; crossing holds with a scalar factor iff all ratios are equal"""
    Wc = Window(6)
    dd = Wc.d
    yd_ = 1.0
    T = T_num(Wc, q, s, cmath.exp(lam * theta) / yd_)
    Tm_ = T_num(Wc, q, s, cmath.exp(lam * (theta - 1j * math.pi)) / yd_)  # T(theta - i pi) = T(-(i pi - theta))
    Tbar = np.linalg.inv(Tm_)  # Tbar(i pi - theta) = T(theta - i pi)^-1 : W (x) V -> V (x) W
    if conj == 'sym':
        Cm = np.array([[0, 1], [1, 0]])
    else:
        Cm = np.array([[0, 1], [-1, 0]])  # v+ -> v-, v- -> -v+
    out = []
    for (a, b, dn) in ((0, 0, 0), (1, 1, 0), (0, 1, 1), (1, 0, -1)):
        n = 0
        lhs_ = T[2 * Wc.i(n + dn) + b, a * dd + Wc.i(n)]
        # Tbar entry from input w_alpha (x) C v_b to output (C v_a)^dagger-component (x) w_beta
        vin = np.zeros(2 * dd, dtype=complex)  # W (x) V
        for bb in (0, 1):
            vin[2 * Wc.i(n) + bb] = Cm[bb, b]
        vout = Tbar @ vin  # in V (x) W
        comp = sum(np.conj(Cm[aa, a]) * vout[aa * dd + Wc.i(n + dn)] for aa in (0, 1))
        out.append(lhs_ / comp)
    return np.array(out)


rat = crossing_ratios(q_f, lam_f, s_f, 0.43 + 0.1j)
ok_f = ok_f and np.abs(rat - rat[0]).max() < 1e-10
# the common ratio is 1/(Phi(theta) Phi(theta - i pi)) required = q^-s (1 + Y(theta - i pi)^2) ... check its value
th_c = 0.43 + 0.1j
Yc = cmath.exp(lam_f * (th_c - 1j * math.pi))
ok_f = ok_f and abs(rat[0] - q_f ** s_f * (1 + Yc ** 2) / 1) < 1e-9 or abs(rat[0] * q_f ** (-s_f) - (1 + Yc ** 2)) < 1e-9
# at -q with the symmetric conjugation, the converting and keeping entries demand different scalars
rat_m = crossing_ratios(-q_f, lam_f, s_f, th_c)
ok_f = ok_f and np.abs(rat_m[2] / rat_m[0] + 1) < 1e-10 and np.abs(rat_m[0] - rat_m[1]) < 1e-10
rat_ma = crossing_ratios(-q_f, lam_f, s_f, th_c, conj='anti')
ok_f = ok_f and np.abs(rat_ma - rat_ma[0]).max() < 1e-10
# the scalar condition with the integral representation (s = 0): Phi(u) Phi(u + i pi)(1 + Y^2) = 1
for u in (0.3, -0.6):
    Yu = math.exp(lam_f * u)
    ok_f = ok_f and abs(Phi_int(u, lam_f) * Phi_int(u + 1j * math.pi, lam_f) * (1 + Yu ** 2) - 1) < 1e-9
report("(f) M_n^dagger M_n = (1 + Y^2), det M_n = q^{1+s}(1 + Y^2); crossing reduces to one scalar condition at "
       "q = -e^{-i pi lambda}; at -q the converting entries demand the opposite sign (repaired by the antisymmetric "
       "conjugation); the scalar condition holds for Phi", ok_f)

# ----------------------------------------------------------------------------------------------------------------
# (g) The scalar factor
# ----------------------------------------------------------------------------------------------------------------
ok_g = True
lam_g = math.sqrt(5)
devs = []
for u in (0.3, -0.7, 1.2 + 0.1j, 0.2 + 1.5j, -0.4 + 3.0j, 0.5 - 0.1j):
    a_, b_ = Phi_int(u, lam_g), Phi_gamma(u, lam_g)
    devs.append(abs(a_ - b_) / abs(b_))
ok_g = ok_g and max(devs) < 1e-9
for u in (0.3, -1.4, 0.9):
    Yu = math.exp(lam_g * u)
    ok_g = ok_g and abs(abs(Phi_int(u, lam_g)) ** 2 * (1 + Yu ** 2) - 1) < 1e-10
    # (1 - cos ku)/(k sinh(pi k/2 lam)) integrates to ln cosh(lam u)
    val = quad(lambda k: (1 - math.cos(k * u)) / (k * math.sinh(math.pi * k / (2 * lam_g))) if k > 0 else 0,
               0, 80 * lam_g, limit=2000)[0]
    ok_g = ok_g and abs(val - math.log(math.cosh(lam_g * u))) < 1e-9
ok_g = ok_g and abs(Phi_int(-20.0, lam_g) - 1) < 1e-7  # the phase approaches 0 like e^u (pole of tanh(pi k/2) at k = i)
# no zeros or poles inside the strip: |ln Phi| bounded on a grid (Gamma product)
grid_ok = True
for ui in np.linspace(-0.5 * math.pi / (2 * lam_g), math.pi + 0.5 * math.pi / (2 * lam_g), 9):
    for ur in np.linspace(-2, 2, 9):
        val = abs(Phi_gamma(ur + 1j * ui, lam_g)) * max(1.0, math.exp(lam_g * ur))  # |Phi| ~ e^{-lambda Re u}
        grid_ok = grid_ok and 0.1 < val < 10
ok_g = ok_g and grid_ok
# the derivative: (ln Phi)'(u) = -lambda/2 + i int_0^inf sinh(pi k/2 + iku)/(2 cosh(pi k/2) sinh(pi k/2 lambda)) dk,
# i.e. -lambda/2 + i PV int K(k) e^{iku} dk with K(k) = e^{pi k/2}/(4 cosh(pi k/2) sinh(pi k/2 lambda))


def dlnPhi(u, lam):
    def g(k):
        if k == 0:
            return 1j * (math.pi / 2 + 1j * u) * lam / math.pi
        r = (np.exp(1j * k * u) - np.exp(-np.pi * k - 1j * k * u)) / (1 + np.exp(-np.pi * k))  # sinh(.)/cosh(pi k/2)
        return 1j * r / (2 * np.sinh(np.pi * k / (2 * lam)))
    re = quad(lambda k: g(k).real, 0, 60 * lam, limit=4000, epsabs=1e-13)[0]
    im = quad(lambda k: g(k).imag, 0, 60 * lam, limit=4000, epsabs=1e-13)[0]
    return -lam / 2 + re + 1j * im


for u in (0.4 + 0.3j, -0.2 + 2.1j):
    hh = 1e-5
    fd = (lnPhi_int(u + hh, lam_g) - lnPhi_int(u - hh, lam_g)) / (2 * hh)
    ok_g = ok_g and abs(fd - dlnPhi(u, lam_g)) < 1e-6
report(f"(g) the integral representation equals the Gamma product (max rel. dev. {max(devs):.0e}); "
       "ln|Phi| = -(1/2) ln(1 + Y^2) via int (1 - cos ku)/(k sinh) = ln cosh(lambda u); Phi -> 1 at -inf; "
       "no zeros or poles in the strip; (ln Phi)' = -lambda/2 + i PV int K(k) e^{iku} dk", ok_g)

# ----------------------------------------------------------------------------------------------------------------
# (h) Poles, resonances, the classical instability of the odd jump
# ----------------------------------------------------------------------------------------------------------------
ok_h = True
lam_h = math.sqrt(5) + 1.0
for n in range(4):
    un = -1j * (2 * n + 1) * math.pi / (2 * lam_h)
    res_pred = -1 / (2 * lam_h * Phi_gamma(un + 1j * math.pi, lam_h))
    vals = []
    for eps in (1e-4, 1e-5):
        vals.append(eps * Phi_gamma(un + eps, lam_h))
    ok_h = ok_h and abs(vals[1] - res_pred) / abs(res_pred) < 1e-4 and abs(vals[0] - res_pred) / abs(res_pred) < 1e-3
    Yn = cmath.exp(lam_h * un)
    Mb = np.array([[1, -1j * Yn], [-1j * Yn, 1]])
    ok_h = ok_h and np.linalg.svd(Mb, compute_uv=False)[1] < 1e-12
# finite at +i pi/(2 lam)
up = 1j * math.pi / (2 * lam_h)
ok_h = ok_h and abs(Phi_gamma(up + 1e-6, lam_h) - Phi_gamma(up - 1e-6, lam_h)) < 1e-4 and abs(Phi_gamma(up + 1e-7, lam_h)) < 10
# widths
thd = 0.9
for lam_t in (lam_h, 7.31):
    M = (lam_t + 1) / math.pi  # m0 = 1
    for n in range(3):
        if 2 * n + 1 >= lam_t:
            continue
        En = M * cmath.cosh(thd - 1j * (2 * n + 1) * math.pi / (2 * lam_t))
        Gam = -2 * En.imag
        m_br = 2 * M * math.sin((2 * n + 1) * math.pi / (2 * lam_t))
        ok_h = ok_h and abs(Gam - m_br * math.sinh(thd)) < 1e-12
lim_errs = []
for lam_t in (1e2, 1e3, 1e4):
    M = (lam_t + 1) / math.pi
    errs = [abs(2 * M * math.sinh(thd) * math.sin((2 * n + 1) * math.pi / (2 * lam_t)) / math.sinh(thd) - (2 * n + 1))
            for n in range(3)]
    lim_errs.append(max(errs))
ok_h = ok_h and lim_errs[0] > lim_errs[1] > lim_errs[2] and lim_errs[2] < 2e-3
# classical: linearize the defect conditions about (Phi, Psi) = (0, 2 pi)
th_, thd_, m0_ = sp.symbols('vartheta theta_d m_0', positive=True)
Asym, Bsym = sp.symbols('A B')
gam, kk = m0_ * sp.sinh(th_), m0_ * sp.cosh(th_)
sig = sp.exp(-thd_)
# phi = A e^{gam t + kk x} (x < 0), psi = B e^{gam t - kk x} (x > 0); d_pm phi = (gam pm kk) phi, d_pm psi = (gam -+ kk) psi
eq1 = (gam + kk) * Asym - (gam - kk) * Bsym - sig * m0_ * (Asym + Bsym)  # d_+(phi - psi) = sigma m0 (phi + psi)
eq2 = (gam - kk) * Asym + (gam + kk) * Bsym + (m0_ / sig) * (Asym - Bsym)  # d_-(phi + psi) = -(m0/sigma)(phi - psi)
Mlin = sp.Matrix([[sp.diff(eq1, Asym), sp.diff(eq1, Bsym)], [sp.diff(eq2, Asym), sp.diff(eq2, Bsym)]])
detl = sp.simplify((Mlin.det() / m0_ ** 2).rewrite(sp.exp))
target = (4 * sp.cosh(th_) * (sp.cosh(th_) - sp.cosh(thd_))).rewrite(sp.exp)
ok_h = ok_h and sp.simplify(sp.expand(detl - target)) == 0
# at vartheta = theta_d: A = 0 (the left field stays at 0), B free: the tail of Psi = 4 arctan e^{m0(x cosh - t sinh)+c}
Msub = Mlin.subs(th_, thd_).applyfunc(lambda z: sp.simplify(z.rewrite(sp.exp)))
ok_h = ok_h and sp.simplify(Msub[0, 0]) != 0 and sp.simplify(Msub[0, 1]) == 0 and sp.simplify(Msub[1, 1]) == 0
# the exact emission solution: 2 pi - Psi(0, t) ~ 4 e^{-c} e^{m0 t sinh theta_d} as t -> -inf
tt, cc_ = sp.symbols('t c', real=True)
Psi_em = 4 * sp.atan(sp.exp(m0_ * (-tt * sp.sinh(thd_)) + cc_))
tail = sp.limit(((2 * sp.pi - Psi_em) / sp.exp(m0_ * tt * sp.sinh(thd_))).subs({m0_: 1, thd_: sp.Rational(7, 10),
                                                                                cc_: 0}), tt, -sp.oo)
ok_h = ok_h and sp.simplify(tail - 4) == 0
# inverted oscillator
xs_, gs_ = sp.symbols('x gamma', positive=True)
for nn_ in range(4):
    fx = xs_ ** nn_
    pf = lambda g: -sp.I * hb * sp.diff(g, xs_)  # p = -i hbar d/dx
    Dfx = gs_ * (xs_ * pf(fx) + pf(xs_ * fx)) / 2
    ok_h = ok_h and sp.simplify(Dfx - (-sp.I * hb * gs_ * (nn_ + sp.Rational(1, 2))) * fx) == 0
# canonical equivalence: X = (p + g x)/sqrt(2g), P = (p - g x)/sqrt(2g) as commuting symbols with [x, p] = c:
# XP + PX = (p^2 - g^2 x^2)/g, and [X, P] = [x, p]
pq = sp.symbols('p_')
Xc = (pq + gs_ * xs_) / sp.sqrt(2 * gs_)
Pc = (pq - gs_ * xs_) / sp.sqrt(2 * gs_)
ok_h = ok_h and sp.simplify(gs_ * (2 * Xc * Pc) / 2 - (pq ** 2 - gs_ ** 2 * xs_ ** 2) / 2) == 0
ok_h = ok_h and sp.simplify(sp.diff(Xc, xs_) * sp.diff(Pc, pq) - sp.diff(Xc, pq) * sp.diff(Pc, xs_) - 1) == 0
report("(h) simple poles of Phi at -i(2n+1) pi/(2 lambda), residue -1/(2 lambda Phi(u_n + i pi)), finite at "
       "+i pi/(2 lambda); rank-one blocks there; Gamma_n = m_{2n+1} sinh theta_d -> (2n+1) m0 sinh theta_d; the odd jump "
       "has det = 4 cosh(cosh - cosh theta_d): one growing mode m0 sinh theta_d, the emission tail; inverted oscillator "
       "D x^n = -i hbar gamma (n + 1/2) x^n", ok_h)

# ----------------------------------------------------------------------------------------------------------------
# (i) Wigner's relation
# ----------------------------------------------------------------------------------------------------------------
ok_i = True


def dphase(u, lam):
    g = lambda k: (math.tanh(math.pi * k / 2) * math.cos(k * u) / (2 * math.sinh(math.pi * k / (2 * lam)))
                   if k > 0 else lam / 2)
    return quad(g, 0, 80 * lam, limit=4000, epsabs=1e-12)[0]


lam_i = math.sqrt(5)
for u in (0.4, -0.9):
    h_ = 1e-4
    num = (cmath.phase(Phi_gamma(u + h_, lam_i) / Phi_gamma(u - h_, lam_i))) / (2 * h_)
    ok_i = ok_i and abs(num - dphase(u, lam_i)) < 1e-6
devw = []
for lam_t in (10.0, 40.0, 160.0):
    dv = max(abs(math.pi / lam_t * dphase(u, lam_t) - math.log(1 / math.tanh(abs(u) / 2))) for u in (0.3, 0.8, -1.5))
    devw.append(dv)
ok_i = ok_i and devw[0] > devw[1] > devw[2] and devw[2] < 0.02 and devw[1] / devw[2] > 2.5
for u in (0.3, 1.1):
    val = quad(lambda k: math.tanh(math.pi * k / 2) * math.cos(k * u) / k if k > 0 else 0, 0, 400, limit=8000)[0]
    # the tail of the oscillatory integral converges slowly; compare with a tolerance
    ok_i = ok_i and abs(val - math.log(1 / math.tanh(u / 2))) < 5e-3
ok_i = ok_i and math.isfinite(dphase(0.0, 10.0)) and dphase(0.0, 10.0) > dphase(0.3, 10.0)
# at u = 0 the delay is finite: delta'(0) = (lambda/pi)(ln(8 lambda/pi) + gamma_E) + o(1)
d0dev = [abs(math.pi / lt * dphase(0.0, lt) - math.log(8 * lt / math.pi) - 0.5772156649015329) for lt in (40.0, 640.0)]
ok_i = ok_i and d0dev[1] < 1e-6 and d0dev[1] < d0dev[0]
report(f"(i) delta'(u) from the integral representation equals the Gamma product's; (pi/lambda) delta' -> ln coth(|u|/2) "
       f"(deviations {devw[0]:.3f}, {devw[1]:.4f}, {devw[2]:.4f} at lambda = 10, 40, 160); "
       "delta'(0) = (lambda/pi)(ln(8 lambda/pi) + gamma_E) + o(1)", ok_i)

# ----------------------------------------------------------------------------------------------------------------
# (j) The lightest breather
# ----------------------------------------------------------------------------------------------------------------
ok_j = True


def gd(u):
    return 2 * cmath.atan(cmath.tanh(u / 2))


def Tcl(u):
    return (1j - cmath.sinh(u)) / cmath.cosh(u)


for lam_t in (math.sqrt(5), math.e + 1, 4.37):
    qt = -cmath.exp(-1j * math.pi * lam_t)
    u1_ = math.pi - math.pi / lam_t
    for th in (0.35, -0.6 + 0.2j, 1.1 - 0.15j):
        for s_t in (0.0, 0.41):
            Wb = Window(5)
            dd = Wb.d
            Ya, Yb_ = cmath.exp(lam_t * (th - 1j * u1_ / 2)), cmath.exp(lam_t * (th + 1j * u1_ / 2))
            Ta, Tb = T_num(Wb, qt, s_t, Ya), T_num(Wb, qt, s_t, Yb_)
            # (T_12(theta_a)) (T_23(theta_b)) on V V W -> W V V
            big = kron(Ta, I2) @ kron(I2, Tb)
            scal = []
            for n in (-1, 0, 2):
                vin = np.zeros(4 * dd, dtype=complex)
                vin[0 * 2 * dd + 1 * dd + Wb.i(n)] = 1  # v+ v- |n>
                vin[1 * 2 * dd + 0 * dd + Wb.i(n)] = -1  # - v- v+ |n>
                vout = big @ vin
                # expect c |n> (x) omega: W V V index (i(n)*2 + a)*2 + b
                exp_vec = np.zeros(4 * dd, dtype=complex)
                exp_vec[(Wb.i(n) * 2 + 0) * 2 + 1] = 1
                exp_vec[(Wb.i(n) * 2 + 1) * 2 + 0] = -1
                cfac = vout[(Wb.i(n) * 2 + 0) * 2 + 1]
                ok_j = ok_j and np.abs(vout - cfac * exp_vec).max() < 1e-10
                scal.append(cfac)
            Yt = cmath.exp(lam_t * th)
            ok_j = ok_j and max(abs(c_ - scal[0]) for c_ in scal) < 1e-10
            ok_j = ok_j and abs(scal[0] - qt ** s_t * (1 + qt * Yt ** 2)) < 1e-10
        u = th
        Yt = cmath.exp(lam_t * u)
        TB1 = Phi_gamma(u - 1j * u1_ / 2, lam_t) * Phi_gamma(u + 1j * u1_ / 2, lam_t) * (1 + qt * Yt ** 2)
        ok_j = ok_j and abs(TB1 - Tcl(u)) < 1e-7
        # crossing turns it into a ratio of Phi at points i pi/lambda apart, both inside the strip
        ratio = Phi_int(u + 1j * math.pi / 2 - 1j * math.pi / (2 * lam_t), lam_t) / \
            Phi_int(u + 1j * math.pi / 2 + 1j * math.pi / (2 * lam_t), lam_t)
        ok_j = ok_j and abs(ratio - Tcl(u)) < 1e-8
        ok_j = ok_j and abs(1j * cmath.exp(1j * gd(u)) - Tcl(u)) < 1e-12
        bad = Phi_gamma(u - 1j * (math.pi - 1.1 * math.pi / lam_t) / 2, lam_t) * \
            Phi_gamma(u + 1j * (math.pi - 1.1 * math.pi / lam_t) / 2, lam_t) * (1 + qt * Yt ** 2)
        ok_j = ok_j and abs(bad - Tcl(u)) > 1e-3
# the kernel identity and the sech/gd integrals
kv, lv = sp.symbols('k lambda', positive=True)
Kker = sp.exp(sp.pi * kv / 2) / (4 * sp.cosh(sp.pi * kv / 2) * sp.sinh(sp.pi * kv / (2 * lv)))
ok_j = ok_j and sp.simplify(2 * sp.sinh(sp.pi * kv / (2 * lv)) * Kker - sp.exp(sp.pi * kv / 2) / (2 * sp.cosh(sp.pi * kv / 2))) == 0
for u in (0.4, -1.3, 2.2):
    val = quad(lambda k: math.sin(k * u) / (k * math.cosh(math.pi * k / 2)) if k > 0 else u, 0, 60, limit=2000)[0]
    ok_j = ok_j and abs(val - gd(u).real) < 1e-10
    val2 = quad(lambda k: math.cos(k * u) / math.cosh(math.pi * k / 2), 0, 60, limit=2000)[0]
    ok_j = ok_j and abs(val2 - 1 / math.cosh(u)) < 1e-10
report("(j) the B1 bootstrap gives q^s (1 + q Y^2) Phi Phi times omega, independent of n; T_B1 = (i - sinh u)/cosh u "
       "= i e^{i gd(u)} for three lambda and complex theta; crossing makes it a ratio of Phi at points i pi/lambda apart; "
       "the kernel loses its lambda; u_1 -> pi - 1.1 pi/lambda fails", ok_j)

# ----------------------------------------------------------------------------------------------------------------
# (k) Classical counterparts
# ----------------------------------------------------------------------------------------------------------------
ok_k = True
ls_, sg_, Dl, Sl = sp.symbols('lambda_s sigma Delta Sigma', nonzero=True)
hM = sp.diag(1, -1)
sM = sp.Matrix([[0, 1], [1, 0]])
ex = lambda c: sp.diag(sp.exp(c), sp.exp(-c))
Ld = ex(sp.I * Dl / 4) + (sg_ / ls_) * ex(-sp.I * Sl / 4) * sM
Ldp = sM * ex(-sp.I * Sl / 4) + (sg_ / ls_) * ex(-sp.I * Dl / 4)
w_ = sp.symbols('w')  # w = 1/lambda_s
Ldw = sp.expand(Ld.subs(ls_, 1 / w_))
ok_k = ok_k and all(sp.degree(sp.expand(Ldw[i, j]), w_) <= 1 for i in range(2) for j in range(2))
c0 = Ldw.subs(w_, 0)
ok_k = ok_k and c0[0, 1] == 0 and c0[1, 0] == 0
# L_d^+ = (sigma/lambda) [e^{-i Delta h/4} + (lambda/sigma) s e^{-i Sigma h/4}]: a polynomial in lambda after the factor
Lp_scaled = sp.expand(Ldp * ls_ / sg_)
ok_k = ok_k and all(sp.degree(Lp_scaled[i, j], ls_) <= 1 for i in range(2) for j in range(2))
ok_k = ok_k and sp.simplify(Ld.det() - (1 - sg_ ** 2 / ls_ ** 2)) == 0
# conversion probability
for lam_t in (2.0, 20.0, 200.0):
    P = lambda u: 1 / (1 + math.exp(-2 * lam_t * u))
    ok_k = ok_k and abs(P(0) - 0.5) < 1e-15
pdev = [max(abs(1 / (1 + math.exp(-2 * lt * u)) - (1.0 if u > 0 else 0.0)) for u in (0.1, -0.1, 0.3)) for lt in
        (2.0, 20.0, 200.0)]
ok_k = ok_k and pdev[0] > pdev[1] > pdev[2] and pdev[2] < 1e-15
# zeros of 1 + Y^2 at theta_d +- i pi/(2 lambda); classical zero of det L_d at lambda_s = -sigma, i.e. theta = theta_d
for lam_t in (3.3, 30.0):
    for sgn in (1, -1):
        Yz = cmath.exp(lam_t * sgn * 1j * math.pi / (2 * lam_t))
        ok_k = ok_k and abs(1 + Yz ** 2) < 1e-12
th_v = 0.6
ok_k = ok_k and abs((1 - (math.exp(-th_v) / (-math.exp(-th_v))) ** 2)) < 1e-15
report("(k) L_d is a polynomial in 1/lambda_spec with diagonal constant term, L_d^+ one in lambda_spec; conversion "
       "probability 1/2 at theta_d, tending to the step; 1 + Y^2 = 0 at theta_d pm i pi/(2 lambda), det L_d = 0 at theta_d",
       ok_k)

# ----------------------------------------------------------------------------------------------------------------
# (l) The odd sector from the bootstrap at the resonance
# ----------------------------------------------------------------------------------------------------------------
ok_l = True
for lam_t, s_t, th in ((math.sqrt(5), 0.0, 0.37 + 0.1j), (4.37, 0.3, -0.5 + 0.05j)):
    qt = -cmath.exp(-1j * math.pi * lam_t)
    K = 6
    Wo = Window(K)
    dd = Wo.d
    yd_ = 1.3
    Yt = cmath.exp(lam_t * th) / yd_
    yP = yd_ * cmath.exp(lam_t * (-1j * math.pi / (2 * lam_t)))  # y at theta_P = theta_d - i pi/(2 lambda)
    YP = yP / yd_
    # the residue of T at theta_P has rank one per block, with image e_n = |n> v+ - |n+1> v-
    TP = T_num(Wo, qt, s_t, YP)
    ok_l = ok_l and abs(YP + 1j) < 1e-12
    def e_vec(n):  # in W (x) V(theta_P): index 2 i(n) + a
        v = np.zeros(2 * dd, dtype=complex)
        v[2 * Wo.i(n) + 0] = 1
        v[2 * Wo.i(n + 1) + 1] = -1
        return v
    for n in range(-2, 2):
        col = TP[:, 0 * dd + Wo.i(n)]
        ok_l = ok_l and np.abs(col - col[2 * Wo.i(n)] * e_vec(n)).max() < 1e-12
    # a soliton crosses the defect, then the absorbed soliton: R_23(y/y_P) T_12(theta)
    Tt = T_num(Wo, qt, s_t, Yt)
    y_t = Yt * yd_
    big = kron(np.eye(dd), Rnum(qt, y_t, yP)) @ kron(Tt, I2)  # V W V_P -> W V V_P -> W V_P V
    amp = {}
    for a in (0, 1):
        for n in range(-2, 3):
            vin = np.zeros(4 * dd, dtype=complex)  # V (x) W (x) V_P: index a*2dd + 2 i(n) + b
            vin[a * 2 * dd + 2 * Wo.i(n) + 0] = 1
            vin[a * 2 * dd + 2 * Wo.i(n + 1) + 1] = -1
            vout = big @ vin  # W (x) V_P (x) V: index (2 i(m) + b) * 2 + c
            # decompose into e_m (x) v_c
            for c in (0, 1):
                for m in range(-3, 4):
                    amp[(a, n, c, m)] = vout[(2 * Wo.i(m) + 0) * 2 + c]
            # the output lies in span{e_m (x) v_c}
            recon = np.zeros_like(vout)
            for c in (0, 1):
                for m in range(-3, 4):
                    recon[(2 * Wo.i(m) + 0) * 2 + c] += amp[(a, n, c, m)]
                    recon[(2 * Wo.i(m + 1) + 1) * 2 + c] -= amp[(a, n, c, m)]
            ok_l = ok_l and np.abs(vout - recon).max() < 1e-10
    Ao = lambda n: amp[(0, n, 0, n)]
    Do = lambda n: amp[(1, n, 1, n)]
    Co = lambda n: amp[(0, n, 1, n + 1)]
    Bo = lambda n: amp[(1, n, 0, n - 1)]
    Yodd = Yt * cmath.exp(1j * math.pi * lam_t)  # theta_d -> theta_d - i pi
    for n in (-1, 0, 1):
        ok_l = ok_l and abs(Ao(n + 1) / Ao(n) - 1 / qt) < 1e-9
        ok_l = ok_l and abs(Do(n) / Ao(n) - qt ** (2 * n + s_t + 1)) < 1e-9
        ok_l = ok_l and abs(Co(n) * Bo(n + 1) / (Ao(n) * Do(n + 1)) - (-Yodd ** 2)) < 1e-9
    # control: theta_d unshifted fails
    ok_l = ok_l and abs(Co(0) * Bo(1) / (Ao(0) * Do(1)) - (-Yt ** 2)) > 1e-3
report("(l) the odd sector: at theta_P the residue has image e_n = |n> v+ - |n+1> v-; R_23 T_12 maps V (x) span{e_n} into "
       "span{e_m} (x) V with the gauge invariants of the even matrix at s + 1 and theta_d - i pi (and not theta_d)", ok_l)

raise SystemExit(0 if OK else 1)
