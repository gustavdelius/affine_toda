"""Checks for Lecture 6 of the defect lecture notes, "The classical algebra of a defect"
(lecture-notes/defects/06-classical-algebra.qmd).

Conventions of the lectures: the dimensionless sine-Gordon field Phi on x < 0 and Psi on x > 0, the action
(1/beta^2) int {int (1/2 Phi_t^2 - 1/2 Phi_x^2 - V(Phi)) dx + (same for Psi) + L_d} dt with V = m0^2 (1 - cos Phi),
d_pm = d_t pm d_x, field equation d_+ d_- Phi = -m0^2 sin Phi. The type-I_- defect has L_d = Phi Psi_t - B,
B = -2 m0 [sigma (cos(Sigma/2) - 1) + sigma^{-1} (cos(Delta/2) - 1)], Sigma = Phi + Psi, Delta = Phi - Psi,
sigma = e^{-theta_d}, P_d = 2 m0 [sigma (cos(Sigma/2) - 1) - sigma^{-1} (cos(Delta/2) - 1)], and the defect conditions
    d_+ Delta = -2 sigma m0 sin(Sigma/2),     d_- Sigma = (2 m0/sigma) sin(Delta/2)      at x = 0.
The 2x2 matrices are h = diag(1, -1), e = [[0, 1], [0, 0]], f = e^T, s = e + f, y = -i e + i f. The Lax pair is
    A_+ = (i/4) d_+Phi h + (m0 lambda/2)(e^{i Phi/2} e + e^{-i Phi/2} f),
    A_- = -(i/4) d_-Phi h - (m0/(2 lambda))(e^{i Phi/2} f + e^{-i Phi/2} e),
a_t = (A_+ + A_-)/2, a_x = (A_+ - A_-)/2, linear problem (d_x + a_x) chi = (d_t + a_t) chi = 0, and the defect matrix
    L_d(lambda) = e^{i Delta h/4} + (sigma/lambda) e^{-i Sigma h/4} s,      chi_Phi(0) = L_d chi_Psi(0).
Canonical brackets {Phi(x), Phi_t(y)} = beta^2 delta(x - y). Written from these definitions, independently of
code/defect_lax_rmatrix.py and defect_classical_an.py. One PASS or FAIL line per check; exit status 1 on
any FAIL. Units m0 = beta = 1 in the numerics.

  (a) Zero curvature: d_+ A_- - d_- A_+ + [A_+, A_-] = -(i/2)(d_+ d_- Phi + m0^2 sin Phi) h; the pair is the book's
      a_1^(1) Lax pair at imaginary coupling (b = i beta_b, sqrt2 beta_b phi = Phi, m = m0/2, H = h/sqrt2); a boost
      by rapidity theta' is compensated by lambda -> e^{theta'} lambda.
  (b) The vacuum: a_x(0) and a_t(0) are multiples of s; at lambda = -e^{-theta} the solutions of the linear problem
      grow and decay like e^{+-(m0/2)(x cosh theta - t sinh theta)} along the eigenvectors (1, +-1) of s, so the ratio
      of the two components is a constant times the soliton's Y.
  (c) Baecklund's transformation is a gauge transformation: with L_d, d_pm L_d = L_d A_pm(Psi) - A_pm(Phi) L_d holds
      whenever d_+ Delta = -2 sigma m0 sin(Sigma/2) and d_- Sigma = (2 m0/sigma) sin(Delta/2) hold (symbolic), so that
      A(Psi) = L_d^{-1} A(Phi) L_d + L_d^{-1} dL_d; it fails for 1.3 sigma in one of the two equations.
  (d) Freezing: for G = G_0 + G_1/lambda, the lambda^1 and lambda^-2 parts of d_t G = G a_t(Psi) - a_t(Phi) G force
      G_0 = e^{i Phi h/4}(p + q s)e^{-i Psi h/4}, G_1 = e^{-i Phi h/4}(p' + q' s)e^{i Psi h/4} (general 2x2 solve);
      the lambda^0 and lambda^-1 parts are: p, q, p', q' constant and four equations, p d_+Delta = -2 m0 q' sin(Sigma/2),
      p' d_-Delta = 2 m0 q sin(Sigma/2), q d_+Sigma = -2 m0 p' sin(Delta/2), q' d_-Sigma = 2 m0 p sin(Delta/2);
      q = p' = 0 gives type I_- with sigma = q'/p, p = q' = 0 gives type I_+ (Psi -> -Psi), three nonzero constants
      force a relation among the defect values, four fix all four first derivatives (rank 4).
  (e) det L_d = 1 - sigma^2/lambda^2; at lambda = -sigma, L_d has rank 1 with kernel (1, e^{i Phi/2}); on the vacuum
      L_d = 1 + (sigma/lambda) s with eigenvalues 1 +- sigma/lambda; their ratio (1 - sigma/lambda)/(1 + sigma/lambda) at
      lambda = -e^{-theta} is Lecture 4's delay z = -coth((theta - theta_d)/2); at lambda = i e^{-theta} it has
      modulus 1 and is Lecture 4's transmission factor (sinh theta_d + i cosh theta)/(cosh theta_d - i sinh theta);
      the gauge e^{i Phi h/4} turns L_d into 1 + (sigma/lambda)(e^{-i Sigma/2} e + e^{i Sigma/2} f).
  (f) The generating function: along the exact passage Phi = 4 arctan Y, Psi = 4 arctan(z Y), the regularized
      monodromy e^{k s x+} T_+ L_d^{-1} T_- e^{-k s x-} (gauge e^{i Phi h/4}, k = (m0/4)(lambda + 1/lambda)) has
      constant diagonal entries in the eigenbasis of s, for real and complex lambda, to 1e-8; with L_d -> 1 or
      sigma -> 1.3 sigma they change by more than 1e-2.
  (g) Charges at lambda -> infinity: the Riccati recursion gives W_1 = (d_+Phi/2m0) y, the order-1/lambda density
      m0/4 - (1/4m0)(1/2 (d_+Phi)^2 + V), a total derivative at order 1/lambda^2, and at order 1/lambda^3 the density
      (1/8m0^3)[(d_+^2 Phi)^2 - (d_+Phi)^4/4 + m0^2 cos Phi (d_+Phi)^2] - d_x(d_+Phi d_+^2 Phi)/(4m0^3) (Euler operator);
      the defect factor gives sigma cos(Sigma/2) at order 1/lambda, so the conserved coefficient is
      sigma - (1/4m0)[beta^2 (E - P) + B - P_d]; at order 1/lambda^3 the defect contribution, with the boundary terms of the
      total derivative, reduces by the first defect condition to
      D_{-3} = (sigma/2m0^2) cos(Sigma/2) d_+Phi d_+Psi + sigma^3 (cos(Sigma/2) - (2/3) cos^3(Sigma/2)).
      Numerically along the passage the spin -1 and spin -3 totals are constant to 1e-10 and fail with L_d -> 1.
  (h) Charges at lambda -> 0: in the gauge e^{-i Phi h/4}, L_d becomes e^{i Delta h/2} + (sigma/lambda) s; the defect
      term at order lambda is cos(Delta/2)/sigma and the bulk density m0/4 - (1/4m0)(1/2 (d_-Phi)^2 + V), so
      E + P + (B + P_d)/beta^2 is conserved; the spin +3 defect contribution is
      D_{+3} = -(1/(2 m0^2 sigma)) cos(Delta/2) d_-Phi d_-Psi + sigma^-3 (cos(Delta/2) - (2/3) cos^3(Delta/2)),
      after the boundary terms of the total derivative d_x(d_-Phi d_-^2 Phi)/(4 m0^3) and the second defect condition;
      numerically both totals are constant along the passage.
  (i) The bulk r-matrix: {a_x(lambda_1; x) (x) a_x(lambda_2; y)} = delta(x - y) [r(lambda_1/lambda_2), a_x(lambda_1) (x) 1
      + 1 (x) a_x(lambda_2)] has the unique solution r(z) = -(beta^2/16(z^2 - 1))[(z^2 + 1) h(x)h + 4 z (e(x)f + f(x)e)]
      among a h(x)h + c (e(x)f + f(x)e); it is the book's r-matrix for h = 2 with beta^2 -> b^2 = -beta^2/2; it solves
      the classical Yang-Baxter equation; r^{t1 t2}(z) = -r(1/z); a lattice transport exp(-Delta a_x) with
      {Phi_n, p_n} = beta^2/Delta satisfies {T (x) T} = c [r, T (x) T] with c -> -1 as Delta -> 0.
  (j) The defect's bracket: L_d = (Phi Psi_t - B)/beta^2 gives {Phi, Psi} = -beta^2; with D = e^{i Delta h/4},
      O = e^{-i Sigma h/4} s: {D_1, O_2} = kappa h(x)h D_1 O_2, {O_1, D_2} = -kappa h(x)h O_1 D_2, kappa = {Phi, Psi}/8,
      [h(x)h, X] = 2 h(x)h X for X = D_1 O_2, O_1 D_2, [Pi, D_1 O_2] = -h(x)h O_1 D_2, [Pi, O_1 D_2] = -h(x)h D_1 O_2
      (Pi = e(x)f + f(x)e); {L_1, L_2} = [r, L_1 L_2] holds exactly for kappa = -beta^2/8, i.e. {Phi, Psi} = -beta^2;
      the r-matrix determined by the defect alone is the bulk one; type I_+ needs {Phi, Psi} = +beta^2, which its
      kinetic term -Phi Psi_t gives; a wrong sign fails.
  (k) Involution: {L^-1 (x) L^-1} = -[r, L^-1 (x) L^-1]; for a chain of three defects with independent canonical pairs,
      {tr M(lambda_1), tr M(lambda_2)} = 0 for M = L^-1 L^-1 L^-1 (and for L L L), and fails for L^-1 L L^-1.
  (l) Two defects: L(sigma_1; Phi_0, Phi_1) L(sigma_2; Phi_1, Phi_12) = L(sigma_2; Phi_0, Phi_2) L(sigma_1; Phi_2, Phi_12)
      when Phi_12 is given by Bianchi's formula with a_k = -sigma_k, for all lambda; a perturbed Phi_12 fails.
  (m) Space as time: for G = int (Phi Psi_t - B) dt, delta G/delta Phi = Phi_x and -delta G/delta Psi = Psi_x are
      the defect conditions; -G generates the canonical map (Phi, -Phi_x) -> (Psi, -Psi_x).
  (n) Exercises and remarks: the type-I_+ matrix s e^{-i Sigma h/4} + (sigma/lambda) e^{-i Delta h/4} satisfies its
      time equation on the type-I_+ conditions; L_d depends on lambda and sigma only through sigma/lambda, which is
      -e^{theta - theta_d} at lambda = -e^{-theta}; L_d(sigma)(Phi, Psi) L_d(-sigma)(Psi, Phi) = det L_d (the inverse
      transformation); the even vacua (0, 0) and (2 pi, 2 pi) give L_d = 1 +- (sigma/lambda) s.
  (o) The a_r^(1) defect in the components of Lecture 5 (real coupling, phi, psi in the hyperplane sum = 0 of R^h,
      f_j = e^{psi_j - phi_{j+1}}, g_j = e^{phi_j - psi_j}): with the Lax pair A_+ = (1/2) d_+phi . H + m lambda
      sum_j e^{alpha_j . phi/2} E_{alpha_j}, A_- = -(1/2) d_-phi . H - (m/lambda) sum_j e^{alpha_j . phi/2} E_{-alpha_j}
      in the defining representation (H . v = diag(v_k), E_{alpha_j} = e_{j,j+1}, E_{-alpha_j} = e_{j+1,j}), zero
      curvature is minus Lecture 5's Toda equation on the diagonal (symbolic, h = 3), and the matrix (L_d)_{kk} = g_k^{1/2},
      (L_d)_{j+1,j} = (sigma/lambda) f_j^{1/2} satisfies d_t L_d = L_d a_t(psi) - a_t(phi) L_d on Lecture 5's light-cone
      conditions d_+(phi_k - psi_k) = -m sigma (f_k - f_{k-1}), d_-(psi_k - phi_{k+1}) = (m/sigma)(g_k - g_{k+1})
      (numeric, h = 2..5, complex lambda); it fails for 1.2 sigma. For h = 2 it is the sinh-Gordon form of L_d.
  (p) The block identities of the freezing calculation, with D = e^{i Delta h/4}, O = e^{-i Sigma h/4} s,
      Ob = e^{i Sigma h/4} s = s e^{-i Sigma h/4} and s_F = e^{i F/2} e + e^{-i F/2} f:
      D s_Psi = s_Phi D, Ob s_Psi = s_Phi Ob, D^-1 s_-Psi = s_-Phi D^-1, O s_-Psi = s_-Phi O,
      O s_Psi - s_Phi O = -2i sin(Sigma/2) h D, D^-1 s_Psi - s_Phi D^-1 = -2i sin(Delta/2) h Ob,
      D s_-Psi - s_-Phi D = 2i sin(Delta/2) h O, Ob s_-Psi - s_-Phi Ob = 2i sin(Sigma/2) h D^-1 (a wrong sign fails);
      A_pm(-Psi) = s A_pm(Psi) s and L_d^+ = L_d(Phi, -Psi) s; L_d = e^{i Delta h/4}(1 + (sigma/lambda) s_-Phi).
  (q) A soliton of rapidity theta carries (1/8m0^3) int q_-+ dx = -(2/3) e^{-+3 theta} (quadrature for four theta;
      -(2/3) e^{-+theta} fails); the odd jump (0, 2 pi) stores D_-+3 = -(2/3) sigma^{+-3} relative to the even
      vacuum, where D_-+3 = sigma^{+-3}/3.
  (r) The a_r^(1) matrix of (o) satisfies d_+ L_d = L_d A_+(psi) - A_+(phi) L_d and the same for d_-, separately, on
      the light-cone conditions (h = 2..5, complex lambda; 1.2 sigma fails), so that it is a gauge transformation, and
      alpha_j . w = psi_j + psi_{j+1} - 2 phi_{j+1} for the book's w = D^T phi - E psi, D = 2(1 - pi)^{-1} (the
      transpose D^T fails for h >= 3).

Runtime: about two minutes.
"""
import cmath
import math

import numpy as np
import sympy as sp
from scipy.integrate import quad
from sympy.calculus.euler import euler_equations
from sympy.physics.quantum import TensorProduct as TP

OK = True


def report(name, ok):
    global OK
    OK = OK and bool(ok)
    print(("PASS " if ok else "FAIL ") + name)


I = sp.I
h = sp.Matrix([[1, 0], [0, -1]])
e = sp.Matrix([[0, 1], [0, 0]])
f = sp.Matrix([[0, 0], [1, 0]])
s = e + f
y = -I * e + I * f
E2 = sp.eye(2)
t, x = sp.symbols('t x', real=True)
lam, m0, sig, beta = sp.symbols('lambda m_0 sigma beta', nonzero=True)


def ex(c):
    return sp.diag(sp.exp(c), sp.exp(-c))


def Xm(F):
    return sp.exp(I * F / 2) * e + sp.exp(-I * F / 2) * f


def Ym(F):
    return sp.exp(I * F / 2) * f + sp.exp(-I * F / 2) * e


def zero_mat(M):
    return all(sp.simplify(sp.expand(z.rewrite(sp.exp))) == 0 for z in M)


# ---------------------------------------------------------------------------------------------------------
# (a) zero curvature, the book's pair, boosts
Phi = sp.Function('Phi')(t, x)
dp = lambda G: sp.diff(G, t) + sp.diff(G, x)
dm = lambda G: sp.diff(G, t) - sp.diff(G, x)


def Ap(F, l=lam):
    return I / 4 * dp(F) * h + m0 * l / 2 * Xm(F)


def Am(F, l=lam):
    return -I / 4 * dm(F) * h - m0 / (2 * l) * Ym(F)


Fc = dp(Am(Phi)) - dm(Ap(Phi)) + Ap(Phi) * Am(Phi) - Am(Phi) * Ap(Phi)
target = -I / 2 * (dp(dm(Phi)) + m0**2 * sp.sin(Phi)) * h
ok_a1 = zero_mat(Fc - target)
# the book's pair for a_1^(1): A_+ = (b/2) d_+phi H + m lam sum c_j e^{b alpha_j phi/2} E_{alpha_j}, alpha_1 = sqrt2,
# alpha_0 = -sqrt2, E_{alpha_1} = e, E_{alpha_0} = f, E_{-alpha_1} = f, E_{-alpha_0} = e, c_j = 1, H = h/sqrt2
bb, phib, mb = sp.symbols('beta_b phi_b m', positive=True)
phiB = sp.Function('phib')(t, x)
b = I * bb
al = sp.sqrt(2)
ApB = b / 2 * dp(phiB) * h / sp.sqrt(2) + mb * lam * (sp.exp(b * al * phiB / 2) * e + sp.exp(-b * al * phiB / 2) * f)
AmB = -b / 2 * dm(phiB) * h / sp.sqrt(2) - mb / lam * (sp.exp(b * al * phiB / 2) * f + sp.exp(-b * al * phiB / 2) * e)
sub = {phiB: Phi / (sp.sqrt(2) * bb), mb: m0 / 2}
ok_a2 = zero_mat((ApB.subs(sub).doit() - Ap(Phi))) and zero_mat((AmB.subs(sub).doit() - Am(Phi)))
# boost: Phi_th(t, x) = Phi(t ch - x sh, x ch - t sh); claim A_pm[Phi_th](lam) = e^{-+th} A_pm[Phi](e^th lam) at x'
thb = sp.Symbol('vartheta', real=True)
Pf = sp.Function('P')
ch_, sh_ = sp.cosh(thb), sp.sinh(thb)
tt, xx = t * ch_ - x * sh_, x * ch_ - t * sh_
Phib = Pf(tt, xx)
u_, v_ = sp.symbols('u v')
P0 = Pf(u_, v_)
dpP = (sp.diff(P0, u_) + sp.diff(P0, v_)).subs({u_: tt, v_: xx})
dmP = (sp.diff(P0, u_) - sp.diff(P0, v_)).subs({u_: tt, v_: xx})
ApOrig = I / 4 * dpP * h + m0 * sp.exp(thb) * lam / 2 * Xm(Phib)
AmOrig = -I / 4 * dmP * h - m0 / (2 * sp.exp(thb) * lam) * Ym(Phib)
ok_a3 = (zero_mat((Ap(Phib) - sp.exp(-thb) * ApOrig).applyfunc(lambda z: sp.simplify(z.doit())))
         and zero_mat((Am(Phib) - sp.exp(thb) * AmOrig).applyfunc(lambda z: sp.simplify(z.doit()))))
report("(a) zero curvature <=> sine-Gordon; the book's a_1^(1) pair at imaginary coupling; boosts rescale lambda",
       ok_a1 and ok_a2 and ok_a3)

# ---------------------------------------------------------------------------------------------------------
# (b) the vacuum
th = sp.Symbol('theta', real=True)
ax0 = (Ap(sp.Integer(0)) - Am(sp.Integer(0))) / 2
at0 = (Ap(sp.Integer(0)) + Am(sp.Integer(0))) / 2
ok_b1 = zero_mat(ax0 - m0 / 4 * (lam + 1 / lam) * s) and zero_mat(at0 - m0 / 4 * (lam - 1 / lam) * s)
kx = (m0 / 4 * (lam + 1 / lam) * x + m0 / 4 * (lam - 1 / lam) * t).subs(lam, -sp.exp(-th))
# chi = exp(-kx s) chi0: along (1, 1) the exponent is -kx, along (1, -1) it is +kx
ok_b2 = sp.simplify((-kx - m0 / 2 * (x * sp.cosh(th) - t * sp.sinh(th))).rewrite(sp.exp)) == 0
# a solution of both equations: chi = exp(-(a_x x + a_t t)) chi0 since they commute
chi = (sp.exp(-kx) * sp.Matrix([1, 1]) * 3 + sp.exp(kx) * sp.Matrix([1, -1]) * 2)
ax0l, at0l = ax0.subs(lam, -sp.exp(-th)), at0.subs(lam, -sp.exp(-th))
ok_b3 = zero_mat(chi.diff(x) + ax0l * chi) and zero_mat(chi.diff(t) + at0l * chi)
ratio = sp.simplify((3 * sp.exp(-kx)) / (2 * sp.exp(kx)))
ok_b4 = sp.simplify(sp.powsimp((ratio / sp.exp(m0 * (x * sp.cosh(th) - t * sp.sinh(th)))).rewrite(sp.exp))
                   - sp.Rational(3, 2)) == 0
report("(b) vacuum: a_x, a_t proportional to s; at lambda = -e^-theta the ratio of components is const * Y",
       ok_b1 and ok_b2 and ok_b3 and ok_b4)

# ---------------------------------------------------------------------------------------------------------
# (c) Baecklund's transformation as a gauge transformation (values and light-cone derivatives as symbols)
P, Q = sp.symbols('Phi Psi')
Pp, Pm, Qp, Qm = sp.symbols('Phi_p Phi_m Psi_p Psi_m')
Pt, Px, Qt, Qx = sp.symbols('Phi_t Phi_x Psi_t Psi_x')


def Ld(l=lam, sg=sig, F=P, G=Q):
    S_, D_ = F + G, F - G
    return ex(I * D_ / 4) + sg / l * ex(-I * S_ / 4) * s


def Apv(F, Fp):
    return I / 4 * Fp * h + m0 * lam / 2 * Xm(F)


def Amv(F, Fm):
    return -I / 4 * Fm * h - m0 / (2 * lam) * Ym(F)


S, Dl = P + Q, P - Q
bt = {Qp: Pp + 2 * sig * m0 * sp.sin(S / 2), Qm: -Pm + 2 * m0 / sig * sp.sin(Dl / 2)}
L0 = Ld()
dLp = L0.diff(P) * Pp + L0.diff(Q) * Qp
dLm = L0.diff(P) * Pm + L0.diff(Q) * Qm
Ep = (dLp - (L0 * Apv(Q, Qp) - Apv(P, Pp) * L0)).subs(bt)
Em = (dLm - (L0 * Amv(Q, Qm) - Amv(P, Pm) * L0)).subs(bt)
ok_c1 = zero_mat(Ep) and zero_mat(Em)
L13 = Ld(sg=sp.Rational(13, 10) * sig)
Ep2 = (L13.diff(P) * Pp + L13.diff(Q) * Qp - (L13 * Apv(Q, Qp) - Apv(P, Pp) * L13)).subs(bt)
ok_c2 = not zero_mat(Ep2)
report("(c) Baecklund's transformation: d_pm L_d = L_d A_pm(Psi) - A_pm(Phi) L_d on the BT; fails for 1.3 sigma",
       ok_c1 and ok_c2)

# ---------------------------------------------------------------------------------------------------------
# (d) freezing: general G = G0 + G1/lambda, time equation at x = 0
def at_v(F, Fx, l=lam):
    return I / 4 * Fx * h + m0 / 4 * (l * Xm(F) - Ym(F) / l)


g = sp.symbols('g0:4')
G0gen = sp.Matrix([[g[0], g[1]], [g[2], g[3]]])
eqs_hi = list(G0gen * Xm(Q) - Xm(P) * G0gen)
sol_hi = sp.solve(eqs_hi, list(g), dict=True)
pp, qq = sp.symbols('p q')
G0form = ex(I * P / 4) * (pp * E2 + qq * s) * ex(-I * Q / 4)
ok_d1 = (len(sol_hi) == 1 and len(G0gen.subs(sol_hi[0]).free_symbols & set(g)) == 2
         and zero_mat(G0form * Xm(Q) - Xm(P) * G0form))
G1gen = sp.Matrix([[g[0], g[1]], [g[2], g[3]]])
sol_lo = sp.solve(list(G1gen * Ym(Q) - Ym(P) * G1gen), list(g), dict=True)
p1, q1 = sp.symbols('p1 q1')
G1form = ex(-I * P / 4) * (p1 * E2 + q1 * s) * ex(I * Q / 4)
ok_d1 = ok_d1 and len(sol_lo) == 1 and zero_mat(G1form * Ym(Q) - Ym(P) * G1form)
# middle orders
ptv, qtv, p1t, q1t = sp.symbols('p_t q_t p1_t q1_t')
Gm_ = G0form + G1form / lam
dGt = (Gm_.diff(P) * Pt + Gm_.diff(Q) * Qt + Gm_.diff(pp) * ptv + Gm_.diff(qq) * qtv
       + Gm_.diff(p1) * p1t + Gm_.diff(q1) * q1t)
Eq = (dGt - (Gm_ * at_v(Q, Qx) - at_v(P, Px) * Gm_)) * lam
coeffs = {}
for i in range(2):
    for j in range(2):
        poly = sp.Poly(sp.expand(Eq[i, j] * lam), lam)
        for (k,), c in poly.terms():
            coeffs[(i, j, k)] = sp.expand(c)
# lambda^1 (k = 3 after multiplying by lambda^2) and lambda^-2 (k = 0) vanish identically
ok_d2 = all(sp.simplify(sp.S(coeffs.get((i, j, k), 0)).rewrite(sp.exp)) == 0 for i in range(2) for j in range(2)
            for k in (0, 3))
Sg_, Dg_ = P + Q, P - Q
dpD = Pt + Px - Qt - Qx
dmD = Pt - Px - Qt + Qx
dpS = Pt + Px + Qt + Qx
dmS = Pt - Px + Qt - Qx
four = [pp * dpD + 2 * m0 * q1 * sp.sin(Sg_ / 2), p1 * dmD - 2 * m0 * qq * sp.sin(Sg_ / 2),
        qq * dpS + 2 * m0 * p1 * sp.sin(Dg_ / 2), q1 * dmS - 2 * m0 * pp * sp.sin(Dg_ / 2)]
mids = [sp.S(coeffs.get((i, j, k), 0)) for i in range(2) for j in range(2) for k in (1, 2)]
mids = [z for z in mids if sp.simplify(z) != 0]
# the eight middle coefficients are an invertible linear combination of the four equations and p_t, q_t, p'_t, q'_t
basis = four + [ptv, qtv, p1t, q1t]
vals_pt = {P: 0.37, Q: -1.21, pp: 1.1, qq: 0.7, p1: 0.4, q1: 1.9, m0: 1.3}
mons = [Pt, Px, Qt, Qx, ptv, qtv, p1t, q1t]


def linvec(z):
    """Coefficients of the linear expression z in the monomials and its constant term, at the numeric point."""
    zn = sp.expand(z.subs(vals_pt))
    const = zn.subs({v: 0 for v in mons})
    return np.array([complex(sp.N(zn.coeff(v))) for v in mons] + [complex(sp.N(const))])


Bv = np.array([linvec(b_) for b_ in basis])          # 8 x 9
Zv = np.array([linvec(z) for z in mids])             # 8 x 9
Cm, res, rk, _ = np.linalg.lstsq(Bv.T, Zv.T, rcond=None)
ok_d3 = len(mids) == 8 and np.max(np.abs(Bv.T @ Cm - Zv.T)) < 1e-12 and np.linalg.matrix_rank(Cm, tol=1e-9) == 8
# type I_-: q = p1 = 0, sigma = q1/p -> the defect conditions
c4 = [z.subs({qq: 0, p1: 0, q1: sig * pp}) for z in four]
solI = sp.solve([c4[0], c4[3]], [Px, Qx], dict=True)[0]
ok_d4 = (sp.simplify(c4[1]) == 0 and sp.simplify(c4[2]) == 0
         and sp.simplify((dpD.subs(solI)) + 2 * sig * m0 * sp.sin(Sg_ / 2)) == 0
         and sp.simplify((dmS.subs(solI)) - 2 * m0 / sig * sp.sin(Dg_ / 2)) == 0)
# the conditions of Lecture 3: Phi_x = Psi_t - B_Phi, Psi_x = Phi_t + B_Psi
Bf = -2 * m0 * (sig * (sp.cos(Sg_ / 2) - 1) + (sp.cos(Dg_ / 2) - 1) / sig)
ok_d4 = ok_d4 and sp.simplify(solI[Px] - (Qt - Bf.diff(P))) == 0 and sp.simplify(solI[Qx] - (Pt + Bf.diff(Q))) == 0
# type I_+: p = q1 = 0, sigma = p1/q -> the I_- conditions with Psi -> -Psi
c4p = [z.subs({pp: 0, q1: 0, p1: sig * qq}) for z in four]
solP = sp.solve([c4p[1], c4p[2]], [Px, Qx], dict=True)[0]
mirror = {Q: -Q, Qt: -Qt}
okP = (sp.simplify(c4p[0]) == 0 and sp.simplify(c4p[3]) == 0
       and sp.simplify(solP[Px] - solI[Px].subs(mirror, simultaneous=True)) == 0
       and sp.simplify(solP[Qx] + solI[Qx].subs(mirror, simultaneous=True)) == 0)
# three nonzero (p1 = 0): the second equation forces q sin(Sigma/2) = 0, a relation among the values
ok_d5 = sp.simplify(four[1].subs(p1, 0) + 2 * m0 * qq * sp.sin(Sg_ / 2)) == 0
# four nonzero: the four equations fix Phi_t, Phi_x, Psi_t, Psi_x
Jac = sp.Matrix([[sp.diff(z, v) for v in [Pt, Px, Qt, Qx]] for z in four])
ok_d5 = ok_d5 and Jac.subs({pp: 1.1, qq: 0.7, p1: 0.4, q1: 1.9}).rank() == 4
report("(d) freezing: extreme orders force the ansatz, middle orders give constants and four equations; "
       "types I_- and I_+ only", ok_d1 and ok_d2 and ok_d3 and ok_d4 and okP and ok_d5)

# ---------------------------------------------------------------------------------------------------------
# (e) det, the zero, the vacuum eigenvalues = delay and transmission
ok_e1 = sp.simplify(L0.det() - (1 - sig**2 / lam**2)) == 0
Lz = L0.subs(lam, -sig)
kern = sp.Matrix([1, sp.exp(I * P / 2)])
ok_e2 = zero_mat(Lz * kern) and Lz.rank(simplify=True) == 1
Lvac = L0.subs({P: 0, Q: 0})
ok_e3 = zero_mat(Lvac - (E2 + sig / lam * s))
ok_e3 = ok_e3 and zero_mat(Lvac * sp.Matrix([1, 1]) - (1 + sig / lam) * sp.Matrix([1, 1]))
thd = sp.Symbol('theta_d', real=True)
zr = ((1 - sig / lam) / (1 + sig / lam)).subs({lam: -sp.exp(-th), sig: sp.exp(-thd)})
ok_e4 = sp.simplify((zr + sp.coth((th - thd) / 2)).rewrite(sp.exp)) == 0
Tr = ((1 - sig / lam) / (1 + sig / lam)).subs({lam: I * sp.exp(-th), sig: sp.exp(-thd)})
Tbook = (sp.sinh(thd) + I * sp.cosh(th)) / (sp.cosh(thd) - I * sp.sinh(th))
ok_e5 = sp.simplify((Tr - Tbook).rewrite(sp.exp)) == 0
ok_e5 = ok_e5 and sp.simplify((Tr * sp.conjugate(Tr)).rewrite(sp.exp)) == 1
Lgauge = sp.simplify(ex(-I * P / 4) * L0 * ex(I * Q / 4))
ok_e6 = zero_mat(Lgauge - (E2 + sig / lam * (sp.exp(-I * S / 2) * e + sp.exp(I * S / 2) * f)))
report("(e) det L_d = 1 - sigma^2/lambda^2, rank 1 at -sigma; vacuum eigenvalue ratio = delay z and "
       "transmission T; gauge form", ok_e1 and ok_e2 and ok_e3 and ok_e4 and ok_e5 and ok_e6)

# ---------------------------------------------------------------------------------------------------------
# numerics for (f)-(h): the exact passage
m0v, thdv, thv, cc = 1.0, 0.2, 0.7, 0.3
sigv = math.exp(-thdv)
zv = -1 / math.tanh((thv - thdv) / 2)
chv, shv = math.cosh(thv), math.sinh(thv)


def fields(tv, xv):
    """Phi (left) and Psi (right) as a whole-line Baecklund pair, with derivatives d_t, d_x."""
    Y = math.exp(m0v * (xv * chv - tv * shv) + cc)
    out = []
    for zz in (1.0, zv):
        w = zz * Y
        F = 4 * math.atan(w)
        dF = 4 * w / (1 + w * w)          # derivative with respect to ln Y
        out.append((F, -m0v * shv * dF, m0v * chv * dF))
    return out


Uh = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
hn = np.diag([1.0, -1.0]).astype(complex)
sn = np.array([[0, 1], [1, 0]], complex)
yn = np.array([[0, -1j], [1j, 0]])
en = np.array([[0, 1], [0, 0]], complex)
fn_ = en.T.copy()


def ax_gauge(F, Ft, Fx, l):
    """a_x in the gauge chi = e^{i Phi h/4} chi~ (all vacua alike)."""
    return (1j / 4 * (Ft + Fx) * hn + m0v * l / 4 * sn
            + m0v / (4 * l) * (np.exp(-1j * F) * en + np.exp(1j * F) * fn_))


def transport(side, tv, l, X=30.0, n=6000):
    """Gauge transport over [-X, 0] (side 0) or [0, X] (side 1), by 4th-order Magnus-like RK steps."""
    xs = np.linspace(-X, 0, n + 1) if side == 0 else np.linspace(0, X, n + 1)
    T = np.eye(2, dtype=complex)
    hstep = xs[1] - xs[0]

    def A(xv):
        F, Ft, Fx = fields(tv, xv)[side]
        return -ax_gauge(F, Ft, Fx, l)
    for k in range(n):
        x0 = xs[k]
        k1 = A(x0) @ T
        k2 = A(x0 + hstep / 2) @ (T + hstep / 2 * k1)
        k3 = A(x0 + hstep / 2) @ (T + hstep / 2 * k2)
        k4 = A(x0 + hstep) @ (T + hstep * k3)
        T = T + hstep / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return T, xs[0], xs[-1]


def expm2(M):
    w, V = np.linalg.eig(M)
    return V @ np.diag(np.exp(w)) @ np.linalg.inv(V)


_cache = {}


def monodromy22(tv, l, mode="L"):
    """The dominant diagonal entry of the regularized monodromy in the eigenbasis of s."""
    if (tv, l) not in _cache:
        _cache[(tv, l)] = (transport(0, tv, l), transport(1, tv, l))
    (Tm, a0, _), (Tp, _, b1) = _cache[(tv, l)]
    F0 = fields(tv, 0.0)
    Sg = F0[0][0] + F0[1][0]
    fac = {"L": 1.0, "wrong": 1.3, "one": 0.0}[mode]
    Lt = np.eye(2) + fac * sigv / l * (np.exp(-1j * Sg / 2) * en + np.exp(1j * Sg / 2) * fn_)
    k = m0v / 4 * (l + 1 / l)
    M = Uh @ Tp @ np.linalg.inv(Lt) @ Tm @ Uh        # in the eigenbasis of s, where the regularization is diagonal
    return np.exp(-k * (b1 - a0)) * M[1, 1]


times = (-6.0, -1.5, 0.0, 1.0, 5.0)
ok_f = True
for l in (2.0, 0.5, 1.3 + 0.7j):
    vals = {mode: np.array([monodromy22(tv, l, mode) for tv in times]) for mode in ("L", "one", "wrong")}
    spreadL = np.max(np.abs(vals["L"] - vals["L"][0])) / abs(vals["L"][0])
    spread1 = np.max(np.abs(vals["one"] - vals["one"][0])) / abs(vals["one"][0])
    spreadw = np.max(np.abs(vals["wrong"] - vals["wrong"][0])) / abs(vals["wrong"][0])
    ok_f = ok_f and spreadL < 1e-8 and spread1 > 1e-2 and spreadw > 1e-2
report("(f) the regularized monodromy has constant diagonal entries along a passage, real and complex lambda; "
       "fails with L_d -> 1 or 1.3 sigma", ok_f)

# ---------------------------------------------------------------------------------------------------------
# (g), (h): charges from the asymptotic expansion
Fx_ = sp.Function('F')(t, x)
NW = 4


def riccati(u, sgn):
    """W_k, Z_k for a^ = eps^-1 (m0/4) h + (i/4) u s + eps (m0/4)(cos F h + sgn sin F y)."""
    A = {-1: m0 / 4 * h, 0: I / 4 * u * s, 1: m0 / 4 * (sp.cos(Fx_) * h + sgn * sp.sin(Fx_) * y)}
    dg = lambda M: sp.diag(M[0, 0], M[1, 1])
    of = lambda M: M - dg(M)
    W = {}
    for k in range(0, NW):
        rhs = sp.zeros(2)
        if k >= 1:
            rhs += -W[k].diff(x)
        rhs += -of(A.get(k, sp.zeros(2)))
        if k >= 2:
            rhs += -(dg(A[1]) * W[k - 1] - W[k - 1] * dg(A[1]))
        for a_ in range(1, k + 2):
            for c_ in range(1, k + 2):
                b_ = k - a_ - c_
                if b_ in A and a_ in W and c_ in W:
                    rhs += W[a_] * of(A[b_]) * W[c_]
        W[k + 1] = sp.Matrix([[0, rhs[0, 1] / (m0 / 2)], [-rhs[1, 0] / (m0 / 2), 0]]).applyfunc(sp.simplify)
    Z = {}
    for k in range(-1, NW):
        M = -dg(A.get(k, sp.zeros(2)))
        for b_ in A:
            c_ = k - b_
            if c_ in W:
                M += -dg(of(A[b_]) * W[c_])
        Z[k] = sp.simplify(M[1, 1])
    return W, Z


uplus = dp(Fx_)
uminus = dm(Fx_)
Winf, Zinf = riccati(uplus, -1)
Wzero, Zzero = riccati(uminus, +1)
Vf = m0**2 * (1 - sp.cos(Fx_))
ok_g1 = (zero_mat(Winf[1] - uplus / (2 * m0) * y)
         and sp.simplify(Zinf[1] - (m0 / 4 - (uplus**2 / 2 + Vf) / (4 * m0))) == 0
         and sp.simplify(Zinf[2] - sp.diff(uplus**2, x) / (8 * m0**2)) == 0)
# order 3: density q/(8 m0^3) - d_x(u v)/(4 m0^3), with u = d_+F, v = d_+^2 F, on shell (F_tt = F_xx - m0^2 sin F)
onshell = lambda z: z.subs(sp.Derivative(Fx_, (t, 2)), sp.diff(Fx_, x, 2) - m0**2 * sp.sin(Fx_))


def jetsub(z):
    """Replace t-derivatives by x-derivatives of F and of p = F_t, on shell."""
    pf = sp.Function('p')(x)
    Fj = sp.Function('Fj')(x)
    z = z.doit()
    rep = {}
    for n in range(4, -1, -1):
        for k in range(0, 5 - n):
            if n == 0 and k == 0:
                continue
            d = sp.Derivative(Fx_, *([(t, n)] if n else []), *([(x, k)] if k else []))
            if n == 0:
                rep[d] = sp.diff(Fj, x, k)
            elif n == 1:
                rep[d] = sp.diff(pf, x, k)
    # reduce t-derivatives of order >= 2 on shell repeatedly
    for _ in range(4):
        z = z.subs(sp.Derivative(Fx_, (t, 2)), sp.diff(Fx_, x, 2) - m0**2 * sp.sin(Fx_)).doit()
        z = z.subs(sp.Derivative(Fx_, (t, 2), (x, 1)), sp.diff(sp.diff(Fx_, x, 2) - m0**2 * sp.sin(Fx_), x)).doit()
        z = z.subs(sp.Derivative(Fx_, (t, 3)), sp.diff(sp.diff(Fx_, x, 2) - m0**2 * sp.sin(Fx_), t)).doit()
    z = z.subs(rep).subs(Fx_, Fj)
    return z, Fj, pf


vplus = dp(uplus)
q3 = vplus**2 - uplus**4 / 4 + m0**2 * sp.cos(Fx_) * uplus**2
diff3 = Zinf[3] - q3 / (8 * m0**3) + sp.diff(uplus * vplus, x) / (4 * m0**3)
dj, Fj, pf = jetsub(diff3)
ok_g2 = sp.simplify(dj) == 0
# the general defect factor ln N_22 to order 3, from W_1, W_2 and the rotated defect matrix
mu = sp.Symbol('mu')
K = 4
aa = sp.symbols('a1:5')
bbs = sp.symbols('b1:5')
cs = sp.symbols('c1:5')
ds = sp.symbols('d1:5')
Sgs = sp.Symbol('Sigma')


def trunc(z):
    z = sp.expand(z)
    return sum(z.coeff(mu, k) * mu**k for k in range(K + 1))


def mm(A, B):
    return (A * B).applyfunc(trunc)


def minv_ser(M):
    det = trunc(M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0])
    di = trunc(sp.series(1 / det, mu, 0, K + 1).removeO())
    return sp.Matrix([[M[1, 1], -M[0, 1]], [-M[1, 0], M[0, 0]]]).applyfunc(lambda z: trunc(z * di))


Wl = sp.Matrix([[1, sum(aa[k - 1] * mu**k for k in range(1, K + 1))], [sum(bbs[k - 1] * mu**k for k in range(1, K + 1)), 1]])
Wr = sp.Matrix([[1, sum(cs[k - 1] * mu**k for k in range(1, K + 1))], [sum(ds[k - 1] * mu**k for k in range(1, K + 1)), 1]])
Uq = sp.Matrix([[1, 1], [1, -1]]) / sp.sqrt(2)
Lhat = (Uq * (E2 + sig * mu * (sp.cos(Sgs / 2) * s + sp.sin(Sgs / 2) * y)) * Uq).applyfunc(sp.simplify)
ok_g3 = zero_mat(Lhat - (E2 + sig * mu * (sp.cos(Sgs / 2) * h - sp.sin(Sgs / 2) * y)))
Nm = mm(mm(minv_ser(Wr), minv_ser(Lhat.applyfunc(trunc))), Wl)
lnN = trunc(sp.series(sp.log(Nm[1, 1]), mu, 0, K + 1).removeO())
d1 = sp.simplify(lnN.coeff(mu, 1))
d3 = sp.expand(lnN.coeff(mu, 3))
ok_g3 = ok_g3 and sp.simplify(d1 - sig * sp.cos(Sgs / 2)) == 0
# substitute W_1 = (u/2m0) y, W_2 = (i/2m0^2) v s on each side (W_2 on shell)
uL, uR, vL, vR = sp.symbols('u_L u_R v_L v_R')
W2check = sp.simplify(jetsub(Winf[2] - I / (2 * m0**2) * vplus * s)[0])
ok_g4 = zero_mat(W2check)
wsub = {aa[0]: -I * uL / (2 * m0), bbs[0]: I * uL / (2 * m0), cs[0]: -I * uR / (2 * m0), ds[0]: I * uR / (2 * m0),
        aa[1]: I * vL / (2 * m0**2), bbs[1]: I * vL / (2 * m0**2), cs[1]: I * vR / (2 * m0**2),
        ds[1]: I * vR / (2 * m0**2)}
d3s = sp.expand(d3.subs(wsub))
# boundary terms of the total derivative: + J(0^-) - J(0^+), J = -u v/(4 m0^3)
D3 = d3s - uL * vL / (4 * m0**3) + uR * vR / (4 * m0**3)
D3 = sp.expand(D3.subs(uR, uL + 2 * sig * m0 * sp.sin(Sgs / 2)))
D3claim = sp.expand((sig / (2 * m0**2) * sp.cos(Sgs / 2) * uL * (uL + 2 * sig * m0 * sp.sin(Sgs / 2))
                     + sig**3 * (sp.cos(Sgs / 2) - sp.Rational(2, 3) * sp.cos(Sgs / 2)**3)))
ok_g5 = sp.simplify(D3 - D3claim) == 0
report("(g) symbolic expansion at lambda -> oo: W_1, the E - P density, a total derivative, the spin -3 density "
       "(Euler operator), the defect contributions sigma cos(Sigma/2) and D_-3", ok_g1 and ok_g2 and ok_g3 and ok_g4 and ok_g5)


# numerics of the charges along the passage
def lam_fun(expr, side):
    """Numeric function of (t, x) for an expression in F(t, x), on the left (Phi) or right (Psi) solution."""
    Y = sp.exp(m0v * (x * sp.cosh(thv) - t * sp.sinh(thv)) + cc)
    G = 4 * sp.atan((1 if side == 0 else zv) * Y)
    return sp.lambdify((t, x), expr.subs(m0, m0v).subs(Fx_, G).doit(), 'numpy')


_lf = {}


def bulk_int(expr, tv, vac=0.0, X=30.0):
    key = sp.srepr(expr)
    if key not in _lf:
        _lf[key] = (lam_fun(expr, 0), lam_fun(expr, 1))
    fl, fr = _lf[key]
    l = quad(lambda xv: complex(fl(tv, xv) - vac).real, -X, 0, limit=400)[0]
    r = quad(lambda xv: complex(fr(tv, xv) - vac).real, 0, X, limit=400)[0]
    return l, r


Vn = lambda F: m0v**2 * (1 - np.cos(F))
dens_m1 = (uplus**2 / 2 + m0**2 * (1 - sp.cos(Fx_)))
dens_p1 = (uminus**2 / 2 + m0**2 * (1 - sp.cos(Fx_)))
dens_m3 = q3 - 0 * Fx_
vminus = dm(uminus)
dens_p3 = vminus**2 - uminus**4 / 4 + m0**2 * sp.cos(Fx_) * uminus**2
upL = lam_fun(uplus, 0)
upR = lam_fun(uplus, 1)
umL = lam_fun(uminus, 0)
umR = lam_fun(uminus, 1)


def defect_contributions(tv, mode="L"):
    F0 = fields(tv, 0.0)
    Sg, Dg = F0[0][0] + F0[1][0], F0[0][0] - F0[1][0]
    sg = sigv if mode == "L" else 0.0
    B = -2 * m0v * (sg * (math.cos(Sg / 2) - 1) + ((math.cos(Dg / 2) - 1) / sigv if mode == "L" else 0.0))
    Pd = 2 * m0v * (sg * (math.cos(Sg / 2) - 1) - ((math.cos(Dg / 2) - 1) / sigv if mode == "L" else 0.0))
    uL_, uR_ = complex(upL(tv, 0.0)).real, complex(upR(tv, 0.0)).real
    wL_, wR_ = complex(umL(tv, 0.0)).real, complex(umR(tv, 0.0)).real
    if mode == "L":
        Dm3 = sigv / (2 * m0v**2) * math.cos(Sg / 2) * uL_ * uR_ + sigv**3 * (math.cos(Sg / 2) - 2 / 3 * math.cos(Sg / 2)**3)
        Dp3 = (-1 / (2 * m0v**2 * sigv) * math.cos(Dg / 2) * wL_ * wR_
               + sigv**-3 * (math.cos(Dg / 2) - 2 / 3 * math.cos(Dg / 2)**3))
    else:
        Dm3 = Dp3 = 0.0
    return B - Pd, B + Pd, Dm3, Dp3


ok_num = True
tot = {k: [] for k in ("m1", "p1", "m3", "p3")}
tot1 = {k: [] for k in ("m1", "p1", "m3", "p3")}
for tv in (-6.0, -2.0, -0.5, 0.0, 0.7, 2.0, 6.0):
    l, r = bulk_int(dens_m1, tv)
    w = defect_contributions(tv)
    w1 = defect_contributions(tv, "none")
    tot["m1"].append(l + r + w[0])
    tot1["m1"].append(l + r + w1[0])
    l, r = bulk_int(dens_p1, tv)
    tot["p1"].append(l + r + w[1])
    tot1["p1"].append(l + r + w1[1])
    l, r = bulk_int(dens_m3, tv)
    tot["m3"].append((l + r) / (8 * m0v**3) + w[2])
    tot1["m3"].append((l + r) / (8 * m0v**3) + w1[2])
    l, r = bulk_int(dens_p3, tv)
    tot["p3"].append((l + r) / (8 * m0v**3) + w[3])
    tot1["p3"].append((l + r) / (8 * m0v**3) + w1[3])
spreads = {k: np.ptp(tot[k]) for k in tot}
spreads1 = {k: np.ptp(tot1[k]) for k in tot1}
ok_num = all(v < 1e-8 for v in spreads.values()) and all(v > 1e-2 for v in spreads1.values())
# the spin -1 total is beta^2 (E - P) + B - P_d = 8 m0 cosh... : its value for one soliton is M e^{-theta} beta^2 = 8 e^-theta
ok_num = ok_num and abs(tot["m1"][0] - 8 * m0v * math.exp(-thv)) < 1e-6 and abs(tot["p1"][0] - 8 * m0v * math.exp(thv)) < 1e-6
# lambda -> 0 side, symbolic: the same steps in the gauge e^{-i Phi h/4}
ok_h1 = (sp.simplify(Zzero[1] - (m0 / 4 - (uminus**2 / 2 + Vf) / (4 * m0))) == 0)
diffp3 = Zzero[3] - dens_p3 / (8 * m0**3) - sp.diff(uminus * vminus, x) / (4 * m0**3)
ok_h1 = ok_h1 and sp.simplify(jetsub(diffp3)[0]) == 0
ok_h1 = ok_h1 and zero_mat(Wzero[1] - uminus / (2 * m0) * y) and zero_mat(jetsub(Wzero[2] + I / (2 * m0**2) * vminus * s)[0])
Lp0 = sp.simplify(ex(I * P / 4) * L0 * ex(-I * Q / 4))
ok_h1 = ok_h1 and zero_mat(Lp0 - (ex(I * Dl / 2) + sig / lam * s))
# rotated: L^' = (sigma/lambda) h (1 + (lambda/sigma) h M0), M0 = U e^{i Delta h/2} U; here mu stands for lambda
Dls = sp.Symbol('Delta')
M0 = (Uq * ex(I * Dls / 2) * Uq).applyfunc(sp.simplify)
inner = (E2 + mu / sig * h * M0).applyfunc(trunc)
Npp = mm(mm(minv_ser(Wr), minv_ser(inner)), h * Wl)
lnNp = trunc(sp.series(sp.log(-Npp[1, 1]), mu, 0, K + 1).removeO())
ok_h1 = ok_h1 and sp.simplify(Npp[1, 1].coeff(mu, 0) + 1) == 0 and sp.simplify(lnNp.coeff(mu, 1) - sp.cos(Dls / 2) / sig) == 0
wL, wR = sp.symbols('w_L w_R')
wsub0 = {aa[0]: -I * wL / (2 * m0), bbs[0]: I * wL / (2 * m0), cs[0]: -I * wR / (2 * m0), ds[0]: I * wR / (2 * m0),
         aa[1]: -I * vL / (2 * m0**2), bbs[1]: -I * vL / (2 * m0**2), cs[1]: -I * vR / (2 * m0**2),
         ds[1]: -I * vR / (2 * m0**2)}
D3p = sp.expand(lnNp.coeff(mu, 3).subs(wsub0)) + wL * vL / (4 * m0**3) - wR * vR / (4 * m0**3)
D3p = sp.expand(D3p.subs(wR, -wL + 2 * m0 / sig * sp.sin(Dls / 2)))
D3pclaim = sp.expand(-1 / (2 * m0**2 * sig) * sp.cos(Dls / 2) * wL * (-wL + 2 * m0 / sig * sp.sin(Dls / 2))
                     + sig**-3 * (sp.cos(Dls / 2) - sp.Rational(2, 3) * sp.cos(Dls / 2)**3))
ok_h1 = ok_h1 and sp.simplify(sp.expand_trig(D3p - D3pclaim)) == 0
report("(g) numerics: the spin -1 and spin -3 totals are constant along the passage (1e-8), equal to 8 e^-theta "
       "for spin -1, and fail without the defect contributions", ok_num)
report("(h) lambda -> 0: L_d -> e^{i Delta h/2} + (sigma/lambda) s, the E + P density and the spin +3 density; "
       "E + P + (B + P_d) and the spin +3 total with D_+3 constant along the passage",
       ok_h1 and spreads["p1"] < 1e-8 and spreads["p3"] < 1e-8)

# ---------------------------------------------------------------------------------------------------------
# (i) the bulk r-matrix
l1, l2, zz_ = sp.symbols('lambda_1 lambda_2 z', nonzero=True)
Pv, Ptv = sp.symbols('Phi Phi_t')


def axv(l):
    return I / 4 * Ptv * h + m0 / 4 * (l * Xm(Pv) + Ym(Pv) / l)


def br_ul(A, B):
    return beta**2 * (TP(A.diff(Pv), B.diff(Ptv)) - TP(A.diff(Ptv), B.diff(Pv)))


alr, gar = sp.symbols('alpha gamma')
Pi = TP(e, f) + TP(f, e)
rg = alr * TP(h, h) + gar * Pi
Br = br_ul(axv(l1), axv(l2))
Cr = rg * (TP(axv(l1), E2) + TP(E2, axv(l2))) - (TP(axv(l1), E2) + TP(E2, axv(l2))) * rg
eqr = [sp.simplify(q) for q in (Br - Cr) if sp.simplify(q) != 0]
solr = sp.solve(eqr, [alr, gar], dict=True)


def rmat(zv_):
    return -beta**2 / (16 * (zv_**2 - 1)) * ((zv_**2 + 1) * TP(h, h) + 4 * zv_ * Pi)


ok_i1 = (len(solr) == 1 and sp.simplify(solr[0][alr].subs(l1, zz_ * l2) + beta**2 * (zz_**2 + 1) / (16 * (zz_**2 - 1))) == 0
         and sp.simplify(solr[0][gar].subs(l1, zz_ * l2) + beta**2 * zz_ / (4 * (zz_**2 - 1))) == 0)
# the book's formula: r = b^2/(4(z^h-1)) [(z^h+1) H(x)H + 2 (z E(x)F + z F(x)E)], H = h/sqrt2, h = 2, b^2 = -beta^2/2
bsq = -beta**2 / 2
rbook = bsq / (4 * (zz_**2 - 1)) * ((zz_**2 + 1) * TP(h, h) / 2 + 2 * (zz_ * TP(e, f) + zz_ * TP(f, e)))
ok_i2 = zero_mat(rbook - rmat(zz_))
# r^{t1 t2}(z) = -r(1/z)
def pt12(M):
    out = sp.zeros(4)
    for i1 in range(2):
        for j1 in range(2):
            for i2 in range(2):
                for j2 in range(2):
                    out[2 * j1 + j2, 2 * i1 + i2] = M[2 * i1 + i2, 2 * j1 + j2]
    return out


ok_i3 = zero_mat(pt12(rmat(zz_)) + rmat(1 / zz_))
# classical Yang-Baxter equation, numerically
def r_num(zv_):
    return np.array(rmat(zv_).subs(beta, 1), dtype=complex)


def emb(M, i, j):
    """Embed a 4x4 two-site operator into sites (i, j) of three."""
    Ms = M.reshape(2, 2, 2, 2)
    out = np.zeros((2,) * 6, complex)
    for a in range(2):
        for b_ in range(2):
            for c in range(2):
                for d in range(2):
                    for k in range(2):
                        idx_in = [None] * 3
                        idx_out = [None] * 3
                        idx_in[i], idx_in[j] = c, d
                        idx_out[i], idx_out[j] = a, b_
                        other = 3 - i - j
                        idx_in[other] = idx_out[other] = k
                        out[tuple(idx_out) + tuple(idx_in)] += Ms[a, b_, c, d]
    return out.reshape(8, 8)


z1, z2, z3 = 1.7 + 0.3j, 0.6 - 0.2j, 1.1 + 0.9j
R12, R13, R23 = emb(r_num(z1 / z2), 0, 1), emb(r_num(z1 / z3), 0, 2), emb(r_num(z2 / z3), 1, 2)
cyb = (R12 @ R13 - R13 @ R12) + (R12 @ R23 - R23 @ R12) + (R13 @ R23 - R23 @ R13)
ok_i4 = np.max(np.abs(cyb)) < 1e-12
# lattice sign: T = exp(-Delta a_x), {Phi, p} = beta^2/Delta with p = Phi_t
def lattice_c(Delta, lv1=1.4, lv2=0.6 + 0.3j):
    Pn, Ptn = 0.83, -0.57
    def T(lv, Pv_, Ptv_):
        A = np.array(axv(lv).subs({Pv: Pv_, Ptv: Ptv_, m0: 1}), dtype=complex)
        return expm2(-Delta * A)
    eps_ = 1e-6
    def dT(lv, which):
        if which == 0:
            return (T(lv, Pn + eps_, Ptn) - T(lv, Pn - eps_, Ptn)) / (2 * eps_)
        return (T(lv, Pn, Ptn + eps_) - T(lv, Pn, Ptn - eps_)) / (2 * eps_)
    br = (np.kron(dT(lv1, 0), dT(lv2, 1)) - np.kron(dT(lv1, 1), dT(lv2, 0))) / Delta
    TT = np.kron(T(lv1, Pn, Ptn), T(lv2, Pn, Ptn))
    rr = r_num(lv1 / lv2)
    com = rr @ TT - TT @ rr
    return np.vdot(com.ravel(), br.ravel()) / np.vdot(com.ravel(), com.ravel())


c1, c2 = lattice_c(0.02), lattice_c(0.01)
ok_i5 = abs(c1 + 1) < 0.05 and abs(c2 + 1) < abs(c1 + 1)
report("(i) bulk r-matrix from the canonical brackets; the book's r for h = 2; r^{t1t2}(z) = -r(1/z); CYBE; "
       "lattice transport bracket -> -[r, T(x)T]", ok_i1 and ok_i2 and ok_i3 and ok_i4 and ok_i5)

# ---------------------------------------------------------------------------------------------------------
# (j) the defect's bracket
# symplectic form of (1/beta^2)(Phi Psi_t - B): Theta = (1/beta^2) Phi dPsi, omega = (1/beta^2) dPhi ^ dPsi,
# {z_I, z_J} = (Omega^{-1})_{IJ} with Omega_{IJ} = omega(d_I, d_J) in the convention {q, p} = 1 for p dq
qv, pvv = sp.symbols('q p_')
# check {q, p} = 1 convention with the general formula {F, G} = dF/dq dG/dp - dF/dp dG/dq, Phi = beta^2 p, Psi = q
PhiE_, PsiE_ = beta**2 * pvv, qv
pb = sp.diff(PhiE_, qv) * sp.diff(PsiE_, pvv) - sp.diff(PhiE_, pvv) * sp.diff(PsiE_, qv)
ok_j1 = sp.simplify(pb + beta**2) == 0
kap = sp.Symbol('kappa')
Dm_ = ex(I * Dl / 4)
Om_ = ex(-I * S / 4) * s


def brPQ(A, B, c):
    return c * (TP(A.diff(P), B.diff(Q)) - TP(A.diff(Q), B.diff(P)))


hh = TP(h, h)
ok_j2 = (zero_mat(brPQ(Dm_, Om_, 8 * kap) - kap * hh * TP(Dm_, Om_))
         and zero_mat(brPQ(Om_, Dm_, 8 * kap) + kap * hh * TP(Om_, Dm_))
         and zero_mat(brPQ(Dm_, Dm_, 1)) and zero_mat(brPQ(Om_, Om_, 1)))
X1, X2 = TP(Dm_, Om_), TP(Om_, Dm_)
ok_j3 = (zero_mat(hh * X1 - X1 * hh - 2 * hh * X1) and zero_mat(hh * X2 - X2 * hh - 2 * hh * X2)
         and zero_mat(Pi * X1 - X1 * Pi + hh * X2) and zero_mat(Pi * X2 - X2 * Pi + hh * X1))
L1, L2 = Ld(l1), Ld(l2)
LHS = brPQ(L1, L2, 8 * kap)
rz = rmat(l1 / l2)
RHS = rz * TP(L1, L2) - TP(L1, L2) * rz
solk = sp.solve([q for q in (LHS - RHS).applyfunc(lambda z: sp.simplify(z.rewrite(sp.exp))) if q != 0], kap, dict=True)
ok_j4 = len(solk) == 1 and sp.simplify(solk[0][kap] + beta**2 / 8) == 0
# the r-matrix determined by the defect alone: alpha h(x)h + gamma Pi with 2 alpha - gamma/z = kappa, 2 alpha - gamma z = -kappa
z_ = sp.Symbol('z')
sw = sp.solve([2 * alr - gar / z_ - kap, 2 * alr - gar * z_ + kap], [alr, gar], dict=True)[0]
ok_j5 = (sp.simplify(sw[alr].subs(kap, -beta**2 / 8) + beta**2 * (z_**2 + 1) / (16 * (z_**2 - 1))) == 0
         and sp.simplify(sw[gar].subs(kap, -beta**2 / 8) + beta**2 * z_ / (4 * (z_**2 - 1))) == 0)
# type I_+: matrix s e^{-i Sigma h/4} + (sigma/lambda) e^{-i Delta h/4}; needs {Phi, Psi} = +beta^2
Lp = lambda l: s * ex(-I * S / 4) + sig / l * ex(-I * Dl / 4)
LHSp = brPQ(Lp(l1), Lp(l2), 8 * kap)
RHSp = rz * TP(Lp(l1), Lp(l2)) - TP(Lp(l1), Lp(l2)) * rz
solkp = sp.solve([q for q in (LHSp - RHSp).applyfunc(lambda z: sp.simplify(z.rewrite(sp.exp))) if q != 0], kap, dict=True)
ok_j6 = len(solkp) == 1 and sp.simplify(solkp[0][kap] - beta**2 / 8) == 0
# the wrong sign fails
ok_j7 = not zero_mat(brPQ(L1, L2, beta**2) - RHS)
report("(j) {Phi, Psi} = -beta^2 from the kinetic term; brackets of D, O; the four commutator identities; "
       "{L_1, L_2} = [r, L_1 L_2] iff {Phi, Psi} = -beta^2; the defect fixes r; type I_+ needs +beta^2",
       ok_j1 and ok_j2 and ok_j3 and ok_j4 and ok_j5 and ok_j6 and ok_j7)

# ---------------------------------------------------------------------------------------------------------
# (k) involution for a chain of defects
rngk = np.random.default_rng(7)
vals3 = rngk.uniform(-2, 2, (3, 2))
sigs3 = [0.7, 1.3, 0.45]
Ld_np = lambda l, sg, Pv_, Qv_: (np.diag([np.exp(1j * (Pv_ - Qv_) / 4), np.exp(-1j * (Pv_ - Qv_) / 4)])
                                 + sg / l * np.diag([np.exp(-1j * (Pv_ + Qv_) / 4), np.exp(1j * (Pv_ + Qv_) / 4)]) @ sn)


def chain(l, which, pert=None):
    M = np.eye(2, dtype=complex)
    for k in range(3):
        Pv_, Qv_ = vals3[k]
        if pert is not None and pert[0] == k:
            if pert[1] == 0:
                Pv_ += pert[2]
            else:
                Qv_ += pert[2]
        Lk = Ld_np(l, sigs3[k], Pv_, Qv_)
        M = M @ (np.linalg.inv(Lk) if which[k] else Lk)
    return M


def trace_bracket(la, lb, which, c=-1.0):
    epsd = 1e-6
    tot_ = 0
    for k in range(3):
        dA = [(np.trace(chain(la, which, (k, j, epsd))) - np.trace(chain(la, which, (k, j, -epsd)))) / (2 * epsd)
              for j in (0, 1)]
        dB = [(np.trace(chain(lb, which, (k, j, epsd))) - np.trace(chain(lb, which, (k, j, -epsd)))) / (2 * epsd)
              for j in (0, 1)]
        tot_ += c * (dA[0] * dB[1] - dA[1] * dB[0])
    return abs(tot_)


la_, lb_ = 0.9 + 0.2j, -1.6 + 0.5j
ok_k1 = (trace_bracket(la_, lb_, (1, 1, 1)) < 1e-7 and trace_bracket(la_, lb_, (0, 0, 0)) < 1e-7
         and trace_bracket(la_, lb_, (1, 0, 1)) > 1e-3)
# {L^-1 (x) L^-1} = -[r, L^-1 (x) L^-1] numerically at a point
Pn, Qn = 0.41, -1.07
def Linv(l, Pv_, Qv_):
    return np.linalg.inv(Ld_np(l, 0.8, Pv_, Qv_))
epsd = 1e-6
dPa = (Linv(la_, Pn + epsd, Qn) - Linv(la_, Pn - epsd, Qn)) / (2 * epsd)
dQa = (Linv(la_, Pn, Qn + epsd) - Linv(la_, Pn, Qn - epsd)) / (2 * epsd)
dPb = (Linv(lb_, Pn + epsd, Qn) - Linv(lb_, Pn - epsd, Qn)) / (2 * epsd)
dQb = (Linv(lb_, Pn, Qn + epsd) - Linv(lb_, Pn, Qn - epsd)) / (2 * epsd)
brI = -1.0 * (np.kron(dPa, dQb) - np.kron(dQa, dPb))
TTi = np.kron(Linv(la_, Pn, Qn), Linv(lb_, Pn, Qn))
rr_ = r_num(la_ / lb_)
ok_k2 = np.max(np.abs(brI + (rr_ @ TTi - TTi @ rr_))) < 1e-8
report("(k) {L^-1, L^-1} = -[r, L^-1 L^-1]; traces of chains of three defects Poisson-commute; mixing L and "
       "L^-1 fails", ok_k1 and ok_k2)

# ---------------------------------------------------------------------------------------------------------
# (l) two defects: Bianchi as a refactorization
rngl = np.random.default_rng(3)
ok_l = True
for trial in range(4):
    s1, s2 = rngl.uniform(0.3, 2.0, 2)
    a1, a2 = -s1, -s2
    P0, P1, P2 = rngl.uniform(-3, 3, 3)
    P12 = P0 + 4 * math.atan((a1 + a2) / (a1 - a2) * math.tan((P1 - P2) / 4))
    for lv in (0.7, -1.3 + 0.4j, 2.2j):
        A_ = Ld_np(lv, s1, P0, P1) @ Ld_np(lv, s2, P1, P12)
        B_ = Ld_np(lv, s2, P0, P2) @ Ld_np(lv, s1, P2, P12)
        ok_l = ok_l and np.max(np.abs(A_ - B_)) < 1e-12
    A_ = Ld_np(0.7, s1, P0, P1) @ Ld_np(0.7, s2, P1, P12 + 0.01)
    B_ = Ld_np(0.7, s2, P0, P2) @ Ld_np(0.7, s1, P2, P12 + 0.01)
    ok_l = ok_l and np.max(np.abs(A_ - B_)) > 1e-4
report("(l) two defects exchange: L(s1; Phi0, Phi1) L(s2; Phi1, Phi12) = L(s2; Phi0, Phi2) L(s1; Phi2, Phi12) "
       "for Bianchi's Phi12; a perturbed Phi12 fails", ok_l)

# ---------------------------------------------------------------------------------------------------------
# (m) space as time
Pf_, Qf_ = sp.Function('Phi0')(t), sp.Function('Psi0')(t)
Bt = -2 * m0 * (sig * (sp.cos((Pf_ + Qf_) / 2) - 1) + (sp.cos((Pf_ - Qf_) / 2) - 1) / sig)
Gd = Pf_ * Qf_.diff(t) - Bt
EL = euler_equations(Gd, [Pf_, Qf_], t)   # dG/dPhi - d/dt dG/dPhi_t = 0 form: lhs = delta G
dG_dPhi, dG_dPsi = EL[0].lhs, EL[1].lhs
# defect conditions of Lecture 3: Phi_x = Psi_t - B_Phi, Psi_x = Phi_t + B_Psi
ok_m = (sp.simplify(dG_dPhi - (Qf_.diff(t) - Bt.diff(Pf_))) == 0
        and sp.simplify(-dG_dPsi - (Pf_.diff(t) + Bt.diff(Qf_))) == 0)
report("(m) space as time: delta G/delta Phi = Phi_x and -delta G/delta Psi = Psi_x; -G generates "
       "(Phi, -Phi_x) -> (Psi, -Psi_x)", ok_m)

# ---------------------------------------------------------------------------------------------------------
# (n) exercises and remarks
Lplus = s * ex(-I * S / 4) + sig / lam * ex(-I * Dl / 4)
# type I_+ conditions: d_+ Sigma = -2 sigma m0 sin(Delta/2), d_- Delta = (2 m0/sigma) sin(Sigma/2)
solp = sp.solve([dpS + 2 * sig * m0 * sp.sin(Dl / 2), dmD - 2 * m0 / sig * sp.sin(S / 2)], [Px, Qx], dict=True)[0]
dLpt = Lplus.diff(P) * Pt + Lplus.diff(Q) * Qt
ok_n1 = zero_mat((dLpt - (Lplus * at_v(Q, Qx) - at_v(P, Px) * Lplus)).subs(solp))
# boost covariance of L_d: L_d(e^th lam; e^th sigma) = L_d(lam; sigma)
ok_n2 = zero_mat(Ld(l=sp.exp(thb) * lam, sg=sp.exp(thb) * sig) - Ld())
# the inverse transformation: L_d(sigma)(Phi, Psi) L_d(-sigma)(Psi, Phi) = det L_d
ok_n3 = zero_mat(Ld() * Ld(sg=-sig, F=Q, G=P) - (1 - sig**2 / lam**2) * E2)
# the two even vacua
ok_n4 = (zero_mat(Ld().subs({P: 0, Q: 0}) - (E2 + sig / lam * s))
         and zero_mat(Ld().subs({P: 2 * sp.pi, Q: 2 * sp.pi}) - (E2 - sig / lam * s)))
report("(n) type-I_+ matrix on its conditions; L_d depends on sigma/lambda only; "
       "L_d(sigma)(Phi, Psi) L_d(-sigma)(Psi, Phi) = det; the even vacua", ok_n1 and ok_n2 and ok_n3 and ok_n4)

# ---------------------------------------------------------------------------------------------------------
# (o) the a_r^(1) defect in components
def ar_lax_symbolic(hh):
    tt_, xx_ = sp.symbols('t x')
    ph = [sp.Function(f'ph{k}')(tt_, xx_) for k in range(hh)]
    dpl = lambda G: sp.diff(G, tt_) + sp.diff(G, xx_)
    dml = lambda G: sp.diff(G, tt_) - sp.diff(G, xx_)
    mm_, ll = sp.symbols('m lambda', positive=True)
    Ej = lambda i, j: sp.Matrix(hh, hh, lambda a, b: 1 if (a, b) == (i % hh, j % hh) else 0)
    ea = [sp.exp((ph[j] - ph[(j + 1) % hh]) / 2) for j in range(hh)]
    Ap_ = sp.diag(*[dpl(p_) / 2 for p_ in ph]) + mm_ * ll * sum((ea[j] * Ej(j, j + 1) for j in range(hh)), sp.zeros(hh))
    Am_ = -sp.diag(*[dml(p_) / 2 for p_ in ph]) - mm_ / ll * sum((ea[j] * Ej(j + 1, j) for j in range(hh)), sp.zeros(hh))
    Fc_ = dpl(Am_) - dml(Ap_) + Ap_ * Am_ - Am_ * Ap_
    toda = [dpl(dml(ph[k])) + mm_**2 * (sp.exp(ph[k] - ph[(k + 1) % hh]) - sp.exp(ph[k - 1] - ph[k])) for k in range(hh)]
    return all(sp.simplify(Fc_[i, j] - (-toda[i] if i == j else 0)) == 0 for i in range(hh) for j in range(hh))


rngo = np.random.default_rng(2)


def ar_time_eq(hh, sg, lv, scale=1.0, mv=1.0):
    phi = rngo.normal(size=hh)
    psi = rngo.normal(size=hh)
    pt = rngo.normal(size=hh)
    qt = rngo.normal(size=hh)
    phi, psi, pt, qt = phi - phi.mean(), psi - psi.mean(), pt - pt.mean(), qt - qt.mean()
    fo = np.exp(psi - np.roll(phi, -1))
    go = np.exp(phi - psi)
    A = np.zeros((2 * hh, 2 * hh))
    bv = np.zeros(2 * hh)
    for k in range(hh):
        A[k, k], A[k, hh + k] = 1, -1
        bv[k] = -mv * sg * (fo[k] - fo[k - 1]) - (pt[k] - qt[k])
        A[hh + k, hh + k], A[hh + k, (k + 1) % hh] = -1, 1
        bv[hh + k] = (mv / sg) * (go[k] - go[(k + 1) % hh]) - (qt[k] - pt[(k + 1) % hh])
    sol = np.linalg.lstsq(A, bv, rcond=None)[0]
    px, qx = sol[:hh] - sol[:hh].mean(), sol[hh:] - sol[hh:].mean()
    ok_sol = np.max(np.abs(A @ np.concatenate([px, qx]) - bv)) < 1e-12
    Eo = lambda i, j: np.outer(np.eye(hh)[i % hh], np.eye(hh)[j % hh])

    def at_ar(F, Fx):
        Ep = sum(np.exp((F[j] - F[(j + 1) % hh]) / 2) * Eo(j, j + 1) for j in range(hh))
        Em = sum(np.exp((F[j] - F[(j + 1) % hh]) / 2) * Eo(j + 1, j) for j in range(hh))
        return 0.5 * np.diag(Fx) + mv / 2 * (lv * Ep - Em / lv)

    def L_ar(Pv_, Qv_):
        ff, gg = np.exp(Qv_ - np.roll(Pv_, -1)), np.exp(Pv_ - Qv_)
        M = np.diag(np.sqrt(gg)).astype(complex)
        for j in range(hh):
            M[(j + 1) % hh, j] += scale * sg / lv * np.sqrt(ff[j])
        return M
    epsd = 1e-6
    Ltd = (L_ar(phi + epsd * pt, psi + epsd * qt) - L_ar(phi - epsd * pt, psi - epsd * qt)) / (2 * epsd)
    return ok_sol, np.max(np.abs(Ltd - (L_ar(phi, psi) @ at_ar(psi, qx) - at_ar(phi, px) @ L_ar(phi, psi))))


ok_o = ar_lax_symbolic(3)
for hh in (2, 3, 4, 5):
    oks, r1 = ar_time_eq(hh, 0.7, 1.3 + 0.4j)
    _, r2 = ar_time_eq(hh, 0.7, 1.3 + 0.4j, scale=1.2)
    ok_o = ok_o and oks and r1 < 1e-8 and r2 > 1e-2
report("(o) a_r^(1) in components: zero curvature is the Toda equation (h = 3); (L_d)_kk = g_k^1/2, "
       "(L_d)_{j+1,j} = (sigma/lambda) f_j^1/2 satisfies the time equation for h = 2..5; fails for 1.2 sigma", ok_o)

# ---------------------------------------------------------------------------------------------------------
# (p) the block identities of the freezing calculation, A_pm(-Psi) = s A_pm(Psi) s, the book's form of L_d
Dd, Ddi = ex(I * Dl / 4), ex(-I * Dl / 4)
Od, Obd = ex(-I * S / 4) * s, s * ex(-I * S / 4)
blocks = [
    Dd * Xm(Q) - Xm(P) * Dd,
    Obd * Xm(Q) - Xm(P) * Obd,
    Ddi * Ym(Q) - Ym(P) * Ddi,
    Od * Ym(Q) - Ym(P) * Od,
    Od * Xm(Q) - Xm(P) * Od + 2 * I * sp.sin(S / 2) * h * Dd,
    Ddi * Xm(Q) - Xm(P) * Ddi + 2 * I * sp.sin(Dl / 2) * h * Obd,
    Dd * Ym(Q) - Ym(P) * Dd - 2 * I * sp.sin(Dl / 2) * h * Od,
    Obd * Ym(Q) - Ym(P) * Obd - 2 * I * sp.sin(S / 2) * h * Ddi,
    Obd - ex(I * S / 4) * s,
]
ok_p1 = all(zero_mat(B_) for B_ in blocks)
ok_p2 = not zero_mat(Od * Xm(Q) - Xm(P) * Od - 2 * I * sp.sin(S / 2) * h * Dd)      # control: wrong sign
ok_p3 = (zero_mat(Apv(-Q, -Qp) - s * Apv(Q, Qp) * s) and zero_mat(Amv(-Q, -Qm) - s * Amv(Q, Qm) * s)
         and zero_mat(s * ex(-I * S / 4) + sig / lam * ex(-I * Dl / 4) - Ld(G=-Q) * s))
ok_p4 = zero_mat(Ld() - ex(I * Dl / 4) * (E2 + sig / lam * Ym(P)))
report("(p) the eight block identities of the freezing calculation (a wrong sign fails); "
       "A_pm(-Psi) = s A_pm(Psi) s and L_d^+ = L_d(Phi, -Psi) s; L_d = e^{i Delta h/4}(1 + (sigma/lambda) s_-Phi)",
       ok_p1 and ok_p2 and ok_p3 and ok_p4)

# ---------------------------------------------------------------------------------------------------------
# (q) the spin -+3 charges of a soliton and the charges stored by the odd jump
uq = sp.Symbol('u', real=True)
thq = sp.Symbol('theta_q', real=True)
Fq = 4 * sp.atan(sp.exp(uq))                       # soliton, u = m0 (x cosh theta - t sinh theta)


def q_soliton(sign):
    # d_pm u = m0 e^{-+theta} (sign = +1 for d_+, -1 for d_-); fields depend on u only
    k = m0 * sp.exp(-sign * thq) * sign
    d1 = k * sp.diff(Fq, uq)
    d2 = k**2 * sp.diff(Fq, uq, 2)
    dens = d2**2 - d1**4 / 4 + m0**2 * sp.cos(Fq) * d1**2
    return sp.lambdify((uq, thq), (dens / (8 * m0**3) / (m0 * sp.cosh(thq))).subs(m0, 1), 'numpy')


ok_q = True
for sign in (1, -1):
    fq = q_soliton(sign)
    for thv_ in (-0.7, 0.0, 0.45, 1.3):
        Qv = quad(lambda uu: fq(uu, thv_), -40, 40, limit=400)[0]
        ok_q = ok_q and abs(Qv + 2 / 3 * math.exp(-3 * sign * thv_)) < 1e-9
        if thv_ != 0.0:          # control: the wrong spin
            ok_q = ok_q and abs(Qv + 2 / 3 * math.exp(-sign * thv_)) > 1e-3
for sgv in (0.4, 1.0, 2.3):
    Dm3 = lambda Sg_: sgv**3 * (math.cos(Sg_ / 2) - 2 / 3 * math.cos(Sg_ / 2)**3)
    Dp3 = lambda Dg_: sgv**-3 * (math.cos(Dg_ / 2) - 2 / 3 * math.cos(Dg_ / 2)**3)
    ok_q = ok_q and abs(Dm3(2 * math.pi) - Dm3(0.0) + 2 / 3 * sgv**3) < 1e-12
    ok_q = ok_q and abs(Dp3(-2 * math.pi) - Dp3(0.0) + 2 / 3 * sgv**-3) < 1e-12
    ok_q = ok_q and abs(Dm3(0.0) - sgv**3 / 3) < 1e-12 and abs(Dp3(0.0) - sgv**-3 / 3) < 1e-12
report("(q) a soliton of rapidity theta carries the spin -+3 charges -(2/3) e^{-+3 theta} (quadrature; the wrong spin "
       "fails); the odd jump (0, 2 pi) stores -(2/3) sigma^{+-3} relative to the even vacuum, where D_-+3 = sigma^{+-3}/3",
       ok_q)

# ---------------------------------------------------------------------------------------------------------
# (r) the a_r^(1) matrix is a gauge transformation in both light-cone directions; the book's form of L_d
rngr = np.random.default_rng(7)


def ar_gauge(hh, scale=1.0, mv=0.9, sg=0.7, lv=0.8 + 0.6j):
    hyp = lambda v: v - v.mean()
    phi, psi = hyp(rngr.normal(size=hh)), hyp(rngr.normal(size=hh))
    pp_phi, pm_phi = hyp(rngr.normal(size=hh)), hyp(rngr.normal(size=hh))
    ff = np.exp(psi - np.roll(phi, -1))
    gg = np.exp(phi - psi)
    pp_psi = pp_phi + mv * sg * (ff - np.roll(ff, 1))                    # d_+(phi_k - psi_k) = -m s (f_k - f_{k-1})
    pm_psi = np.roll(pm_phi, -1) + mv / sg * (gg - np.roll(gg, -1))      # d_-(psi_k - phi_{k+1}) = (m/s)(g_k - g_{k+1})

    def Apl(F, dF):
        M = np.diag(dF / 2).astype(complex)
        for j in range(hh):
            M[j, (j + 1) % hh] += mv * lv * np.exp((F[j] - F[(j + 1) % hh]) / 2)
        return M

    def Amn(F, dF):
        M = -np.diag(dF / 2).astype(complex)
        for j in range(hh):
            M[(j + 1) % hh, j] += -mv / lv * np.exp((F[j] - F[(j + 1) % hh]) / 2)
        return M

    def Lar(F, G):
        M = np.diag(np.exp((F - G) / 2)).astype(complex)
        for j in range(hh):
            M[(j + 1) % hh, j] += scale * sg / lv * np.exp((G[j] - F[(j + 1) % hh]) / 2)
        return M

    def dL(dF, dG, eps=1e-6):
        return (Lar(phi + eps * dF, psi + eps * dG) - Lar(phi - eps * dF, psi - eps * dG)) / (2 * eps)
    rp = dL(pp_phi, pp_psi) - (Lar(phi, psi) @ Apl(psi, pp_psi) - Apl(phi, pp_phi) @ Lar(phi, psi))
    rm = dL(pm_phi, pm_psi) - (Lar(phi, psi) @ Amn(psi, pm_psi) - Amn(phi, pm_phi) @ Lar(phi, psi))
    return max(np.abs(rp).max(), np.abs(rm).max())


def ar_w(hh, transpose=False):
    # orthonormal basis of the hyperplane, pi with (pi v)_k = v_{k+1}, D = 2 (1 - pi)^{-1}, E = 1 - D, w = D^T phi - E psi
    Bh = np.linalg.qr(np.vstack([np.ones(hh), rngr.normal(size=(hh - 1, hh))]).T)[0][:, 1:]
    piM = np.roll(np.eye(hh), 1, axis=1)
    Dh = 2 * np.linalg.inv(np.eye(hh - 1) - Bh.T @ piM @ Bh)
    if transpose:
        Dh = Dh.T
    Eh = np.eye(hh - 1) - Dh
    phi, psi = Bh @ rngr.normal(size=hh - 1), Bh @ rngr.normal(size=hh - 1)
    w = Bh @ (Dh.T @ (Bh.T @ phi) - Eh @ (Bh.T @ psi))
    lhs = np.array([w[j] - w[(j + 1) % hh] for j in range(hh)])
    rhs = np.array([psi[j] + psi[(j + 1) % hh] - 2 * phi[(j + 1) % hh] for j in range(hh)])
    return np.abs(lhs - rhs).max()


ok_r = True
for hh in (2, 3, 4, 5):
    ok_r = ok_r and ar_gauge(hh) < 1e-8 and ar_gauge(hh, scale=1.2) > 1e-2
    ok_r = ok_r and ar_w(hh) < 1e-12
    if hh > 2:
        ok_r = ok_r and ar_w(hh, transpose=True) > 1e-3
report("(r) a_r^(1): d_pm L_d = L_d A_pm(psi) - A_pm(phi) L_d separately on the light-cone conditions, h = 2..5 "
       "(fails for 1.2 sigma); alpha_j . w = psi_j + psi_{j+1} - 2 phi_{j+1} for w = D^T phi - E psi (fails for D^T)",
       ok_r)

raise SystemExit(0 if OK else 1)
