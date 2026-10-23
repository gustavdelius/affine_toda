"""Checks for Lecture 3 of the defect lecture notes, "Freezing it at a point"
(lecture-notes/defects/03-freezing.qmd).

Conventions of the lecture: Phi on x<0, Psi on x>0, both dimensionless (Phi = beta_SG phi); the action is
    S = (1/beta^2) int dt { int_{-inf}^0 l(Phi) dx + int_0^inf l(Psi) dx + [ D Phi Psi_t - B(Phi,Psi) ]_{x=0} },
    l(Phi) = Phi_t^2/2 - Phi_x^2/2 - V(Phi),   sine-Gordon V = m0^2 (1 - cos Phi);
energy E = (1/beta^2) int (Phi_t^2/2 + Phi_x^2/2 + V), momentum P = -(1/beta^2) int Phi_t Phi_x (a soliton has
M sinh theta, M = 8 m0/beta^2); d_pm = d_t pm d_x; s = Phi + Psi, d = Phi - Psi. Baecklund transformation of
Lecture 2:
    d_+ (Phi - Psi) = 2 a m0 sin((Phi+Psi)/2),   d_- (Phi + Psi) = -(2 m0/a) sin((Phi-Psi)/2).
Written from these definitions, independently of code/defect_classical_sg.py.
One PASS or FAIL line per check; exit status 1 on any FAIL.

  (a) Varying the action: bulk boundary terms -Phi_x dPhi (left) and +Psi_x dPsi (right); a kinetic term
      A(Phi,Psi) Phi_t + C(Phi,Psi) Psi_t enters the conditions only through w = C_Phi - A_Psi:
      Phi_x = w Psi_t - B_Phi, Psi_x = w Phi_t + B_Psi.
  (b) The frozen Baecklund transformation is the defect condition with D = 1 and
      B = 2 m0 (a cos(s/2) + cos(d/2)/a), a = -sigma: beta^2 times the book's defect potential, up to a constant.
  (c) Energy: Phi_t Phi_x - Psi_t Psi_x = -dB/dt on the defect conditions, for any w and any B.
  (d) Momentum flux: dP/dt = [T^xx]_{0+} - [T^xx]_{0-}; on the conditions it is
      (1-w^2)(Psi_t^2-Phi_t^2)/2 + w(Phi_t B_Psi + Psi_t B_Phi) + B_Psi^2/2 - B_Phi^2/2 + V(Phi) - V(Psi).
  (e) A continuous field with a point potential U: dP/dt = U' Phi_x(0-) + U'^2/2, not a total derivative;
      for the delta-defect of Lecture 1 the time-averaged momentum fluxes differ by -k^2 |R|^2.
  (f) d'Alembert in field space: B = f(s) + g(d) gives P_d = D (g - f) with the required gradient and
      B_Phi^2/2 - B_Psi^2/2 = 2 f'(s) g'(d); (d_s^2 - d_d^2) [V(Phi) - V(Psi)] = 0 for any V.
  (g) Solving 2 f'(s) g'(d) = V(Phi) - V(Psi) for sine-Gordon: V(Phi) - V(Psi) = 2 m0^2 sin(s/2) sin(d/2); with
      f = A cos(s/2) + A' sin(s/2), g = C cos(d/2) + C' sin(d/2) the only solutions are A' = C' = 0, A C = 4 m0^2,
      the frozen Baecklund potential; V(Phi) - V(Psi) is not of product form for double sine-Gordon or phi^4.
  (h) The exponential solutions reproduce Bowcock-Corrigan-Zambon's list (their 2003 eqs. (2.10)-(2.11)),
      and the free massive case their 2004 eq. (2.7).
  (i) Unfreezing: imposed for all x, the defect conditions with D = 1 give
      Phi_tt - Phi_xx = (B_PhiPhi - B_PsiPsi) Psi_t - d_Phi W, Psi_tt - Psi_xx = (B_PhiPhi - B_PsiPsi) Phi_t + d_Psi W,
      W = (B_Phi^2 - B_Psi^2)/2; so the momentum conditions make them a Baecklund transformation
      (sine-Gordon, and Liouville next to the free wave equation).
  (j) Conservation along exact solutions: E + B/beta^2 and P + P_d/beta^2 are constant on the emission solution
      (0, B_{-sigma} 0) and on the Bianchi pair (Phi_1, Phi_12), whose defect residual is checked too.
  (k) The light-cone split: B + P_d = -4 m0 (cos(d/2) - 1)/sigma and B - P_d = -4 m0 sigma (cos(s/2) - 1);
      a static jump of 2 pi stores M (cosh theta_d, sinh theta_d); even jumps are minima of B, odd ones maxima.
  (l) No reflection from the fluxes: the time-averaged energy and momentum fluxes of an incident plus reflected
      wave are -omega k (1-|R|^2)/2 and k^2 (1+|R|^2)/2; the linearized frozen Baecklund conditions give R = 0,
      |T| = 1, and the book's transmission factor (eq-shg-linear-transmission).
  (m) Type I_+: Psi -> -Psi maps the D = 1 defect to the D = -1 defect with sigma and 1/sigma exchanged.
  (n) The massless limit: at m0 = 0 the conditions say Psi_L = Phi_L, Psi_R = -Phi_R; in the scaling limit
      m0 -> 0 with mu = m0/sigma fixed a right-mover is transmitted with (mu + 2 i omega)/(mu - 2 i omega),
      the limit of the massive transmission factor.
  (o) Spin 3: the frozen Baecklund potential satisfies the conditions (2.6)-(2.9) of Bowcock-Corrigan-Zambon 2003
      with lambda^2 = -1/4; for massless free fields any chiral potential g(Phi - Psi) conserves momentum, while
      their (2.9) admits exponentials but not (Phi - Psi)^4.
  (p) With such a g the right-movers obey u' + w' = g'(u - w - c), and Q_4 = int (d_- field)^4 changes at the
      rate -8 g' X' (g'^2 + X'^2) (exercise).
  (q) Exercises: the quadratic part of the frozen Baecklund potential is the free defect with lambda = sigma;
      sigma -> -sigma is Psi -> Psi + 2 pi; B >= 0 for sigma > 0; freezing Psi = psi_0 gives a boundary potential
      -M cos((Phi - Phi_0)/2) + const, and every (M > 0, Phi_0) arises.

Runtime: about ten seconds.
"""
import numpy as np
import sympy as sp
from scipy.integrate import quad

OK = True


def report(name, ok):
    global OK
    OK = OK and bool(ok)
    print(("PASS " if ok else "FAIL ") + name)


rng = np.random.default_rng(3)
P, Q = sp.symbols('Phi Psi')                      # defect values Phi(0,t), Psi(0,t)
Pt, Qt, Px, Qx = sp.symbols('Phi_t Psi_t Phi_x Psi_x')
m0 = sp.Symbol('m0', positive=True)
sig = sp.Symbol('sigma', nonzero=True)
w = sp.Symbol('w')
s, d = sp.symbols('s d')
Bf = sp.Function('B')(P, Q)
Vf = sp.Function('V')
ff, gf = sp.Function('f'), sp.Function('g')
tt, xx = sp.symbols('t x', real=True)
A_, B_ = (P + Q) / 2, (P - Q) / 2                 # s/2 and d/2


def V_sg(f):
    return m0**2 * (1 - sp.cos(f))


def zero(expr):
    return sp.simplify(sp.expand(expr)) == 0


# ---------------------------------------------------------------------------------------------------------
# (a) varying the action
X = sp.Symbol('X')
Phi_test = sp.sin(X) + X**2 / 3
eta = sp.exp(X) * (1 + X)                       # test variation, nonzero at x = 0
L_ = 2
# left: int_{-L}^0 -Phi_x eta_x = -[Phi_x eta]_{0} + [Phi_x eta]_{-L} + int Phi_xx eta
lhs = sp.integrate(-sp.diff(Phi_test, X) * sp.diff(eta, X), (X, -L_, 0))
rhs = -(sp.diff(Phi_test, X) * eta).subs(X, 0) + (sp.diff(Phi_test, X) * eta).subs(X, -L_) \
    + sp.integrate(sp.diff(Phi_test, X, 2) * eta, (X, -L_, 0))
# right: int_0^L -Psi_x eta_x = +[Psi_x eta]_{0} - [Psi_x eta]_{L} + int Psi_xx eta
lhs_r = sp.integrate(-sp.diff(Phi_test, X) * sp.diff(eta, X), (X, 0, L_))
rhs_r = (sp.diff(Phi_test, X) * eta).subs(X, 0) - (sp.diff(Phi_test, X) * eta).subs(X, L_) \
    + sp.integrate(sp.diff(Phi_test, X, 2) * eta, (X, 0, L_))
ok_bulk = abs(sp.N(lhs - rhs)) < 1e-12 and abs(sp.N(lhs_r - rhs_r)) < 1e-12

T = sp.Symbol('t')
p_t, q_t = sp.Function('p')(T), sp.Function('q')(T)
Af, Cf, Bt = sp.Function('A'), sp.Function('C'), sp.Function('B')
Ld = Af(p_t, q_t) * sp.diff(p_t, T) + Cf(p_t, q_t) * sp.diff(q_t, T) - Bt(p_t, q_t)


def el(L, fn):
    """d L/d fn - d/dt (d L/d fn_t): the coefficient of delta fn in the variation of int L dt"""
    return sp.diff(L, fn) - sp.diff(sp.diff(L, sp.diff(fn, T)), T)


wpq = sp.diff(Cf(p_t, q_t), p_t) - sp.diff(Af(p_t, q_t), q_t)
# coefficient of dPhi: -Phi_x + el_p = 0 ; coefficient of dPsi: +Psi_x + el_q = 0
ok_a1 = zero(el(Ld, p_t) - (wpq * sp.diff(q_t, T) - sp.diff(Bt(p_t, q_t), p_t)))
ok_a2 = zero(-el(Ld, q_t) - (wpq * sp.diff(p_t, T) + sp.diff(Bt(p_t, q_t), q_t)))
report("(a) bulk boundary terms -Phi_x dPhi and +Psi_x dPsi; the kinetic term enters only through "
       "w = C_Phi - A_Psi: Phi_x = w Psi_t - B_Phi, Psi_x = w Phi_t + B_Psi", ok_bulk and ok_a1 and ok_a2)

# ---------------------------------------------------------------------------------------------------------
# (b) the frozen Baecklund transformation
a = sp.Symbol('a', nonzero=True)
Bbt = 2 * m0 * (a * sp.cos(A_) + sp.cos(B_) / a)
cond_x = {Px: Qt - sp.diff(Bbt, P), Qx: Pt + sp.diff(Bbt, Q)}
bt1 = (Pt + Px) - (Qt + Qx) - 2 * a * m0 * sp.sin(A_)          # d_+(Phi - Psi) - 2 a m0 sin(s/2)
bt2 = (Pt - Px) + (Qt - Qx) + 2 * m0 / a * sp.sin(B_)          # d_-(Phi + Psi) + (2 m0/a) sin(d/2)
ok_b1 = zero(bt1.subs(cond_x)) and zero(bt2.subs(cond_x))
book_Bd = -2 * m0 * (sig * (sp.cos(A_) - 1) + (sp.cos(B_) - 1) / sig)     # beta^2 x eq-sg-defect-potential
ok_b2 = zero(Bbt.subs(a, -sig) - book_Bd + 2 * m0 * (sig + 1 / sig))
report("(b) the Baecklund equations at x = 0 are the defect conditions with D = 1 and "
       "B = 2 m0 (a cos(s/2) + cos(d/2)/a); with a = -sigma this is beta^2 times the book's defect potential",
       ok_b1 and ok_b2)

# ---------------------------------------------------------------------------------------------------------
# (c) energy is always conserved
Bp, Bq = sp.diff(Bf, P), sp.diff(Bf, Q)
gen_x = {Px: w * Qt - Bp, Qx: w * Pt + Bq}
report("(c) dE/dt = [Phi_t Phi_x]_{0-} - [Psi_t Psi_x]_{0+} = -dB/dt for every w and every B",
       zero((Pt * Px - Qt * Qx).subs(gen_x) + Bp * Pt + Bq * Qt))

# ---------------------------------------------------------------------------------------------------------
# (d) the momentum flux
ph = sp.Function('phi')(tt, xx)
ddt = sp.diff(-sp.diff(ph, tt) * sp.diff(ph, xx), tt).subs(sp.diff(ph, tt, 2), sp.diff(ph, xx, 2) - sp.diff(Vf(ph), ph))
Txx = sp.diff(ph, tt)**2 / 2 + sp.diff(ph, xx)**2 / 2 - Vf(ph)
ok_d1 = zero(ddt + sp.diff(Txx, xx))
flux_P = (Qt**2 / 2 + Qx**2 / 2 - Vf(Q) - Pt**2 / 2 - Px**2 / 2 + Vf(P)).subs(gen_x)
claimed = (1 - w**2) * (Qt**2 - Pt**2) / 2 + w * (Pt * Bq + Qt * Bp) + Bq**2 / 2 - Bp**2 / 2 + Vf(P) - Vf(Q)
report("(d) d_t(-Phi_t Phi_x) = -d_x T^xx in the bulk; on the conditions dP/dt = (1-w^2)(Psi_t^2-Phi_t^2)/2 "
       "+ w(Phi_t B_Psi + Psi_t B_Phi) + B_Psi^2/2 - B_Phi^2/2 + V(Phi) - V(Psi)", ok_d1 and zero(flux_P - claimed))

# ---------------------------------------------------------------------------------------------------------
# (e) a continuous field with a point potential U
Uf = sp.Function('U')(P)
Pxm = sp.Symbol('Phi_x_minus')
Pxp = Pxm + sp.diff(Uf, P)                     # Phi_x(0+) - Phi_x(0-) = U'(Phi(0)), Phi continuous
flux_cont = (Pt**2 / 2 + Pxp**2 / 2 - Vf(P)) - (Pt**2 / 2 + Pxm**2 / 2 - Vf(P))
ok_e1 = zero(flux_cont - (sp.diff(Uf, P) * Pxm + sp.diff(Uf, P)**2 / 2))
mm, gg, kk = 1.3, 0.7, 0.9
om = np.sqrt(kk**2 + mm**2)
Rd = gg / (2j * kk - gg)
Td = 2j * kk / (2j * kk - gg)


def avg_flux(amp_in, amp_out, k, m, x0):
    """time averages of Phi_t Phi_x and of T^xx = Phi_t^2/2 + Phi_x^2/2 - m^2 Phi^2/2 at x0 for
    Phi = Re[(amp_in e^{ikx} + amp_out e^{-ikx}) e^{-i omega t}]"""
    omg_ = np.sqrt(k**2 + m**2)
    ts = np.linspace(0, 2 * np.pi / omg_, 6001)[:-1]
    u = amp_in * np.exp(1j * k * x0) + amp_out * np.exp(-1j * k * x0)
    ux = 1j * k * (amp_in * np.exp(1j * k * x0) - amp_out * np.exp(-1j * k * x0))
    f_ = np.real(u * np.exp(-1j * omg_ * ts))
    ft_ = np.real(-1j * omg_ * u * np.exp(-1j * omg_ * ts))
    fx_ = np.real(ux * np.exp(-1j * omg_ * ts))
    return np.mean(ft_ * fx_), np.mean(ft_**2 / 2 + fx_**2 / 2 - m**2 * f_**2 / 2)


diff_flux = avg_flux(Td, 0, kk, mm, 0.0)[1] - avg_flux(1, Rd, kk, mm, 0.0)[1]
ok_e2 = abs(diff_flux + kk**2 * abs(Rd)**2) < 1e-10
report("(e) continuous field with a point potential U: dP/dt = U' Phi_x(0-) + U'^2/2, not a total derivative; "
       "delta-defect: averaged momentum fluxes differ by -k^2 |R|^2", ok_e1 and ok_e2)

# ---------------------------------------------------------------------------------------------------------
# (f) d'Alembert in field space
Bdal = ff(P + Q) + gf(P - Q)
ok_f1 = zero(sp.diff(Bdal, P, 2) - sp.diff(Bdal, Q, 2))
ok_f2 = True
for Dv in (1, -1):
    Pd_ = Dv * (gf(P - Q) - ff(P + Q))
    # needed: D (Phi_t B_Psi + Psi_t B_Phi) = -dP_d/dt, i.e. d_Phi P_d = -D B_Psi and d_Psi P_d = -D B_Phi
    ok_f2 = ok_f2 and zero(sp.diff(Pd_, P) + Dv * sp.diff(Bdal, Q)) and zero(sp.diff(Pd_, Q) + Dv * sp.diff(Bdal, P))
ok_f3 = zero(sp.diff(Bdal, P)**2 / 2 - sp.diff(Bdal, Q)**2 / 2
             - 2 * sp.diff(ff(P + Q), P) * sp.diff(gf(P - Q), P))
rhs_fe = Vf((s + d) / 2) - Vf((s - d) / 2)
ok_f4 = zero(sp.diff(rhs_fe, s, 2) - sp.diff(rhs_fe, d, 2))
report("(f) B = f(s) + g(d) solves B_PhiPhi = B_PsiPsi, P_d = D(g - f) has the required gradient for D = +-1, "
       "B_Phi^2/2 - B_Psi^2/2 = 2 f'(s) g'(d), and (d_s^2 - d_d^2)[V(Phi) - V(Psi)] = 0 for any V",
       ok_f1 and ok_f2 and ok_f3 and ok_f4)

# ---------------------------------------------------------------------------------------------------------
# (g) solving the functional equation for sine-Gordon
rhs_sg = sp.expand_trig(sp.expand(V_sg((s + d) / 2) - V_sg((s - d) / 2)))
ok_g0 = zero(rhs_sg - 2 * m0**2 * sp.sin(s / 2) * sp.sin(d / 2))
A0, A1, C0, C1 = sp.symbols('A Ap C Cp')
fs = A0 * sp.cos(s / 2) + A1 * sp.sin(s / 2)
gd = C0 * sp.cos(d / 2) + C1 * sp.sin(d / 2)
ss_, cs_, sd_, cd_ = sp.symbols('ss cs sd cd')
diffg = sp.expand(2 * sp.diff(fs, s) * sp.diff(gd, d) - 2 * m0**2 * sp.sin(s / 2) * sp.sin(d / 2))
diffg = diffg.subs({sp.sin(s / 2): ss_, sp.cos(s / 2): cs_, sp.sin(d / 2): sd_, sp.cos(d / 2): cd_})
# the four products ss*sd, ss*cd, cs*sd, cs*cd are linearly independent functions of (s, d)
eqs = sp.Poly(diffg, ss_, cs_, sd_, cd_).coeffs()
sols = sp.solve(eqs, [A1, C1, C0], dict=True)
ok_g1 = len(sols) == 1 and sols[0][A1] == 0 and sols[0][C1] == 0 and zero(sols[0][C0] - 4 * m0**2 / A0)
B_sol = (-2 * m0 * sig * sp.cos(s / 2) - 2 * m0 / sig * sp.cos(d / 2)).subs({s: P + Q, d: P - Q})
ok_g2 = zero(B_sol - Bbt.subs(a, -sig))


def rank_ratio(Vexpr, var, npts=7):
    """second singular value / first of the matrix R(s_i, d_j), R = V((s+d)/2) - V((s-d)/2):
    zero exactly when R(s, d) = F(s) G(d)"""
    Rfun = sp.lambdify((s, d), Vexpr.subs(var, (s + d) / 2) - Vexpr.subs(var, (s - d) / 2), 'numpy')
    S_, D_ = rng.normal(size=npts), rng.normal(size=npts)
    Mx = np.array([[Rfun(si, dj) for dj in D_] for si in S_], dtype=float)
    sv = np.linalg.svd(Mx, compute_uv=False)
    return sv[1] / sv[0]


yv = sp.Symbol('y')
r_sg = rank_ratio(1 - sp.cos(yv), yv)
r_shg = rank_ratio(sp.cosh(yv) - 1, yv)
r_dsg = rank_ratio((1 - sp.cos(yv)) + sp.Rational(3, 10) * (1 - sp.cos(2 * yv)), yv)
r_phi4 = rank_ratio(yv**2 / 2 + yv**4 / 4, yv)
ok_g3 = r_sg < 1e-12 and r_shg < 1e-12 and r_dsg > 1e-3 and r_phi4 > 1e-3
report("(g) sine-Gordon: V(Phi)-V(Psi) = 2 m0^2 sin(s/2) sin(d/2); the functional equation forces A' = C' = 0 and "
       "A C = 4 m0^2, the frozen Baecklund potential with a = -sigma; product form holds for sine- and sinh-Gordon "
       f"(rank ratios {r_sg:.0e}, {r_shg:.0e}) but not for double sine-Gordon ({r_dsg:.1e}) or phi^4 ({r_phi4:.1e})",
       ok_g0 and ok_g1 and ok_g2 and ok_g3)

# ---------------------------------------------------------------------------------------------------------
# (h) the general exponential solution and BCZ's list; the free massive case
al, alp, ga, gap, ka = sp.symbols('alpha alphap gamma gammap kappa', nonzero=True)
F_ = al * sp.exp(ka * s) + alp * sp.exp(-ka * s)
G_ = ga * sp.exp(ka * d) + gap * sp.exp(-ka * d)
V1 = 2 * al * ga * sp.exp(2 * ka * P) + 2 * alp * gap * sp.exp(-2 * ka * P)
V2 = -2 * al * gap * sp.exp(2 * ka * Q) - 2 * alp * ga * sp.exp(-2 * ka * Q)
ok_h1 = zero(sp.expand(2 * F_ * G_).subs({s: P + Q, d: P - Q}) - (V1 - V2))
# BCZ 2003 (2.10)-(2.11): B = a e^{l(f1+f2)} + b e^{l(f1-f2)} + c e^{-l(f1-f2)} + d e^{-l(f1+f2)},
# V_k = A_k e^{2 l f_k} + B_k e^{-2 l f_k}, A1 = 2 l^2 a b, A2 = 2 l^2 a c, B1 = 2 l^2 c d, B2 = 2 l^2 b d
fB = al / ka * sp.exp(ka * s) - alp / ka * sp.exp(-ka * s)           # f with f' = F
gB = ga / ka * sp.exp(ka * d) - gap / ka * sp.exp(-ka * d)           # g with g' = G
ok_h2 = zero(sp.diff(fB, s) - F_) and zero(sp.diff(gB, d) - G_)
aB, bB, cB, dB = al / ka, ga / ka, -gap / ka, -alp / ka
ok_h2 = ok_h2 and zero((fB + gB).subs({s: P + Q, d: P - Q}) - (aB * sp.exp(ka * (P + Q)) + bB * sp.exp(ka * (P - Q))
                                                               + cB * sp.exp(-ka * (P - Q)) + dB * sp.exp(-ka * (P + Q))))
ok_h2 = ok_h2 and all(zero(l_ - r_) for l_, r_ in [
    (2 * ka**2 * aB * bB, V1.coeff(sp.exp(2 * ka * P))), (2 * ka**2 * aB * cB, V2.coeff(sp.exp(2 * ka * Q))),
    (2 * ka**2 * cB * dB, V1.coeff(sp.exp(-2 * ka * P))), (2 * ka**2 * bB * dB, V2.coeff(sp.exp(-2 * ka * Q)))])
# free massive (kappa = 0): BCZ 2004 eq. (2.7): B = (m lam/4)(phi+psi)^2 + (m/(4 lam))(phi-psi)^2
mfr, lam = sp.symbols('m lambda', positive=True)
Bfree = mfr * lam / 4 * (P + Q)**2 + mfr / (4 * lam) * (P - Q)**2
ok_h3 = zero(sp.diff(Bfree, P)**2 / 2 - sp.diff(Bfree, Q)**2 / 2 - (mfr**2 * P**2 / 2 - mfr**2 * Q**2 / 2)) \
    and zero(sp.diff(Bfree, P, 2) - sp.diff(Bfree, Q, 2))
report("(h) exponential solutions: 2 F G = V1(Phi) - V2(Psi) with BCZ 2003 (2.10)-(2.11) coefficients; the free "
       "massive defect of BCZ 2004 (2.7) satisfies the momentum conditions", ok_h1 and ok_h2 and ok_h3)

# ---------------------------------------------------------------------------------------------------------
# (i) unfreezing: the defect conditions imposed everywhere
phL, psL = sp.Function('Phi')(tt, xx), sp.Function('Psi')(tt, xx)


def unfreeze(Bexpr):
    """Phi_tt - Phi_xx and Psi_tt - Psi_xx implied by Phi_x = Psi_t - B_Phi, Psi_x = Phi_t + B_Psi for all x"""
    Bp_ = sp.diff(Bexpr, P).subs({P: phL, Q: psL})
    Bq_ = sp.diff(Bexpr, Q).subs({P: phL, Q: psL})
    sysx = {sp.diff(phL, xx): sp.diff(psL, tt) - Bp_, sp.diff(psL, xx): sp.diff(phL, tt) + Bq_}
    Phixx = sp.diff(sp.diff(phL, tt) + Bq_, tt) - sp.diff(Bp_, xx)      # d_x(Psi_t - B_Phi) = d_t Psi_x - d_x B_Phi
    Psixx = sp.diff(sp.diff(psL, tt) - Bp_, tt) + sp.diff(Bq_, xx)      # d_x(Phi_t + B_Psi) = d_t Phi_x + d_x B_Psi
    e1 = sp.expand((sp.diff(phL, tt, 2) - Phixx).subs(sysx))
    e2 = sp.expand((sp.diff(psL, tt, 2) - Psixx).subs(sysx))
    return e1, e2


e1, e2 = unfreeze(Bf)
Wf = (Bp**2 - Bq**2) / 2
back = {P: phL, Q: psL}
gap_ = (sp.diff(Bf, P, 2) - sp.diff(Bf, Q, 2))
ok_i1 = zero(e1 - ((gap_ * sp.Symbol('PSIT')) - sp.diff(Wf, P)).subs(back).subs(sp.Symbol('PSIT'), sp.diff(psL, tt))) \
    and zero(e2 - ((gap_ * sp.Symbol('PHIT')) + sp.diff(Wf, Q)).subs(back).subs(sp.Symbol('PHIT'), sp.diff(phL, tt)))
e1s, e2s = unfreeze(Bbt.subs(a, -sig))
ok_i2 = zero(sp.expand_trig(e1s + m0**2 * sp.sin(phL))) and zero(sp.expand_trig(e2s + m0**2 * sp.sin(psL)))
BL = (al * sp.exp(ka * (P + Q)) + ga * sp.exp(ka * (P - Q))) / ka      # alpha' = gamma' = 0
e1l, e2l = unfreeze(BL)
ok_i3 = zero(e1l + sp.diff(2 * al * ga * sp.exp(2 * ka * phL), phL)) and zero(e2l)
report("(i) imposed for all x, the conditions give Phi_tt - Phi_xx = (B_PhiPhi - B_PsiPsi) Psi_t - d_Phi W and "
       "Psi_tt - Psi_xx = (B_PhiPhi - B_PsiPsi) Phi_t + d_Psi W; sine-Gordon on both sides for the frozen Baecklund "
       "potential, Liouville next to the wave equation for the Liouville-free potential", ok_i1 and ok_i2 and ok_i3)

# ---------------------------------------------------------------------------------------------------------
# (j) conservation along exact solutions (numerical quadrature, analytic derivatives)
mv = 1.0
Mv = 8 * mv                                   # beta^2 = 1
thd = 0.6
sg_v = np.exp(-thd)


def kink_expr(theta, c):
    return sp.exp(mv * (xx * sp.cosh(theta) - tt * sp.sinh(theta)) + c)


def lamb(expr):
    return (sp.lambdify((tt, xx), expr, 'numpy'), sp.lambdify((tt, xx), sp.diff(expr, tt), 'numpy'),
            sp.lambdify((tt, xx), sp.diff(expr, xx), 'numpy'))


def bulk_EP(F3, tv, x0, x1):
    F, Ft, Fx = F3
    e = lambda xv: Ft(tv, xv)**2 / 2 + Fx(tv, xv)**2 / 2 + mv**2 * (1 - np.cos(F(tv, xv)))
    p = lambda xv: -Ft(tv, xv) * Fx(tv, xv)
    return (quad(e, x0, x1, limit=500, epsabs=1e-13, epsrel=1e-13)[0],
            quad(p, x0, x1, limit=500, epsabs=1e-13, epsrel=1e-13)[0])


def defect_BP(ph0, ps0, sgv):
    Bv = -2 * mv * (sgv * (np.cos((ph0 + ps0) / 2) - 1) + (np.cos((ph0 - ps0) / 2) - 1) / sgv)
    Pv = 2 * mv * (sgv * (np.cos((ph0 + ps0) / 2) - 1) - (np.cos((ph0 - ps0) / 2) - 1) / sgv)
    return Bv, Pv


def totals(L3, R3, sgv, tv, Lbox=45.0):
    El, Pl = bulk_EP(L3, tv, -Lbox, 0.0)
    Er, Pr = bulk_EP(R3, tv, 0.0, Lbox)
    Bv, Pv = defect_BP(L3[0](tv, 0.0), R3[0](tv, 0.0), sgv)
    return El + Er + Bv, Pl + Pr + Pv, Bv, Pv


def defect_residual(L3, R3, sgv, tv):
    ph0, pt_, px_ = (fn(tv, 0.0) for fn in L3)
    ps0, qt_, qx_ = (fn(tv, 0.0) for fn in R3)
    Bp_ = mv * (sgv * np.sin((ph0 + ps0) / 2) + np.sin((ph0 - ps0) / 2) / sgv)
    Bq_ = mv * (sgv * np.sin((ph0 + ps0) / 2) - np.sin((ph0 - ps0) / 2) / sgv)
    return max(abs(px_ - (qt_ - Bp_)), abs(qx_ - (pt_ + Bq_)))


zero3 = (lambda tv, xv: 0.0 * xv, lambda tv, xv: 0.0 * xv, lambda tv, xv: 0.0 * xv)
emit = lamb(4 * sp.atan(kink_expr(thd, 0)))
times = (-30.0, -1.0, 0.5, 3.0, 30.0)
vals = [totals(zero3, emit, sg_v, tv) for tv in times]
res0 = max(defect_residual(zero3, emit, sg_v, tv) for tv in times)
Es_ = np.array([v_[0] for v_ in vals])
Ps_ = np.array([v_[1] for v_ in vals])
ok_j1 = res0 < 1e-12 and np.ptp(Es_) < 1e-8 and np.ptp(Ps_) < 1e-8 \
    and abs(Es_[0] - Mv * np.cosh(thd)) < 1e-8 and abs(Ps_[0] - Mv * np.sinh(thd)) < 1e-8 \
    and abs(vals[0][2] - Mv * np.cosh(thd)) < 1e-8 and abs(vals[-1][2]) < 1e-8
# Bianchi pair: Phi = Phi_1 (rapidity th1), Psi = Phi_12 = B_{a2} Phi_1 with a2 = -sigma
th1 = -0.3
a1v, a2v = -np.exp(-th1), -sg_v
coef = (a1v + a2v) / (a1v - a2v)
Y1 = kink_expr(th1, 0.4)
Y2 = kink_expr(thd, -0.7)
left1 = lamb(4 * sp.atan(Y1))
right1 = lamb(4 * sp.atan(coef * (Y1 - Y2) / (1 + Y1 * Y2)))
times1 = (-10.0, -2.0, 0.0, 1.5, 9.0)
res1 = max(defect_residual(left1, right1, sg_v, tv) for tv in times1)
vals1 = [totals(left1, right1, sg_v, tv) for tv in times1]
ok_j2 = res1 < 1e-12 and np.ptp([v_[0] for v_ in vals1]) < 1e-8 and np.ptp([v_[1] for v_ in vals1]) < 1e-8
sE1, sP1 = np.ptp([v_[0] for v_ in vals1]), np.ptp([v_[1] for v_ in vals1])
report(f"(j) E + B and P + P_d constant on the emission solution (spreads {np.ptp(Es_):.0e}, {np.ptp(Ps_):.0e}; "
       f"B goes from M cosh theta_d to 0) and on the Bianchi pair (Phi_1, Phi_12) (spreads {sE1:.0e}, {sP1:.0e}); "
       f"defect residuals {max(res0, res1):.0e}", ok_j1 and ok_j2)

# ---------------------------------------------------------------------------------------------------------
# (k) light-cone split and stored energy-momentum
Bd = -2 * m0 * (sig * (sp.cos(A_) - 1) + (sp.cos(B_) - 1) / sig)
Pd = 2 * m0 * (sig * (sp.cos(A_) - 1) - (sp.cos(B_) - 1) / sig)
ok_k1 = zero(Bd + Pd + 4 * m0 / sig * (sp.cos(B_) - 1)) and zero(Bd - Pd + 4 * m0 * sig * (sp.cos(A_) - 1))
# P_d is D (g - f) up to a constant, with f = -2 m0 sigma cos(s/2), g = -2 m0 cos(d/2)/sigma
ok_k1 = ok_k1 and zero(Pd - ((-2 * m0 / sig * sp.cos(B_)) - (-2 * m0 * sig * sp.cos(A_))) - 2 * m0 * (1 / sig - sig))
thd_s = sp.Symbol('theta_d', real=True)
jump = {P: 0, Q: 2 * sp.pi}
ok_k2 = zero((Bd.subs(jump) - 8 * m0 * sp.cosh(thd_s)).subs(sig, sp.exp(-thd_s)).rewrite(sp.exp)) \
    and zero((Pd.subs(jump) - 8 * m0 * sp.sinh(thd_s)).subs(sig, sp.exp(-thd_s)).rewrite(sp.exp))
ok_k3 = True
for n_ in range(-2, 3):
    pt_ = {P: 0, Q: 2 * sp.pi * n_}
    grad = [sp.diff(Bd, vv).subs(pt_) for vv in (P, Q)]
    H = sp.hessian(Bd, (P, Q)).subs(pt_).subs({sig: sp.Rational(7, 10), m0: 1})
    ev = [float(sp.re(e_)) for e_ in H.eigenvals()]
    ok_k3 = ok_k3 and all(zero(g_) for g_ in grad) and \
        (all(e_ > 0 for e_ in ev) if n_ % 2 == 0 else all(e_ < 0 for e_ in ev))
report("(k) B + P_d = -4 m0 (cos(d/2) - 1)/sigma, B - P_d = -4 m0 sigma (cos(s/2) - 1); a static jump of 2 pi "
       "stores M (cosh theta_d, sinh theta_d); even jumps are minima of B, odd jumps maxima",
       ok_k1 and ok_k2 and ok_k3)

# ---------------------------------------------------------------------------------------------------------
# (l) no reflection from the fluxes; the linearized frozen Baecklund conditions
Rc = complex(*rng.normal(size=2)) * 0.4
kv, mv2 = 0.8, 1.1
omv = np.sqrt(kv**2 + mv2**2)
avgE, avgP = avg_flux(1, Rc, kv, mv2, 0.37)
ok_l1 = abs(avgE + omv * kv * (1 - abs(Rc)**2) / 2) < 1e-10 and abs(avgP - kv**2 * (1 + abs(Rc)**2) / 2) < 1e-10
avgE_t, avgP_t = avg_flux(Rc, 0, kv, mv2, -0.5)          # transmitted wave of amplitude Rc
ok_l1 = ok_l1 and abs(avgE_t + omv * kv * abs(Rc)**2 / 2) < 1e-10 and abs(avgP_t - kv**2 * abs(Rc)**2 / 2) < 1e-10
ok_l2 = True
for th, thd_v in [(0.43, -0.2), (-0.7, 0.9), (1.4, 1.1)]:
    sgv = np.exp(-thd_v)
    E_, p_ = mv * np.cosh(th), mv * np.sinh(th)
    Rs, Ts = sp.symbols('R T')
    B2 = mv * sgv * (P + Q)**2 / 4 + mv / sgv * (P - Q)**2 / 4      # quadratic part of the frozen BT potential
    # Phi = e^{-iEt}(e^{ipx} + R e^{-ipx}), Psi = e^{-iEt} T e^{ipx}, at x = 0
    phi0, psi0 = 1 + Rs, Ts
    eq1 = sp.I * p_ * (1 - Rs) - (-sp.I * E_ * Ts - sp.diff(B2, P).subs({P: phi0, Q: psi0}))
    eq2 = sp.I * p_ * Ts - (-sp.I * E_ * (1 + Rs) + sp.diff(B2, Q).subs({P: phi0, Q: psi0}))
    sol = sp.solve([eq1, eq2], [Rs, Ts], dict=True)[0]
    Tnum = complex(sol[Ts])
    Tbook = (np.sinh(thd_v) + 1j * np.cosh(th)) / (np.cosh(thd_v) - 1j * np.sinh(th))
    ok_l2 = ok_l2 and abs(complex(sol[Rs])) < 1e-12 and abs(abs(Tnum) - 1) < 1e-12 and abs(Tnum - Tbook) < 1e-12
report("(l) averaged energy flux -omega k (1-|R|^2)/2 and momentum flux k^2 (1+|R|^2)/2 on the left; the "
       "linearized frozen Baecklund conditions give R = 0, |T| = 1 and the book's eq-shg-linear-transmission",
       ok_l1 and ok_l2)

# ---------------------------------------------------------------------------------------------------------
# (m) type I_+ via Psi -> -Psi
B_I = -2 * m0 * (sig * sp.cos(A_) + sp.cos(B_) / sig)
Qn = sp.Symbol('Psi_prime')
B_In = B_I.subs(Q, Qn)
# D = 1 conditions for (Phi, Psi'), then Psi' = -Psi (so Psi'_t = -Psi_t, Psi'_x = -Psi_x)
c1_ = (Px - (-Qt - sp.diff(B_In, P))).subs(Qn, -Q)
c2_ = (-Qx - (Pt + sp.diff(B_In, Qn))).subs(Qn, -Q)
B_plus = -2 * m0 * (sig * sp.cos(B_) + sp.cos(A_) / sig)             # sigma and 1/sigma exchanged
e1p = Px - (-Qt - sp.diff(B_plus, P))                                # D = -1: Phi_x = -Psi_t - B_Phi
e2p = Qx - (-Pt + sp.diff(B_plus, Q))                                # D = -1: Psi_x = -Phi_t + B_Psi
report("(m) Psi -> -Psi maps the D = 1 defect with sigma cos(s/2) + cos(d/2)/sigma to the D = -1 defect with "
       "sigma and 1/sigma exchanged (type I_+)",
       zero(c1_ - e1p) and zero(c2_ + e2p) and zero(B_In.subs(Qn, -Q) - B_plus))

# ---------------------------------------------------------------------------------------------------------
# (n) the massless limit and the scaling limit
FL, FR = sp.Function('FL'), sp.Function('FR')
ph_full = FL(tt + xx) + FR(tt - xx)
ps_glued = FL(tt + xx) - FR(tt - xx)
cA = lambda ph_, ps_: sp.diff(ph_, xx) - sp.diff(ps_, tt)           # Phi_x - Psi_t   (m0 = 0)
cB = lambda ph_, ps_: sp.diff(ps_, xx) - sp.diff(ph_, tt)           # Psi_x - Phi_t
ok_n1 = zero(cA(ph_full, ps_glued)) and zero(cB(ph_full, ps_glued))
GL, GR = sp.Function('GL'), sp.Function('GR')
ps_full = GL(tt + xx) + GR(tt - xx)
dplus = lambda f_: sp.diff(f_, tt) + sp.diff(f_, xx)
dminus = lambda f_: sp.diff(f_, tt) - sp.diff(f_, xx)
ok_n1 = ok_n1 and zero(cA(ph_full, ps_full) - cB(ph_full, ps_full) - dplus(ph_full - ps_full)) \
    and zero(cA(ph_full, ps_full) + cB(ph_full, ps_full) + dminus(ph_full + ps_full))
# d_+ kills right-movers, d_- kills left-movers
ok_n1 = ok_n1 and zero(dplus(FR(tt - xx))) and zero(dminus(FL(tt + xx)))
mu, omg = sp.symbols('mu omega', positive=True)
Bsc = -2 * mu * sp.cos(B_)
dm_sum = ((Pt - Px) + (Qt - Qx)).subs({Px: Qt - sp.diff(Bsc, P), Qx: Pt + sp.diff(Bsc, Q)})
ok_n2 = zero(dm_sum - 2 * mu * sp.sin(B_))
Tm = sp.Symbol('Tm')
# right-movers u = e^{-i omega t}, w = Tm e^{-i omega t} at x = 0: d_- acts as 2 d/dt; 2 mu sin(d/2) ~ mu (u - w)
Tm_sol = sp.solve(sp.Eq(2 * (-sp.I * omg) * (1 + Tm), mu * (1 - Tm)), Tm)[0]
ok_n2 = ok_n2 and zero(Tm_sol - (mu + 2 * sp.I * omg) / (mu - 2 * sp.I * omg))
ok_n3 = True
for m0v in (1e-3, 1e-5):
    om_v, mu_v = 0.7, 1.9
    thv, thdv = np.log(2 * om_v / m0v), np.log(mu_v / m0v)
    Tmass = (np.sinh(thdv) + 1j * np.cosh(thv)) / (np.cosh(thdv) - 1j * np.sinh(thv))
    ok_n3 = ok_n3 and abs(Tmass - (mu_v + 2j * om_v) / (mu_v - 2j * om_v)) < 10 * m0v
report("(n) at m0 = 0 the conditions are d_+(Phi-Psi) = 0, d_-(Phi+Psi) = 0, i.e. Psi_L = Phi_L, Psi_R = -Phi_R; "
       "in the scaling limit right-movers are transmitted with (mu + 2 i omega)/(mu - 2 i omega), the limit of "
       "the massive factor", ok_n1 and ok_n2 and ok_n3)

# ---------------------------------------------------------------------------------------------------------
# (o) spin 3: BCZ 2003 (2.6)-(2.9) with lambda1 = lambda2 = lambda, lambda^2 = -1/4 (V''' = 4 lambda^2 V')
l2 = sp.Rational(-1, 4)
ok_o0 = zero(sp.diff(V_sg(P), P, 3) - 4 * l2 * sp.diff(V_sg(P), P))
Bn = B_I
conds = [2 * (l2 * sp.diff(Bn, P)**2 - l2 * sp.diff(Bn, Q)**2) - sp.diff(V_sg(P), P, 2) + sp.diff(V_sg(Q), Q, 2),
         sp.diff(Bn, Q, 2, P) - l2 * sp.diff(Bn, P), sp.diff(Bn, P, 2, Q) - l2 * sp.diff(Bn, Q),
         sp.diff(Bn, Q, 3) - l2 * sp.diff(Bn, Q), sp.diff(Bn, P, 3) - l2 * sp.diff(Bn, P)]
ok_o1 = all(zero(sp.expand_trig(c_)) for c_ in conds)
gm = gf(P - Q)
ok_o2 = zero(sp.diff(gm, P, 2) - sp.diff(gm, Q, 2)) and zero(sp.diff(gm, P)**2 / 2 - sp.diff(gm, Q)**2 / 2)
lam_s = sp.Symbol('lambda', nonzero=True)
ok_o3 = not zero(sp.diff((P - Q)**4, Q, 3) - lam_s**2 * sp.diff((P - Q)**4, Q)) \
    and zero(sp.diff(sp.exp(lam_s * (P - Q)), Q, 3) - lam_s**2 * sp.diff(sp.exp(lam_s * (P - Q)), Q))
report("(o) the frozen Baecklund potential satisfies BCZ's spin-3 conditions (2.6)-(2.9) with lambda^2 = -1/4; "
       "massless free fields: any g(Phi - Psi) conserves momentum, but (2.9) admits e^{lambda(Phi-Psi)}, "
       "not (Phi - Psi)^4", ok_o0 and ok_o1 and ok_o2 and ok_o3)

# massless free with B = g(Phi - Psi): chiral right-mover law and the quartic charge
cc = sp.Symbol('c')
FLf, FRf, GRf = sp.Function('FL'), sp.Function('FR'), sp.Function('GR')
phm = FLf(tt + xx) + FRf(tt - xx)
psm = FLf(tt + xx) + cc + GRf(tt - xx)                 # left-movers glued, Psi_L = Phi_L + c
Bg = gf(P - Q)
condA = sp.diff(phm, xx) - (sp.diff(psm, tt) - sp.diff(Bg, P).subs({P: phm, Q: psm}))
condB = sp.diff(psm, xx) - (sp.diff(phm, tt) + sp.diff(Bg, Q).subs({P: phm, Q: psm}))
# at x = 0: condA - condB and condA + condB in terms of u = FR(t), w = GR(t)
u_, w_ = sp.Function('u')(tt), sp.Function('w')(tt)
at0 = lambda e_: sp.simplify(e_.subs(xx, 0).doit().subs({FRf(tt): u_, GRf(tt): w_}).doit())
diffAB = at0(condA - condB)
sumAB = at0(condA + condB)
# expected: difference vanishes identically; sum: -2(u' + w') + 2 g'(u - w - c) = 0
Xs = u_ - w_ - cc
gp = sp.Function('gp')
sum_expected = -2 * (sp.diff(u_, tt) + sp.diff(w_, tt)) + 2 * sp.Subs(sp.diff(gf(yv), yv), yv, Xs).doit()
ok_o4 = zero(diffAB) and zero(sp.simplify(sumAB - sum_expected))
# quartic charge: 16 (w'^4 - u'^4) with u' + w' = g', u' - w' = X'
gpv, Xd = sp.symbols('gp Xdot')
ud, wd = (gpv + Xd) / 2, (gpv - Xd) / 2
ok_o5 = zero(16 * (wd**4 - ud**4) + 8 * gpv * Xd * (gpv**2 + Xd**2))
report("(p) massless free fields with B = g(Phi - Psi): left-movers glued, right-movers obey u' + w' = g'(u - w - c); "
       "dQ_4/dt = 16 (w'^4 - u'^4) = -8 g' X' (g'^2 + X'^2)", ok_o4 and ok_o5)

# ---------------------------------------------------------------------------------------------------------
# (q) exercises: free limit, sigma -> -sigma, positivity, freezing Psi to a constant
Sg, Dl = sp.symbols('Sigma Delta')
Bsd = -2 * m0 * (sig * (sp.cos(Sg / 2) - 1) + (sp.cos(Dl / 2) - 1) / sig)
quad_part = sp.series(Bsd, Sg, 0, 3).removeO()
quad_part = sp.series(quad_part, Dl, 0, 3).removeO()
ok_p1 = zero(quad_part - (m0 * sig / 4 * Sg**2 + m0 / (4 * sig) * Dl**2))
ok_p2 = zero(Bd.subs(Q, Q + 2 * sp.pi) - Bd.subs(sig, -sig) - (Bd.subs(Q, Q + 2 * sp.pi) - Bd.subs(sig, -sig)).subs({P: 0, Q: 0}))
Bdn = sp.lambdify((P, Q, sig, m0), Bd, 'numpy')
grid = rng.uniform(-10, 10, size=(400, 2))
ok_p3 = all(Bdn(p0, q0, sv, 1.0) >= -1e-12 for (p0, q0) in grid for sv in (0.3, 1.0, 2.5))
ps0 = sp.Symbol('psi0', real=True)
sg_r = sp.Symbol('sigma_r', positive=True)
cplx = sg_r * sp.exp(sp.I * ps0 / 2) + sp.exp(-sp.I * ps0 / 2) / sg_r
lhs_p = sg_r * sp.cos((P + ps0) / 2) + sp.cos((P - ps0) / 2) / sg_r
# the identity sigma cos((Phi+psi0)/2) + cos((Phi-psi0)/2)/sigma = Re(c e^{i Phi/2}), and surjectivity onto c != 0
ok_p4 = True
lhs_f = sp.lambdify((P, ps0, sg_r), lhs_p, 'numpy')
c_f = sp.lambdify((ps0, sg_r), cplx, 'numpy')
for _ in range(10):
    Pv_, psv, sgv_ = rng.normal(), rng.uniform(-6, 6), rng.uniform(0.2, 3)
    ok_p4 = ok_p4 and abs(lhs_f(Pv_, psv, sgv_) - (c_f(psv, sgv_) * np.exp(1j * Pv_ / 2)).real) < 1e-12
from scipy.optimize import brentq
ok_p5 = True
for _ in range(10):
    r_, chi_ = rng.uniform(0.05, 6), rng.uniform(-np.pi, np.pi)
    Xc, Yc = r_ * np.cos(chi_), r_ * np.sin(chi_)
    # ellipse X = a cos(psi0/2), Y = b sin(psi0/2), a = sigma + 1/sigma, b = sigma - 1/sigma; find a > 2
    h_ = lambda a_: Xc**2 / a_**2 + Yc**2 / (a_**2 - 4) - 1
    a_sol = brentq(h_, 2 + 1e-12, 1e6, xtol=1e-15, rtol=1e-15)
    sv = (a_sol + np.sqrt(a_sol**2 - 4)) / 2
    psv = 2 * np.arctan2(Yc / (sv - 1 / sv), Xc / (sv + 1 / sv))
    err_ = abs(c_f(psv, sv) - (Xc + 1j * Yc)) / r_
    ok_p5 = ok_p5 and err_ < 1e-9
report("(q) quadratic part of B is the free defect with lambda = sigma; sigma -> -sigma is Psi -> Psi + 2 pi; "
       "B >= 0 for sigma > 0; B(Phi, psi_0) = Re(c e^{i Phi/2}) up to the factor and constants, and every c != 0 "
       "arises from some sigma > 0 and psi_0", ok_p1 and ok_p2 and ok_p3 and ok_p4 and ok_p5)

raise SystemExit(0 if OK else 1)
