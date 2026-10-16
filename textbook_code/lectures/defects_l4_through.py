"""Checks for Lecture 4 of the defect lecture notes, "What gets through"
(lecture-notes/defects/04-what-gets-through.qmd).

Conventions of the lectures: Phi on x<0, Psi on x>0, dimensionless (Phi = beta_SG phi), V = m0^2 (1 - cos Phi);
E = (1/beta^2) int (Phi_t^2/2 + Phi_x^2/2 + V), P = -(1/beta^2) int Phi_t Phi_x, M = 8 m0/beta^2; d_pm = d_t pm d_x.
The type-I_- defect with sigma = e^{-theta_d} has the conditions, at x = 0,
    Phi_x = Psi_t - B_Phi,  Psi_x = Phi_t + B_Psi,   B = -2 m0 [sigma (cos((Phi+Psi)/2) - 1) + (cos((Phi-Psi)/2) - 1)/sigma],
equivalently d_+(Phi - Psi) = -2 sigma m0 sin((Phi+Psi)/2), d_-(Phi + Psi) = (2 m0/sigma) sin((Phi-Psi)/2).
A soliton of rapidity theta is Phi = 4 arctan Y, Y = exp(m0 (x cosh theta - t sinh theta) + c); Y -> -Y is the
antisoliton. Written from these definitions, independently of textbook_code/defect_classical_sg.py.
One PASS or FAIL line per check; exit status 1 on any FAIL. Units m0 = beta = 1 (M = 8) in the numerics.

  (a) The substitution Psi = 4 arctan(z Y): sin((Phi +- Psi)/2) = 2 (1 +- z) Y (1 -+ z Y^2)/((1+Y^2)(1+z^2 Y^2)),
      d_+ ln Y = m0 e^{-theta}, d_- ln Y = -m0 e^{theta}; the two conditions reduce to e^{-theta}(1-z) = -sigma (1+z)
      and sigma e^{theta} (1+z) = -(1-z), with the common solution z = (e^{-theta}+sigma)/(e^{-theta}-sigma)
      = -coth((theta-theta_d)/2), the Bianchi coefficient (a + a_d)/(a - a_d), a = -e^{-theta}, a_d = -sigma.
  (b) The pair solves the defect conditions for all t (residuals, complex c, both sides of theta_d), and is a
      Baecklund pair on the whole line; it is the limit of Lecture 2's (Phi_1, Phi_12) as soliton 2 is pushed away.
  (c) Kinematics: the static states (2 pi m, 2 pi n) store (0, 0) or M (cosh theta_d, sinh theta_d); after
      absorption by an odd state the outgoing E^2 - P^2 is 2 M^2 (1 - cosh(theta - theta_d)) <= 0; the reverse
      triangle inequality for future timelike vectors.
  (d) Reading z: |z| > 1 for real theta != theta_d; the outgoing centre measured on the exact solution sits at
      x = t tanh theta - (c + ln|z|)/(m0 cosh theta); time delay ln|z|/(m0 sinh theta); for theta > theta_d the
      outgoing field is an antikink, the jump grows by 4 pi and B returns to 0; E + B/beta^2 and P + P_d/beta^2
      are constant along a converted passage.
  (e) The limits: z -> 1 as theta_d -> +infinity and z -> -1 as theta_d -> -infinity; for large sigma the
      conditions force Psi = -Phi with Psi_x = -Phi_x (charge conjugation line).
  (f) Absorption: (soliton of rapidity theta_d, 2 pi) solves the conditions; defect values go from (2 pi, 2 pi)
      to (0, 2 pi); E + B and P + P_d constant along it; for theta_d < 0, (soliton, 2 pi) emits to the left and
      (0, soliton) absorbs from the right.
  (g) Capture: on the exact solution the defect values pass half-way (Phi(0) = pi, Psi(0) = +-pi) at times
      that differ by exactly ln|z|/(m0 sinh theta); in between, B/beta^2 rises to within O(theta - theta_d) of
      M cosh theta_d.
  (h) Several solitons: the two-soliton field with Y_k -> z_k Y_k on the right solves the conditions (soliton-
      antisoliton and soliton-soliton, one converted, complex c); final displacements measured on the exact
      solution are the sums of the defect and collision displacements in both orders of events.
  (i) A breather passes with Y_1 -> z(theta + i alpha) Y_1, Y_2 -> z(theta - i alpha) Y_2;
      |z(u + i alpha)|^2 = (cosh u + cos alpha)/(cosh u - cos alpha) > 1; its envelope is delayed.
  (j) Linear waves: allowing a reflected wave, the linearized conditions give R = 0 and
      T = (e^{-theta} + i sigma)/(e^{-theta} - i sigma) = z(theta - theta_d + i pi/2), the book's T^cl;
      |T| = 1, T(theta_d) = i, limits +-1; a wave from the right has bar R = 0, bar T(theta) = T(i pi - theta),
      T(theta) bar T(-theta) = 1 = T(theta) T(theta + i pi), T(theta) T(-theta) != 1; the soliton at theta + i pi/2
      is the plane wave; pole at theta_d - i pi/2, zero at theta_d + i pi/2.
  (k) One function: |z(u)| = 1 exactly on Im u = +-pi/2, > 1 for 0 < |Im u| < pi/2, < 1 for pi/2 < Im u < pi;
      z(u + i pi) = 1/z(u) = -tanh(u/2); and the odd-sector pair (Phi, 2 pi + 4 arctan(z_odd Y)) solves the
      conditions with z_odd = 1/z.
  (l) Half a soliton: z_odd^2 = A(theta - theta_d) = tanh^2, z^2 = 1/A; T_odd^2 = (sinh u + i)/(sinh u - i), the
      transmission of a linear wave through a kink (from the exact fluctuation (tanh(m0 x) - i k/m0) e^{ikx});
      T^2 = 1/T_kink; with the constants -M/2 (cosh, sinh) theta_d the odd and even states lie on the mass shell
      of mass M/2.
  (m) Two defects at x = 0 and x = L: Psi_1 = 4 arctan(z_1 Y), Psi_2 = 4 arctan(z_1 z_2 Y) satisfy both sets of
      conditions; two identical defects never convert.
  (n) Exercises: type I_+ has z_+ = -z; an antisoliton has the same z; a soliton from the right (theta < 0) is
      delayed, converted iff theta > theta_d; Wigner: M sinh theta Delta t = (8/beta^2) ln|z|; theta_d < 0
      converts every soliton from the left (NOT gate).
  (o) Odd sector near theta_d: the right defect value passes half-way before the left one, by
      ln(1/|z_odd|)/(m0 sinh theta): the stored soliton leaves before the incoming one is swallowed.

Runtime: a few seconds.
"""
import cmath
import math

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar

OK = True


def report(name, ok):
    global OK
    OK = OK and bool(ok)
    print(("PASS " if ok else "FAIL ") + name)


rng = np.random.default_rng(4)
tt, xx = sp.symbols('t x', real=True)
th, thd_s = sp.symbols('theta theta_d', real=True)
sg_s = sp.Symbol('sigma', positive=True)
z_s, Y_s = sp.symbols('z Y')
m0s = sp.Symbol('m0', positive=True)

# ---------------------------------------------------------------------------------------------------------
# (a) the substitution, symbolically and at random complex points
ok_a = True
for _ in range(6):
    Yv = complex(rng.normal(), rng.normal())
    zv = complex(rng.normal(), rng.normal())
    Ph, Ps = 4 * cmath.atan(Yv), 4 * cmath.atan(zv * Yv)
    den = (1 + Yv**2) * (1 + zv**2 * Yv**2)
    ok_a &= abs(cmath.sin((Ph + Ps) / 2) - 2 * (1 + zv) * Yv * (1 - zv * Yv**2) / den) < 1e-12
    ok_a &= abs(cmath.sin((Ph - Ps) / 2) - 2 * (1 - zv) * Yv * (1 + zv * Yv**2) / den) < 1e-12
Ylog = m0s * (xx * sp.cosh(th) - tt * sp.sinh(th))
ok_a &= sp.simplify(((sp.diff(Ylog, tt) + sp.diff(Ylog, xx)) - m0s * sp.exp(-th)).rewrite(sp.exp)) == 0
ok_a &= sp.simplify(((sp.diff(Ylog, tt) - sp.diff(Ylog, xx)) + m0s * sp.exp(th)).rewrite(sp.exp)) == 0
# d_+(Phi - Psi) with Phi = 4 atan Y, Psi = 4 atan(zY), d_+ Y = m0 e^{-theta} Y:
lhs1 = 4 * m0s * sp.exp(-th) * Y_s * (1 / (1 + Y_s**2) - z_s / (1 + z_s**2 * Y_s**2))
rhs1 = -2 * sg_s * m0s * 2 * (1 + z_s) * Y_s * (1 - z_s * Y_s**2) / ((1 + Y_s**2) * (1 + z_s**2 * Y_s**2))
lhs2 = -4 * m0s * sp.exp(th) * Y_s * (1 / (1 + Y_s**2) + z_s / (1 + z_s**2 * Y_s**2))
rhs2 = (2 * m0s / sg_s) * 2 * (1 - z_s) * Y_s * (1 + z_s * Y_s**2) / ((1 + Y_s**2) * (1 + z_s**2 * Y_s**2))
common = 4 * m0s * Y_s / ((1 + Y_s**2) * (1 + z_s**2 * Y_s**2))
red1 = sp.simplify((lhs1 - rhs1) / (common * (1 - z_s * Y_s**2)))
red2 = sp.simplify((lhs2 - rhs2) / (common * (1 + z_s * Y_s**2)))
ok_a &= sp.simplify(red1 - (sp.exp(-th) * (1 - z_s) + sg_s * (1 + z_s))) == 0
ok_a &= sp.simplify(red2 - (-sp.exp(th) * (1 + z_s) - (1 - z_s) / sg_s)) == 0
z1 = sp.solve(red1, z_s)
z2 = sp.solve(red2, z_s)
zform = (sp.exp(-th) + sg_s) / (sp.exp(-th) - sg_s)
ok_a &= len(z1) == 1 and len(z2) == 1 and sp.simplify(z1[0] - zform) == 0 and sp.simplify(z2[0] - zform) == 0
zcoth = (-sp.coth((th - thd_s) / 2)).rewrite(sp.exp)
ok_a &= sp.simplify(zform.subs(sg_s, sp.exp(-thd_s)) - zcoth) == 0
a_, ad_ = -sp.exp(-th), -sg_s
ok_a &= sp.simplify((a_ + ad_) / (a_ - ad_) - zform) == 0
report("(a) substituting Psi = 4 arctan(zY): both conditions reduce to linear equations for z with the same "
       "solution z = (e^-theta + sigma)/(e^-theta - sigma) = -coth((theta - theta_d)/2), the Bianchi coefficient",
       ok_a)

# ---------------------------------------------------------------------------------------------------------
# numerical machinery
m0 = 1.0
M = 8.0


def zfun(theta, thd):
    return (cmath.exp(-theta) + cmath.exp(-thd)) / (cmath.exp(-theta) - cmath.exp(-thd))


def Yexpr(theta, c):
    return sp.exp(m0 * (xx * sp.cosh(theta) - tt * sp.sinh(theta)) + c)


def lamb(F):
    """Phi = 4 atan F: return Phi, Phi_t, Phi_x as numeric functions of (t, x)."""
    Phi = 4 * sp.atan(F)
    fs = [sp.lambdify((tt, xx), e, 'numpy') for e in (Phi, sp.diff(Phi, tt), sp.diff(Phi, xx))]
    return fs


def residual(L3, R3, sgv, tv, x0=0.0):
    ph0, pt_, px_ = (complex(fn(tv, x0)) for fn in L3)
    ps0, qt_, qx_ = (complex(fn(tv, x0)) for fn in R3)
    Bp = m0 * (sgv * cmath.sin((ph0 + ps0) / 2) + cmath.sin((ph0 - ps0) / 2) / sgv)
    Bq = m0 * (sgv * cmath.sin((ph0 + ps0) / 2) - cmath.sin((ph0 - ps0) / 2) / sgv)
    return max(abs(px_ - (qt_ - Bp)), abs(qx_ - (pt_ + Bq)))


def defect_BP(ph0, ps0, sgv):
    Bv = -2 * m0 * (sgv * (math.cos((ph0 + ps0) / 2) - 1) + (math.cos((ph0 - ps0) / 2) - 1) / sgv)
    Pv = 2 * m0 * (sgv * (math.cos((ph0 + ps0) / 2) - 1) - (math.cos((ph0 - ps0) / 2) - 1) / sgv)
    return Bv, Pv


def bulk_EP(F3, tv, x0, x1):
    F, Ft, Fx = F3
    e = lambda xv: Ft(tv, xv)**2 / 2 + Fx(tv, xv)**2 / 2 + m0**2 * (1 - np.cos(F(tv, xv)))
    p = lambda xv: -Ft(tv, xv) * Fx(tv, xv)
    return (quad(e, x0, x1, limit=500, epsabs=1e-12, epsrel=1e-12)[0],
            quad(p, x0, x1, limit=500, epsabs=1e-12, epsrel=1e-12)[0])


def totals(L3, R3, sgv, tv, Lbox=50.0):
    El, Pl = bulk_EP(L3, tv, -Lbox, 0.0)
    Er, Pr = bulk_EP(R3, tv, 0.0, Lbox)
    Bv, Pv = defect_BP(float(L3[0](tv, 0.0)), float(R3[0](tv, 0.0)), sgv)
    return El + Er + Bv, Pl + Pr + Pv


# ---------------------------------------------------------------------------------------------------------
# (b) the pair is an exact solution; Baecklund pair on the line; limit of the Bianchi pair
thd = 0.9
sgv = math.exp(-thd)
ok_b = True
worst = 0.0
for theta in (0.3, 1.7, -0.6):
    zv = zfun(theta, thd).real
    c = complex(rng.normal(), rng.normal())
    L3 = lamb(Yexpr(theta, c))
    R3 = lamb(zv * Yexpr(theta, c))
    for tv in (-4.0, -0.7, 0.0, 1.3, 5.0):
        worst = max(worst, residual(L3, R3, sgv, tv))
        # Baecklund pair at other points of the line
        worst = max(worst, residual(L3, R3, sgv, tv, x0=float(rng.uniform(-3, 3))))
ok_b &= worst < 1e-12
# Bianchi: tan(Phi_12/4) = B (Y1 - Y2)/(1 + Y1 Y2), B = (a1+a2)/(a1-a2), a2 = -sigma; push soliton 2 to +infinity
theta = 0.4
Bc = zfun(theta, thd).real
pushed = 0.0
for c2 in (-30.0, -50.0):
    for (tv, xv) in ((0.0, 0.0), (-2.0, 1.0), (3.0, -1.5)):
        Y1 = math.exp(m0 * (xv * math.cosh(theta) - tv * math.sinh(theta)) + 0.3)
        Y2 = math.exp(m0 * (xv * math.cosh(thd) - tv * math.sinh(thd)) + c2)
        pushed = max(pushed, abs(Bc * (Y1 - Y2) / (1 + Y1 * Y2) - Bc * Y1) / abs(Bc * Y1))
ok_b &= pushed < 1e-10
report(f"(b) (4 arctan Y, 4 arctan(zY)) solves the defect conditions for all t and is a Baecklund pair on the "
       f"line (residual {worst:.0e}); it is Lecture 2's (Phi_1, Phi_12) with soliton 2 pushed away "
       f"(relative difference {pushed:.0e})", ok_b)

# ---------------------------------------------------------------------------------------------------------
# (c) kinematics
ok_c = True
for mm in range(-2, 3):
    for nn in range(-2, 3):
        Bv, Pv = defect_BP(2 * math.pi * mm, 2 * math.pi * nn, sgv)
        if (nn - mm) % 2 == 0:
            ok_c &= abs(Bv) < 1e-12 and abs(Pv) < 1e-12
        else:
            ok_c &= abs(Bv - M * math.cosh(thd)) < 1e-12 and abs(Pv - M * math.sinh(thd)) < 1e-12
Ms = sp.Symbol('M', positive=True)
pout2 = (Ms * (sp.cosh(th) - sp.cosh(thd_s)))**2 - (Ms * (sp.sinh(th) - sp.sinh(thd_s)))**2
ok_c &= sp.simplify((pout2 - 2 * Ms**2 * (1 - sp.cosh(th - thd_s))).rewrite(sp.exp)) == 0
# reverse triangle inequality: invariant mass of a sum of future timelike vectors >= sum of masses
for _ in range(200):
    ms_ = rng.uniform(0.1, 2, size=3)
    ths = rng.normal(size=3) * 2
    E_ = sum(m * math.cosh(t) for m, t in zip(ms_, ths))
    P_ = sum(m * math.sinh(t) for m, t in zip(ms_, ths))
    ok_c &= E_**2 - P_**2 >= ms_.sum()**2 - 1e-9
report("(c) static states store (0, 0) (even jumps) or M (cosh, sinh) theta_d (odd jumps); after an odd final state "
       "the outgoing E^2 - P^2 = 2 M^2 (1 - cosh(theta - theta_d)) <= 0; invariant mass of a sum >= sum of masses",
       ok_c)

# ---------------------------------------------------------------------------------------------------------
# (d) reading z: delay, conversion, bookkeeping, conservation along a converted passage
ok_d = True
for theta in np.linspace(-3, 3, 61):
    if abs(theta - thd) > 1e-6:
        ok_d &= abs(zfun(theta, thd)) > 1
ok_d &= zfun(0.3, thd).real > 0 and zfun(1.7, thd).real < 0
for theta, c in ((0.3, 0.2), (1.7, -0.4)):
    zv = zfun(theta, thd).real
    Psi = lambda tv, xv: 4 * math.atan(zv * math.exp(m0 * (xv * math.cosh(theta) - tv * math.sinh(theta)) + c))
    target = math.pi if zv > 0 else -math.pi
    for tv in (8.0, 14.0):
        xpred = tv * math.tanh(theta) - (c + math.log(abs(zv))) / (m0 * math.cosh(theta))
        xc = brentq(lambda xv: Psi(tv, xv) - target, xpred - 5, xpred + 5, xtol=1e-14)
        ok_d &= abs(xc - xpred) < 1e-10
    # orientation of the outgoing field: kink for z > 0, antikink for z < 0
    slope = Psi(10.0, 10.0 * math.tanh(theta) + 0.01) - Psi(10.0, 10.0 * math.tanh(theta) - 0.01)
    ok_d &= (slope > 0) == (zv > 0)
# time delay = displacement / velocity
theta = 0.3
zv = zfun(theta, thd).real
ok_d &= abs((math.log(zv) / (m0 * math.cosh(theta))) / math.tanh(theta) - math.log(zv) / (m0 * math.sinh(theta))) \
    < 1e-14
# bookkeeping for the converted passage, Psi shifted by 4 pi: defect values (2 pi, 2 pi) -> (0, 4 pi)
theta, c = 1.7, -0.4
zv = zfun(theta, thd).real
L3 = lamb(Yexpr(theta, c))
R3n = lamb(zv * Yexpr(theta, c))
R3 = [lambda tv, xv, f=R3n[0]: f(tv, xv) + 4 * math.pi, R3n[1], R3n[2]]
early = (float(L3[0](-60.0, 0.0)), float(R3[0](-60.0, 0.0)))
late = (float(L3[0](60.0, 0.0)), float(R3[0](60.0, 0.0)))
ok_d &= abs(early[0] - 2 * math.pi) < 1e-9 and abs(early[1] - 2 * math.pi) < 1e-9
ok_d &= abs(late[0]) < 1e-9 and abs(late[1] - 4 * math.pi) < 1e-9
ok_d &= abs(defect_BP(*early, sgv)[0]) < 1e-9 and abs(defect_BP(*late, sgv)[0]) < 1e-9
tots = [totals(L3, R3, sgv, tv) for tv in (-12.0, -1.0, 0.0, 0.8, 12.0)]
sE = np.ptp([v[0] for v in tots])
sP = np.ptp([v[1] for v in tots])
ok_d &= sE < 1e-8 and sP < 1e-8 and abs(tots[0][0] - M * math.cosh(theta)) < 1e-8
report(f"(d) |z| > 1; outgoing centre at x = t tanh theta - (c + ln|z|)/(m0 cosh theta); antikink iff theta > "
       f"theta_d; converted passage takes the defect from (2 pi, 2 pi) to (0, 4 pi) with B = 0 at both ends; "
       f"E + B and P + P_d constant along it (spreads {sE:.0e}, {sP:.0e})", ok_d)

# ---------------------------------------------------------------------------------------------------------
# (e) limits of the defect parameter
ok_e = abs(zfun(0.5, 30.0) - 1) < 1e-12 and abs(zfun(0.5, -30.0) + 1) < 1e-12
# large sigma: Psi = -Phi at x = 0 up to O(1/sigma); with Psi = -Phi the sum of the conditions gives
# Phi_x + Psi_x = d_t(Phi + Psi) - (2 m0/sigma) sin((Phi - Psi)/2) -> 0
P_, Q_ = sp.symbols('P Q')
Bsym = -2 * m0s * (sg_s * (sp.cos((P_ + Q_) / 2) - 1) + (sp.cos((P_ - Q_) / 2) - 1) / sg_s)
sum_cond = sp.simplify((sp.diff(Bsym, Q_) - sp.diff(Bsym, P_)))
ok_e &= sp.simplify(sum_cond + 2 * m0s / sg_s * sp.sin((P_ - Q_) / 2)) == 0
# the difference of the conditions: d_+(Phi - Psi) = -2 sigma m0 sin((Phi+Psi)/2) forces sin((Phi+Psi)/2) = O(1/sigma)
diff_cond = sp.simplify(-sp.diff(Bsym, P_) - sp.diff(Bsym, Q_))
ok_e &= sp.simplify(diff_cond + 2 * sg_s * m0s * sp.sin((P_ + Q_) / 2)) == 0
report("(e) z -> 1 as theta_d -> +infinity and z -> -1 as theta_d -> -infinity; for large sigma the conditions "
       "force sin((Phi+Psi)/2) -> 0 and Phi_x + Psi_x -> d_t(Phi + Psi): the charge-conjugation line", ok_e)

# ---------------------------------------------------------------------------------------------------------
# (f) absorption and emission
Yd = Yexpr(thd_s, 0)
Phi_abs = 4 * sp.atan(Yd)
Psi_abs = 2 * sp.pi
r1 = sp.diff(Phi_abs, xx) - (sp.diff(Psi_abs, tt) - sp.diff(Bsym, P_).subs({P_: Phi_abs, Q_: Psi_abs}))
r2 = sp.diff(Psi_abs, xx) - (sp.diff(Phi_abs, tt) + sp.diff(Bsym, Q_).subs({P_: Phi_abs, Q_: Psi_abs}))
subsd = {sg_s: sp.exp(-thd_s), m0s: 1}
ok_f = True
for _ in range(5):
    vals = {thd_s: rng.normal(), tt: rng.normal(), xx: rng.normal()}
    ok_f &= abs(complex(sp.N(r1.subs(subsd).subs(vals)))) < 1e-12 and abs(complex(sp.N(r2.subs(subsd).subs(vals)))) < 1e-12
for thd_v in (0.6, -0.5):
    sgd = math.exp(-thd_v)
    L3 = lamb(Yexpr(thd_v, 0.3))
    R3 = [lambda tv, xv: 2 * math.pi + 0 * xv, lambda tv, xv: 0 * xv, lambda tv, xv: 0 * xv]
    tots = [totals(L3, R3, sgd, tv) for tv in (-15.0, -1.0, 0.0, 2.0, 15.0)]
    ok_f &= np.ptp([v[0] for v in tots]) < 1e-8 and np.ptp([v[1] for v in tots]) < 1e-8
    early = float(L3[0](-80.0, 0.0))
    late = float(L3[0](80.0, 0.0))
    if thd_v > 0:      # absorbed from the left: (2 pi, 2 pi) -> (0, 2 pi)
        ok_f &= abs(early - 2 * math.pi) < 1e-9 and abs(late) < 1e-9
        ok_f &= abs(defect_BP(late, 2 * math.pi, sgd)[0] - M * math.cosh(thd_v)) < 1e-9
    else:              # theta_d < 0: (0, 2 pi) -> (2 pi, 2 pi), a soliton leaves to the left
        ok_f &= abs(early) < 1e-9 and abs(late - 2 * math.pi) < 1e-9
        ok_f &= abs(tots[0][0] - M * math.cosh(thd_v)) < 1e-8
# for theta_d < 0, (0, soliton theta_d) absorbs from the right: Psi(0) goes 0 -> 2 pi
thd_v = -0.5
Psi0 = lambda tv: 4 * math.atan(math.exp(-m0 * tv * math.sinh(thd_v)))
ok_f &= abs(Psi0(-80.0)) < 1e-9 and abs(Psi0(80.0) - 2 * math.pi) < 1e-9
report("(f) (soliton of rapidity theta_d, 2 pi) solves the conditions; it absorbs from the left for theta_d > 0 "
       "(defect values (2 pi, 2 pi) -> (0, 2 pi), stored M cosh theta_d) and emits to the left for theta_d < 0; "
       "E + B and P + P_d constant; (0, soliton) absorbs from the right for theta_d < 0", ok_f)

# ---------------------------------------------------------------------------------------------------------
# (g) capture near theta_d
thd_g = 0.8
sgg = math.exp(-thd_g)
ok_g = True
for eps in (5e-4, -5e-4, 2e-2):
    theta = thd_g + eps
    zv = zfun(theta, thd_g).real
    c = 0.0
    Phi0 = lambda tv: 4 * math.atan(math.exp(-m0 * tv * math.sinh(theta) + c))
    Psi0 = lambda tv: 4 * math.atan(zv * math.exp(-m0 * tv * math.sinh(theta) + c))
    t1 = brentq(lambda tv: Phi0(tv) - math.pi, -50, 50, xtol=1e-14)
    tgt = math.pi if zv > 0 else -math.pi
    t2 = brentq(lambda tv: Psi0(tv) - tgt, -50, 80, xtol=1e-14)
    ok_g &= abs((t2 - t1) - math.log(abs(zv)) / (m0 * math.sinh(theta))) < 1e-10
    tm = 0.5 * (t1 + t2)
    Bmid = defect_BP(Phi0(tm), Psi0(tm), sgg)[0]
    if abs(eps) < 1e-3:
        ok_g &= abs(Bmid / (M * math.cosh(thd_g)) - 1) < 2 * abs(eps)
capt = math.log(abs(zfun(thd_g + 5e-4, thd_g))) / (m0 * math.sinh(thd_g + 5e-4))
report(f"(g) the defect values pass half-way at times differing by exactly ln|z|/(m0 sinh theta) (capture time "
       f"{capt:.1f}/m0 for theta - theta_d = 5e-4); in between B/beta^2 is close to M cosh theta_d", ok_g)

# ---------------------------------------------------------------------------------------------------------
# (h) several solitons
def two_sol(Y1, Y2, th1, th2):
    Bc = -1 / sp.tanh(sp.Rational(1, 2) * (th1 - th2))
    return Bc * (Y1 - Y2) / (1 + Y1 * Y2)


ok_h = True
worst = 0.0
for (th1, th2, sign2) in ((1.5, 0.4, 1), (1.5, 0.4, -1), (0.2, -0.7, 1)):
    c1, c2 = complex(0.3, 0.2), complex(-0.5, -0.1)
    z1, z2 = zfun(th1, thd).real, zfun(th2, thd).real
    Y1, Y2 = Yexpr(th1, c1), sign2 * Yexpr(th2, c2)
    L3 = lamb(two_sol(Y1, Y2, th1, th2))
    R3 = lamb(two_sol(z1 * Y1, z2 * Y2, th1, th2))
    for tv in (-3.0, -0.5, 0.0, 0.7, 4.0):
        worst = max(worst, residual(L3, R3, sgv, tv))
ok_h &= worst < 1e-11


# classical factorization: two histories, measured on the exact solution
def centre(F, tv, xguess):
    """position of the maximum of |dPhi/dx| for Phi = 4 atan F near xguess (well-separated solitons)."""
    Fx = sp.lambdify(xx, sp.diff(4 * sp.atan(F.subs(tt, tv)), xx), 'numpy')
    res = minimize_scalar(lambda xv: -abs(Fx(xv)), bracket=(xguess - 1.0, xguess, xguess + 1.0), tol=1e-12)
    return res.x


th1, th2 = 1.5, 0.4                 # soliton 1 fast (converted), soliton 2 slow (passes)
z1, z2 = zfun(th1, thd).real, zfun(th2, thd).real
A12 = math.tanh((th1 - th2) / 2) ** 2
shifts, sides = [], []
for c2 in (-6.0, 6.0):              # the collision happens on the right (c2 < 0) or on the left (c2 > 0)
    Y1, Y2 = Yexpr(th1, 0.0), Yexpr(th2, c2)
    F_left = two_sol(Y1, Y2, th1, th2)
    F_right = two_sol(z1 * Y1, z2 * Y2, th1, th2)
    t_in, t_out = -40.0, 40.0
    x_in = centre(F_left, t_in, t_in * math.tanh(th1) + 0.5 * math.log(A12) / (m0 * math.cosh(th1)))
    x_out = centre(F_right, t_out, t_out * math.tanh(th1)
                   - (math.log(abs(z1)) + 0.5 * math.log(A12)) / (m0 * math.cosh(th1)))
    shifts.append((x_out - t_out * math.tanh(th1)) - (x_in - t_in * math.tanh(th1)))
    t_coll = -c2 / (m0 * math.cosh(th2) * (math.tanh(th1) - math.tanh(th2)))
    sides.append(t_coll * math.tanh(th1))
pred = -(math.log(abs(z1)) + math.log(A12)) / (m0 * math.cosh(th1))
ok_h &= sides[0] > 0 > sides[1]
ok_h &= abs(shifts[0] - shifts[1]) < 1e-6 and abs(shifts[0] - pred) < 1e-6
report(f"(h) two solitons pass with Y_k -> z_k Y_k (residual {worst:.0e}); the final displacement of the fast "
       f"soliton is -(ln|z_1| + ln A)/(m0 cosh theta_1) whether the collision happens before or after the defect "
       f"(measured {shifts[0]:.6f}, {shifts[1]:.6f})", ok_h)

# ---------------------------------------------------------------------------------------------------------
# (i) a breather
theta, alpha = 0.5, 0.7
zp, zm = zfun(theta + 1j * alpha, thd), zfun(theta - 1j * alpha, thd)
Y1 = sp.exp(m0 * (xx * sp.cosh(theta + sp.I * alpha) - tt * sp.sinh(theta + sp.I * alpha)))
Y2 = sp.exp(m0 * (xx * sp.cosh(theta - sp.I * alpha) - tt * sp.sinh(theta - sp.I * alpha)))
Bb = -1 / sp.tanh(sp.I * alpha)
L3 = lamb(Bb * (Y1 - Y2) / (1 + Y1 * Y2))
R3 = lamb(Bb * (zp * Y1 - zm * Y2) / (1 + zp * zm * Y1 * Y2))
worst = max(residual(L3, R3, sgv, tv) for tv in (-3.0, 0.0, 0.4, 2.5))
ok_i = worst < 1e-11 and abs(zm - zp.conjugate()) < 1e-14
ok_i &= abs(L3[0](0.3, -0.4).imag) < 1e-12 and abs(R3[0](0.3, 0.4).imag) < 1e-12
for u in np.linspace(-3, 3, 13):
    for al in (0.2, 0.9, 1.4):
        zz = -1 / cmath.tanh((u + 1j * al) / 2)
        ok_i &= abs(abs(zz)**2 - (math.cosh(u) + math.cos(al)) / (math.cosh(u) - math.cos(al))) < 1e-10
        ok_i &= abs(zz) > 1
report(f"(i) a breather passes with Y_1 -> z(theta + i alpha) Y_1, Y_2 -> z(theta - i alpha) Y_2, real on both "
       f"sides (residual {worst:.0e}); |z(u + i alpha)|^2 = (cosh u + cos alpha)/(cosh u - cos alpha) > 1", ok_i)

# ---------------------------------------------------------------------------------------------------------
# (j) linear waves
R_, T_ = sp.symbols('R T')
u_ = sp.exp(-th)                        # e^{-theta}
# incident from the left: Phi = e^{i(px - wt)} + R e^{-i(px + wt)}, Psi = T e^{i(px - wt)}
# d_+ e^{i(px-wt)} = -i m0 e^{-theta}, d_- e^{i(px-wt)} = -i m0 e^{theta};
# d_+ e^{-i(px+wt)} = -i m0 e^{theta}, d_- e^{-i(px+wt)} = -i m0 e^{-theta}
I = sp.I
eq1 = sp.Eq(-I * m0s * u_ * (1 - T_) - I * m0s / u_ * R_, -sg_s * m0s * (1 + R_ + T_))
eq2 = sp.Eq(-I * m0s / u_ * (1 + T_) - I * m0s * u_ * R_, (m0s / sg_s) * (1 + R_ - T_))
sol = sp.solve([eq1, eq2], [R_, T_], dict=True)[0]
Tform = (u_ + I * sg_s) / (u_ - I * sg_s)
ok_j = sp.simplify(sol[R_]) == 0 and sp.simplify(sol[T_] - Tform) == 0
# the plane-wave derivatives used above
pw = sp.exp(I * m0s * (xx * sp.sinh(th) - tt * sp.cosh(th)))
ok_j &= sp.simplify((sp.diff(pw, tt) + sp.diff(pw, xx)) / pw + I * m0s * sp.exp(-th)) == 0
ok_j &= sp.simplify((sp.diff(pw, tt) - sp.diff(pw, xx)) / pw + I * m0s * sp.exp(th)) == 0
pw2 = sp.exp(-I * m0s * (xx * sp.sinh(th) + tt * sp.cosh(th)))
ok_j &= sp.simplify((sp.diff(pw2, tt) + sp.diff(pw2, xx)) / pw2 + I * m0s * sp.exp(th)) == 0
ok_j &= sp.simplify((sp.diff(pw2, tt) - sp.diff(pw2, xx)) / pw2 + I * m0s * sp.exp(-th)) == 0
# the soliton at theta + i pi/2 is the plane wave
Ywave = sp.lambdify((th, tt, xx), Yexpr(th + I * sp.pi / 2, 0), 'numpy')
pwn = sp.lambdify((th, tt, xx), pw.subs(m0s, 1), 'numpy')
ok_j &= max(abs(Ywave(*v) - pwn(*v)) for v in rng.normal(size=(5, 3))) < 1e-12
Tnum = lambda theta, thd_: (cmath.exp(-theta) + 1j * math.exp(-thd_)) / (cmath.exp(-theta) - 1j * math.exp(-thd_))
Tbook = lambda theta, thd_: (cmath.sinh(thd_) + 1j * cmath.cosh(theta)) / (cmath.cosh(thd_) - 1j * cmath.sinh(theta))
for theta in np.linspace(-2.5, 2.5, 11):
    ok_j &= abs(Tnum(theta, thd) - zfun(theta + 1j * math.pi / 2, thd)) < 1e-12
    ok_j &= abs(Tnum(theta, thd) - Tbook(theta, thd)) < 1e-12
    ok_j &= abs(abs(Tnum(theta, thd)) - 1) < 1e-12
    ok_j &= abs(Tnum(theta, thd) * Tnum(theta + 1j * math.pi, thd) - 1) < 1e-12
ok_j &= abs(Tnum(thd, thd) - 1j) < 1e-12
ok_j &= abs(Tnum(0.3, 40.0) - 1) < 1e-12 and abs(Tnum(0.3, -40.0) + 1) < 1e-12
ok_j &= abs(Tnum(-40.0, thd) - 1) < 1e-12 and abs(Tnum(40.0, thd) + 1) < 1e-12
ok_j &= abs(Tnum(0.7, thd) * Tnum(-0.7, thd) - 1) > 0.1
# wave from the right: Psi = e^{-i(px + wt)} + bR e^{i(px - wt)}, Phi = bT e^{-i(px + wt)}
bR, bT = sp.symbols('bR bT')
eq1r = sp.Eq(-I * m0s / u_ * bT - (-I * m0s / u_ - I * m0s * u_ * bR), -sg_s * m0s * (bT + 1 + bR))
eq2r = sp.Eq(-I * m0s * u_ * bT + (-I * m0s * u_ - I * m0s / u_ * bR), (m0s / sg_s) * (bT - 1 - bR))
solr = sp.solve([eq1r, eq2r], [bR, bT], dict=True)[0]
ok_j &= sp.simplify(solr[bR]) == 0
bTf = sp.lambdify((th, sg_s, m0s), solr[bT], 'mpmath')
for theta in (-1.1, 0.2, 0.9, 1.6):
    v = complex(bTf(theta, sgv, 1.0))
    ok_j &= abs(v - Tnum(1j * math.pi - theta, thd)) < 1e-12
    ok_j &= abs(Tnum(theta, thd) * complex(bTf(-theta, sgv, 1.0)) - 1) < 1e-12
# pole and zero
ok_j &= abs(Tnum(thd + 1j * math.pi / 2, thd)) < 1e-12 and abs(1 / Tnum(thd - 1j * math.pi / 2 + 1e-9, thd)) < 1e-6
report("(j) linear waves: R = 0 and T = (e^-theta + i sigma)/(e^-theta - i sigma) = z(theta - theta_d + i pi/2), "
       "the book's T^cl; |T| = 1, T(theta_d) = i, limits; from the right bar R = 0, bar T(theta) = T(i pi - theta), "
       "unitarity and crossing, no parity; the soliton at theta + i pi/2 is the plane wave; pole and zero", ok_j)

# ---------------------------------------------------------------------------------------------------------
# (k) one function on the complex plane
ok_k = True
for re in np.linspace(-3, 3, 13):
    for im, rel in ((math.pi / 2, 0), (-math.pi / 2, 0), (0.4, 1), (-1.2, 1), (2.0, -1), (2.9, -1)):
        zz = -1 / cmath.tanh((re + 1j * im) / 2)
        if rel == 0:
            ok_k &= abs(abs(zz) - 1) < 1e-12
        elif rel == 1:
            ok_k &= abs(zz) > 1
        else:
            ok_k &= abs(zz) < 1
        ok_k &= abs(-1 / cmath.tanh((re + 1j * im + 1j * math.pi) / 2) - 1 / zz) < 1e-12
        ok_k &= abs(1 / zz + cmath.tanh((re + 1j * im) / 2)) < 1e-12
# odd-sector solution: (Phi, 2 pi + 4 atan(z_odd Y)) with z_odd = 1/z
worst = 0.0
for theta in (0.3, 1.7):
    zo = 1 / zfun(theta, thd).real
    L3 = lamb(Yexpr(theta, 0.2))
    R3n = lamb(zo * Yexpr(theta, 0.2))
    R3 = [lambda tv, xv, f=R3n[0]: f(tv, xv) + 2 * math.pi, R3n[1], R3n[2]]
    worst = max(worst, max(residual(L3, R3, sgv, tv) for tv in (-2.0, 0.0, 1.5)))
ok_k &= worst < 1e-12
report(f"(k) |z(u)| = 1 on Im u = +-pi/2, > 1 inside, < 1 for pi/2 < Im u < pi; z(u + i pi) = 1/z(u) = -tanh(u/2); "
       f"the odd sector passes solitons with z_odd = 1/z (residual {worst:.0e})", ok_k)

# ---------------------------------------------------------------------------------------------------------
# (l) half a soliton
ok_l = True
for u in np.linspace(-3, 3, 13):
    if abs(u) < 1e-9:
        continue
    zz = -1 / math.tanh(u / 2)
    A = math.tanh(u / 2) ** 2
    ok_l &= abs((1 / zz) ** 2 - A) < 1e-12 and abs(zz**2 - 1 / A) < 1e-9 * zz**2
    Todd = (cmath.exp(-u) - 1j) / (cmath.exp(-u) + 1j)       # sigma -> -sigma at theta_d = 0
    Tev = (cmath.exp(-u) + 1j) / (cmath.exp(-u) - 1j)
    Tk = (math.sinh(u) + 1j) / (math.sinh(u) - 1j)
    ok_l &= abs(Todd**2 - Tk) < 1e-12 and abs(Tev**2 - 1 / Tk) < 1e-12
    ok_l &= abs(Todd - Tev.conjugate()) < 1e-12
# the transmission factor of a static kink, from the exact fluctuation
k_ = sp.Symbol('k', positive=True)
fl = (sp.tanh(m0s * xx) - I * k_ / m0s) * sp.exp(I * k_ * xx)
cosK = 1 - 2 / sp.cosh(m0s * xx) ** 2                     # cos(4 arctan e^{m0 x})
ode = -sp.diff(fl, xx, 2) + m0s**2 * cosK * fl - (k_**2 + m0s**2) * fl
ok_l &= sp.simplify(ode.rewrite(sp.exp)) == 0
ok_l &= abs(math.cos(4 * math.atan(math.exp(0.37))) - (1 - 2 / math.cosh(0.37) ** 2)) < 1e-14
Tk_from = ((1 - I * k_ / m0s) / (-1 - I * k_ / m0s))
for u in (-1.3, 0.4, 2.2):
    kv = math.sinh(u)
    ok_l &= abs(complex(Tk_from.subs({k_: kv, m0s: 1})) - (math.sinh(u) + 1j) / (math.sinh(u) - 1j)) < 1e-12
# energy with BCZ's constants
for thd_v in (-0.7, 0.4):
    sgd = math.exp(-thd_v)
    Bo, Po = defect_BP(0.0, 2 * math.pi, sgd)
    Be, Pe = defect_BP(0.0, 0.0, sgd)
    shift = (M / 2 * math.cosh(thd_v), M / 2 * math.sinh(thd_v))
    ok_l &= abs(Bo - shift[0] - M / 2 * math.cosh(thd_v)) < 1e-12 and abs(Po - shift[1] - M / 2 * math.sinh(thd_v)) < 1e-12
    ok_l &= abs((Be - shift[0]) - M / 2 * cmath.cosh(thd_v + 1j * math.pi)) < 1e-12
    ok_l &= abs((Pe - shift[1]) - M / 2 * cmath.sinh(thd_v + 1j * math.pi)) < 1e-12
report("(l) z_odd^2 = A(theta - theta_d), z^2 = 1/A; T_odd^2 = (sinh u + i)/(sinh u - i), the transmission of a wave "
       "through a kink (exact fluctuation), T^2 its inverse; with BCZ's constants odd and even states have "
       "(M/2)(cosh, sinh) of theta_d and theta_d + i pi", ok_l)

# ---------------------------------------------------------------------------------------------------------
# (m) two defects
ok_m = True
Lsep = 2.5
thA, thB = 0.9, -0.4
sA, sB = math.exp(-thA), math.exp(-thB)
worst = 0.0
for theta in (0.3, 1.2):
    zA, zB = zfun(theta, thA).real, zfun(theta, thB).real
    L3 = lamb(Yexpr(theta, 0.1))
    M3 = lamb(zA * Yexpr(theta, 0.1))
    R3 = lamb(zA * zB * Yexpr(theta, 0.1))
    for tv in (-1.0, 0.5, 3.0):
        worst = max(worst, residual(L3, M3, sA, tv, 0.0), residual(M3, R3, sB, tv, Lsep))
ok_m &= worst < 1e-12
for theta in np.linspace(-2, 2, 9):
    if abs(theta - thA) > 1e-6:
        ok_m &= zfun(theta, thA).real ** 2 > 0
report(f"(m) two defects in a row multiply the factors: Psi_1 = 4 arctan(z_1 Y), Psi_2 = 4 arctan(z_1 z_2 Y) "
       f"(residual {worst:.0e}); identical defects never convert", ok_m)

# ---------------------------------------------------------------------------------------------------------
# (n) exercises
ok_n = True
# type I_+: Phi_x = -Psi_t - B+_Phi, Psi_x = -Phi_t + B+_Psi with B+ = B(sigma <-> 1/sigma); solution 4 atan(-zY)
worst = 0.0
for theta in (0.3, 1.7):
    zp_ = -zfun(theta, thd).real
    L3 = lamb(Yexpr(theta, 0.2))
    R3 = lamb(zp_ * Yexpr(theta, 0.2))
    for tv in (-1.5, 0.0, 2.0):
        ph0, pt_, px_ = (float(f(tv, 0.0)) for f in L3)
        ps0, qt_, qx_ = (float(f(tv, 0.0)) for f in R3)
        sp_ = 1 / sgv                                       # sigma and 1/sigma exchanged in B
        Bp = m0 * (sp_ * math.sin((ph0 + ps0) / 2) + math.sin((ph0 - ps0) / 2) / sp_)
        Bq = m0 * (sp_ * math.sin((ph0 + ps0) / 2) - math.sin((ph0 - ps0) / 2) / sp_)
        worst = max(worst, abs(px_ - (-qt_ - Bp)), abs(qx_ - (-pt_ + Bq)))
# antisoliton: Y -> -Y on both sides, same z
for theta in (0.3, 1.7):
    zv = zfun(theta, thd).real
    L3 = lamb(-Yexpr(theta, 0.2))
    R3 = lamb(-zv * Yexpr(theta, 0.2))
    worst = max(worst, max(residual(L3, R3, sgv, tv) for tv in (-1.0, 0.0, 1.0)))
ok_n &= worst < 1e-12
# soliton from the right, theta < 0: incoming on the right has factor z, outgoing on the left factor 1: relative 1/z
for theta in (-0.3, -1.4):
    for thd_v in (0.9, -0.8):
        zv = zfun(theta, thd_v).real
        dx_out = math.log(abs(zv)) / (m0 * math.cosh(theta))   # outgoing (left) minus incoming (right) centre
        ok_n &= dx_out > 0                                     # behind, for a left mover: delayed
        ok_n &= (zv < 0) == (theta > thd_v)
# Wigner: M sinh theta * Delta t = (8/beta^2) ln|z|
theta = 0.4
ok_n &= abs(M * math.sinh(theta) * math.log(abs(zfun(theta, thd))) / (m0 * math.sinh(theta))
            - 8 * math.log(abs(zfun(theta, thd)))) < 1e-12
# NOT gate: theta_d < 0 converts every soliton from the left
ok_n &= all(zfun(theta, -0.3).real < 0 for theta in np.linspace(0.01, 5, 50))
report(f"(n) type I_+ passes with -z, an antisoliton with z (residual {worst:.0e}); a soliton from the right is "
       f"delayed and converted iff theta > theta_d; Wigner's phase slope (8/beta^2) ln|z|; theta_d < 0 is a NOT gate",
       ok_n)

# ---------------------------------------------------------------------------------------------------------
# (o) odd sector near theta_d: emission before absorption
ok_o = True
for eps in (1e-3, -1e-3):
    theta = thd_g + eps
    zo = 1 / zfun(theta, thd_g).real
    Phi0 = lambda tv: 4 * math.atan(math.exp(-m0 * tv * math.sinh(theta)))
    Psi0 = lambda tv: 2 * math.pi + 4 * math.atan(zo * math.exp(-m0 * tv * math.sinh(theta)))
    t1 = brentq(lambda tv: Phi0(tv) - math.pi, -60, 60, xtol=1e-14)
    tgt = 3 * math.pi if zo > 0 else math.pi
    t2 = brentq(lambda tv: Psi0(tv) - tgt, -80, 60, xtol=1e-14)
    ok_o &= t2 < t1 and abs((t1 - t2) - math.log(1 / abs(zo)) / (m0 * math.sinh(theta))) < 1e-10
    # between the two times the jump is even: B close to 0
    tm = 0.5 * (t1 + t2)
    ok_o &= abs(defect_BP(Phi0(tm), Psi0(tm), sgg)[0]) < 0.05 * M
report("(o) odd sector, theta near theta_d: the right defect value passes half-way first, by ln(1/|z_odd|)/(m0 sinh "
       "theta), and in between the jump is even: the stored soliton leaves before the incoming one is swallowed", ok_o)

raise SystemExit(0 if OK else 1)
