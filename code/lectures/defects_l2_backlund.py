"""Checks for Lecture 2 of the defect lecture notes, "Bäcklund's machine"
(lecture-notes/defects/02-backlund-machine.qmd).

Conventions of the lecture: Phi = beta_SG phi, d_pm = d_t pm d_x, sine-Gordon d_+ d_- Phi = -m0^2 sin Phi;
asymptotic coordinates u = m0 (x+t)/2, v = m0 (x-t)/2; Baecklund transformation Psi = B_a Phi:
    d_+ (Phi - Psi) = 2 a m0 sin((Phi+Psi)/2),   d_- (Phi + Psi) = -(2 m0/a) sin((Phi-Psi)/2).
Written from these definitions, independently of code/defect_classical_sg.py.
One PASS or FAIL line per check; exit status 1 on any FAIL.

  (a) Chebyshev nets: the metric du^2 + 2 cos(w) du dv + dv^2 has Gaussian curvature -w_uv/sin(w) (Brioschi).
  (b) The dictionary: d_+ = m0 d_u, d_- = -m0 d_v, so sine-Gordon is Bour's equation w_uv = sin w with w = Phi.
  (c) The static kink is the pseudosphere: first and second fundamental forms of
      (sech s cos t, sech s sin t, s - tanh s), s = u+v, t = u-v, are those of w = 4 arctan e^{u+v};
      the tractrix has unit tangent length; the kink energy density is (4 m0^2/beta^2) times the squared radius.
  (d) The Baecklund system is consistent exactly on solutions: cross-differentiation gives sine-Gordon for both.
  (e) From the vacuum: d_+ ln Y = -a m0, d_- ln Y = m0/a; a = -e^{-theta} gives the kink of rapidity theta;
      a = +e^{-theta} gives the antikink on the vacuum 2 pi; B_{-a} undoes B_a.
  (f) Boosts: a boost by rapidity vartheta maps a Baecklund pair with parameter a to one with a e^{-vartheta}.
  (g) Bianchi permutability: the algebraic formula follows from the four d_+ equations; the Phi_12 it defines
      satisfies all four first-order systems; from the vacuum it is the two-soliton solution with
      A = tanh^2((theta1-theta2)/2).
  (h) The collision shift: soliton 1 is displaced forward by -ln A/(m0 cosh theta1), measured on the exact solution.
  (i) Complex-conjugate parameters give the breather 4 arctan(cot(alpha) sin(m0 t sin alpha)/cosh(m0 x cos alpha)),
      a solution of sine-Gordon, with energy 2M cos(alpha).
  (j) The geometric construction: r' = r + (sin s/sin w)[sin((w-w')/2) r_u + sin((w+w')/2) r_v], tan(s/2) = |a|,
      maps the pseudosphere to a surface of curvature -1 with asymptotic angle w' = B_a w; the segment has
      length sin(sigma), is tangent to both surfaces, and the normals meet at the angle sigma;
      cos(sigma) = tanh(theta) and sin(sigma) = 1/cosh(theta).
  (k) Pushing the second soliton to x = +infinity in Bianchi's formula leaves the first one displaced by the
      Bianchi coefficient (a1+a2)/(a1-a2) (the delay of Lecture 4).
  (l) Bianchi's cube (exercise): Phi_123 from Phi_1, Phi_12, Phi_13 solves sine-Gordon.

Runtime: about a minute (the geometric check uses 30-digit numerical derivatives).
"""
import math

import mpmath as mp
import numpy as np
import sympy as sp

OK = True


def report(name, ok):
    global OK
    OK = OK and bool(ok)
    print(("PASS " if ok else "FAIL ") + name)


rng = np.random.default_rng(2)
t, x, u, v = sp.symbols('t x u v', real=True)
m0 = sp.Symbol('m0', positive=True)
a = sp.Symbol('a', nonzero=True)


def dp(f):
    return sp.diff(f, t) + sp.diff(f, x)


def dm(f):
    return sp.diff(f, t) - sp.diff(f, x)


def is_zero(expr, subs_list):
    for s in subs_list:
        if abs(complex(sp.N(expr.subs(s)))) > 1e-9:
            return False
    return True


def rand_pts(n=4, extra=None):
    pts = []
    for _ in range(n):
        d = {t: rng.normal(), x: rng.normal(), m0: 0.5 + rng.random()}
        if extra:
            d.update(extra())
        pts.append(d)
    return pts


# ---------------------------------------------------------------------------------------------
# (a) Chebyshev nets
# ---------------------------------------------------------------------------------------------
w = sp.Function('w')(u, v)
E, F, G = sp.Integer(1), sp.cos(w), sp.Integer(1)
Eu, Ev = sp.diff(E, u), sp.diff(E, v)
Fu, Fv = sp.diff(F, u), sp.diff(F, v)
Gu, Gv = sp.diff(G, u), sp.diff(G, v)
Evv, Euv, Guu = sp.diff(E, v, 2), sp.diff(F, u, v), sp.diff(G, u, 2)
M1 = sp.Matrix([[-Evv / 2 + Fu.diff(v) - Guu / 2, Eu / 2, Fu - Ev / 2],
                [Fv - Gu / 2, E, F],
                [Gv / 2, F, G]])
M2 = sp.Matrix([[0, Ev / 2, Gu / 2], [Ev / 2, E, F], [Gu / 2, F, G]])
K_brioschi = sp.simplify((M1.det() - M2.det()) / (E * G - F ** 2) ** 2)
report("(a) Chebyshev net du^2+2cos(w)dudv+dv^2 has curvature K = -w_uv/sin w",
       sp.simplify(K_brioschi + sp.diff(w, u, v) / sp.sin(w)) == 0)

# ---------------------------------------------------------------------------------------------
# (b) the dictionary to field theory
# ---------------------------------------------------------------------------------------------
W = sp.Function('W')
uu, vv = m0 * (x + t) / 2, m0 * (x - t) / 2
Phi_uv = W(uu, vv)
dpW = sp.simplify(dp(Phi_uv) - m0 * sp.Subs(sp.Derivative(W(u, v), u), (u, v), (uu, vv)).doit())
dmW = sp.simplify(dm(Phi_uv) + m0 * sp.Subs(sp.Derivative(W(u, v), v), (u, v), (uu, vv)).doit())
lhs = sp.simplify(dp(dm(Phi_uv)) + m0 ** 2 * sp.Subs(sp.Derivative(W(u, v), u, v), (u, v), (uu, vv)).doit())
report("(b) d_+ = m0 d_u, d_- = -m0 d_v, d_+d_- = -m0^2 d_u d_v", dpW == 0 and dmW == 0 and lhs == 0)

# ---------------------------------------------------------------------------------------------
# (c) the pseudosphere
# ---------------------------------------------------------------------------------------------
s_, tau_ = u + v, u - v
r = sp.Matrix([sp.sech(s_) * sp.cos(tau_), sp.sech(s_) * sp.sin(tau_), s_ - sp.tanh(s_)])
om = 4 * sp.atan(sp.exp(s_))
ru, rv = r.diff(u), r.diff(v)
N = ru.cross(rv)
N = N / sp.sqrt(N.dot(N))
I_E, I_F, I_G = ru.dot(ru), ru.dot(rv), rv.dot(rv)
II_L, II_M, II_N = r.diff(u, 2).dot(N), r.diff(u, v).dot(N), r.diff(v, 2).dot(N)
pts = [{u: rng.normal(), v: rng.normal()} for _ in range(5)]
pts = [p for p in pts if abs(p[u] + p[v]) > 0.05]      # away from the rim, where N flips
ok = is_zero(I_E - 1, pts) and is_zero(I_G - 1, pts) and is_zero(I_F - sp.cos(om), pts)
ok &= is_zero(II_L, pts) and is_zero(II_N, pts) and is_zero(sp.Abs(II_M) - sp.Abs(sp.sin(om)), pts)
report("(c) pseudosphere: I = du^2+2cos(w)dudv+dv^2, II = +-2 sin(w) dudv with w = 4 arctan e^{u+v}", ok)
sv = sp.Symbol('s', positive=True)
rr, zz = sp.sech(sv), sv - sp.tanh(sv)
tang = sp.Matrix([sp.diff(rr, sv), sp.diff(zz, sv)])
lam = rr / (-tang[0])                    # move along the tangent until the radius vanishes
seg = sp.simplify(sp.sqrt((lam * tang[0]) ** 2 + (lam * tang[1]) ** 2))
report("(c) the tractrix: the tangent segment from the curve to the axis has length 1",
       is_zero(seg - 1, [{sv: 0.3}, {sv: 1.7}, {sv: 3.1}]))
beta = sp.Symbol('beta', positive=True)
kink = 4 * sp.atan(sp.exp(m0 * x))
dens = (sp.Rational(1, 2) * sp.diff(kink, x) ** 2 + m0 ** 2 * (1 - sp.cos(kink))) / beta ** 2
ok = is_zero(dens - 4 * m0 ** 2 / beta ** 2 * sp.sech(m0 * x) ** 2, [{x: 0.4, m0: 1.3, beta: 0.7}, {x: -1.1, m0: 0.6, beta: 1.9}])
# int (4 m0^2/beta^2) sech^2(m0 x) dx = (4 m0/beta^2) [tanh] = 8 m0/beta^2
ok &= sp.simplify(sp.diff(sp.tanh(sv), sv) - sp.sech(sv) ** 2) == 0
ok &= sp.limit(4 * sp.tanh(sv), sv, sp.oo) - sp.limit(4 * sp.tanh(-sv), sv, sp.oo) == 8
report("(c) kink energy density = (4 m0^2/beta^2) sech^2(m0 x) = (4m0^2/beta^2) radius^2; mass 8 m0/beta^2", ok)

# ---------------------------------------------------------------------------------------------
# (d) the Baecklund system is consistent exactly on solutions
# ---------------------------------------------------------------------------------------------
P, Q = sp.Function('P')(t, x), sp.Function('Q')(t, x)
bt1 = 2 * a * m0 * sp.sin((P + Q) / 2)        # = d_+(P - Q)
bt2 = -2 * m0 / a * sp.sin((P - Q) / 2)       # = d_-(P + Q)
# d_- of the first, using the second for d_-(P+Q); d_+ of the second, using the first for d_+(P-Q)
A_, B_ = (P + Q) / 2, (P - Q) / 2
dm_first = 2 * a * m0 * sp.cos(A_) * bt2 / 2          # d_+d_-(P-Q)
dp_second = -2 * m0 / a * sp.cos(B_) * bt1 / 2        # d_+d_-(P+Q)
ddP = sp.simplify((dm_first + dp_second) / 2)
ddQ = sp.simplify((dp_second - dm_first) / 2)
report("(d) cross-differentiation: d_+d_- P = -m0^2 sin P and d_+d_- Q = -m0^2 sin Q",
       sp.simplify(sp.expand_trig(ddP + m0 ** 2 * sp.sin(P))) == 0
       and sp.simplify(sp.expand_trig(ddQ + m0 ** 2 * sp.sin(Q))) == 0)

# ---------------------------------------------------------------------------------------------
# (e) from the vacuum
# ---------------------------------------------------------------------------------------------
th = sp.Symbol('theta', real=True)
c0 = sp.Symbol('c', real=True)


def soliton_from_vacuum(aval):
    # d_+ lnY = -a m0, d_- lnY = m0/a  =>  lnY = A x + B t + c
    Ax = -m0 / 2 * (aval + 1 / aval)
    Bt = m0 / 2 * (1 / aval - aval)
    return 4 * sp.atan(sp.exp(Ax * x + Bt * t + c0))


def bt_residual(Pf, Qf, aval):
    r1 = dp(Pf - Qf) - 2 * aval * m0 * sp.sin((Pf + Qf) / 2)
    r2 = dm(Pf + Qf) + 2 * m0 / aval * sp.sin((Pf - Qf) / 2)
    return r1, r2


def sg_residual(f):
    return dp(dm(f)) + m0 ** 2 * sp.sin(f)


Psi = soliton_from_vacuum(-sp.exp(-th))
pts = rand_pts(extra=lambda: {th: rng.normal(), c0: rng.normal()})
r1, r2 = bt_residual(sp.Integer(0), Psi, -sp.exp(-th))
ok = is_zero(r1, pts) and is_zero(r2, pts) and is_zero(sg_residual(Psi), pts)
kink_th = 4 * sp.atan(sp.exp(m0 * (x * sp.cosh(th) - t * sp.sinh(th)) + c0))
ok &= is_zero(Psi - kink_th, pts)
report("(e) a = -e^{-theta}: B_a(0) = 4 arctan exp(m0(x cosh(theta) - t sinh(theta)) + c), a kink of rapidity theta", ok)
Psi2 = soliton_from_vacuum(sp.exp(-th))
r1, r2 = bt_residual(sp.Integer(0), Psi2, sp.exp(-th))
ok = is_zero(r1, pts) and is_zero(r2, pts)
ok &= is_zero(Psi2 - (2 * sp.pi - 4 * sp.atan(sp.exp(m0 * (x * sp.cosh(th) - t * sp.sinh(th)) - c0))), pts)
report("(e) a = +e^{-theta}: an antikink of rapidity theta on the vacuum 2 pi", ok)
# B_{-a} undoes B_a: if (P,Q) satisfy B_a then (Q,P) satisfy B_{-a}
r1a, r2a = bt_residual(Psi, sp.Integer(0), sp.exp(-th))
report("(e) (Psi, 0) is a Baecklund pair with parameter -a: B_{-a} undoes B_a",
       is_zero(r1a, pts) and is_zero(r2a, pts))

# ---------------------------------------------------------------------------------------------
# (f) boosts rescale the parameter
# ---------------------------------------------------------------------------------------------
vt = sp.Symbol('vartheta', real=True)
boost = {t: t * sp.cosh(vt) - x * sp.sinh(vt), x: x * sp.cosh(vt) - t * sp.sinh(vt)}
Pb = Psi.subs(boost, simultaneous=True)         # the kink boosted by vartheta
Zb = sp.Integer(0)
r1, r2 = bt_residual(Zb, Pb, -sp.exp(-th) * sp.exp(-vt))
pts = rand_pts(extra=lambda: {th: rng.normal(), c0: rng.normal(), vt: rng.normal()})
ok = is_zero(r1, pts) and is_zero(r2, pts)
ok &= is_zero(Pb - kink_th.subs(th, th + vt), pts)
report("(f) a boost by vartheta maps a Baecklund pair with a to one with a e^{-vartheta}; theta -> theta + vartheta", ok)

# ---------------------------------------------------------------------------------------------
# (g) Bianchi permutability
# ---------------------------------------------------------------------------------------------
P0, P1, P2, P12 = sp.symbols('P0 P1 P2 P12', real=True)
a1s, a2s = sp.symbols('a1 a2', nonzero=True)
# sum of the d_+ equations along the two paths, equated: a1[s(P0+P1) - s(P2+P12)] = a2[s(P0+P2) - s(P1+P12)]
lhs = a1s * (sp.sin((P0 + P1) / 2) - sp.sin((P2 + P12) / 2)) - a2s * (sp.sin((P0 + P2) / 2) - sp.sin((P1 + P12) / 2))
Xq, Yq = (P0 - P12) / 4, (P1 - P2) / 4
fact = 2 * sp.cos((P0 + P1 + P2 + P12) / 4) * (a1s * sp.sin(Xq + Yq) - a2s * sp.sin(Xq - Yq))
pts = [{P0: rng.normal(), P1: rng.normal(), P2: rng.normal(), P12: rng.normal(), a1s: rng.normal(), a2s: rng.normal()}
       for _ in range(6)]
ok = is_zero(lhs - fact, pts)
# a1 sin(X+Y) = a2 sin(X-Y)  <=>  tan X = -(a1+a2)/(a1-a2) tan Y, i.e. tan((P12-P0)/4) = (a1+a2)/(a1-a2) tan((P1-P2)/4)
for pt in pts:
    Yv = complex(sp.N(Yq.subs(pt))).real
    a1v, a2v = pt[a1s], pt[a2s]
    Xv = math.atan(-(a1v + a2v) / (a1v - a2v) * math.tan(Yv))
    ok &= abs(a1v * math.sin(Xv + Yv) - a2v * math.sin(Xv - Yv)) < 1e-12
report("(g) the four d_+ equations give a1 sin(X+Y) = a2 sin(X-Y), i.e. tan((P12-P0)/4) = (a1+a2)/(a1-a2) tan((P1-P2)/4)", ok)


def num_solution(a1, a2, xi1, xi2):
    """Phi_1, Phi_2 = B_{a_k}(0) with complex constants, and Phi_12 by Bianchi; as numpy-callable functions."""
    def Yk(ak, xik, T, X):
        Ax = -0.5 * (ak + 1 / ak)
        Bt = 0.5 * (1 / ak - ak)
        return np.exp(Ax * X + Bt * T + xik)
    def phis(T, X):
        Y1, Y2 = Yk(a1, xi1, T, X), Yk(a2, xi2, T, X)
        c = (a1 + a2) / (a1 - a2)
        P1n, P2n = 4 * np.arctan(Y1), 4 * np.arctan(Y2)
        P12n = 4 * np.arctan(c * (Y1 - Y2) / (1 + Y1 * Y2))
        return P1n, P2n, P12n
    return phis


def check_bt_numeric(Pf, Qf, aval, T, X, h=1e-4):
    """Residuals of B_a for numerically given P, Q (functions of T, X), m0 = 1, central differences."""
    def d_pm(f, sgn):
        return (f(T + h, X) - f(T - h, X)) / (2 * h) + sgn * (f(T, X + h) - f(T, X - h)) / (2 * h)
    r1 = d_pm(lambda a_, b_: Pf(a_, b_) - Qf(a_, b_), 1) - 2 * aval * np.sin((Pf(T, X) + Qf(T, X)) / 2)
    r2 = d_pm(lambda a_, b_: Pf(a_, b_) + Qf(a_, b_), -1) + 2 / aval * np.sin((Pf(T, X) - Qf(T, X)) / 2)
    return max(abs(r1), abs(r2))


th1n, th2n = 0.9, -0.4
a1n, a2n = -math.exp(-th1n), -math.exp(-th2n)
phis = num_solution(a1n, a2n, 0.3, -0.2)
ok = True
for _ in range(6):
    T, X = rng.normal(), rng.normal()
    P1f = lambda a_, b_: phis(a_, b_)[0]
    P2f = lambda a_, b_: phis(a_, b_)[1]
    P12f = lambda a_, b_: phis(a_, b_)[2]
    zero = lambda a_, b_: 0.0
    ok &= check_bt_numeric(zero, P1f, a1n, T, X) < 1e-6
    ok &= check_bt_numeric(zero, P2f, a2n, T, X) < 1e-6
    ok &= check_bt_numeric(P1f, P12f, a2n, T, X) < 1e-6
    ok &= check_bt_numeric(P2f, P12f, a1n, T, X) < 1e-6
report("(g) Phi_12 from the Bianchi formula satisfies all four first-order systems (B_a1, B_a2 on both sides)", ok)
# compare with the Hirota form tan(Phi/4) = (X1+X2)/(1 - A X1 X2), A = tanh^2(theta12/2)
tt = math.tanh((th1n - th2n) / 2)
ok = True
for _ in range(6):
    T, X = rng.normal() * 2, rng.normal() * 2
    Y1 = math.exp(-0.5 * (a1n + 1 / a1n) * X + 0.5 * (1 / a1n - a1n) * T + 0.3)
    Y2 = math.exp(-0.5 * (a2n + 1 / a2n) * X + 0.5 * (1 / a2n - a2n) * T - 0.2)
    X1, X2 = -Y1 / tt, Y2 / tt
    hir = (X1 + X2) / (1 - tt ** 2 * X1 * X2)
    bia = ((a1n + a2n) / (a1n - a2n)) * (Y1 - Y2) / (1 + Y1 * Y2)
    ok &= abs(hir - bia) < 1e-12 * (1 + abs(bia))
ok &= abs(((a1n - a2n) / (a1n + a2n)) ** 2 - tt ** 2) < 1e-14
report("(g) from the vacuum: tan(Phi_12/4) = (X1+X2)/(1 - A X1 X2) with A = ((a1-a2)/(a1+a2))^2 = tanh^2(theta12/2)", ok)

# topological charge of Phi_12: with Y1, Y2 > 0 it is a soliton and an antisoliton (charge 0); flipping the sign
# of Y2 (xi2 -> xi2 + i pi) gives two solitons (charge 2).
Xs = np.linspace(-40, 40, 200001)
def charge(sign2):
    Y1 = np.exp(-0.5 * (a1n + 1 / a1n) * Xs + 0.3)
    Y2 = sign2 * np.exp(-0.5 * (a2n + 1 / a2n) * Xs - 0.2)
    cc = (a1n + a2n) / (a1n - a2n)
    ph = np.unwrap(4 * np.arctan2(cc * (Y1 - Y2), 1 + Y1 * Y2), period=2 * np.pi)
    return (ph[-1] - ph[0]) / (2 * np.pi)
q_pos, q_neg = charge(1.0), charge(-1.0)
report(f"(g) charge of Phi_12 from the vacuum: {q_pos:+.3f} with Y1, Y2 > 0, {q_neg:+.3f} with Y2 -> -Y2",
       abs(q_pos) < 1e-6 and abs(abs(q_neg) - 2) < 1e-6)

# ---------------------------------------------------------------------------------------------
# (h) the collision shift
# ---------------------------------------------------------------------------------------------
# two kinks, theta1 = 0.9 overtaking theta2 = -0.4 (head on).  Track the zero of the soliton-1 profile:
# far from soliton 2, Phi_12 near soliton 1 is a kink centred where |X1| = 1 (before) or A|X1| = 1 (after).
A = tt ** 2
def centre_of_1(T):
    # position where the soliton-1 factor crosses the midpoint: solve |X1| (or A|X1|) = 1 via the exact solution
    Xs = np.linspace(-40, 40, 400001)
    Y1 = np.exp(-0.5 * (a1n + 1 / a1n) * Xs + 0.5 * (1 / a1n - a1n) * T + 0.3)
    Y2 = np.exp(-0.5 * (a2n + 1 / a2n) * Xs + 0.5 * (1 / a2n - a2n) * T - 0.2)
    c = (a1n + a2n) / (a1n - a2n)
    val = np.unwrap(4 * np.arctan2(c * (Y1 - Y2), 1 + Y1 * Y2) * 1.0)
    # the gradient peak of the field away from soliton 2 marks soliton 1
    grad = np.abs(np.gradient(val, Xs))
    x2 = (-0.2 - 0.5 * (1 / a2n - a2n) * T) / (-0.5 * (a2n + 1 / a2n))   # centre of soliton 2: Y2 = 1
    mask = np.abs(Xs - x2) > 8
    return Xs[mask][np.argmax(grad[mask])]
v1 = math.tanh(th1n)
Tin, Tout = -25.0, 25.0
x_in, x_out = centre_of_1(Tin), centre_of_1(Tout)
shift = (x_out - v1 * Tout) - (x_in - v1 * Tin)
predict = -math.log(A) / math.cosh(th1n)
report(f"(h) soliton 1 is displaced forward by -ln A/(m0 cosh theta1): measured {shift:.4f}, predicted {predict:.4f}",
       abs(shift - predict) < 2e-3)

# ---------------------------------------------------------------------------------------------
# (i) the breather
# ---------------------------------------------------------------------------------------------
al = sp.Symbol('alpha', positive=True)
a1c = -sp.exp(-sp.I * al)
a2c = -sp.exp(sp.I * al)
Y1c = sp.exp(-(a1c + 1 / a1c) / 2 * m0 * x + (1 / a1c - a1c) / 2 * m0 * t)
Y2c = sp.exp(-(a2c + 1 / a2c) / 2 * m0 * x + (1 / a2c - a2c) / 2 * m0 * t)
tan_q = (a1c + a2c) / (a1c - a2c) * (Y1c - Y2c) / (1 + Y1c * Y2c)
breather_tan = sp.cot(al) * sp.sin(m0 * t * sp.sin(al)) / sp.cosh(m0 * x * sp.cos(al))
pts = [{t: rng.normal(), x: rng.normal(), m0: 0.5 + rng.random(), al: 0.2 + 1.2 * rng.random()} for _ in range(5)]
ok = is_zero(tan_q - breather_tan, pts)
br = 4 * sp.atan(breather_tan)
ok &= is_zero(sg_residual(br), pts)
report("(i) a2 = conj(a1) = -e^{i alpha}: Bianchi gives the breather 4 arctan(cot a sin(m0 t sin a)/cosh(m0 x cos a))", ok)
# energy of the breather: (1/beta^2) int [1/2 Phi_t^2 + 1/2 Phi_x^2 + m0^2 (1 - cos Phi)] dx = 2 M cos(alpha), M = 8 m0/beta^2
def breather_energy(alpha, m=1.0, T=0.37):
    Xs = np.linspace(-60, 60, 120001)
    sa, ca = math.sin(alpha), math.cos(alpha)
    def f(T_, X_):
        return 4 * np.arctan(ca / sa * np.sin(m * T_ * sa) / np.cosh(m * X_ * ca))
    h = 1e-5
    Pt = (f(T + h, Xs) - f(T - h, Xs)) / (2 * h)
    Px = np.gradient(f(T, Xs), Xs)
    dens = 0.5 * Pt ** 2 + 0.5 * Px ** 2 + m ** 2 * (1 - np.cos(f(T, Xs)))
    return np.trapz(dens, Xs)
ok = all(abs(breather_energy(al_) - 2 * 8 * math.cos(al_)) < 1e-3 for al_ in (0.4, 0.9, 1.3))
report("(i) breather energy = 2 M cos(alpha) (beta = m0 = 1)", ok)

# ---------------------------------------------------------------------------------------------
# (j) the geometric construction of Baecklund
# ---------------------------------------------------------------------------------------------
mp.mp.dps = 30
thg = mp.mpf('0.6')
ag1, ag2 = mp.mpf(-1), -mp.e ** (-thg)
cg = (ag1 + ag2) / (ag1 - ag2)


def Yg(ak, U, V, xi=0):
    return mp.e ** (-ak * U - V / ak + xi)          # d_u lnY = -a, d_v lnY = -1/a  (m0 = 1)


def omg(U, V):
    return 4 * mp.atan(Yg(ag1, U, V))


def omgp(U, V):
    Pa, Pb = 4 * mp.atan(Yg(ag1, U, V)), 4 * mp.atan(Yg(ag2, U, V, mp.mpf(1) / 3))
    return 4 * mp.atan(cg * mp.tan((Pa - Pb) / 4))


def rg(U, V):
    s, ta = U + V, U - V
    return mp.matrix([mp.sech(s) * mp.cos(ta), mp.sech(s) * mp.sin(ta), s - mp.tanh(s)])


def dvec(fun, U, V, n):
    return mp.matrix([mp.diff(lambda p_, q_: fun(p_, q_)[i], (U, V), n) for i in range(3)])


def cross(p_, q_):
    return mp.matrix([p_[1] * q_[2] - p_[2] * q_[1], p_[2] * q_[0] - p_[0] * q_[2], p_[0] * q_[1] - p_[1] * q_[0]])


def dot(p_, q_):
    return sum(p_[i] * q_[i] for i in range(3))


def forms(fun, U, V):
    Ru, Rv = dvec(fun, U, V, (1, 0)), dvec(fun, U, V, (0, 1))
    Nn = cross(Ru, Rv)
    Nn = Nn / mp.sqrt(dot(Nn, Nn))
    return (dot(Ru, Ru), dot(Ru, Rv), dot(Rv, Rv), dot(dvec(fun, U, V, (2, 0)), Nn),
            dot(dvec(fun, U, V, (1, 1)), Nn), dot(dvec(fun, U, V, (0, 2)), Nn), Nn)


sig = 2 * mp.atan(abs(ag2))


def rgp(U, V):
    w_, wp_ = omg(U, V), omgp(U, V)
    Ru, Rv = dvec(rg, U, V, (1, 0)), dvec(rg, U, V, (0, 1))
    return rg(U, V) + (mp.sin(sig) / mp.sin(w_)) * (mp.sin((w_ - wp_) / 2) * Ru + mp.sin((w_ + wp_) / 2) * Rv)


ok = True
for (U0, V0) in ((mp.mpf('0.37'), mp.mpf('-0.81')), (mp.mpf('-0.6'), mp.mpf('0.25'))):
    # Baecklund pair (omega -> omega') in u, v: ((w-w')/2)_u = a sin((w+w')/2), ((w+w')/2)_v = (1/a) sin((w-w')/2)
    bu = mp.diff(lambda p_, q_: (omg(p_, q_) - omgp(p_, q_)) / 2, (U0, V0), (1, 0)) - ag2 * mp.sin((omg(U0, V0) + omgp(U0, V0)) / 2)
    bv = mp.diff(lambda p_, q_: (omg(p_, q_) + omgp(p_, q_)) / 2, (U0, V0), (0, 1)) - mp.sin((omg(U0, V0) - omgp(U0, V0)) / 2) / ag2
    ok &= abs(bu) < 1e-20 and abs(bv) < 1e-20
    E1, F1, G1, L1, M1_, N1, Nv1 = forms(rg, U0, V0)
    E2, F2, G2, L2, M2_, N2, Nv2 = forms(rgp, U0, V0)
    K2 = (L2 * N2 - M2_ ** 2) / (E2 * G2 - F2 ** 2)
    dd = rgp(U0, V0) - rg(U0, V0)
    ok &= abs(K2 + 1) < 1e-15
    ok &= abs(E2 - 1) < 1e-15 and abs(G2 - 1) < 1e-15 and abs(F2 - mp.cos(omgp(U0, V0))) < 1e-15
    ok &= abs(mp.sqrt(dot(dd, dd)) - mp.sin(sig)) < 1e-15
    ok &= abs(dot(dd, Nv1)) < 1e-15 and abs(dot(dd, Nv2)) < 1e-15
    ok &= abs(abs(dot(Nv1, Nv2)) - mp.cos(sig)) < 1e-15
ok &= abs(mp.cos(sig) - mp.tanh(thg)) < 1e-25 and abs(mp.sin(sig) - 1 / mp.cosh(thg)) < 1e-25
report("(j) Baecklund's construction: new surface has K = -1 and angle w' = B_a w; segment of length sin(sigma) tangent "
       "to both; normals at angle sigma; cos(sigma) = tanh(theta), sin(sigma) = 1/cosh(theta)", ok)

# ---------------------------------------------------------------------------------------------
# (k) pushing the second soliton to infinity leaves the first one displaced by the Bianchi coefficient
# ---------------------------------------------------------------------------------------------
# xi2 -> -infinity moves soliton 2 to x = +infinity (Y2 -> 0 in any bounded region): tan(Phi_12/4) -> c Y1.
ok = True
cc = (a1n + a2n) / (a1n - a2n)
for T, X in ((0.3, -1.2), (-0.7, 0.4), (1.1, 2.0)):
    Y1 = math.exp(-0.5 * (a1n + 1 / a1n) * X + 0.5 * (1 / a1n - a1n) * T + 0.3)
    Y2 = math.exp(-0.5 * (a2n + 1 / a2n) * X + 0.5 * (1 / a2n - a2n) * T - 60.0)
    ok &= abs(cc * (Y1 - Y2) / (1 + Y1 * Y2) - cc * Y1) < 1e-12
report("(k) with soliton 2 pushed to x = +infinity, Phi_12 = 4 arctan(c Y1): soliton 1 displaced by c = (a1+a2)/(a1-a2)", ok)

# ---------------------------------------------------------------------------------------------
# (l) Bianchi's cube: three solitons (exercise)
# ---------------------------------------------------------------------------------------------
ths = (1.1, 0.2, -0.7)
aa = [-math.exp(-q) for q in ths]
cs = (0.4, -0.3, 0.9)
def Yn(k, T, X):
    ak = aa[k]
    return math.exp(-0.5 * (ak + 1 / ak) * X + 0.5 * (1 / ak - ak) * T + cs[k])
def tri(T, X):
    P = [4 * math.atan(Yn(k, T, X)) for k in range(3)]
    def pair(i, j):   # Phi_ij from Phi_i, Phi_j over the vacuum
        return 4 * math.atan((aa[i] + aa[j]) / (aa[i] - aa[j]) * math.tan((P[i] - P[j]) / 4))
    P12, P13 = pair(0, 1), pair(0, 2)
    return P[0] + 4 * math.atan((aa[1] + aa[2]) / (aa[1] - aa[2]) * math.tan((P12 - P13) / 4))
hh = 1e-3
res = []
for _ in range(12):
    T, X = rng.uniform(-3, 3), rng.uniform(-3, 3)
    f0 = tri(T, X)
    near = lambda d: f0 + ((tri(*d) - f0 + math.pi) % (2 * math.pi) - math.pi)   # unwrap locally
    ftt = (near((T + hh, X)) - 2 * f0 + near((T - hh, X))) / hh ** 2
    fxx = (near((T, X + hh)) - 2 * f0 + near((T, X - hh))) / hh ** 2
    res.append(abs(ftt - fxx + math.sin(f0)))
report(f"(l) Bianchi's cube gives a three-soliton solution of sine-Gordon (max residual {max(res):.1e})", max(res) < 1e-4)

raise SystemExit(0 if OK else 1)
