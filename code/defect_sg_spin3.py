"""The spin -3 and +3 charges of the sine-Gordon jump defect in closed form (Part VIII, ch. 38: sec-defect-lax).

Written from the book's definitions, without defect_common.py or defect_classical_sg.py.
Conventions (ch. 37): beta = beta_SG, m0, sine-Gordon V = (m0^2/beta^2)(1 - cos beta phi), field equation
phi_tt - phi_xx = -(m0^2/beta) sin beta phi, phi on x < 0, psi on x > 0, d_+- = d_t +- d_x.  Type I_- defect potential
(eq-sg-defect-potential) B = -(2 m0/beta^2)[sigma (cos u - 1) + sigma^-1 (cos v - 1)], u = beta(phi + psi)/2, v = beta(phi - psi)/2,
conditions (eq-sg-defect-conditions) phi_x = psi_t - B_phi, psi_x = phi_t + B_psi at x = 0, sigma = exp(-theta_d).
Spin -3 density q_-(phi) = beta^2 [(d_+^2 phi)^2 - (beta^2/4)(d_+ phi)^4 + m0^2 cos(beta phi)(d_+ phi)^2]
(eq-spin3-current continued to sine-Gordon), spin +3 density q_+ the same with d_-.  Window terms
  D_-3 = (sigma beta^2/(2 m0^2)) cos u d_+phi d_+psi + sigma^3 (cos u - (2/3) cos^3 u),
  D_+3 = -(beta^2/(2 m0^2 sigma)) cos v d_-phi d_-psi + sigma^-3 (cos v - (2/3) cos^3 v),
and Q_-+3 = (1/(8 m0^3)) [int_{x<0} q_-+(phi) + int_{x>0} q_-+(psi)] + D_-+3.

Checks (book ids they support):
 1. Bulk: d_- T + d_+ X = 0 on shell, T = (d_+^2 phi)^2 - (beta^2/4)(d_+ phi)^4, X = m0^2 cos(beta phi)(d_+ phi)^2, so q_- has
    the current j_- = beta^2 (X - T); and the same with d_+ <-> d_- (j_+ = beta^2 (T - X)).  Exact, in jet variables.
 2. Window, spin -3: on the type I_- conditions, (j_-(psi) - j_-(phi))/(8 m0^3) + dD_-3/dt = 0 at x = 0, identically in
    (phi, psi, phi_t, psi_t, phi_tt, psi_tt), sigma symbolic.  Controls: D_-3 with sigma -> 1.1 sigma fails; the spin -3 flux
    with the cos u half of B replaced by cos(0.9 u) is not a total time derivative (nonzero Euler operator).
 3. Window, spin +3: the same for D_+3 and j_+ (exactly), with the same two controls.
 4. Soliton charges: for phi = (4/beta) arctan exp(m0(x cosh theta - t sinh theta)), (1/(8 m0^3)) int q_-+ dx
    = -(2/3) e^{-+3 theta} (quadrature, several theta), and beta^2 (E -+ P)/(8 m0) = e^{-+theta}.
 5. Passage and absorption: along the exact passage (phi, psi) = (soliton, soliton with Y -> zY), z = z_I-(theta)
    (eq-sg-delay), and along the absorption (soliton, 2 pi/beta) at theta = theta_d, Q_-3 and Q_+3 are constant in t and equal
    -(2/3) e^{-+3 theta} + (1/3) sigma^{+-3}; without the window term they are not constant.
 6. Stored charge: D_-+3 at the odd jump (0, 2 pi/beta) minus D_-+3 at the even vacuum (0, 0) is -(2/3) sigma^{+-3}
    = -(2/3) e^{-+3 theta_d}, the charge of a soliton of rapidity theta_d (check 4).  With the constant of D_-+3 as given, D_-+3
    is half the charge of a soliton of rapidity theta_d at the odd jump and of rapidity theta_d + i pi at the vacuum
    (the pattern of prp-half-soliton (1)).
 7. Prior art: Caudrelier's defect term [caudrelier2008, eq. (3.28)], (alpha/3) cos((vt + v)/2)(vt_x^2 + vt_x v_x + v_x^2 + 2 alpha^2)
    with his condition (vt - v)_x = 2 alpha sin((vt + v)/2) (eq. (3.29), x a light-cone coordinate), becomes (2 alpha^3/sigma^3) D_-3
    under vt = beta phi, v = beta psi, d_+ = kappa d_x, alpha = -sigma m0/kappa.
Runtime: about 30 seconds.
"""
import numpy as np
import sympy as sp
from scipy.integrate import quad


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


be, m0, sig = sp.symbols('beta m0 sigma', positive=True)

# ---------- 1. bulk conservation law in jet variables
NJ = 5
U = sp.symbols('u0:%d' % NJ)        # u_n = d_+^n phi
V0 = sp.symbols('v0')               # d_- phi


def Dp(F):
    return sum(U[n + 1]*sp.diff(F, U[n]) for n in range(NJ - 1))


Fsrc = -(m0**2/be)*sp.sin(be*U[0])  # d_+ d_- phi
dm_u = [V0, Fsrc]
for n in range(2, NJ - 1):
    dm_u.append(Dp(dm_u[-1]))


def Dm(F):
    return sum(dm_u[n]*sp.diff(F, U[n]) for n in range(NJ - 2))


T = U[2]**2 - be**2/4*U[1]**4
X = m0**2*sp.cos(be*U[0])*U[1]**2
ok1 = sp.simplify(Dm(T) + Dp(X)) == 0
ctrl1 = sp.simplify(Dm(U[2]**2 + be**2/4*U[1]**4) + Dp(X)) != 0   # the sinh-Gordon sign fails here
report('1. d_-[(d_+^2 phi)^2 - (beta^2/4)(d_+ phi)^4] + d_+[m0^2 cos(beta phi)(d_+ phi)^2] = 0 on shell (and d_+ <-> d_-, '
       'by the same algebra); control: the sinh-Gordon sign fails', ok1 and ctrl1)

# ---------- 2-3. the window terms
p, s, pt, st, ptt, stt = sp.symbols('p s p_t s_t p_tt s_tt')


def Bpot(pp, ss, sg, cu=1):
    u, v = be*(pp + ss)/2, be*(pp - ss)/2
    return -(2*m0/be**2)*(sg*(sp.cos(cu*u) - 1) + (sp.cos(v) - 1)/sg)


def ddt(F):
    return (sp.diff(F, p)*pt + sp.diff(F, s)*st + sp.diff(F, pt)*ptt + sp.diff(F, st)*stt)


def boundary_jets(sg, cu=1):
    B = Bpot(p, s, sg, cu)
    px = st - sp.diff(B, p)
    sx = pt + sp.diff(B, s)
    pxx = ptt + (m0**2/be)*sp.sin(be*p)
    sxx = stt + (m0**2/be)*sp.sin(be*s)
    ptx, stx = ddt(px), ddt(sx)
    return dict(p1p=pt + px, p2p=ptt + 2*ptx + pxx, s1p=st + sx, s2p=stt + 2*stx + sxx,
                p1m=pt - px, p2m=ptt - 2*ptx + pxx, s1m=st - sx, s2m=stt - 2*stx + sxx)


def jminus(f0, f1, f2):
    return be**2*(m0**2*sp.cos(be*f0)*f1**2 - f2**2 + be**2/4*f1**4)


def jplus(f0, f1, f2):
    return be**2*(f2**2 - be**2/4*f1**4 - m0**2*sp.cos(be*f0)*f1**2)


def Dm3(pp, ss, d1p, d1s, sg):
    u = be*(pp + ss)/2
    return sg*be**2/(2*m0**2)*sp.cos(u)*d1p*d1s + sg**3*(sp.cos(u) - sp.Rational(2, 3)*sp.cos(u)**3)


def Dp3(pp, ss, d1p, d1s, sg):
    v = be*(pp - ss)/2
    return -be**2/(2*m0**2*sg)*sp.cos(v)*d1p*d1s + sg**-3*(sp.cos(v) - sp.Rational(2, 3)*sp.cos(v)**3)


def residual(chir, sgD=None, cu=1):
    J = boundary_jets(sig, cu)
    sgD = sig if sgD is None else sgD
    if chir == '-':
        flux = (jminus(s, J['s1p'], J['s2p']) - jminus(p, J['p1p'], J['p2p']))/(8*m0**3)
        DD = Dm3(p, s, J['p1p'], J['s1p'], sgD)
    else:
        flux = (jplus(s, J['s1m'], J['s2m']) - jplus(p, J['p1m'], J['p2m']))/(8*m0**3)
        DD = Dp3(p, s, J['p1m'], J['s1m'], sgD)
    return flux, flux + ddt(DD)


def is_zero_numeric(expr, ntest=6, tol=1e-11):
    f = sp.lambdify((p, s, pt, st, ptt, stt, be, m0, sig), expr, 'numpy')
    rng = np.random.default_rng(1)
    worst = 0.0
    for _ in range(ntest):
        vals = rng.normal(size=6)
        bb, mm, sg = 0.5 + rng.random(), 0.5 + rng.random(), 0.3 + rng.random()
        worst = max(worst, abs(f(*vals, bb, mm, sg)))
    return worst < tol


def euler_P(flux):
    """Largest value at random points of the Euler operator in t, with respect to phi(0, t), of a function of
    (phi, psi) and their first two time derivatives.  A total time derivative has zero Euler operator."""
    p3, p4, s3, s4 = sp.symbols('p_3 p_4 s_3 s_4')
    jp, js = [p, pt, ptt, p3, p4], [s, st, stt, s3, s4]

    def tot(G):
        return sum(jp[k + 1]*sp.diff(G, jp[k]) + js[k + 1]*sp.diff(G, js[k]) for k in range(4))
    E = 0
    for k in range(3):
        term = sp.diff(flux, jp[k])
        for _ in range(k):
            term = -tot(term)
        E += term
    fE = sp.lambdify(jp + js + [be, m0, sig], E, 'numpy')
    rng = np.random.default_rng(2)
    return max(abs(fE(*rng.normal(size=10), 0.8, 1.1, 0.6)) for _ in range(4))


def exact_zero(expr):
    """Exact test: rewrite the trigonometric functions as exponentials and simplify."""
    return sp.simplify(sp.expand(expr.rewrite(sp.exp))) == 0


ok = {}
for chir in '-+':
    flux, res = residual(chir)
    ok_exact = exact_zero(res)
    _, res_bad = residual(chir, sgD=sp.Float(1.1)*sig)
    ctrl_a = not is_zero_numeric(res_bad)
    flux_bad, _ = residual(chir, cu=sp.Float(0.9))
    ok[chir] = ok_exact and ctrl_a and euler_P(flux) < 1e-10 and euler_P(flux_bad) > 1e-3
report('2. spin -3: (j_-(psi) - j_-(phi))/(8 m0^3) + dD_-3/dt = 0 on the type I_- conditions, sigma symbolic; controls: '
       'sigma -> 1.1 sigma in D_-3 fails, cos u -> cos 0.9u in B leaves a nonzero Euler operator', ok['-'])
report('3. spin +3: (j_+(psi) - j_+(phi))/(8 m0^3) + dD_+3/dt = 0, with the same controls', ok['+'])

# ---------- 4-6. numerics on solitons
bv, mv = 0.7, 1.3                     # generic beta, m0
xs, ts = sp.symbols('x t', real=True)
th, cY = sp.symbols('theta cY', real=True)
Ysym = sp.exp(mv*(xs*sp.cosh(th) - ts*sp.sinh(th)))
fsol = 4/bv*sp.atan(cY*Ysym)


def dens(fexpr, chir):
    d = (lambda F: sp.diff(F, ts) + sp.diff(F, xs)) if chir == '-' else (lambda F: sp.diff(F, ts) - sp.diff(F, xs))
    d1, d2 = d(fexpr), d(d(fexpr))
    q = bv**2*(d2**2 - bv**2/4*d1**4 + mv**2*sp.cos(bv*fexpr)*d1**2)
    return sp.lambdify((xs, ts, th, cY), q, 'numpy'), sp.lambdify((xs, ts, th, cY), d1, 'numpy'), \
        sp.lambdify((xs, ts, th, cY), fexpr, 'numpy')


qm, d1m, fval = dens(fsol, '-')
qp, d1p_, _ = dens(fsol, '+')
ok4 = True
for thv in [0.0, 0.4, -0.9, 1.3]:
    Lx = 40/(mv*np.cosh(thv))
    Qm = quad(lambda x: qm(x, 0.0, thv, 1.0), -Lx, Lx, limit=400, points=[0])[0]/(8*mv**3)
    Qp = quad(lambda x: qp(x, 0.0, thv, 1.0), -Lx, Lx, limit=400, points=[0])[0]/(8*mv**3)
    ok4 &= abs(Qm + 2/3*np.exp(-3*thv)) < 1e-8 and abs(Qp + 2/3*np.exp(3*thv)) < 1e-8
    # energy and momentum: E -+ P = int (1/2)(d_-+ phi)^2 + V
    Vf = lambda x: mv**2/bv**2*(1 - np.cos(bv*fval(x, 0.0, thv, 1.0)))
    Em = quad(lambda x: 0.5*d1m(x, 0.0, thv, 1.0)**2 + Vf(x), -Lx, Lx, limit=400)[0]
    Ep = quad(lambda x: 0.5*d1p_(x, 0.0, thv, 1.0)**2 + Vf(x), -Lx, Lx, limit=400)[0]
    ok4 &= abs(bv**2*Em/(8*mv) - np.exp(-thv)) < 1e-8 and abs(bv**2*Ep/(8*mv) - np.exp(thv)) < 1e-8
report('4. a soliton of rapidity theta has (1/(8 m0^3)) int q_-+ = -(2/3) e^{-+3 theta} and beta^2 (E -+ P)/(8 m0) = e^{-+theta}',
       ok4)


def check_history(thv, sgv, z_right, psi_const=None):
    """Q_-3, Q_+3 (with and without window terms) along a solution; returns spreads and means."""
    out = {}
    Lx = 40/(mv*np.cosh(thv))
    for chir, qf, d1f, Dfun in [('-', qm, d1m, Dm3), ('+', qp, d1p_, Dp3)]:
        Dl = sp.lambdify((p, s, sp.Symbol('a'), sp.Symbol('c'), be, m0, sig),
                         Dfun(p, s, sp.Symbol('a'), sp.Symbol('c'), sig), 'numpy')
        tot, bulk = [], []
        for tv in np.linspace(-6, 6, 9):
            Ql = quad(lambda x: qf(x, tv, thv, 1.0), -Lx + tv*np.tanh(thv), 0, limit=400)[0]
            if psi_const is None:
                Qr = quad(lambda x: qf(x, tv, thv, z_right), 0, Lx + tv*np.tanh(thv), limit=400)[0]
                sv, d1s = fval(0.0, tv, thv, z_right), d1f(0.0, tv, thv, z_right)
            else:
                Qr, sv, d1s = 0.0, psi_const, 0.0
            pv, d1pv = fval(0.0, tv, thv, 1.0), d1f(0.0, tv, thv, 1.0)
            b3 = (Ql + Qr)/(8*mv**3)
            bulk.append(b3)
            tot.append(b3 + Dl(pv, sv, d1pv, d1s, bv, mv, sgv))
        out[chir] = (np.ptp(tot), np.mean(tot), np.ptp(bulk))
    return out


f_x = sp.lambdify((xs, ts, th, cY), sp.diff(fsol, xs), 'numpy')
f_t = sp.lambdify((xs, ts, th, cY), sp.diff(fsol, ts), 'numpy')
Bphi = sp.lambdify((p, s, be, m0, sig), sp.diff(Bpot(p, s, sig), p), 'numpy')
Bpsi = sp.lambdify((p, s, be, m0, sig), sp.diff(Bpot(p, s, sig), s), 'numpy')


def cond_resid(thv, sgv, z_right, psi_const=None):
    """Residual of the type I_- conditions phi_x = psi_t - B_phi, psi_x = phi_t + B_psi at x = 0."""
    worst = 0.0
    for tv in np.linspace(-6, 6, 13):
        pv, px_, pt_ = fval(0.0, tv, thv, 1.0), f_x(0.0, tv, thv, 1.0), f_t(0.0, tv, thv, 1.0)
        if psi_const is None:
            sv, sx_, st_ = fval(0.0, tv, thv, z_right), f_x(0.0, tv, thv, z_right), f_t(0.0, tv, thv, z_right)
        else:
            sv, sx_, st_ = psi_const, 0.0, 0.0
        worst = max(worst, abs(px_ - (st_ - Bphi(pv, sv, bv, mv, sgv))), abs(sx_ - (pt_ + Bpsi(pv, sv, bv, mv, sgv))))
    return worst


ok5 = True
# passage: theta != theta_d, right field the soliton with Y -> zY
for thv, thd in [(0.5, 1.2), (1.1, 0.4)]:
    sgv = np.exp(-thd)
    zv = (np.exp(-thv) + sgv)/(np.exp(-thv) - sgv)
    ok5 &= cond_resid(thv, sgv, zv) < 1e-12 and cond_resid(thv, 1.1*sgv, zv) > 1e-3
    out = check_history(thv, sgv, zv)
    ok5 &= (out['-'][0] < 1e-8 and abs(out['-'][1] - (-2/3*np.exp(-3*thv) + sgv**3/3)) < 1e-8 and out['-'][2] > 1e-3
            and out['+'][0] < 1e-7 and abs(out['+'][1] - (-2/3*np.exp(3*thv) + sgv**-3/3)) < 1e-7 and out['+'][2] > 1e-3)
# absorption at theta = theta_d > 0: (soliton, 2 pi/beta)
thd = 0.6
sgv = np.exp(-thd)
ok5 &= cond_resid(thd, sgv, None, psi_const=2*np.pi/bv) < 1e-12
out = check_history(thd, sgv, None, psi_const=2*np.pi/bv)
ok5 &= (out['-'][0] < 1e-8 and abs(out['-'][1] + sgv**3/3) < 1e-8 and out['-'][2] > 1e-3
        and out['+'][0] < 1e-7 and abs(out['+'][1] + sgv**-3/3) < 1e-7 and out['+'][2] > 1e-3)
report('5. Q_-3 and Q_+3 are constant along exact passages (theta < theta_d and theta > theta_d) and along the absorption '
       'at theta = theta_d, with the predicted values; the bulk parts alone are not constant', ok5)

# 6. stored charge
a_, c_ = sp.symbols('a c')
odd = [Dm3(0, 2*sp.pi/be, 0, 0, sig) - Dm3(0, 0, 0, 0, sig), Dp3(0, 2*sp.pi/be, 0, 0, sig) - Dp3(0, 0, 0, 0, sig)]
ok6 = (sp.simplify(odd[0] + sp.Rational(2, 3)*sig**3) == 0 and sp.simplify(odd[1] + sp.Rational(2, 3)*sig**-3) == 0
       and sp.simplify(Dm3(0, 0, 0, 0, sig) - sig**3/3) == 0)
# half-soliton form: D at the odd jump is half the charge -(2/3) e^{-+3 theta} at theta = theta_d, and at the vacuum half the
# charge at theta_d + i pi (sigma = e^{-theta_d})
thd_ = sp.symbols('theta_d')
half = lambda thv, sgn: -sp.Rational(1, 3)*sp.exp(sgn*3*thv)
sub = {sig: sp.exp(-thd_)}
ok6 &= all(sp.simplify(x) == 0 for x in [
    Dm3(0, 2*sp.pi/be, 0, 0, sig).subs(sub) - half(thd_, -1), Dp3(0, 2*sp.pi/be, 0, 0, sig).subs(sub) - half(thd_, 1),
    Dm3(0, 0, 0, 0, sig).subs(sub) - half(thd_ + sp.I*sp.pi, -1), Dp3(0, 0, 0, 0, sig).subs(sub) - half(thd_ + sp.I*sp.pi, 1)])
report('6. the odd jump stores D_-+3(0, 2pi/beta) - D_-+3(0, 0) = -(2/3) sigma^{+-3} = -(2/3) e^{-+3 theta_d}, the spin -+3 '
       'charge of a soliton of rapidity theta_d; D_-+3 is half the charge of a soliton of rapidity theta_d at the odd jump and '
       'of rapidity theta_d + i pi at the vacuum', ok6)

# 7. Caudrelier's (3.28) after his defect condition (3.29)
kap = sp.symbols('kappa', nonzero=True)
P1, S1 = sp.symbols('P1 S1')          # d_+ phi, d_+ psi
alpha = -sig*m0/kap
uu = be*(p + s)/2
vt_x, v_x = be*P1/kap, be*S1/kap      # vt = beta phi, v = beta psi, d_x(Caudrelier) = d_+/kappa
DC = alpha/3*sp.cos(uu)*(vt_x**2 + vt_x*v_x + v_x**2 + 2*alpha**2)
# his (3.29): (vt - v)_x = 2 alpha sin((vt + v)/2), i.e. d_+(phi - psi) = -2 sigma m0 sin(u)/beta, the book's condition
cond_ok = sp.simplify((vt_x - v_x).subs(S1, P1 + 2*sig*m0*sp.sin(uu)/be) - 2*alpha*sp.sin(uu)) == 0
DC_on = DC.subs(S1, P1 + 2*sig*m0*sp.sin(uu)/be)
ours = (2*alpha**3/sig**3)*Dm3(p, s, P1, P1 + 2*sig*m0*sp.sin(uu)/be, sig)
ok7 = cond_ok and sp.simplify(sp.expand(sp.expand_trig(DC_on - ours))) == 0
report('7. Caudrelier (3.28), on his condition (3.29), is (2 alpha^3/sigma^3) D_-3: the closed form is his, in light-cone '
       'variables', ok7)
