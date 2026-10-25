"""A sine-Gordon defect fused onto a Dirichlet boundary, the defect UV-IR heuristic, folding, and chirality.

Book ids: prp-defect-fusion-k (3) (ch. 43, sec-defect-boundary-fusion); prp-sg-fused-spin-three,
cnj-defect-parameter (ch. 38); sec-defect-folding (ch. 39); prp-diagonal-borel, "Reading the result" (ch. 43).

Written without defect_common.py, from the book's definitions: sigma(x, u) and R_0 of (eq-bsg-k-full) (ch. 30),
Phi of (eq-sgd-scalar) and (eq-sgd-scalar-integral) (ch. 38), the transmission matrix (eq-sg-transmission-matrix),
the spin-3 densities and defect terms of prp-sg-spin-three (eq-dc-spin-three, ch. 39), and the linear factors
(eq-an-linear-transmission, ch. 40).  Conventions: y = exp(lambda theta), y_d = exp(lambda theta_d), u = -i theta in
the reflection matrices, Y = y/y_d, beta = beta_SG, sigma = exp(-theta_d^cl), partial_pm = partial_t pm partial_x.

Checks:
 1. Malmstén step for sigma: the four Gamma functions of the l-th factor give the signed sum
    -4 e^{-t/2} e^{-(2l-1) lam t} sinh(lam t/2) sinh(xt/pi + lam t/2); the sum over l gives the kernel of
    ln sigma(x,u) = ln(cos x/cos(x+lam u)) - int sinh(lam u t/pi) sinh(xt/pi+lam t/2)/(sinh(t/2) cosh(lam t/2)) dt/t
    (sympy).  The integral equals the Gamma product at complex (x, u) in the strip |Re x| + |Re lam u| < pi/2, and
    sigma(i vartheta, u) -> 1 as vartheta -> infinity (the Dirichlet limit used for K_D), also in the split form of the
    proof, ln cosh(lam d) - (1/2) ln cosh(lam(theta-d)) - (1/2) ln cosh(lam(theta+d)) + i int (bounded amplitude) cos(kd),
    against the Gamma product; the raw sin(k theta) sin(kd) part tends to lam theta, so it is not a Riemann-Lebesgue term.
 2. The two integrals of the proof at x = i lam theta_d, u = -i theta (k = lam t/pi), against the Gamma products of
    sigma and of Phi; the integral of (1 - cos kv)/(k sinh(pi k/2lam)) is ln cosh(lam v); the elementary remainder
    lam theta + ln cosh(lam theta_d) - ln cosh(lam(theta+theta_d)) is ln((1+y_d^-2)/(1+y^-2 y_d^-2)) (sympy).
 3. The scalar identity sigma(i lam theta_d, -i theta) = (1+y_d^-2)/(1+y^-2 y_d^-2) Phi(theta-theta_d)/Phi(-theta-theta_d)
    with both sides from their Gamma products, three couplings, real and complex rapidities.  Control: with
    (lam+1) theta_d in sigma it fails.
 4. The matrix part (sympy, identically in y, y_d, e^{i eta}, q^n, q^s): on the sector {v+ |n>, v- |n+1>},
    T(-theta)^-1 K_D T(theta) = c (Lambda D_n)^-1 M_GZ (Lambda D_n), Lambda = diag(1, -e^{i eta}),
    D_n = diag(q^-n, q^{n+1+s}), at vartheta = lam theta_d; c = 1 is exactly the scalar identity of check 3.
    Controls: vartheta = -lam theta_d and vartheta = theta_d fail.
 5. Spin-3 charges (sympy, exact): q_{-3}, q_{+3} are conserved with the currents j_{-3}, j_{+3} on shell; the defect
    flux equals -d/dt (D_{-3} + D_{+3}) on the type-I_- conditions (re-derives the conventions of prp-sg-spin-three);
    at a Dirichlet end the flux of q_{-3}+q_{+3} is -(beta/m0) sin(beta phi_D) d/dt psi_x.
 6. The fused boundary (prp-sg-fused-spin-three): with phi_x = -dB_d/dphi(phi, phi_D) at x = 0 and
    psi_x = phi_t + dB_d/dpsi(phi, phi_D), Q = int (q_{-3}+q_{+3}) + D_{-3} + D_{+3} - (beta/m0) sin(beta phi_D) psi_x
    is conserved identically, sigma and phi_D symbolic.  The steps of the proof by hand: with f = phi(t,0), v = f_t,
    a = f_tt, F = phi_x(f), the flux is -(2b^2/m0^3) F' v a + (b^4/4m0^3) F v^3 + U' v; F'' = -b^2 F/4; the boundary
    terms are -(b^2/m0^3) F' v^2 + U(f), the defect terms alone having the linear term (b/m0) sin(b phi_D) v.
    Controls: without the last term, and with q_{-3} - q_{+3} (and D_{-3} - D_{+3}), it is not.
 7. The boost heuristic for cnj-defect-parameter: h = beta^2/8pi = 1/(lam+1); theta_d^cl = (1-h) theta_d is
    lam theta_d = (lam+1) theta_d^cl; at lam = -2/B it is theta_d = (1 - B/2) theta_d^cl, and with
    B = 2b^2/(8pi+b^2) its first correction is -b^2/8pi (Bajnok-Simon (17), (23)).  Bajnok-Holló-Watts (101)-(102):
    -i(b+3)pi/5 = (6/5)(theta_b - i pi) with theta_b = i pi(3-b)/6, and 6/5 = 1 - h for h = -1/5.  Konik-LeClair
    Sec. 5.1: d(g) = 1 - betahat^2/2 with betahat^2 = beta^2/4pi is 1 - h.
 8. Folding (sec-defect-folding): rho pi rho = pi^-1 and rho maps the omega^k eigenvectors of pi to omega^-k ones, and
    D^T on mode k equals D on mode h-k (h = 3..6), so the mirror image transmits particle k with T_{h-k}.
    T_k/T_{h-k} = (1 - x^2 e^{2i gamma})/(1 - x^2 e^{-2i gamma}); the pair folds iff sigma_2 = +-sigma_1
    (controls: four other ratios fail); -sigma_1 is transparent; +sigma_1 gives 1/f_{k/h}(u) and the folded soliton
    delay (cosh u + sin gamma)/(cosh u - sin gamma); for h = 3 this is the Tzitzeica z_1 of ch. 39; for 2k = h the
    ratio is identically 1 (no condition) and the delay is coth^2(u/2).
 9. Chirality (prp-diagonal-borel, reading): at m0 = 0 the type-I_- conditions give psi_L = phi_L, psi_R = -phi_R and
    type I_+ gives psi_L = -phi_L, psi_R = phi_R; a left-mover f(t+x) carries E - P = (1/2)(d_+ phi)^2 and no E + P;
    on V(theta) the f_1 act with 1/y and the e_1 with y; the type-I_- matrix part tends to the twist
    diag(q^-n, q^{n+s}) as theta - theta_d -> -infinity and to the off-diagonal conversion as theta - theta_d -> +infinity.
Runtime: about two and a half minutes.
"""
import mpmath as mp
import numpy as np
import sympy as sp

mp.mp.dps = 25
PI = mp.pi
LAMS = [mp.sqrt(2), mp.e, 2 * mp.sqrt(3)]   # generic couplings: q is not a root of unity


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), mp.mpf(1))


# ---------------------------------------------------------------- scalar factors from the book's Gamma products

def lnsigma_gamma(x, u, lam):
    """ln sigma(x, u) of (eq-bsg-k-full): cos x/cos(x + lam u) times prod_l F_l(u)/F_l(-u)."""
    w = lam * u / PI

    def term(l):
        a1 = mp.mpf(1) / 2 + x / PI + (2 * l - 1) * lam
        a2 = mp.mpf(1) / 2 - x / PI + (2 * l - 1) * lam
        a3 = mp.mpf(1) / 2 - x / PI + (2 * l - 2) * lam
        a4 = mp.mpf(1) / 2 + x / PI + 2 * l * lam
        d = lambda a: mp.loggamma(a - w) - mp.loggamma(a + w)
        return d(a1) + d(a2) - d(a3) - d(a4)
    return mp.log(mp.cos(x) / mp.cos(x + lam * u)) + mp.nsum(term, [1, mp.inf])


def lnphi_gamma(th, lam):
    """ln Phi(theta) of (eq-sgd-scalar): prefactor times lim r_N/rbar_N, ratio taken factor by factor."""
    z = 1j * lam * th / (2 * PI)

    def term(k):
        args = [(k * lam + mp.mpf(1) / 4, 1), ((k + 1) * lam + mp.mpf(3) / 4, 1),
                ((k + mp.mpf(1) / 2) * lam + mp.mpf(1) / 4, -1), ((k + mp.mpf(1) / 2) * lam + mp.mpf(3) / 4, -1)]
        return sum(s * (mp.loggamma(a - z) - mp.loggamma(a + z)) for a, s in args)
    return 1j * PI * (1 + lam) / 4 - mp.log(1 + 1j * mp.exp(lam * th)) + mp.nsum(term, [0, mp.inf])


def lnsigma_int(x, u, lam):
    """The integral representation derived in the proof of prp-defect-fusion-k (3)."""
    def f(t):
        if t == 0:
            return lam * u / PI * 2 * (x / PI + lam / 2)
        return (mp.sinh(lam * u * t / PI) * mp.sinh(x * t / PI + lam * t / 2)
                / (mp.sinh(t / 2) * mp.cosh(lam * t / 2)) / t)
    return mp.log(mp.cos(x) / mp.cos(x + lam * u)) - mp.quad(f, [0, 1, 10, mp.inf])


# ---------------------------------------------------------------- 1. the Malmstén step and the integral of ln sigma

def check_sigma_integral():
    t, l, lam, X = sp.symbols('t l lambda X', positive=True)
    e = sp.exp
    signed = (e(-t / 2) * (e(-X) * e(-(2 * l - 1) * lam * t) + e(X) * e(-(2 * l - 1) * lam * t)
                           - e(X) * e(-(2 * l - 2) * lam * t) - e(-X) * e(-2 * l * lam * t)))
    closed = -4 * e(-t / 2) * e(-(2 * l - 1) * lam * t) * sp.sinh(lam * t / 2) * sp.sinh(X + lam * t / 2)
    ok1 = sp.simplify((signed - closed).rewrite(sp.exp)) == 0
    geo = e(-lam * t) / (1 - e(-2 * lam * t))                       # sum_{l>=1} e^{-(2l-1) lam t}, geometric
    ok2 = sp.simplify((geo - 1 / (2 * sp.sinh(lam * t))).rewrite(sp.exp)) == 0
    # kernel: 2 sinh(wt)/(1-e^-t) * (-e^{-t/2} sinh(X + lam t/2)/cosh(lam t/2)) = -sinh(wt) sinh(..)/(sinh(t/2) cosh(..))
    w = sp.symbols('w')
    lhs = 2 * sp.sinh(w * t) / (1 - e(-t)) * (-4 * e(-t / 2) * sp.sinh(lam * t / 2) * sp.sinh(X + lam * t / 2)
                                              / (2 * sp.sinh(lam * t)))
    rhs = -sp.sinh(w * t) * sp.sinh(X + lam * t / 2) / (sp.sinh(t / 2) * sp.cosh(lam * t / 2))
    ok3 = sp.simplify((lhs - rhs).rewrite(sp.exp)) == 0
    errs = []
    for lam_n in LAMS[:2]:
        for x, u in [(mp.mpf('0.4') + 0.7j, mp.mpf('0.1') - 0.3j), (-0.2 + 0.3j, 0.05 + 0.9j), (0.9j, -0.6j)]:
            errs.append(rel(mp.exp(lnsigma_int(x, u, lam_n)), mp.exp(lnsigma_gamma(x, u, lam_n))))
    dev = [abs(mp.exp(lnsigma_gamma(1j * vt, 0.3j, LAMS[0])) - 1) for vt in (10, 20)]   # Dirichlet limit
    # the Dirichlet limit as proved: at x = i lam d, u = -i theta, real theta, ln sigma = ln cosh(lam d)
    # - (1/2) ln cosh(lam(theta-d)) - (1/2) ln cosh(lam(theta+d)) + i int sin(k theta) tanh(pi k/2) cos(kd)/(k sinh(pi k/2lam)),
    # whose amplitude is bounded at k = 0; the raw sin(k theta) sin(kd)/(k sinh) part tends to lam theta, not to 0
    lam, th = LAMS[0], mp.mpf('0.43')
    kq = lambda f: mp.quad(f, mp.linspace(0, 60, 121) + [mp.inf])
    rl, raw = [], []
    for d in (1, 2, 4, 8):
        osc = kq(lambda k: mp.sin(k * th) * mp.tanh(PI * k / 2) * mp.cos(k * d) / (k * mp.sinh(PI * k / (2 * lam))))
        r1 = (mp.log(mp.cosh(lam * d)) - (mp.log(mp.cosh(lam * (th - d))) + mp.log(mp.cosh(lam * (th + d)))) / 2
              + 1j * osc)
        errs.append(rel(mp.exp(r1), mp.exp(lnsigma_gamma(1j * lam * d, -1j * th, lam))))
        rl.append(abs(r1))
        raw.append(kq(lambda k: mp.sin(k * th) * mp.sin(k * d) / (k * mp.sinh(PI * k / (2 * lam)))))
    ok4 = rl[-1] < rl[0] / 100 and abs(raw[-1] - lam * th) < 1e-6
    report('Malmsten step for sigma (sympy) and the integral of ln sigma = Gamma product, max rel err %.1e; '
           'sigma(i vartheta, u) - 1 = %.0e, %.0e at vartheta = 10, 20; split form of the Dirichlet limit: |ln sigma| '
           '= %.0e at d = 8, raw part -> lam theta' % (float(max(errs)), dev[0], dev[1], rl[-1]),
           ok1 and ok2 and ok3 and ok4 and max(errs) < 1e-12 and dev[1] < dev[0] / 100)


# ---------------------------------------------------------------- 2. the steps of the proof of the scalar identity

def check_proof_steps():
    errs = []
    for lam in LAMS:
        for th, thd in [(mp.mpf('0.37'), mp.mpf('-0.52')), (mp.mpf('-1.1'), mp.mpf('0.23')), (mp.mpf('0.8'), mp.mpf('0.6'))]:
            kint = lambda f: mp.quad(f, [0, 1, 10, mp.inf])   # tanh-sinh: the endpoint k = 0 is never sampled
            kern = lambda k, s: (mp.sin(k * th) * mp.sinh(PI * k / 2 + s * 1j * k * thd)
                                 / (mp.sinh(PI * k / (2 * lam)) * mp.cosh(PI * k / 2) * k))
            lnsig = mp.log(mp.cosh(lam * thd) / mp.cosh(lam * (th - thd))) + 1j * kint(lambda k: kern(k, 1))
            errs.append(rel(mp.exp(lnsig), mp.exp(lnsigma_gamma(1j * lam * thd, -1j * th, lam))))
            lnphir = -lam * th + 1j * kint(lambda k: kern(k, -1))
            errs.append(rel(mp.exp(lnphir), mp.exp(lnphi_gamma(th - thd, lam) - lnphi_gamma(-th - thd, lam))))
            for v in (th - thd, th + thd):
                errs.append(abs(kint(lambda k: (1 - mp.cos(k * v)) / (k * mp.sinh(PI * k / (2 * lam))))
                                - mp.log(mp.cosh(lam * v))))
    lam, th, thd = sp.symbols('lambda theta theta_d', real=True)
    y, yd = sp.exp(lam * th), sp.exp(lam * thd)
    rem = lam * th + sp.log(sp.cosh(lam * thd)) - sp.log(sp.cosh(lam * (th + thd)))
    ok = sp.simplify(sp.exp(rem).rewrite(sp.exp) - (1 + yd**-2) / (1 + y**-2 * yd**-2)) == 0
    report('proof steps: the two integrals, the ln cosh integral, the elementary remainder, max err %.1e'
           % float(max(errs)), ok and max(errs) < 1e-12)


# ---------------------------------------------------------------- 3. the scalar identity from the Gamma products

def check_scalar_identity():
    errs, ctrl = [], []
    for lam in LAMS:
        for th, thd in [(0.37, -0.52), (-1.1, 0.23), (0.8, 0.6), (0.3 + 0.05j, -0.2), (0.5, 0.4 + 0.3j)]:
            th, thd = mp.mpc(th), mp.mpc(thd)
            y, yd = mp.exp(lam * th), mp.exp(lam * thd)
            rhs = (1 + yd**-2) / (1 + y**-2 * yd**-2) * mp.exp(lnphi_gamma(th - thd, lam) - lnphi_gamma(-th - thd, lam))
            errs.append(rel(mp.exp(lnsigma_gamma(1j * lam * thd, -1j * th, lam)), rhs))
            ctrl.append(rel(mp.exp(lnsigma_gamma(1j * (lam + 1) * thd, -1j * th, lam)), rhs))
    report('scalar identity of prp-defect-fusion-k (3) from the Gamma products, max rel err %.1e; '
           'control (lam+1) theta_d min err %.1e' % (float(max(errs)), float(min(ctrl))),
           max(errs) < 1e-12 and min(ctrl) > 1e-3)


# ---------------------------------------------------------------- 4. the matrix part

def check_matrix_part():
    y, yd, E, Qn, Qs, PhiR, Sig = sp.symbols('y y_d E Q_n Q_s PhiR Sigma', nonzero=True)
    A = lambda Y: sp.Matrix([[1, -sp.I * Y], [-sp.I * Y, 1]])
    Dn = sp.diag(1 / Qn, Qn * Qs * sp.Symbol('q'))           # diag(q^-n, q^{n+1+s}); q^{n+1+s} = Q_n q Q_s
    Y, Yp = y / yd, 1 / (y * yd)
    kp, km = y * E + 1 / (y * E), y / E + E / y
    # K_f / (R_0 sigma(eta,u)) with PhiR = Phi(theta-theta_d)/Phi(-theta-theta_d); cos eta = (E + 1/E)/2
    cos_eta = (E + 1 / E) / 2
    Kf = PhiR * (Dn.inv() * A(Yp).inv() * sp.diag(kp, km) * A(Y) * Dn) / (2 * cos_eta)
    res = {}
    for name, vt in [('lam theta_d', yd), ('-lam theta_d', 1 / yd), ('theta_d', sp.Symbol('z_d'))]:
        ch, sh = (vt + 1 / vt) / 2, (vt - 1 / vt) / 2
        cl, sl = (y + 1 / y) / 2, (y - 1 / y) / (2 * sp.I)
        ce, se = cos_eta, (E - 1 / E) / (2 * sp.I)
        M = sp.Matrix([[cl * ce * ch - sl * se * sh, -sl * cl], [-sl * cl, cl * ce * ch + sl * se * sh]])
        G = Sig / (ce * ch) * M                                 # (eq-bsg-k-full) / (R_0 sigma(eta,u))
        Lam = sp.diag(1, -E)
        target = (Lam * Dn).inv() * G * (Lam * Dn)
        ratios = [sp.simplify(Kf[i, j] / target[i, j]) for i in range(2) for j in range(2)]
        res[name] = (all(sp.simplify(r - ratios[0]) == 0 for r in ratios), ratios[0])
    ok_main, ratio = res['lam theta_d']
    expected = PhiR * (1 + yd**-2) / ((1 + y**-2 * yd**-2) * Sig)
    ok_scalar = sp.simplify(ratio - expected) == 0
    report('matrix part: K_f = (Lambda D_n)^-1 GZ (Lambda D_n) at vartheta = lam theta_d, common factor 1 iff the '
           'scalar identity; controls -lam theta_d, theta_d fail',
           ok_main and ok_scalar and not res['-lam theta_d'][0] and not res['theta_d'][0])


# ---------------------------------------------------------------- 5, 6. spin-3 charges

b, m0, sg, phD = sp.symbols('beta m_0 sigma phi_D', positive=True)


def q3(d1, d2, f, sign):
    """q_{-3} (sign=-1, d = partial_+) or q_{+3} (sign=+1, d = partial_-), prp-sg-spin-three."""
    return b**2 / (8 * m0**3) * (d2**2 - b**2 / 4 * d1**4 + m0**2 * sp.cos(b * f) * d1**2)


def j3(d1, d2, f, sign):
    """j_{-3} of the proof of prp-sg-spin-three; j_{+3} by parity (partial_+ <-> partial_-, j -> -j)."""
    jm = b**2 / (8 * m0**3) * (m0**2 * sp.cos(b * f) * d1**2 - d2**2 + b**2 / 4 * d1**4)
    return jm if sign < 0 else -jm


def check_bulk_currents():
    xp, xm = sp.symbols('x_p x_m')
    F = sp.Function('F')(xp, xm)
    eom = -m0**2 / (4 * b) * sp.sin(b * F)                         # F_{x_p x_m}: d_+ d_- = 4 d_{xp} d_{xm}

    def reduce(expr):
        for _ in range(8):
            ders = [d for d in expr.atoms(sp.Derivative)
                    if d.expr == F and dict(d.variable_count).get(xp, 0) and dict(d.variable_count).get(xm, 0)]
            if not ders:
                return expr
            sub = {}
            for d in ders:
                vc = dict(d.variable_count)
                g = eom
                if vc[xp] > 1:
                    g = sp.diff(g, xp, vc[xp] - 1)
                if vc[xm] > 1:
                    g = sp.diff(g, xm, vc[xm] - 1)
                sub[d] = g
            expr = expr.subs(sub)
        return expr
    dp = lambda g: 2 * sp.diff(g, xp)
    dm = lambda g: 2 * sp.diff(g, xm)
    dt = lambda g: sp.diff(g, xp) + sp.diff(g, xm)
    dx = lambda g: sp.diff(g, xp) - sp.diff(g, xm)
    ok = True
    for sign, d in [(-1, dp), (1, dm)]:
        q, j = q3(d(F), d(d(F)), F, sign), j3(d(F), d(d(F)), F, sign)
        ok &= sp.simplify(reduce(sp.expand(dt(q) + dx(j)))) == 0
    return ok


def is_zero(expr, angles):
    """Exact zero test for expressions in cos, sin of the fields: with beta a/2 -> A for each field a in angles, every
    exponential e^{i(nA + mB + c)} becomes U_A^n U_B^m e^{ic}, and the numerator of the rational function must vanish."""
    A = sp.symbols('A0:%d' % len(angles))
    Us = sp.symbols('U0:%d' % len(angles), nonzero=True)
    e = expr.subs({a: 2 * Ai / b for a, Ai in zip(angles, A)}).rewrite(sp.exp)

    def to_u(z):
        arg = sp.expand(z.args[0] / sp.I)
        p = sp.Poly(arg, *A)
        cs = [p.coeff_monomial(Ai) for Ai in A]
        out = sp.exp(sp.I * sp.expand(arg - sum(c * Ai for c, Ai in zip(cs, A))))
        for c, Ui in zip(cs, Us):
            out *= Ui**c
        return out
    e = e.replace(lambda z: isinstance(z, sp.exp) and z.args[0].has(*A), to_u)
    return sp.expand(sp.numer(sp.together(sp.expand(e)))) == 0


def boundary_jets():
    """phi(0,t) = f with time derivatives f1, f2, f3; D_t is the total time derivative on these jets."""
    f, f1, f2, f3, g, g1, g2, g3 = sp.symbols('f f1 f2 f3 g g1 g2 g3')
    Dt = lambda e, jets: sum(sp.diff(e, a) * a_next for a, a_next in jets)
    return (f, f1, f2, f3), (g, g1, g2, g3), Dt


def Bd(ph, ps):
    """(eq-sg-defect-potential), type I_-."""
    return -2 * m0 / b**2 * (sg * (sp.cos(b * (ph + ps) / 2) - 1) + (sp.cos(b * (ph - ps) / 2) - 1) / sg)


def Dterms(ph, ps, dpp, dpm, dsp, dsm):
    """D_{-3}, D_{+3} of (eq-dc-spin-three) with c_pm = cos(beta(phi pm psi)/2)."""
    cp, cm = sp.cos(b * (ph + ps) / 2), sp.cos(b * (ph - ps) / 2)
    Dm = sg * b**2 / (2 * m0**2) * cp * dpp * dsp + sg**3 * (cp - sp.Rational(2, 3) * cp**3)
    Dp = -b**2 / (2 * m0**2 * sg) * cm * dpm * dsm + (cm - sp.Rational(2, 3) * cm**3) / sg**3
    return Dm, Dp


def second_derivs(f, f1, f2, fx, fx_t):
    """partial_pm^2 phi at x = 0 from phi_tt = f2, phi_tx = d/dt phi_x, phi_xx = phi_tt + (m0^2/beta) sin(beta phi)."""
    fxx = f2 + m0**2 / b * sp.sin(b * f)
    return f2 + 2 * fx_t + fxx, f2 - 2 * fx_t + fxx


def check_defect_flux():
    (f, f1, f2, f3), (g, g1, g2, g3), Dt = boundary_jets()
    jets = [(f, f1), (f1, f2), (f2, f3), (g, g1), (g1, g2), (g2, g3)]
    PH, PS = sp.symbols('PH PS')
    Bphi = sp.diff(Bd(PH, PS), PH).subs({PH: f, PS: g})
    Bpsi = sp.diff(Bd(PH, PS), PS).subs({PH: f, PS: g})
    phx, psx = g1 - Bphi, f1 + Bpsi                                # (eq-sg-defect-conditions), D = 1
    php, phm, psp, psm = f1 + phx, f1 - phx, g1 + psx, g1 - psx
    p2p, p2m = second_derivs(f, f1, f2, phx, Dt(phx, jets))
    s2p, s2m = second_derivs(g, g1, g2, psx, Dt(psx, jets))
    Dm, Dp = Dterms(f, g, php, phm, psp, psm)
    ok = True
    for sign, (a1, a2, c1, c2), D in [(-1, (php, p2p, psp, s2p), Dm), (1, (phm, p2m, psm, s2m), Dp)]:
        flux = -j3(a1, a2, f, sign) + j3(c1, c2, g, sign)
        ok &= is_zero(flux + Dt(D, jets), [f, g])
        ok &= not is_zero(flux + Dt(D.subs(sg, sg * sp.Rational(11, 10)), jets), [f, g])   # control
    # Dirichlet end: psi_t = 0, psi = phi_D; flux of q_-3 + q_+3 is -(beta/m0) sin(beta phi_D) psi_tx
    psx_, pstx = sp.symbols('psi_x psi_tx')
    s2p_d, s2m_d = 2 * pstx + m0**2 / b * sp.sin(b * phD), -2 * pstx + m0**2 / b * sp.sin(b * phD)
    jD = j3(psx_, s2p_d, phD, -1) + j3(-psx_, s2m_d, phD, 1)
    okD = sp.simplify(jD + b / m0 * sp.sin(b * phD) * pstx) == 0
    report('spin-3 currents conserved on shell; defect flux = -d/dt D_{-+3} (type I_-); Dirichlet flux '
           '= -(beta/m0) sin(beta phi_D) d_t psi_x', check_bulk_currents() and ok and okD)


def check_fused_charge():
    (f, f1, f2, f3), _, Dt = boundary_jets()
    jets = [(f, f1), (f1, f2), (f2, f3)]
    PH, PS = sp.symbols('PH PS')
    phx = -sp.diff(Bd(PH, PS), PH).subs({PH: f, PS: phD})         # fused boundary condition
    psx = f1 + sp.diff(Bd(PH, PS), PS).subs({PH: f, PS: phD})     # second defect condition at psi = phi_D
    php, phm = f1 + phx, f1 - phx
    p2p, p2m = second_derivs(f, f1, f2, phx, Dt(phx, jets))
    Dm, Dp = Dterms(f, phD, php, phm, psx, -psx)
    jsum = j3(php, p2p, f, -1) + j3(phm, p2m, f, 1)
    jdiff = j3(php, p2p, f, -1) - j3(phm, p2m, f, 1)
    extra = -b / m0 * sp.sin(b * phD) * psx
    zero = lambda j, Bt: is_zero(-j + Dt(Bt, jets), [f, phD])
    ok = zero(jsum, Dm + Dp + extra)
    no_extra = not zero(jsum, Dm + Dp)
    diff = not any(zero(jdiff, Dm - Dp + c * psx) for c in (0, -b / m0 * sp.sin(b * phD), b / m0 * sp.sin(b * phD)))
    # the steps of the proof by hand: v = f1, a = f2, F = phi_x(f)
    v, a, F = f1, f2, phx
    Fp = sp.diff(F, f)
    bracket = b**2 / (2 * m0) * F * sp.cos(b * f) - b / m0 * Fp * sp.sin(b * f) + b**4 / (4 * m0**3) * F**3
    steps = is_zero(jsum - (-2 * b**2 / m0**3 * Fp * v * a + b**4 / (4 * m0**3) * F * v**3 + bracket * v), [f, phD])
    steps &= is_zero(sp.diff(F, f, 2) + b**2 / 4 * F, [f, phD])
    H = sp.expand(Dm + Dp + extra)
    steps &= is_zero(H.coeff(v, 2) + b**2 / m0**3 * Fp, [f, phD]) and is_zero(H.coeff(v, 1), [f, phD])
    steps &= is_zero(sp.diff(H.coeff(v, 0), f) - bracket, [f, phD])
    lin = is_zero(sp.expand(Dm + Dp).coeff(v, 1) - b / m0 * sp.sin(b * phD), [f, phD])   # what the Dirichlet term cancels
    report('fused boundary conserves int(q_-3 + q_+3) + D_-3 + D_+3 - (beta/m0) sin(beta phi_D) psi_x; steps of the '
           'proof (flux in v, a; F\'\' = -beta^2 F/4; boundary term = -(beta^2/m0^3) F\' v^2 + U(f), U\' = bracket); '
           'controls (no Dirichlet term; q_-3 - q_+3) fail', ok and no_extra and diff and steps and lin)


# ---------------------------------------------------------------- 7. the boost heuristic

def check_boost_heuristic():
    beta, lam, th, B, bb, bLY = sp.symbols('beta lambda theta_d B b b_LY', positive=True)
    h = beta**2 / (8 * sp.pi)
    lam_of_beta = 8 * sp.pi / beta**2 - 1
    ok1 = sp.simplify(h - 1 / (lam_of_beta + 1)) == 0
    hl = 1 / (lam + 1)
    ok2 = sp.simplify((1 - hl) * th * (lam + 1) - lam * th) == 0               # theta_cl (lam+1) = lam theta_d
    ok3 = sp.simplify((1 / (1 - hl)).subs(lam, -2 / B) - (1 - B / 2)) == 0      # theta_d = (1 - B/2) theta_cl
    Bb = 2 * bb**2 / (8 * sp.pi + bb**2)
    ser = sp.series(1 - Bb / 2, bb, 0, 4).removeO()
    ok4 = sp.simplify(ser - (1 - bb**2 / (8 * sp.pi))) == 0
    thb = sp.I * sp.pi * (3 - bLY) / 6
    ok5 = sp.simplify(-sp.I * (bLY + 3) * sp.pi / 5 - sp.Rational(6, 5) * (thb - sp.I * sp.pi)) == 0
    ok6 = 1 - sp.Rational(-1, 5) == sp.Rational(6, 5)
    betahat2 = beta**2 / (4 * sp.pi)
    ok7 = sp.simplify((1 - betahat2 / 2) - (1 - h)) == 0
    report('boost heuristic: h = 1/(lam+1), sinh-Gordon (1-B/2) and one loop, BHW 6/5 = 1-h, KL d(g) = 1-h',
           all([ok1, ok2, ok3, ok4, ok5, ok6, ok7]))


# ---------------------------------------------------------------- 8. folding

def check_folding():
    ok_lin = True
    for h in range(3, 7):
        P = np.roll(np.eye(h), 1, axis=0)                          # pi: e_j -> e_{j+1}
        R = np.eye(h)[[(-j) % h for j in range(h)]]                 # rho: e_j -> e_{-j}
        ok_lin &= np.allclose(R @ P @ R, P.T)
        om = np.exp(2j * np.pi / h)
        for k in range(1, h):
            v = om ** (-k * np.arange(h))                          # P v = om^k v
            ok_lin &= np.allclose(P @ v, om**k * v) and np.allclose(P @ (R @ v), om**(-k) * (R @ v))
            Dk, DTk, Dhk = 2 / (1 - om**k), 2 / (1 - om**-k), 2 / (1 - om**(h - k))
            ok_lin &= np.isclose(DTk, Dhk) and np.isclose(Dk, 1 + 1j / np.tan(np.pi * k / h))
    x, g, u = sp.symbols('x gamma u')
    E = sp.exp(sp.I * g)
    Tk = lambda xx: (1 + xx * E) / (1 + xx / E)
    Thk = lambda xx: ((1 + xx * E) / (1 + xx / E)).subs(g, sp.pi - g)
    ok_hk = sp.simplify(Thk(x) - (1 - x / E) / (1 - x * E)) == 0
    ratio = sp.simplify(Tk(x) / Thk(x))
    ok_ratio = sp.simplify(ratio - (1 - x**2 * E**2) / (1 - x**2 / E**2)) == 0
    x1 = sp.Symbol('x1')
    num = sp.factor(sp.numer(sp.together(ratio.subs(x, x1) - ratio)))
    ok_pm = sp.simplify(num.subs(x1, x)) == 0 and sp.simplify(num.subs(x1, -x)) == 0
    ctrl = all(sp.simplify(num.subs({x1: c * x, g: sp.pi / 3, x: sp.Rational(7, 10)})) != 0
               for c in (2, sp.Rational(1, 2), sp.I, -3))
    ok_transp = sp.simplify(Tk(x) * Thk(-x) - 1) == 0
    prod = (Tk(x) * Thk(x)).subs(x, sp.exp(u))
    ok_f = sp.simplify((prod - (sp.sinh(u) - sp.I * sp.sin(g)) / (sp.sinh(u) + sp.I * sp.sin(g))).rewrite(sp.exp)) == 0
    zz = (Tk(x) * Thk(x)).subs(x, sp.exp(u - sp.I * sp.pi / 2))
    ok_z = sp.simplify((zz - (sp.cosh(u) + sp.sin(g)) / (sp.cosh(u) - sp.sin(g))).rewrite(sp.exp)) == 0
    z1 = sp.coth(u / 2 - sp.I * sp.pi / 12) * sp.coth(u / 2 + sp.I * sp.pi / 12)
    vals = [complex(sp.N((zz.subs(g, sp.pi / 3) - z1).subs(u, uu))) for uu in (0.37, -1.2 + 0.4j, 2.1 - 0.3j)]
    ok_tz = max(abs(v) for v in vals) < 1e-12
    # 2k = h (gamma = pi/2): no condition, and the delay is coth^2(u/2), with a pole at u = 0
    ok_h2 = sp.simplify(ratio.subs(g, sp.pi / 2) - 1) == 0
    ok_h2 &= sp.simplify((zz.subs(g, sp.pi / 2) - sp.coth(u / 2)**2).rewrite(sp.exp)) == 0
    report('folding: mirror transmits k as h-k; T_k/T_{h-k} Mobius in x^2; folds iff sigma_2 = +-sigma_1 '
           '(4 controls fail); transparency; 1/f_{k/h}; soliton delay; Tzitzeica z_1; 2k = h',
           all([ok_lin, ok_hk, ok_ratio, ok_pm, ctrl, ok_transp, ok_f, ok_z, ok_tz, ok_h2]))


# ---------------------------------------------------------------- 9. chirality

def check_chirality():
    # at m0 = 0, x = 0: phi_x = phi_L' - phi_R', phi_t = phi_L' + phi_R', and the same for psi
    pl, pr, sl, sr = sp.symbols('pl pr sl sr')
    ok_types = True
    for D, want in [(1, {sl: pl, sr: -pr}), (-1, {sl: -pl, sr: pr})]:
        sol = sp.solve([(pl - pr) - D * (sl + sr), (sl - sr) - D * (pl + pr)], [sl, sr], dict=True)[0]
        ok_types &= all(sp.simplify(sol[k] - v) == 0 for k, v in want.items())
    t, x = sp.symbols('t x')
    lm = sp.Function('F')(t + x)
    E_ = (sp.diff(lm, t)**2 + sp.diff(lm, x)**2) / 2
    P_ = -sp.diff(lm, t) * sp.diff(lm, x)                         # P of prp-sg-defect-momentum
    ok_ep = (sp.simplify(E_ + P_) == 0
             and sp.simplify(E_ - P_ - (sp.diff(lm, t) + sp.diff(lm, x))**2 / 2) == 0)
    lam, th, thd, zeta = sp.symbols('lambda theta theta_d zeta', real=True)
    y = sp.exp(lam * th)                                           # e_1 v_- = y v_+, f_1 v_+ = v_-/y
    ok_spin = (sp.simplify(y.subs(th, th + zeta) / y - sp.exp(lam * zeta)) == 0
               and sp.simplify((1 / y).subs(th, th + zeta) * y - sp.exp(-lam * zeta)) == 0)
    Y = sp.exp(lam * (th - thd))
    qn, qs = sp.symbols('q_n q_s', nonzero=True)
    blk = sp.Matrix([[1 / qn, -sp.I * Y * qn * qs], [-sp.I * Y / qn, qn * qs]])   # A(Y) D_n of the sector
    ok_twist = blk.subs(Y, 0) == sp.diag(1 / qn, qn * qs)
    ok_conv = all(sp.limit((blk[i, i] / blk[i, 1 - i]).subs({lam: 1, thd: 0}), th, sp.oo) == 0 for i in range(2))
    report('chirality: massless I_- passes left-movers, I_+ right-movers; left-movers carry E - P only; '
           'f_1 ~ 1/y, e_1 ~ y; the I_- block is a twist for theta << theta_d, converting for theta >> theta_d',
           all([ok_types, ok_ep, ok_spin, ok_twist, ok_conv]))


if __name__ == '__main__':
    check_sigma_integral()
    check_proof_steps()
    check_scalar_identity()
    check_matrix_part()
    check_defect_flux()
    check_fused_charge()
    check_boost_heuristic()
    check_folding()
    check_chirality()
