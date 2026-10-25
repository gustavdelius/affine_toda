"""The sine-Gordon defect scalar factor as an integral, its poles, the lightest breather, and the two Weyl pairs.

Book ids: eq-sgd-scalar, prp-sgd-scalar-integral, eq-sgd-scalar-integral, prp-sg-breather-transmission (ch. 38);
prp-sgd-pole-ladder, prp-sgd-odd-jump (ch. 38, sec-sg-defect-spectrum); prp-dq-two-weyl-pairs (ch. 41).

Written without defect_common.py.  Conventions of the book: q = -exp(-i pi lambda) = exp(-8 pi^2 i/beta_SG^2),
Part VI coproduct Delta(f_j) = f_j (x) k_j^-1 + 1 (x) f_j, y = exp(lambda theta), y_d = exp(lambda theta_d),
u = theta - theta_d, Y = exp(lambda u).  Phi_G is (-@eq-sgd-scalar) as printed, evaluated here as the convergent
product of the ratios r_k/rbar_k (each ratio is 1 + O(k^-2)) with mpmath.nprod-style extrapolation of the log-sum.

Checks:
 1. The integral (eq-sgd-scalar-integral) equals the Gamma product at real and complex u, for three couplings,
    up to Im u = pi/2 + pi/(2 lambda).                                        [prp-sgd-scalar-integral]
 2. The proof: Malmstén's integral for log Gamma turns log(r/rbar) into Nemes's kernel
    sin(ku) sinh(pi k (1+1/lambda)/2)/(k sinh(pi k/lambda) cosh(pi k/2)); it equals the Gamma ratio, and it equals the
    imaginary part of the integral of check 1 plus the phase of the prefactor (arctan identity).
 3. Unitarity |Phi|^2 (1+Y^2) = 1 and crossing Phi(u) Phi(u+i pi)(1+Y^2) = 1 of the integral at complex u;
    Phi -> 1 as u -> -infinity; no zeros or poles on a grid of the strip -pi/2lam < Im u < pi + pi/2lam.
 4. Poles (prp-sgd-pole-ladder): eps Phi_G(u_n + eps) -> -1/(2 lambda Phi(u_n + i pi)) at u_n = -i(2n+1)pi/(2 lambda),
    every n < lambda, with 1/eps scaling; finite at the midpoints u = -i n pi/lambda (eps Phi -> 0 linearly).
 5. The origin of the poles: for even n the pole sits in the Gamma ratio r/rbar, for odd n in the prefactor
    1/(1 + iY) (the Gamma ratio is finite there).  So the poles with -pi/2 < Im u < 0 number #{n: 2n+1 < lambda},
    about lambda/2, while the Gamma functions alone give about lambda/4 (Nemes, arXiv:0909.3268, eq. (67)).
 6. Widths: Gamma_n = -2 Im M cosh(theta_d + u_n) = m_{2n+1} sinh theta_d, m_k = 2M sin(k pi/2lambda);
    Gamma_n/Gamma_0 -> 2n + 1 as lambda -> infinity, with deviations that fall like 1/lambda^2.
 7. Classical odd jump (prp-sgd-odd-jump): the type-I_- defect conditions (eq-sg-defect-conditions, with B_d of
    eq-sg-defect-potential) linearized about (phi, psi) = (0, 2 pi/beta) give, for phi = A e^{gt + kx}, psi - 2pi/beta
    = B e^{gt - kx}, k^2 = g^2 + m0^2, the determinant proportional to k (k - m0 cosh theta_d): modes
    g = +-m0 sinh theta_d, k = m0 cosh theta_d; the growing one has A = 0 for theta_d != 0.  At theta_d = 0 (g = 0,
    k = m0) the whole matching matrix vanishes, so A and B are independent static modes.  Control: about (0, 0) the determinant is
    proportional to k (k + m0 cosh theta_d), with no localized mode.
 8. The emitted soliton of rapidity theta_d (prp-sg-delay (2)) starts as the growing mode: (psi - 2 pi/beta) e^{-gt}
    at x = 0 tends to a constant as t -> -infinity; the pair (0, soliton) solves the full defect conditions.
 9. Inverted oscillator: H = (p^2 - g^2 x^2)/2 = g(XP + PX)/2 with [X, P] = [x, p], and H X^n = -i hbar g (n+1/2) X^n.
10. The book's transmission matrix (eq-sg-transmission-matrix) intertwines U_q(b^-) on V(theta) (x) W; the vector
    v+ v- - v- v+ of V(theta - i u_1/2) (x) V(theta + i u_1/2) is invariant; T_12 T_23 maps it, times |n>, to
    q^s (1 + q Y^2) |n> (x) omega for every n (sympy).  Control: at y1/y2 = -q the |n+-1> terms survive.
11. Lightest breather: Phi_G(u - iu_1/2) Phi_G(u + iu_1/2)(1 + qY^2) = T_cl = (i - sinh u)/cosh u, both printed forms
    of eq-shg-linear-transmission, three couplings, complex u.  Control: u_1 -> pi - 1.1 pi/lambda fails.
12. The proof of 11: crossing gives Phi(v - eps)/Phi(v + eps), v = u + i pi/2, eps = i pi/(2 lambda); the integral
    gives log T_B1 = i pi/2 + i int_0^inf sin(ku)/(k cosh(pi k/2)) dk = i pi/2 + i gd(u).
13. Two Weyl pairs (prp-dq-two-weyl-pairs): with [phi(0), psi(0)] = -i hbar realized on functions (phi = x,
    psi = i hbar d/dx), U = e^{-b(phi-psi)/2sqrt2}, W = e^{b phi/sqrt2} have U W = qhat W U, qhat = e^{i hbar b^2/4};
    U_s = e^{-sqrt2 pi (phi-psi)/(hbar b)}, X = e^{2 sqrt2 pi phi/(hbar b)} have U_s X = q_s X U_s,
    q_s = e^{4 pi^2 i/(hbar b^2)}; cross pairs commute up to -1; U^2, W^2 commute with U_s, X; at b = i beta/sqrt2,
    hbar = 1: qhat = e^{-i beta^2/8}, q_s = -e^{-i pi lambda}; U_s = q^{n+s/2} on the jump 2 pi(2n+s)/beta, X shifts n by 1.
14. Faddeev's normalization (math/9912078, (14)-(18)): his (16) follows from (14); with his printed (18) the two Weyl
    pairs do not commute for generic hbar; with hbar alpha beta~ = -hbar alpha~ beta = 1 they commute, and
    q = e^{i pi tau} gives q~ = e^{-i pi/tau}, as in his introduction and (22) (ch. 41 notes).
Runtime: about 5 minutes.
"""
import mpmath as mp
import numpy as np
import sympy as sp

mp.mp.dps = 20
PI = mp.pi


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


# ---------------------------------------------------------------- the scalar factor as printed (eq-sgd-scalar)

def log_ratio_k(k, z, lam):
    """log of r_k(z)/rbar_k(z), z = i lam theta/(2 pi): the k-th factor of r over that of rbar."""
    a1, a2 = k*lam + mp.mpf(1)/4, (k + 1)*lam + mp.mpf(3)/4
    b1, b2 = (k + mp.mpf(1)/2)*lam + mp.mpf(1)/4, (k + mp.mpf(1)/2)*lam + mp.mpf(3)/4
    num = mp.gamma(a1 - z)*mp.gamma(a2 - z)*mp.gamma(b1 + z)*mp.gamma(b2 + z)
    den = mp.gamma(b1 - z)*mp.gamma(b2 - z)*mp.gamma(a1 + z)*mp.gamma(a2 + z)
    return mp.log(num/den)


def gamma_ratio(u, lam, K=12):
    """r/rbar of (-@eq-sgd-scalar): the first K factors exactly, the tail by nsum (Richardson on a 1/k^2 tail)."""
    z = 1j*lam*u/(2*PI)
    head = mp.fsum(log_ratio_k(k, z, lam) for k in range(K))
    tail = mp.nsum(lambda k: log_ratio_k(k, z, lam), [K, mp.inf])
    return mp.exp(head + tail)


def phi_gamma(u, lam):
    return mp.exp(1j*PI*(1 + lam)/4)/(1 + 1j*mp.exp(lam*u))*gamma_ratio(u, lam)


# ---------------------------------------------------------------- the integral (eq-sgd-scalar-integral)

def log_phi_int(u, lam):
    """-lam u/2 + i pi lam/4 - ln 2/2 + int_0^inf [cosh(pi k/2 + iku)/cosh(pi k/2) - 1] dk/(2k sinh(pi k/2lam))."""
    u = mp.mpc(u)

    def f(k):
        # cosh(a + b)/cosh(a) - 1, a = pi k/2, b = iku: for k < 1 as 2 sinh^2(b/2) + tanh(a) sinh(b), which has no
        # cancellation at small k; for k >= 1 as (e^b + e^{-2a-b})/(1 + e^{-2a}) - 1, which has none at large k
        if k == 0:
            return lam*(1j*PI*u - u**2)/(2*PI)
        b = 1j*k*u
        if k < 1:
            num = 2*mp.sinh(b/2)**2 + mp.tanh(PI*k/2)*mp.sinh(b)
        else:
            num = (mp.exp(b) + mp.exp(-PI*k - b))/(1 + mp.exp(-PI*k)) - 1
        return num/(2*k*mp.sinh(PI*k/(2*lam)))
    # the integrand decays like e^{-rk} (the three exponentials of the proof of prp-sgd-scalar-integral); cut off
    # where e^{-rk} < 10^-26, in steps shorter than a period of e^{iku}
    r = min(PI/(2*lam), PI/(2*lam) + mp.im(u), PI + PI/(2*lam) - mp.im(u))
    K = 60/r
    steps = int(mp.ceil(K*max(1, abs(mp.re(u)))/2)) + 1
    return -lam*u/2 + 1j*PI*lam/4 - mp.log(2)/2 + mp.quad(f, mp.linspace(0, K, steps + 1))


def phi_int(u, lam):
    return mp.exp(log_phi_int(u, lam))


LAMS = [mp.sqrt(2) + 1, mp.e + 1, 3*mp.sqrt(2)]       # generic couplings (lambda irrational), all > 1


def rel(a, b):
    return abs(a - b)/abs(b)


def check_integral_equals_gamma():
    errs = []
    for lam in LAMS:
        top = PI/2 + PI/(2*lam)
        for u in [mp.mpf('0.37'), mp.mpf('-1.3'), mp.mpf('0.8') + 0.3j*PI/(2*lam), mp.mpf('-0.4') - 0.7j*PI/(2*lam),
                  mp.mpf('0.2') + 0.9j*top, mp.mpf('1.1') + 0.5j*PI]:
            errs.append(rel(phi_int(u, lam), phi_gamma(u, lam)))
    report(f'1. integral representation = Gamma product (3 couplings, 6 points each, real and complex u up to '
           f'Im u = 0.9(pi/2 + pi/2lam)); max rel. error {float(max(errs)):.1e}', max(errs) < 1e-12)


def check_malmsten_step():
    """The Malmstén kernel equals log(r/rbar) for real u, and the crossing form of the integral."""
    errs_gamma, errs_identity = [], []
    for lam in LAMS[:2]:
        for u in [mp.mpf('0.37'), mp.mpf('-1.3'), mp.mpf('2.2')]:
            nemes = mp.quad(lambda k: mp.sin(k*u)*mp.sinh(PI*k*(1 + 1/lam)/2)/(k*mp.sinh(PI*k/lam)*mp.cosh(PI*k/2)),
                            [0, 1, 4, 16, mp.inf])
            errs_gamma.append(abs(mp.exp(1j*nemes) - gamma_ratio(u, lam)))
            # phase of the prefactor + Nemes kernel = pi lam/4 + the tanh kernel of the crossing form
            tanh_kernel = mp.quad(lambda k: mp.tanh(PI*k/2)*mp.sin(k*u)/(2*k*mp.sinh(PI*k/(2*lam))),
                                  [0, 1, 4, 16, mp.inf])
            lhs = PI*(1 + lam)/4 - mp.atan(mp.exp(lam*u)) + nemes
            errs_identity.append(abs(lhs - (PI*lam/4 + tanh_kernel)))
    report(f'2. Malmsten step: exp(i x Nemes kernel) = r/rbar (err {float(max(errs_gamma)):.1e}); '
           f'prefactor phase + Nemes kernel = pi lam/4 + tanh kernel (err {float(max(errs_identity)):.1e})',
           max(errs_gamma) < 1e-15 and max(errs_identity) < 1e-15)


def check_unitarity_crossing_strip():
    lam = LAMS[0]
    errs = []
    for u in [mp.mpf('0.3'), mp.mpf('-0.9'), mp.mpf('0.5') + 0.2j*PI/(2*lam), mp.mpf('-0.2') - 0.6j*PI/(2*lam)]:
        Y2 = mp.exp(2*lam*u)
        errs.append(abs(phi_int(u, lam)*phi_int(u + 1j*PI, lam)*(1 + Y2) - 1))
        if mp.im(u) == 0:
            errs.append(abs(abs(phi_int(u, lam))**2*(1 + Y2) - 1))
    lim = abs(phi_int(-30, lam) - 1)
    # grid of the strip: |log Phi| finite and smooth -- no zeros or poles (numerical support for the analytic proof)
    lo, hi = -PI/(2*lam), PI + PI/(2*lam)
    vals = [abs(phi_int(x + 1j*(lo + (hi - lo)*t), lam)) for x in [-2, -0.5, 0, 0.5, 2]
            for t in [0.1, 0.3, 0.5, 0.7, 0.9]]
    ok_grid = all(mp.mpf('1e-6') < v < mp.mpf('1e6') for v in vals)
    report(f'3. unitarity and crossing of the integral (err {float(max(errs)):.1e}); Phi(-30) - 1 = {float(lim):.1e}; '
           f'|Phi| in [{float(min(vals)):.2g}, {float(max(vals)):.2g}] on a 5x5 grid of the strip',
           max(errs) < 1e-12 and lim < 1e-12 and ok_grid)


def poles_of(lam):
    return [-1j*(2*n + 1)*PI/(2*lam) for n in range(int(mp.floor(lam)) + 1) if n < lam]


def check_pole_ladder():
    worst, scal = 0, 0
    for lam in LAMS:
        for n, un in enumerate(poles_of(lam)):
            res_pred = -1/(2*lam*phi_int(un + 1j*PI, lam))
            e1, e2 = mp.mpf('1e-8'), mp.mpf('1e-9')
            r1 = e1*phi_gamma(un + e1, lam)
            r2 = e2*phi_gamma(un + e2, lam)
            res = (e1*r2 - e2*r1)/(e1 - e2)               # eps*Phi = Res + c eps + O(eps^2): remove the linear term
            worst = max(worst, rel(res, res_pred))
            scal = max(scal, abs(r1 - r2)/abs(r2))        # a simple pole: eps*Phi converges, error O(eps)
    # midpoints between the poles: finite, so eps*Phi -> 0 linearly
    lam = LAMS[2]
    mids = [-1j*n*PI/lam for n in range(1, 4)]
    lin = all(abs(mp.mpf('1e-8')*phi_gamma(m + mp.mpf('1e-8'), lam)) < mp.mpf('1e-7') for m in mids)
    report(f'4. simple poles of the Gamma product exactly at u_n = -i(2n+1)pi/2lam, n < lambda (3 couplings, '
           f'{sum(len(poles_of(l)) for l in LAMS)} poles), residues -1/(2 lam Phi(u_n + i pi)) to {float(worst):.1e}; '
           f'finite between them', worst < 1e-10 and scal < 1e-6 and lin)


def check_pole_origin():
    lam = mp.mpf(21) + mp.sqrt(2)/3                     # lambda ~ 21.47: 11 poles with -pi/2 < Im u < 0
    in_gamma, in_prefactor = [], []
    e1, e2 = mp.mpf('1e-9'), mp.mpf('1e-10')
    for n, un in enumerate(poles_of(lam)):
        if mp.im(un) <= -PI/2:
            continue
        # a pole of r/rbar keeps eps*(r/rbar) fixed as eps -> eps/10; a finite value makes it fall tenfold
        ratio = abs(e1*gamma_ratio(un + e1, lam))/abs(e2*gamma_ratio(un + e2, lam))
        (in_gamma if abs(ratio - 1) < 0.01 else in_prefactor).append(n)
        if not (abs(ratio - 1) < 0.01 or abs(ratio - 10) < 0.1):
            in_gamma.append(None)                       # neither: fail
    upper = [n for n in range(int(lam) + 1) if 2*n + 1 < lam]
    expect_gamma = [n for n in upper if n % 2 == 0]
    expect_pref = [n for n in upper if n % 2 == 1]
    report(f'5. lambda = {float(lam):.3f}: {len(upper)} poles with -pi/2 < Im u < 0 (#(2n+1 < lambda) ~ lambda/2); '
           f'{len(in_gamma)} in the Gamma ratio (n even, ~lambda/4 as Nemes counted), {len(in_prefactor)} in 1/(1+iY)',
           in_gamma == expect_gamma and in_prefactor == expect_pref)


def check_widths():
    th_d = mp.mpf('0.7')
    ok = True
    devs = []
    for lam in [mp.mpf(10) + mp.sqrt(2), mp.mpf(100) + mp.sqrt(2), mp.mpf(1000) + mp.sqrt(2)]:
        M = mp.mpf(1)
        m = lambda k: 2*M*mp.sin(k*PI/(2*lam))
        for n in range(4):
            un = -1j*(2*n + 1)*PI/(2*lam)
            gam = -2*mp.im(M*mp.cosh(th_d + un))
            ok &= abs(gam - m(2*n + 1)*mp.sinh(th_d)) < mp.mpf('1e-15')
        gam0 = m(1)*mp.sinh(th_d)
        devs.append(max(abs(m(2*n + 1)*mp.sinh(th_d)/gam0 - (2*n + 1)) for n in range(4)))
    falls = devs[1] < devs[0]/50 and devs[2] < devs[1]/50
    report(f'6. widths Gamma_n = m_(2n+1) sinh theta_d exactly; max |Gamma_n/Gamma_0 - (2n+1)|, n <= 3, at lambda ~ '
           f'10, 100, 1000: {", ".join(f"{float(d):.1e}" for d in devs)}', ok and falls)


# ---------------------------------------------------------------- classical odd jump

def linearized_determinant(odd):
    t, x, beta, m0, thd, A, B, g, k = sp.symbols('t x beta m0 theta_d A B gamma kappa')
    sigma = sp.exp(-thd)
    ph, ps = sp.symbols('phi psi')
    Bd = -2*m0/beta**2*(sigma*(sp.cos(beta*(ph + ps)/2) - 1) + (sp.cos(beta*(ph - ps)/2) - 1)/sigma)
    psi0 = 2*sp.pi/beta if odd else 0
    phi_f = A*sp.exp(g*t + k*x)
    psi_f = psi0 + B*sp.exp(g*t - k*x)
    # eq-sg-defect-conditions, D = 1: phi_x = psi_t - dB/dphi, psi_x = phi_t + dB/dpsi, at x = 0; linear order in A, B
    eqs = []
    for lhs, rhs in [(sp.diff(phi_f, x), sp.diff(psi_f, t) - sp.diff(Bd, ph)),
                     (sp.diff(psi_f, x), sp.diff(phi_f, t) + sp.diff(Bd, ps))]:
        e = (lhs - rhs).subs({ph: phi_f, ps: psi_f}).subs(x, 0)
        eps = sp.Symbol('eps')
        e = sp.series(e.subs({A: eps*A, B: eps*B}), eps, 0, 2).removeO().coeff(eps)
        eqs.append(sp.expand(sp.simplify(e*sp.exp(-g*t))))
    Mtx = sp.Matrix([[sp.diff(e, A), sp.diff(e, B)] for e in eqs])
    det = sp.simplify(Mtx.det().subs(g**2, k**2 - m0**2))
    return det, Mtx, (g, k, m0, thd, beta)


def check_odd_jump():
    det, Mtx, (g, k, m0, thd, beta) = linearized_determinant(odd=True)
    target = k*(k - m0*sp.cosh(thd))
    ratio = sp.simplify((det/target).rewrite(sp.exp))
    const = ratio.free_symbols <= {beta} and ratio != 0
    grow = Mtx.subs({g: m0*sp.sinh(thd), k: m0*sp.cosh(thd)})
    # growing mode: the first equation forces A = 0 (theta_d != 0), the second is then empty
    a_zero = sp.simplify(grow[0, 1].rewrite(sp.exp)) == 0 and sp.simplify(grow[1, 1].rewrite(sp.exp)) == 0
    a_zero = a_zero and sp.simplify(grow[1, 0].rewrite(sp.exp) + 2*m0*sp.sinh(thd).rewrite(sp.exp)) == 0
    # theta_d = 0: gamma = 0, kappa = m0, and the whole matching matrix vanishes (two independent static zero modes)
    static = all(sp.simplify(e) == 0 for e in Mtx.subs({thd: 0, g: 0, k: m0}))
    det_e, _, (g2, k2, m02, thd2, _) = linearized_determinant(odd=False)
    ratio_e = sp.simplify((det_e/(k2*(k2 + m02*sp.cosh(thd2)))).rewrite(sp.exp))
    control = ratio_e.free_symbols <= {beta} and ratio_e != 0
    report(f'7. odd jump: determinant = {sp.factor(ratio)} * kappa (kappa - m0 cosh theta_d); growing mode '
           f'gamma = m0 sinh theta_d has A = 0 for theta_d != 0; at theta_d = 0 the matching matrix vanishes (A, B '
           f'independent static modes); control (0, 0): kappa (kappa + m0 cosh theta_d), no localized mode',
           bool(const and a_zero and static and control))


def check_emission_tail():
    beta, m0, thd = mp.mpf('1.3'), mp.mpf(1), mp.mpf('0.6')
    sigma = mp.exp(-thd)

    def psi(t, x):
        # prp-sg-backlund with theta = theta_d and Im xi = pi/2: psi -> 2pi/beta as t -> -infinity at x = 0
        E = mp.exp(m0*(x*mp.cosh(thd) - t*mp.sinh(thd)) + 1j*PI/2)
        return mp.re(2j/beta*mp.log((1 - E)/(1 + E)))      # = (4/beta) arctan(-iE): 2pi/beta before, 0 after

    vals = [(psi(t, 0) - 2*PI/beta)*mp.exp(-m0*mp.sinh(thd)*t) for t in [-20, -30, -40]]
    const = abs(vals[2] - vals[1]) < mp.mpf('1e-8')*abs(vals[1])
    # the pair (0, soliton) solves eq-sg-defect-conditions (type I_-) at x = 0
    res = 0
    for t in [-1, 0, 0.7]:
        psi_t = mp.diff(lambda tt: psi(tt, 0), t)
        psi_x = mp.diff(lambda xx: psi(t, xx), 0)
        dB_dphi = -2*m0/beta**2*(-sigma*mp.sin(beta*psi(t, 0)/2)*beta/2 - mp.sin(-beta*psi(t, 0)/2)*beta/(2*sigma))
        dB_dpsi = -2*m0/beta**2*(-sigma*mp.sin(beta*psi(t, 0)/2)*beta/2 + mp.sin(-beta*psi(t, 0)/2)*beta/(2*sigma))
        res = max(res, abs(0 - (psi_t - dB_dphi)), abs(psi_x - (0 + dB_dpsi)))
    report(f'8. the emitted soliton of rapidity theta_d starts as the growing mode: (psi - 2pi/beta) e^(-gamma t) -> '
           f'{float(vals[2]):.6f}; (0, soliton) solves the defect conditions (residual {float(res):.1e})',
           const and res < 1e-8)


def check_inverted_oscillator():
    gam, hb = sp.symbols('gamma hbar', positive=True)
    x, p = sp.symbols('x p', commutative=False)
    X, n = sp.symbols('X n')
    # canonical change X = (p + g x)/sqrt(2g), P = (p - g x)/sqrt(2g): [X, P] = ([p, x] g - g [x, p])/(2g) = [x, p]
    c_xp = sp.Symbol('c')                               # [x, p] = c, [p, x] = -c
    ok1 = sp.simplify((-gam*(-c_xp) + gam*c_xp)/(2*gam) - c_xp) == 0     # [p + gx, p - gx]/2g
    # XP + PX, operator ordering kept: the cross terms cancel without using [x, p]
    sym = sp.expand((p + gam*x)*(p - gam*x) + (p - gam*x)*(p + gam*x))/(2*gam)
    ok2 = sp.expand(gam/2*sym - (p**2 - gam**2*x**2)/2) == 0
    # P = -i hbar d/dX: H X^n = (g/2)(X P + P X) X^n = -i hbar g (n + 1/2) X^n
    f = X**n
    Hf = gam/2*(X*(-sp.I*hb*sp.diff(f, X)) + (-sp.I*hb)*sp.diff(X*f, X))
    ok3 = sp.simplify(Hf/f + sp.I*hb*gam*(n + sp.Rational(1, 2))) == 0
    report('9. inverted oscillator: [X, P] = [x, p]; H = (p^2 - g^2 x^2)/2 = g(XP + PX)/2 (noncommutative); '
           'H X^n = -i hbar g (n + 1/2) X^n', ok1 and ok2 and ok3)


# ---------------------------------------------------------------- the bootstrap of the lightest breather

def soliton_rep(y, q):
    """V(theta) in Part VI conventions, basis (v+, v-): e_1 v- = y v+, f_1 v+ = v-/y, e_0 v+ = y v-, f_0 v- = v+/y."""
    k1 = sp.diag(q, 1/q)
    f1 = sp.Matrix([[0, 0], [1/y, 0]])
    f0 = sp.Matrix([[0, 1/y], [0, 0]])
    return {'k1': k1, 'k0': k1.inv(), 'f1': f1, 'f0': f0}


def check_bootstrap_algebra():
    q, Yb, s, mu, yd = sp.symbols('q Y_2 s mu y_d', nonzero=True)
    N = 3                                  # window n = -N..N of W, equations on the interior
    ns = list(range(-N, N + 1))
    dW = len(ns)

    def T_matrix(Y):
        """eq-sg-transmission-matrix (scalar factor dropped): V (x) W -> W (x) V, basis index (a, n) -> (n, a)."""
        Tm = sp.zeros(2*dW, 2*dW)
        idx_in = lambda a, i: a*dW + i                    # V (x) W
        idx_out = lambda i, a: i*2 + a                    # W (x) V
        for i, n in enumerate(ns):
            Tm[idx_out(i, 0), idx_in(0, i)] += q**(-n)
            if i + 1 < dW:
                Tm[idx_out(i + 1, 1), idx_in(0, i)] += -sp.I*Y*q**(-n)
            Tm[idx_out(i, 1), idx_in(1, i)] += q**(n + s)
            if i - 1 >= 0:
                Tm[idx_out(i - 1, 0), idx_in(1, i)] += -sp.I*Y*q**(n + s)
        return Tm

    # intertwining with U_q(b^-): Delta(f_j) = f_j (x) k_j^-1 + 1 (x) f_j, Delta(k) = k (x) k; W: f_1|n> = c|n-1>,
    # f_0|n> = c|n+1>, k_1|n> = q^(2n+s)|n>, c = mu/y_d with mu = i/(q - 1/q); Y = y/y_d
    y = sp.Symbol('y', nonzero=True)
    c = mu/yd
    V = soliton_rep(y, q)
    k1W = sp.diag(*[q**(2*n + s) for n in ns])
    f1W = sp.zeros(dW, dW)
    f0W = sp.zeros(dW, dW)
    for i in range(dW):
        if i - 1 >= 0:
            f1W[i - 1, i] = c
        if i + 1 < dW:
            f0W[i + 1, i] = c
    W = {'k1': k1W, 'k0': k1W.inv(), 'f1': f1W, 'f0': f0W}
    kp = sp.kronecker_product
    Tm = T_matrix(y/yd)
    interior_in = [a*dW + i for a in range(2) for i in range(1, dW - 1)]
    interior_out = [i*2 + a for i in range(2, dW - 2) for a in range(2)]
    res = []
    for j in ['1', '0']:
        lhs = Tm*(kp(V['f' + j], W['k' + j].inv()) + kp(sp.eye(2), W['f' + j]))
        rhs = (kp(W['f' + j], V['k' + j].inv()) + kp(sp.eye(dW), V['f' + j]))*Tm
        D = (lhs - rhs).subs(mu, sp.I/(q - 1/q))
        res += [sp.simplify(D[r, cc]) for r in interior_out for cc in interior_in]
    ok_int = all(r == 0 for r in res)

    # invariance of omega in V(y1) (x) V(y2) at y1 = q y2
    y2 = sp.Symbol('y2', nonzero=True)
    V1, V2 = soliton_rep(q*y2, q), soliton_rep(y2, q)
    omega = sp.Matrix([0, 1, -1, 0])
    inv = all(sp.simplify(e) == 0 for j in ['1', '0']
              for e in (kp(V1['f' + j], V2['k' + j].inv()) + kp(sp.eye(2), V2['f' + j]))*omega)
    inv &= sp.simplify(kp(V1['k1'], V2['k1'])*omega - omega) == sp.zeros(4, 1)

    def bootstrap(ratio):
        """T_12(theta_1) T_23(theta_2) on omega (x) |n>, n = 0, with Y_1 = ratio * Y_2; returns the vector in W (x) V (x) V."""
        Y_1, Y_2 = ratio*Yb, Yb
        T1, T2 = T_matrix(Y_1), T_matrix(Y_2)
        vec = sp.zeros(2*2*dW, 1)                # V (x) V (x) W, index (a, b, i)
        i0 = ns.index(0)
        vec[(0*2 + 1)*dW + i0] = 1
        vec[(1*2 + 0)*dW + i0] = -1
        # T_23: V (x) [V (x) W] -> V (x) [W (x) V]
        P23 = kp(sp.eye(2), T2)
        v1 = P23*vec                             # index (a, i, b) = a*(dW*2) + i*2 + b
        # T_12 acts on (a, i) -> (i, a): reorder (a, i, b) into blocks
        out = sp.zeros(dW*2*2, 1)                # W (x) V (x) V, index (i, a, b)
        for b in range(2):
            sub = sp.Matrix([v1[a*(dW*2) + i*2 + b] for a in range(2) for i in range(dW)])
            img = T1*sub                         # index (i, a) = i*2 + a
            for i in range(dW):
                for a in range(2):
                    out[(i*2 + a)*2 + b] += img[i*2 + a]
        return out, i0

    out, i0 = bootstrap(q)
    expect = sp.zeros(dW*4, 1)
    expect[(i0*2 + 0)*2 + 1] = 1
    expect[(i0*2 + 1)*2 + 0] = -1
    Ysq = q*Yb**2                                # Y^2 = Y_1 Y_2
    ok_b = sp.simplify(out - q**s*(1 + q*Ysq)*expect) == sp.zeros(dW*4, 1)
    out_c, _ = bootstrap(-q)
    survives = any(sp.simplify(out_c[(i*2 + a)*2 + b2]) != 0 for i in [i0 - 1, i0 + 1] for a in range(2)
                   for b2 in range(2))
    report('10. T intertwines U_q(b^-) on V (x) W (window, sympy); omega is invariant at y1 = q y2; '
           'T_12 T_23 (omega (x) |n>) = q^s (1 + q Y^2) |n> (x) omega; control y1 = -q y2 leaves |n+-1> terms',
           ok_int and inv and ok_b and survives)


def t_cl_forms(u, th, thd):
    a = (1j - mp.sinh(u))/mp.cosh(u)
    b = -1j*mp.sinh(u/2 - 1j*PI/4)/mp.sinh(u/2 + 1j*PI/4)
    c = (mp.sinh(thd) + 1j*mp.cosh(th))/(mp.cosh(thd) - 1j*mp.sinh(th))
    return a, b, c


def check_breather():
    errs, ctrl = [], []
    thd = mp.mpf('0.45')
    for lam in LAMS:
        q = -mp.exp(-1j*PI*lam)
        u1 = PI - PI/lam
        for u in [mp.mpf('0.6'), mp.mpf('-1.1'), mp.mpf('0.3') + 0.2j, mp.mpf('-0.5') - 0.15j]:
            tb1 = phi_gamma(u - 1j*u1/2, lam)*phi_gamma(u + 1j*u1/2, lam)*(1 + q*mp.exp(2*lam*u))
            a, b, c = t_cl_forms(u, u + thd, thd)
            errs.append(max(rel(tb1, a), rel(b, a), rel(c, a)))
        u = mp.mpf('0.6')
        u1c = PI - 1.1*PI/lam
        tbc = phi_gamma(u - 1j*u1c/2, lam)*phi_gamma(u + 1j*u1c/2, lam)*(1 + q*mp.exp(2*lam*u))
        ctrl.append(rel(tbc, t_cl_forms(u, u + thd, thd)[0]))
    report(f'11. T_B1 = Phi(u - iu_1/2) Phi(u + iu_1/2)(1 + qY^2) = (i - sinh u)/cosh u = both forms of T_cl '
           f'(3 couplings, complex u; err {float(max(errs)):.1e}); control u_1 -> pi - 1.1pi/lam: err '
           f'{float(min(ctrl)):.2f}', max(errs) < 1e-12 and min(ctrl) > 1e-2)


def check_breather_proof():
    errs = []
    for lam in LAMS[:2]:
        q = -mp.exp(-1j*PI*lam)
        u1 = PI - PI/lam
        for u in [mp.mpf('0.6'), mp.mpf('-1.1')]:
            v, eps = u + 1j*PI/2, 1j*PI/(2*lam)
            direct = phi_gamma(u - 1j*u1/2, lam)*phi_gamma(u + 1j*u1/2, lam)*(1 + q*mp.exp(2*lam*u))
            via_crossing = phi_int(v - eps, lam)/phi_int(v + eps, lam)
            gd = mp.quad(lambda k: mp.sin(k*u)/(k*mp.cosh(PI*k/2)), [0, 1, 4, 16, mp.inf])
            closed = mp.exp(1j*PI/2 + 1j*gd)
            errs += [rel(via_crossing, direct), rel(closed, direct), abs(gd - 2*mp.atan(mp.tanh(u/2)))]
    report(f'12. proof of 11: = Phi(v - eps)/Phi(v + eps) by crossing; = exp(i pi/2 + i int sin(ku)/(k cosh(pi k/2)) '
           f'dk) = i e^(i gd u) (err {float(max(errs)):.1e})', max(errs) < 1e-12)


# ---------------------------------------------------------------- two Weyl pairs in the defect values (ch. 41)

def check_two_weyl_pairs():
    x, hb = sp.symbols('x hbar', positive=True)
    b, s_, n_ = sp.symbols('b s n')
    f = sp.Function('f')

    def expo(a1, a2):
        """e^{a1 phi + a2 psi} on functions, phi = x, psi = i hbar d/dx, so [phi, psi] = -i hbar."""
        return lambda g: sp.exp(sp.I*hb*a1*a2/2)*sp.exp(a1*x)*g.subs(x, x + sp.I*hb*a2)

    # it is the exponential: d/ds F_s = (a1 phi + a2 psi) F_s
    a1, a2, s = sp.symbols('a1 a2 s')
    F = sp.exp(sp.I*hb*s**2*a1*a2/2)*sp.exp(s*a1*x)*f(x + sp.I*hb*s*a2)
    gen = a1*x*F + a2*sp.I*hb*sp.diff(F, x)
    ok_exp = sp.simplify(sp.diff(F, s) - gen) == 0
    ok_ccr = sp.simplify(x*sp.I*hb*sp.diff(f(x), x) - sp.I*hb*sp.diff(x*f(x), x) + sp.I*hb*f(x)) == 0

    r2 = sp.sqrt(2)
    U = (-b/(2*r2), b/(2*r2))                   # e^{-b(phi - psi)/2sqrt2}
    Wp = (b/r2, 0)                               # e^{b phi/sqrt2}
    Us = (-r2*sp.pi/(hb*b), r2*sp.pi/(hb*b))     # e^{-sqrt2 pi (phi - psi)/(hbar b)}
    Xp = (2*r2*sp.pi/(hb*b), 0)                  # e^{2 sqrt2 pi phi/(hbar b)}
    sq = lambda e: (2*e[0], 2*e[1])

    def phase(A, B):
        """e^A e^B = phase * e^B e^A, computed on f."""
        eA, eB = expo(*A), expo(*B)
        lhs = eA(eB(f(x)))
        rhs = eB(eA(f(x)))
        return sp.simplify(lhs/rhs)

    qhat = phase(U, Wp)
    qs = phase(Us, Xp)
    ok_q = sp.simplify(qhat - sp.exp(sp.I*hb*b**2/4)) == 0 and sp.simplify(qs - sp.exp(4*sp.pi**2*sp.I/(hb*b**2))) == 0
    ok_signs = (sp.simplify(phase(U, Xp) + 1) == 0 and sp.simplify(phase(Us, Wp) + 1) == 0
                and sp.simplify(phase(U, Us) - 1) == 0 and sp.simplify(phase(Wp, Xp) - 1) == 0)
    ok_sq = sp.simplify(phase(sq(U), Xp) - 1) == 0 and sp.simplify(phase(Us, sq(Wp)) - 1) == 0
    # the sine-Gordon values: b = i beta/sqrt2, hbar = 1; lambda = 8 pi/beta^2 - 1
    beta = sp.Symbol('beta', positive=True)
    lam = 8*sp.pi/beta**2 - 1
    qh_sg = sp.simplify(qhat.subs({b: sp.I*beta/r2, hb: 1}))
    qs_sg = sp.simplify(qs.subs({b: sp.I*beta/r2, hb: 1}))
    ok_sg = sp.simplify(qh_sg - sp.exp(-sp.I*beta**2/8)) == 0
    for bv in [1.3, 2.9]:                                   # q_s = -e^{-i pi lambda}, numerically at two couplings
        ok_sg &= abs(complex((qs_sg + sp.exp(-sp.I*sp.pi*lam)).subs(beta, bv).evalf())) < 1e-12
    # U_s on a defect state of jump psi - phi = 2 pi (2n + s)/beta: e^{-sqrt2 pi (phi - psi)/b} = q^{n + s/2}
    expo_eig = Us[0].subs({b: sp.I*beta/r2, hb: 1})*(-2*sp.pi*(2*n_ + s_)/beta)
    ok_eig = sp.simplify(expo_eig - (n_ + s_/2)*(-8*sp.pi**2*sp.I/beta**2)) == 0
    # X_+ lowers phi - psi by 4 pi/beta: [phi - psi, X] = -(4 pi/beta) X, i.e. n -> n + 1
    shift = sp.simplify(Xp[0].subs({b: sp.I*beta/r2, hb: 1})*(-sp.I))       # [phi - psi, a phi] = a [phi - psi, phi]
    ok_shift = sp.simplify(shift + 4*sp.pi/beta) == 0
    report('13. two Weyl pairs from [phi, psi] = -i hbar (realized on functions): qhat = e^(i hbar b^2/4), '
           'q_s = e^(4 pi^2 i/(hbar b^2)); cross pairs commute up to -1, squares exactly; at b = i beta/sqrt2, hbar = 1, '
           'qhat = e^(-i beta^2/8), q_s = -e^(-i pi lambda), U_s = q^(n+s/2) on |n>, X shifts n by 1',
           ok_exp and ok_ccr and ok_q and ok_signs and ok_sq and ok_sg and ok_eig and ok_shift)


def check_faddeev_normalization():
    """Faddeev math/9912078: [Q, P] = 2 pi i hbar (14), u = e^{alpha P}, v = e^{beta Q} (15), uv = q^2 vu (16)."""
    al, be, alt, bet, hb = sp.symbols('alpha beta alpha_t beta_t hbar', nonzero=True)
    PQ = -2*sp.pi*sp.I*hb                                # [P, Q]
    weyl = lambda a, b: sp.exp(a*b*PQ)                   # e^{aP} e^{bQ} = weyl(a, b) e^{bQ} e^{aP}
    q = sp.exp(-sp.I*sp.pi*hb*al*be)                     # (16): q^2 = weyl(alpha, beta), ln q = pi alpha beta hbar/i
    ok16 = sp.simplify(weyl(al, be) - q**2) == 0
    # printed (18), alpha beta_t = alpha_t beta = hbar: the pairs do not commute for generic hbar
    printed = weyl(al, bet).subs(bet, hb/al)
    fails = abs(complex(printed.subs(hb, 0.7).evalf()) - 1) > 1e-3
    # corrected: hbar alpha beta_t = 1, hbar alpha_t beta = -1: they commute, and q = e^{i pi tau} gives qt = e^{-i pi/tau}
    sub = {bet: 1/(hb*al), alt: -1/(hb*be)}
    commute = sp.simplify(weyl(al, bet).subs(sub) - 1) == 0 and sp.simplify(weyl(alt, be).subs(sub) - 1) == 0
    tau = -hb*al*be                                      # q = e^{i pi tau}
    log_qt = (-sp.I*sp.pi*hb*alt*bet).subs(sub)         # ln q~ = pi alpha~ beta~ hbar/i
    ok_qt = sp.simplify(log_qt - (-sp.I*sp.pi/tau)) == 0
    report('14. Faddeev (14)-(16) consistent; printed (18) alpha beta~ = alpha~ beta = hbar leaves a phase '
           'e^(-2 pi i hbar^2); hbar alpha beta~ = -hbar alpha~ beta = 1 commutes and gives q~ = e^(-i pi/tau)',
           ok16 and fails and commute and ok_qt)


if __name__ == '__main__':
    check_integral_equals_gamma()
    check_malmsten_step()
    check_unitarity_crossing_strip()
    check_pole_ladder()
    check_pole_origin()
    check_widths()
    check_odd_jump()
    check_emission_tail()
    check_inverted_oscillator()
    check_bootstrap_algebra()
    check_breather()
    check_breather_proof()
    check_two_weyl_pairs()
    check_faddeev_normalization()
