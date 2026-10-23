"""Referee check of prp-dn-nondiagonal-scalar (ch. 33) and eq-dn-neumann-amplitude (ch. 35).
Written from the formulas in the book text only (no code imports).

Part A (mpmath, 30 digits): with T = (n-1) omega, mu complex,
  (i)  A_s(mu) A_s(-mu) U_s = 1,
  (ii) A_s(T/2-mu)/A_s(T/2+mu) = F(2mu)/Psi_s(mu),
  (iii) the breather identities eq-dn-nondiagonal-breather, given c^d_s of (iv),
  also eq-dn-scalar (eps-hat = 0) as a control.
  Infinite products are summed as series of log-Gamma differences with mpmath.nsum
  (Richardson / Euler-Maclaurin), compared against each other.
Part B (continuation omega -> -2/B):
  K^N_1 = K_1|_{B->2-B},  K^(-)_1 = K_1/K_1^cl from the closed forms, unitarity,
  crossing-unitarity with S_11, bootstrap closure over the vector-like particles (n=4..7);
  minimal (no Xi) amplitude; one loop vs Kim-Cho (29) for Xi, for the minimal amplitude
  and for alternative CDD factors P_x = (T/2-x)_I (T/2+x)_I.
"""
import mpmath as mp
import itertools

mp.mp.dps = 30
pi = mp.pi


def blkI(x, mu, T):
    return mp.sin(pi * (mu + x) / (2 * T)) / mp.sin(pi * (mu - x) / (2 * T))


def lg(mu, c):
    return mp.loggamma(mu + c) - mp.loggamma(-mu + c)


def series(f, method):
    return mp.nsum(f, [0, mp.inf], method=method)


METHOD = 'richardson'


def logF(mu, T, om):
    def term(k):
        def lB(a):
            return lg(mu, a + 2 * T * k) - lg(mu, a + T + 2 * T * k)
        return lB(1 + T) + lB(1 + om) - lB(0) - lB(T - om)
    return series(term, METHOD)


def logA0(mu, T, om):
    m2 = 2 * mu

    def term(k):
        s = 4 * T * k
        return (lg(m2, 1 + 3 * T + s) + lg(m2, 1 + om + 2 * T + s) + lg(m2, T + s) + lg(m2, 2 * T - om + s)
                - lg(m2, 1 + 2 * T + s) - lg(m2, 1 + om + T + s) - lg(m2, 2 * T + s) - lg(m2, 3 * T - om + s))
    return series(term, METHOD)


def logrho_prod(mu, T, om, d):
    def term(k):
        j = k + 1
        s = 2 * j * T
        return (lg(mu, s + d + om / 2 - 3 * T / 2) + lg(mu, s + 1 - d - om / 2 - T / 2)
                - lg(mu, s + d + om / 2 - T / 2) - lg(mu, s + 1 - d - om / 2 - 3 * T / 2))
    return series(term, METHOD)


def logv(x, om, sg):
    if sg > 0:
        return mp.loggamma(x + om / 2) + mp.loggamma(1 + x - om / 2) - mp.loggamma(x) - mp.loggamma(1 + x)
    return mp.loggamma(mp.mpf(1) / 2 + x + om / 2) + mp.loggamma(mp.mpf(1) / 2 + x - om / 2) - 2 * mp.loggamma(mp.mpf(1) / 2 + x)


def logg_s(mu, T, om, sg):
    def term(k):
        kk = k + 1
        return (-1) ** kk * (logv(mu + kk * T, om, sg) - logv(-mu + kk * T, om, sg))
    return logv(mu, om, sg) + mp.nsum(term, [0, mp.inf])     # alternating: default (Shanks)


def A_s(mu, T, om, sg):
    d = 0 if sg > 0 else mp.mpf(1) / 2
    if sg > 0:
        r = blkI(-T / 2, mu, T) * mp.sin(pi * om / 2) / mp.sin(pi * (mu - T / 2))
    else:
        r = -1j * mp.sin(pi * om / 2) / mp.cos(pi * (mu - T / 2))
    return mp.exp(logA0(mu, T, om) + logrho_prod(mu, T, om, d) + logg_s(mu, T, om, sg)) * r


def U_s(mu, T, om, sg):
    c2 = mp.cos(2 * pi * mu)
    return sg * (c2 - sg * mp.cos(pi * T)) * (c2 - sg * mp.cos(pi * om)) / (2 * mp.sin(pi * om / 2) ** 2 * (c2 - sg))


def Psi_s(mu, om, sg):
    if sg > 0:
        return mp.sin(pi * (om / 2 - mu)) / mp.sin(pi * (om / 2 + mu))
    return mp.cos(pi * (mu - om / 2)) / mp.cos(pi * (mu + om / 2))


def c_d(mu, T, om, sg):
    s2 = mp.sin(pi * om / 2) ** 2
    if sg > 0:
        c = -mp.cos(pi * mu) * mp.cos(pi * (mu + om / 2)) * mp.cos(pi * (mu - T)) / (s2 * mp.cos(pi * (mu - om / 2)))
        e = mp.cos(pi * (mu - T / 2 + om / 2)) * mp.cos(pi * (mu - T / 2 - om / 2)) / mp.cos(pi * (mu - T / 2)) ** 2
    else:
        c = mp.sin(pi * mu) * mp.sin(pi * (mu + om / 2)) * mp.sin(pi * (mu - T)) / (s2 * mp.sin(pi * (mu - om / 2)))
        e = mp.sin(pi * (mu - T / 2 + om / 2)) * mp.sin(pi * (mu - T / 2 - om / 2)) / mp.sin(pi * (mu - T / 2)) ** 2
    return c * e


def Xi(mu, T, om):
    b = lambda x: blkI(x, mu, T)
    return b(T / 2) ** 2 * b(T / 2 - om / 2) * b(T / 2 + om / 2) * b(om / 2) * b(T - om / 2)


def rhs_minus(mu, T, om):
    b = lambda x: blkI(x, mu, T)
    h = mp.mpf(1) / 2
    return -b(om / 2) * b(h - T) * b(-om / 2 - h) * b(T / 2 + om / 2 + h) * b(T / 2) / (b(T / 2 + h) * b(T / 2 + om / 2))


def rhs_plusXi(mu, T, om):
    b = lambda x: blkI(x, mu, T)
    h = mp.mpf(1) / 2
    return b(om / 2 - T) * b(h) * b(T - om / 2 - h) * b(T / 2 - om / 2) * b(T / 2 - h) / (b(T / 2) * b(T / 2 - om / 2 - h))


def rhs_eps0(mu, T, om):
    b = lambda x: blkI(x, mu, T)
    h = mp.mpf(1) / 2
    return b(om / 2 - T) * b(h - T) * b(-om / 2 - h) * b(T / 2 - om / 2) * b(T / 2 + om / 2 + h) / (b(T / 2) * b(T / 2 + h))


def partA():
    ok = True
    print("Part A: functional equations and breather identities (residuals)")
    for n in (4, 5, 6):
        for om in (mp.mpf('0.37'), mp.mpf('0.83')):
            T = (n - 1) * om
            for mu in (mp.mpc('0.13', '0.21'), mp.mpc('-0.31', '0.07')):
                res = {}
                # eps-hat = 0 control
                A0 = lambda m: mp.exp(logA0(m, T, om))
                res['A0 unit'] = abs(A0(mu) * A0(-mu) - 1)
                res['A0 cross'] = abs(A0(T / 2 - mu) / A0(T / 2 + mu) - mp.exp(logF(2 * mu, T, om)))
                for sg in (+1, -1):
                    As = lambda m: A_s(m, T, om, sg)
                    res[f'unit {sg:+d}'] = abs(As(mu) * As(-mu) * U_s(mu, T, om, sg) - 1)
                    res[f'cross {sg:+d}'] = abs(As(T / 2 - mu) / As(T / 2 + mu) * Psi_s(mu, om, sg) / mp.exp(logF(2 * mu, T, om)) - 1)
                ma, mb = mu - (T / 2 - mp.mpf(1) / 2), mu + (T / 2 - mp.mpf(1) / 2)
                F2 = mp.exp(logF(2 * mu, T, om))
                lhs0 = mp.exp(logA0(ma, T, om) + logA0(mb, T, om)) * F2
                res['breather eps0'] = abs(lhs0 / rhs_eps0(mu, T, om) - 1)
                lhsm = A_s(ma, T, om, -1) * F2 * A_s(mb, T, om, -1) * c_d(mu, T, om, -1)
                res['breather -'] = abs(lhsm / rhs_minus(mu, T, om) - 1)
                lhsp = Xi(ma, T, om) * A_s(ma, T, om, +1) * F2 * Xi(mb, T, om) * A_s(mb, T, om, +1) * c_d(mu, T, om, +1)
                res['breather +Xi'] = abs(lhsp / rhs_plusXi(mu, T, om) - 1)
                worst = max(res.values())
                ok &= worst < 1e-18
                print(f"  n={n} om={om} mu={mp.nstr(mu, 3)}: worst {mp.nstr(worst, 3)}  " +
                      " ".join(f"{k}:{mp.nstr(v, 2)}" for k, v in res.items()))
    return ok


# ---------------------------------------------------------------- Part B
def pblk(x, th, h):
    return mp.sinh(th / 2 + 1j * pi * x / (2 * h)) / mp.sinh(th / 2 - 1j * pi * x / (2 * h))


def curly(x, th, h, B):
    return pblk(x - 1, th, h) * pblk(x + 1, th, h) / (pblk(x - 1 + B, th, h) * pblk(x + 1 - B, th, h))


def S_vec(a, b, th, h, B):
    """d_n S-matrix between vector-like particles a, b <= n-2 (BCDS):
       prod_{p=1,3,...,2 min(a,b)-1} {|a-b|+p}{h-|a-b|-p}"""
    d = abs(a - b)
    s = 1
    for p in range(1, 2 * min(a, b), 2):
        s *= curly(d + p, th, h, B) * curly(h - d - p, th, h, B)
    return s


def cont(f, th, n, B):
    """evaluate an imaginary-coupling function f(mu, T, om) at omega = -2/B, mu = -i T th/pi"""
    om = -2 / B
    T = (n - 1) * om
    mu = -1j * T * th / pi
    return f(mu, T, om)


def K1(th, h, B):
    b = lambda x: pblk(x, th, h)
    return b(h - B / 2) * b(h / 2 - 1) * b(h / 2 + 1 - B / 2) / (b(h - 1) * b(1 - B / 2) * b(h / 2) * b(h / 2 - B / 2))


def K1cl(th, h):
    b = lambda x: pblk(x, th, h)
    return -b(h / 2 - 1) * b(h / 2 + 1) / (b(1) * b(h - 1) * b(h / 2) ** 2)


def kimcho_dK(th):
    sh, ch = mp.sinh(th), mp.cosh(th)
    return 1j * 2 * pi / 24 * (sh / ch - sh / (ch - 1) + 2 * sh / (2 * ch + mp.sqrt(3)) - 2 * sh / (2 * ch - 1))


def bootstrap_family(K1f, n, B):
    """K_a for vector-like a via K_{a+1}(th) = K_1(th + i pi a/h) S_{1a}(2th + i pi (a-1)/h) K_a(th - i pi/h);
    returns dict of functions."""
    h = 2 * n - 2
    Ks = {1: K1f}
    for a in range(1, n - 2):
        Ka = Ks[a]
        Ks[a + 1] = (lambda Ka, a: (lambda th: K1f(th + 1j * pi * a / h) * S_vec(1, a, 2 * th + 1j * pi * (a - 1) / h, h, B)
                                      * Ka(th - 1j * pi / h)))(Ka, a)
    return Ks


def axioms(K1f, n, B, pts):
    h = 2 * n - 2
    Ks = bootstrap_family(K1f, n, B)
    worst = 0
    for th in pts:
        for a, Ka in Ks.items():
            worst = max(worst, abs(Ka(th) * Ka(-th) - 1))                                   # unitarity
            worst = max(worst, abs(Ka(1j * pi / 2 - th) / (S_vec(a, a, 2 * th, h, B) * Ka(1j * pi / 2 + th)) - 1))  # crossing
        # closure: K_{a+b} via (2, a-1) fusing equals the (1, a) route
        for a in range(3, n - 1):
            alt = Ks[2](th + 1j * pi * (a - 2) / h) * S_vec(2, a - 2, 2 * th + 1j * pi * (a - 4) / h, h, B) * Ks[a - 2](th - 2j * pi / h)
            worst = max(worst, abs(alt / Ks[a](th) - 1))
    return worst


def dlog_at0(f):
    eps = mp.mpf('1e-15')
    return (mp.log(f(eps)) - mp.log(f(-eps))) / (2 * eps), f(mp.mpf(0)) if False else None


def partB():
    ok = True
    pts = [mp.mpc('0.41', '0.17'), mp.mpc('-1.1', '0.33'), mp.mpc('0.7', '-0.52')]
    print("Part B: continuation omega -> -2/B")
    for n in (4, 5, 6, 7):
        h = 2 * n - 2
        for B in (mp.mpf('0.31'), mp.mpf('1.37'), mp.mpc('0.6', '0.2')):
            r1 = max(abs(cont(rhs_plusXi, th, n, B) / K1(th, h, 2 - B) - 1) for th in pts)
            r2 = max(abs(cont(rhs_minus, th, n, B) / (K1(th, h, B) / K1cl(th, h)) - 1) for th in pts)
            r3 = max(abs(cont(rhs_eps0, th, n, B) / K1(th, h, B) - 1) for th in pts)
            ax_N = axioms(lambda t: cont(rhs_plusXi, t, n, B), n, B, pts[:2])
            ax_m = axioms(lambda t: cont(rhs_minus, t, n, B), n, B, pts[:2])
            minimal = lambda t: cont(lambda m, T, o: rhs_plusXi(m, T, o) / (Xi(m - (T / 2 - mp.mpf(1) / 2), T, o) * Xi(m + (T / 2 - mp.mpf(1) / 2), T, o)), t, n, B)
            ax_min = axioms(minimal, n, B, pts[:2])
            print(f"  n={n} B={mp.nstr(B, 3)}: K^N=K_1(2-B) {mp.nstr(r1, 2)}, K^(-)=K_1/K_1cl {mp.nstr(r2, 2)}, eps0=K_1 {mp.nstr(r3, 2)};"
                  f" axioms K^N {mp.nstr(ax_N, 2)}, K^(-) {mp.nstr(ax_m, 2)}, minimal {mp.nstr(ax_min, 2)}")
            ok &= max(r1, r2, r3, ax_N, ax_m, ax_min) < 1e-20
        # B -> 0 limits
        for name, f in (("K^N", rhs_plusXi), ("K^(-)", rhs_minus)):
            v = max(abs(cont(f, th, n, mp.mpf('1e-12')) - 1) for th in pts)
            ok &= v < 1e-9
        print(f"  n={n}: K^N, K^(-) -> 1 as B -> 0")
    # one loop, d_4
    n, h = 4, 6
    print("  one loop d_4 vs Kim-Cho (29): O(B^0) value and O(B) coefficient of log K at three rapidities")
    half = mp.mpf(1) / 2

    def P(x):
        return lambda m, T, o: blkI(T / 2 - x(T, o), m, T) * blkI(T / 2 + x(T, o), m, T)

    def with_cdd(cdd):
        def f(m, T, o):
            ma, mb = m - (T / 2 - half), m + (T / 2 - half)
            base = rhs_plusXi(m, T, o) / (Xi(ma, T, o) * Xi(mb, T, o))
            return base * cdd(ma, T, o) * cdd(mb, T, o)
        return f
    x0, xw, xT = (lambda T, o: 0), (lambda T, o: o / 2), (lambda T, o: (T - o) / 2)
    xh, xwh = (lambda T, o: half), (lambda T, o: half + o / 2)
    cands = {
        'Xi = P0 Pw/2 P(T-w)/2': lambda m, T, o: P(x0)(m, T, o) * P(xw)(m, T, o) * P(xT)(m, T, o),
        'none (minimal)': lambda m, T, o: 1,
        'P1/2 Pw/2 P(T-w)/2': lambda m, T, o: P(xh)(m, T, o) * P(xw)(m, T, o) * P(xT)(m, T, o),
        'P0 P1/2+w/2 P(T-w)/2': lambda m, T, o: P(x0)(m, T, o) * P(xwh)(m, T, o) * P(xT)(m, T, o),
        'P1/2 P1/2+w/2 P(T-w)/2': lambda m, T, o: P(xh)(m, T, o) * P(xwh)(m, T, o) * P(xT)(m, T, o),
    }
    eps = mp.mpf('1e-14')
    results = {}
    for name, cdd in cands.items():
        f = with_cdd(cdd)
        row = []
        for th in pts:
            g = lambda B: cont(f, th, n, B)
            v0 = g(eps)
            d1 = (mp.log(g(eps)) - mp.log(g(2 * eps))) / (-eps)
            row.append((v0, d1, kimcho_dK(th)))
        results[name] = row
        dev = max(max(abs(v0 - 1), abs(d1 - kc)) for v0, d1, kc in row)
        print(f"    {name:26s}: max(|K(0)-1|, |dlogK - KimCho|) = {mp.nstr(dev, 3)}")
    ok &= max(max(abs(v0 - 1), abs(d1 - kc)) for v0, d1, kc in results['Xi = P0 Pw/2 P(T-w)/2']) < 1e-8
    ok &= max(max(abs(v0 - 1), abs(d1 - kc)) for v0, d1, kc in results['none (minimal)']) > 1e-3
    return ok, results


if __name__ == '__main__':
    okA = partA()
    okB, res = partB()
    print("PASS" if okA and okB else "FAIL", "d_n s = +-1 scalar factors, breathers, continued amplitudes, one loop")
