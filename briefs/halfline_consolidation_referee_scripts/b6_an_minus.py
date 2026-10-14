"""Referee check of the a_n s = -1 scalar factor and amplitudes (ch. 33 prp-an-nondiagonal-scalar,
eq-an-nondiagonal-scalar-minus, eq-an-minus-breather; ch. 35 "The third matrix"), from the formulas
in the text only.
 1. A^(-) satisfies A(mu)A(-mu)U_- = 1 and A(T/2-mu)/A(T/2+mu) = F_11(2mu)/Psi_-(mu), n = 2, 3, 4;
 2. A^(-)(mu_a) F_11(2mu) A^(-)(mu_b) c_-(mu) = -(w/2)_I (1/2-T)_I (-w/2-1/2)_I, n = 1..4;
 3. continuation of that closed form = -(1)(h-B/2)(-1+B/2) = -(1)(h-1) K_1;
 4. K^(-)_a = -(a)(h-a) K_a: unitarity, crossing-unitarity, bootstrap at every fusing (n <= 7);
 5. duals -(a)(h-a) K^N_a = CDRS (5.16) K^0_a = prod_k <k-1/2>/<h-k+1/2>, <x> = (x+1/2)/(x-1/2+B/2);
    n = 2: CDRS (5.7) -(1)(2+B/2)/(B/2).
"""
import mpmath as mp

mp.mp.dps = 30
pi = mp.pi
half = mp.mpf(1) / 2


def lg(x, c):
    return mp.loggamma(x + c) - mp.loggamma(-x + c)


def blkI(x, mu, T):
    return mp.sin(pi * (mu + x) / (2 * T)) / mp.sin(pi * (mu - x) / (2 * T))


def loga1(mu, n, om):
    h = n + 1; T = h * om / 2
    m2 = 2 * mu
    def term(k):
        kk = k + 1; s = 4 * kk * T      # G_{2k}: + 2k h omega = 4 k T
        return (lg(m2, s - 3 * h * om / 2 + 1) + lg(m2, s - h * om / 2)
                - lg(m2, s - 3 * n * om / 2 - om / 2 + 1) - lg(m2, s - n * om / 2 - 3 * om / 2))
    return mp.nsum(term, [0, mp.inf])


def logGprod(mu, n, om, d):
    h = n + 1; T = h * om / 2
    def term(k):
        j = k + 1; s = 2 * j * T
        return (lg(mu, s + d + om / 2 - 3 * T / 2) + lg(mu, s + 1 - d - om / 2 - T / 2)
                - lg(mu, s + d + om / 2 - T / 2) - lg(mu, s + 1 - d - om / 2 - 3 * T / 2))
    return mp.nsum(term, [0, mp.inf])


def Aminus(mu, n, om):
    T = (n + 1) * om / 2
    pref = -1j * blkI(-T / 2, mu, T) * mp.sin(pi * om / 2) / mp.cos(pi * (mu - T / 2))
    return pref * mp.exp(logGprod(mu, n, om, half) + loga1(mu, n, om))


def F11(mu, n, om):
    h = n + 1; T = h * om / 2
    def term(k):
        j = k + 1; s = 2 * j * T
        return lg(mu, s - (h - 1) * om + 1) + lg(mu, s - om) - lg(mu, s - 2 * T + 1) - lg(mu, s)
    return -mp.exp(mp.nsum(term, [0, mp.inf]))


def check12():
    ok = True
    for n in (1, 2, 3, 4):
        for om in (mp.mpf('0.43'), mp.mpf('0.91')):
            T = (n + 1) * om / 2
            for mu in (mp.mpc('0.12', '0.23'), mp.mpc('-0.3', '0.11')):
                r = {}
                if n >= 2:
                    U = (-mp.cos(2 * pi * mu) - mp.cos(pi * T)) / (2 * mp.sin(pi * om / 2) ** 2)
                    Psi = mp.cos(pi * (mu - om / 2)) / mp.cos(pi * (mu + om / 2))
                    r['unit'] = abs(Aminus(mu, n, om) * Aminus(-mu, n, om) * U - 1)
                    r['cross'] = abs(Aminus(T / 2 - mu, n, om) / Aminus(T / 2 + mu, n, om) * Psi / F11(2 * mu, n, om) - 1)
                ma, mb = mu - (T / 2 - half), mu + (T / 2 - half)
                c = mp.sin(pi * mu) * mp.sin(pi * (mu + om / 2)) * mp.sin(pi * (mu - T)) / (mp.sin(pi * om / 2) ** 2 * mp.sin(pi * (mu - om / 2)))
                lhs = Aminus(ma, n, om) * F11(2 * mu, n, om) * Aminus(mb, n, om) * c
                rhs = -blkI(om / 2, mu, T) * blkI(half - T, mu, T) * blkI(-om / 2 - half, mu, T)
                r['breather'] = abs(lhs / rhs - 1)
                print(f"  n={n} om={om} mu={mp.nstr(mu, 3)}: " + ", ".join(f"{k} {mp.nstr(v, 2)}" for k, v in r.items()))
                ok &= max(r.values()) < 1e-20
    return ok


def pb(x, th, h):
    return mp.sinh(th / 2 + 1j * pi * x / (2 * h)) / mp.sinh(th / 2 - 1j * pi * x / (2 * h))


def curly(x, th, h, B):
    return pb(x - 1, th, h) * pb(x + 1, th, h) / (pb(x - 1 + B, th, h) * pb(x + 1 - B, th, h))


def San(a, b, th, h, B):
    r = 1
    for p in range(abs(a - b) + 1, a + b, 2):
        r *= curly(p, th, h, B)
    return r


def Ku(a, th, h, B):
    r = 1
    for c in range(1, a + 1):
        r *= pb(c - 1, th, h) * pb(h + 1 - c - B / 2, th, h) / (pb(h - c, th, h) * pb(c - B / 2, th, h))
    return r


def Km(a, th, h, B):
    return -pb(a, th, h) * pb(h - a, th, h) * Ku(a, th, h, B)


def angle(h, x):
    return x


def check345():
    ok = True
    pts = [mp.mpc('0.41', '0.17'), mp.mpc('-1.1', '0.33')]
    for n in range(2, 8):
        h = n + 1
        for B in (mp.mpf('0.37'), mp.mpc('1.2', '0.3')):
            # 3: continuation of the closed form
            def cont(th):
                om = -2 / B; T = h * om / 2; mu = -1j * T * th / pi
                return -blkI(om / 2, mu, T) * blkI(half - T, mu, T) * blkI(-om / 2 - half, mu, T)
            r3 = max(abs(cont(t) / Km(1, t, h, B) - 1) for t in pts)
            r3b = max(abs(cont(t) / (-pb(1, t, h) * pb(h - B / 2, t, h) * pb(-1 + B / 2, t, h)) - 1) for t in pts)
            # 4: axioms for K^(-)_a = -(a)(h-a) K_a
            r4 = 0
            for t in pts:
                for a in range(1, n + 1):
                    r4 = max(r4, abs(Km(a, t, h, B) * Km(a, -t, h, B) - 1))
                    r4 = max(r4, abs(Km(a, 1j * pi / 2 - t, h, B) / (San(a, h - a, 2 * t, h, B) * Km(h - a, 1j * pi / 2 + t, h, B)) - 1))
                    r4 = max(r4, abs(Km(a, t, h, B) / Km(h - a, t, h, B) - 1))
                # bootstrap K_{a+b}(th) = K_a(th + i pi b/h) S_ab(2th + i pi (b-a)/h) K_b(th - i pi a/h), a+b <= n
                for a in range(1, n + 1):
                    for b in range(1, n + 1 - a):
                        lhs = Km(a + b, t, h, B)
                        rhs = Km(a, t + 1j * pi * b / h, h, B) * San(a, b, 2 * t + 1j * pi * (b - a) / h, h, B) * Km(b, t - 1j * pi * a / h, h, B)
                        r4 = max(r4, abs(lhs / rhs - 1))
            # 5: duals vs CDRS (5.16)
            r5 = 0
            for t in pts:
                for a in range(1, n + 1):
                    dual = -pb(a, t, h) * pb(h - a, t, h) * Ku(a, t, h, 2 - B)
                    ang = lambda x: pb(x + half, t, h) / pb(x - half + B / 2, t, h)
                    cdrs = 1
                    for k in range(1, a + 1):
                        cdrs *= ang(k - half) / ang(h - k + half)
                    r5 = max(r5, abs(dual / cdrs - 1))
            print(f"  n={n} B={mp.nstr(B, 3)}: continuation {mp.nstr(max(r3, r3b), 2)}, axioms {mp.nstr(r4, 2)}, dual = CDRS (5.16) {mp.nstr(r5, 2)}")
            ok &= max(r3, r3b, r4, r5) < 1e-20
            # B -> 0 limit
            ok &= max(abs(Km(a, pts[0], h, mp.mpf('1e-15')) - 1) for a in range(1, n + 1)) < 1e-10
    # n = 2: CDRS (5.7)
    B = mp.mpf('0.37'); h = 3
    r = max(abs(-pb(1, t, h) * pb(2, t, h) * Ku(1, t, h, 2 - B) / (-pb(1, t, h) * pb(2 + B / 2, t, h) / pb(B / 2, t, h)) - 1) for t in pts)
    print(f"  n=2: dual of K^(-)_1 = CDRS (5.7) -(1)(2+B/2)/(B/2): {mp.nstr(r, 2)}")
    ok &= r < 1e-20
    return ok


if __name__ == '__main__':
    ok1 = check12()
    ok2 = check345()
    print("PASS" if ok1 and ok2 else "FAIL", "a_n s = -1: scalar factor, breather, continued amplitudes, axioms, duals = CDRS (5.16)")
