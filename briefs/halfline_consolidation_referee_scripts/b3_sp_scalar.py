"""Referee check of prp-sp-scalar, eq-sp-scalar-eqs (ch. 33) and eq-sp-breather (ch. 35).

Matrices: the book's library boundary_common (a_n vector and conjugate representations, R-check
solver), as allowed by the brief; the de Vega - Gonzalez-Ruiz K's are built from the formulas of
sec-preserving-conjugating. Scalar functions: mpmath, written from the formulas in the text.
Checks
 1. K_V solves the braid reflection equation (sanity of the K convention);
 2. R-check_{V Vbar} is 1 on every transmission v_i (x) vbar_j -> vbar_j (x) v_i, i != j;
 3. crossing matrix parts: k_B(th) = Phi_1 R_{VVbar}(2th) k_A(-th), k_A(th) = Phi_2 R_{VbarV}(2th) k_B(-th),
    Phi_1 = -(-q)^{1-l} s(mu), Phi_2 = -(-q)^{l-1} s(mu)  (independent of eta);
 4. eq-sp-scalar solves eq-sp-scalar-eqs at complex mu (n,l) = (2,1), (3,2), (3,0);
 5. breather: the singlet in V(th_a) (x) Vbar(th_b) is mapped by (1 x a K_V(th_a)) F(2mu) R (1 x b~ K_Vbar(th_b))
    to K_{1->n}(mu) times the singlet in Vbar(-th_b) (x) V(-th_a) (normalizations as stated in ch. 35).
"""
import os
import sys
import numpy as np
import mpmath as mp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code'))
import boundary_common as bc   # noqa: E402

mp.mp.dps = 25
pi = np.pi


def KV(n, om, th, eta, l, bar=False):
    h = n + 1; T = h * om / 2
    d = []
    for j in range(h):
        up = j < l
        a = 2 * j if up else 2 * j - h
        sg = (+1 if up else -1)
        if not bar:
            d.append(np.exp(a * om * th) * np.sinh(eta + sg * T * th))
        else:
            etab = eta + 1j * pi * om * (l - h / 2)
            d.append(np.exp(-a * om * th) * np.sinh(etab - sg * T * th))
    return np.diag(d)


def rep(n, q, th, om, conj=False):
    return bc.an_vector(n, q, np.exp(om * th), conj=conj)


def matrix_checks(n, l, om, eta):
    q = bc.q_phys(om); h = n + 1; T = h * om / 2
    out = {}
    # 1. reflection equation (braid form) with R on V (x) V
    def Rf(z):
        R, _ = bc.rcheck(bc.an_vector(n, q, z), bc.an_vector(n, q, 1.0))
        return R
    Kf = lambda y: KV(n, om, np.log(y) / om, eta, l)
    res, _ = bc.re_preserving(Rf, Kf, np.exp(om * (0.3 + 0.2j)), np.exp(om * (-0.5 + 0.1j)))
    out['RE'] = res
    # 2. transmissions of R_{V Vbar}
    th = 0.27 + 0.19j
    R, _ = bc.rcheck(rep(n, q, th, om), rep(n, q, -th, om, conj=True))
    N = h
    tr = [R[j * N + i, i * N + j] for i in range(N) for j in range(N) if i != j]
    out['transmissions'] = max(abs(t - 1) for t in tr)
    # 3. crossing matrix parts
    mu = -1j * T * th / pi
    s = np.sin(pi * (2 * mu - T)) / np.sin(pi * (2 * mu + T))
    def kA(t):
        K = KV(n, om, 1j * pi / 2 - t, eta, l, bar=True)
        v = np.zeros(N * N, complex)
        for i in range(N):
            v[i * N + i] = K[i, i]
        return v
    def kB(t):
        K = KV(n, om, 1j * pi / 2 - t, eta, l)
        v = np.zeros(N * N, complex)
        for i in range(N):
            v[i * N + i] = K[i, i]
        return v
    R1, _ = bc.rcheck(rep(n, q, th, om), rep(n, q, -th, om, conj=True))       # V(th) x Vbar(-th) -> Vbar(-th) x V(th)
    R2, _ = bc.rcheck(rep(n, q, th, om, conj=True), rep(n, q, -th, om))       # Vbar(th) x V(-th) -> V(-th) x Vbar(th)
    lhs1, rhs1 = kB(th), R1 @ kA(-th)
    lhs2, rhs2 = kA(th), R2 @ kB(-th)
    i0 = np.argmax(abs(rhs1)); P1 = lhs1[i0] / rhs1[i0]
    i0 = np.argmax(abs(rhs2)); P2 = lhs2[i0] / rhs2[i0]
    out['cross prop'] = max(np.abs(lhs1 - P1 * rhs1).max() / np.abs(lhs1).max(), np.abs(lhs2 - P2 * rhs2).max() / np.abs(lhs2).max())
    out['Phi1'] = abs(P1 / (-(-q) ** (1 - l) * s) - 1)
    out['Phi2'] = abs(P2 / (-(-q) ** (l - 1) * s) - 1)
    return out


# ------------------------------------------------------------------ scalar functions (mpmath)
def lg(x, c):
    return mp.loggamma(x + c) - mp.loggamma(-x + c)


def logcalG(x, c, cp, T, step=2):
    """log of calG_{c,c'}(x) with G_j(x,a) = Gamma(x + step*j*T + a)/Gamma(-x + step*j*T + a)"""
    def term(k):
        j = k + 1
        s = step * j * T
        return lg(x, s + c) + lg(x, s + 1 - c) - lg(x, s + cp) - lg(x, s + 1 - cp)
    return lg(x, c) - lg(x, cp) + mp.nsum(term, [0, mp.inf])


def blkI(x, mu, T):
    return mp.sin(mp.pi * (mu + x) / (2 * T)) / mp.sin(mp.pi * (mu - x) / (2 * T))


def F11(mu, n, om):
    h = n + 1; T = h * om / 2
    def term(k):
        j = k + 1; s = 2 * j * T
        return lg(mu, s - (h - 1) * om + 1) + lg(mu, s - om) - lg(mu, s - 2 * T + 1) - lg(mu, s)
    return -mp.exp(mp.nsum(term, [0, mp.inf]))


def ST(mu, n, om):
    return mp.sin(mp.pi * mu) / mp.sin(mp.pi * (om - mu)) * F11(mu, n, om)


def scalars(n, l, om, eta):
    h = n + 1; T = h * om / 2
    eps = eta / (1j * mp.pi)
    epsb = eps + om * l - T
    rho = lambda mu: blkI(om / 2, mu, T) * mp.exp(logcalG(2 * mu, 2 * T + om, 0, T, step=4))
    def a(mu):
        th = 1j * mp.pi * mu / T
        return -1j * rho(mu) * mp.exp((1 - l) * om * th) * mp.exp(logcalG(mu, eps + 2 * T, epsb + T, T)) / mp.sin(mp.pi * (eps + mu))
    def b(mu):
        th = 1j * mp.pi * mu / T
        return -1j * rho(mu) * mp.exp(-(1 - l) * om * th) * mp.exp(logcalG(mu, epsb + 2 * T, eps + T, T)) / mp.sin(mp.pi * (epsb + mu))
    return a, b, T, eps, epsb


def scalar_checks(n, l, om, eta, mu):
    a, b, T, eps, epsb = scalars(n, l, om, eta)
    q = -mp.exp(-1j * mp.pi * om)
    etab = eta + 1j * mp.pi * om * (l - mp.mpf(h_of(n)) / 2)
    U = mp.sinh(eta + 1j * mp.pi * mu) * mp.sinh(eta - 1j * mp.pi * mu)
    Ub = mp.sinh(etab + 1j * mp.pi * mu) * mp.sinh(etab - 1j * mp.pi * mu)
    s = mp.sin(mp.pi * (2 * mu - T)) / mp.sin(mp.pi * (2 * mu + T))
    F = ST(T - 2 * mu, n, om)
    r = {}
    r['unit a'] = abs(a(mu) * a(-mu) * U - 1)
    r['unit b'] = abs(b(mu) * b(-mu) * Ub - 1)
    r['cross 1'] = abs(-(-q) ** (1 - l) * s * a(T / 2 - mu) / (F * b(T / 2 + mu)) - 1)
    r['cross 2'] = abs(-(-q) ** (l - 1) * s * b(T / 2 - mu) / (F * a(T / 2 + mu)) - 1)
    return r


def h_of(n):
    return n + 1


def singlet(A, B):
    M = []
    for j in bc.nodes_of(A):
        M.append(bc.cop(A, B, 'e', j)); M.append(bc.cop(A, B, 'f', j))
        M.append(bc.cop(A, B, 'k', j) - np.eye(A['dim'] * B['dim']))
    _, s, Vh = np.linalg.svd(np.vstack(M))
    assert s[-1] / s[0] < 1e-10 and s[-2] / s[0] > 1e-6
    return Vh[-1].conj()


def breather_check(n, l, om, eta, th):
    q = bc.q_phys(om); h = n + 1; T = h * om / 2; N = h
    mu = -1j * T * th / pi
    tha = th - 1j * pi / 2 + 1j * pi / (2 * T); thb = th + 1j * pi / 2 - 1j * pi / (2 * T)
    mua, mub = mu - (T / 2 - 0.5), mu + (T / 2 - 0.5)
    A = rep(n, q, tha, om); Bb = rep(n, q, thb, om, conj=True)
    sv = singlet(A, Bb); sv = sv / sv[0]                         # v_1 (x) vbar_1  (hw of V, lowest of Vbar)
    R, _ = bc.rcheck(rep(n, q, tha, om), rep(n, q, -thb, om, conj=True))
    w = bc.one_k(KV(n, om, tha, eta, l), N) @ R @ bc.one_k(KV(n, om, thb, eta, l, bar=True), N) @ sv
    sv2 = singlet(rep(n, q, -thb, om, conj=True), rep(n, q, -tha, om))
    i2 = n * N + n                                               # vbar_{n+1} (hw of Vbar) (x) v_{n+1} (lowest of V)
    sv2 = sv2 / sv2[i2]
    c = w[i2]
    prop = np.abs(w - c * sv2).max() / np.abs(w).max()
    a, b, T_, eps, epsb = scalars(n, l, mp.mpf(om), mp.mpc(eta))
    mu_m, mua_m, mub_m = mp.mpc(mu), mp.mpc(mua), mp.mpc(mub)
    full = a(mua_m) * ST(T_ - 2 * mu_m, n, mp.mpf(om)) * b(mub_m) * mp.mpc(c)
    # eq-sp-breather
    d = (T_ - om) / 2
    hf = mp.mpf(1) / 2
    K0 = (blkI(T_ / 2, mu_m, T_) * blkI(T_ / 2 - hf, mu_m, T_) * blkI(d, mu_m, T_) * blkI(hf - d, mu_m, T_)
          * blkI(hf + d, mu_m, T_) * blkI(T_ - hf - d, mu_m, T_))
    Aa = eps + (T_ - 1) / 2; Ab = Aa + om * l - T_
    thm = 1j * mp.pi * mu_m / T_
    E = (mp.sinh(thm) - 1j * mp.sin(mp.pi * Aa / T_)) / (mp.sinh(thm) + 1j * mp.sin(mp.pi * Ab / T_))
    K1n = mp.exp(-2j * mp.pi * l / h) * K0 * E
    return prop, abs(full / K1n - 1)


if __name__ == '__main__':
    ok = True
    om = 0.41; eta = 0.3 + 0.7j
    for (n, l) in [(2, 1), (3, 2), (3, 0), (2, 0)]:
        r = matrix_checks(n, l, om, eta)
        print(f"n={n} l={l}: " + ", ".join(f"{k} {v:.1e}" for k, v in r.items()))
        ok &= max(r.values()) < 1e-9
    for (n, l) in [(2, 1), (3, 2), (3, 1)]:
        for mu in (mp.mpc('0.13', '0.21'), mp.mpc('-0.27', '0.09')):
            r = scalar_checks(n, l, mp.mpf(om), mp.mpc(eta), mu)
            print(f"scalar n={n} l={l} mu={mp.nstr(mu, 3)}: " + ", ".join(f"{k} {mp.nstr(v, 2)}" for k, v in r.items()))
            ok &= max(r.values()) < 1e-20
    for (n, l) in [(2, 1), (3, 2), (2, 0), (3, 1)]:
        for th in (0.31 + 0.17j, -0.2 + 0.4j):
            prop, res = breather_check(n, l, om, eta, th)
            print(f"breather n={n} l={l} th={th}: singlet proportionality {prop:.1e}, |result/K_(1->n) - 1| {mp.nstr(res, 3)}")
            ok &= prop < 1e-9 and res < 1e-9
    print("PASS" if ok else "FAIL", "soliton-preserving crossing functions, scalar factors and breather amplitude")
