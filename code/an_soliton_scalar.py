"""Scalar factor of the a_{h-1}^(1) soliton S-matrix: chapter 21 (sec-scalar-factor,
sec-an-example), h = 2, 3, 4, at the generic coupling omega = 1 + sqrt(2)/2.

Variables of chapter 21: t = T theta/(i pi), x = e^{2 pi i t}, q = -e^{-i pi omega},
2T = h omega.  Species-1 solitons, R_11 = P_sym + rho P_antisym.

- tau(t) = (x - 1)/(q x - 1/q) = -sin(pi t)/sin(pi(t - omega)) is the transmission entry
  of R_11 (charges e_i, e_k with i != k), divided by the identical-charge entry. It is the
  same for every h, because Jimbo's sl_h R-matrix restricted to e_i x e_k, e_k x e_i is the
  six-vertex one. It is ch. 7's b(z) = q(z-1)/(z-q^2) with q -> 1/q (coproduct of Part VI),
  and its pole sits at t = omega.
- Book: F = c f of sec-scalar-factor with the single zero a_1 = omega, c = sin pi(t - omega),
  the regularisation r = 0.371 of the text, and the constant fixed by F(t)F(-t) = 1.
- Hollowood, hep-th/9203076, eq. (4.5): S_{ik->ki}, the transmission amplitude S_T. His n is
  h, his lambda is omega (his pole tower (4.7) is t = omega - p), and i theta/2pi = -t/(2T).
  His identical-charge amplitude is Phi_11 = S_T/tau.

The test that separates them: crossing particle 1 of the transmission process e_i e_k (i != k)
gives antisoliton(i) + soliton(k), whose charge sector in V_bar x V is one-dimensional. So
unitarity there is a scalar condition, S_T(T - t) S_T(T + t) = 1. The c f form satisfies
f(T - t) = f(t) instead, and then S_T(T-t)S_T(T+t) = (sin^2 pi T - sin^2 pi t)/(sin^2 pi omega
- sin^2 pi t), which is 1 only for T = +-omega mod 1, i.e. for h = 2 at generic omega.
Semiclassically, ch. 10 (sec-multi-solitons) gives -i dln S/dtheta -> (2h/beta^2) ln A_11.
"""
import numpy as np
import mpmath as mp
from scipy.special import loggamma, psi

OM = 1 + np.sqrt(2) / 2           # generic: q is not a root of unity
R_OFF = 0.371                      # regularisation offset of sec-scalar-factor
TS = np.array([0.31 + 0.47j, -0.2 + 0.9j, 0.77 - 0.35j])


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


def rich(fun, n):
    """Richardson extrapolation of a truncated sum with an O(1/n) tail."""
    a, b = fun(n), fun(2 * n)
    return 2 * b - a


def tau(t, om):
    return -np.sin(np.pi * t) / np.sin(np.pi * (t - om))


def dlntau(t, om):
    return np.pi / np.tan(np.pi * t) - np.pi / np.tan(np.pi * (t - om))


# ---- Hollowood (4.5), z = i theta/2pi = -t/(2T), n lambda = 2T, 1/n = 1/h ----
def _holl_args(t, h, om, J):
    T = h * om / 2
    z = (-np.atleast_1d(t) / (2 * T))[None, :]
    j = np.arange(1, J + 1)[:, None]
    u, v, r = j / (2 * T), (j - 1) / (2 * T), 1 / h
    plus = [1 + z + u, 1 + z + v, -z + r + u, -z - r + v]
    minus = [-z + u, -z + v, 1 + z + r + u, 1 + z - r + v]
    return plus, minus


def ln_ST_holl(t, h, om, J=20000):
    def part(J):
        plus, minus = _holl_args(t, h, om, J)
        return (sum(loggamma(a) for a in plus) - sum(loggamma(a) for a in minus)).sum(axis=0)
    return rich(part, J)


def dln_ST_holl(t, h, om, J=20000):      # d/dt
    T = h * om / 2

    def part(J):
        plus, minus = _holl_args(t, h, om, J)
        s = [1, 1, -1, -1]                 # d/dz of the argument
        num = sum(si * psi(a) for si, a in zip(s, plus))
        den = sum(si * psi(a) for si, a in zip([-1, -1, 1, 1], minus))
        return (num - den).sum(axis=0) * (-1 / (2 * T))
    return rich(part, J)


# ---- book: F = c f, f1(t) = Gamma(t - a) Gamma(1 + t + a)/pi, a = omega ----
def _L1(s, a):
    s = np.asarray(s, dtype=complex)
    return loggamma(s - a) + loggamma(1 + s + a) - np.log(np.pi)


def _dL1(s, a):
    return psi(s - a) + psi(1 + s + a)


def ln_f(t, h, om, M=20000):
    """ln f up to a t-independent constant (the regularised product of the text)."""
    T = h * om / 2

    def part(M):
        m = np.arange(1, M + 1)[:, None]
        tt = np.atleast_1d(t)[None, :]
        out = 0
        for sgn, sh, A in ((1, 1, 2 * T * (m - 1)), (1, -1, 2 * T * (m - .5)),
                           (-1, 1, 2 * T * (m - .5)), (-1, -1, 2 * T * m)):
            out = out + sgn * (_L1(sh * tt + A, om) - _L1(A + R_OFF, om))
        return out.sum(axis=0)
    return rich(part, M)


def dln_F(t, h, om, M=20000):            # d/dt ln(c f)
    T = h * om / 2

    def part(M):
        m = np.arange(1, M + 1)[:, None]
        tt = np.atleast_1d(t)[None, :]
        s = (_dL1(tt + 2 * T * (m - 1), om) - _dL1(-tt + 2 * T * (m - .5), om)
             - _dL1(tt + 2 * T * (m - .5), om) + _dL1(-tt + 2 * T * m, om))
        return s.sum(axis=0)
    return rich(part, M) + np.pi / np.tan(np.pi * (np.atleast_1d(t) - om))


def ln_F(t, h, om):
    """ln(c f) with the constant fixed by F(t)F(-t) = 1."""
    c = lambda s: np.log(np.sin(np.pi * (s - om)))
    t0 = np.array([0.23 + 0.31j])
    k2 = -(c(t0) + c(-t0) + ln_f(t0, h, om) + ln_f(-t0, h, om))[0]
    return c(np.atleast_1d(t)) + ln_f(t, h, om) + k2 / 2


# 0. Jimbo's sl_3 R-matrix (x_book = x^2): Yang-Baxter, and its transmission ratio is tau
def jimbo(x, q, N=3):
    E = lambda i, j: np.outer(np.eye(N)[i], np.eye(N)[j])
    R = np.zeros((N * N, N * N), complex)
    for i in range(N):
        for j in range(N):
            if i == j:
                R += (q * x - 1 / (q * x)) * np.kron(E(i, i), E(i, i))
            else:
                R += (x - 1 / x) * np.kron(E(i, i), E(j, j))
                R += (q - 1 / q) * (x if i < j else 1 / x) * np.kron(E(i, j), E(j, i))
    return R


qg = -np.exp(-1j * np.pi * OM)
I3 = np.eye(3)
P23 = np.zeros((27, 27))
for a in range(3):
    for b in range(3):
        for c in range(3):
            P23[a * 9 + c * 3 + b, a * 9 + b * 3 + c] = 1
x1, x2, x3 = 1.3 + 0.4j, 0.7 - 0.2j, 1.1 + 0.9j
R12 = lambda u: np.kron(jimbo(u, qg), I3)
R23 = lambda u: np.kron(I3, jimbo(u, qg))
R13 = lambda u: P23 @ R12(u) @ P23
ybe = np.abs(R12(x1 / x2) @ R13(x1 / x3) @ R23(x2 / x3) - R23(x2 / x3) @ R13(x1 / x3) @ R12(x1 / x2)).max()
xs = np.sqrt(np.exp(2j * np.pi * TS))
Rd = [jimbo(u, qg) for u in xs]
ratio = np.array([R[0 * 3 + 1, 0 * 3 + 1] / R[0, 0] for R in Rd])
report(f'Jimbo sl_3 R-matrix: Yang-Baxter [{ybe:.1e}]; transmission/identical entry = tau = '
       f'-sin(pi t)/sin(pi(t - omega)), pole at t = omega [{np.abs(ratio - tau(TS, OM)).max():.1e}]',
       ybe < 1e-12 and np.abs(ratio - tau(TS, OM)).max() < 1e-12)

# 1. calibration at h = 2 against Zamolodchikov (eq-sg-s0, eq-sg-smatrix), lambda = T = omega
xi = np.pi / OM


def lnS0(t):
    th = 1j * np.pi * t / OM
    g = lambda k: mp.sinh((mp.pi - xi) * k / 2) * mp.sin(th * k) / (k * mp.sinh(xi * k / 2) * mp.cosh(mp.pi * k / 2))
    return complex(mp.log(-1) - 1j * mp.quad(g, [0, 1, 5, 20, mp.inf]))


S0 = np.exp([lnS0(t) for t in TS])
eps = 1e-5
dS0 = np.array([(lnS0(t + eps) - lnS0(t - eps)) / (2 * eps) for t in TS])
STz = S0 * tau(TS, OM)
d1 = np.abs(dln_F(TS, 2, OM) - dS0).max()
d2 = np.abs(np.exp(ln_ST_holl(TS, 2, OM)) / STz - 1).max()
report(f'h = 2: c f has the log-derivative of S_0 [{d1:.1e}], and Hollowood (4.5) equals S_T = '
       f'S_0 sinh(lambda theta)/sinh(lambda(i pi - theta)) [{d2:.1e}]', d1 < 1e-6 and d2 < 1e-6)

for h in (2, 3, 4):
    T = h * OM / 2
    # 2. the book's c f: unitarity and f(T - t) = f(t)
    u = ln_F(TS, h, OM) + ln_F(-TS, h, OM)
    sym = ln_f(T - TS, h, OM) - ln_f(TS, h, OM)
    e2 = max(np.abs(np.exp(u) - 1).max(), np.abs(np.exp(sym) - 1).max())
    report(f'h = {h}: c f satisfies F(t)F(-t) = 1 and f(T - t) = f(t) [{e2:.1e}]', e2 < 1e-5)
    # 3. Hollowood: unitarity of Phi_11 = S_T/tau, and of the soliton-antisoliton transmission
    phi = np.exp(ln_ST_holl(TS, h, OM) + ln_ST_holl(-TS, h, OM)) / (tau(TS, OM) * tau(-TS, OM))
    cu = np.exp(ln_ST_holl(T - TS, h, OM) + ln_ST_holl(T + TS, h, OM))
    e3 = max(np.abs(phi - 1).max(), np.abs(cu - 1).max())
    report(f'h = {h}: Hollowood (4.5): Phi_11(t)Phi_11(-t) = 1 and S_T(T - t)S_T(T + t) = 1 [{e3:.1e}]', e3 < 1e-6)
    # 4. the book's c f in the soliton-antisoliton channel
    cf = np.exp(ln_F(T - TS, h, OM) + ln_F(T + TS, h, OM)) * tau(T - TS, OM) * tau(T + TS, OM)
    pred = (np.sin(np.pi * T) ** 2 - np.sin(np.pi * TS) ** 2) / (np.sin(np.pi * OM) ** 2 - np.sin(np.pi * TS) ** 2)
    dev = np.abs(cf - 1).max()
    ok = np.abs(cf / pred - 1).max() < 1e-5 and (dev < 1e-5 if h == 2 else dev > 1e-2)
    report(f'h = {h}: c f gives S_T(T - t)S_T(T + t) = (sin^2 pi T - sin^2 pi t)/(sin^2 pi omega - sin^2 pi t) '
           f'[{np.abs(cf / pred - 1).max():.1e}], so |S_T(T - t)S_T(T + t) - 1| = {dev:.1e}: '
           + ('unitary' if h == 2 else 'not unitary for h >= 3'), ok)
    # 5. Hollowood's bound-state pole: simple pole of S_T at t = omega, Phi_11 regular there
    if h > 2:
        tp = OM + np.array([1e-4, 1e-5])
        st = np.exp(ln_ST_holl(tp, h, OM))
        ph = st / tau(tp, OM)
        ok5 = abs(st[1] * 1e-5 / (st[0] * 1e-4) - 1) < 1e-3 and abs(ph[1] / ph[0] - 1) < 1e-3
        report(f'h = {h}: Hollowood S_T has a simple pole at t = omega (theta = 2 pi i/h), the top of '
               f'the tower t = omega - p, and Phi_11 = S_T/tau is regular there', ok5)

# 6. semiclassics: beta^2/(2h) (-i dln S_T/dtheta) -> ln A_11 (sec-multi-solitons), 4pi/beta^2 = omega + 1
th = np.array([0.5, 1.0, 2.0])
for h in (2, 3, 4):
    lnA = np.log((np.cosh(th) - 1) / (np.cosh(th) - np.cos(2 * np.pi / h)))
    errH, valF = [], []
    for om in (20 * np.sqrt(2), 40 * np.sqrt(2)):
        T = h * om / 2
        t = T * th / (1j * np.pi)
        scale = (T / (1j * np.pi)) * (-1j) * 2 * np.pi / (h * (om + 1))
        vH = (scale * dln_ST_holl(t, h, om, J=int(100 * 2 * T))).real
        vF = (scale * (dln_F(t, h, om, M=5000) + dlntau(t, om))).real
        errH.append(np.abs(vH / lnA - 1).max())
        valF.append(vF)
    okH = errH[1] < 0.6 * errH[0] and errH[1] < 0.03
    relF = np.abs(valF[1] / lnA - 1)
    okF = relF.max() < 0.03 if h == 2 else relF.min() > 0.3
    report(f'h = {h}: semiclassical phase: Hollowood -> ln A_11 (relative error {errH[0]:.3f}, {errH[1]:.3f} at '
           f'omega = 28.3, 56.6); c f gives {np.round(valF[1], 3)} against ln A_11 = {np.round(lnA, 3)}', okH and okF)
