"""Numerical fluctuation operator A = -d^2/dx^2 + M(x) about a static Hirota (multi-)soliton,
built from polynomial tau functions, plus exact Hirota fluctuation solutions evaluated on the line.

M(x) = sum_j n_j alpha_j alpha_j^T exp(i alpha_j . phi_s),  exp(i alpha_j.phi_s) = prod_{l~j} tau_l / tau_j^2
(beta = 1, m = 1, simply laced).  Folded theories are treated as invariant sectors of the parent.
"""
import numpy as np
import mpmath as mp
from lie import Algebra, delta_mp
from hirota import MpCtx, soliton_deg, fluct_deg


class Background:
    def __init__(self, A, t, mu, xi=None):
        """t: coefficient lists (mp or complex) of tau_j in E = exp(mu x + xi)"""
        self.A = A
        self.t = [np.array([complex(c) for c in tj]) for tj in t]
        self.tmp = t
        self.mu = float(mu)
        self.D = np.array([len(tj) - 1 for tj in t])
        if xi is None:
            xi = self.choose_xi()
        self.xi = xi

    def taus_scaled(self, x, xi=None):
        """returns tau_j(x) * exp(-D_j * max(mu x + Re xi, 0)) (scaled to avoid overflow), as array (r+1, len(x))"""
        xi = self.xi if xi is None else xi
        x = np.atleast_1d(x)
        L = self.mu * x + xi.real
        s = np.maximum(L, 0)
        out = []
        for j, tj in enumerate(self.t):
            acc = np.zeros(len(x), complex)
            for p, c in enumerate(tj):
                acc += c * np.exp(p * (L + 1j * xi.imag) - self.D[j] * s)
            out.append(acc)
        return np.array(out), s

    def choose_xi(self, re=0.0):
        """choose Im xi maximizing min_{x,j} |tau_j| / sum_p |t_p E^p| (scale-free distance from a zero)"""
        best, bt = -1, None
        xs = np.linspace(-12, 12, 2401) / self.mu
        L = self.mu * xs + re
        absT = []
        for j, tj in enumerate(self.t):
            absT.append(sum(abs(c) * np.exp(p * L - self.D[j] * np.maximum(L, 0)) for p, c in enumerate(tj)))
        absT = np.array(absT)
        for th in np.linspace(0, 2 * np.pi, 721)[:-1]:
            T, s = self.taus_scaled(xs, re + 1j * th)
            m = np.min(np.abs(T) / absT)
            if m > best:
                best, bt = m, th
        self.min_tau = best
        return re + 1j * bt

    def expo(self, x):
        """exp(i alpha_j . phi_s) for each j: array (r+1, len x)"""
        T, s = self.taus_scaled(x)
        A = self.A
        out = []
        for j in range(A.r + 1):
            pr = np.ones(T.shape[1], complex)
            for l in A.nbrs[j]:
                pr = pr * T[l]
            out.append(pr / T[j] ** 2)   # scale factors cancel because sum_{l~j} D_l = 2 D_j
        return np.array(out)

    def M(self, x):
        E = self.expo(x)
        A = self.A
        al = A.alpha
        x = np.atleast_1d(x)
        out = np.zeros((len(x), A.r, A.r), complex)
        for j in range(A.r + 1):
            out += A.n[j] * E[j][:, None, None] * np.outer(al[j], al[j])[None]
        return out

    def phi(self, x):
        """phi_s = i sum_j alpha_j ln tau_j (log of scaled tau plus D_j s); branch continuous along x"""
        T, s = self.taus_scaled(x)
        lg = np.log(T) + self.D[:, None] * s[None]
        lg = np.unwrap(lg.imag, axis=1) * 1j + lg.real
        return 1j * (self.A.alpha.T @ lg)


def single_soliton(A, a, dps=30, xi=None):
    ctx = MpCtx(dps)
    d, lam, _ = delta_mp(A, a, dps=dps)
    mu = mp.sqrt(lam)
    t, dr, res = soliton_deg(A, d, mu, ctx=ctx)
    assert dr < mp.mpf(10) ** (-dps + 8) and res < mp.mpf(10) ** (-dps + 8), (dr, res)
    return Background(A, t, mu, xi), t, mu


class Channels:
    """high-precision first-order data for all species"""
    def __init__(self, A, dps=30):
        self.d, self.lam = {}, {}
        for b in A.mass:
            d, lam, _ = delta_mp(A, b, dps=dps)
            self.d[b], self.lam[b] = d, lam


def exact_solution(A, t, mu, ch, b, z, xi, x, dps=30, normalize_plus=False):
    """Hirota fluctuation solution in channel b at z = ik/m_b on grid x: returns (u, u', X, s).
    u(x) = sum_j alpha_j f_j(x), f_j = e^{ikx} S_j(E)/T_j(E), E = exp(mu x + xi).
    Normalized so that u -> e^{ikx} u_b(raw) at -infinity, with raw u_b = sum_j alpha_j delta^b_j."""
    ctx = MpCtx(dps)
    s, dr, res = fluct_deg(A, t, mu, ch.d[b], ch.lam[b], z, ctx=ctx)
    D = [len(tj) - 1 for tj in t]
    w = np.array([complex(s[j][D[j]] / t[j][D[j]]) for j in range(A.r + 1)])
    d0 = np.array([complex(c) for c in ch.d[b]])
    X = complex(np.vdot(d0, w) / np.vdot(d0, d0))
    mb = float(mp.sqrt(ch.lam[b]))
    ik = complex(z) * mb
    muf = float(mu)
    x = np.atleast_1d(x)
    L = muf * x + xi.real
    sc = np.maximum(L, 0)
    f = np.zeros((A.r + 1, len(x)), complex)
    fp = np.zeros((A.r + 1, len(x)), complex)
    for j in range(A.r + 1):
        Tj = np.zeros(len(x), complex); Tpj = np.zeros(len(x), complex)
        Sj = np.zeros(len(x), complex); Spj = np.zeros(len(x), complex)
        for p, c in enumerate(t[j]):
            e = np.exp(p * (L + 1j * xi.imag) - D[j] * sc)
            Tj += complex(c) * e; Tpj += complex(c) * p * muf * e
        for q, c in enumerate(s[j]):
            e = np.exp(q * (L + 1j * xi.imag) - D[j] * sc)
            Sj += complex(c) * e; Spj += complex(c) * q * muf * e
        ph = np.exp(ik * x)
        f[j] = ph * Sj / Tj
        fp[j] = ph * (ik * Sj / Tj + (Spj * Tj - Sj * Tpj) / Tj ** 2)
    u = A.alpha.T @ f
    up = A.alpha.T @ fp
    return u, up, X, (dr, res)
