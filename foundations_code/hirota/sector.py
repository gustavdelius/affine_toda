"""A fluctuation sector: background soliton (polynomial tau functions of a simply-laced parent), a set of
channels (first-order data, invariant for folded theories) and the projector onto the field subspace.
Provides: exact Hirota solutions psi_b^-(x,lambda), psi_b^+(x,lambda); Wronskian matrix; Fourier-grid operator;
determinant ratio; spectrum."""
import numpy as np
import mpmath as mp
from lie import Algebra, delta_mp
from hirota import MpCtx, soliton_deg, fluct_deg
from fluctop import Background


class Sector:
    def __init__(self, A, t, mu, channels, P=None, xi=None, dps=30, label=''):
        """channels: list of (name, delta (mp list), lam_b (mp))"""
        self.A, self.t, self.mu = A, t, mu
        self.bg = Background(A, t, mu, xi)
        self.xi = self.bg.xi
        self.channels = channels
        self.P = np.eye(A.r) if P is None else P
        self.rf = self.P.shape[1]
        self.dps = dps
        self.label = label
        self.D = [len(tj) - 1 for tj in t]
        self.raw = {nm: A.alpha.T @ np.array([complex(c) for c in d]) for nm, d, lb in channels}
        self.mb2 = {nm: float(lb) for nm, d, lb in channels}

    # ---------------------------------------------------------------- exact solutions
    def coeffs(self, nm, z):
        d, lb = [(dd, l) for n2, dd, l in self.channels if n2 == nm][0]
        s, dr, res = fluct_deg(self.A, self.t, self.mu, d, lb, z, ctx=MpCtx(self.dps))
        D = self.D
        w = np.array([complex(s[j][D[j]] / self.t[j][D[j]]) for j in range(self.A.r + 1)])
        d0 = np.array([complex(c) for c in d])
        X = complex(np.vdot(d0, w) / np.vdot(d0, d0))
        return s, X, float(dr), float(res)

    def evaluate(self, s, z, nm, x):
        A = self.A
        mb = np.sqrt(self.mb2[nm])
        ik = complex(z) * mb
        muf = float(self.mu)
        xi = self.xi
        x = np.atleast_1d(x)
        L = muf * x + xi.real
        sc = np.maximum(L, 0)
        f = np.zeros((A.r + 1, len(x)), complex)
        fp = np.zeros_like(f)
        for j in range(A.r + 1):
            Tj = np.zeros(len(x), complex); Tpj = np.zeros_like(Tj)
            Sj = np.zeros_like(Tj); Spj = np.zeros_like(Tj)
            for p, c in enumerate(self.t[j]):
                e = np.exp(p * (L + 1j * xi.imag) - self.D[j] * sc)
                Tj += complex(c) * e; Tpj += complex(c) * p * muf * e
            for q, c in enumerate(s[j]):
                e = np.exp(q * (L + 1j * xi.imag) - self.D[j] * sc)
                Sj += complex(c) * e; Spj += complex(c) * q * muf * e
            ph = np.exp(ik * x)
            f[j] = ph * Sj / Tj
            fp[j] = ph * (ik * Sj / Tj + (Spj * Tj - Sj * Tpj) / Tj ** 2)
        return A.alpha.T @ f, A.alpha.T @ fp

    def kappa(self, nm, lam):
        k = np.sqrt(complex(self.mb2[nm] - lam))
        return k if k.real >= 0 else -k

    def wronskian_matrix(self, lam, xs):
        names = [c[0] for c in self.channels]
        psim, psip, Xs = {}, {}, {}
        for nm in names:
            kap = self.kappa(nm, lam)
            mb = np.sqrt(self.mb2[nm])
            zm = mp.mpc(kap / mb)
            zp = mp.mpc(-kap / mb)
            s, X, _, _ = self.coeffs(nm, zm)
            psim[nm] = self.evaluate(s, zm, nm, xs)
            Xs[nm] = X
            s2, X2, _, _ = self.coeffs(nm, zp)
            u, up = self.evaluate(s2, zp, nm, xs)
            psip[nm] = (u / X2, up / X2)
        W = np.zeros((len(names), len(names), len(xs)), complex)
        for i, b in enumerate(names):
            for j, d in enumerate(names):
                u, up = psim[b]
                v, vp = psip[d]
                W[i, j] = np.sum(u * vp - up * v, axis=0)
        return names, W, Xs

    # ---------------------------------------------------------------- grid operator
    def grid(self, N, Lp):
        x = -Lp / 2 + Lp * np.arange(N) / N
        k = 2 * np.pi * np.fft.fftfreq(N, d=Lp / N)
        F = np.fft.fft(np.eye(N), axis=0)
        D2 = np.real(np.fft.ifft(-(k ** 2)[:, None] * F, axis=0))
        Mx = self.bg.M(x)
        Mf = np.einsum('ia,xij,jb->xab', self.P, Mx, self.P)
        M0 = self.P.T @ self.A.M0 @ self.P
        rf = self.rf
        H = np.kron(-D2, np.eye(rf)).astype(complex)
        H0 = H.copy()
        for i in range(N):
            H[i * rf:(i + 1) * rf, i * rf:(i + 1) * rf] += Mf[i]
            H0[i * rf:(i + 1) * rf, i * rf:(i + 1) * rf] += M0
        trV = np.sum(np.trace(Mf - M0[None], axis1=1, axis2=2)) * (Lp / N)
        kmax = np.pi * N / Lp
        return H, H0, trV, kmax, x

    def Xprod(self, lam, Xfuncs):
        """prod_b X_b(lambda) using given closed forms Xfuncs[nm](z)"""
        p = 1.0 + 0j
        for nm, d, lb in self.channels:
            kap = self.kappa(nm, lam)
            p *= Xfuncs[nm](kap / np.sqrt(self.mb2[nm]))
        return p
