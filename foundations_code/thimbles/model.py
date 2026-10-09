"""Reduced model: static lattice energy of imaginary-coupling a_n^(1) Toda on a chain.

Variables u_x = beta*phi_x in C^r, x = 1..L; boundary u_0 = 0, u_{L+1} = 2 pi lam.
    S[u] = sum_{x=0}^{L} (1/2)(u_{x+1}-u_x).(u_{x+1}-u_x) + kappa sum_{x=1}^{L} W(u_x),
    W(u) = -sum_j n_j (exp(i alpha_j.u) - 1).
Weight exp(-S/hbar), hbar = g beta^2, kappa = z beta^2. Measure N d^{rL}u with
N = (2 pi hbar)^{-r(L+1)/2}, so that Z = heat kernel at kappa = 0.
A phase eps can be given: thimbles are computed for exp(-e^{-i eps} S/|hbar|).
Thimble geometry depends only on (kappa, L, lam, eps), not on |hbar|.
"""
import numpy as np
from lie import An


class Chain:
    def __init__(self, g, L, kappa, lam, eps=0.0):
        self.g = g
        self.r = g.r
        self.L = L
        self.n = g.r * L
        self.kappa = kappa
        self.lam = np.asarray(lam, dtype=float)
        self.uL = 2 * np.pi * self.lam
        self.alpha = g.alpha  # (h, r)
        self.kac = g.kac
        self.eps = eps
        self.phase = np.exp(-1j * eps)
        self.Kmin = 0.5 * np.dot(self.uL, self.uL) / (L + 1)  # min of kinetic term on R
        self.Imax = kappa * L * np.sum(self.kac)  # bound on |Im S| on R

    # --- shapes: U is (..., L, r) complex
    def pad(self, U):
        sh = U.shape[:-2]
        z0 = np.zeros(sh + (1, self.r), dtype=complex)
        zL = np.broadcast_to(self.uL.astype(complex), sh + (1, self.r))
        return np.concatenate([z0, U, zL], axis=-2)

    def action(self, U):
        P = self.pad(U)
        D = P[..., 1:, :] - P[..., :-1, :]
        kin = 0.5 * np.sum(D * D, axis=(-1, -2))
        E = np.exp(1j * np.einsum('jr,...xr->...xj', self.alpha, U))
        pot = -self.kappa * np.sum((E - 1) * self.kac, axis=(-1, -2))
        return kin + pot

    def grad(self, U):
        P = self.pad(U)
        lap = 2 * U - P[..., :-2, :] - P[..., 2:, :]
        E = np.exp(1j * np.einsum('jr,...xr->...xj', self.alpha, U))
        gW = -1j * np.einsum('...xj,jr->...xr', E * self.kac, self.alpha)
        return lap + self.kappa * gW

    def hess(self, U):
        """Full Hessian (n x n) at a single configuration U (L, r)."""
        L, r = self.L, self.r
        H = np.zeros((L, r, L, r), dtype=complex)
        E = np.exp(1j * (U @ self.alpha.T))  # (L, h)
        for x in range(L):
            H[x, :, x, :] = 2 * np.eye(r) + self.kappa * np.einsum('j,ja,jb->ab', E[x] * self.kac, self.alpha, self.alpha)
            if x + 1 < L:
                H[x, :, x + 1, :] = -np.eye(r)
                H[x + 1, :, x, :] = -np.eye(r)
        return H.reshape(L * r, L * r)

    def derivs34(self, U):
        """Site-local cubic and quartic coefficients: E (L,h) to build S3, S4."""
        return np.exp(1j * (U @ self.alpha.T))

    # --- Newton
    def newton(self, U0, maxit=60, tol=1e-12):
        U = U0.copy()
        for it in range(maxit):
            G = self.grad(U).reshape(-1)
            nG = np.linalg.norm(G)
            if not np.isfinite(nG):
                return U, False
            if nG < tol:
                return U, True
            H = self.hess(U)
            try:
                d = np.linalg.solve(H, -G)
            except np.linalg.LinAlgError:
                return U, False
            # damping
            step = 1.0
            nd = np.linalg.norm(d)
            if nd > 2.0:
                step = 2.0 / nd
            U = U + step * d.reshape(U.shape)
            if np.max(np.abs(U.imag)) > 30:
                return U, False
        G = self.grad(U).reshape(-1)
        return U, np.linalg.norm(G) < 1e-9

    # --- thimble tangent spaces
    def tangent(self, U):
        """Return mu (n,), V (n x n complex): columns v_k with conj(Ht v_k) = mu_k v_k, mu_k > 0,
        Ht = e^{-i eps} H. J tangent = real span of V (orientation V), K tangent = i V."""
        Ht = self.phase * self.hess(U)
        A, B = Ht.real, Ht.imag
        M = np.block([[A, -B], [-B, -A]])
        w, X = np.linalg.eigh(M)
        idx = np.argsort(-w)[: self.n]
        mu = w[idx]
        Vr = X[:, idx]
        V = Vr[: self.n] + 1j * Vr[self.n:]
        return mu, V

    def one_loop(self, U, hbar, V=None, mu=None):
        """Gaussian thimble integral (with measure normalization N), orientation V. Returns complex."""
        if V is None:
            mu, V = self.tangent(U)
        n = self.n
        N = (2 * np.pi * hbar) ** (-self.r * (self.L + 1) / 2)
        S = self.action(U)
        # integral of exp(-e^{-i eps} S / hbar): Hessian scaled by phase; mu are eigenvalues for Ht
        return N * np.exp(-self.phase * S / hbar) * (2 * np.pi * hbar) ** (n / 2) * np.linalg.det(V) / np.sqrt(np.prod(mu))

    def two_loop_c1(self, U):
        """First correction c1 in Z_sigma = Z_1loop (1 + hbar c1 + ...), for action e^{-i eps} S."""
        L, r = self.L, self.r
        Ht = self.phase * self.hess(U)
        G = np.linalg.inv(Ht).reshape(L, r, L, r)
        E = np.exp(1j * (U @ self.alpha.T)) * self.kac  # (L,h)
        k = self.kappa * self.phase
        al = self.alpha
        # third derivative of kappa W at site x: i kappa sum_j E_j al_a al_b al_c
        # fourth: - kappa sum_j E_j al^4
        Gxx = np.einsum('xayb->xyab', G)  # (L,L,r,r)
        # contractions a_j^T G_xy a_k : (L,L,h,h)
        AG = np.einsum('ja,xyab,kb->xyjk', al, Gxx, al)
        diagAG = np.array([[AG[x, x, j, j] for j in range(al.shape[0])] for x in range(L)])  # (L,h)
        # quartic: -1/8 f_ijkl G_ij G_kl ; f = -k E al^4
        quart = -0.125 * np.sum(-k * E * diagAG ** 2)
        # cubic coefficient g_x,j = i k E_xj
        gc = 1j * k * E  # (L,h)
        # dumbbell: 1/8 sum g g (a G_xx a)(a_j G_xy a_k)(a_k G_yy a_k)
        dumb = 0.125 * np.einsum('xj,yk,xj,xyjk,yk->', gc, gc, diagAG, AG, diagAG)
        sun = (1.0 / 12) * np.einsum('xj,yk,xyjk->', gc, gc, AG ** 3)
        return quart + dumb + sun


def lin_interp(ch):
    L = ch.L
    xs = np.arange(1, L + 1) / (L + 1)
    return np.outer(xs, ch.uL).astype(complex)


def find_critical_points(ch, nstart=4000, seed=0, re_spread=np.pi, im_spread=1.5, extra_starts=None,
                         tol_dup=1e-6, verbose=False):
    rng = np.random.default_rng(seed)
    base = lin_interp(ch)
    found = []
    starts = []
    for k in range(nstart):
        U0 = base + rng.uniform(-re_spread, re_spread, base.shape) + 1j * rng.normal(0, im_spread, base.shape)
        starts.append(U0)
    if extra_starts is not None:
        starts.extend(extra_starts)
    hits = []
    for U0 in starts:
        U, ok = ch.newton(U0)
        if not ok:
            continue
        # dedupe
        dup = False
        for i, V in enumerate(found):
            if np.max(np.abs(V - U)) < tol_dup:
                hits[i] += 1
                dup = True
                break
        if not dup:
            found.append(U)
            hits.append(1)
    return found, hits
