"""Lie-algebraic data for a_n^(1), roots of squared length 2, in R^n coordinates.

The epsilon basis e_1..e_{n+1} of R^{n+1} is projected onto the hyperplane sum=0
with an orthonormal basis B (rows), so all vectors live in R^n.
Affine roots alpha_0..alpha_n with alpha_i = e_i - e_{i+1} (i=1..n), alpha_0 = e_{n+1} - e_1.
Kac labels n_j = 1. Weights of Lambda^a: sum_{i in S} e_i - (a/h) sum e, |S| = a.
"""
import itertools
import numpy as np


def hyperplane_basis(N):
    # orthonormal basis of {x in R^N : sum x = 0}, rows
    M = np.eye(N) - np.ones((N, N)) / N
    U, s, Vt = np.linalg.svd(M)
    B = Vt[: N - 1]
    return B


class An:
    def __init__(self, n):
        self.n = n
        self.r = n
        self.h = n + 1
        N = n + 1
        self.B = hyperplane_basis(N)
        E = np.eye(N)
        roots = [E[N - 1] - E[0]]  # alpha_0
        for i in range(n):
            roots.append(E[i] - E[i + 1])
        self.alpha_eps = np.array(roots)  # (h, N), affine roots in eps basis
        self.alpha = self.alpha_eps @ self.B.T  # (h, r)
        self.kac = np.ones(self.h)
        # fundamental weights lambda_a (a=1..n) in eps basis
        self.eps = E

    def to_r(self, v_eps):
        return np.asarray(v_eps) @ self.B.T

    def weight_from_subset(self, S):
        N = self.n + 1
        v = np.zeros(N)
        for i in S:
            v[i] += 1
        v -= len(S) / N
        return self.to_r(v)

    def fundamental_weights_rep(self, a):
        """All weights of Lambda^a (as dict label->vector in R^r); label = subset (1-based)."""
        out = {}
        for S in itertools.combinations(range(self.n + 1), a):
            out[tuple(i + 1 for i in S)] = self.weight_from_subset(S)
        return out

    def eps_coords(self, lam):
        """Return lambda in eps coordinates (sum zero)."""
        return np.asarray(lam) @ self.B

    def label(self, lam, tol=1e-6):
        """Try to identify lambda as sum_i c_i e_i - mean. Return integer eps vector (normalized min 0)."""
        v = self.eps_coords(lam)
        # v = c - mean(c) with c integer; find c by shifting so that entries are integers
        # choose shift t such that v + t is integer: t = frac
        t = -v[0] + np.round(v[0])
        best = None
        for k in range(self.n + 1):
            t = np.round(v[k]) - v[k]
            c = v + t
            if np.allclose(c, np.round(c), atol=tol):
                c = np.round(c).astype(int)
                c = c - c.min()
                best = tuple(c)
                break
        return best

    def mass_vector(self):
        # classical masses (lattice units, kappa=1): eigenvalues of sum alpha alpha^T = 4 sin^2(pi a/h)
        return np.array([2 * np.sin(np.pi * a / self.h) for a in range(1, self.n + 1)])


if __name__ == "__main__":
    for n in (2, 3):
        g = An(n)
        print("a_%d" % n, "alpha.alpha^T diag:", np.round(g.alpha @ g.alpha.T, 6))
        M2 = (g.alpha.T * g.kac) @ g.alpha
        print("mass^2 eigenvalues", np.round(np.linalg.eigvalsh(M2), 6), "expected", np.round(g.mass_vector() ** 2, 6))
        for a in range(1, n + 1):
            ws = g.fundamental_weights_rep(a)
            print(" Lambda^%d:" % a, {k: np.round(np.dot(v, v), 4) for k, v in ws.items()})
