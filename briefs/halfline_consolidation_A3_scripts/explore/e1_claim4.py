"""Claim 4 exploration: sigma = dagger, imaginary coupling, book's a_t (boundary_classical.py).

Boundary data at x=0: Re phi(0) = v fixed, Im phi(0) free, d_x Re phi free, d_x Im phi = 0.
K a_t(lam) + a_t(1/lam)^dagger K = 0 (lam real).
"""
import itertools
import numpy as np
from booklib import load_classical

C = load_classical()
algebra, a_t = C['algebra'], C['a_t']
beta = C['beta']
rng = np.random.default_rng(11)


def eps_S(S, h):
    v = np.zeros(h); v[list(S)] = 1.0
    return v - len(S) / h


def sample(n, v, lam, dIm=0.0):
    roots, Es, Hs, nj, sub = algebra(f'a{n}(1)')
    w = sub @ rng.normal(size=sub.shape[1])          # Im phi(0)
    u = sub @ rng.normal(size=sub.shape[1])          # d_x Re phi(0)
    s = sub @ rng.normal(size=sub.shape[1])          # d_x Im phi(0) (should be 0)
    phi = v + 1j * w
    dx = u + 1j * dIm * s
    A = a_t(lam, phi, dx, roots, Es, Hs, nj, coupling='imag')
    B = a_t(1 / lam, phi, dx, roots, Es, Hs, nj, coupling='imag').conj().T
    return A, B


def resid(K, A, B):
    return np.linalg.norm(K @ A + B @ K) / (np.linalg.norm(K) * (np.linalg.norm(A) + np.linalg.norm(B)))


def nullspace(n, v, lam, nsamp=6, dIm=0.0):
    rows = []
    h = n + 1
    for _ in range(nsamp):
        A, B = sample(n, v, lam, dIm)
        cols = []
        for p in range(h):
            for q in range(h):
                K = np.zeros((h, h), complex); K[p, q] = 1
                cols.append((K @ A + B @ K).ravel())
        rows.append(np.array(cols).T)
    M = np.vstack(rows)
    U, sv, Vh = np.linalg.svd(M)
    nul = int((sv < 1e-10 * sv[0]).sum())
    return nul, sv, Vh[-1].conj().reshape(h, h)


for n in (1, 2, 3):
    h = n + 1
    print(f'=== a_{n}^(1)')
    for S in itertools.chain.from_iterable(itertools.combinations(range(1, h), k) for k in range(h)):
        nu = eps_S(S, h)
        shift = eps_S(tuple(rng.choice(h, size=rng.integers(0, h), replace=False)), h)  # random weight
        v = (np.pi / beta) * nu + (2 * np.pi / beta) * shift
        D = np.diag(np.exp(1j * np.pi * nu))
        r1 = max(resid(np.eye(h), *sample(n, v, lam)) for lam in (0.6, 1.7))
        rD = max(resid(D, *sample(n, v, lam)) for lam in (0.6, 1.7))
        rDbad = max(resid(D, *sample(n, v, lam, dIm=1.0)) for lam in (0.6, 1.7))
        nul, sv, K = nullspace(n, v, 0.6)
        Kn = K / K[0, 0]
        dK = np.linalg.norm(Kn - D / D[0, 0])
        # perturbed boundary value: not in (pi/beta) coweights
        pert = rng.normal(size=h); pert -= pert.mean()
        nul_p, sv_p, _ = nullspace(n, v + 0.05 * pert, 0.6)
        nul_d, _, _ = nullspace(n, v, 0.6, dIm=1.0)
        print(f'S={S!s:10} K=1: {r1:.1e}  K=D: {rD:.1e}  (D, no Neumann on Im: {rDbad:.1e})  '
              f'general K: nullity {nul} |K/K00 - D/D00| = {dK:.1e} sv_min/sv_max {sv[-1]/sv[0]:.1e}; '
              f'perturbed v: nullity {nul_p} (sv_min/sv_max {sv_p[-1]/sv_p[0]:.1e}); d_x Im != 0: nullity {nul_d}')
