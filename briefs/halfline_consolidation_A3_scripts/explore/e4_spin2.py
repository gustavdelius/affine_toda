"""Claim 5, independent of Sklyanin: the spin +-2 charges of a_n^(1) Toda (n >= 2) at a Dirichlet wall.

Book conventions: d_+- = d_t +- d_x, EOM d_+ d_- phi = -grad V, V = (m^2/beta^2) sum_j (e^{beta alpha_j.phi} - 1) (real coupling).
Ansatz (as in Corrigan-Dorey-Rietdijk-Sasaki hep-th/9404108, sec. 2): T = A(p,p,p) + B(q,p), p = d_+ phi, q = d_+^2 phi,
d_- T = d_+ Theta with Theta = -B(p, grad V)  <=>  3 A(p,p,grad V) = 2 B(p, Hess V p) for all p, phi.
Mirror current: Tbar = A(pb,pb,pb) + B(qb,pb), Thetabar = -B(pb, grad V), pb = d_- phi, d_+ Tbar = d_- Thetabar.
On x <= 0: d/dt int (T - Theta) = (T + Theta)|_0, d/dt int (Tbar - Thetabar) = -(Tbar + Thetabar)|_0.
Dirichlet phi(0,t) = phi_0: with u = d_x phi(0,t): p = u, pb = -u, q = 2 udot + grad V, qb = -2 udot + grad V.
"""
import itertools
import numpy as np
from booklib import load_classical, load_images

C = load_classical(); I = load_images()
rng = np.random.default_rng(3)
beta, m = 0.83, 1.27


def spin2_current(n, bt=beta):
    roots, Es, Hs, nj, sub = C['algebra'](f'a{n}(1)')
    r = [sub.T @ a for a in roots]                     # roots in an orthonormal basis of the field space (n-dim)
    sym3 = [c for c in itertools.combinations_with_replacement(range(n), 3)]
    asym = [c for c in itertools.combinations(range(n), 2)]
    nu = len(sym3) + len(asym)
    def unpack(vec):
        A = np.zeros((n, n, n), complex); B = np.zeros((n, n), complex)
        for v, (i, j, k) in zip(vec[:len(sym3)], sym3):
            for P in set(itertools.permutations((i, j, k))):
                A[P] = v
        for v, (i, j) in zip(vec[len(sym3):], asym):
            B[i, j] = v; B[j, i] = -v
        return A, B
    rows = []
    for col in range(nu):
        e = np.zeros(nu); e[col] = 1
        A, B = unpack(e)
        eqs = []
        for rj in r:                                   # (3/beta) A(.,.,r_j) - (B r_j) r_j^T - r_j (B r_j)^T = 0
            M = (3 / bt) * np.einsum('xyc,c->xy', A, rj) - np.outer(B @ rj, rj) - np.outer(rj, B @ rj)
            eqs.extend(M[np.triu_indices(n)])
        rows.append(eqs)
    Mat = np.array(rows).T
    U, sv, Vh = np.linalg.svd(Mat)
    nul = int((sv < 1e-10 * sv[0]).sum()) + max(0, nu - len(sv))
    sol = Vh[-1].conj()
    A, B = unpack(sol)
    return nul, A, B, r, sub


def gradV(phi, r, bt=beta, mm=m):
    return (mm**2 / bt) * sum(rj * np.exp(bt * (rj @ phi)) for rj in r)


def hessV(phi, r, bt=beta, mm=m):
    return mm**2 * sum(np.outer(rj, rj) * np.exp(bt * (rj @ phi)) for rj in r)


for n in (1, 2, 3, 4):
    nul, A, B, r, sub = spin2_current(n)
    if n == 1:
        print(f'a_1: dimension of the spin-2 current space = {nul} (sinh-Gordon has no spin-2 charge)')
        continue
    # identity check at random jets
    res = 0
    for _ in range(20):
        phi = rng.normal(size=n); p = rng.normal(size=n)
        res = max(res, abs(3 * np.einsum('abc,a,b,c->', A, p, p, gradV(phi, r)) - 2 * p @ B @ (hessV(phi, r) @ p)))
    # Dirichlet fluxes
    u = rng.normal(size=n); ud = rng.normal(size=n); phi0 = rng.normal(size=n)
    A3 = np.einsum('abc,a,b,c->', A, u, u, u)
    Fplus = 2 * A3 - 4 * u @ B @ gradV(phi0, r)
    Fminus = 4 * ud @ B @ u
    print(f'a_{n}: spin-2 current space dim {nul}; jet identity residual {res:.1e}; |A| = {np.linalg.norm(A):.3f}, |B| = {np.linalg.norm(B):.3f}; '
          f'at a Dirichlet wall: energy-like flux 2A(u,u,u) - 4B(u,gradV) = {Fplus:.3f} (no udot), momentum-like flux 4B(udot,u) = {Fminus:.3f} '
          f'(exact iff B symmetric: |B + B^T| = {np.linalg.norm(B + B.T):.1e}, |B - B^T| = {np.linalg.norm(B - B.T):.2f})')

# independent check of the conservation law on an exact (complexified) solution: the imaginary-coupling solitons of the book
# solve the real-coupling equations with beta_real = i beta_imag (V_imag(beta) = V_real(i beta)).
Sol = I['Sol']; bi = I['beta']
n = 2
nul, A, B, r, sub = spin2_current(n, bt=1j * bi)
sol = Sol(n, [1, 2], [0.5, -0.3], [0.2 + 0.4j, -0.6 + 0.9j])
x0s = np.linspace(-1.5, 1.5, 7); t0 = 0.2; hh = 1e-3


def phi_at(x, t):
    xg = np.concatenate([np.linspace(-30, x - 0.01, 3000), [x]])
    ph, _, _ = sol.phi(xg, t)
    return sub.T @ ph[:, -1]                          # n-dim components


def deriv_jet(x, t):
    """p = d_+ phi, q = d_+^2 phi, pb = d_- phi, qb = d_-^2 phi by 4th-order finite differences of the exact solution."""
    f = lambda a, b: phi_at(x + a, t + b)
    d = {}
    st = [(-2, 1 / 12), (-1, -8 / 12), (1, 8 / 12), (2, -1 / 12)]
    # first derivatives along the light-cone directions (t +- x)
    for s, key in [(1, 'p'), (-1, 'pb')]:
        d[key] = sum(c * f(s * k * hh, k * hh) for k, c in st) / hh   # d/d tau of phi(x + s tau, t + tau) = d_t + s d_x
    st2 = [(-2, -1 / 12), (-1, 16 / 12), (0, -30 / 12), (1, 16 / 12), (2, -1 / 12)]
    for s, key in [(1, 'q'), (-1, 'qb')]:
        d[key] = sum(c * f(s * k * hh, k * hh) for k, c in st2) / hh**2
    return d


def TTh(phi, d, key_p, key_q):
    p, q = d[key_p], d[key_q]
    T = np.einsum('abc,a,b,c->', A, p, p, p) + q @ B @ p
    Th = -p @ B @ gradV(phi, r, 1j * bi, I['m'])
    return T, Th


# conservation law d_- T = d_+ Theta checked by differencing T and Theta along the light-cone directions
worst = 0; scale = 0
for x0 in x0s:
    def TT(a, b):
        ph = phi_at(x0 + a, t0 + b); dj = deriv_jet(x0 + a, t0 + b)
        return TTh(ph, dj, 'p', 'q')
    k = 2e-3
    dmT = (TT(-k, k)[0] - TT(k, -k)[0]) / (2 * k)      # d_- = d_t - d_x: move (x - k, t + k) vs (x + k, t - k) ... d/ds phi(x - s, t + s)
    dpTh = (TT(k, k)[1] - TT(-k, -k)[1]) / (2 * k)     # d_+ = d_t + d_x
    worst = max(worst, abs(dmT - dpTh)); scale = max(scale, abs(dmT))
print(f'a_2 spin-2 conservation law on an exact 2-soliton solution (beta_real = i beta): |d_- T - d_+ Theta| = {worst:.1e} (|d_- T| ~ {scale:.1e})')
