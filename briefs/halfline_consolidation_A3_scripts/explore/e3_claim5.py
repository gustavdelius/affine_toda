"""Claim 5 exploration: Dirichlet phi(0,t) = phi_0 (all components), field-independent K, book's a_t.

K a_t(lam) + sigma(a_t(1/lam)) K = 0 for all d_x phi(0) (free), phi(0) = phi_0 fixed.
sigma: 'T' X^T, 'S' -X (Sklyanin antipode), 'dag' X^dagger, 'nconj' -conj(X); optionally P sigma(X) P^{-1}.
"""
import itertools
import numpy as np
from booklib import load_classical

C = load_classical()
algebra, a_t, eij = C['algebra'], C['a_t'], C['eij']
rng = np.random.default_rng(5)

SIG = {'T': lambda X: X.T, 'S': lambda X: -X, 'dag': lambda X: X.conj().T, 'nconj': lambda X: -X.conj()}


def nullity(rep, phi0, sig, lam, coupling='real', P=None, nsamp=6, cplx_dx=False):
    roots, Es, Hs, nj, sub = rep
    N = Es[0].shape[0]
    s0 = SIG[sig]
    s = (lambda X: P @ s0(X) @ np.linalg.inv(P)) if P is not None else s0
    rows = []
    for _ in range(nsamp):
        dx = sub @ rng.normal(size=sub.shape[1])
        if cplx_dx:
            dx = dx + 1j * (sub @ rng.normal(size=sub.shape[1]))
        A = a_t(lam, phi0, dx, roots, Es, Hs, nj, coupling)
        B = s(a_t(1 / lam, phi0, dx, roots, Es, Hs, nj, coupling))
        cols = []
        for p in range(N):
            for q in range(N):
                K = np.zeros((N, N), complex); K[p, q] = 1
                cols.append((K @ A + B @ K).ravel())
        rows.append(np.array(cols).T)
    M = np.vstack(rows)
    sv = np.linalg.svd(M, compute_uv=False)
    nul = int((sv < 1e-10 * sv[0]).sum())
    if nul == 0:
        return 0, sv[-1] / sv[0]
    U, sv, Vh = np.linalg.svd(M)
    Kr = sum(rng.normal() * Vh[-1 - k].conj() for k in range(nul)).reshape(N, N)
    svK = np.linalg.svd(Kr, compute_uv=False)
    return (nul if svK[-1] > 1e-7 * svK[0] else -nul), sv[-1] / sv[0]


# ---------- representations
def defining(n):
    return algebra(f'a{n}(1)')


def adjoint(n):
    """adjoint of sl(h) in a Hilbert-Schmidt-orthonormal basis, so rho(X^T) = rho(X)^T."""
    roots, Es, Hs, nj, sub = algebra(f'a{n}(1)')
    h = n + 1
    basis = [eij(h, i, k) for i in range(h) for k in range(h) if i != k]
    for k in range(1, h):                                   # orthonormal traceless diagonals
        d = np.zeros(h); d[:k] = 1; d[k] = -k
        basis.append(np.diag(d / np.sqrt(k * (k + 1))))
    Bm = np.array([b.ravel() for b in basis])               # rows: orthonormal basis vectors
    def rho(X):
        return np.array([[np.vdot(Bm[i], (X @ basis[j] - basis[j] @ X).ravel()) for j in range(len(basis))] for i in range(len(basis))])
    return roots, [rho(E) for E in Es], [rho(H) for H in Hs], nj, sub


def wedge2(n):
    roots, Es, Hs, nj, sub = algebra(f'a{n}(1)')
    h = n + 1
    pairs = list(itertools.combinations(range(h), 2))
    vecs = []
    for i, k in pairs:
        v = np.zeros((h, h)); v[i, k] = 1 / np.sqrt(2); v[k, i] = -1 / np.sqrt(2); vecs.append(v.ravel())
    Pm = np.array(vecs)
    def rho(X):
        G = np.kron(X, np.eye(h)) + np.kron(np.eye(h), X)
        return Pm @ G @ Pm.T
    return roots, [rho(E) for E in Es], [rho(H) for H in Hs], nj, sub


def v_plus_vbar(n):
    roots, Es, Hs, nj, sub = algebra(f'a{n}(1)')
    h = n + 1; Z = np.zeros((h, h))
    rho = lambda X: np.block([[X, Z], [Z, -X.T]])
    return roots, [rho(E) for E in Es], [rho(H) for H in Hs], nj, sub


def check_rep(rep):
    roots, Es, Hs, nj, sub = rep
    ok = all(np.allclose((E @ E.T - E.T @ E) @ E - E @ (E @ E.T - E.T @ E), 2 * E) for E in Es)   # [h_j, E_j] = 2 E_j
    for w in sub.T:
        pH = sum(w[k] * Hs[k] for k in range(len(Hs)))
        ok &= all(np.allclose(pH @ E - E @ pH, (a @ w) * E) for a, E in zip(roots, Es))
    return ok


def cyclic_P(h, k=1):
    P = np.zeros((h, h))
    for i in range(h):
        P[(i + k) % h, i] = 1
    return P


def reflect_P(h):
    """realizes the diagram reflection on the defining rep up to transposition: J_ik = (-1)^i delta_{i, h-1-k}"""
    J = np.zeros((h, h))
    for i in range(h):
        J[i, h - 1 - i] = (-1) ** i
    return J


lam_gen = 0.83
print('=== real coupling, defining representation, generic phi_0 and phi_0 = 0')
for n in (1, 2, 3, 4):
    rep = defining(n); roots, Es, Hs, nj, sub = rep; h = n + 1
    for label, phi0 in [('generic', sub @ rng.normal(size=sub.shape[1])), ('zero', np.zeros(h))]:
        res = {s: nullity(rep, phi0, s, lam_gen) for s in ('T', 'S')}
        resP = {f'{s}P{k}': nullity(rep, phi0, s, lam_gen, P=cyclic_P(h, k))[0] for s in ('T', 'S') for k in range(1, h)}
        resJ = {f'{s}J': nullity(rep, phi0, s, lam_gen, P=reflect_P(h))[0] for s in ('T', 'S')}
        resAnti = {s: nullity(rep, phi0, s, lam_gen)[0] for s in ('dag', 'nconj')}
        print(f'a_{n} phi0 {label}: T {res["T"][0]} (svmin {res["T"][1]:.1e}), S {res["S"][0]} (svmin {res["S"][1]:.1e}); '
              f'permuted {sorted(set(resP.values()))}, reflected {resJ}, antilinear on real data {resAnti}')

print('=== mechanism: lam^(2h) = 1 gives a K (cycle condition closes)')
for n in (2, 3):
    rep = defining(n); roots, Es, Hs, nj, sub = rep; h = n + 1
    phi0 = sub @ rng.normal(size=sub.shape[1])
    for lam in (np.exp(1j * np.pi / h), np.exp(1j * np.pi / (2 * h)), 1.0, -1.0):
        print(f'a_{n} lam={np.round(lam, 3)} lam^(2h)={np.round(lam**(2*h), 3)}: T {nullity(rep, phi0, "T", lam)[0]}, S {nullity(rep, phi0, "S", lam)[0]}')

print('=== other representations (real coupling, generic phi_0)')
reps = {'a2 adjoint (8)': adjoint(2), 'a3 Lambda^2 (6)': wedge2(3), 'a2 V+Vbar (6)': v_plus_vbar(2), 'a3 adjoint (15)': adjoint(3),
        'a1 adjoint (3)': adjoint(1)}
d4 = algebra('d4(1)')
reps['d4 vector (8)'] = d4
for name, rep in reps.items():
    roots, Es, Hs, nj, sub = rep
    phi0 = sub @ rng.normal(size=sub.shape[1])
    print(f'{name}: rep ok {check_rep(rep)}; T {nullity(rep, phi0, "T", lam_gen)}, S {nullity(rep, phi0, "S", lam_gen)}')

print('=== imaginary coupling, complex fields, full complex Dirichlet phi(0) = phi_0 (complex), d_x phi complex free')
for n in (1, 2, 3):
    rep = defining(n); roots, Es, Hs, nj, sub = rep
    phi0 = sub @ (rng.normal(size=sub.shape[1]) + 1j * rng.normal(size=sub.shape[1]))
    out = {s: nullity(rep, phi0, s, lam_gen, coupling='imag', cplx_dx=True)[0] for s in ('T', 'S', 'dag', 'nconj')}
    print(f'a_{n}: {out}')
