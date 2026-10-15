"""Checks of prp-krein-regularization and the *Checks* paragraph after it (chapter 48): the light-cone (brickwork)
lattice built from Jimbo's U_q(sl_n^) R-matrix is exactly Krein-unitary.

G(x) = S R(x)/(xq - 1/(xq)), R as in the book (and in krein_bethe.py), q = e^{i mu}.  On 2M sites (periodic),
U = U_e U_o, U_o = prod of G over bonds (1,2),(3,4),..., U_e over (2,3),...,(2M,1); x = e^{2 Lambda}.
Pi: j -> 2M+1-j, Pi_s: j -> 2-j (mod 2M), C: colour reversal on every site, K: complex conjugation.

One PASS/FAIL line per check; exit status non-zero if one fails.
 1. part 1 for n = 2, 3, 4: G^T = G, G(x) G(1/x) = 1, conj G(x) = S G(x)^{-1} S (random mu, x);
 2. parts 2-4 for n = 3 on 6 sites: U^dag Pi U = Pi; Pi_s conj(U) Pi_s = U^{-1}; C Pi U = U C Pi;
 3. for n = 3, neither Pi nor C alone commutes with U;
 4. n = 2, on 6 and 8 sites: every eigenvalue of U on the unit circle;
 5. n = 3, on 6 and 8 sites: the spectrum of each colour sector is invariant under lambda -> 1/conj(lambda); sectors
    with at most two colours lie on the unit circle; some sector with all three colours has eigenvalues off it;
 6. the off-circle pairs appear as soon as Lambda > 0: n = 3, 6 sites, Lambda = 0.01, 0.03, 0.1, 0.3;
 7. Euclidean weights: 6 sites, n = 3, mu = 1, x = e^{i(pi-mu)/2} (|x| = 1); the ratio of Tr V^m, V = U_e U_o, to the
    same trace with the phases x^{+-1} of the weights removed (every entry of G replaced by its modulus) is
    0.78, 0.66, 0.49, 0.18 for m = 1, 2, 4, 8.
Runtime: 24 min measured on a heavily loaded machine; the 8-site sectors dominate."""
import contextlib
import io
import itertools
import sys

import numpy as np

with contextlib.redirect_stdout(io.StringIO()):  # krein_bethe prints its own checks on import
    from krein_bethe import R as Rt              # normalized R(x)/(xq - 1/(xq))

ok = True


def check(name, cond):
    global ok
    ok &= bool(cond)
    print(("PASS " if cond else "FAIL ") + name)


def swap(n):
    P = np.zeros((n * n, n * n))
    for i in range(n):
        for j in range(n):
            P[i * n + j, j * n + i] = 1
    return P


def gate(n, x, q):
    return swap(n) @ Rt(n, x, q)


def apply2(T, M, i, j, n, L):
    """Left-multiply the columns of T (vectors on (C^n)^L) by the two-site matrix M on factors (i, j)."""
    T = T.reshape((n,) * L + (-1,))
    T = np.moveaxis(T, [i, j], [0, 1])
    T = np.tensordot(M.reshape(n, n, n, n), T, axes=([2, 3], [0, 1]))
    return np.moveaxis(T, [0, 1], [i, j]).reshape(n**L, -1)


def evolve(T, Gm, n, L):
    """U T with U = U_e U_o on L = 2M sites (0-indexed bonds (0,1),(2,3),... then (1,2),...,(L-1,0))."""
    for b in range(0, L, 2):
        T = apply2(T, Gm, b, b + 1, n, L)
    for b in range(1, L, 2):
        T = apply2(T, Gm, b, (b + 1) % L, n, L)
    return T


def site_perm(perm, n, L):
    """Operator |a_0 ... a_{L-1}> -> |a'> with a'_{perm[k]} = a_k."""
    D = n**L
    P = np.zeros((D, D))
    for idx in itertools.product(range(n), repeat=L):
        new = [0] * L
        for k in range(L):
            new[perm[k]] = idx[k]
        P[np.ravel_multi_index(new, [n] * L), np.ravel_multi_index(idx, [n] * L)] = 1
    return P


def sectors(n, L):
    out = {}
    for idx in itertools.product(range(n), repeat=L):
        out.setdefault(tuple(sorted(idx)), []).append(np.ravel_multi_index(idx, [n] * L))
    return out


def sector_spectra(n, L, x, q):
    Gm = gate(n, x, q)
    res = {}
    for key, ids in sectors(n, L).items():
        T = np.zeros((n**L, len(ids)), dtype=complex)
        T[ids, range(len(ids))] = 1
        res[key] = np.linalg.eigvals(evolve(T, Gm, n, L)[ids, :])
    return res


def paired(ev, tol=1e-8):
    return max(min(abs(e - 1 / np.conj(f)) for f in ev) for e in ev) < tol


rng = np.random.default_rng(5)

# 1. part 1
for n in (2, 3, 4):
    mu = rng.uniform(0.3, 2.8)
    q, x = np.exp(1j * mu), np.exp(rng.normal())
    G, Gi, S = gate(n, x, q), gate(n, 1 / x, q), swap(n)
    e = max(np.abs(G.T - G).max(), np.abs(G @ Gi - np.eye(n * n)).max(),
            np.abs(G.conj() - S @ np.linalg.inv(G) @ S).max())
    check(f"part 1, n={n}, mu={mu:.3f}: G^T = G, G(x)G(1/x) = 1, conj G = S G^-1 S (max error {e:.1e})", e < 1e-12)

# 2-3. parts 2-4 on 6 sites, n = 3
n, L = 3, 6
mu = 1.3
q, x = np.exp(1j * mu), np.exp(2 * 0.4)
U = evolve(np.eye(n**L, dtype=complex), gate(n, x, q), n, L)
Pi = site_perm([L - 1 - k for k in range(L)], n, L)
Pis = site_perm([(-k) % L for k in range(L)], n, L)
C1 = np.eye(n)[::-1]
CL = C1
for _ in range(L - 1):
    CL = np.kron(CL, C1)
Uinv = np.linalg.inv(U)
e2 = np.abs(U.conj().T @ Pi @ U - Pi).max()
e3 = np.abs(Pis @ U.conj() @ Pis - Uinv).max()
e4 = np.abs(CL @ Pi @ U - U @ CL @ Pi).max()
check(f"part 2, n=3, 6 sites: U^dag Pi U = Pi (error {e2:.1e})", e2 < 1e-10)
check(f"part 3: Pi_s conj(U) Pi_s = U^-1 (error {e3:.1e})", e3 < 1e-10)
check(f"part 4: C Pi commutes with U (error {e4:.1e})", e4 < 1e-10)
cP, cC = np.abs(Pi @ U - U @ Pi).max(), np.abs(CL @ U - U @ CL).max()
check(f"n=3: neither Pi ({cP:.2f}) nor C ({cC:.2f}) alone commutes with U", cP > 1e-3 and cC > 1e-3)
del U, Uinv, Pi, Pis, CL

# 4. n = 2 on the unit circle
for L in (6, 8):
    sp = sector_spectra(2, L, np.exp(2 * 0.4), np.exp(1j * 1.3))
    dev = max(np.abs(np.abs(ev) - 1).max() for ev in sp.values())
    check(f"n=2, {L} sites: every eigenvalue on the unit circle (max ||l|-1| {dev:.1e})", dev < 1e-10)

# 5. n = 3 by colour content
for L in (6, 8):
    for mu, Lam in ((1.3, 0.4), (2.4, 0.8)):
        sp = sector_spectra(3, L, np.exp(2 * Lam), np.exp(1j * mu))
        pair = all(paired(ev) for ev in sp.values())
        two = max(np.abs(np.abs(ev) - 1).max() for k, ev in sp.items() if len(set(k)) <= 2)
        three = max(np.abs(np.abs(ev) - 1).max() for k, ev in sp.items() if len(set(k)) == 3)
        check(f"n=3, {L} sites, mu={mu}, Lambda={Lam}: sectors paired under l -> 1/conj l; two-colour sectors on the "
              f"circle (max dev {two:.1e}); three-colour sectors off it (max dev {three:.3f})",
              pair and two < 1e-10 and three > 1e-4)

# 6. as soon as Lambda > 0
devs = []
for Lam in (0.01, 0.03, 0.1, 0.3):
    sp = sector_spectra(3, 6, np.exp(2 * Lam), np.exp(1j * 1.3))
    devs.append(max(np.abs(np.abs(ev) - 1).max() for k, ev in sp.items() if len(set(k)) == 3))
check("n=3, 6 sites, mu=1.3: off-circle pairs at Lambda = 0.01, 0.03, 0.1, 0.3 (max dev "
      + ", ".join(f"{d:.1e}" for d in devs) + ")", min(devs) > 1e-8)

# 7. Euclidean weights
n, L, mu = 3, 6, 1.0
q, x = np.exp(1j * mu), np.exp(1j * (np.pi - mu) / 2)
G = gate(n, x, q)
V = evolve(np.eye(n**L, dtype=complex), G, n, L)
Vabs = evolve(np.eye(n**L, dtype=complex), np.abs(G), n, L)
ratios = []
for m in (1, 2, 4, 8):
    t, tabs = np.trace(np.linalg.matrix_power(V, m)), np.trace(np.linalg.matrix_power(Vabs, m))
    ratios.append(t / tabs)
print("   ratios Tr V^m / Tr |V|^m, m = 1, 2, 4, 8:", ", ".join(f"{r:.4f}" for r in ratios))
check("Euclidean weights: ratios real and equal to 0.78, 0.66, 0.49, 0.18 to two digits",
      all(abs(r.imag) < 1e-9 for r in ratios)
      and np.allclose(np.round([r.real for r in ratios], 2), [0.78, 0.66, 0.49, 0.18]))
sys.exit(0 if ok else 1)
