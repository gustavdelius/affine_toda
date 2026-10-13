"""Classical boundary conditions of affine Toda theory on the half-line x <= 0 (Part VII, ch. 30).

Book conventions: S = int dt [ int_{-inf}^0 L dx - B(phi(0,t)) ], so d_x phi|_0 = -grad B, with
  B_real = (2m/beta^2) sum_j C_j sqrt(2 n_j/alpha_j^2) (e^{beta alpha_j.phi/2} - 1),  B_imag = B_real(beta -> i beta).
Lax pair eq-lax-pair, a_t = (A_+ + A_-)/2, a_x = (A_+ - A_-)/2.

1. Sklyanin/BCDR compatibility K a_t(lam) + a_t(1/lam)^T K = 0 at x=0 (hep-th/9809140 (2.13), (4.1)) with
   d_x phi = -grad B_real substituted, solved for a field-independent K in a faithful representation.
   An invertible K exists exactly for the classified C (tbl-durham), for a_1^(1), a_2^(1), a_3^(1),
   d_4^(1), c_2^(1) (= b_2^(1)), c_3^(1), b_3^(1), g_2^(1), a_2^(2), a_4^(2), a_6^(2), d_3^(2) (folded reps);
   none exists for the normalization without n_j or with an overall factor 1/sqrt2 (normalization controls).
   Delius 98a (1.5) (C_j ~ sqrt(n_j alpha_j^2/2)) is integrable for a_n, c_n, a_2n^(2), d_n+1^(2), not for b_3, d_4, g_2.
   Corrections to BCDR App. D: c_2^(1) = b_2^(1) also allows (C_0, C_2 = +-1, C_1 free); d_3^(2) = a_3^(2) obeys
   the union of the d_n^(2) and a_{2n-1}^(2) rules; g_2^(1) does NOT allow C_2 (short) free -- only C = 0 or all +-1.
2. Soliton-preserving condition (sigma = hermitian conjugation, K = 1): compatibility <=> Dirichlet on
   Re phi with values in (2 pi/beta) x coweight lattice, Neumann on Im phi (imaginary coupling).
3. Theta (eq-theta-data) on the half-line: B_imag(-conj phi) = conj B_imag(phi) for all phi iff C_j real.
4. Dictionary: Delius 98a (eps = -C), Delius 98b (A_i = -C_i), Delius-Gandenberger (C = C),
   Corrigan-Delius (eps_0, eps_1) = (C_0, C_1), Ghoshal-Zamolodchikov (M_0, phi_0), Saleur-Skorik-Warner
   (zeta, eta), Baseilhac-Delius-George oscillator (C_0, C_1) = (cos p^, cos q^).
5. Boundary oscillator: energy of field + oscillator is conserved (symbolic).
6. eq-boundary-k-equation solved directly (referee pass 2026-10-11): e_6^(1) in the 27; the reductions K -> K^{-1} (C -> -C) and
   K -> D K D^{-1} (coweight sign flips); f_4^(1) and d_4^(3) by folding the 27 (d_4^(3): BCDR App. D branch C_0 = C_1 = 0,
   C_2 free has NO solution); a_5^(2), a_7^(2).
7. Sklyanin's formalism (a_1^(1), a_2^(1)): r(z) = (beta^2/4) Delius hep-th/9809140 (3.4) from the canonical brackets;
   r^(t1 t2)(z) = -r(1/z); the Durham K solves the reflection equation (2.7), K_- = K_+(-1/lam) solves (2.8) and the x_-
   condition; K_C(lam)^T ~ K_-C(-1/lam); {tau, tau} = 0 from {T, T} = [r, T T]; sigma = dagger with K = 1.
"""
import itertools
import numpy as np
import sympy as sp


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


rng = np.random.default_rng(7)

# ---------------------------------------------------------------- 1. representations


def eij(N, i, j):
    M = np.zeros((N, N)); M[i, j] = 1.0
    return M


def an_rep(n):
    """Defining rep of a_n^(1), cyclic labels alpha_j = e_j - e_{j+1} (j mod h), E_j = e_{j,j+1}."""
    h = n + 1
    roots = [np.eye(h)[j] - np.eye(h)[(j + 1) % h] for j in range(h)]
    Es = [eij(h, j, (j + 1) % h) for j in range(h)]
    Hs = [eij(h, k, k) for k in range(h)]
    return roots, Es, Hs


def d4_rep():
    """Vector rep of so(8). Basis order e1..e4, e-1..e-4. alpha_0=-e1-e2, alpha_1=e1-e2, alpha_2=e2-e3,
    alpha_3=e3-e4, alpha_4=e3+e4; Kac labels (1,1,2,1,1)."""
    N = 8
    p = lambda i: i - 1; mi = lambda i: 3 + i
    E_minus = lambda i, j: eij(N, p(i), p(j)) - eij(N, mi(j), mi(i))          # e_i - e_j
    E_plus = lambda i, j: eij(N, p(i), mi(j)) - eij(N, p(j), mi(i))           # e_i + e_j
    e = np.eye(4)
    roots = [-e[0] - e[1], e[0] - e[1], e[1] - e[2], e[2] - e[3], e[2] + e[3]]
    Es = [E_plus(1, 2).T, E_minus(1, 2), E_minus(2, 3), E_minus(3, 4), E_plus(3, 4)]
    Hs = [eij(N, p(k), p(k)) - eij(N, mi(k), mi(k)) for k in range(1, 5)]
    return roots, Es, Hs, [1, 1, 2, 1, 1]


def fold(roots, Es, orbits):
    """Folded nodes: root = orbit average, generator = sum of the orbit's generators (unnormalized)."""
    return [sum(roots[j] for j in O) / len(O) for O in orbits], [sum(Es[j] for j in O) for O in orbits]


def chevalley_normalize(roots, Es, Hs, sub=None):
    """Rescale E_j so that [E_j, E_j^T] = 2 alpha_j.H/alpha_j^2; return normalized E's (asserts it works).
    [phi.H, E_j] = (alpha_j.phi) E_j is checked for phi in the field subspace sub."""
    out = []
    for a, E in zip(roots, Es):
        aH = sum(a[k] * Hs[k] for k in range(len(Hs)))
        target = 2 * aH / (a @ a)
        comm = E @ E.T - E.T @ E
        s = np.vdot(comm.ravel(), target.ravel()) / np.vdot(comm.ravel(), comm.ravel())
        assert s > 0 and np.allclose(s * comm, target), 'not Chevalley-normalizable'
        out.append(np.sqrt(s) * E)
        sub_ = np.eye(len(Hs)) if sub is None else sub
        for w in sub_.T:
            pH = sum(w[k] * Hs[k] for k in range(len(Hs)))
            assert np.allclose(pH @ E - E @ pH, (a @ w) * E)
    return out


def algebra(name):
    """Return (roots, normalized E_j, H_k, Kac labels n_j, basis of the field subspace)."""
    if name.startswith('a') and name.endswith('(1)'):
        n = int(name[1:-3]); roots, Es, Hs = an_rep(n)
        h = n + 1; sub = np.linalg.svd(np.ones((1, h)))[2][1:].T          # sum-zero subspace
        return roots, chevalley_normalize(roots, Es, Hs), Hs, [1] * h, sub
    if name == 'd4(1)':
        roots, Es, Hs, nj = d4_rep()
        return roots, chevalley_normalize(roots, Es, Hs), Hs, nj, np.eye(4)
    if name in ('b3(1)', 'd3(2)', 'g2(1)'):  # foldings of d_4^(1) (vector rep of so(8)); d4 nodes (0,1,2,3,4)
        roots, Es, Hs, nj = d4_rep()
        orbits, perm, labels = {'b3(1)': ([[0], [1], [2], [3, 4]], {0: 0, 1: 1, 2: 2, 3: 4, 4: 3}, [1, 1, 2, 2]),
                                'd3(2)': ([[0, 1], [2], [3, 4]], {0: 1, 1: 0, 2: 2, 3: 4, 4: 3}, [1, 1, 1]),
                                'g2(1)': ([[0], [2], [1, 3, 4]], {0: 0, 2: 2, 1: 3, 3: 4, 4: 1}, [1, 2, 3])}[name]
        r, E = fold(roots, Es, orbits)
        sub = invariant_subspace(roots, perm, sumzero=False)
        return r, chevalley_normalize(r, E, Hs, sub), Hs, labels, sub
    if name == 'c3(1)':                      # fold a_5^(1): fix 0, 3; swap 1 <-> 5, 2 <-> 4
        roots, Es, Hs = an_rep(5)
        r, E = fold(roots, Es, [[0], [1, 5], [2, 4], [3]])
        sub = invariant_subspace(roots, {0: 0, 1: 5, 5: 1, 2: 4, 4: 2, 3: 3})
        return r, chevalley_normalize(r, E, Hs, sub), Hs, [1, 2, 2, 1], sub
    if name == 'c2(1)':                      # fold a_3^(1): fix 0, 2; swap 1 <-> 3   (book sec-affine-folding)
        roots, Es, Hs = an_rep(3)
        r, E = fold(roots, Es, [[0], [1, 3], [2]])
        sub = invariant_subspace(roots, {1: 3, 3: 1, 0: 0, 2: 2})
        return r, chevalley_normalize(r, E, Hs, sub), Hs, [1, 2, 1], sub
    if name == 'a2(2)':                      # fold a_2^(1): fix 1; swap 0 <-> 2.  nodes (short, long)
        roots, Es, Hs = an_rep(2)
        r, E = fold(roots, Es, [[0, 2], [1]])
        sub = invariant_subspace(roots, {0: 2, 2: 0, 1: 1})
        return r, chevalley_normalize(r, E, Hs, sub), Hs, [2, 1], sub
    if name in ('a4(2)', 'a6(2)'):           # fold a_{2n}^(1) by the reflection fixing node n: nodes long..shortest
        n = int(name[1]) // 2; roots, Es, Hs = an_rep(2 * n)
        orbits = [[n]] + [[n - k, n + k] for k in range(1, n + 1)]
        perm = {j: (2 * n - j) % (2 * n + 1) for j in range(2 * n + 1)}
        r, E = fold(roots, Es, orbits)
        sub = invariant_subspace(roots, perm)
        nj = [1] + [2] * n
        return r, chevalley_normalize(r, E, Hs, sub), Hs, nj, sub
    raise ValueError(name)


def invariant_subspace(roots, perm, sumzero=True):
    """Field subspace with alpha_j.phi = alpha_perm(j).phi (and sum-zero: a_n reps live in R^h)."""
    h = len(roots[0])
    rows = [roots[j] - roots[perm[j]] for j in perm] + ([np.ones(h)] if sumzero else [])
    u, s, vt = np.linalg.svd(np.array(rows))
    rank = int((s > 1e-9).sum())
    return vt[rank:].T


# check the labels: sum_j n_j alpha_j = 0 for every algebra used
ok = True
for name in ['a1(1)', 'a2(1)', 'a3(1)', 'd4(1)', 'c2(1)', 'c3(1)', 'b3(1)', 'g2(1)', 'a2(2)', 'a4(2)', 'a6(2)', 'd3(2)']:
    roots, Es, Hs, nj, sub = algebra(name)
    ok &= np.allclose(sum(n * a for n, a in zip(nj, roots)), 0)
    ok &= all(np.allclose(sub.T @ (a - sub @ (sub.T @ a)), 0) for a in roots)
report('representations: Chevalley-normalized generators, sum_j n_j alpha_j = 0 (a1,a2,a3,d4,c2,c3,b3,g2,a2^(2),a4^(2),a6^(2),d3^(2))', ok)

# ---------------------------------------------------------------- 2. K-matrix compatibility (Sklyanin / BCDR (3.14))
beta, m = 0.83, 1.27


def a_t(lam, phi, dxphi, roots, Es, Hs, nj, coupling='real'):
    """time component of the book's Lax pair at one point: a_t = (A_+ + A_-)/2."""
    b = beta if coupling == 'real' else 1j * beta
    out = (b / 2) * sum(dxphi[k] * Hs[k] for k in range(len(Hs)))
    for a, E, n in zip(roots, Es, nj):
        c = np.sqrt(n * (a @ a) / 2)
        out = out + (m / 2) * c * np.exp(b * (a @ phi) / 2) * (lam * E - E.T / lam)
    return out


def gradB(phi, C, roots, nj, norm, coupling='real'):
    b = beta if coupling == 'real' else 1j * beta
    # B = (2m/b^2) sum_j C_j N_j (e^{b alpha_j.phi/2}-1)  ->  grad B = (m/b) sum_j C_j N_j alpha_j e^{...}
    return (m / b) * sum(Cj * norm(n, a) * a * np.exp(b * (a @ phi) / 2) for Cj, n, a in zip(C, nj, roots))


book_norm = lambda n, a: np.sqrt(2 * n / (a @ a))


def K_nullity(name, C, norm=book_norm, coupling='real', nsamp=4):
    """dimension of {K : K a_t(lam) + a_t(1/lam)^T K = 0 for all sampled fields}, with d_x phi = -grad B.
    Returns 0 if there is no solution, minus the dimension if no solution is invertible."""
    roots, Es, Hs, nj, sub = algebra(name)
    N = Es[0].shape[0]
    lam = 0.7 + 0.4j
    rows = []
    for s in range(nsamp):
        phi = sub @ (rng.normal(size=sub.shape[1]) + 0.3j * rng.normal(size=sub.shape[1]))
        dx = -gradB(phi, C, roots, nj, norm, coupling)
        A = a_t(lam, phi, dx, roots, Es, Hs, nj, coupling)
        Bm = a_t(1 / lam, phi, dx, roots, Es, Hs, nj, coupling).T
        cols = []
        for p in range(N):
            for q in range(N):
                K = eij(N, p, q)
                cols.append((K @ A + Bm @ K).ravel())
        rows.append(np.array(cols).T)
    M = np.vstack(rows)
    U, sv, Vh = np.linalg.svd(M)
    nul = int((sv < 1e-9 * sv[0]).sum())
    if nul == 0:
        return 0
    # an invertible K must exist in the solution space (reducible reps: K is free on trivial summands)
    v = sum(rng.normal() * Vh[-1 - k].conj() for k in range(nul)).reshape(N, N)
    svK = np.linalg.svd(v, compute_uv=False)
    return nul if svK[-1] > 1e-6 * svK[0] else -nul


signs = lambda r: [np.array(s, float) for s in itertools.product([1, -1], repeat=r)]

cases = {
    'a1(1)': ([np.zeros(2), np.ones(2), np.array([0.37, -1.9]), np.array([2.1 + 0.4j, 0.3 - 1j]), np.array([0, 1.5])], []),
    'a2(1)': ([np.zeros(3)] + signs(3), [np.array([1, 1, 0.]), np.array([.5, .5, .5]), np.array([1, -1, 0.3]), np.array([0, 0, 1.])]),
    'a3(1)': ([np.zeros(4)] + signs(4), [np.array([1, 1, 1, 0.]), np.array([1, 0, 1, 0.]), np.array([.9, 1, 1, 1])]),
    'd4(1)': ([np.zeros(5)] + signs(5), [np.array([1, 1, 0, 1, 1.]), np.array([0, 0, 1, 0, 0.]),
                                          np.array([1, 1, np.sqrt(2), 1, 1]), np.array([1, 1, 1 / np.sqrt(2), 1, 1])]),
    # c_2^(1) = b_2^(1), nodes (0, 1, 2) = (long, short, long): BCDR c_n rule AND b_n rule (long = +-1, short free)
    'c2(1)': ([np.zeros(3)] + signs(3) + [np.array([0.4, 0, -1.7]), np.array([1, 0, 0.]), np.array([2 + 1j, 0, 0.3]),
               np.array([1, 0.5, 1]), np.array([1, 0.5, -1]), np.array([-1, 2.3 + 1j, 1])],
              [np.array([1, 0.5, 0.9]), np.array([1, 1, 0.]), np.array([0.3, 1, -1]), np.array([0, 1, 0.])]),
    # c_3^(1), nodes (0,1,2,3) = (long, short, short, long)
    'c3(1)': ([np.zeros(4), np.ones(4), np.array([1, -1, 1, -1.]), np.array([0.4, 0, 0, 2.]), np.array([1, 0, 0, 1.])],
              [np.array([1, 0.5, 1, 1]), np.array([1, 1, 0.5, 1]), np.array([1, 0.5, 0.5, 1]), np.array([0.3, 1, 1, 1])]),
    # b_3^(1), nodes (0,1,2,3) = (long, long, long, short)
    'b3(1)': ([np.zeros(4), np.ones(4), np.array([1, -1, -1, 0.3]), np.array([-1, 1, -1, 2 + 1j]), np.array([1, 1, 1, 0.])],
              [np.array([1, 1, 0.5, 1]), np.array([0, 0, 0, 1.]), np.array([0.5, 1, 1, 1]), np.array([1, 1, 0, 0.3])]),
    # g_2^(1), nodes (0, 2, {1,3,4}) = (long, long, short).  BCDR App. D: 'C_0, C_1 = +-1, C_2 arbitrary';
    # the full solve finds only C_2 = +-1 (the bad list contains BCDR's C_2 = 0.3, 0, 2-i)
    'g2(1)': ([np.zeros(3)] + signs(3),
              [np.array([1, -1, 0.3]), np.array([-1, -1, 0.]), np.array([1, 1, 2 - 1j]), np.array([1, 1, 0.5]),
               np.array([1, 0.5, 1]), np.array([0, 0, 1.]), np.array([0, 1, 1.]), np.array([0.5, 1, 0.2])]),
    # a_2^(2) nodes (0, 1) = (short, long)
    'a2(2)': ([np.zeros(2), np.array([0.6, 1]), np.array([-2.2, -1]), np.array([0, 0.7]), np.array([1, 1.])],
              [np.array([0.5, 0.5]), np.array([1, 0.5])]),
    # a_4^(2) nodes (0, 1, 2) = (long, middle, shortest)
    'a4(2)': ([np.zeros(3), np.array([1, -1, 0.4]), np.array([-1, -1, 2.3]), np.array([1, 0.6, 0]),
               np.array([-1, -1.7, 0]), np.array([0.8, 0, 0])],
              [np.array([0.5, 1, 0.3]), np.array([1, 0.5, 0.3]), np.array([0.5, 0, 0.3]), np.array([1, 0, 0.3]),
               np.array([0, 1, 0.3]), np.array([0, 0, 1.])]),
    # d_3^(2) = a_3^(2), nodes ({0,1}, 2, {3,4}) = (short, long, short): BCDR d_n^(2) rule (middle +-1, ends free)
    # AND a_{2n-1}^(2) rule at n=2 (shorts 0, long free)
    'd3(2)': ([np.zeros(3), np.ones(3), np.array([0.3, 1, -2.]), np.array([1 + 1j, -1, 0.]), np.array([0, 1, 0.]),
               np.array([0, 0.5, 0.]), np.array([0, 2 + 1j, 0.])],
              [np.array([1, 0.5, 1]), np.array([0.3, 0, 1]), np.array([1, 0, 0.]), np.array([0.3, 0.5, 0.])]),
    # a_6^(2) nodes (0,1,2,3) = (long, middle, middle, shortest)
    'a6(2)': ([np.zeros(4), np.array([1, -1, 1, 0.4]), np.array([0.8, 0, 0, 0]), np.array([1, -1, 1, 0])],
              [np.array([1, 0.6, 0, 0]), np.array([1, -1, 0.6, 0]), np.array([1, 1, 0, 0]), np.array([0.5, 1, 1, 0.3])]),
}
for name, (good, bad) in cases.items():
    ng = [K_nullity(name, C) for C in good]
    nb = [K_nullity(name, C) for C in bad]
    d0 = ng[0]          # C = 0: K = 1 and the commutant of the representation (1 if irreducible)
    report(f'{name}: invertible K exists for the {len(good)} classified C, not for {len(bad)} others  '
           f'[nullities: good {sorted(set(ng))}, bad {sorted(set(nb))}]', all(x == d0 for x in ng) and all(x <= 0 for x in nb))

# imaginary coupling: same statement with B_imag and the beta -> i beta Lax pair
ok = all(K_nullity('a2(1)', C, coupling='imag') == 1 for C in signs(3)) and K_nullity('a2(1)', np.array([1, 1, 0.]), coupling='imag') == 0
ok &= all(K_nullity('d4(1)', C, coupling='imag') == 1 for C in signs(5)[:6])
report('imaginary coupling (B_imag, beta -> i beta): same classification (a2, d4)', ok)

# normalization controls on d_4^(1) (n_2 = 2) and a_2^(1)
no_n = lambda n, a: np.sqrt(2 / (a @ a))                    # drop n_j
half = lambda n, a: np.sqrt(2 * n / (a @ a)) / np.sqrt(2)   # overall factor 1/sqrt2 (normalization control; BCDR (1.4)
                                                             # prints |A_i| = 2 sqrt(n_i), i.e. our |C_i| = 1)
ok = all(K_nullity('d4(1)', C, norm=no_n) == 0 for C in signs(5)[:8])
ok &= K_nullity('d4(1)', np.ones(5), norm=half) == 0 and K_nullity('a2(1)', np.ones(3), norm=half) == 0
report('normalization controls: without the factor sqrt(n_j), or with an overall factor 1/sqrt2, no K exists for |C_j| = 1 (d4, a2)', ok)

# Delius 98a (1.5): C_j proportional to sqrt(n_j alpha_j^2/2) (single eps, largest |C_j| = 1)
d15 = {}
for name in ['a2(1)', 'c3(1)', 'a2(2)', 'a4(2)', 'a6(2)', 'd3(2)', 'b3(1)', 'd4(1)', 'g2(1)']:
    roots, Es, Hs, nj, sub = algebra(name)
    Cd = np.sqrt(np.array(nj) * np.array([a @ a for a in roots]) / 2)
    d15[name] = K_nullity(name, Cd / Cd.max())
report('Delius 98a (1.5) (C_j ~ sqrt(n_j alpha_j^2/2)): integrable for a2, c3, a2^(2), a4^(2), a6^(2), d3^(2); not for b3, d4, g2  '
       f'[nullities {d15}]', all(d15[k] > 0 for k in ['a2(1)', 'c3(1)', 'a2(2)', 'a4(2)', 'a6(2)', 'd3(2)'])
       and all(d15[k] <= 0 for k in ['b3(1)', 'd4(1)', 'g2(1)']))

# phi = 0 is a solution of the boundary condition iff sum_j C_j sqrt(2n_j/alpha_j^2) alpha_j = 0
roots, Es, Hs, nj, sub = algebra('d4(1)')
vac = [s for s in signs(5) if np.allclose(sum(c * book_norm(n, a) * a for c, n, a in zip(s, nj, roots)), 0)]
roots2, _, _, nj2, _ = algebra('a3(1)')
vac2 = [s for s in signs(4) if np.allclose(sum(c * a for c, a in zip(s, roots2)), 0)]
report(f'phi=0 satisfies the C=+-1 boundary condition: d4 for {len(vac)} of 32 sign choices, a3 for {len(vac2)} of 16 (C = +-(1,1,1,1))', len(vac) == 0 and len(vac2) == 2)

# ---------------------------------------------------------------- 3. soliton-preserving condition (sigma = dagger, K = 1)
lam_s = sp.symbols('lambda', positive=True)
bs, ms = sp.symbols('beta m', positive=True)
u = sp.symbols('u0:3', real=True); v = sp.symbols('v0:3', real=True)
ux = sp.symbols('ux0:3', real=True); vx = sp.symbols('vx0:3', real=True)
h3 = 3
al = [sp.Matrix([1 if k == j else (-1 if k == (j + 1) % h3 else 0) for k in range(h3)]) for j in range(h3)]
Esym = [sp.Matrix(3, 3, lambda r, c: 1 if (r, c) == (j, (j + 1) % 3) else 0) for j in range(3)]
phi = sp.Matrix([u[k] + sp.I * v[k] for k in range(3)])
dphi = sp.Matrix([ux[k] + sp.I * vx[k] for k in range(3)])


def at_sym(lmb, coupling):
    b = bs if coupling == 'real' else sp.I * bs
    out = (b / 2) * sp.diag(*dphi)
    for j in range(3):
        out += (ms / 2) * sp.exp(b * (al[j].T * phi)[0] / 2) * (lmb * Esym[j] - Esym[j].T / lmb)
    return out


for coupling, cond in [('imag', 'Dirichlet on Re phi in (2pi/beta) coweights, Neumann on Im phi'),
                       ('real', 'Neumann on Re phi, Dirichlet on Im phi in (2pi/beta) coweights')]:
    S = at_sym(lam_s, coupling) + at_sym(1 / lam_s, coupling).H       # K = 1, sigma = hermitian conjugation
    S = sp.simplify(S.applyfunc(lambda z: sp.expand_complex(z)))
    diag_part = [sp.simplify(S[k, k]) for k in range(3)]
    off = [sp.simplify(S[j, (j + 1) % 3] / lam_s) for j in range(3)]
    b = bs if coupling == 'real' else sp.I * bs
    # expected: diagonal = beta*d_x Re phi (real) or -beta*d_x Im phi (imag); off = i m Im(e^{b alpha.phi/2})
    exp_diag = [bs * ux[k] if coupling == 'real' else -bs * vx[k] for k in range(3)]
    exp_off = [sp.I * ms * sp.im(sp.expand_complex(sp.exp(b * (al[j].T * phi)[0] / 2))) for j in range(3)]
    ok = all(sp.simplify(d - e) == 0 for d, e in zip(diag_part, exp_diag)) and all(sp.simplify(o - e) == 0 for o, e in zip(off, exp_off))
    # zero set of Im e^{i beta alpha_j.phi/2} = e^{-beta alpha.v/2} sin(beta alpha_j.u/2): beta alpha_j.Re phi in 2 pi Z
    report(f'soliton-preserving K=1, {coupling} coupling: K a_t + a_t(1/lam)^dagger = 0 <=> {cond}', ok)

# a numerical illustration of the 'iff': pick Re phi = (2 pi/beta) * coweight, check vanishing; shift by half -> fails with K = 1 (it holds with K = e^{i pi nu.H}: boundary_dirichlet.py)
bv, mv = 0.9, 1.0
w1 = np.array([2, -1, -1]) / 3.0                          # fundamental coweight of a_2 (alpha_1.w1 = 1 in cyclic labels)
alv = [np.array(a, dtype=float).ravel() for a in [[1, -1, 0], [0, 1, -1], [-1, 0, 1]]]
w1 = np.array([2 / 3, -1 / 3, -1 / 3])
def sp_res(Rephi, Imphi, dReal, dImag):
    phi_ = Rephi + 1j * Imphi
    res = [-bv * dImag[k] for k in range(3)] + [np.imag(np.exp(1j * bv * (a @ phi_) / 2)) for a in alv]
    return max(abs(np.array(res)))
ok = sp_res(2 * np.pi / bv * w1, np.array([0.3, -0.1, -0.2]), np.zeros(3), np.zeros(3)) < 1e-14
ok &= sp_res(np.pi / bv * w1, np.array([0.3, -0.1, -0.2]), np.zeros(3), np.zeros(3)) > 0.1
report('soliton-preserving: Re phi(0) = (2pi/beta) w passes, (pi/beta) w (half a coweight) fails with K_+ = 1', ok)

# ---------------------------------------------------------------- 4. Theta on the half-line
Cs = sp.symbols('C0:3'); Cb = [sp.conjugate(c) for c in Cs]
def Bimag(ph, C):
    return -(2 * ms / bs**2) * sum(C[j] * (sp.exp(sp.I * bs * (al[j].T * ph)[0] / 2) - 1) for j in range(3))
phis = sp.Matrix([u[k] + sp.I * v[k] for k in range(3)])
diff = sp.expand(sp.conjugate(Bimag(phis, Cs)) - Bimag(-phis.conjugate(), Cs))
# diff = -(2m/beta^2) sum_j (conj C_j - C_j)(e^{-i beta alpha_j.conj(phi)/2} - 1): vanishes for all phi iff C_j real
expect = -(2 * ms / bs**2) * sum((Cb[j] - Cs[j]) * (sp.exp(-sp.I * bs * (al[j].T * phis.conjugate())[0] / 2) - 1) for j in range(3))
report('Theta: conj B_imag(phi) - B_imag(-conj phi) = -(2m/beta^2) sum_j (conj C_j - C_j)(e^{-i beta alpha_j.conj phi/2}-1)  => Theta-invariant iff C_j real',
       sp.simplify(diff - sp.expand(expect)) == 0)

# ---------------------------------------------------------------- 5. dictionary (a_1^(1) and a_n^(1))
x = sp.symbols('x', real=True)
f = sp.symbols('f')          # phi(0) for a_1^(1): alpha_1 = sqrt2, alpha_0 = -sqrt2, n_j = 1
C0, C1 = sp.symbols('C_0 C_1')
r2 = sp.sqrt(2)
def a1_B(coupling):
    b = bs if coupling == 'real' else sp.I * bs
    return (2 * ms / b**2) * (C1 * (sp.exp(b * r2 * f / 2) - 1) + C0 * (sp.exp(-b * r2 * f / 2) - 1))
dxphi_real = -sp.diff(a1_B('real'), f)
dxphi_imag = -sp.diff(a1_B('imag'), f)

# Corrigan-Delius hep-th/9909145 (2.4): d_x phi = (sqrt2 m/beta)(eps0 e^{-beta phi/sqrt2} - eps1 e^{beta phi/sqrt2})
e0, e1 = sp.symbols('epsilon_0 epsilon_1')
cd = (r2 * ms / bs) * (e0 * sp.exp(-bs * f / r2) - e1 * sp.exp(bs * f / r2))
report('Corrigan-Delius (2.4): (eps_0, eps_1) = (C_0, C_1), same phi, beta, m',
       sp.simplify(dxphi_real.subs({C0: e0, C1: e1}) - cd) == 0)

# Ghoshal-Zamolodchikov (5.8): B = -M cos(beta_SG (phi - phi0)/2), beta_SG = sqrt2 beta, m_0 = 2m, same phi
M0, ph0 = sp.symbols('M_0 phi_0', real=True)
bSG, m0 = r2 * bs, 2 * ms
gz = -M0 * sp.cos(bSG * (f - ph0) / 2)
Cgz = {C0: bSG**2 * M0 / (4 * m0) * sp.exp(sp.I * bSG * ph0 / 2), C1: bSG**2 * M0 / (4 * m0) * sp.exp(-sp.I * bSG * ph0 / 2)}
report('Ghoshal-Zamolodchikov (5.8): C_{0,1} = (beta_SG^2 M_0/4 m_0) e^{+-i beta_SG phi_0/2}  (B_imag - B_GZ = const)',
       sp.simplify(sp.expand(sp.diff(a1_B('imag').subs(Cgz) - gz, f).rewrite(sp.exp))) == 0)

# Saleur-Skorik-Warner (2.14): mu = M cos(phi0/2) = 2 cosh zeta cos eta, nu = M sin(phi0/2) = 2 sinh zeta sin eta,
# their M = M_0 beta_SG^2/(2 m_0) (field on x>=0, units m_0=1), phi0_SSW = beta_SG phi_0
ze, et = sp.symbols('zeta eta', real=True)
MSSW = 2 * sp.sqrt((sp.cosh(ze) * sp.cos(et))**2 + (sp.sinh(ze) * sp.sin(et))**2)
C0ssw = (MSSW / 2) * (sp.cosh(ze) * sp.cos(et) + sp.I * sp.sinh(ze) * sp.sin(et)) / (MSSW / 2)
ok = sp.simplify(sp.expand_complex(sp.cosh(ze + sp.I * et)) - (sp.cosh(ze) * sp.cos(et) + sp.I * sp.sinh(ze) * sp.sin(et))) == 0
# with C_0 = (mu + i nu)/2, C_1 = (mu - i nu)/2 [from GZ map with M_SSW e^{i phi0/2}/2]:
ok &= sp.simplify((2 * sp.cosh(ze) * sp.cos(et) + sp.I * 2 * sp.sinh(ze) * sp.sin(et)) / 2 - sp.expand_complex(sp.cosh(ze + sp.I * et))) == 0
report('Saleur-Skorik-Warner (2.14): C_0 = cosh(zeta + i eta), C_1 = cosh(zeta - i eta); C=1 <-> (0,0), Neumann <-> (0, pi/2)',
       ok and sp.cosh(sp.I * sp.pi / 2).rewrite(sp.cos) == 0)

# Delius 98a (1.5) for a_n^(1): phi_D = i beta phi (imaginary coupling), m = 1: d_x phi_D = eps sum_i alpha_i e^{alpha_i.phi_D/2}
# book: d_x phi = (i m/beta) sum_j C_j alpha_j e^{i beta alpha_j.phi/2}  =>  d_x phi_D = -m sum_j C_j alpha_j e^{alpha_j.phi_D/2}
pD = sp.symbols('pD0:3')
lhs_book = [sp.I * bs * (sp.I * ms / bs) * sum(Cs[j] * al[j][k] * sp.exp((al[j].T * sp.Matrix(pD))[0] / 2) for j in range(3)) for k in range(3)]
eps = sp.symbols('epsilon')
delius = [ms * eps * sum(al[j][k] * sp.exp((al[j].T * sp.Matrix(pD))[0] / 2) for j in range(3)) for k in range(3)]
report('Delius 98a (1.5) (a_n, |alpha|^2=2): eps = -C  (phi_D = i beta phi); Delius 98b (1.2): A_i = -C_i',
       all(sp.simplify(lhs_book[k].subs({Cs[0]: -eps, Cs[1]: -eps, Cs[2]: -eps}) - delius[k]) == 0 for k in range(3)))

# Delius-Gandenberger hep-th/9904002 (1.2): beta d_x phi + m sum_i C_i alpha_i e^{beta alpha_i.phi/2} = 0, same real field
phr = sp.symbols('p0:3', real=True)
book_real = [-(ms / bs) * sum(Cs[j] * al[j][k] * sp.exp(bs * (al[j].T * sp.Matrix(phr))[0] / 2) for j in range(3)) for k in range(3)]
dg = [-(ms / bs) * sum(Cs[j] * al[j][k] * sp.exp(bs * (al[j].T * sp.Matrix(phr))[0] / 2) for j in range(3)) for k in range(3)]
report('Delius-Gandenberger (1.2): C_i(DG) = C_i(book) for a_n^(1); "(+...+)" = C = +1 = Delius-98a eps = -1',
       all(sp.simplify(a - b) == 0 for a, b in zip(book_real, dg)))

# Baseilhac-Delius-George nlin/0201007 (2): H_osc = -(2 m_0/beta_SG^2)(e^{-i b phi/2} cos p^ + e^{i b phi/2} cos q^)
ph_, qh_ = sp.symbols('phat qhat', real=True)
Hosc = -(2 * m0 / bSG**2) * (sp.exp(-sp.I * bSG * f / 2) * sp.cos(ph_) + sp.exp(sp.I * bSG * f / 2) * sp.cos(qh_))
report('Baseilhac-Delius-George (2): H_osc = B_imag with C_0 = cos p^, C_1 = cos q^ (beta_SG = sqrt2 beta, m_0 = 2m)',
       sp.simplify(sp.diff(a1_B('imag').subs({C0: sp.cos(ph_), C1: sp.cos(qh_)}) - Hosc, f)) == 0)

# ---------------------------------------------------------------- 6. boundary oscillator: energy conservation
# dE/dt = phi_t(0) phi_x(0) + dH_osc/dt with phi_x(0) = -dH_osc/dphi and {p^, q^} = beta_SG^2/4:
#   dp^/dt = -(beta_SG^2/4) dH/dq^,  dq^/dt = +(beta_SG^2/4) dH/dp^   (BDG (9), (10))
pt = sp.symbols('phi_t')
pdot = -(bSG**2 / 4) * sp.diff(Hosc, qh_)
qdot = (bSG**2 / 4) * sp.diff(Hosc, ph_)
ok = sp.simplify(pdot - (-(m0 / 2) * sp.exp(sp.I * bSG * f / 2) * sp.sin(qh_))) == 0
ok &= sp.simplify(qdot - ((m0 / 2) * sp.exp(-sp.I * bSG * f / 2) * sp.sin(ph_))) == 0
dEdt = pt * (-sp.diff(Hosc, f)) + sp.diff(Hosc, f) * pt + sp.diff(Hosc, ph_) * pdot + sp.diff(Hosc, qh_) * qdot
report('BDG oscillator: equations (9), (10) reproduced; d/dt (E_field + H_osc) = 0 with BC d_x phi = -dH_osc/dphi', ok and sp.simplify(dEdt) == 0)
# BC (7) of BDG in book variables
bc7 = -(sp.I * m0 / bSG) * (sp.exp(-sp.I * bSG * f / 2) * sp.cos(ph_) - sp.exp(sp.I * bSG * f / 2) * sp.cos(qh_))
report('BDG (7): d_x phi(0) = -(i m_0/beta_SG)(e^{-i beta_SG phi/2} cos p^ - e^{i beta_SG phi/2} cos q^) = -dH_osc/dphi',
       sp.simplify(bc7 + sp.diff(Hosc, f)) == 0)

# ---------------------------------------------------------------- 7. eq-boundary-k-equation solved directly: exceptional and twisted rows
# K(-C_j h_j + lam E_j - F_j/lam) = (C_j h_j + lam E_j - F_j/lam) K,  F_j = E_j^T, h_j = [E_j, F_j], [h_j, E_j] = 2 E_j.
# e_6^(1) in the 27 (minuscule, Bourbaki labels, affine node E_0 = lowering by the highest root); f_4^(1) and d_4^(3) by
# folding the 27 (Z_2: 1<->6, 3<->5; Z_3: arms {0,1,6}, {2,3,5}, node 4 fixed); a_5^(2), a_7^(2) in the 2n-dim rep of the
# twisted loop algebra of sl(2n).  (Referee pass 2026-10-11; e_6^(2) in the 56 of e_7 was not finished there.)
import scipy.sparse as sps


def kn_normalize(Es):
    out = []
    for E in Es:
        h_ = E @ E.T - E.T @ E
        c = np.vdot((h_ @ E - E @ h_).ravel(), E.ravel()).real / np.vdot(E.ravel(), E.ravel()).real
        assert c > 0 and np.allclose(h_ @ E - E @ h_, c * E)
        out.append(E * np.sqrt(2 / c))
    return out


def kn_cartan(Es):
    """A_ij from [h_i, E_j] = A_ij E_j; asserts [E_i, F_j] = 0 for i != j."""
    A = np.zeros((len(Es), len(Es)))
    for i, Ei in enumerate(Es):
        h_ = Ei @ Ei.T - Ei.T @ Ei
        for j, Ej in enumerate(Es):
            c = h_ @ Ej - Ej @ h_
            A[i, j] = np.vdot(c.ravel(), Ej.ravel()).real / np.vdot(Ej.ravel(), Ej.ravel()).real
            assert np.allclose(c, A[i, j] * Ej)
            if i != j:
                assert np.allclose(Ei @ Ej.T - Ej.T @ Ei, 0)
    return np.round(A, 9)


def kn_null(Es, C, lam=0.83 + 0.37j):
    """nullity of eq-boundary-k-equation (all nodes); minus the nullity if no solution is invertible."""
    N = Es[0].shape[0]; I = sps.identity(N, format='csr'); G = 0
    for E, c in zip(Es, C):
        h_ = E @ E.T - E.T @ E
        X = -c * h_ + lam * E - E.T / lam; Y = c * h_ + lam * E - E.T / lam
        L = sps.kron(sps.csr_matrix(X.T), I) - sps.kron(I, sps.csr_matrix(Y))      # column-major vec(K X - Y K)
        G = G + (L.conj().T @ L)
    w, V = np.linalg.eigh(G.toarray())
    nul = int((w < 1e-10 * w[-1]).sum())
    if nul == 0:
        return 0
    Kr = (V[:, :nul] @ rng.normal(size=nul)).reshape(N, N, order='F')
    sv = np.linalg.svd(Kr, compute_uv=False)
    return nul if sv[-1] > 1e-7 * sv[0] else -nul


def kn_K(Es, C, lam):
    """the (unique up to scale) solution K for an irreducible rep."""
    N = Es[0].shape[0]; I = np.eye(N); G = 0
    for E, c in zip(Es, C):
        h_ = E @ E.T - E.T @ E
        L = np.kron((-c * h_ + lam * E - E.T / lam).T, I) - np.kron(I, c * h_ + lam * E - E.T / lam)
        G = G + L.conj().T @ L
    return np.linalg.eigh(G)[1][:, 0].reshape(N, N, order='F')


def e6_27():
    A = 2 * np.eye(6)
    for i, j in [(1, 3), (3, 4), (4, 5), (5, 6), (2, 4)]:
        A[i - 1, j - 1] = A[j - 1, i - 1] = -1
    W = [(1, 0, 0, 0, 0, 0)]; idx = {W[0]: 0}; k = 0                       # weights in Dynkin labels
    while k < len(W):
        mu = np.array(W[k])
        for i in range(6):
            if mu[i] == 1:
                nu = tuple((mu - A[i]).astype(int))
                if nu not in idx:
                    idx[nu] = len(W); W.append(nu)
        k += 1
    N = len(W); Es = []
    for i in range(6):
        E = np.zeros((N, N))
        for nu, bb in idx.items():
            if nu[i] == -1:
                E[idx[tuple((np.array(nu) + A[i]).astype(int))], bb] = 1
        Es.append(E)
    F = {tuple(np.eye(6, dtype=int)[i]): Es[i].T for i in range(6)}; front = list(F)
    while front:                                                             # lowering operators F_beta, beta > 0
        new = []
        for bt in front:
            for i in range(6):
                c = tuple(np.array(bt) + np.eye(6, dtype=int)[i])
                if c not in F:
                    M = Es[i].T @ F[bt] - F[bt] @ Es[i].T
                    if np.linalg.norm(M) > 1e-9:
                        F[c] = M; new.append(c)
        front = new
    theta = max(F, key=sum)
    rootc = np.array([np.linalg.solve(A, np.array(w, float)) for w in W])   # weights in simple-root coordinates
    return kn_normalize([F[theta]] + Es), np.array((1,) + theta), rootc


def summands(Es):
    """invariant subspaces of a reducible rep from a generic element of the commutant."""
    N = Es[0].shape[0]; I = np.eye(N)
    L = np.vstack([np.kron(g.T, I) - np.kron(I, g) for g in Es + [E.T for E in Es]])
    w, V = np.linalg.eigh(L.T @ L)
    Z = (V[:, w < 1e-9 * w[-1]] @ rng.normal(size=int((w < 1e-9 * w[-1]).sum()))).reshape(N, N, order='F')
    ev, U = np.linalg.eigh(Z + Z.T)
    groups = []
    for k in range(N):
        if groups and abs(ev[k] - ev[groups[-1][0]]) < 1e-6:
            groups[-1].append(k)
        else:
            groups.append([k])
    return [U[:, g] for g in groups]


E6, n6, rootc = e6_27()
A6 = kn_cartan(E6)
ok = np.allclose(A6 @ n6, 0) and np.allclose(n6, [1, 1, 2, 2, 3, 2, 1])
# uniform +1, 5 random sign patterns covering both classes prod_j C_j^{n_j} = +-1, and C = 0 admit an invertible K;
# off-list values do not.  (With the reductions below this covers all 128 sign patterns.)
pats = [np.ones(7)] + [np.array(p, float) for p in rng.choice([1, -1], size=(5, 7))]
pats[1][0] = -np.prod(pats[1][1:] ** n6[1:])                          # make sure both classes occur
classes = sorted(set(int(np.prod(p ** n6)) for p in pats))
ng = [kn_null(E6, p) for p in pats] + [kn_null(E6, np.zeros(7))]
bad6 = [np.array([1, 1, 1, 0.5, 1, 1, 1]), np.array([1, 1, 1, 1, 0, 1, 1.]), np.array([0, 0, 0, 0, 0, 0, 1.]), np.array([1, -1, 1, 1, 1, 1, 1.3])]
nb = [kn_null(E6, p) for p in bad6]
report(f'e6^(1) (27, affine Cartan matrix verified): invertible K for uniform +1, 5 random sign patterns (classes prod C^n = {classes}) and C = 0; '
       f'none for 4 off-list C  [nullities: good {sorted(set(ng))}, bad {sorted(set(nb))}]',
       ok and classes == [-1, 1] and all(x == 1 for x in ng) and all(x <= 0 for x in nb))

# reductions: K -> K^{-1} maps C -> -C;  K -> D K D^{-1}, D = e^{i pi mu.H}, flips C_j -> (-1)^{alpha_j.mu} C_j
lam = 0.71 + 0.2j
C = pats[2]; K = kn_K(E6, C, lam)
res = lambda Kx, Cx: max(np.linalg.norm(Kx @ (-c * (E @ E.T - E.T @ E) + lam * E - E.T / lam) - (c * (E @ E.T - E.T @ E) + lam * E - E.T / lam) @ Kx)
                         for E, c in zip(E6, Cx)) / np.linalg.norm(Kx)
ok = res(K, C) < 1e-10 and res(np.linalg.inv(K), -C) < 1e-10
for k in range(6):                                  # mu = fundamental coweight k+1: alpha_j.mu = delta_jk, alpha_0.mu = -theta_k
    D = np.diag(np.exp(1j * np.pi * rootc[:, k]))
    Cf = C * np.array([(-1) ** n6[k + 1]] + [(-1) ** (j == k) for j in range(6)])
    ok &= res(D @ K @ np.linalg.inv(D), Cf) < 1e-10 and res(D @ K @ np.linalg.inv(D), C) > 1e-3
report('reductions (e6, 27): K^{-1} solves the equation for -C; D K D^{-1} with D = e^{i pi mu.H} for C_j -> (-1)^{alpha_j.mu} C_j '
       '(6 fundamental coweights); so for simply-laced rows only the classes of prod_j C_j^{n_j} need checking', ok)

# f_4^(1): Z_2 fold, nodes (0, 2, 4, {3,5}, {1,6}) = (C_0, ..., C_4), 3 and 4 short
F4 = kn_normalize([E6[0], E6[2], E6[4], E6[3] + E6[5], E6[1] + E6[6]])
ok = np.allclose(kn_cartan(F4), [[2, -1, 0, 0, 0], [-1, 2, -1, 0, 0], [0, -1, 2, -1, 0], [0, 0, -2, 2, -1], [0, 0, 0, -1, 2]])
good = [np.ones(5), np.array([1, -1, 1, -1, -1.]), np.array([-1, -1, 1, 1, -1.]), np.array([1, 1, 1, 0, 0.]), np.array([-1, 1, -1, 0, 0.]),
        np.array([1, -1, -1, 0, 0.]), np.zeros(5)]
bad = [np.array([1, 1, 1, 0.5, 1]), np.array([1, 1, 1, 1, 0.]), np.array([1, 1, 1, 0, 1.]), np.array([1, 1, 0.4, 0, 0]), np.array([0, 0, 1, 0, 0.]),
       np.array([1, 1, 1, 0, 0.6])]
ng = [kn_null(F4, c) for c in good]; nb = [kn_null(F4, c) for c in bad]
report(f'f4^(1) (27 = 26 + 1 of the folded algebra): invertible K for all +-1 and for (C_0, C_1, C_2) = +-1, C_3 = C_4 = 0; none for 6 off-list C  '
       f'[nullities: good {sorted(set(ng))}, bad {sorted(set(nb))}]', ok and all(x > 0 for x in ng) and all(x <= 0 for x in nb))

# d_4^(3): Z_3 fold, nodes ({0,1,6}, {2,3,5}, 4) = (C_0, C_1, C_2), C_2 long.  BCDR App. D also lists C_0 = C_1 = 0, C_2 free: wrong.
D43 = kn_normalize([E6[0] + E6[1] + E6[6], E6[2] + E6[3] + E6[5], E6[4]])
ok = np.allclose(kn_cartan(D43), [[2, -1, 0], [-1, 2, -3], [0, -1, 2]])
good = [np.array(p, float) for p in itertools.product([1, -1], repeat=3)] + [np.zeros(3)]
bcdr = [np.array([0, 0, c]) for c in (0.7, 2 + 1j, 1, -1)]
bad = [np.array([1, 1, 0.5]), np.array([0, 1, 1.]), np.array([1, 0, 1.]), np.array([0.5, 0, 0]), np.array([1, 1, 0.])]
info = []
for P in summands(D43):
    R = [P.T @ E @ P for E in D43]
    if all(np.allclose(r_, 0) for r_ in R):
        continue                                                            # trivial summands
    R = kn_normalize(R)
    ng = [kn_null(R, c) for c in good]; nbc = [kn_null(R, c) for c in bcdr]; nb = [kn_null(R, c) for c in bad]
    info.append(f'dim {P.shape[1]}: good {sorted(set(ng))}, C_0=C_1=0 {sorted(set(nbc))}, bad {sorted(set(nb))}')
    ok &= all(x > 0 for x in ng) and all(x == 0 for x in nbc) and all(x <= 0 for x in nb)
report('d4^(3) (faithful summands of the folded 27): invertible K for all 8 sign patterns and C = 0; NO solution for C_0 = C_1 = 0 with '
       f'C_2 = 0.7, 2+i, +1, -1 (corrects BCDR App. D) or 5 other off-list C  [{"; ".join(info)}]', ok and len(info) == 2)

# a_{2n-1}^(2), n = 3, 4: sl(2n) with sigma(X) = -J X^T J^{-1}; nodes (0, 1, ..., n), node 0 in g_1, node n the long end
for n in (3, 4):
    N = 2 * n; J = np.zeros((N, N))
    for i in range(n):
        J[i, N - 1 - i] = 1; J[N - 1 - i, i] = -1
    sig = lambda X: -J @ X.T @ np.linalg.inv(J)
    Es = [(eij(N, i, i + 1) + sig(eij(N, i, i + 1))) / 2 for i in range(n)]
    E0 = (eij(N, N - 1, 1) - sig(eij(N, N - 1, 1))) / 2
    A2 = kn_normalize([E0] + Es); Ac = kn_cartan(A2)
    ok = Ac[n - 1, n] == -2 and Ac[n, n - 1] == -1 and Ac[0, 2] == -1 and Ac[1, 2] == -1 and Ac[0, 1] == 0
    good = [np.ones(n + 1)] + [rng.choice([1., -1.], n + 1) for _ in range(4)] + [np.r_[np.zeros(n), 0.7], np.r_[np.zeros(n), 2 + 1j], np.zeros(n + 1)]
    bad = [np.r_[np.ones(n), 0.5], np.r_[0.5, np.ones(n)], np.r_[np.ones(n), 0.], np.r_[np.zeros(n - 1), 1, 1.], np.r_[1, np.zeros(n - 1), 1.]]
    ng = [kn_null(A2, c) for c in good]; nb = [kn_null(A2, c) for c in bad]
    report(f'a{2 * n - 1}^(2) ({N}-dim): invertible K for all +-1 (uniform and 4 random) and for 0 except the long end (0.7, 2+i); none for 5 off-list C  '
           f'[nullities: good {sorted(set(ng))}, bad {sorted(set(nb))}]', ok and all(x == 1 for x in ng) and all(x <= 0 for x in nb))

# ---------------------------------------------------------------- 8. Sklyanin's formalism: r-matrix and reflection equations (a_1^(1), a_2^(1))
# Defining rep of a_{N-1}^(1): E_0 = e_{N-1,0}, E_j = e_{j-1,j}; F_j(lam) = lam E_j + E_j^T/lam (exponential part of a_x).


def rm_setup(N):
    E = [eij(N, N - 1, 0)] + [eij(N, j - 1, j) for j in range(1, N)]
    Hs = []
    for k in range(1, N):
        d = np.zeros(N); d[:k] = 1; d[k] = -k; Hs.append(np.diag(d / np.sqrt(k * (k + 1))))
    return E, [e @ e.T - e.T @ e for e in E], Hs


def r_toda(N, z, bt=beta):
    """r(z) = beta^2/(4(z^h-1)) [(z^h+1) sum_a H_a(x)H_a + 2 sum_{alpha>0} (z^{ht} E_alpha(x)E_-alpha + z^{h-ht} E_-alpha(x)E_alpha)]
    (alpha^2 = 2), i.e. (beta^2/4) x Delius hep-th/9809140 (3.4)."""
    E, h_, Hs = rm_setup(N)
    r = (z**N + 1) * sum(np.kron(H, H) for H in Hs)
    for i in range(N):
        for k in range(i + 1, N):
            r = r + 2 * (z**(k - i) * np.kron(eij(N, i, k), eij(N, k, i)) + z**(N - k + i) * np.kron(eij(N, k, i), eij(N, i, k)))
    return bt**2 / (4 * (z**N - 1)) * r


def pt(X, N, which):
    """partial transpose in End(V (x) V): which = 1, 2 or 12."""
    X4 = X.reshape(N, N, N, N)
    perm = {1: (2, 1, 0, 3), 2: (0, 3, 2, 1), 12: (2, 3, 0, 1)}[which]
    return X4.transpose(perm).reshape(N * N, N * N)


def RE_plus(K, rf, l1, l2, N):           # eq-classical-reflection, sigma = transposition
    I = np.eye(N); K1, K2 = np.kron(K(l1), I), np.kron(I, K(l2))
    return K1 @ K2 @ rf(l1 / l2) + pt(rf(l2 / l1), N, 12) @ K1 @ K2 + K1 @ pt(rf(l1 * l2), N, 2) @ K2 + K2 @ pt(rf(1 / (l1 * l2)), N, 1) @ K1


def RE_minus(K, rf, l1, l2, N):          # the K_- equation of prp-sklyanin (b)
    I = np.eye(N); K1, K2 = np.kron(K(l1), I), np.kron(I, K(l2))
    return K1 @ K2 @ pt(rf(l1 / l2), N, 12) + rf(l2 / l1) @ K1 @ K2 + K1 @ pt(rf(l1 * l2), N, 1) @ K2 + K2 @ pt(rf(1 / (l1 * l2)), N, 2) @ K1


def durham_K(N, C, lam):
    E, h_, _ = rm_setup(N)
    K = kn_K(E, C, lam)
    return K / np.linalg.det(K) ** (1 / N)


def tau_bracket(N, rf, Kp, Km, l1, l2):
    """{tau(l1), tau(l2)}, tau(l) = tr(K_+(l) T(l) K_-(1/l) T(1/l)^T), by the chain rule from {T(a) (x), T(b)} = [r(a/b), T(a) (x) T(b)]
    alone, with random T's.  The book's transport matrix (d_x T = -a_x T) has the opposite sign, {T (x) T} = -[r, T (x) T]
    (prp-sklyanin; defect_lax_rmatrix.py check 3); the sign is harmless here, because the reflection equations and the
    involution are homogeneous in r."""
    lams = [l1, 1 / l1, l2, 1 / l2]
    Ts = [rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N)) for _ in range(4)]
    def grads(lam, T, Tp):
        return Kp(lam).T @ Tp @ Km(1 / lam).T, Kp(lam) @ T @ Km(1 / lam)
    g1, g2 = grads(l1, Ts[0], Ts[1]), grads(l2, Ts[2], Ts[3]); tot = 0
    for a_, ga in zip([0, 1], g1):
        for b_, gb in zip([2, 3], g2):
            r = rf(lams[a_] / lams[b_]); TT = np.kron(Ts[a_], Ts[b_])
            tot += np.einsum('ij,kl,ikjl->', ga, gb, (r @ TT - TT @ r).reshape(N, N, N, N))
    return tot


pairs = [(0.7, 1.3), (1.9, 0.45), (0.8 + 0.3j, 1.1 - 0.2j), (2.2, -0.6)]
for N in (2, 3):
    E, h_, Hs = rm_setup(N); I = np.eye(N); rf = lambda z, N=N: r_toda(N, z)
    # r is the r-matrix of the canonical brackets: [r(l1/l2), F_j(l1)(x)1 + 1(x)F_j(l2)] = (beta^2/4)(F_j(l1)(x)h_j - h_j(x)F_j(l2)),
    # [r, H(x)1 + 1(x)H] = 0  (from {phi_a(x), pi_b(y)} = delta_ab delta(x-y))
    br = 0
    for l1, l2 in pairs[:3]:
        r = rf(l1 / l2)
        for j in range(N):
            X = np.kron(l1 * E[j] + E[j].T / l1, I) + np.kron(I, l2 * E[j] + E[j].T / l2)
            tgt = (beta**2 / 4) * (np.kron(l1 * E[j] + E[j].T / l1, h_[j]) - np.kron(h_[j], l2 * E[j] + E[j].T / l2))
            br = max(br, np.linalg.norm(r @ X - X @ r - tgt))
        br = max(br, max(np.linalg.norm(r @ (np.kron(H, I) + np.kron(I, H)) - (np.kron(H, I) + np.kron(I, H)) @ r) for H in Hs))
    sym = max(np.linalg.norm(pt(rf(z), N, 12) + rf(1 / z)) for z in (0.37, 1.8, 0.5 + 0.4j))
    report(f'a{N - 1}^(1): r(z) = (beta^2/4) Delius (3.4) is the r-matrix of the canonical brackets of a_x, and r^(t1 t2)(z) = -r(1/z)  '
           f'[residuals {br:.1e}, {sym:.1e}]', br < 1e-12 and sym < 1e-12)
    Cs = [(1, 1), (0.3 + 0.2j, -1.7), (2.5, 0.4j), (0, 0), (-1, 1)] if N == 2 else [tuple(p) for p in itertools.product([1, -1], repeat=3)] + [(0, 0, 0)]
    rp = rm = xm = tr = 0
    for C in Cs:
        C = np.array(C, complex)
        Kf = lambda lam, C=C: durham_K(N, C, lam)
        Km = lambda lam, Kf=Kf: Kf(-1 / lam)
        rp = max(rp, max(np.linalg.norm(RE_plus(Kf, rf, *p, N)) for p in pairs) / np.linalg.norm(np.kron(Kf(0.7), Kf(1.3))))
        rm = max(rm, max(np.linalg.norm(RE_minus(Km, rf, *p, N)) for p in pairs) / np.linalg.norm(np.kron(Kf(0.7), Kf(1.3))))
        # K_C(-1/lam) solves the x_- condition K_-(lam) a_t(lam)^T + a_t(1/lam) K_-(lam) = 0 for the same equation d_x phi = -grad B_C
        for lam in (0.83, 1.7 + 0.3j):
            at = lambda l, j: -C[j] * h_[j] + l * E[j] - E[j].T / l
            xm = max(xm, max(np.linalg.norm(Km(lam) @ at(lam, j).T + at(1 / lam, j) @ Km(lam)) for j in range(N)) / np.linalg.norm(Km(lam)))
            A_, B_ = Kf(lam).T, durham_K(N, -C, -1 / lam)
            tr = max(tr, np.linalg.norm(A_ - np.vdot(B_.ravel(), A_.ravel()) / np.vdot(B_.ravel(), B_.ravel()) * B_) / np.linalg.norm(A_))
    lab = 'generic complex C_0, C_1' if N == 2 else 'all 8 sign patterns'
    report(f'a{N - 1}^(1), {lab} and Neumann: the Durham K solves eq-classical-reflection; K_- = K_+(-1/lam) solves the K_- equation and the '
           f'x_- condition for d_x phi = -grad B_C; K_C(lam)^T ~ K_-C(-1/lam)  [rel. residuals {rp:.1e}, {rm:.1e}, {xm:.1e}, {tr:.1e}]',
           max(rp, rm, xm, tr) < 1e-12)
    # involution from the bracket of T alone (prp-sklyanin (b)), different C at the two ends; control: random K_+
    C1, C2 = np.array(Cs[0], complex), np.array(Cs[1], complex)
    Kp = lambda lam: durham_K(N, C1, lam); Km2 = lambda lam: durham_K(N, C2, -1 / lam)
    tb = max(abs(tau_bracket(N, rf, Kp, Km2, *p)) for p in pairs[:3])
    Kr = rng.normal(size=(N, N)); tbad = abs(tau_bracket(N, rf, lambda l: Kr, Km2, 0.7, 1.3))
    report(f'a{N - 1}^(1): {{tau(l1), tau(l2)}} = 0 from {{T (x), T}} = [r, T (x) T] with K_+ = K_C, K_- = K_C\'(-1/lam)  [{tb:.1e}; random K_+: {tbad:.1e}]',
           tb < 1e-10 and tbad > 1e-3)
    # sigma = hermitian conjugation, K_+ = 1, real spectral parameters: r is real, so the reflection equation reduces to
    # r(z) + r^(t1 t2)(1/z) = 0 and r^(t2)(w) + r^(t1)(1/w) = 0
    re = max(max(np.linalg.norm(rf(l1 / l2) + pt(rf(l2 / l1), N, 12)), np.linalg.norm(pt(rf(l1 * l2), N, 2) + pt(rf(1 / (l1 * l2)), N, 1)))
             for l1, l2 in [(0.7, 1.3), (1.9, 0.45), (2.2, -0.6)])
    report(f'a{N - 1}^(1), sigma = dagger, K_+ = 1, real lambda: r real, r(z) + r^(t1 t2)(1/z) = 0 and r^(t2)(w) + r^(t1)(1/w) = 0  [{re:.1e}]',
           re < 1e-12 and np.abs(rf(0.7 / 1.3).imag).max() == 0)
