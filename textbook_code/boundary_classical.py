"""Classical boundary conditions of affine Toda theory on the half-line x <= 0 (Part VII, ch. 30).

Book conventions: S = int dt [ int_{-inf}^0 L dx - B(phi(0,t)) ], so d_x phi|_0 = -grad B, with
  B_real = (2m/beta^2) sum_j C_j sqrt(2 n_j/alpha_j^2) (e^{beta alpha_j.phi/2} - 1),  B_imag = B_real(beta -> i beta).
Lax pair eq-lax-pair, a_t = (A_+ + A_-)/2, a_x = (A_+ - A_-)/2.

1. Sklyanin/BCDR compatibility K a_t(lam) + a_t(1/lam)^T K = 0 at x=0 (hep-th/9809140 (2.13), (4.1)) with
   d_x phi = -grad B_real substituted, solved for a field-independent K in a faithful representation.
   An invertible K exists exactly for the classified C (tbl-durham), for a_1^(1), a_2^(1), a_3^(1),
   d_4^(1), c_2^(1) (= b_2^(1)), c_3^(1), b_3^(1), g_2^(1), a_2^(2), a_4^(2), a_6^(2), d_3^(2) (folded reps);
   none exists for the normalizations without n_j or for the literal reading |A_i| = sqrt(2 n_i) of BCDR (1.4).
   Corrections to BCDR App. D: c_2^(1) = b_2^(1) also allows (C_0, C_2 = +-1, C_1 free); d_3^(2) = a_3^(2) obeys
   the union of the d_n^(2) and a_{2n-1}^(2) rules; g_2^(1) does NOT allow C_2 (short) free -- only C = 0 or all +-1.
2. Soliton-preserving condition (sigma = hermitian conjugation, K = 1): compatibility <=> Dirichlet on
   Re phi with values in (2 pi/beta) x coweight lattice, Neumann on Im phi (imaginary coupling).
3. Theta (eq-theta-data) on the half-line: B_imag(-conj phi) = conj B_imag(phi) for all phi iff C_j real.
4. Dictionary: Delius 98a (eps = -C), Delius 98b (A_i = -C_i), Delius-Gandenberger (C = C),
   Corrigan-Delius (eps_0, eps_1) = (C_0, C_1), Ghoshal-Zamolodchikov (M_0, phi_0), Saleur-Skorik-Warner
   (zeta, eta), Baseilhac-Delius-George oscillator (C_0, C_1) = (cos p^, cos q^).
5. Boundary oscillator: energy of field + oscillator is conserved (symbolic).
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
lit = lambda n, a: np.sqrt(2 * n / (a @ a)) / np.sqrt(2)    # literal |A_i| = sqrt(2 n_i) reading of BCDR (1.4)
ok = all(K_nullity('d4(1)', C, norm=no_n) == 0 for C in signs(5)[:8])
ok &= K_nullity('d4(1)', np.ones(5), norm=lit) == 0 and K_nullity('a2(1)', np.ones(3), norm=lit) == 0
report('controls: without the factor sqrt(n_j), or with |A_i| = sqrt(2 n_i), no K exists for |C_j| = 1 (d4, a2)', ok)

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

# a numerical illustration of the 'iff': pick Re phi = (2 pi/beta) * coweight, check vanishing; shift by half -> fails
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
report('soliton-preserving: Re phi(0) = (2pi/beta) w passes, (pi/beta) w (half a coweight) fails', ok)

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
