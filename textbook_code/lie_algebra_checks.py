"""Checks of Lie-algebraic statements of Part II (chapters 5-8).

1. Coxeter plane (sec-weyl-coxeter, fig-root-systems): in B_3 a short-root orbit and a long-root
   orbit of the bicoloured Coxeter element project to h-gons of the same radius; in E_8 the eight
   radii are distinct and proportional to the Perron-Frobenius eigenvector.
2. Twisted foldings (tbl-foldings, sec-affine-folding): orbit averages, coprime labels, the factor
   k, dual labels n_j^vee = k n_j alpha_j^2/2 as the null vector of the transposed Cartan matrix,
   and h, h^vee, for every twisted row, including both parents of a_{2n}^{(2)}, whose long node 0
   has n_0^vee = 2 (sec-amplitude-and-gradation).
3. The two parents of a_{2n}^{(2)}: the restricted mass matrices are the coprime-label mass matrix
   at m_F^2 = m_A^2 = 2 m_D^2, with masses^2 4 m_F^2 sin^2(pi a/(2n+1)).
4. Gradation (sec-gradations): sum_j n_j^vee s_j = 4 pi k h/beta^2 - h^vee for X^{(k)}.
   Duality of the self-dual a_{2n}^{(2)} (sec-affine-folding): rescaled coroots reproduce the roots,
   beta^vee = 8 pi/beta, B = (beta^2/4pi)/(1+beta^2/8pi) goes to 2 - B; for a_2^{(2)} the
   Bullough-Dodd S-matrix in this normalization.
5. Diagram symmetries (sec-affine-dynkin): the rotation of the a_n^{(1)} cycle multiplies each mass
   eigenvector by a phase; the reflection maps species a to n+1-a.
6. S^2 shift (sec-qg-crossing): V(x)^** = V(x t), found by solving for the intertwiner, with
   t = q^{h^vee (alpha_0,alpha_0)} in the normalization of sec-uqghat: q^4 for spin 1/2 of
   U_q(sl_2^), q^8 for the spinor of U_q(d_3^{(2)}), q^3 = q^{h^vee} for the 3-dimensional
   representation of U_q(a_2^{(2)}) with alpha_0^2 = 1 (q^6 if instead node 0 is the long root).
7. Twisted tensor-product graph (sec-tpg): the spinor R-matrix of U_q(d_3^{(2)}) from Jimbo's
   equations has eigenvalues 1, <1/2>, <1/2><1>_+ with q_i = q on the long roots, that is
   1, <1>, <1><2>_+ in the q of sec-uqghat.
"""
import itertools
from math import gcd
from functools import reduce

import numpy as np
import sympy as sp

np.set_printoptions(precision=6, suppress=True)


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


# ------------------------------------------------------------------ root systems
def cartan(roots):
    return np.array([[2 * a @ b / (a @ a) for b in roots] for a in roots])


def closure(simple):
    """All roots generated from the simple roots by simple reflections."""
    roots = [np.array(a, float) for a in simple]
    seen = {tuple(np.round(a, 9)) for a in roots}
    frontier = list(roots)
    while frontier:
        new = []
        for v in frontier:
            for a in simple:
                a = np.array(a, float)
                w = v - 2 * (a @ v) / (a @ a) * a
                key = tuple(np.round(w, 9))
                if key not in seen:
                    seen.add(key)
                    roots.append(w)
                    new.append(w)
        frontier = new
    return roots


def simple_from_cartan(A):
    """Simple roots of a simply-laced algebra (all squared lengths 2) from its Cartan matrix."""
    return list(np.linalg.cholesky(np.array(A, float)))


def affine_data(simple):
    """alpha_0 = -theta and Kac labels (n_0=1, n_i = coefficients of the highest root)."""
    S = np.array(simple)
    roots = closure(simple)
    coeffs = [np.linalg.solve(S.T, r) for r in roots]
    top = max(range(len(roots)), key=lambda i: coeffs[i].sum())
    n = np.rint(coeffs[top]).astype(int)
    return [-roots[top]] + list(S), [1] + list(n)


def dynkin_cartan(r, edges):
    A = 2 * np.eye(r)
    for i, j in edges:
        A[i, j] = A[j, i] = -1
    return A


def bicoloured_radii(simple):
    """Radii of the projections of the roots onto the Coxeter plane, per orbit, with root lengths."""
    S = [np.array(a, float) for a in simple]
    r = len(S)
    A = cartan(S)
    colour = [None] * r
    colour[0] = 0
    changed = True
    while changed:
        changed = False
        for i in range(r):
            for j in range(r):
                if A[i, j] != 0 and i != j and colour[i] is not None and colour[j] is None:
                    colour[j] = 1 - colour[i]
                    changed = True

    def refl(a):
        return np.eye(len(a)) - 2 * np.outer(a, a) / (a @ a)

    w = np.eye(len(S[0]))
    for c in (0, 1):
        for i in range(r):
            if colour[i] == c:
                w = w @ refl(S[i])
    roots = closure(S)
    ev, V = np.linalg.eig(w)
    h = len(roots) // r
    k = np.argmin(abs(ev - np.exp(2j * np.pi / h)))
    v = V[:, k]
    e1, e2 = v.real / np.linalg.norm(v.real), v.imag / np.linalg.norm(v.imag)
    left = {tuple(np.round(a, 9)) for a in roots}
    orbits = []
    for a in roots:
        if tuple(np.round(a, 9)) not in left:
            continue
        orb, b = [], a
        while tuple(np.round(b, 9)) in left:
            left.discard(tuple(np.round(b, 9)))
            orb.append(b)
            b = w @ b
        rad = [np.hypot(b @ e1, b @ e2) for b in orb]
        orbits.append((len(orb), round(float(a @ a), 9), float(np.mean(rad)), float(np.ptp(rad))))
    return h, orbits, w


# 1. Coxeter plane radii -------------------------------------------------------
e = np.eye(3)
h, orbs, _ = bicoloured_radii([e[0] - e[1], e[1] - e[2], e[2]])
print('B3: h =', h, ' orbits (size, length^2, radius):',
      [(o[0], o[1], round(o[2], 6)) for o in orbs])
short = [o[2] for o in orbs if abs(o[1] - 1) < 1e-9]
long_ = [o[2] for o in orbs if abs(o[1] - 2) < 1e-9]
report('B3: every orbit is a regular h-gon (constant radius, h roots)',
       all(o[0] == h and o[3] < 1e-12 for o in orbs))
report('B3: a short-root orbit and a long-root orbit have the same Coxeter-plane radius',
       any(abs(s - l) < 1e-12 for s in short for l in long_))

E8 = dynkin_cartan(8, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)])
S8 = simple_from_cartan(E8)
h, orbs, _ = bicoloured_radii(S8)
rad = sorted(o[2] for o in orbs)
pf = np.sort(np.abs(np.linalg.eigh(2 * np.eye(8) - E8)[1][:, 0]))
print('E8: h =', h, ' radius ratios:', np.array(rad) / rad[0])
report('E8: eight orbits of 30 roots with distinct radii',
       len(orbs) == 8 and all(o[0] == 30 for o in orbs) and min(np.diff(rad)) > 1e-3)
report('E8: radii proportional to the Perron-Frobenius eigenvector',
       np.allclose(np.array(rad) / rad[0], pf / pf[0], atol=1e-12))


# 2. twisted foldings ----------------------------------------------------------
def fold(aff, labels, tau, k):
    """Orbit averages, raw labels |O| n_O, coprime labels, squared lengths, and checks."""
    N = len(aff)
    A = cartan(aff)
    ok_aut = all(abs(A[tau[i], tau[j]] - A[i, j]) < 1e-12 for i in range(N) for j in range(N))
    orbits, left = [], set(range(N))
    while left:
        i = min(left)
        orb = [i]
        while tau[orb[-1]] != i:
            orb.append(tau[orb[-1]])
        left -= set(orb)
        orbits.append(orb)
    # the orbit average is the projection on the fixed subspace only if tau preserves
    # the inner products, which ok_aut checks; labels are constant on orbits
    ok_lab = all(len({labels[j] for j in o}) == 1 for o in orbits)
    alO = [np.mean([aff[j] for j in o], axis=0) for o in orbits]
    raw = [len(o) * labels[o[0]] for o in orbits]
    c = reduce(gcd, raw)
    prim = [x // c for x in raw]
    L = [float(a @ a) for a in alO]
    AO = cartan(alO)
    null = np.allclose(sum(p * a for p, a in zip(prim, alO)), 0)
    ndual = [k * p * l / 2 for p, l in zip(prim, L)]
    ok_dual = (np.allclose(AO.T @ np.array(ndual), 0) and np.allclose(ndual, np.rint(ndual))
               and reduce(gcd, [int(round(x)) for x in ndual]) == 1)
    return dict(ok=ok_aut and ok_lab and null and ok_dual and abs(max(L) - 2) < 1e-9, orbits=orbits,
                raw=raw, c=c, prim=prim, L=L, ndual=[int(round(x)) for x in ndual],
                h=sum(prim), hdual=int(round(sum(ndual))), alO=alO)


def d_affine(N):
    """d_N^(1): node 0 = -e1-e2, 1..N-1 = e_i - e_{i+1}, N = e_{N-1}+e_N."""
    E = np.eye(N)
    simple = [E[i] - E[i + 1] for i in range(N - 1)] + [E[N - 2] + E[N - 1]]
    aff = [-E[0] - E[1]] + simple
    lab = [1, 1] + [2] * (N - 3) + [1, 1]
    return aff, lab


def a_affine(N):
    """a_{N-1}^(1) as a cycle: alpha_j = e_j - e_{j+1 mod N}, in R^N (restricted to sum x = 0)."""
    E = np.eye(N)
    return [E[j] - E[(j + 1) % N] for j in range(N)], [1] * N


def to_hyperplane(vecs):
    """Coordinates in an orthonormal basis of the hyperplane sum x = 0."""
    N = len(vecs[0])
    Q = np.linalg.qr(np.vstack([np.ones(N), np.eye(N)[:-1]]).T)[0][:, 1:]
    return [Q.T @ v for v in vecs]


rows = []
for n in (2, 3, 4):                                    # a_{2n-1}^(2) <- d_{2n}^(1)
    N = 2 * n
    aff, lab = d_affine(N)
    tau = [N, N - 1] + [N - i for i in range(2, N - 1)] + [1, 0]
    rows.append((f'a_{2*n-1}^(2) <- d_{N}^(1)', fold(aff, lab, tau, 2), 2 * n - 1, 2 * n))
for n in (2, 3, 4):                                    # d_{n+1}^(2) <- d_{n+2}^(1)
    N = n + 2
    aff, lab = d_affine(N)
    tau = [1, 0] + list(range(2, N - 1)) + [N, N - 1]
    rows.append((f'd_{n+1}^(2) <- d_{N}^(1)', fold(aff, lab, tau, 2), n + 1, 2 * n))
E7 = dynkin_cartan(7, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (2, 6)])   # chain 0-1-2-3-4-5, leg 6 at 2
aff, lab = affine_data(simple_from_cartan(E7))
# affine node attaches to the end of the long chain opposite to the length-3 arm
A7 = cartan(aff)
ends = [i for i in range(1, 8) if sum(abs(A7[i]) > 0.5) == 2]
chain = [0]
while len(chain) < 8:
    nb = [j for j in range(8) if j not in chain and abs(A7[chain[-1], j]) > 0.5 and j != 7]
    if not nb:
        break
    chain.append(nb[0])
tau = list(range(8))
for i, j in zip(chain, reversed(chain)):
    tau[i] = j
rows.append(('e_6^(2) <- e_7^(1)', fold(aff, lab, tau, 2), 9, 12))
E6 = dynkin_cartan(6, [(0, 1), (1, 2), (2, 3), (3, 4), (2, 5)])         # centre 2, legs (1,0), (3,4), (5)
aff, lab = affine_data(simple_from_cartan(E6))
A6 = cartan(aff)
centre = 3                                                               # simple index 2 -> affine index 3
legs = []
for start in [j for j in range(7) if abs(A6[centre, j]) > 0.5 and j != centre]:
    leg = [start]
    nxt = [j for j in range(7) if abs(A6[start, j]) > 0.5 and j not in (start, centre)]
    leg += nxt
    legs.append(leg)
legs.sort(key=lambda l: 0 if 0 in l else 1)                              # the leg of the affine node first
tau = list(range(7))
for a, b in zip(legs, legs[1:] + legs[:1]):
    tau[a[0]], tau[a[1]] = b[0], b[1]
rows.append(('d_4^(3) <- e_6^(1)', fold(aff, lab, tau, 3), 4, 6))

for name, f, h, hd in rows:
    print(f'{name}: orbits {f["orbits"]}, raw labels {f["raw"]}, labels {f["prim"]}, '
          f'lengths^2 {np.round(f["L"], 4).tolist()}, dual labels {f["ndual"]}')
    k = 3 if 'd_4^(3)' in name else 2
    report(f'{name}: automorphism, coprime labels, common factor k={k}, (h, h^vee) = ({h}, {hd})',
           f['ok'] and f['c'] == k and f['h'] == h and f['hdual'] == hd)

# both parents of a_{2n}^(2)
for n in (1, 2, 3):
    affA, labA = a_affine(2 * n + 1)
    affA = to_hyperplane(affA)
    tauA = [(-j) % (2 * n + 1) for j in range(2 * n + 1)]
    fA = fold(affA, labA, tauA, 2)
    N = 2 * n + 2
    affD, labD = d_affine(N)
    # ends 0,1 at node 2 and N-1,N at node N-2; 0 -> N-1 -> 1 -> N -> 0, chain reversed
    tauD = list(range(N + 1))
    tauD[0], tauD[N - 1], tauD[1], tauD[N] = N - 1, 1, N, 0
    for i in range(2, N - 1):
        tauD[i] = N - i
    fD = fold(affD, labD, tauD, 2)
    okA = fA['ok'] and fA['c'] == 1 and fA['h'] == 2 * n + 1 and fA['hdual'] == 2 * n + 1
    okD = fD['ok'] and fD['c'] == 2 and fD['h'] == 2 * n + 1 and fD['hdual'] == 2 * n + 1
    print(f'a_{2*n}^(2) <- a_{2*n}^(1): raw {fA["raw"]}, lengths^2 {np.round(fA["L"], 4).tolist()};'
          f' <- d_{N}^(1): raw {fD["raw"]}, lengths^2 {np.round(fD["L"], 4).tolist()}')
    report(f'a_{2*n}^(2) from a_{2*n}^(1) (labels coprime) and from d_{N}^(1) (common factor 2)', okA and okD)
    report(f'a_{2*n}^(2): dual labels {fA["ndual"]}, so the long node 0 has n_0^vee = 2',
           fA['ndual'] == [2] * n + [1] and abs(fA['L'][0] - 2) < 1e-9)

    # 3. masses: restrict the parent mass matrix to the invariant subspace
    def restricted(aff, lab, tau):
        M = np.array([aff[j] - aff[tau[j]] for j in range(len(aff))])
        _, s, Vt = np.linalg.svd(M)
        B = Vt[np.sum(s > 1e-10):].T                    # orthonormal basis of the invariant subspace
        M2 = sum(l * np.outer(a, a) for l, a in zip(lab, aff))
        return np.sort(np.linalg.eigvalsh(B.T @ M2 @ B))
    mA = restricted(affA, labA, tauA)                 # in units of m_A^2
    mD = restricted(affD, labD, tauD)                 # in units of m_D^2
    # coprime-label a_{2n}^(2) mass matrix from the folded roots of the A route
    M2F = sum(p * np.outer(a, a) for p, a in zip(fA['prim'], fA['alO']))
    mF = np.sort(np.linalg.eigvalsh(M2F))
    mF = mF[mF > 1e-9]
    closed = np.sort([4 * np.sin(np.pi * a / (2 * n + 1)) ** 2 for a in range(1, n + 1)])
    print(f'  masses^2: coprime {mF}, A parent {mA}, D parent {mD}')
    report(f'a_{2*n}^(2): masses^2 4 m_F^2 sin^2(pi a/{2*n+1}); A parent at m_F = m_A; D parent at m_F^2 = 2 m_D^2',
           np.allclose(mF, closed) and np.allclose(mA, closed) and np.allclose(mD, 2 * closed))

# 4. gradation of twisted algebras ----------------------------------------------
b = sp.symbols('beta', positive=True)
ok = True
for name, f, h, hd in rows:
    k = 3 if 'd_4^(3)' in name else 2
    L = [sp.nsimplify(round(x, 9)) for x in f['L']]
    s = [8 * sp.pi / (b**2 * l) - 1 for l in L]
    nd = [sp.nsimplify(k * p * l / 2) for p, l in zip(f['prim'], L)]
    tot = sp.simplify(sum(x * y for x, y in zip(nd, s)) - (4 * sp.pi * k * h / b**2 - hd))
    ok &= tot == 0
report('twisted X^(k): sum_j n_j^vee s_j = 4 pi k h/beta^2 - h^vee, i.e. 4 pi h/beta_lit^2 - h^vee', ok)

# 4b. strong-weak duality of the self-dual a_{2n}^(2) ----------------------------
ok = True
for n in (1, 2, 3):
    affA, labA = a_affine(2 * n + 1)
    affA = to_hyperplane(affA)
    fA = fold(affA, labA, [(-j) % (2 * n + 1) for j in range(2 * n + 1)], 2)
    co = [2 * a / (a @ a) / 2 for a in fA['alO']]          # coroots, rescaled by 1/2
    G, Gco = [np.array([[a @ b for b in v] for a in v]) for v in (fA['alO'], co)]
    ok &= np.allclose(Gco, G[::-1, ::-1])                 # the same roots, chain reversed
report('a_{2n}^(2) (n=1,2,3): coroots rescaled by 1/2 are its roots with the chain reversed, '
       'so beta^vee = 2 (4 pi/beta) = 8 pi/beta', ok)
Bbook = (b**2 / (4 * sp.pi)) / (1 + b**2 / (8 * sp.pi))
blit = b / sp.sqrt(2)                                     # longest roots of squared length 4
Blit = (blit**2 / (2 * sp.pi)) / (1 + blit**2 / (4 * sp.pi))
report('a_{2n}^(2): B = (beta_lit^2/2pi)/(1+beta_lit^2/4pi) = (beta^2/4pi)/(1+beta^2/8pi), '
       'and B(8 pi/beta) = 2 - B(beta)',
       sp.simplify(Bbook - Blit) == 0 and sp.simplify(Bbook.subs(b, 8 * sp.pi / b) - (2 - Bbook)) == 0)
# a_2^(2), the Bullough-Dodd model: action (1/16pi)(d phi)^2 + mu (e^{sqrt2 b phi} + 2 e^{-b phi/sqrt2}) and
# S-matrix of O. Alekseev, arXiv:1106.4758, eqs. (3.1), (3.5); its b <-> 1/b is the duality above.
bb = sp.symbols('b', positive=True)
expo = sp.sqrt(2) * bb * sp.sqrt(8 * sp.pi)               # phi = sqrt(8 pi) x canonical field
ok_norm = sp.simplify((expo / sp.sqrt(2)).subs(bb, b / sp.sqrt(8 * sp.pi)) - b) == 0  # e^{sqrt2 beta phi}: b^2 = beta^2/8pi


def fA_(A, th):
    return np.tanh((th + 1j * np.pi * A) / 2) / np.tanh((th - 1j * np.pi * A) / 2)


def S_bd(th, bv):
    Q = bv + 1 / bv
    t = lambda a: np.tanh((th + a) / 2)
    return (t(2j * np.pi / 3) * t(-2j * np.pi / (3 * bv * Q)) * t(-2j * np.pi * bv / (3 * Q))
            / (t(-2j * np.pi / 3) * t(2j * np.pi / (3 * Q * bv)) * t(2j * np.pi * bv / (3 * Q))))


ok_bd = True
for beta_v in (np.sqrt(2.0), np.e, 5.3):
    bv = beta_v / np.sqrt(8 * np.pi)
    Bv = float(Bbook.subs(b, beta_v))
    for th in (0.37, 1.9 + 0.2j, -0.8):
        ok_bd &= abs(S_bd(th, bv) - fA_(2 / 3, th) * fA_(-Bv / 3, th) * fA_((Bv - 2) / 3, th)) < 1e-12
        ok_bd &= abs(S_bd(th, bv) - S_bd(th, 1 / bv)) < 1e-12
report('a_2^(2) (Bullough-Dodd): b^2 = beta^2/8pi; S = f_{2/3} f_{-B/3} f_{(B-2)/3} with the B above; '
       'invariant under b -> 1/b, i.e. beta -> 8 pi/beta', ok_norm and ok_bd)

# 5. diagram symmetries of a_n^(1) acting on the particles ------------------------
ok_rot = ok_ref = True
for N in range(3, 9):
    aff, lab = a_affine(N)
    # the mass eigenvector of species a has components exp(2 pi i a j/N) on the roots alpha_j
    for a in range(1, N):
        v = np.exp(2j * np.pi * a * np.arange(N) / N)
        rot = np.roll(v, -1)                           # rotation j -> j+1 of the affine diagram
        ref = v[(-np.arange(N)) % N]                   # reflection j -> -j
        ok_rot &= np.allclose(rot, np.exp(2j * np.pi * a / N) * v)
        ok_ref &= np.allclose(ref, np.exp(2j * np.pi * (N - a) * np.arange(N) / N))
report('a_n^(1): rotation of the cycle acts on species a by the phase e^{2 pi i a/(n+1)}', ok_rot)
report('a_n^(1): reflection of the cycle maps species a to species n+1-a', ok_ref)


# 6. S^2 shift -------------------------------------------------------------------
def qint(n, q):
    return (q**n - q**-n) / (q - q**-1)


def qbinom(n, s, q):
    num = np.prod([qint(n - i, q) for i in range(s)]) if s else 1
    den = np.prod([qint(i + 1, q) for i in range(s)]) if s else 1
    return num / den


def check_relations(rep, A, qi):
    e, f, k = rep
    r = len(e)
    err = 0
    for i in range(r):
        for j in range(r):
            err = max(err, np.abs(k[i] @ e[j] @ np.linalg.inv(k[i]) - qi[i] ** A[i][j] * e[j]).max())
            err = max(err, np.abs(k[i] @ f[j] @ np.linalg.inv(k[i]) - qi[i] ** -A[i][j] * f[j]).max())
            rhs = (k[i] - np.linalg.inv(k[i])) / (qi[i] - 1 / qi[i]) if i == j else 0 * k[i]
            err = max(err, np.abs(e[i] @ f[j] - f[j] @ e[i] - rhs).max())
            if i != j:
                m = int(round(1 - A[i][j]))
                for g in (e, f):
                    ser = sum((-1) ** s * qbinom(m, s, qi[i]) * np.linalg.matrix_power(g[i], m - s)
                              @ g[j] @ np.linalg.matrix_power(g[i], s) for s in range(m + 1))
                    err = max(err, np.abs(ser).max())
    return err


def at(rep, x):
    """Homogeneous gradation: e_0 -> x e_0, f_0 -> f_0/x."""
    e, f, k = rep
    return ([x * e[0]] + e[1:], [f[0] / x] + f[1:], k)


def dual(rep):
    """pi^*(a) = pi(S(a))^T with S(e) = -e k^{-1}, S(f) = -k f, S(k) = k^{-1} (eq-antipode)."""
    e, f, k = rep
    ki = [np.linalg.inv(m) for m in k]
    return ([(-a @ b).T for a, b in zip(e, ki)], [(-b @ a).T for a, b in zip(f, k)], [m.T for m in ki])


def shifts(rep, x, cands, target=None):
    """Values t for which some M != 0 satisfies M rho(g) = pi_{x t}(g) M for all generators."""
    found = []
    d = rep[0][0].shape[0]
    for t in cands:
        tgt = at(target if target is not None else base, x * t)
        rows_ = []
        for g1, g2 in zip(sum(rep, []), sum(tgt, [])):
            rows_.append(np.kron(np.eye(d), g1.T) - np.kron(g2, np.eye(d)))
        s = np.linalg.svd(np.vstack(rows_), compute_uv=False)
        if s[-1] < 1e-9 * s[0]:
            found.append(t)
    return found


def shift_test(name, rep0, A, qi, p, expect_dd, expect_d=None):
    global base
    base = rep0
    err = check_relations(at(rep0, 0.7 + 0.2j), A, qi)
    x = 0.83 + 0.31j
    V = at(rep0, x)
    cands = [sgn * p**m for m in range(-30, 31) for sgn in (1, -1, 1j, -1j)]
    dd = shifts(dual(dual(V)), x, cands)
    d = shifts(dual(V), x, cands)
    ok = err < 1e-10 and len(dd) == 1 and abs(dd[0] - expect_dd) < 1e-9
    if expect_d is not None:
        ok &= all(abs(t**2 - expect_dd) < 1e-9 for t in d) and len(d) == 1 and abs(d[0] - expect_d) < 1e-9
    print(f'{name}: relations {err:.1e}; V** = V(x t) for t = {np.round(dd, 6)}; '
          f'V* = V(x x_c) for x_c = {np.round(d, 6)}')
    return ok, dd, d


q = 1.1 * np.exp(1j * np.sqrt(2))                    # generic, not a root of unity

# U_q(sl_2^), spin 1/2, eq-evaluation-rep: alpha_0^2 = 2, h^vee = 2
E_ = np.array([[0, 1], [0, 0]], complex)
F_ = E_.T.copy()
K_ = np.diag([q, 1 / q])
rep = ([F_, E_], [E_, F_], [np.linalg.inv(K_), K_])
A2 = [[2, -2], [-2, 2]]
report('U_q(sl_2^) spin 1/2: V** = V(x q^4), V* = V(x q^2)',
       shift_test('sl2 spin 1/2', rep, A2, [q, q], q, q**4, q**2)[0])

# U_q(d_3^(2)) spinor (chapter 23 realization, sec-uqghat normalization: q_i = q on short roots)
n = 2
sp_, sm_, I2 = np.array([[0, 1], [0, 0]], complex), np.array([[0, 0], [1, 0]], complex), np.eye(2)


def site(op, i):
    mats = [I2] * n
    mats[i] = op
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


Ws = np.array([[0.5 if bit == 0 else -0.5 for bit in bits] for bits in itertools.product([0, 1], repeat=n)])
e = [site(sm_, 0)] + [site(sp_, i - 1) @ site(sm_, i) for i in range(1, n)] + [site(sp_, n - 1)]
f = [m.T.copy() for m in e]
k = ([np.diag(q ** (-2 * Ws[:, 0])).astype(complex)]
     + [np.diag(q ** (2 * (Ws[:, i - 1] - Ws[:, i]))).astype(complex) for i in range(1, n)]
     + [np.diag(q ** (2 * Ws[:, n - 1])).astype(complex)])
qd = [q] + [q**2] * (n - 1) + [q]
Ad = [[2, -2, 0], [-1, 2, -1], [0, -2, 2]]
spinor = (e, f, k)
report('U_q(d_3^(2)) spinor: V** = V(x q^{2h^vee}) = V(x q^8), V* = V(-x q^4)',
       shift_test('d3^(2) spinor', spinor, Ad, qd, q, q**8, -q**4)[0])

# U_q(a_2^(2)), roots of squared lengths 1 (short, q_s = q^{1/2}) and 4 (long, q_l = q^2), h^vee = 3
p = np.sqrt(q)                                       # q_s
qs, ql = p, p**4
# (i) alpha_0 short, alpha_0^2 = 1 (sec-uqghat): horizontal U_{q_l}(sl_2) on 2 + 1, basis u+, u-, w
e1 = np.zeros((3, 3), complex); e1[0, 1] = 1
e0 = np.zeros((3, 3), complex); e0[2, 0] = 1; e0[1, 2] = 1
f0 = np.zeros((3, 3), complex); f0[0, 2] = qint(2, qs); f0[2, 1] = qint(2, qs)
k1 = np.diag([ql, 1 / ql, 1]); k0 = np.diag([1 / qs**2, qs**2, 1])
repK = ([e0, e1], [f0, e1.T.copy()], [k0, k1])
AK = [[2, -4], [-1, 2]]
_, ddK, dK = shift_test('a2^(2), alpha_0^2 = 1', repK, AK, [qs, ql], p, q**3, None)
base = repK
pm = shifts(at(repK, -0.83 - 0.31j), -0.83 - 0.31j, [-1])     # is V(x) = V(-x)?
okK = (check_relations(at(repK, 0.7), AK, [qs, ql]) < 1e-10
       and any(abs(t - q**3) < 1e-9 for t in ddK) and all(abs(t - q**6) > 1e-6 for t in ddK)
       and all(min(abs(t - q**3), abs(t + q**3)) < 1e-9 for t in ddK) and pm == [-1]
       and all(any(abs(c**2 - t) < 1e-9 for t in ddK) for c in dK) and len(dK) > 0)
print(f'  a2^(2), alpha_0^2 = 1: V(x) = V(-x): {pm == [-1]}; x_c^2 = {np.round([c**2 for c in dK], 6)}')
# (ii) alpha_0 long, alpha_0^2 = 4: horizontal U_{q_s}(sl_2) on spin 1, e_0 ~ F^2
Es = np.zeros((3, 3), complex); Es[0, 1] = qint(1, qs) * qint(2, qs); Es[1, 2] = qint(2, qs) * qint(1, qs)
Fs = np.zeros((3, 3), complex); Fs[1, 0] = 1; Fs[2, 1] = 1
Ks = np.diag([qs**2, 1, qs**-2])
k0 = np.linalg.matrix_power(np.linalg.inv(Ks), 2)
c0 = ((k0 - np.linalg.inv(k0)) / (ql - 1 / ql))[2, 2] / (Fs @ Fs @ Es @ Es - Es @ Es @ Fs @ Fs)[2, 2]
repB = ([Fs @ Fs, Es], [c0 * Es @ Es, Fs], [k0, Ks])
AB = [[2, -1], [-4, 2]]
okB = shift_test('a2^(2), alpha_0^2 = 4', repB, AB, [ql, qs], p, q**6, -q**3)[0]
report('U_q(a_2^(2)) with alpha_0^2 = 1: V** = V(x q^{h^vee}) = V(x q^3), not V(x q^{2h^vee}); '
       'here V(x) = V(-x), so also V(-x q^3), and x_c^2 = +-q^3', okK)
report('U_q(a_2^(2)) with node 0 the long root (alpha_0^2 = 4): V** = V(x q^{2h^vee}) = V(x q^6), V* = V(-x q^3)', okB)

# 7. twisted tensor-product graph: spinor of U_q(d_3^(2)) --------------------------
def coproduct(r1, r2):
    """eq-coproduct: Delta(e) = e x k + 1 x e, Delta(f) = f x 1 + k^{-1} x f, Delta(k) = k x k."""
    e1_, f1_, k1_ = r1
    e2_, f2_, k2_ = r2
    I = np.eye(e1_[0].shape[0])
    return ([np.kron(a, kk) + np.kron(I, c) for a, kk, c in zip(e1_, k2_, e2_)],
            [np.kron(a, I) + np.kron(np.linalg.inv(kk), c) for a, kk, c in zip(f1_, k1_, f2_)],
            [np.kron(a, c) for a, c in zip(k1_, k2_)])


def rcheck(rep, x, y):
    """Solve R Delta_{x,y}(g) = Delta_{y,x}(g) R for R: V(x) x V(y) -> V(y) x V(x)."""
    D1 = coproduct(at(rep, x), at(rep, y))
    D2 = coproduct(at(rep, y), at(rep, x))
    d = D1[0][0].shape[0]
    rows_ = [np.kron(np.eye(d), g1.T) - np.kron(g2, np.eye(d)) for g1, g2 in zip(sum(D1, []), sum(D2, []))]
    _, s, Vt = np.linalg.svd(np.vstack(rows_))
    return Vt[-1].conj().reshape(d, d), s[-2] / s[0], s[-1] / s[0]


x, y = 1.3 + 0.4j, 0.9 - 0.2j
z = x / y
R, s2, s1 = rcheck(spinor, x, y)
R = R / R[0, 0]                                      # 1 on the highest-weight vector
ev = np.linalg.eigvals(R)
groups = []
for v in ev:
    for g in groups:
        if abs(g[0] - v) < 1e-7:
            g.append(v)
            break
    else:
        groups.append([v])
found = sorted([(len(g), np.mean(g)) for g in groups], key=lambda t: -t[0])


def br(a, sgn, qq):
    return (1 + sgn * z * qq ** (2 * a)) / (z + sgn * qq ** (2 * a))


book_q = [1, br(0.5, -1, q**2), br(0.5, -1, q**2) * br(1, +1, q**2)]   # q_i = q on long roots: q -> q^2
kac_q = [1, br(1, -1, q), br(1, -1, q) * br(2, +1, q)]
half = [1, br(0.5, -1, q), br(0.5, -1, q) * br(1, +1, q)]
mult = [g[0] for g in found]
vals = [g[1] for g in found]
print(f'd3^(2) spinor: solution space 1-dimensional (singular values {s2:.2e}, {s1:.1e}); '
      f'eigenvalue multiplicities {mult}')
report('d_3^(2) spinor R-matrix: multiplicities 10, 5, 1 (Lambda^2, Lambda^1, Lambda^0 of C^5)',
       mult == [10, 5, 1] and s1 < 1e-12 and s2 > 1e-6)
report('d_3^(2) spinor: rho = 1, <1/2>, <1/2><1>_+ with q_i = q on long roots = 1, <1>, <1><2>_+ in the q of sec-uqghat',
       np.allclose(vals, book_q, atol=1e-8) and np.allclose(vals, kac_q, atol=1e-8))
report('d_3^(2) spinor: the factors <1/2>, <1/2><1>_+ do not hold in the q of sec-uqghat',
       not np.allclose(vals, half, atol=1e-6))
