"""Theta-breaking by colour content for N fundamental solitons of a_{n-1}^(1) (prp-krein-fundamental, chapter 48).

T_1 = R_1N(x_1/x_N) ... R_12(x_1/x_2) on (C^n)^N, with Jimbo's R normalized as in krein_bethe.py, q = exp(i mu),
0 < mu < pi, and random real rapidities (x_k = exp(1.5 g_k), g_k standard normal).  T_1 conserves the colour content,
so the sectors are the weight spaces.  They are grouped by their pattern of colour multiplicities m_c, e.g. (2,1,1)
for contents like (0,0,1,2), and for N >= 5 one sector of each pattern is computed.  The deviation of a sector is
max ||s|-1| over the eigenvalues s of T_1 on it, maximized over the rapidity samples.  Complex conjugation of q
conjugates T_1, so mu in (0, pi) covers every |q| = 1.

At q = -1 the normalized R is +1 on equal colours and -1 on different ones, so T_1 is diagonal, with eigenvalue
(-1)^(N - m_c) when soliton 1 has colour c.  Rule observed below: near q = -1 a sector is unbroken exactly when each of
these two eigenspaces involves at most two colours, i.e. at most two colours have even multiplicity and at most two
have odd multiplicity.

Checks the claims of the Numerics paragraph after prp-krein-fundamental (exit status non-zero if one fails):
 1. sectors with at most two colours are unbroken (part 3 of the proposition), in every case below;
 2. N = 3, n = 3 and 4: the sector with three colours is broken at every mu sampled in [0.3, 3] (deviation 0.08-0.49);
 3. N = 4, n = 3 and 4: the sectors with three colours are broken for mu < 2 pi/3 and unbroken for 2 pi/3 < mu < pi,
    i.e. for Re q < -1/2; for n = 4 the sector with four colours is broken at every mu sampled;
 4. the rule above at mu = 2.9 and 3.05 for every pattern with at least three colours, n = 3 with N <= 7 and n = 4
    with N <= 6 (for N = 6, n = 4, the patterns (2,2,1,1) and (3,1,1,1));
 5. N = 5, 6: the patterns (3,1,1), (2,2,2) and (2,1,1,1) are broken at every mu sampled;
 6. the onset depends on the pattern: (2,2,1) is unbroken from mu = 2.1 on, (4,1,1) only from between 2.45 and 2.55.
Runtime: about two minutes.
"""
import contextlib
import io
import itertools
import sys

import numpy as np

with contextlib.redirect_stdout(io.StringIO()):  # krein_bethe prints its own checks on import
    from krein_bethe import R

MUS = [0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.0, 2.09, 2.1, 2.2, 2.5, 2.8, 3.0]
UNBROKEN, BROKEN = 1e-8, 1e-3


def apply2(T, M, i, j, n, N):
    """Left-multiply T by the two-site matrix M (indexed (a_i a_j, b_i b_j)) acting on factors i, j."""
    T = T.reshape((n,) * N + (-1,))
    T = np.moveaxis(T, [i, j], [0, 1])
    T = np.tensordot(M.reshape(n, n, n, n), T, axes=([2, 3], [0, 1]))
    return np.moveaxis(T, [0, 1], [i, j]).reshape(n**N, -1)


def sector_basis(n, N):
    sectors = {}
    for idx in itertools.product(range(n), repeat=N):
        sectors.setdefault(tuple(sorted(idx)), []).append(np.ravel_multi_index(idx, [n] * N))
    return sectors


def deviation(n, N, ids, x, q):
    """max ||s|-1| over the eigenvalues of T_1 on the sector spanned by the basis vectors ids."""
    T = np.zeros((n**N, len(ids)), dtype=complex)
    T[ids, range(len(ids))] = 1
    for j in range(1, N):
        T = apply2(T, R(n, x[0] / x[j], q), 0, j, n, N)
    return np.abs(np.abs(np.linalg.eigvals(T[ids, :])) - 1).max()


def pattern(content):
    return tuple(sorted((content.count(c) for c in set(content)), reverse=True))


def scan(n, N, samples, mus=MUS, seed=1, only=None):
    """Largest deviation per colour pattern; only = list of patterns to compute (one sector each), default all."""
    sectors = sector_basis(n, N)
    if only is not None:
        sectors = {k: v for k, v in sectors.items() if pattern(k) in only}
        sectors = {pattern(k): (k, v) for k, v in sorted(sectors.items())}
        sectors = {k: v for k, v in sectors.values()}
    rng = np.random.default_rng(seed)
    table = {}
    for mu in mus:
        q = np.exp(1j * mu)
        worst = {}
        for _ in range(samples):
            x = np.exp(1.5 * rng.normal(size=N))
            for key, ids in sectors.items():
                worst[pattern(key)] = max(worst.get(pattern(key), 0.0), deviation(n, N, ids, x, q))
        table[mu] = worst
    pats = sorted({p for w in table.values() for p in w}, key=lambda p: (len(p), p))
    print(f"\nn={n}, N={N} ({samples} rapidity samples per mu): max ||s|-1| by colour multiplicities")
    print("  mu     Re q   " + "".join(f"{str(p):>14}" for p in pats))
    for mu, w in table.items():
        print(f"  {mu:5.2f} {np.cos(mu):+.3f} " + "".join(f"{w[p]:14.1e}" for p in pats))
    return table


ok = True


def check(name, cond):
    global ok
    ok &= bool(cond)
    print(("PASS " if cond else "FAIL ") + name)


def rule(p, N):
    """Predicted unbroken near q = -1: each parity class of N - m_c has at most two colours."""
    odd = sum((N - m) % 2 for m in p)
    return max(odd, len(p) - odd) <= 2


THIRD = 2 * np.pi / 3
NEAR = [2.9, 3.05]
t33 = scan(3, 3, 40)
t43 = scan(4, 3, 40)
t34 = scan(3, 4, 40)
t44 = scan(4, 4, 15)
t35 = scan(3, 5, 10, mus=sorted(MUS + [2.45, 2.55, 3.05]))
t36 = scan(3, 6, 6, mus=sorted(MUS + [2.45, 2.55, 3.05]))
t37 = scan(3, 7, 3, mus=NEAR)
t45 = scan(4, 5, 4, mus=sorted(MUS + [3.05]), only=[(2, 1, 1, 1), (2, 2, 1), (3, 1, 1)])
t46 = scan(4, 6, 3, mus=[2.0] + NEAR, only=[(2, 2, 1, 1), (3, 1, 1, 1)])

print()
runs = {(3, 3): t33, (4, 3): t43, (3, 4): t34, (4, 4): t44, (3, 5): t35, (3, 6): t36, (3, 7): t37, (4, 5): t45,
        (4, 6): t46}
for (n, N), t in runs.items():
    few = [d for w in t.values() for p, d in w.items() if len(p) <= 2]
    if few:
        check(f"n={n} N={N}: sectors with at most two colours unbroken (max dev {max(few):.1e})", max(few) < UNBROKEN)
d33 = [t33[mu][(1, 1, 1)] for mu in MUS]
print(f"n=3 N=3: three-colour deviation from {min(d33):.2f} to {max(d33):.2f} for 0.3 <= mu <= 3")
check("N=3, n=3 and 4: three-colour sector broken at every mu sampled",
      min(d33) > BROKEN and all(t43[mu][(1, 1, 1)] > BROKEN for mu in MUS))
check("N=4, n=3 and 4: three-colour sectors broken for mu < 2pi/3 and unbroken above",
      all((t[mu][(2, 1, 1)] < UNBROKEN) == (mu > THIRD) for t in (t34, t44) for mu in MUS)
      and all(t[mu][(2, 1, 1)] > BROKEN for t in (t34, t44) for mu in MUS if mu < THIRD))
check("N=4, n=4: four-colour sector broken at every mu sampled", all(t44[mu][(1, 1, 1, 1)] > BROKEN for mu in MUS))
for (n, N), t in runs.items():
    for p in sorted({p for mu in NEAR if mu in t for p in t[mu] if len(p) >= 3}):
        devs = [t[mu][p] for mu in NEAR if mu in t]
        want = rule(p, N)
        check(f"n={n} N={N} {p}: near q=-1 {'unbroken' if want else 'broken'} as the rule predicts "
              f"({', '.join(f'{d:.1e}' for d in devs)})",
              all(d < UNBROKEN for d in devs) if want else all(d > BROKEN for d in devs))
check("N=5, 6: (3,1,1), (2,2,2) and (2,1,1,1) broken at every mu sampled",
      all(d > BROKEN for t, p in ((t35, (3, 1, 1)), (t36, (2, 2, 2)), (t45, (2, 1, 1, 1))) for d in
          (w[p] for w in t.values())))
check("onsets: (2,2,1) unbroken from mu = 2.1 on and broken at 2.0; (4,1,1) broken at 2.45, unbroken from 2.55 on",
      t35[2.0][(2, 2, 1)] > BROKEN and all(t35[mu][(2, 2, 1)] < UNBROKEN for mu in t35 if mu >= 2.1)
      and t36[2.45][(4, 1, 1)] > BROKEN and all(t36[mu][(4, 1, 1)] < UNBROKEN for mu in t36 if mu >= 2.55))
sys.exit(0 if ok else 1)
