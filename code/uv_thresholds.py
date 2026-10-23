"""First collapse threshold of sec-collapse-thresholds (chapter 46) for the non-simply-laced affine algebras.

beta^2_UV / pi = min_S 8(|S|-1)/D(S), D(S) = sum_k |alpha_k|^2 - |sum_k alpha_k|^2 > 0, over multisets S of the
simple affine roots projected to the finite space, longest roots of squared length 2.  Clusters are enumerated up to
size K = max(12, 2 + size of the full neutral cluster).  The projected roots are, with e_i
orthonormal (and rescaled so that the longest roots have norm 2):

  b_r^(1)       alpha_0 = -(e1+e2), e_i - e_(i+1), e_r                 labels 1,1,2,...,2
  c_r^(1)       alpha_0 = -2 e1,    e_i - e_(i+1), 2 e_r               labels 1,2,...,2,1
  g_2^(1), f_4^(1)  alpha_0 = -highest root                            labels 1,3,2 and 1,2,3,4,2
  a_2r^(2)      alpha_0 = -e1,      e_i - e_(i+1), 2 e_r   (BC_r)      labels 2,...,2,1
  a_(2r-1)^(2)  alpha_0 = -(e1+e2), e_i - e_(i+1), 2 e_r   (C_r)       labels 1,1,2,...,2,1
  d_(r+1)^(2)   alpha_0 = -e1,      e_i - e_(i+1), e_r     (B_r)       labels 1,...,1
  e_6^(2), d_4^(3)  alpha_0 = -highest short root of F_4, G_2          labels 1,1,2,3,2 and 1,2,1

Checks (exit status non-zero if one fails): the labels annihilate the projected roots; beta^2_UV = 4 pi for every
algebra listed (ranks up to 5), attained by a longest root and a neighbour; for b_r^(1) the full neutral cluster also
gives 4 pi; and, as in the proof in the text, every root that is not longest has |alpha|^2 <= 1, the longest roots
have integer inner products with all simple roots.  Runtime: a few seconds.
"""
import itertools
import sys

import numpy as np


def e(n, i):
    v = np.zeros(n)
    v[i] = 1
    return v


def classical(kind, r):
    """Projected simple affine roots (alpha_0 first) and labels, before normalization."""
    chain = [e(r, i) - e(r, i + 1) for i in range(r - 1)]
    if kind == "b1":
        return [-(e(r, 0) + e(r, 1))] + chain + [e(r, r - 1)], [1, 1] + [2] * (r - 1)
    if kind == "c1":
        return [-2 * e(r, 0)] + chain + [2 * e(r, r - 1)], [1] + [2] * (r - 1) + [1]
    if kind == "a2r_2":
        return [-e(r, 0)] + chain + [2 * e(r, r - 1)], [2] * r + [1]
    if kind == "a2r-1_2":
        return [-(e(r, 0) + e(r, 1))] + chain + [2 * e(r, r - 1)], [1, 1] + [2] * (r - 2) + [1]
    if kind == "dr+1_2":
        return [-e(r, 0)] + chain + [e(r, r - 1)], [1] * (r + 1)
    a1, a2 = np.array([1.0, -1, 0]), np.array([-2.0, 1, 1])                                # G_2: short, long
    f4 = [e(4, 1) - e(4, 2), e(4, 2) - e(4, 3), e(4, 3), (e(4, 0) - e(4, 1) - e(4, 2) - e(4, 3)) / 2]
    if kind == "g2_1":
        return [-(3 * a1 + 2 * a2), a1, a2], [1, 3, 2]
    if kind == "d4_3":
        return [-(2 * a1 + a2), a1, a2], [1, 2, 1]
    if kind == "f4_1":
        return [-(e(4, 0) + e(4, 1))] + f4, [1, 2, 3, 4, 2]
    if kind == "e6_2":
        return [-e(4, 0)] + f4, [1, 1, 2, 3, 2]


def threshold(G, K):
    """min over multisets of size 2..K of 8(|S|-1)/D(S), and a minimizing multiplicity vector."""
    best, arg = np.inf, None
    for size in range(2, K + 1):
        for m in itertools.combinations_with_replacement(range(len(G)), size):
            v = np.bincount(m, minlength=len(G))
            D = np.diag(G) @ v - v @ G @ v
            if D > 1e-9 and 8 * (size - 1) / D < best - 1e-9:
                best, arg = 8 * (size - 1) / D, v
    return best, arg


ok = True


def check(name, cond):
    global ok
    ok &= bool(cond)
    print(("PASS " if cond else "FAIL ") + name)


cases = [("b1", r) for r in (3, 4, 5)] + [("c1", r) for r in (2, 3, 4, 5)] + [("g2_1", 2), ("f4_1", 4)] \
    + [("a2r_2", r) for r in (1, 2, 3, 4, 5)] + [("a2r-1_2", r) for r in (3, 4, 5)] \
    + [("dr+1_2", r) for r in (2, 3, 4, 5)] + [("e6_2", 4), ("d4_3", 2)]
for kind, r in cases:
    roots, labels = classical(kind, r)
    roots = np.array(roots, dtype=float)
    G = roots @ roots.T
    G *= 2 / G.diagonal().max()
    name = f"{kind} rank {r}"
    check(f"{name}: labels annihilate the projected roots", np.allclose(np.array(labels) @ roots, 0))
    longest = np.isclose(G.diagonal(), 2)
    check(f"{name}: shorter roots have norm <= 1, longest roots integral against all", np.all(G.diagonal()[~longest] <= 1 + 1e-12)
          and np.allclose(G[longest], np.round(G[longest])))
    K = max(12, sum(labels) + 2)
    t, v = threshold(G, K)
    pair = [(i, j) for i in np.flatnonzero(longest) for j in range(len(G)) if j != i and np.isclose(G[i, j], -1)]
    print(f"   {name}: beta^2_UV = {t:.6f} pi (clusters up to {K}), attained by multiplicities {tuple(v)}")
    check(f"{name}: beta^2_UV = 4 pi, attained by a longest root and a neighbour", np.isclose(t, 4) and len(pair) > 0)
    if kind == "b1":
        n = np.array(labels)
        check(f"{name}: the full neutral cluster also gives 4 pi",
              np.isclose(8 * (n.sum() - 1) / (G.diagonal() @ n), 4))
sys.exit(0 if ok else 1)
