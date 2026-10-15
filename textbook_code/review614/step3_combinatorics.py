"""Checks of Step 3 of the proof of thm-nn-distances (chapter 47): the cluster bound on products of m-th-neighbour
distances.  One PASS/FAIL line per check; exit status non-zero if one fails.

1. Kruskal dendrograms and the hook-length sum (Step 3.3 and 3.5).  For random labelled trees tau with
   n = 3, ..., 8 vertices (5 trees for each n, fixed seed), every ordering of the n-1 edges is enumerated and
   Kruskal's algorithm is run on it.  Checks: the orderings giving a dendrogram T are exactly the linear
   extensions of T's tree order; their number is (n-1)!/prod_j h_j (hook-length formula); and
   sum_{T ~ tau} prod_j 1/h_j = 1 exactly (rational arithmetic).
2. The nested-scale integral of Step 3.4.  int_{R(T)} prod_j s_j^{1 - gamma c_j} ds = prod_j 1/E_j,
   E_j = 2 h_j - gamma C_j, over R(T) = {s_j <= s_parent(j), s <= 1}, computed exactly with sympy:
   on the 4-leaf caterpillar dendrogram with m = 2, gamma = 6/5 (both sides 25/24), on the balanced 4-leaf
   dendrogram (5/24), and on random dendrograms with up to 6 leaves, m = 1, 2, 3 and rational gamma < 2m/(m+1).
3. The bound tilde delta_k >= s_{i(k)} of Step 3.2.  For random point sets in the plane (n = 4..12, m = 1, 2, 3,
   fixed seed), with the Euclidean minimum spanning tree, its single-linkage dendrogram, and i(k) the lowest
   ancestor of k with at least m+1 leaves, the distance from u_k to its m-th nearest other point is at least
   the length of edge i(k).
Runtime: about 1 min (2 min on a loaded machine)."""
import itertools
import math
import random
import sys
from fractions import Fraction

import numpy as np
import sympy as sp

ok = True


def check(name, cond):
    global ok
    ok &= bool(cond)
    print(("PASS " if cond else "FAIL ") + name)


def random_tree(n, rnd):
    """Random labelled tree on n vertices from a random Pruefer sequence; returns the edge list."""
    if n == 2:
        return [(0, 1)]
    seq = [rnd.randrange(n) for _ in range(n - 2)]
    degree = [1] * n
    for v in seq:
        degree[v] += 1
    edges = []
    for v in seq:
        leaf = min(i for i in range(n) if degree[i] == 1)
        edges.append((leaf, v))
        degree[leaf] -= 1
        degree[v] -= 1
    u, w = [i for i in range(n) if degree[i] == 1]
    edges.append((u, w))
    return edges


def kruskal_dendrogram(n, edges, order):
    """Single-linkage dendrogram: internal node = edge; parent[j] = edge whose merger absorbs the cluster of j.
    Returns a canonical tuple of parents (None for the root)."""
    comp = list(range(n))              # union-find over vertices
    top = {}                           # root vertex of a cluster -> last edge merged into it
    parent = {}

    def find(v):
        while comp[v] != v:
            comp[v] = comp[comp[v]]
            v = comp[v]
        return v

    for j in order:
        a, b = edges[j]
        ra, rb = find(a), find(b)
        for r in (ra, rb):
            if r in top:
                parent[top[r]] = j
        comp[ra] = rb
        top.pop(ra, None)
        top[rb] = j
    parent[order[-1]] = None
    return tuple(parent[j] for j in range(len(edges)))


def subtree_sizes(parent):
    """h_j = number of internal nodes in the subtree of j."""
    h = [1] * len(parent)
    for j in range(len(parent)):
        p = parent[j]
        while p is not None:
            h[p] += 1
            p = parent[p]
    return h


def is_linear_extension(order, parent):
    pos = {j: i for i, j in enumerate(order)}
    return all(parent[j] is None or pos[j] < pos[parent[j]] for j in range(len(parent)))


# ---- 1. Kruskal dendrograms and hook lengths ----
rnd = random.Random(7)
ok1 = ok2 = ok3 = True
ntrees = 0
for n in range(3, 9):
    for _ in range(5):
        edges = random_tree(n, rnd)
        ntrees += 1
        groups = {}
        for order in itertools.permutations(range(n - 1)):
            T = kruskal_dendrogram(n, edges, order)
            groups.setdefault(T, []).append(order)
        total = Fraction(0)
        for T, orders in groups.items():
            h = subtree_sizes(T)
            ext = [o for o in itertools.permutations(range(n - 1)) if is_linear_extension(o, T)]
            ok1 &= sorted(ext) == sorted(orders)
            ok2 &= len(orders) * math.prod(h) == math.factorial(n - 1)
            total += Fraction(1, math.prod(h))
        ok3 &= total == 1
check(f"{ntrees} random trees with 3-8 vertices: the orderings giving each Kruskal dendrogram T are exactly the "
      f"linear extensions of T", ok1)
check("their number is (n-1)!/prod_j h_j (hook-length formula)", ok2)
check("sum over dendrograms T ~ tau of prod_j 1/h_j = 1 exactly", ok3)


# ---- 2. nested-scale integral ----
def leaves_below(parent, nleaf_of):
    """number of leaves in the subtree of each internal node; nleaf_of[j] = leaves attached directly to j."""
    L = list(nleaf_of)
    for j in range(len(parent)):
        p = parent[j]
        while p is not None:
            L[p] += nleaf_of[j]
            p = parent[p]
    return L


def nested_integral(parent, nleaf_of, m, gamma):
    """Exact int_{R(T)} prod_j s_j^{1-gamma c_j} ds and prod_j 1/E_j."""
    k = len(parent)
    h = subtree_sizes(parent)
    L = leaves_below(parent, nleaf_of)
    # i(k) for every leaf: lowest ancestor (starting at the node it hangs from) with >= m+1 leaves
    c = [0] * k
    for j in range(k):
        for _ in range(nleaf_of[j]):
            q = j
            while q is not None and L[q] < m + 1:
                q = parent[q]
            if q is not None:
                c[q] += 1
    C = [0] * k
    for j in range(k):
        p = j
        while p is not None:
            C[p] += c[j]
            p = parent[p]
    E = [2 * h[j] - gamma * C[j] for j in range(k)]
    s = sp.symbols(f"s0:{k}", positive=True)
    children = {j: [i for i in range(k) if parent[i] == j] for j in range(k)}

    def F(j, t):                       # integral over the subtree of j with s_j <= t
        integrand = s[j] ** (1 - gamma * c[j])
        for ch in children[j]:
            integrand *= F(ch, s[j])
        return sp.integrate(integrand, (s[j], 0, t))

    root = [j for j in range(k) if parent[j] is None][0]
    lhs = sp.nsimplify(sp.simplify(F(root, 1)))
    rhs = math.prod(sp.Rational(1) / sp.nsimplify(e) for e in E)
    return sp.simplify(lhs), sp.simplify(rhs), E


g = sp.Rational(6, 5)
# caterpillar: node 0 = (12), node 1 = ((12)3), node 2 = root; leaves attached: 2, 1, 1
lhs, rhs, E = nested_integral((1, 2, None), (2, 1, 1), 2, g)
check(f"4-leaf caterpillar, m=2, gamma=6/5: integral = {lhs}, prod 1/E_j = {rhs} (E = {E})",
      lhs == rhs == sp.Rational(25, 24))
lhs, rhs, E = nested_integral((2, 2, None), (2, 2, 0), 2, g)
check(f"4-leaf balanced, m=2, gamma=6/5: integral = {lhs}, prod 1/E_j = {rhs}", lhs == rhs == sp.Rational(5, 24))


def random_dendrogram(nleaf, rnd):
    """Random binary dendrogram: merge random pairs of clusters; returns parents and leaves attached."""
    clusters = [("leaf", i) for i in range(nleaf)]
    parent, attached = [], []
    while len(clusters) > 1:
        a, b = rnd.sample(range(len(clusters)), 2)
        j = len(parent)
        parent.append(None)
        attached.append(0)
        for x in (clusters[a], clusters[b]):
            if x[0] == "leaf":
                attached[j] += 1
            else:
                parent[x[1]] = j
        clusters = [x for i, x in enumerate(clusters) if i not in (a, b)] + [("node", j)]
    return tuple(parent), tuple(attached)


rnd = random.Random(3)
good, cnt = True, 0
for nleaf in (3, 4, 5, 6):
    for m in (1, 2, 3):
        if nleaf < m + 1:
            continue
        for _ in range(3):
            P, A = random_dendrogram(nleaf, rnd)
            gmax = Fraction(2 * m, m + 1)
            gam = sp.Rational(rnd.randrange(1, 20), 20) * sp.Rational(gmax.numerator, gmax.denominator)
            lhs, rhs, E = nested_integral(P, A, m, gam)
            if min(E) <= 0:
                continue
            good &= sp.simplify(lhs - rhs) == 0
            cnt += 1
check(f"{cnt} random dendrograms with 3-6 leaves, m=1,2,3, rational gamma < 2m/(m+1): integral = prod 1/E_j exactly",
      good)


# ---- 3. tilde delta_k >= s_{i(k)} ----
def emst(pts):
    n = len(pts)
    D = np.linalg.norm(pts[:, None] - pts[None], axis=-1)
    in_tree, edges = {0}, []
    while len(in_tree) < n:
        best = min(((D[i, j], i, j) for i in in_tree for j in range(n) if j not in in_tree))
        edges.append((best[1], best[2], best[0]))
        in_tree.add(best[2])
    return D, edges


rng = np.random.default_rng(11)
good, cnt = True, 0
for _ in range(400):
    n = int(rng.integers(4, 13))
    pts = rng.random((n, 2)) * rng.choice([0.5, 1, 3])
    D, edges = emst(pts)
    order = sorted(range(n - 1), key=lambda j: edges[j][2])
    for m in (1, 2, 3):
        if n < m + 1:
            continue
        # single linkage: track for each point the size of its cluster and find the first merger reaching m+1
        comp = list(range(n))
        size = [1] * n

        def find(v):
            while comp[v] != v:
                v = comp[v]
            return v

        sk = [None] * n
        for j in order:
            a, b, length = edges[j]
            ra, rb = find(a), find(b)
            new = size[ra] + size[rb]
            if new >= m + 1:
                for k in range(n):
                    if sk[k] is None and find(k) in (ra, rb):
                        sk[k] = length
            comp[ra] = rb
            size[rb] = new
        for k in range(n):
            dm = np.sort(np.delete(D[k], k))[m - 1]
            good &= dm >= sk[k] * (1 - 1e-12)
            cnt += 1
check(f"{cnt} (configuration, point, m) samples: m-th-neighbour distance >= length of edge i(k)", good)
sys.exit(0 if ok else 1)
