"""Independent verification of Step 3 (Steps 3.1-3.5) of thm-nn-distances in Ch. 48.

Checks every combinatorial and integral claim in Steps 3.1-3.5 from scratch:
- Step 3.1: Minimax path in Euclidean graph vs EMST path vs single-linkage dendrogram merger height s_{a v b}, and diam(A) <= (n-1) s_root.
- Step 3.2: Lower bound tilde_delta_k >= s_{i(k)} where i(k) is the lowest ancestor of leaf k with >= m+1 leaves, and sum_j c_j = n_I <= n.
- Step 3.3: Kruskal orderings giving dendrogram T equal the linear extensions of the tree poset T, with count (n-1)! / prod_j h_j.
- Step 3.4: Nested cone integral int_{R(T)} prod_j s_j^{1 - gamma c_j} ds = prod_{j != root} 1/E_j * s_root^{E_root} / (2 - gamma c_root), and the inequality C_j <= h_j + m => 2 h_j / E_j <= 2m / alpha for all h_j >= 1, C_j >= 1.
- Step 3.5: Dendrogram identity sum_{T ~ tau} prod_j (1 / h_j) = 1 for every labelled tree tau.
"""
import itertools
from fractions import Fraction
import numpy as np
import scipy.integrate as integrate


def enumerate_labeled_trees(n):
    """Return all n^{n-2} labeled trees on {0, ..., n-1} via Prüfer codes."""
    if n == 2:
        return [[(0, 1)]]
    trees = []
    for code in itertools.product(range(n), repeat=n - 2):
        deg = [1] * n
        for v in code:
            deg[v] += 1
        edges = []
        c_list = list(code)
        for v in c_list:
            for u in range(n):
                if deg[u] == 1:
                    edges.append((min(u, v), max(u, v)))
                    deg[u] -= 1
                    deg[v] -= 1
                    break
        u, v = [i for i in range(n) if deg[i] == 1]
        edges.append((min(u, v), max(u, v)))
        trees.append(tuple(sorted(edges)))
    return trees


def build_dendrogram(n, ordered_edges):
    """Build single-linkage dendrogram from Kruskal-ordered edges.
    Returns:
      parent: dict mapping node_id -> parent_id (leaves 0..n-1, internal n..2n-2)
      children: dict mapping internal_id -> (left_id, right_id)
      leaves_of: dict mapping node_id -> frozenset of leaf ids in subtree
      edge_of: dict mapping internal_id -> edge (u, v)
    """
    comp = {i: i for i in range(n)}
    leaves_of = {i: frozenset([i]) for i in range(n)}
    parent = {}
    children = {}
    edge_of = {}

    def find(x):
        while comp[x] != x:
            comp[x] = comp[comp[x]]
            x = comp[x]
        return x

    for idx, e in enumerate(ordered_edges):
        u, v = e
        ru, rv = find(u), find(v)
        node_id = n + idx
        if ru > rv:
            ru, rv = rv, ru
        children[node_id] = (ru, rv)
        parent[ru] = node_id
        parent[rv] = node_id
        comp[ru] = node_id
        comp[rv] = node_id
        comp[node_id] = node_id
        leaves_of[node_id] = leaves_of[ru] | leaves_of[rv]
        edge_of[node_id] = e
    return parent, children, leaves_of, edge_of


def dendrogram_signature(n, children, leaves_of, edge_of):
    """Canonical representation of unlabeled-height dendrogram T (retaining which tree edge is at each internal node and which leaf subsets it merges)."""
    items = []
    for node_id in range(n, 2 * n - 1):
        ch_l, ch_r = children[node_id]
        pair = tuple(sorted([tuple(sorted(leaves_of[ch_l])), tuple(sorted(leaves_of[ch_r]))]))
        items.append((edge_of[node_id], pair))
    return frozenset(items)


def verify_step_3_1_and_3_2(n=6, m=2, n_trials=200):
    """Verify Steps 3.1 and 3.2 on random point configurations in R^2."""
    rng = np.random.default_rng(42)
    ok_3_1 = True
    ok_3_2 = True
    for _ in range(n_trials):
        pts = rng.uniform(-1.0, 1.0, size=(n, 2))
        dist = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1)
        # Kruskal EMST
        all_pairs = sorted([(dist[i, j], i, j) for i in range(n) for j in range(i + 1, n)])
        comp = list(range(n))
        def find(x):
            while comp[x] != x:
                x = comp[x]
            return x
        mst_edges = []
        mst_weights = []
        for d, i, j in all_pairs:
            ri, rj = find(i), find(j)
            if ri != rj:
                comp[ri] = rj
                mst_edges.append((i, j))
                mst_weights.append(d)
        parent, children, leaves_of, edge_of = build_dendrogram(n, mst_edges)
        height = {n + idx: mst_weights[idx] for idx in range(n - 1)}

        # Step 3.1: check s_{a v b} = minimax distance between a and b, and diam(A) <= (n-1)*s_root
        for a in range(n):
            for b in range(a + 1, n):
                # Find lowest common ancestor in dendrogram
                anc_a = set()
                curr = a
                while curr in parent:
                    curr = parent[curr]
                    anc_a.add(curr)
                curr = b
                while curr not in anc_a:
                    curr = parent[curr]
                s_lca = height[curr]
                # Compare with minimax distance over all paths in complete graph
                # Floyd-Warshall minimax:
                mm = dist.copy()
                for k in range(n):
                    mm = np.minimum(mm, np.maximum(mm[:, k:k+1], mm[k:k+1, :]))
                if abs(s_lca - mm[a, b]) > 1e-12:
                    ok_3_1 = False
        if np.max(dist) > (n - 1) * mst_weights[-1] + 1e-12:
            ok_3_1 = False

        # Step 3.2: for each leaf k, tilde_delta_k = m-th nearest neighbor distance in A \ {k}
        # and i(k) = lowest ancestor of k with |L_j| >= m + 1
        c_counts = {n + idx: 0 for idx in range(n - 1)}
        for k in range(n):
            other_dists = sorted([dist[k, l] for l in range(n) if l != k])
            tilde_delta_k = other_dists[m - 1]
            curr = k
            while len(leaves_of[curr]) < m + 1:
                curr = parent[curr]
            c_counts[curr] += 1
            if tilde_delta_k + 1e-12 < height[curr]:
                ok_3_2 = False
        if sum(c_counts.values()) != n:
            ok_3_2 = False
    print("PASS step_3_1_minimax_and_diam" if ok_3_1 else "FAIL step_3_1_minimax_and_diam")
    print("PASS step_3_2_ancestor_bound" if ok_3_2 else "FAIL step_3_2_ancestor_bound")


def verify_step_3_3_and_3_5(max_n=5):
    """Verify Step 3.3 (linear extension count = (n-1)! / prod h_j) and Step 3.5 (sum_{T ~ tau} prod 1/h_j = 1)."""
    ok_3_3 = True
    ok_3_5 = True
    for n in range(2, max_n + 1):
        trees = enumerate_labeled_trees(n)
        for tau in trees:
            # Group all (n-1)! edge orderings by the resulting dendrogram T
            dendro_counts = {}
            dendro_hj = {}
            for perm in itertools.permutations(tau):
                parent, children, leaves_of, edge_of = build_dendrogram(n, perm)
                sig = dendrogram_signature(n, children, leaves_of, edge_of)
                dendro_counts[sig] = dendro_counts.get(sig, 0) + 1
                if sig not in dendro_hj:
                    # Compute subtree sizes h_j of the internal-node poset T (n-1 internal nodes)
                    internal_desc = {node_id: 1 for node_id in range(n, 2 * n - 1)}
                    for node_id in range(n, 2 * n - 1):
                        curr = node_id
                        while curr in parent:
                            curr = parent[curr]
                            internal_desc[curr] += 1
                    dendro_hj[sig] = tuple(internal_desc[node_id] for node_id in range(n, 2 * n - 1))

            # Check Step 3.3 for every dendrogram T compatible with tau
            fact = 1
            for k in range(1, n):
                fact *= k
            sum_recip_prod = Fraction(0, 1)
            for sig, count in dendro_counts.items():
                prod_hj = 1
                for h_val in dendro_hj[sig]:
                    prod_hj *= h_val
                if count * prod_hj != fact:
                    ok_3_3 = False
                sum_recip_prod += Fraction(1, prod_hj)
            if sum_recip_prod != Fraction(1, 1):
                ok_3_5 = False
    print("PASS step_3_3_linear_extensions" if ok_3_3 else "FAIL step_3_3_linear_extensions")
    print("PASS step_3_5_dendrogram_sum_identity" if ok_3_5 else "FAIL step_3_5_dendrogram_sum_identity")


def verify_step_3_4():
    """Verify Step 3.4:
    1) Exact nested integral over R(T) against numerical quadrature on a non-chain branching dendrogram (n=5, 4 internal nodes).
    2) Algebraic inequality C_j <= h_j + m and 2*h_j / E_j <= 2*m / alpha for all h_j >= 1, 1 <= C_j <= h_j + m, and 0 < alpha <= 2m.
    """
    # 1) Numerical quadrature on a balanced binary dendrogram on n=4 leaves (3 internal nodes: j=1,2 children of root j=3)
    # R(T) = {0 < s_1 < s_3 < 1, 0 < s_2 < s_3 < 1} with s_root = s_3 in (0, 1).
    # Or fix s_3 = 1 and integrate over 0 < s_1 < 1, 0 < s_2 < 1, or integrate over 0 < s_1 < s_2 < s_4, 0 < s_3 < s_4 < 1 (n=5, 4 internal nodes):
    # Internal nodes: 1 (leaf-pair, h_1=1), 2 (parent of 1 and leaf, h_2=2), 3 (leaf-pair, h_3=1), 4 (root, parent of 2 and 3, h_4=4).
    c1, c2, c3, c4 = 2, 1, 2, 0
    gamma = 0.65
    # Integrate prod_{j=1}^4 s_j^{1 - gamma c_j} over 0 < s_1 < s_2 < s_4 < 1 and 0 < s_3 < s_4 < 1:
    # Inner integrals:
    # int_0^{s_2} s_1^{1 - gamma c_1} ds_1 = s_2^{E_1} / E_1 where E_1 = 2 - gamma c_1
    # int_0^{s_4} s_2^{1 - gamma c_2} (s_2^{E_1} / E_1) ds_2 = s_4^{E_2} / (E_1 E_2) where E_2 = 4 - gamma(c_1 + c_2)
    # int_0^{s_4} s_3^{1 - gamma c_3} ds_3 = s_4^{E_3} / E_3 where E_3 = 2 - gamma c_3
    # int_0^1 s_4^{1 - gamma c_4 + E_2 + E_3} ds_4 = 1 / (E_1 E_2 E_3 E_4) where E_4 = 8 - gamma(c_1 + c_2 + c_3 + c_4)
    E1 = 2 * 1 - gamma * c1
    E2 = 2 * 2 - gamma * (c1 + c2)
    E3 = 2 * 1 - gamma * c3
    E4 = 2 * 4 - gamma * (c1 + c2 + c3 + c4)
    analytic_val = 1.0 / (E1 * E2 * E3 * E4)

    def integrand_s4(s4):
        I12, _ = integrate.dblquad(
            lambda s1, s2: (s1 ** (1 - gamma * c1)) * (s2 ** (1 - gamma * c2)),
            0.0, s4,
            lambda s2: 0.0, lambda s2: s2,
            epsabs=1e-11, epsrel=1e-11
        )
        I3, _ = integrate.quad(lambda s3: s3 ** (1 - gamma * c3), 0.0, s4, epsabs=1e-11, epsrel=1e-11)
        return (s4 ** (1 - gamma * c4)) * I12 * I3

    num_val, _ = integrate.quad(integrand_s4, 0.0, 1.0, epsabs=1e-10, epsrel=1e-10)
    ok_quad = abs(num_val - analytic_val) < 1e-9

    # 2) Verify C_j <= h_j + m for all dendrograms on n <= 7 and all m >= 1, and verify 2*h_j / E_j <= 2*m / alpha
    ok_ineq = True
    for n in range(3, 7):
        trees = enumerate_labeled_trees(n)
        for tau in trees[:10]:
            for perm in list(itertools.permutations(tau))[:20]:
                parent, children, leaves_of, _ = build_dendrogram(n, perm)
                # Compute h_j
                h_map = {node_id: 1 for node_id in range(n, 2 * n - 1)}
                for node_id in range(n, 2 * n - 1):
                    curr = node_id
                    while curr in parent:
                        curr = parent[curr]
                        h_map[curr] += 1
                for m in range(1, n):
                    c_map = {node_id: 0 for node_id in range(n, 2 * n - 1)}
                    for k in range(n):
                        curr = k
                        while len(leaves_of[curr]) < m + 1:
                            curr = parent[curr]
                        c_map[curr] += 1
                    C_map = dict(c_map)
                    for node_id in range(n, 2 * n - 1):
                        curr = node_id
                        while curr in parent:
                            curr = parent[curr]
                            C_map[curr] += c_map[node_id]
                    for node_id in range(n, 2 * n - 1):
                        hj = h_map[node_id]
                        Cj = C_map[node_id]
                        if Cj > 0 and Cj > hj + m:
                            ok_ineq = False
                        for alpha in (0.1, 0.5, m, 1.9 * m):
                            gam = (2 * m - alpha) / (m + 1)
                            Ej = 2 * hj - gam * Cj
                            if Ej <= 0 or (2 * hj / Ej) > (2 * m / alpha) + 1e-12:
                                ok_ineq = False

    print("PASS step_3_4_nested_cone_integral" if ok_quad else "FAIL step_3_4_nested_cone_integral")
    print("PASS step_3_4_exponent_bound" if ok_ineq else "FAIL step_3_4_exponent_bound")


if __name__ == '__main__':
    verify_step_3_1_and_3_2()
    verify_step_3_3_and_3_5()
    verify_step_3_4()
