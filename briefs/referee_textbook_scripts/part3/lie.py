"""Lie-algebraic data for simply-laced affine algebras, in the normalization of
atft_foundations.qmd Section 2: |alpha|^2 = 2, sum_j n_j alpha_j = 0, M0 = m^2 sum_j n_j a_j a_j^T, m = 1.

Nodes are labelled 0..r, node 0 affine (alpha_0 = -theta).  Bourbaki labelling of the finite nodes.
Species are labelled by finite nodes 1..r (particle a <-> node a, via the Perron-Frobenius vector),
with degenerate masses resolved by the characters of the centre (diagram automorphisms acting on
first-order Hirota data by scalars), in the convention delta_{z(j)} = chi_z(a) delta_j,
chi_z(a) = exp(-2 pi i (C^{-1})_{z(0),a}), which reproduces delta_j = omega^{ja} for a_n.
"""
import itertools
import numpy as np


def finite_edges(typ, r):
    if typ == 'a':
        return [(i, i + 1) for i in range(1, r)]
    if typ == 'd':
        e = [(i, i + 1) for i in range(1, r - 1)]
        e.append((r - 2, r))
        return e
    if typ == 'e':
        e = [(1, 3), (3, 4), (2, 4)] + [(i, i + 1) for i in range(4, r)]
        return e
    raise ValueError(typ)


def highest_root_marks(typ, r):
    if typ == 'a':
        return [1] * r
    if typ == 'd':
        return [1] + [2] * (r - 3) + [1, 1]
    if typ == 'e':
        return {6: [1, 2, 2, 3, 2, 1], 7: [2, 2, 3, 4, 3, 2, 1], 8: [2, 3, 4, 6, 5, 4, 3, 2]}[r]
    raise ValueError


def affine_attach(typ, r):
    """finite node(s) joined to node 0"""
    if typ == 'a':
        return [1, r] if r > 1 else [1]
    if typ == 'd':
        return [2]
    if typ == 'e':
        return {6: [2], 7: [1], 8: [8]}[r]


def coxeter_number(typ, r):
    return {'a': r + 1, 'd': 2 * r - 2, 'e': {6: 12, 7: 18, 8: 30}.get(r)}[typ]


class Algebra:
    def __init__(self, typ, r):
        self.typ, self.r = typ, r
        self.name = f"{typ}{r}"
        self.h = coxeter_number(typ, r)
        C = 2 * np.eye(r)
        for (i, j) in finite_edges(typ, r):
            C[i - 1, j - 1] = C[j - 1, i - 1] = -1
        self.Cfin = C
        L = np.linalg.cholesky(C)          # rows = simple roots
        marks = np.array(highest_root_marks(typ, r), float)
        theta = marks @ L
        self.alpha = np.vstack([-theta, L])  # (r+1) x r
        self.n = np.concatenate([[1.0], marks])
        A = self.alpha
        self.C = np.round(A @ A.T)           # affine Cartan matrix (simply laced), exact integers
        assert np.allclose(self.C @ self.n, 0)
        assert np.allclose(np.diag(self.C), 2)
        self.nbrs = [[k for k in range(r + 1) if k != j and abs(self.C[j, k]) > 0.5] for j in range(r + 1)]
        self.bondmult = [[int(round(-self.C[j, k])) for k in self.nbrs[j]] for j in range(r + 1)]
        self.N = np.diag(self.n)
        self.NC = self.N @ self.C
        self.M0 = (A.T * self.n) @ A
        self.Cinv = np.linalg.inv(C)
        self._automorphisms()
        self._species()

    # ------------------------------------------------------------------ symmetries
    def _automorphisms(self):
        r = self.r
        adj = (np.abs(self.C) > 0.5) & ~np.eye(r + 1, dtype=bool)
        autos = []
        # brute force via backtracking on permutations preserving adjacency
        nodes = list(range(r + 1))
        deg = adj.sum(1)

        def bt(perm):
            j = len(perm)
            if j == r + 1:
                autos.append(tuple(perm))
                return
            for c in nodes:
                if c in perm or deg[c] != deg[j] or self.n[c] != self.n[j]:
                    continue
                ok = all(adj[j, i] == adj[c, perm[i]] for i in range(j))
                if ok:
                    bt(perm + [c])
        bt([])
        self.autos = autos
        # linear action on h: R alpha_j = alpha_{g(j)}
        self.autoR = []
        for g in autos:
            # solve R from finite simple roots: R alpha_i = alpha_{g(i)} for i = 1..r
            Afin = self.alpha[1:]
            Ag = self.alpha[[g[i] for i in range(1, r + 1)]]
            R = np.linalg.solve(Afin, Ag).T   # R @ alpha_i = alpha_g(i)
            assert np.allclose(R @ self.alpha.T, self.alpha[list(g)].T)
            self.autoR.append(R)
        # centre: automorphisms whose linear part lies in the Weyl group
        self.centre = [g for g, R in zip(autos, self.autoR) if self._in_weyl(R)]

    def _in_weyl(self, R):
        r = self.r
        rho = self.Cinv.sum(0) @ self.alpha[1:]  # sum of fundamental weights (in root coords -> vector)
        # fundamental weights: omega_i = sum_j Cinv[i,j] alpha_j
        v = R @ rho
        w = np.eye(r)
        for _ in range(10000):
            done = True
            for i in range(r):
                ai = self.alpha[i + 1]
                if v @ ai < -1e-9:
                    s = np.eye(r) - np.outer(ai, ai)
                    v = s @ v
                    w = s @ w
                    done = False
            if done:
                break
        return np.allclose(w @ R, np.eye(r), atol=1e-8)

    # ------------------------------------------------------------------ species
    def _species(self):
        r, h = self.r, self.h
        ev, V = np.linalg.eig(self.NC)
        idx = np.argsort(ev.real)
        ev, V = ev[idx].real, V[:, idx]
        assert abs(ev[0]) < 1e-9
        ev, V = ev[1:], V[:, 1:]
        # Perron-Frobenius masses -> node labels
        w, U = np.linalg.eigh(self.Cfin)
        pf = np.abs(U[:, 0])
        # group eigenvalues
        groups = []
        for i, e in enumerate(ev):
            for g in groups:
                if abs(ev[g[0]] - e) < 1e-8:
                    g.append(i)
                    break
            else:
                groups.append([i])
        masses = np.sqrt(ev)
        scale = None
        # scale: lightest mass <-> smallest PF component
        scale = masses.min() / pf.min()
        node_mass = pf * scale
        self.mass_of_node = {a + 1: node_mass[a] for a in range(r)}
        species = {}
        for g in groups:
            m = masses[g[0]]
            nodes = [a + 1 for a in range(r) if abs(node_mass[a] - m) < 1e-7]
            assert len(nodes) == len(g), (self.name, m, nodes, g)
            sub = V[:, g]
            if len(g) == 1:
                species[nodes[0]] = sub[:, 0]
                continue
            # diagonalize centre action on the degenerate eigenspace
            vecs = self._split_by_centre(sub, nodes)
            species.update(vecs)
        # normalize delta_0 = 1 when possible (delta_0 = 1 for all species: check)
        for a in species:
            d = species[a]
            d = d / d[0]
            species[a] = d
        self.delta = species
        self.mass = {a: float(np.sqrt(np.real((self.NC @ species[a])[0] / species[a][0]))) for a in species}
        for a in species:
            assert np.allclose(self.NC @ species[a], self.mass[a] ** 2 * species[a], atol=1e-8)
        # channel vectors u_b = sum_j alpha_j delta_j, normalized u_b^T u_bbar = 1
        self.conj = {}
        raw = {a: self.alpha.T @ species[a] for a in species}
        for a in species:
            best = None
            for b in species:
                if abs(self.mass[a] - self.mass[b]) < 1e-8:
                    if abs(raw[a] @ raw[b]) > 1e-8:
                        best = b
            self.conj[a] = best
        self.u = {}
        for a in species:
            b = self.conj[a]
            s = raw[a] @ raw[b]
            self.u[a] = raw[a] / np.sqrt(s)
        for a in species:
            for b in species:
                val = self.u[a] @ self.u[b]
                tgt = 1.0 if b == self.conj[a] else 0.0
                assert abs(val - tgt) < 1e-8, (a, b, val)

    def chi(self, z, a):
        """character of centre element z on species a: exp(-2 pi i Cinv[z(0), a])"""
        i = z[0]
        if i == 0:
            return 1.0 + 0j
        return np.exp(-2j * np.pi * self.Cinv[i - 1, a - 1])

    def _split_by_centre(self, sub, nodes):
        # find combination diagonalizing all centre permutations, with characters chi(z, a)
        out = {}
        for a in nodes:
            # projector onto common eigenspace with eigenvalues chi(z,a)
            P = np.eye(sub.shape[0], dtype=complex)
            vec = None
            # build vector: sum over centre of chi^{-1} * permuted
            # acting on delta: (g.delta)_j = delta_{g(j)}
            base = sub
            cands = []
            for k in range(sub.shape[1]):
                d = sub[:, k].astype(complex)
                acc = np.zeros_like(d)
                for z in self.centre:
                    pd = d[list(z)]
                    acc += pd / self.chi(z, a)
                if np.linalg.norm(acc) > 1e-8:
                    cands.append(acc)
            assert cands, (self.name, a)
            v = max(cands, key=np.linalg.norm)
            # check eigen-property
            for z in self.centre:
                assert np.allclose(v[list(z)], self.chi(z, a) * v, atol=1e-8), (self.name, a, z)
            out[a] = v
        # check independence
        M = np.array([out[a] for a in nodes]).T
        assert np.linalg.matrix_rank(M, tol=1e-8) == len(nodes), (self.name, nodes)
        return out

    def summary(self):
        s = [f"{self.name}: h={self.h}, n={self.n.astype(int).tolist()}, |centre|={len(self.centre)}, |Aut|={len(self.autos)}"]
        for a in sorted(self.mass):
            s.append(f"  species {a}: m^2={self.mass[a]**2:.6f} m={self.mass[a]:.6f} conj={self.conj[a]}")
        return "\n".join(s)


if __name__ == '__main__':
    for typ, r in [('a', 3), ('a', 4), ('d', 4), ('d', 5), ('d', 6), ('d', 7), ('e', 6), ('e', 7), ('e', 8)]:
        A = Algebra(typ, r)
        print(A.summary())
        print("  centre:", A.centre)


def delta_mp(alg, a, dps=40, iters=6):
    """high-precision first-order data for species a (inverse iteration in mpmath, preserving
    the centre character); returns (delta list of mpc, m_a^2 as mpf)"""
    import mpmath as mp
    mp.mp.dps = dps
    r1 = alg.r + 1
    NC = mp.matrix([[int(round(alg.NC[i, j])) for j in range(r1)] for i in range(r1)])
    v = mp.matrix([mp.mpc(complex(x)) for x in alg.delta[a]])
    lam = mp.mpf(alg.mass[a] ** 2)
    for it in range(iters):
        # Rayleigh-type update using left eigvec: N^{-1} symmetrizes: N^{-1/2} (NC) N^{1/2} symmetric
        # use generalized Rayleigh quotient with weight N^{-1}: (v^T C v)/(v^T N^{-1} v)
        Cv = mp.matrix([sum(int(round(alg.C[i, j])) * v[j] for j in range(r1)) for i in range(r1)])
        num = sum(mp.conj(v[i]) * Cv[i] for i in range(r1))
        den = sum(mp.conj(v[i]) * v[i] / int(round(alg.n[i])) for i in range(r1))
        lam = num / den
        shift = lam + mp.mpf(10) ** (-dps // 2)
        A = NC - shift * mp.eye(r1)
        v = mp.lu_solve(A, v)
        v = v / v[0]
    # character projection with exact phases
    acc = [mp.mpc(0)] * r1
    for zz in alg.centre:
        i0 = zz[0]
        chi = mp.mpc(1) if i0 == 0 else mp.expj(-2 * mp.pi * mp.mpf(alg.Cinv[i0 - 1, a - 1]))
        # Cinv entries are rationals with small denominators; refine
        if i0 != 0:
            from fractions import Fraction
            fr = Fraction(float(alg.Cinv[i0 - 1, a - 1])).limit_denominator(12)
            chi = mp.expj(-2 * mp.pi * mp.mpf(fr.numerator) / fr.denominator)
        for j in range(r1):
            acc[j] += v[zz[j]] / chi
    acc = [x / acc[0] for x in acc]
    Cv = [sum(int(round(alg.C[i, j])) * acc[j] for j in range(r1)) for i in range(r1)]
    num = sum(mp.conj(acc[i]) * Cv[i] for i in range(r1))
    den = sum(mp.conj(acc[i]) * acc[i] / int(round(alg.n[i])) for i in range(r1))
    lam = num / den
    resid = max(abs(sum(NC[i, j] * acc[j] for j in range(r1)) - lam * acc[i]) for i in range(r1))
    return acc, mp.re(lam), resid
