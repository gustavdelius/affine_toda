"""Small root-system library. Normalization: longest roots have length^2 2.
Cartan matrix a_ij = 2 (alpha_i, alpha_j)/(alpha_i, alpha_i)  (book eq-cartan-matrix)."""
import numpy as np
import itertools

def gram(typ, n):
    """Gram matrix (alpha_i, alpha_j) of simple roots, Bourbaki numbering (0-based)."""
    G = np.zeros((n, n))
    def link(i, j, v):
        G[i, j] = G[j, i] = v
    if typ in 'ADE':
        for i in range(n): G[i, i] = 2
        if typ == 'A':
            for i in range(n-1): link(i, i+1, -1)
        elif typ == 'D':
            for i in range(n-2): link(i, i+1, -1)
            link(n-3, n-1, -1)
        elif typ == 'E':   # Bourbaki: 1-3-4-5-6(-7-8), 2 attached to 4
            # 0-based: chain 0-2-3-4-5-..., node 1 attached to node 3
            chain = [0, 2, 3] + list(range(4, n))
            for a, b in zip(chain, chain[1:]): link(a, b, -1)
            link(1, 3, -1)
    elif typ == 'B':   # alpha_n short (len^2 1)
        for i in range(n-1): G[i, i] = 2
        G[n-1, n-1] = 1
        for i in range(n-1): link(i, i+1, -1)
    elif typ == 'C':   # alpha_1..alpha_{n-1} short (1), alpha_n long (2)
        for i in range(n-1): G[i, i] = 1
        G[n-1, n-1] = 2
        for i in range(n-2): link(i, i+1, -0.5)
        link(n-2, n-1, -1)
    elif typ == 'F':   # alpha1, alpha2 long; alpha3, alpha4 short
        G = np.array([[2, -1, 0, 0], [-1, 2, -1, 0], [0, -1, 1, -0.5], [0, 0, -0.5, 1]], float)
    elif typ == 'G':   # alpha1 short (2/3), alpha2 long (2)
        G = np.array([[2/3, -1], [-1, 2]], float)
    return G

def cartan_from_gram(G):
    return np.array([[2*G[i, j]/G[i, i] for j in range(len(G))] for i in range(len(G))])

def simple_roots(G):
    return np.linalg.cholesky(G)          # rows alpha_i, alpha_i . alpha_j = G_ij

def refl(a):
    return np.eye(len(a)) - 2*np.outer(a, a)/(a @ a)

def all_roots(al):
    roots = {tuple(np.round(a, 10)): a for a in al}
    frontier = list(al)
    while frontier:
        new = []
        for v in frontier:
            for a in al:
                w = v - 2*(v @ a)/(a @ a)*a
                t = tuple(np.round(w, 10))
                if t not in roots:
                    roots[t] = w; new.append(w)
        frontier = new
    return list(roots.values())

class RootSystem:
    def __init__(self, typ, n):
        self.typ, self.n = typ, n
        self.G = gram(typ, n); self.A = cartan_from_gram(self.G)
        self.al = simple_roots(self.G)
        self.roots = all_roots(list(self.al))
        coords = [np.linalg.solve(self.al.T, r) for r in self.roots]
        self.coords = [np.round(c).astype(int) for c in coords]
        assert all(np.allclose(c, ci) for c, ci in zip(coords, self.coords))
        self.pos = [r for r, c in zip(self.roots, self.coords) if c.sum() > 0]
        k = int(np.argmax([c.sum() for c in self.coords]))
        self.theta = self.roots[k]; self.marks = self.coords[k]
        self.len2 = np.diag(self.G)
        self.comarks = self.marks*self.len2/2
        self.h = 1 + self.marks.sum(); self.hv = 1 + self.comarks.sum()
        self.lam = np.linalg.inv(self.A) @ (self.al/ (self.len2[:, None]/2))  # placeholder, recomputed below
        # fundamental weights: alpha_i^vee . lambda_j = delta_ij
        alv = 2*self.al/self.len2[:, None]
        self.alv = alv
        self.lam = np.linalg.solve(alv, np.eye(n)).T     # rows lambda_j with alv_i . lambda_j = delta_ij
        self.colv = np.linalg.solve(self.al, np.eye(n)).T  # rows lambda_j^vee with alpha_i . lambda_j^vee = delta_ij

    def coxeter(self, order=None):
        order = range(self.n) if order is None else order
        w = np.eye(self.n)
        for i in order: w = w @ refl(self.al[i])
        return w

    def bicolour(self):
        col = {0: 0}; stack = [0]
        while stack:
            i = stack.pop()
            for j in range(self.n):
                if j != i and abs(self.A[i, j]) > 1e-12 and j not in col:
                    col[j] = 1 - col[i]; stack.append(j)
        return [i for i in range(self.n) if col[i] == 0], [i for i in range(self.n) if col[i] == 1]

def order_of(w, maxp=200):
    P = np.eye(len(w))
    for p in range(1, maxp):
        P = P @ w
        if np.allclose(P, np.eye(len(w)), atol=1e-8): return p
    return None

def exponents_of(w, h):
    ev = np.linalg.eigvals(w)
    m = np.sort(np.round((np.angle(ev) % (2*np.pi))*h/(2*np.pi), 6))
    return m
