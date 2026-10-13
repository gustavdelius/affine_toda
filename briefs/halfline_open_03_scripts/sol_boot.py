"""Soliton bootstrap V^(x n) -> Vbar for the s = +-1 reflection matrices, relative to the diagonal K = 1.

M = full reflection of n solitons V(th_1) (x) ... (x) V(th_n) (boundary on the right): process particles
k = n..1; move incoming k to the right end through the already reflected ones (R-check V (x) Vbar), reflect with
K at the right end.  Final order: Vbar(-th_n) (x) ... (x) Vbar(-th_1).
iota_in: Vbar(th) -> (x)_k V(th_k), iota_out: V(-th) -> (x)_k Vbar(-th_{n+1-k}).
M iota_in = c * iota_out K_W(th), with K_W = K^T (s = +-1) or 1 (diagonal).
P(mu) = c^(s)/c^(d) is independent of the normalization of the embeddings and of the bulk scalar factors.
Bootstrap for A^(s) = A_d * Y:  Y(mu) = P(mu) prod_k Y(mu_k),  mu_k = mu + T - k om.
"""
import itertools
import numpy as np
import scipy.sparse as sps
from boundary_common import q_phys, an_vector
from nk_common import fast_rcheck, k_book


def tensor_gen(reps, kind, j):
    """Iterated Part VI coproduct of generator (kind, j) on reps[0] (x) ... (x) reps[-1] (sparse)."""
    n = len(reps)
    tot = None
    for p in range(n):
        mats = []
        for r, R in enumerate(reps):
            if kind == 'e':
                M = R[('k', j)] if r < p else (R[('e', j)] if r == p else np.eye(R['dim']))
            else:   # f: f at p, k^-1 to the right
                M = np.eye(R['dim']) if r < p else (R[('f', j)] if r == p else np.linalg.inv(R[('k', j)]))
            mats.append(sps.csr_matrix(M))
        T = mats[0]
        for M in mats[1:]:
            T = sps.kron(T, M, format='csr')
        tot = T if tot is None else tot + T
    return tot


def embed(rep_small, reps, wt_small, wt_big):
    """Intertwiner X: rep_small -> (x) reps, weight restricted (weights as integer vectors mod (1,...,1))."""
    nodes = sorted(k[1] for k in rep_small if isinstance(k, tuple) and k[0] == 'e')
    d_s, D = rep_small['dim'], int(np.prod([R['dim'] for R in reps]))
    norm = lambda w: tuple(np.array(w) - np.array(w)[0])
    ws = [norm(w) for w in wt_small]
    wb = [norm(w) for w in wt_big]
    cols = [(o, i) for i in range(d_s) for o in range(D) if wb[o] == ws[i]]
    idx = {c: a for a, c in enumerate(cols)}
    rows, cs, vals = [], [], []
    nr = 0
    for j in nodes:
        for kind in 'ef':
            G = tensor_gen(reps, kind, j).tocsc()
            g = rep_small[(kind, j)]
            # residual[r, i] = sum_o G[r,o] X[o,i] - sum_i' X[r,i'] g[i',i]
            for a, (o, i) in enumerate(cols):
                col = G[:, o]
                for r, v in zip(col.indices, col.data):
                    rows.append(nr + r*d_s + i); cs.append(a); vals.append(v)
                for i2 in range(d_s):
                    if g[i, i2] != 0:      # X[o,i] contributes to residual[o, i2] with -g[i, i2]
                        rows.append(nr + o*d_s + i2); cs.append(a); vals.append(-g[i, i2])
            nr += D*d_s
    Mx = sps.csr_matrix((vals, (rows, cs)), shape=(nr, len(cols)))
    Gm = (Mx.conj().T @ Mx).toarray()
    ev, U = np.linalg.eigh(Gm)
    assert ev[0] < 1e-12*ev[-1] and ev[1] > 1e-8*ev[-1], (ev[:3]/ev[-1])
    v = U[:, 0]
    X = np.zeros((D, d_s), complex)
    for (o, i), val in zip(cols, v):
        X[o, i] = val
    return X


def apply_two(vecs, M, p, n, N):
    """Apply N^2 x N^2 matrix M on tensor factors (p, p+1) of vecs (N^n x m)."""
    m = vecs.shape[1]
    t = vecs.reshape((N**p, N*N, N**(n - p - 2), m))
    t = np.einsum('ab,ibjm->iajm', M, t)
    return t.reshape(N**n, m)


def apply_one(vecs, K, p, n, N):
    m = vecs.shape[1]
    t = vecs.reshape((N**p, N, N**(n - p - 1), m))
    t = np.einsum('ab,ibjm->iajm', K, t)
    return t.reshape(N**n, m)


def boot_coef(n, om, th, s, order=+1):
    """Returns (c, residual) with M iota_in = c iota_out K_W.  s = 0: diagonal K = 1."""
    q = q_phys(om); N = n + 1; T = N*om/2
    mu = -1j*T*th/np.pi
    mus = [mu + T - k*om for k in range(1, n + 1)]
    if order < 0:
        mus = mus[::-1]
    ths = [1j*np.pi*m/T for m in mus]
    Kv = (lambda t: np.eye(N)) if s == 0 else (lambda t: k_book(n, om, t, s))
    Kw = (lambda t: np.eye(N)) if s == 0 else (lambda t: k_book(n, om, t, s).T)
    E = np.eye(N, dtype=int)
    Vs = [an_vector(n, q, np.exp(om*t)) for t in ths]
    wt_in = [-E[i] for i in range(N)]
    wt_big_in = [sum(E[i] for i in tup) for tup in itertools.product(range(N), repeat=n)]
    io_in = embed(an_vector(n, q, np.exp(om*th), True), Vs, wt_in, wt_big_in)
    Ws = [an_vector(n, q, np.exp(-om*t), True) for t in ths[::-1]]
    wt_out = [E[i] for i in range(N)]
    wt_big_out = [-sum(E[i] for i in tup) for tup in itertools.product(range(N), repeat=n)]
    io_out = embed(an_vector(n, q, np.exp(-om*th)), Ws, wt_out, wt_big_out)
    # reflection
    v = io_in.copy()
    # position contents: list of ('in', k) or ('out', k)
    pos = [('in', k) for k in range(n)]
    for k in range(n - 1, -1, -1):
        p = k
        while p < n - 1:
            # swap incoming k at p with reflected (out, j) at p+1
            j = pos[p + 1][1]
            R = fast_rcheck(an_vector(n, q, np.exp(om*ths[k])), an_vector(n, q, np.exp(-om*ths[j]), True))
            v = apply_two(v, R, p, n, N)
            pos[p], pos[p + 1] = pos[p + 1], pos[p]
            p += 1
        v = apply_one(v, Kv(ths[k]), n - 1, n, N)
        pos[n - 1] = ('out', k)
    target = io_out @ Kw(th)
    c = np.vdot(target.ravel(), v.ravel())/np.vdot(target.ravel(), target.ravel())
    return c, np.linalg.norm(v - c*target)/np.linalg.norm(v)


if __name__ == '__main__':
    import sys
    om = 0.37
    for n in (2, 3):
        for th in (0.31 + 0.07j,):
            for s in (0, 1, -1):
                try:
                    print(n, s, boot_coef(n, om, th, s))
                except AssertionError as e:
                    print(n, s, 'embedding failed', e)
