import numpy as np
from dtw import rep, weights
def coprod(a, b, kind, i):
    # Delta(e)=e(x)1 + k(x)e ; Delta(f)=f(x)k^-1 + 1(x)f   (a, b = rep dicts)
    ea, fa, ka, _ = a; eb, fb, kb, _ = b; Na, Nb = ka[0].shape[0], kb[0].shape[0]
    if kind == 'e': return np.kron(ea[i], np.eye(Nb)) + np.kron(ka[i], eb[i])
    return np.kron(fa[i], np.linalg.inv(kb[i])) + np.kron(np.eye(Na), fb[i])
def intertwiner(A, B, nodes, wA, wB):
    """Solve  Rc * Delta_{A,B}(g) = Delta_{B,A}(g) * Rc  for Rc: A(x)B -> B(x)A, weight-preserving."""
    Na, Nb = wA.shape[0], wB.shape[0]; D = Na*Nb
    win = (wA[:, None, :] + wB[None, :, :]).reshape(D, -1); wout = (wB[:, None, :] + wA[None, :, :]).reshape(D, -1)
    pairs = [(r, c) for r in range(D) for c in range(D) if np.allclose(wout[r], win[c])]
    col = {p: t for t, p in enumerate(pairs)}; M = len(pairs)
    rows = []
    for kind in ('e', 'f'):
        for i in nodes:
            L = coprod(A, B, kind, i); Rm = coprod(B, A, kind, i)   # Rc L - Rm Rc = 0
            # entry (r, c): sum_m Rc[r,m] L[m,c] - sum_m Rm[r,m] Rc[m,c]
            for r in range(D):
                for c in range(D):
                    row = {}
                    for m in np.nonzero(np.abs(L[:, c]) > 1e-14)[0]:
                        if (r, m) in col: row[col[(r, m)]] = row.get(col[(r, m)], 0) + L[m, c]
                    for m in np.nonzero(np.abs(Rm[r, :]) > 1e-14)[0]:
                        if (m, c) in col: row[col[(m, c)]] = row.get(col[(m, c)], 0) - Rm[r, m]
                    if row:
                        v = np.zeros(M, complex)
                        for t, val in row.items(): v[t] = val
                        rows.append(v)
    Amat = np.array(rows); _, s, Vh = np.linalg.svd(Amat)
    null = np.sum(s < 1e-9*s[0]) + max(0, M - len(s))
    sol = Vh[-1].conj(); Rc = np.zeros((D, D), complex)
    for (r, c), t in col.items(): Rc[r, c] = sol[t]
    return Rc, null, s[-2:]/s[0], M
if __name__ == "__main__":
    q = 0.7
    for n in [2, 3]:
        W = weights(n); N = 2*n+2
        out = []
        for x in [1.9, 0.37+0.8j]:
            Rc, null, sv, M = intertwiner(rep(n, q, x), rep(n, q, 1.0), range(n+1), W, W)
            out.append((x, null, sv))
        print(f"n={n}: unknowns (weight-preserving entries) = {M}; null-space dimension at generic x: {[o[1] for o in out]} "
              f"(smallest singular values {[np.round(o[2], 14).tolist() for o in out]})")
