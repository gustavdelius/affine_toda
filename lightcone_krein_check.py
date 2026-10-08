"""
Light-cone (brickwork) regularization of the U_q(sl_n^) vertex model:
checks of the Krein structure, discrete symmetries, complex quasi-energies,
and the Euclidean phase problem.  Companion to the discussion of open problem 8.

Conventions
  Jimbo R-matrix, homogeneous gradation, q = e^{i mu} (|q| = 1, critical regime),
  multiplicative spectral parameter x real  (Minkowski / light-cone).
  Gate  G(x) = P R(x) / (x q - 1/(x q)).
  Light-cone cutoff rapidity theta_cut enters as x = exp(mu * n * theta_cut / pi)
  (Zinn-Justin's rescaling Theta = mu h theta / 2pi, gate argument 2 Theta).
  Brickwork  U = U_even U_odd  on L sites (L even, periodic).
Doc normalization (a_{n-1}^{(1)}):  beta^2 = 4 (pi - mu).
"""
import itertools
import numpy as np
import scipy.sparse as sp


def gate(n, x, q, strip_phase=False):
    d = n * n
    R = np.zeros((d, d), complex)
    for i in range(n):
        for j in range(n):
            if i == j:
                R[i*n+i, i*n+i] = x*q - 1/(x*q)
            else:
                R[i*n+j, i*n+j] = x - 1/x
                R[i*n+j, j*n+i] = (q - 1/q) * (abs(x) if strip_phase else (x if i < j else 1/x))
    P = np.zeros((d, d))
    for i in range(n):
        for j in range(n):
            P[i*n+j, j*n+i] = 1
    return P @ R / (x*q - 1/(x*q)), P


def site_perm(n, L, f):
    """permutation operator moving the content of site j to site f(j)"""
    D = n**L
    rows, cols = [], []
    for c in itertools.product(range(n), repeat=L):
        c2 = [0]*L
        for j in range(L):
            c2[f(j)] = c[j]
        rows.append(np.ravel_multi_index(c2, [n]*L)); cols.append(np.ravel_multi_index(c, [n]*L))
    return sp.csr_matrix((np.ones(D), (rows, cols)), shape=(D, D))


def bond(G, j, k, n, L):
    Gt = G.reshape(n, n, n, n); rows, cols, vals = [], [], []
    for c in itertools.product(range(n), repeat=L):
        a = np.ravel_multi_index(c, [n]*L)
        for o1 in range(n):
            for o2 in range(n):
                v = Gt[o1, o2, c[j], c[k]]
                if v != 0:
                    c2 = list(c); c2[j], c2[k] = o1, o2
                    rows.append(np.ravel_multi_index(c2, [n]*L)); cols.append(a); vals.append(v)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n**L,)*2)


def brickwork(n, L, x, q, strip_phase=False):
    G, _ = gate(n, x, q, strip_phase)
    Uo = Ue = sp.identity(n**L, dtype=complex, format="csr")
    for j in range(0, L, 2):
        Uo = bond(G, j, j+1, n, L) @ Uo
    for j in range(1, L, 2):
        Ue = bond(G, j, (j+1) % L, n, L) @ Ue
    return (Ue @ Uo).tocsr()


def sectors(n, L):
    s = {}
    for c in itertools.product(range(n), repeat=L):
        s.setdefault(tuple(np.bincount(c, minlength=n)), []).append(np.ravel_multi_index(c, [n]*L))
    return {k: np.array(v) for k, v in s.items()}


if __name__ == "__main__":
    mu = 1.0; q = np.exp(1j*mu)

    # 1. gate identities (any n): G^T = G,  G(x)G(1/x) = 1,  conj G(x) = S G(x)^{-1} S
    for n in (2, 3, 4):
        G, S = gate(n, np.exp(0.8), q); Gi = np.linalg.inv(G)
        print(f"n={n}: |G^T-G|={np.abs(G.T-G).max():.0e}  "
              f"|G(x)G(1/x)-1|={np.abs(G@gate(n, np.exp(-0.8), q)[0]-np.eye(n*n)).max():.0e}  "
              f"|conj G - S G^-1 S|={np.abs(G.conj()-S@Gi@S).max():.0e}")

    # 2. brickwork: Krein unitarity and symmetries (sl_3, L=6)
    n, L = 3, 6
    U = brickwork(n, L, np.exp(mu*n*1.0/np.pi), q).toarray()
    Pi = site_perm(n, L, lambda j: L-1-j).toarray()        # bond-centred reflection
    Pi_s = site_perm(n, L, lambda j: (-j) % L).toarray()     # site-centred reflection
    C = np.zeros_like(Pi)
    for c in itertools.product(range(n), repeat=L):
        C[np.ravel_multi_index(tuple(n-1-t for t in c), [n]*L), np.ravel_multi_index(c, [n]*L)] = 1
    Ui = np.linalg.inv(U)
    print("Krein:      |U^dag Pi U - Pi|            =", f"{np.abs(U.conj().T@Pi@U-Pi).max():.1e}")
    print("antilinear: |Pi_s conj(U) Pi_s - U^-1|   =", f"{np.abs(Pi_s@U.conj()@Pi_s-Ui).max():.1e}")
    print("unitary:    |[C Pi, U]| =", f"{np.abs(C@Pi@U-U@C@Pi).max():.1e}",
          "  |[Pi, U]| =", f"{np.abs(Pi@U-U@Pi).max():.2f}", "  |[C, U]| =", f"{np.abs(C@U-U@C).max():.2f}")

    # 3. complex quasi-energies (|lambda| != 1) by sector
    for n, L in [(2, 8), (3, 6), (3, 8)]:
        secs = sectors(n, L)
        reps = sorted({tuple(sorted(k, reverse=True)) for k in secs}, reverse=True)
        for m, th in [(1.0, 2.0), (2.4, 2.0)]:
            Uf = brickwork(n, L, np.exp(m*n*th/np.pi), np.exp(1j*m))
            out = []
            for r in reps:
                ev = np.linalg.eigvals(Uf[secs[r]][:, secs[r]].toarray())
                out.append(f"{tuple(int(t) for t in r)}:{int((np.abs(np.abs(ev)-1) > 1e-6).sum())}/{len(ev)}")
            print(f"sl_{n} L={L} mu={m} (beta^2/2pi={4*(np.pi-m)/(2*np.pi):.2f}) theta_cut={th}: off-circle", " ".join(out))

    # 4. Euclidean phase problem: x = e^{i phi}.  W = same model with the phases e^{+-i phi}
    #    of the exchange weights removed (all entries of W are real and >= 0).
    n, L, phi = 3, 6, (np.pi - mu)/2
    V = brickwork(n, L, np.exp(1j*phi), q).toarray()
    W = brickwork(n, L, np.exp(1j*phi), q, strip_phase=True).toarray()
    print(f"Euclidean sl_3, L={L}: W real and >= 0: {np.abs(W.imag).max() < 1e-12 and W.real.min() > -1e-12}")
    ev = np.linalg.eigvals(V); lead = ev[np.argmax(np.abs(ev))]
    print(f"  leading eigenvalue of V: {lead:.4f}")
    for m in (1, 2, 4, 8):
        a, b = np.trace(np.linalg.matrix_power(V, m)), np.trace(np.linalg.matrix_power(W, m))
        print(f"  m={m}: Tr V^m = {a.real:.4g}   average phase Tr V^m / Tr W^m = {(a/b).real:.3f}")
