"""One-loop particle masses of the folded (real-coupling) Toda theory, in the same normalization as the
soliton computation: V_real = (m^2/beta^2) sum_j n_j (e^{beta alpha_j.phi} - 1) of the PARENT, restricted to
the sigma-invariant subspace. Bubble diagram as in soliton_rmatrix_code/smatrix/oneloop_ratio.py
(tadpoles vanish by normal ordering). Returns masses and d log m_a / d beta^2 (real coupling)."""
import numpy as np
from scipy.integrate import quad
from toda import affine, perm_matrix

def parent_roots(alg, r):
    C, n, _ = affine(alg, r)
    G = C[1:, 1:]
    L = np.linalg.cholesky(G)
    roots = np.zeros((r+1, r))
    roots[1:] = L
    roots[0] = -sum(n[i]*roots[i] for i in range(1, r+1))/n[0]
    assert np.allclose(roots @ roots.T, C)
    return roots, n

def J(ma2, mb, mc):
    return quad(lambda x: 1/(x*mb**2 + (1 - x)*mc**2 - x*(1 - x)*ma2), 0, 1, limit=200)[0]/(4*np.pi)

def folded_particles(alg, r, perm=None):
    roots, n = parent_roots(alg, r)
    N = r + 1
    if perm is None:
        W = np.eye(r)
    else:
        # induced linear map on root space: alpha_j -> alpha_{perm[j]}
        A = roots[1:]                         # basis
        B = np.array([roots[perm[j]] for j in range(1, N)])
        S = np.linalg.solve(A, B).T           # S alpha_j = alpha_{perm j}  (columns)
        S = (B.T @ np.linalg.inv(A.T))
        assert np.allclose(np.array([S @ roots[j] for j in range(N)]), roots[perm]), "not an automorphism"
        # group average projector
        G = [np.eye(r)]; g = S.copy()
        while not np.allclose(g, np.eye(r)):
            G.append(g); g = g @ S
        Pi = sum(G)/len(G)
        u, s, vt = np.linalg.svd(Pi)
        W = u[:, s > 0.5]                     # orthonormal basis of invariant subspace
    Ae = roots @ W                            # projected roots in an orthonormal basis of W
    M2 = sum(n[j]*np.outer(Ae[j], Ae[j]) for j in range(N))
    m2, V = np.linalg.eigh(M2)
    Ae = Ae @ V
    C3 = np.einsum('i,ia,ib,ic->abc', n, Ae, Ae, Ae)
    m = np.sqrt(m2); d = len(m2)
    # self-energy matrix within each degenerate mass block, then diagonalize
    dm2 = np.zeros(d)
    a = 0
    while a < d:
        blk = [a]
        while blk[-1] + 1 < d and abs(m2[blk[-1] + 1] - m2[a]) < 1e-9: blk.append(blk[-1] + 1)
        S = np.zeros((len(blk), len(blk)))
        for i, ai in enumerate(blk):
            for k, ak in enumerate(blk):
                S[i, k] = -0.5*sum(C3[ai, b, c]*C3[ak, b, c]*J(m2[a], m[b], m[c]) for b in range(d) for c in range(d)
                                   if abs(C3[ai, b, c]) > 1e-9 or abs(C3[ak, b, c]) > 1e-9)
        ev = np.linalg.eigvalsh((S + S.T)/2)
        dm2[blk] = np.sort(ev)
        a = blk[-1] + 1
    return m, dm2/(2*m2)        # masses, d log m_a / d beta^2 (real coupling, m and beta = 1 units)

if __name__ == "__main__":
    pi = np.pi
    for r in [2, 3, 4, 5]:
        m, dl = folded_particles('a', r)
        h = r + 1
        print(f"a_{r}: masses", np.round(m, 5), " dlog m/dbeta^2:", np.round(dl, 6), " -(1/8h)cot(pi/h) =", round(-1/(8*h)/np.tan(pi/h), 6))
    m, dl = folded_particles('d', 5); print("d_5:", np.round(m, 5), np.round(dl, 6), round(-1/(8*8)/np.tan(pi/8), 6))
    m, dl = folded_particles('e', 6); print("e_6:", np.round(m, 5), np.round(dl, 6), round(-1/(8*12)/np.tan(pi/12), 6))
