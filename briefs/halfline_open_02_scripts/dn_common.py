"""Common helpers for open problem 2 (d_n vector soliton reflection), book conventions.
Imports the book's boundary_common.py unchanged from the repository (read only)."""
import sys
sys.path.insert(0, '/home/gustav/Git/affine_toda/textbook_code')
import numpy as np
import scipy.sparse as sps
from boundary_common import (q_phys, dn_vector, nodes_of, cop, kron, k_intertwiner, intertwiners, dual,
                             prop_residual, _mk)


def weights(rep):
    ks = [np.diag(rep[('k', j)]) for j in nodes_of(rep)]
    return [tuple(k[i] for k in ks) for i in range(rep['dim'])]


def _cop_sparse(A, B, kind, j):
    IA, IB = sps.identity(A['dim'], format='csr'), sps.identity(B['dim'], format='csr')
    c = sps.csr_matrix
    if kind == 'e':
        return sps.kron(c(A[('e', j)]), IB) + sps.kron(c(A[('k', j)]), c(B[('e', j)]))
    if kind == 'f':
        return sps.kron(c(A[('f', j)]), c(np.linalg.inv(B[('k', j)]))) + sps.kron(IA, c(B[('f', j)]))
    return sps.kron(c(A[('k', j)]), c(B[('k', j)]))


def fast_rcheck(A, B, norm=True):
    """R-check A(x)B -> B(x)A (Part VI coproduct), weight-restricted; normalized to 1 on hw(x)hw."""
    dA, dB = A['dim'], B['dim']
    d = dA*dB
    wA, wB = weights(A), weights(B)
    rd = lambda v: tuple(np.round(np.array(v), 6) + 0)
    win = [rd(np.array(wA[a])*np.array(wB[b])) for a in range(dA) for b in range(dB)]
    wout = [rd(np.array(wB[b])*np.array(wA[a])) for b in range(dB) for a in range(dA)]
    cols = [(o, i) for i in range(d) for o in range(d) if wout[o] == win[i]]
    colidx = [o + d*i for (o, i) in cols]
    blocks = []
    Id = sps.identity(d, format='csr')
    for j in nodes_of(A):
        for kind in 'ef':
            P = _cop_sparse(A, B, kind, j)
            Q = _cop_sparse(B, A, kind, j)
            M = sps.kron(P.T, Id) - sps.kron(Id, Q)
            blocks.append(M.tocsc()[:, colidx])
    M = sps.vstack(blocks)
    G = (M.conj().T @ M).toarray()
    ev, U = np.linalg.eigh(G)
    assert ev[0] < 1e-14*ev[-1] and ev[1] > 1e-9*ev[-1], (ev[:3]/ev[-1])
    v = U[:, 0]
    X = np.zeros((d, d), complex)
    for (o, i), val in zip(cols, v):
        X[o, i] = val
    if norm:
        i_in = A['hw']*dB + B['hw']
        i_out = B['hw']*dA + A['hw']
        X = X/X[i_out, i_in]
    return X


def Kdn(n, q, y, eps):
    """book intertwiner K: V(y) -> V(1/y), coideal c_j = 1/q, eps_j."""
    null, rel = k_intertwiner(dn_vector(n, q, y), dn_vector(n, q, 1/y), [1/q]*(n + 1), eps)
    assert len(null) == 1, len(null)
    return null[0]


def basis_labels(n):
    """index -> signed label i in +-1..+-n"""
    N = 2*n
    return [i + 1 for i in range(n)] + [-(n - k) for k in range(n)]   # index n+k -> v_{-(n-k)}
