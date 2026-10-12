"""Helpers for the s = +-1 scalar-factor computation (scratch).  Book conventions (boundary_common.py)."""
import numpy as np
import scipy.sparse as sps
from boundary_common import q_phys, an_vector, nodes_of


def weights(rep):
    """Diagonal of log_q k_j for each node -> tuple weight label per basis vector (using k matrices)."""
    d = rep['dim']
    ks = [np.diag(rep[('k', j)]) for j in nodes_of(rep)]
    # label by rounded complex values of the k eigenvalues
    return [tuple(k[i] for k in ks) for i in range(d)]


def _cop_sparse(A, B, kind, j):
    IA, IB = sps.identity(A['dim'], format='csr'), sps.identity(B['dim'], format='csr')
    c = sps.csr_matrix
    if kind == 'e':
        return sps.kron(c(A[('e', j)]), IB) + sps.kron(c(A[('k', j)]), c(B[('e', j)]))
    if kind == 'f':
        return sps.kron(c(A[('f', j)]), c(np.linalg.inv(B[('k', j)]))) + sps.kron(IA, c(B[('f', j)]))
    return sps.kron(c(A[('k', j)]), c(B[('k', j)]))


def fast_rcheck(A, B):
    """R-check A(x)B -> B(x)A, weight-restricted dense solve; normalized to 1 on hw(x)hw."""
    dA, dB = A['dim'], B['dim']
    d = dA*dB
    wA, wB = weights(A), weights(B)
    rd = lambda v: tuple(np.round(np.array(v), 6) + 0)
    win = [rd(np.array(wA[a])*np.array(wB[b])) for a in range(dA) for b in range(dB)]
    wout = [rd(np.array(wB[b])*np.array(wA[a])) for b in range(dB) for a in range(dA)]
    cols = [(o, i) for i in range(d) for o in range(d) if wout[o] == win[i]]
    colidx = [o + d*i for (o, i) in cols]          # column-major vec index of X[o, i]
    blocks = []
    Id = sps.identity(d, format='csr')
    for j in nodes_of(A):
        for kind in 'ef':
            P = _cop_sparse(A, B, kind, j)          # on A(x)B
            Q = _cop_sparse(B, A, kind, j)          # on B(x)A
            M = sps.kron(P.T, Id) - sps.kron(Id, Q)
            blocks.append(M.tocsc()[:, colidx])
    M = sps.vstack(blocks).toarray()
    _, s, Vh = np.linalg.svd(M, full_matrices=False)
    v = Vh[-1].conj()
    assert s[-1]/s[0] < 1e-9 and s[-2]/s[0] > 1e-7, (s[-1]/s[0], s[-2]/s[0])
    X = np.zeros((d, d), complex)
    for (o, i), val in zip(cols, v):
        X[o, i] = val
    i_in = A['hw']*dB + B['hw']
    i_out = B['hw']*dA + A['hw']
    return X/X[i_out, i_in]


def k_closed(n, y, s, qq):
    N = n + 1
    w = np.sqrt(-qq + 0j)
    es = 1/np.sqrt((1 - qq)*(1 - 1/qq))
    X = w*y
    K = np.zeros((N, N), complex)
    for i in range(N):
        for j in range(N):
            if i == j:
                K[i, j] = (w*X**(-N/2) - s/w*X**(N/2))/(es*(1/qq - qq))
            elif j > i:
                K[i, j] = X**(j - i - N/2)
            else:
                K[i, j] = s*X**(j - i + N/2)
    return K


def k_book(n, om, th, s):
    """eq-an-k-matrix with the book's branch X^p = exp(p om (th - i pi/2)), w = (-q)^(1/2) = exp(-i pi om/2),
    eps_* (q^-1 - q) = w - w^-1 (eps_* = 1/(w + w^-1)).  Rapidity th (complex)."""
    N = n + 1
    L = om*(th - 0.5j*np.pi)
    Xp = lambda p: np.exp(p*L)
    w = np.exp(-0.5j*np.pi*om)
    K = np.zeros((N, N), complex)
    for i in range(N):
        for j in range(N):
            if i == j:
                K[i, j] = (w*Xp(-N/2) - s/w*Xp(N/2))/(w - 1/w)
            elif j > i:
                K[i, j] = Xp(j - i - N/2)
            else:
                K[i, j] = s*Xp(j - i + N/2)
    return K
