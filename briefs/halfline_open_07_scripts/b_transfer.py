"""Problem 7(b): double-row (boundary Bethe-Yang) transfer matrices of N fundamental a_n solitons on an interval
with soliton-conjugating boundaries.  Boundary codes: 'U' uniform (eps = 0, K = 1), 'N' Neumann (s = +1), 'M' (s = -1).
Right wall: K^R_V(theta) = A K(theta): V(y) -> Vbar(1/y); K^R_Vbar = A K(theta)^T (report of problem 3).
Left wall = image of the right wall under the parity map omega (e <-> f, k -> k^-1), which maps V(y) to Vbar(1/y)
in the code bases: K^L_{abar}(-theta) = K^R_a(theta) as matrices ('omega' mirror).  The alternative 'swap'
(K^L_{abar}(-theta) = K^R_{abar}(theta)) is kept only to test it.
Scalars: bulk hw-normalized R-checks (c_R = 1, scalar phases omitted), boundary A_min (|A_min|^2 cK = 1)."""
import numpy as np
from p07_common import Rc, k_book, A_min, Cmat
bar = {'V': 'W', 'W': 'V'}

def Kmat(code, a, theta, n, om):
    N = n + 1
    if code == 'U':
        return np.eye(N, dtype=complex)
    s = 1 if code == 'N' else -1
    K = A_min(n, theta, s, om)*k_book(n, theta, s, om)
    return K if a == 'V' else K.T

def embed2(M, p, N, d):          # two-site matrix on positions p, p+1 of d^N
    return np.kron(np.kron(np.eye(d**p), M), np.eye(d**(N - p - 2)))

def embed1(M, p, N, d):
    return np.kron(np.kron(np.eye(d**p), M), np.eye(d**(N - p - 1)))

def Tdr(j, types, th, left, right, n, om, mirror='omega'):
    """Round trip of particle j (0-based).  types[k] in {'V','W'}; th real or complex rapidities."""
    d = n + 1; Np = len(types); y = np.exp(om*np.asarray(th, complex))
    M = np.eye(d**Np, dtype=complex)
    order = list(range(Np))                 # particle at each position
    yy = {k: y[k] for k in range(Np)}; tt = dict(enumerate(types))
    def swap(p):                             # R-check at positions p, p+1
        nonlocal M
        a, b = order[p], order[p + 1]
        M = embed2(Rc(n, tt[a], tt[b], yy[a]/yy[b], om), p, Np, d) @ M
        order[p], order[p + 1] = b, a
    p = j
    while p < Np - 1:                        # 1. move right through j+1..N-1
        swap(p); p += 1
    M = embed1(Kmat(right, tt[j], th[j], n, om), Np - 1, Np, d) @ M    # 2. right wall
    tt[j] = bar[tt[j]]; yy[j] = 1/y[j]
    while p > 0:                             # 3. move left through everybody
        swap(p - 1); p -= 1
    a_in = tt[j]                             # 4. left wall: incoming abar at -theta
    KL = Kmat(left, bar[a_in], th[j], n, om) if mirror == 'omega' else Kmat(left, a_in, th[j], n, om)
    M = embed1(KL, 0, Np, d) @ M
    tt[j] = bar[a_in]; yy[j] = y[j]
    while p < j:                             # 5. move right through 0..j-1
        swap(p); p += 1
    assert order == list(range(Np)) and tt[j] == types[j]
    return M

def Zshift(n, s=1):
    N = n + 1; Z = np.zeros((N, N), complex)
    for i in range(N):
        Z[(i + 1) % N, i] = 1
    Z[0, N - 1] = s
    return Z

def kron_all(Ms):
    out = np.eye(1)
    for M in Ms:
        out = np.kron(out, M)
    return out
