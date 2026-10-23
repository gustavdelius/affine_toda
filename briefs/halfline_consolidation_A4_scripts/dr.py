"""Independent double-row transfer matrices for N a_2 solitons on an interval (book conventions only)."""
import sys
sys.path.insert(0, '/home/gustav/Git/affine_toda/code')
import numpy as np
from boundary_common import q_phys, an_vector, rcheck

n, d = 2, 3


def kmat(om, th, s):                       # eq-an-k-matrix, analytic branch (an_nondiagonal_scalar.kmat)
    N = n + 1
    L = om*(th - 0.5j*np.pi)
    w = np.exp(-0.5j*np.pi*om)
    q = q_phys(om)
    eps = 1/(2*np.cos(np.pi*om/2))
    K = np.empty((N, N), complex)
    for i in range(N):
        for j in range(N):
            K[i, j] = ((w*np.exp(-N/2*L) - s/w*np.exp(N/2*L))/(eps*(1/q - q)) if i == j
                       else np.exp((j - i - N/2)*L) if j > i else s*np.exp((j - i + N/2)*L))
    return K


def wall(code, typ, om, th):
    """Physical right-wall amplitude for an incoming particle of type typ ('V' or 'W') at real rapidity th,
    up to a unimodular scalar: K/sqrt(c_K) on V, its transpose on Vbar; 1 for the uniform wall."""
    if code == 'U':
        return np.eye(d, dtype=complex)
    s = 1 if code == 'N' else -1
    T = (n + 1)*om/2
    cK = (np.cosh(2*T*th) - s*np.cos(np.pi*T))/(2*np.sin(np.pi*om/2)**2)
    K = kmat(om, th, s)/np.sqrt(cK)
    return K if typ == 'V' else K.T


def rep(typ, om, th):
    return an_vector(n, q_phys(om), np.exp(om*th), typ == 'W')


def op2(M, p, N):
    return np.kron(np.kron(np.eye(d**p), M), np.eye(d**(N - p - 2)))


def op1(M, p, N):
    return np.kron(np.kron(np.eye(d**p), M), np.eye(d**(N - p - 1)))


bar = {'V': 'W', 'W': 'V'}


def transfer(j, types, ths, left, right, om):
    N = len(types)
    state = [(types[k], ths[k], k) for k in range(N)]     # (type, rapidity, label) at each position
    M = np.eye(d**N, dtype=complex)

    def swap(p):
        nonlocal M
        (ta, xa, la), (tb, xb, lb) = state[p], state[p + 1]
        M = op2(rcheck(rep(ta, om, xa), rep(tb, om, xb))[0], p, N) @ M
        state[p], state[p + 1] = state[p + 1], state[p]

    p = j
    while p < N - 1:
        swap(p); p += 1
    t, x, l = state[p]
    M = op1(wall(right, t, om, x), p, N) @ M            # right wall: t(x) -> tbar(-x)
    state[p] = (bar[t], -x, l)
    while p > 0:
        swap(p - 1); p -= 1
    t, x, l = state[p]                                   # incoming t at -theta_j; mirror of the right wall:
    M = op1(wall(left, bar[t], om, -x), p, N) @ M       # matrix on t(-theta) = right-wall matrix on tbar(theta)
    state[p] = (bar[t], -x, l)
    while p < j:
        swap(p); p += 1
    assert [s[2] for s in state] == list(range(N)) and state[j][0] == types[j] and abs(state[j][1] - ths[j]) < 1e-14
    return M


def zshift(s):
    Z = np.zeros((d, d), complex)
    for i in range(d):
        Z[(i + 1) % d, i] = 1
    Z[0, d - 1] = s
    return Z


rng = np.random.default_rng(7)


def sector_dev(types, ths, left, right, om, zeta=1):
    """max ||Lambda| - 1| over the joint eigenvalues of the T_j in the Z-sector zeta (Z = twisted shift on V, its
    inverse transpose... chosen so that it commutes; returned with the commutator residuals)."""
    N = len(types)
    s = -1 if 'M' in (left + right) else 1
    Ts = [transfer(j, types, ths, left, right, om) for j in range(N)]
    Z1 = zshift(s)
    zs = {'V': Z1, 'W': Z1}
    Z = np.eye(1)
    for t in types:
        Z = np.kron(Z, zs[t])
    comm = max([np.linalg.norm(Ts[0] @ T - T @ Ts[0]) for T in Ts] + [np.linalg.norm(Z @ T - T @ Z) for T in Ts])
    G = sum(rng.normal()*T for T in Ts) + 0.37*Z
    _, P = np.linalg.eig(G)
    Pi = np.linalg.inv(P)
    zd = np.diag(Pi @ Z @ P)
    sel = np.abs(zd - zeta*s**N) < 1e-6
    dev = max(np.abs(np.abs(np.diag(Pi @ T @ P)[sel]) - 1).max() for T in Ts)
    return dev, comm, sel.sum()
