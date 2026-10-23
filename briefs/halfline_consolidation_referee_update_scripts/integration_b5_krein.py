import os
"""B5: two a_2 solitons on an interval, Krein structure (sec-boundary-krein, ch. 35).  Own transfer matrices; book helpers
used: boundary_common.an_vector, rcheck, q_phys (R-check normalized to 1 on highest weights).
Walls: U (K = 1) and N (eq-an-k-matrix, s = +1, divided by sqrt(c_K), c_K = (cosh 2T theta - cos pi T)/(2 sin^2 pi omega/2));
on Vbar the matrix is K^T.  Mirror rule: the left wall's matrix on tbar at -theta is the right wall's matrix on t at theta.
Sector zeta = 1 of Z (x) Z, Z the cyclic shift; deviation = max ||Lambda| - 1| over eigenvalues of T_1, T_2 there.
"""
import sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code'))
from boundary_common import an_vector, rcheck, q_phys

n, d = 2, 3


def Kneu(om, th):
    N = n + 1
    X = np.exp(om*(th - 0.5j*np.pi))
    w = np.exp(-0.5j*np.pi*om)
    q = q_phys(om)
    eps = 1/(2*np.cos(np.pi*om/2))
    K = np.empty((N, N), complex)
    for i in range(N):
        for j in range(N):
            if i == j:
                K[i, j] = (w*X**(-N/2) - X**(N/2)/w)/(eps*(1/q - q))
            elif j > i:
                K[i, j] = X**(j - i - N/2)
            else:
                K[i, j] = X**(j - i + N/2)
    T = N*om/2
    cK = (np.cosh(2*T*th) - np.cos(np.pi*T))/(2*np.sin(np.pi*om/2)**2)
    return K/np.sqrt(cK)


def Kwall(code, t, om, th):
    if code == 'U':
        return np.eye(d, dtype=complex)
    K = Kneu(om, th)
    return K if t == 'V' else K.T


def R(t1, x1, t2, x2, om):
    q = q_phys(om)
    return rcheck(an_vector(n, q, np.exp(om*x1), t1 == 'W'), an_vector(n, q, np.exp(om*x2), t2 == 'W'))[0]


def T_two(j, types, ths, left, right, om, mirror=True):
    """Double-row transfer matrix of particle j (0 or 1) for two particles, ordered left to right."""
    bar = {'V': 'W', 'W': 'V'}
    I = np.eye(d)
    st = [(types[0], ths[0]), (types[1], ths[1])]
    M = np.eye(d*d, dtype=complex)
    if j == 0:
        # particle 0 crosses particle 1 (moving right)
        M = R(st[0][0], st[0][1], st[1][0], st[1][1], om) @ M
        st = [st[1], st[0]]
        t, x = st[1]
        M = np.kron(I, Kwall(right, t, om, x)) @ M
        st[1] = (bar[t], -x)
        M = R(st[0][0], st[0][1], st[1][0], st[1][1], om) @ M
        st = [st[1], st[0]]
        t, x = st[0]
        M = np.kron(Kwall(left, bar[t] if mirror else t, om, -x), I) @ M
        st[0] = (bar[t], -x)
    else:
        t, x = st[1]
        M = np.kron(I, Kwall(right, t, om, x)) @ M
        st[1] = (bar[t], -x)
        M = R(st[0][0], st[0][1], st[1][0], st[1][1], om) @ M
        st = [st[1], st[0]]
        t, x = st[0]
        M = np.kron(Kwall(left, bar[t] if mirror else t, om, -x), I) @ M
        st[0] = (bar[t], -x)
        M = R(st[0][0], st[0][1], st[1][0], st[1][1], om) @ M
        st = [st[1], st[0]]
    return M


Z1 = np.roll(np.eye(3), 1, axis=0)
ZZ = np.kron(Z1, Z1)
Pz = (np.eye(9) + ZZ + ZZ@ZZ)/3
Uz = np.linalg.svd(Pz)[0][:, :3]          # orthonormal basis of the zeta = 1 sector


def deviation(types, ths, walls, om):
    Ts = [T_two(j, types, ths, walls[0], walls[1], om) for j in (0, 1)]
    out = 0
    for T in Ts:
        Ts_ = Uz.conj().T @ T @ Uz
        out = max(out, np.abs(np.abs(np.linalg.eigvals(Ts_)) - 1).max())
    return out


if __name__ == '__main__':
    # (a) commutativity and the mirror rule
    for walls in ('UN', 'NN', 'NU'):
        for types in ('VV', 'VW'):
            om = 0.37
            ths = [0.7, 0.25]
            for mirror in (True, False):
                T0 = T_two(0, list(types), ths, walls[0], walls[1], om, mirror)
                T1 = T_two(1, list(types), ths, walls[0], walls[1], om, mirror)
                c = np.linalg.norm(T0@T1 - T1@T0)
                zc = max(np.linalg.norm(ZZ@T - T@ZZ) for T in (T0, T1))
                print(f'(a) walls {walls} {types} mirror={mirror}: ||[T1,T2]|| = {c:.1e}, ||[ZZ,T]|| = {zc:.1e}')

    # (b) VV thresholds at (theta_1, theta_2) = (12, 1e-3)
    def bisect(f, lo, hi, it=18):
        for _ in range(it):
            mid = (lo + hi)/2
            if f(mid):
                hi = mid
            else:
                lo = mid
        return (lo + hi)/2
    for walls in ('UN', 'NN'):
        thr = bisect(lambda om: deviation(['V', 'V'], [12.0, 1e-3], walls, om) > 1e-6, 0.4, 0.9)
        print(f'(b) VV {walls}: threshold at theta = (12, 1e-3): {thr:.5f}')

    # (c) V Vbar: grid of slow rapidities
    grid = [(t1, t2) for t1 in (12.0, 20.0) for t2 in (1e-3, 0.3, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0)]
    for walls in ('UN', 'NN'):
        def broken(om, pts):
            return max(max(deviation(['V', 'W'], [a, b], walls, om), deviation(['V', 'W'], [b, a], walls, om)) for a, b in pts) > 1e-6
        thr_grid = bisect(lambda om: broken(om, grid), 0.4, 0.9, 14)
        thr_zero = bisect(lambda om: broken(om, [(12.0, 1e-3)]), 0.4, 0.9, 14)
        # where does the onset occur just above the grid threshold?
        om = thr_grid + 0.003
        devs = {(a, b): deviation(['V', 'W'], [a, b], walls, om) for a, b in grid}
        worst = max(devs, key=devs.get)
        print(f'(c) V Vbar {walls}: threshold over the grid {thr_grid:.4f} (largest deviation at {thr_grid + 0.003:.4f}: '
              f'theta = {worst}, {devs[worst]:.3f}); with theta_2 -> 0 only: {thr_zero:.4f}')
