"""Krein structure on the half-line for a_n^(1) vector solitons (Part VII, ch. 35, sec-boundary-krein; half-line open
problem 7(a), (b); and the soliton-antisoliton case of prp-krein-fundamental, ch. 40).

Book conventions: physical q = -exp(-i pi omega) (|q| = 1), principal gradation y = exp(omega theta), V and Vbar =
an_vector, R-check = rcheck (normalized to 1 on highest weights), coideal b_j = e_j + c_j f_j k_j + eps_j k_j with
c_j = q^{-1}. C is the colour reversal v_i -> v_{n+2-i} on V and on Vbar; nu is the diagram reflection j -> n+1-j
(fixing node 0); tau is the antilinear anti-automorphism tau(e_j) = e_(nu j), tau(f_j) = f_(nu j), tau(k_j) = k_(nu j).
Checks (omega = 0.37, so that k_closed's branch is the analytic one):
 1. For real y: pi(a)^dagger = C pi(tau a) C on V(y) and Vbar(y), for every generator.
 2. For real rapidities: R-check_ab^dagger (C x C) R-check_ab = C x C for ab = VV, VVbar, VbarV, VbarVbar (n = 2, 3).
 3. tau(b_j) = b_(nu j) with eps_j replaced by conj(eps_j); so B_eps is tau-stable iff eps_(nu j) = conj(eps_j).
 4. The s = +-1 matrices (eq-an-k-matrix): K^dagger C K = c_K C with c_K = (cosh 2T theta - s cos pi T)/(2 sin^2(pi omega/2))
    >= 0, T = (n+1) omega/2, and C K C = K^T (n = 2, 3).
 5. One soliton between the uniform wall (K = 1) and the Neumann wall (s = +1), n = 2: K/sqrt(c_K) (the modulus of the
    physical scalar factor is 1/sqrt(c_K)) has one eigenvalue of modulus 1 (the Z_3 sector zeta = 1) and a pair
    (rho, 1/rho), rho != 1 (zeta = exp(+-2 pi i/3)), with rho = |y zeta - w|/|y zeta - conj w|, y = exp(omega theta), w = exp(-i pi omega/2).
 6. Z_(n+1) sectors: [Z, K] = 0, [Z x Z, R-check] = 0, C Z C = Z^(-1).
 7. Two and three solitons of a_2 on an interval (open problem 7(b)): the double-row transfer matrices commute (and
    fail to with the unmirrored left wall); in the sector zeta = 1 the levels are unimodular for uniform walls on both
    sides, and for UN and NN walls below the coupling thresholds (VV: 0.69-0.71 UN, 0.53-0.54 NN; VW: 0.65-0.68 UN,
    0.49-0.51 NN), not above; three solitons are broken.
Runtime: about 2 minutes.
"""
import itertools
import numpy as np
from boundary_common import an_vector, rcheck, q_phys, coideal_gen, report

OM = 0.37
Q = q_phys(OM)


def k_closed(n, y, s, qq=Q):
    """eq-an-k-matrix (as in boundary_coideal.k_closed; valid on the analytic branch for omega < 1)."""
    N = n + 1
    w = np.sqrt(-qq + 0j)
    es = 1/np.sqrt((1 - qq)*(1 - 1/qq))
    X = w*y
    K = np.zeros((N, N), complex)
    for i in range(N):
        for j in range(N):
            K[i, j] = ((w*X**(-N/2) - s/w*X**(N/2))/(es*(1/qq - qq)) if i == j
                       else X**(j - i - N/2) if j > i else s*X**(j - i + N/2))
    return K


# ---------------------------------------------------------------- 1. the antilinear anti-automorphism tau
dev1 = 0
for n in (2, 3, 4):
    N = n + 1
    C = np.fliplr(np.eye(N))
    nu = lambda j: (N - j) % N
    for conj in (False, True):
        rep = an_vector(n, Q, np.exp(OM*0.7), conj=conj)
        for j in range(N):
            for kind in 'efk':
                dev1 = max(dev1, np.linalg.norm(rep[(kind, j)].conj().T - C @ rep[(kind, nu(j))] @ C))
report(f'pi(a)^dagger = C pi(tau a) C on V(y), Vbar(y), y real, all generators, n = 2, 3, 4 [{dev1:.0e}]', dev1 < 1e-14)

# ---------------------------------------------------------------- 2. Krein-unitarity of every R-check
dev2 = 0
for n in (2, 3):
    N = n + 1
    CC = np.kron(np.fliplr(np.eye(N)), np.fliplr(np.eye(N)))
    for a in (False, True):
        for b in (False, True):
            for th in (0.4, 1.2):
                R = rcheck(an_vector(n, Q, np.exp(OM*th), conj=a), an_vector(n, Q, 1.0, conj=b))[0]
                dev2 = max(dev2, np.linalg.norm(R.conj().T @ CC @ R - CC))
report(f'R-check^dagger (C x C) R-check = C x C for VV, VVbar, VbarV, VbarVbar at real rapidity, untwisted C, '
       f'n = 2, 3 [{dev2:.0e}]', dev2 < 1e-12)

# ---------------------------------------------------------------- 3. tau-stability of the coideal
rng = np.random.default_rng(1)
dev3, bad3 = 0, np.inf
for n in (2, 3):
    N = n + 1
    C = np.fliplr(np.eye(N))
    nu = lambda j: (N - j) % N
    rep = an_vector(n, Q, np.exp(OM*0.7))
    eps = rng.normal(size=N) + 1j*rng.normal(size=N)
    eps_real = eps.copy()
    for j in range(N):                       # impose eps_(nu j) = conj(eps_j)
        if nu(j) == j:
            eps_real[j] = eps[j].real
        elif j < nu(j):
            eps_real[nu(j)] = np.conj(eps_real[j])
    for e, out in ((eps_real, 'ok'), (eps, 'bad')):
        for j in range(N):
            lhs = C @ coideal_gen(rep, j, 1/Q, e[j]).conj().T @ C
            rhs = coideal_gen(rep, nu(j), 1/Q, np.conj(e[j]))      # tau(b_j) = b_(nu j)(conj eps_j)
            dev3 = max(dev3, np.linalg.norm(lhs - rhs))
        if out == 'bad':
            bad3 = min(bad3, max(abs(e[nu(j)] - np.conj(e[j])) for j in range(N)))
report(f'C pi(b_j)^dagger C = pi(b_(nu j)) with eps -> conj(eps): tau(b_j) = b_(nu j)(conj eps_j) [{dev3:.0e}]; '
       f'B_eps is tau-stable iff eps_(nu j) = conj(eps_j)', dev3 < 1e-14 and bad3 > 1e-3)

# ---------------------------------------------------------------- 4. K^dagger C K = c_K C for s = +-1
dev4 = 0
for n in (2, 3):
    N, T = n + 1, (n + 1)*OM/2
    C = np.fliplr(np.eye(N))
    for s in (1, -1):
        for th in (0.3, 0.9, -0.6):
            K = k_closed(n, np.exp(OM*th), s)
            cK = (np.cosh(2*T*th) - s*np.cos(np.pi*T))/(2*np.sin(np.pi*OM/2)**2)
            dev4 = max(dev4, np.linalg.norm(K.conj().T @ C @ K - cK*C)/abs(cK), np.linalg.norm(C @ K @ C - K.T),
                       0 if cK > 0 else 1)
report(f'K^dagger C K = c_K C, c_K = (cosh 2T theta - s cos pi T)/(2 sin^2(pi omega/2)) > 0, and C K C = K^T '
       f'(s = +-1, n = 2, 3) [{dev4:.0e}]', dev4 < 1e-13)

# ---------------------------------------------------------------- 5. one soliton between uniform and Neumann walls
n, N, T = 2, 3, 3*OM/2
th = 0.55
K = k_closed(n, np.exp(OM*th), 1)
cK = (np.cosh(2*T*th) - np.cos(np.pi*T))/(2*np.sin(np.pi*OM/2)**2)
mods = np.sort(np.abs(np.linalg.eigvals(K/np.sqrt(cK))))
w = np.exp(-1j*np.pi*OM/2)
yy = np.exp(OM*th)
rho = abs(yy*np.exp(2j*np.pi/3) - w)/abs(yy*np.exp(2j*np.pi/3) - np.conj(w))
ok5 = abs(mods[1] - 1) < 1e-12 and abs(mods[0]*mods[2] - 1) < 1e-12 and abs(mods[2] - 1) > 1e-3 \
    and min(abs(mods[2] - rho), abs(mods[2] - 1/rho)) < 1e-12
report(f'one soliton, uniform/Neumann walls, a_2: eigenvalue moduli of K/sqrt(c_K) = {mods.round(6).tolist()} '
       f'(1 and a pair r, 1/r = rho^(+-1), rho = {rho:.6f})', ok5)

# ---------------------------------------------------------------- 6. the Z_(n+1) sectors
# Z: twisted cyclic shift, Z v_i = v_(i+1), Z v_n = s v_0, Z^(n+1) = s. K (s = +-1) commutes with Z, R-check commutes
# with Z x Z, and C Z C = Z^(-1): the Krein form pairs the Z-sector zeta with conj(zeta) = 1/zeta.
dev6 = 0
for n in (2, 3):
    N = n + 1
    C = np.fliplr(np.eye(N))
    for s in (1, -1):
        Z = np.zeros((N, N), complex)
        for i in range(N - 1):
            Z[i + 1, i] = 1
        Z[0, N - 1] = s
        K = k_closed(n, np.exp(OM*0.45), s)
        dev6 = max(dev6, np.linalg.norm(Z @ K - K @ Z), np.linalg.norm(C @ Z @ C - np.linalg.inv(Z)))
    Z1 = np.roll(np.eye(N), 1, axis=0)
    R = rcheck(an_vector(n, Q, np.exp(OM*0.4)), an_vector(n, Q, 1.0))[0]
    dev6 = max(dev6, np.linalg.norm(np.kron(Z1, Z1) @ R - R @ np.kron(Z1, Z1)))
report(f'Z_(n+1) sectors: [Z, K] = 0 with the twisted shift (Z^(n+1) = s), [Z x Z, R-check] = 0, C Z C = Z^(-1) '
       f'(n = 2, 3) [{dev6:.0e}]', dev6 < 1e-12)


# ---------------------------------------------------------------- 7. two and three solitons on the interval (a_2)
# Double-row transfer matrices (sec-boundary-bethe-yang): particle j crosses j+1..N-1, reflects at the right wall,
# crosses everybody, reflects at the left wall (the mirror image: its matrix on tbar at -theta is the right wall's
# matrix on t at theta) and returns. Walls: U (K = 1) and N (s = +1, K/sqrt(c_K), on the analytic branch, so valid
# for omega > 1 too); bulk R-check normalized on highest weights. The scalar phases are unimodular at real rapidity
# and are omitted. Deviation = max ||Lambda| - 1| over the joint eigenvalues of the T_j in the sector zeta = 1.
def kmat_an(n, om, th, s):
    N = n + 1
    L = om*(th - 0.5j*np.pi)
    w = np.exp(-0.5j*np.pi*om)
    q = q_phys(om)
    es = 1/(2*np.cos(np.pi*om/2))
    K = np.empty((N, N), complex)
    for i in range(N):
        for j in range(N):
            K[i, j] = ((w*np.exp(-N/2*L) - s/w*np.exp(N/2*L))/(es*(1/q - q)) if i == j
                       else np.exp((j - i - N/2)*L) if j > i else s*np.exp((j - i + N/2)*L))
    return K


def wall(code, typ, om, th, n=2):
    if code == 'U':
        return np.eye(n + 1, dtype=complex)
    T = (n + 1)*om/2
    cK = (np.cosh(2*T*th) - np.cos(np.pi*T))/(2*np.sin(np.pi*om/2)**2)
    K = kmat_an(n, om, th, 1)/np.sqrt(cK)
    return K if typ == 'V' else K.T


def transfer(j, types, ths, left, right, om, mirror=True, n=2):
    d, N = n + 1, len(types)
    bar = {'V': 'W', 'W': 'V'}
    rep = lambda t, x: an_vector(n, q_phys(om), np.exp(om*x), t == 'W')
    op = lambda M, p, k: np.kron(np.kron(np.eye(d**p), M), np.eye(d**(N - p - k)))
    st = [(types[k], ths[k]) for k in range(N)]
    M = np.eye(d**N, dtype=complex)

    def swap(p):
        nonlocal M
        M = op(rcheck(rep(*st[p]), rep(*st[p + 1]))[0], p, 2) @ M
        st[p], st[p + 1] = st[p + 1], st[p]
    p = j
    while p < N - 1:
        swap(p); p += 1
    t, x = st[p]
    M = op(wall(right, t, om, x), p, 1) @ M
    st[p] = (bar[t], -x)
    while p > 0:
        swap(p - 1); p -= 1
    t, x = st[p]
    M = op(wall(left, bar[t] if mirror else t, om, -x), p, 1) @ M
    st[p] = (bar[t], -x)
    while p < j:
        swap(p); p += 1
    return M


rng7 = np.random.default_rng(7)


def sector(types, ths, walls, om, mirror=True):
    Ts = [transfer(j, types, ths, walls[0], walls[1], om, mirror) for j in range(len(types))]
    Z1 = np.roll(np.eye(3), 1, axis=0).astype(complex)
    Z = Z1
    for _ in types[1:]:
        Z = np.kron(Z, Z1)
    comm = max(np.linalg.norm(Ts[0] @ T - T @ Ts[0]) + np.linalg.norm(Z @ T - T @ Z) for T in Ts)
    _, P = np.linalg.eig(sum(rng7.normal()*T for T in Ts) + 0.37*Z)
    Pi = np.linalg.inv(P)
    sel = np.abs(np.diag(Pi @ Z @ P) - 1) < 1e-6
    return max(np.abs(np.abs(np.diag(Pi @ T @ P)[sel]) - 1).max() for T in Ts), comm


t2s = np.geomspace(1e-4, 3, 8)
dev = lambda types, walls, om: max(sector(list(types), [12.0, t2], walls, om)[0] for t2 in t2s)
comm7 = max(sector(list(t), [0.7, 0.25], w, om)[1] for t in ('VV', 'VW') for w in ('UN', 'NN') for om in (0.37, 0.83))
bad7 = sector(['V', 'V'], [0.7, 0.25], 'NN', 0.37, mirror=False)[1]
report(f'a_2 double-row transfer matrices: [T_1, T_2] = 0 and [Z x Z, T_j] = 0 for UN, NN walls, VV and VW '
       f'[{comm7:.0e}]; with the unmirrored left wall [T_1, T_2] = {bad7:.2f} (control)', comm7 < 1e-12 and bad7 > 0.1)
uu = max(dev(t, 'UU', om) for t in ('VV', 'VW') for om in (0.37, 0.83, 1.43))
rows, ok7 = [], uu < 1e-6
for t, w, lo, hi in (('VV', 'UN', 0.69, 0.71), ('VW', 'UN', 0.65, 0.68), ('VV', 'NN', 0.53, 0.54), ('VW', 'NN', 0.49, 0.51)):
    a, b = dev(t, w, lo), dev(t, w, hi)
    ok7 &= a < 1e-6 and b > 1e-2
    rows.append(f'{w} {t}: {a:.0e} at {lo}, {b:.2f} at {hi}')
report(f'two solitons, sector zeta = 1, theta_1 = 12, theta_2 in [1e-4, 3]: UU unbroken at omega = 0.37, 0.83, 1.43 '
       f'[{uu:.0e}]; thresholds ' + '; '.join(rows), ok7)
g3 = [8.0, 2.0, 0.5, 0.05]
vvv = min(max(sector(['V']*3, list(c), w, om)[0] for c in itertools.combinations(g3, 3))
          for w in ('UU', 'UN', 'NN') for om in (0.3, 0.83))
report(f'three solitons VVV, sector zeta = 1: broken for UU, UN, NN at omega = 0.3, 0.83 (min over cases of the max '
       f'deviation {vvv:.2f})', vvv > 0.1)
