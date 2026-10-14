"""Referee check of ch. 35 sec-boundary-krein, "Two and three solitons" (a_2^(1)), and of the
mirror convention for the left wall.  Own construction of the double-row transfer matrices in
braid form; representations and R-check from the book's library boundary_common (allowed).

T_j (particle j of multiplet mu_j at rapidity theta_j):
  move j right through j+1..N with R_{mu_j mu_k}(theta_j - theta_k);
  right wall K_R(theta_j): mu_j(theta_j) -> mubar_j(-theta_j);
  move j left through all others with R_{mu_k mubar_j}(theta_k + theta_j);
  left wall: mubar_j(-theta_j) -> mu_j(theta_j), with the MIRROR convention
     K_L on mubar at -theta = (the right wall's matrix on mu at theta)   [book]
  or the ALTERNATIVE  K_L = right wall's matrix on mubar at theta;
  move j right through 1..j-1 with R_{mu_j mu_k}(theta_j - theta_k).
Walls: U (K = 1) and N (eq-an-k-matrix, s = +1), each matrix part normalized to be Krein-unitary,
K^dag C K = C, C = reversal (prp-boundary-krein, part 2). R-check normalized to 1 on highest weights.
Scalar factors are phases at real rapidity and are dropped.
"""
import os
import sys
import itertools
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'textbook_code'))
import boundary_common as bc   # noqa: E402

pi = np.pi
n = 2; N3 = n + 1
Cmat = np.eye(N3)[::-1]


def rep(mu, om, th):
    q = bc.q_phys(om)
    return bc.an_vector(n, q, np.exp(om * th), conj=(mu == 'W'))


_Rcache = {}


def R(mu, nu, om, t1, t2):
    key = (mu, nu, om, t1, t2)
    if key not in _Rcache:
        Rm, _ = bc.rcheck(rep(mu, om, t1), rep(nu, om, t2))
        _Rcache[key] = Rm
    return _Rcache[key]


def KN(om, th):
    """eq-an-k-matrix, s = +1: V(th) -> Vbar(-th)"""
    q = bc.q_phys(om); s = 1
    X = (-q) ** 0.5 * np.exp(om * th)
    es = 1 / (2 * np.cos(pi * om / 2))
    K = np.zeros((N3, N3), complex)
    for i in range(N3):
        for j in range(N3):
            if j > i:
                K[i, j] = X ** (j - i - N3 / 2)
            elif j < i:
                K[i, j] = s * X ** (j - i + N3 / 2)
            else:
                K[i, i] = ((-q) ** 0.5 * X ** (-N3 / 2) - s * (-q) ** (-0.5) * X ** (N3 / 2)) / (es * (1 / q - q))
    return K


def Kwall(kind, mu, om, th):
    """right-wall matrix on multiplet mu at rapidity th, Krein-normalized"""
    K = np.eye(N3, dtype=complex) if kind == 'U' else KN(om, th)
    if mu == 'W':
        K = Cmat @ K @ Cmat
    cK = (K.conj().T @ Cmat @ K)[0, -1].real / Cmat[0, -1]
    return K / np.sqrt(cK)


def op_on(slot, A, dims):
    mats = [np.eye(d) for d in dims]
    mats[slot] = A
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def op_adj(slot, A, dims):
    """A acts on slots slot, slot+1"""
    left = int(np.prod(dims[:slot])) if slot > 0 else 1
    right = int(np.prod(dims[slot + 2:])) if slot + 2 < len(dims) else 1
    return np.kron(np.kron(np.eye(left), A), np.eye(right))


def transfer(j, mus, ths, om, left, right, convention='mirror'):
    Np = len(mus)
    dims = [N3] * Np
    order = list(range(Np))                         # particle at each slot
    lab = {k: (mus[k], ths[k]) for k in range(Np)}  # current multiplet and rapidity of each particle
    M = np.eye(N3 ** Np, dtype=complex)

    def swap(pos):                                  # particles at pos, pos+1 exchange places
        nonlocal M
        a, b = order[pos], order[pos + 1]
        (ma, ta), (mb, tb) = lab[a], lab[b]
        M = op_adj(pos, R(ma, mb, om, ta, tb), dims) @ M
        order[pos], order[pos + 1] = b, a
    pos = order.index(j)
    while pos < Np - 1:
        swap(pos); pos += 1
    mu, th = lab[j]
    M = op_on(Np - 1, Kwall(right, mu, om, th), dims) @ M
    lab[j] = ('W' if mu == 'V' else 'V', -th)
    while pos > 0:
        swap(pos - 1); pos -= 1
    mub, mth = lab[j]
    orig = mus[j]
    if convention == 'mirror':
        KL = Kwall(left, orig, om, -mth)          # right wall's matrix on mu at theta (= -(-theta))
    else:
        KL = Kwall(left, mub, om, -mth)           # right wall's matrix on mubar at theta
    M = op_on(0, KL, dims) @ M
    lab[j] = (orig, -mth)
    while order.index(j) < j:
        swap(order.index(j))
    assert order == list(range(Np))
    return M


def Zop(mus):
    Z = np.zeros((N3, N3));
    for i in range(N3):
        Z[(i + 1) % N3, i] = 1                     # s = +1 twist
    mats = [Z if m == 'V' else Z for m in mus]
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def sector_eigs(Ts, Zf, zeta=1.0):
    w, U = np.linalg.eig(Zf)
    P = U[:, np.abs(w - zeta) < 1e-8]
    Q, _ = np.linalg.qr(P)
    Tr = [Q.conj().T @ T @ Q for T in Ts]
    rng = np.random.default_rng(0)
    mix = sum(rng.normal() * T for T in Tr)
    ev, W = np.linalg.eig(mix)
    lams = np.array([[np.vdot(W[:, i], T @ W[:, i]) / np.vdot(W[:, i], W[:, i]) for T in Tr] for i in range(len(ev))])
    return lams


def broken(mus, om, left, right, ths):
    Ts = [transfer(j, mus, ths, om, left, right) for j in range(len(mus))]
    lam = sector_eigs(Ts, Zop(mus))
    return np.max(np.abs(np.abs(lam) - 1))


def main():
    ok = True
    om = 0.61
    # sanity: Krein unitarity of R and K, Z commutes
    Rv = R('V', 'W', om, 0.4, -0.3)
    CC = np.kron(Cmat, Cmat)
    print(f"R_VW Krein residual {np.abs(Rv.conj().T @ CC @ Rv - CC).max():.1e};"
          f" K_N Krein residual {np.abs(Kwall('N', 'V', om, 0.37).conj().T @ Cmat @ Kwall('N', 'V', om, 0.37) - Cmat).max():.1e}")
    # mirror convention: commutation of T_1, T_2 (and T_3)
    for mus in (('V', 'V'), ('V', 'W'), ('V', 'V', 'W')):
        ths = [0.83, 0.31, -0.2][:len(mus)]
        for lw, rw in (('U', 'N'), ('N', 'N'), ('N', 'U')):
            res = {}
            for conv in ('mirror', 'other'):
                Ts = [transfer(j, mus, ths, om, lw, rw, conv) for j in range(len(mus))]
                res[conv] = max(np.abs(Ts[a] @ Ts[b] - Ts[b] @ Ts[a]).max() for a in range(len(Ts)) for b in range(a))
            Zf = Zop(mus)
            Ts = [transfer(j, mus, ths, om, lw, rw) for j in range(len(mus))]
            zc = max(np.abs(T @ Zf - Zf @ T).max() for T in Ts)
            CN = Cmat
            for _ in mus[1:]:
                CN = np.kron(CN, Cmat)
            kr = max(np.abs(T.conj().T @ CN @ T - CN).max() for T in Ts)
            print(f"{''.join(mus)} walls {lw}{rw}: [T_i,T_j] mirror {res['mirror']:.1e}, other {res['other']:.1e};"
                  f" [T, Z] {zc:.1e}; Krein {kr:.1e}")
            ok &= res['mirror'] < 1e-9 and res['other'] > 1e-4 and zc < 1e-9 and kr < 1e-9
    # thresholds (two solitons): limit theta_1 large, theta_2 -> 0
    def maxbreak(mus, lw, rw, om):
        return max(broken(mus, om, lw, rw, [t1, t2]) for t1 in (10.0, 12.0) for t2 in (1e-4, 1e-3, 1e-2))

    targets = {(('V', 'V'), 'U', 'N'): 0.698, (('V', 'W'), 'U', 'N'): 0.666,
               (('V', 'V'), 'N', 'N'): 0.536, (('V', 'W'), 'N', 'N'): 0.501}
    for (mus, lw, rw), tgt in targets.items():
        lo, hi = 0.30, 1.0
        assert maxbreak(mus, lw, rw, lo) < 1e-8 and maxbreak(mus, lw, rw, hi) > 1e-8
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            if maxbreak(mus, lw, rw, mid) > 1e-8:
                hi = mid
            else:
                lo = mid
        print(f"threshold {''.join(mus)} {lw}{rw}: omega = {0.5 * (lo + hi):.4f} (book {tgt})")
        ok &= abs(0.5 * (lo + hi) - tgt) < 2e-3
    # below threshold: a coarse full grid is unbroken; UU unbroken at all omega sampled
    grid = [1e-4, 0.05, 0.4, 1.3, 3.0, 7.0, 12.0]
    for mus in (('V', 'V'), ('V', 'W')):
        for om_ in (0.15, 0.45, 0.9, 1.3):
            b = max(broken(mus, om_, 'U', 'U', [t1, t2]) for t1 in grid for t2 in grid if t1 != t2)
            ok &= b < 1e-8
        print(f"UU {''.join(mus)}: unbroken on the grid at omega = 0.15, 0.45, 0.9, 1.3")
        b = max(broken(mus, 0.45, 'U', 'N', [t1, t2]) for t1 in grid for t2 in grid if t1 != t2)
        print(f"UN {''.join(mus)} at omega = 0.45 (below threshold), full grid: max ||Lambda|-1| = {b:.1e}")
        ok &= b < 1e-8
    # three solitons
    for lw, rw in (('U', 'U'), ('U', 'N'), ('N', 'N')):
        out = []
        for om_ in (0.2, 0.5, 0.9):
            b = max(broken(('V', 'V', 'V'), om_, lw, rw, list(t)) for t in ((0.9, 0.4, 0.1), (2.0, 0.7, 0.05), (5.0, 1.1, 0.3)))
            out.append(b)
        print(f"three solitons VVV walls {lw}{rw}: max ||Lambda|-1| at omega 0.2, 0.5, 0.9 = {np.round(out, 4)}")
        ok &= min(out) > 1e-6
    print("PASS" if ok else "FAIL", "boundary Krein: mirror convention, two-soliton thresholds, three solitons broken (a_2)")


if __name__ == '__main__':
    main()
