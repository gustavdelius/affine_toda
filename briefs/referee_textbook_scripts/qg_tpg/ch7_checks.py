"""Chapter 7 checks: central element, quasitriangularity, Hermitian form on spin 1.

1. sec-uqghat: is prod_i k_i^{n_i} (Kac labels) or prod_i k_i^{n_i^vee} (dual labels,
   n_i^vee = n_i alpha_i^2/2, tbl-untwisted) central?  Exact exponent bookkeeping
   plus explicit matrices for the c_2^(1) and b_2^(1) vector representations.
2. sec-universal-r: for R = q^{H(x)H/2} sum_n q^{n(n-1)/2}(q-q^-1)^n/[n]! E^n(x)F^n and
   the coproduct eq-coproduct, check on spin 1/2 and spin 1:
   (Delta(x)1)R = R13 R23, (1(x)Delta)R = R13 R12, YBE R12R13R23 = R23R13R12,
   (S(x)1)R = R^-1, (S(x)S)R = R.
3. exr-hermitian-spin1: invariant form <a v, w> = <v, a^* w>, E^*=F, K^*=K^-1, q=e^{i gamma}.
"""
import numpy as np
from fractions import Fraction as Fr
from qaff import qnum, qfact


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


# ---------- 1. central element ----------
# conjugation of e_j by prod_i k_i^{m_i} multiplies it by q^{sum_i m_i d_i a_ij}, d_i = alpha_i^2/2
cases = {
    'c_2^(1)': ([[2, -1, 0], [-2, 2, -2], [0, -1, 2]], [Fr(1), Fr(1, 2), Fr(1)], [1, 2, 1]),
    'b_2^(1)': ([[2, 0, -1], [0, 2, -1], [-2, -2, 2]], [Fr(1), Fr(1), Fr(1, 2)], [1, 1, 2]),
    'g_2^(1)': ([[2, -1, 0], [-1, 2, -1], [0, -3, 2]], [Fr(1), Fr(1), Fr(1, 3)], [1, 2, 3]),
    'd_3^(2)': ([[2, -2, 0], [-1, 2, -1], [0, -2, 2]], [Fr(1, 2), Fr(1), Fr(1, 2)], [1, 1, 1]),
    'a_2^(2)': ([[2, -4], [-1, 2]], [Fr(1, 4), Fr(1)], [2, 1]),
}
for name, (A, d, n) in cases.items():
    r = len(n)
    nv = [n[i]*d[i] for i in range(r)]
    # check n is the null vector: sum_i n_i alpha_i = 0  <=> sum_i n_i d_i a_ij = 0 for all j
    shift = lambda m: [sum(m[i]*d[i]*A[i][j] for i in range(r)) for j in range(r)]
    print(f'{name}: n={n}, n^vee={[str(x) for x in nv]}, q-exponents for prod k^n: {[str(s) for s in shift(n)]}, '
          f'for prod k^(n^vee): {[str(s) for s in shift(nv)]}')
    report(f'{name}: prod k_i^(n_i) central', all(s == 0 for s in shift(n)))
    report(f'{name}: prod k_i^(n_i^vee) NOT central (book text claims it is)', not all(s == 0 for s in shift(nv)) or n == nv)

# explicit c_2^(1) vector rep: compute prod k^n and prod k^(n^vee)
t = 0.83+0.41j; q = t*t
h = [np.diag([-1, 0, 0, 1]), np.diag([1, -1, 1, -1]), np.diag([0, 1, -1, 0])]
qi = [q, t, q]
k = [np.diag(qi[i]**np.diag(h[i]).astype(complex)) for i in range(3)]
Kn = k[0]@k[1]@k[1]@k[2]; Knv = k[0]@k[1]@k[2]
report('c_2^(1) vector rep: k0 k1^2 k2 = 1', np.allclose(Kn, np.eye(4)))
report('c_2^(1) vector rep: k0 k1 k2 is not a scalar', not np.allclose(Knv, Knv[0, 0]*np.eye(4)))
print('   k0 k1 k2 =', np.round(np.diag(Knv), 4), ' (q^{1/2} =', np.round(t, 4), ')')

# ---------- 2. quasitriangularity ----------


def rep(twoj, q):
    dd = twoj+1
    E, F = np.zeros((dd, dd), complex), np.zeros((dd, dd), complex)
    for kk in range(dd-1):
        F[kk+1, kk] = 1
    for kk in range(1, dd):
        E[kk-1, kk] = qnum(kk, q)*qnum(twoj+1-kk, q)
    H = np.array([twoj-2*kk for kk in range(dd)], float)
    return E, F, np.diag(q**H), H


def kr(*ms):
    out = ms[0]
    for m in ms[1:]:
        out = np.kron(out, m)
    return out


def qHH(Ha, Hb, q, sign=1):
    return np.diag([q**(sign*a*b/2) for a in Ha for b in Hb])


def coef(n, q):
    return q**(n*(n-1)/2)*(q-1/q)**n/qfact(n, q)


for twoj in (1, 2):
    E, F, K, H = rep(twoj, q); d = twoj+1; I = np.eye(d); Ki = np.linalg.inv(K)
    mp = np.linalg.matrix_power
    R = qHH(H, H, q) @ sum(coef(n, q)*np.kron(mp(E, n), mp(F, n)) for n in range(d))
    # embeddings into V(x)V(x)V
    D = np.diag
    H1, H2, H3 = [np.kron(np.kron(a, b), c) for a, b, c in ((D(H), I, I), (I, D(H), I), (I, I, D(H)))]
    expo = lambda X, Y: np.diag(q**(np.diag(X)*np.diag(Y)/2))
    E1, E2 = kr(E, I, I), kr(I, E, I); F2, F3 = kr(I, F, I), kr(I, I, F)
    K2 = kr(I, K, I)
    R12 = kr(R, I)
    R23 = kr(I, R)
    # R13 = sum over terms: q^{H1 H3/2} sum c_n E1^n F3^n
    E1n = lambda n: mp(E1, n)
    R13 = expo(H1, H3) @ sum(coef(n, q)*E1n(n)@mp(F3, n) for n in range(d))
    # (Delta(x)1)R = q^{(H1+H2)H3/2} sum c_n (E1 K2 + E2)^n F3^n
    DE = kr(E, K, I) + kr(I, E, I)
    lhs1 = expo(H1+H2, H3) @ sum(coef(n, q)*mp(DE, n)@mp(F3, n) for n in range(2*d))
    report(f'spin {twoj}/2: (Delta(x)1)R = R13 R23', np.allclose(lhs1, R13@R23))
    report(f'spin {twoj}/2: (Delta(x)1)R != R23 R13 (order matters)', not np.allclose(lhs1, R23@R13))
    # (1(x)Delta)R = q^{H1(H2+H3)/2} sum c_n E1^n (F2 + K2^-1 F3)^n
    DF = kr(I, F, I) + kr(I, Ki, F)
    lhs2 = expo(H1, H2+H3) @ sum(coef(n, q)*E1n(n)@mp(DF, n) for n in range(2*d))
    report(f'spin {twoj}/2: (1(x)Delta)R = R13 R12', np.allclose(lhs2, R13@R12))
    report(f'spin {twoj}/2: (1(x)Delta)R != R12 R13 (order matters)', not np.allclose(lhs2, R12@R13))
    report(f'spin {twoj}/2: YBE R12 R13 R23 = R23 R13 R12', np.allclose(R12@R13@R23, R23@R13@R12))
    SE = -E@Ki; SF = -K@F
    S1R = sum(coef(n, q)*np.kron(mp(SE, n), I)@qHH(H, H, q, -1)@np.kron(I, mp(F, n)) for n in range(d))
    report(f'spin {twoj}/2: (S(x)1)R = R^-1', np.allclose(S1R, np.linalg.inv(R)))
    SSR = sum(coef(n, q)*np.kron(mp(SE, n), mp(SF, n)) for n in range(d)) @ qHH(H, H, q)
    report(f'spin {twoj}/2: (S(x)S)R = R', np.allclose(SSR, R))

# evaluation rep eq-evaluation-rep: q-Serre (a_01=a_10=-2) for spin 1/2..2
for twoj in (1, 2, 3, 4):
    E, F, K, H = rep(twoj, q); mp = np.linalg.matrix_power; xx = 0.37-0.8j
    e = [xx*F, E]; f = [E/xx, F]
    ok = True
    for (i, j) in ((0, 1), (1, 0)):
        for X in (e, f):
            S = sum((-1)**s_*qfact(3, q)/(qfact(s_, q)*qfact(3-s_, q))*mp(X[i], 3-s_)@X[j]@mp(X[i], s_) for s_ in range(4))
            ok &= np.allclose(S, 0)
    report(f'spin {twoj}/2 evaluation rep satisfies both q-Serre relations', ok)

# ---------- 3. exr-hermitian-spin1 ----------
for gamma in (0.4, 1.2, 2.0, 2.9):
    qq = np.exp(1j*gamma)
    E, F, K, H = rep(2, qq)
    # Gram matrix G with <v,w> = v^dag G w ; invariance <a v,w> = <v,a^* w>  <=> a^dag G = G a^*
    # solve linear equations for Hermitian G
    ops = [(E, F), (F, E), (K, np.linalg.inv(K))]
    rows = []
    for a, astar in ops:
        # a^dag G - G a^* = 0, vectorised (row-major): (a^dag (x) I - I (x) astar^T) vec(G)
        rows.append(np.kron(a.conj().T, np.eye(3)) - np.kron(np.eye(3), astar.T))
    L = np.vstack(rows); u, s, vh = np.linalg.svd(L)
    G = vh[-1].conj().reshape(3, 3); G = G/G[0, 0]
    ok_unique = s[-2] > 1e-6 and s[-1] < 1e-10
    want = np.diag([1, 2*np.cos(gamma), 4*np.cos(gamma)**2])
    ev = np.linalg.eigvalsh((G+G.conj().T)/2)
    sig = (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()))
    report(f'spin-1 form at gamma={gamma}: unique, = diag(1,[2],[2]^2), signature {sig}', ok_unique and np.allclose(G, want))
