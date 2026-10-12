"""Chapter 40 (sec-q-operators): defect operators on the circle, Q-operators, the TQ relation and the twist limit.

Ids checked: prp-defect-transfer, cnj-defect-tq, prp-defect-twist-limit (and prp-twist-defect of ch. 36 as the limit).

Setting: N sine-Gordon solitons with rapidities theta_j (an inhomogeneous spin-1/2 chain, y_j = e^{lambda theta_j}),
lambda = 1.37, q = q_phys, Part VI conventions.  R_aj(y/y_j) = P Rcheck (boundary_common.rcheck), the twisted bulk
transfer matrix t_w(y) = Tr_a[w^{h_a} R_aN ... R_a1], and the defect operator with defect space W,
      D_w(y_W) = Tr_W[ w^{h_W} L_N^-1 ... L_1^-1 ],     L_j = borel_l_operator(sign, y_j, W) on V_j (x) W,
so that X = L_aW is the intertwiner of the train argument (X R_aj L_j^-1 = L_j^-1 R_aj X is eq-triangle-rtt).
Fock modules (U_q(b^+), BLZ) give Q-operators; Weyl modules (U_q(b^-), type I_-) give the defect operators.

Checks:
 1. Fock, |w| = 0.8: the twisted trace over the window n = 0..K is stable under K = 40 -> 80 (difference ~ |w|^2K).
 2. Fock: [Q(y_F), Q(y_F')] = 0, [Q(y_F), t_w(y)] = 0 for two y; controls: the twist w' != w, and the product
    L_N ... L_1 without inverses, do not commute with t_w.
 3. TQ relation (Fock), no fitted constants:
      t_w(i y_F) Q_s(y_F) = q^{N/2} Q_{s+1}(q y_F) + prod_j q^{-1/2} (q z_j^2 + q^-1)/(z_j^2 + 1) Q_{s-1}(y_F/q),  z_j = y_j/y_F;
    control: the shifts q^{+-2} fail.  (N = 4.)
 4. Weyl (b^-): in the sector of soliton charge m = sum sigma_j, <n|M|n>_W = q^{mn} <0|M|0>_W (M = the monodromy), so
    the window trace sum_{|n|<=K} w^{2n+s} <n|M|n> is a sector-dependent scalar times D_W := <0|M|0>; for |w| != 1
    this scalar grows without bound as K increases, the normalized operator does not change.
 5. Weyl: [D_W(y_d), D_W(y_d')] = 0; in the sector m, D_W commutes with t_w for w = q^{-m/2} (the support of the
    formal sum sum_n (w^2 q^m)^n), and not with the untwisted t for m != 0.
 6. cnj-defect-tq (Weyl), no fitted constants, every sector m with w_m = q^{-m/2}:
      t_{w_m}(i y_d) D_W(y_d) = q^{-N/2} D_W(q y_d) + prod_j q^{1/2} (q z_j^2 + q^-1)/(z_j^2 + 1) D_W(y_d/q),  z_j = y_j/y_d;
    checked for N = 2, 3, 4; for N = 4 the family {D_W(y_d)} spans 3 dimensions in the sector m = 0 (dimension 6), so
    the test is not empty (for N <= 3 it is weaker); control: the shifts q^{+-2} fail.
 7. prp-defect-twist-limit: (a) the full T (Phi included, eq-sg-transmission-matrix) tends to T_inf:
    v_+-|n> -> q^{s/2} q^{-+(2n+s)/2} |n> v_+- as theta - theta_d -> -infinity, with error ~ e^{min(1,lambda)(theta-theta_d)};
    (b) the defect Bethe-Yang transfer matrix T_1 = L_1 S_1N ... S_12 of sec-defect-circle tends, on each |n>, to
    q^{s/2} times the Bloch-twisted bulk T_1 of prp-twist-defect with e^{i alpha_n} = q^{-(2n+s)/2};
    (c) as theta - theta_d -> +infinity, T/(Phi Y) tends to the charge-conjugating limit v_+|n> -> -i q^-n |n+1> v_-,
    v_-|n> -> -i q^{n+s} |n-1> v_+; (d) D_W tends to q^{ms/2}, a function of the charge alone.
Windows: Fock n = 0..40 (80 for check 1), Weyl n = -6..6 (D_W) and -4..4 (T_1).  Runtime about 10 seconds.
"""
import numpy as np
from boundary_common import report, q_phys, sl2_spin, rcheck, prop_residual
from defect_common import fock_module_b, borel_l_operator, flip_VW, sg_tmatrix, tmatrix_from_elements

LAM = 1.37
q = q_phys(LAM)
qh = np.sqrt(q + 0j)
kap = q - 1/q
kron = np.kron
P22 = flip_VW(2, 2)


def Rmat(z):
    """R = P Rcheck on V(y1) (x) V(y2), z = y1/y2, normalized to 1 on v_+ (x) v_+."""
    return P22 @ rcheck(sl2_spin(1, q, z), sl2_spin(1, q, 1.0))[0]


_perm = {}


def embed(L, j, N, dW, dV=2):
    """L on V_j (x) W  ->  operator on V^{(x)N} (x) W (W last); dense, for small spaces."""
    key = (j, N, dW)
    if key not in _perm:
        shape = [dV]*N + [dW]; Dt = dV**N*dW
        idx = np.arange(Dt).reshape(shape)
        order = [k for k in range(N) if k != j] + [j, N]
        P = np.zeros((Dt, Dt)); P[np.arange(Dt), idx.transpose(order).ravel()] = 1
        _perm[key] = P
    P = _perm[key]
    return P.T @ kron(np.eye(dV**(N - 1)), L) @ P


def apply_local(L, j, N, dW, M):
    """(L acting on V_j (x) W) @ M, for M with rows indexed by V^{(x)N} (x) W (W last), without forming the embedding."""
    X = M.shape[1]
    T = np.moveaxis(M.reshape([2]*N + [dW, X]), [j, N], [0, 1])
    T = np.tensordot(L.reshape(2, dW, 2, dW), T, axes=([2, 3], [0, 1]))
    return np.moveaxis(T, [0, 1], [j, N]).reshape(M.shape)


def ptrace(M, D, dW, w):
    return np.einsum('awbw,w->ab', M.reshape(D, dW, D, dW), w)


def bulk_t(y, ys, w):
    """t_w(y) = Tr_a[w^{h_a} R_aN(y/y_N) ... R_a1(y/y_1)] on V^{(x)N}."""
    N = len(ys); D = 2**N
    M = np.eye(2*D, dtype=complex)
    for j in range(N):
        M = apply_local(P22 @ Rmat(y/ys[j]) @ P22.T, j, N, 2, M)      # R_aj with a stored last
    return ptrace(M, D, 2, np.array([w, 1/w]))


def monodromy(sign, kind, yW, ys, s=0.0, K=40, inverse=True):
    N = len(ys); D = 2**N
    W = fock_module_b(sign, q, yW, K, s, kind); dW = W['dim']
    M = np.eye(D*dW, dtype=complex)
    for j in range(N):
        L = borel_l_operator(sign, ys[j], W)
        M = apply_local(np.linalg.inv(L) if inverse else L, j, N, dW, M)
    return M, W


def Q_fock(yF, ys, w, s=0.0, K=40, inverse=True):
    M, W = monodromy(1, 'fock', yF, ys, s, K, inverse)
    return ptrace(M, 2**len(ys), W['dim'], w**W['h'])


def D_weyl(yd, ys, s=0.0, K=6):
    M, W = monodromy(-1, 'weyl', yd, ys, s, K)
    D, dW = 2**len(ys), W['dim']
    i0 = list(W['n']).index(0)
    return M.reshape(D, dW, D, dW)[:, i0, :, i0], M, W


def charges(N):
    return np.array([sum(1 if (k >> (N - 1 - b)) & 1 == 0 else -1 for b in range(N)) for k in range(2**N)])


def comm(A, B):
    return np.abs(A@B - B@A).max()/(np.abs(A).max()*np.abs(B).max())


ys4 = np.exp(LAM*np.array([0.31, -0.42, 0.77, 0.05]))
ys3 = ys4[:3]
mag4 = charges(4)
w = 0.8*np.exp(0.4j)
yF = 1.3*np.exp(0.15j)

# ---------------------------------------------------------------- 1. Fock: convergence of the twisted trace
Q40, Q80 = Q_fock(yF, ys4, w, K=40), Q_fock(yF, ys4, w, K=80)
dK = np.abs(Q40 - Q80).max()/np.abs(Q80).max()
report(f'Fock Q-operator (N = 4, |w| = 0.8): twisted trace over n = 0..K stable under K = 40 -> 80 [{dK:.0e}, '
       f'|w|^80 = {abs(w)**80:.0e}]', dK < 1e-6)
w = 0.5*np.exp(0.4j)                      # checks 2-3: truncation error |w|^80 ~ 1e-24
Q40 = Q_fock(yF, ys4, w, K=40)

# ---------------------------------------------------------------- 2. Fock: commutation
Q2 = Q_fock(0.6*np.exp(-0.5j), ys4, w)
c_QQ = comm(Q40, Q2)
c_Qt = max(comm(Q40, bulk_t(y, ys4, w)) for y in (0.7*np.exp(0.2j), 2.3))
c_bad1 = comm(Q40, bulk_t(0.7*np.exp(0.2j), ys4, 0.5))
c_bad2 = comm(Q_fock(yF, ys4, w, inverse=False), bulk_t(0.7*np.exp(0.2j), ys4, w))
report(f'prp-defect-transfer (Fock, Q-operators): [Q(y_F), Q(y_F\')] = 0 [{c_QQ:.0e}], [Q, t_w(y)] = 0 [{c_Qt:.0e}]; '
       f'controls: other twist [{c_bad1:.2f}], L_N...L_1 without inverses [{c_bad2:.2f}]',
       c_QQ < 1e-12 and c_Qt < 1e-12 and c_bad1 > 1e-3 and c_bad2 > 1e-3)


# ---------------------------------------------------------------- 3. TQ relation (Fock), predicted coefficients
def tq_fock(shift, ys, s=0.0):
    z = ys/yF
    a = qh**len(ys)
    d = np.prod((q*z**2 + 1/q)/(qh*(z**2 + 1)))
    lhs = bulk_t(1j*yF, ys, w) @ Q_fock(yF, ys, w, s)
    rhs = a*Q_fock(shift*yF, ys, w, s + 1) + d*Q_fock(yF/shift, ys, w, s - 1)
    return np.abs(lhs - rhs).max()/np.abs(lhs).max()


r_tq = tq_fock(q, ys4)
r_tq_bad = tq_fock(q**2, ys4)
report(f'TQ relation for the Fock Q-operator (N = 4): t_w(i y_F) Q_s(y_F) = q^(N/2) Q_s+1(q y_F) + prod_j q^(-1/2)(q z_j^2 + '
       f'q^-1)/(z_j^2 + 1) Q_s-1(y_F/q) with no fitted constant [{r_tq:.0e}]; control: shifts q^(+-2) [{r_tq_bad:.2f}]',
       r_tq < 1e-10 and r_tq_bad > 1e-2)

# ---------------------------------------------------------------- 4. Weyl: factorization of the trace
yd = 1.3*np.exp(0.1j)
S = 0.3
DW, M, Wm = D_weyl(yd, ys4, S)
D16, dW = 16, Wm['dim']
M4 = M.reshape(D16, dW, D16, dW)
fac = 0
for k, n in enumerate(Wm['n']):
    if abs(n) <= 3:
        fac = max(fac, np.abs(M4[:, k, :, k] - DW*q**(np.outer(np.ones(D16), mag4)*n)).max())
w2 = 0.8
tr = {}
for Kw in (2, 3):
    sel = np.abs(Wm['n']) <= Kw
    tr[Kw] = sum(w2**(2*n + S)*M4[:, k, :, k] for k, n in enumerate(Wm['n']) if sel[k])
m0 = mag4 == 2
grow = np.abs(tr[3][np.ix_(m0, m0)]).max()/np.abs(tr[2][np.ix_(m0, m0)]).max()
norm_dev = prop_residual(tr[3][np.ix_(m0, m0)], DW[np.ix_(m0, m0)])[0]
report(f'Weyl defect operator: in the sector of charge m, <n|M|n> = q^(mn) <0|M|0> [{fac:.0e}]; the window trace with '
       f'|w| = 0.8 changes by a factor {grow:.2f} from |n| <= 2 to |n| <= 3 (no convergence), while the normalized operator '
       f'is D_W = <0|M|0> in every window [{norm_dev:.0e}]', fac < 1e-12 and grow > 1.3 and norm_dev < 1e-12)

# ---------------------------------------------------------------- 5. Weyl: commutation
DW2 = D_weyl(0.6*np.exp(-0.3j), ys4, S)[0]
c_DD = comm(DW, DW2)
c_Dt, c_un = 0, 0
for m in sorted(set(mag4)):
    sel = mag4 == m
    for y in (0.7*np.exp(0.2j), 2.3):
        A = DW[np.ix_(sel, sel)]
        c_Dt = max(c_Dt, comm(A, bulk_t(y, ys4, qh**(-m))[np.ix_(sel, sel)]))
        if m != 0 and sel.sum() > 1:
            c_un = max(c_un, comm(A, bulk_t(y, ys4, 1.0)[np.ix_(sel, sel)]))
report(f'prp-defect-transfer (Weyl, type I_-): [D_W(y_d), D_W(y_d\')] = 0 [{c_DD:.0e}]; in the sector m, D_W commutes with '
       f't_w for w = q^(-m/2) [{c_Dt:.0e}], not with the untwisted t when m != 0 [{c_un:.2f}]',
       c_DD < 1e-12 and c_Dt < 1e-12 and c_un > 1e-3)


# ---------------------------------------------------------------- 6. cnj-defect-tq (Weyl), predicted coefficients
def tq_weyl(shift, ys, s=S):
    N = len(ys); mg = charges(N)
    z = ys/yd
    a = qh**(-N)
    d = np.prod(qh*(q*z**2 + 1/q)/(z**2 + 1))
    D0, Dp, Dm = D_weyl(yd, ys, s)[0], D_weyl(shift*yd, ys, s)[0], D_weyl(yd/shift, ys, s)[0]
    worst = 0
    for m in sorted(set(mg)):
        sel = mg == m
        t = bulk_t(1j*yd, ys, qh**(-m))
        lhs = (t @ D0)[np.ix_(sel, sel)]
        rhs = a*Dp[np.ix_(sel, sel)] + d*Dm[np.ix_(sel, sel)]
        worst = max(worst, np.abs(lhs - rhs).max()/np.abs(lhs).max())
    return worst


fam = np.array([D_weyl(yd*np.exp(0.3*k + 0.2j*k), ys4, S)[0][np.ix_(mag4 == 0, mag4 == 0)].ravel() for k in range(5)])
sv = np.linalg.svd(fam, compute_uv=False)
rank0 = int((sv/sv[0] > 1e-10).sum())
r_w, r_w_bad = max(tq_weyl(q, ys4[:k]) for k in (2, 3, 4)), tq_weyl(q**2, ys4)
report(f'cnj-defect-tq (Weyl, N = 2, 3, 4): in every sector m, t_w(i y_d) D_W(y_d) = q^(-N/2) D_W(q y_d) + prod_j q^(1/2)'
       f'(q z_j^2 + q^-1)/(z_j^2 + 1) D_W(y_d/q), w = q^(-m/2), no fitted constant [{r_w:.0e}]; family rank {rank0} in the '
       f'6-dimensional sector m = 0; control: shifts q^(+-2) [{r_w_bad:.2f}]', r_w < 1e-10 and rank0 == 3 and r_w_bad > 1e-2)

# ---------------------------------------------------------------- 7. prp-defect-twist-limit
TD, NWL = 0.31, 4


def T_inf(W, s):
    V_ = sl2_spin(1, q, 1.0)
    return tmatrix_from_elements(lambda n: q**-n, lambda n: 0, lambda n: 0, lambda n: q**(n + s), V_, W)


errs = []
for s_ in (0.0, S):
    for Dl in (-3.0, -6.0):
        T, V_, W_ = sg_tmatrix(TD + Dl, LAM, q, TD, NWL, s=s_)
        errs.append(np.abs(T - T_inf(W_, s_)).max())
rate = errs[0]/errs[1]
report(f'prp-defect-twist-limit (a): T(theta) -> diag(q^-n, q^(n+s)) = q^(s/2) q^(-+(2n+s)/2) as theta - theta_d -> -inf: '
       f'errors {errs[0]:.1e}, {errs[1]:.1e} at theta - theta_d = -3, -6 (ratio {rate:.1f}, e^(3 min(1, lambda)) = '
       f'{np.exp(3*min(1, LAM)):.1f})', errs[1] < 1e-2 and abs(rate/np.exp(3*min(1, LAM)) - 1) < 0.2)

# (b) defect Bethe-Yang transfer matrix T_1 = L_1 S_1N ... S_12 on V^{(x)3} (x) W
th3 = np.log(ys3)/LAM


def T1_defect(theta_d, s_, limit=False):
    N = 3; D = 8
    T, V_, W_ = sg_tmatrix(th3[0], LAM, q, theta_d, NWL, s=s_)
    if limit:
        T = T_inf(W_, s_)
    dW = W_['dim']
    L1 = flip_VW(2, dW).T @ T                                       # L = P T on V_1 (x) W
    M = np.eye(D*dW, dtype=complex)
    for k in range(1, N):                                           # S_12 first, then ..., S_1N
        Sk = Rmat(ys3[0]/ys3[k])
        # S_1k on V_1 (x) V_k inside V^3 (x) W
        idx = np.arange(D*dW).reshape([2, 2, 2, dW])
        order = [m for m in range(N) if m not in (0, k)] + [0, k, N]
        P = np.zeros((D*dW, D*dW)); P[np.arange(D*dW), idx.transpose(order).ravel()] = 1
        oth = [m for m in range(N) if m not in (0, k)]
        big = kron(np.eye(2**len(oth)), kron(Sk, np.eye(dW)))
        M = (P.T @ big @ P) @ M
    return embed(L1, 0, N, dW) @ M, W_


errsT = []
for Dd in (3.0, 6.0):
    A, W_ = T1_defect(th3[0] + Dd, S)
    B, _ = T1_defect(th3[0] + Dd, S, limit=True)
    sel = np.tile(np.abs(W_['n']) <= NWL - 1, 8)
    errsT.append(np.abs((A - B)[np.ix_(sel, sel)]).max()/np.abs(B[np.ix_(sel, sel)]).max())
# the limit is the Bloch twist of prp-twist-defect on each |n>: q^{s/2} e^{i alpha_n h_1} T_1^bulk
B, W_ = T1_defect(10.0, S, limit=True)
dW = W_['dim']
T1b = np.eye(8, dtype=complex)
for k in (1, 2):
    idx = np.arange(8).reshape([2, 2, 2]); oth = [m for m in range(3) if m not in (0, k)]
    order = oth + [0, k]
    P = np.zeros((8, 8)); P[np.arange(8), idx.transpose(order).ravel()] = 1
    T1b = (P.T @ kron(np.eye(2), Rmat(ys3[0]/ys3[k])) @ P) @ T1b
h1 = np.repeat([1.0, -1.0], 4)
B4 = B.reshape(8, dW, 8, dW)
dev_tw = 0
for kk, n in enumerate(W_['n']):
    tw = q**(S/2)*np.diag(qh**(-(2*n + S)*h1)) @ T1b
    dev_tw = max(dev_tw, np.abs(B4[:, kk, :, kk] - tw).max())
report(f'prp-defect-twist-limit (b): the defect Bethe-Yang transfer matrix T_1 = L_1 S_13 S_12 (sec-defect-circle) converges '
       f'[{errsT[0]:.0e} -> {errsT[1]:.0e} for theta_d - theta_1 = 3 -> 6] to q^(s/2) e^(i alpha_n h_1) T_1^bulk on each |n>, '
       f'e^(i alpha_n) = q^(-(2n+s)/2): the twist defect of prp-twist-defect [{dev_tw:.0e}]',
       errsT[1] < errsT[0]/10 and dev_tw < 1e-12)

# (c) theta - theta_d -> +inf: charge conjugation dressed by the shift
errC = []
for Dl in (3.0, 6.0):
    T, V_, W_ = sg_tmatrix(TD + Dl, LAM, q, TD, NWL, s=S)
    Tc = tmatrix_from_elements(lambda n: 0, lambda n: -1j*q**(n + S), lambda n: -1j*q**-n, lambda n: 0, V_, W_)
    errC.append(prop_residual(T, Tc)[0])
# (d) the Weyl defect operator in the limit: <0|M|0> with L^-1 -> diag  =>  q^{ms/2}
ys_far = ys4*np.exp(-LAM*8.0)
DWf = D_weyl(yd, ys_far, S)[0]
dev_d = np.abs(DWf - np.diag(q**(mag4*S/2))).max()
report(f'prp-defect-twist-limit (c), (d): as theta - theta_d -> +inf, T -> charge conjugation with the shift |n> -> |n+-1> '
       f'[{errC[0]:.0e} -> {errC[1]:.0e}]; as theta_j - theta_d -> -inf the Weyl defect operator D_W -> q^(ms/2), a function '
       f'of the charge alone [{dev_d:.0e}]', errC[1] < errC[0]/10 and dev_d < 1e-4)
