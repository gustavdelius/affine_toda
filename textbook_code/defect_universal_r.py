"""Chapter 40 (sec-borel, sec-t-intertwiner, sec-universal-r-affine): q-oscillators, Borel L-operators and the
Khoroshkin-Tolstoy universal R-matrix of U_q(sl2^) in the conventions of Part VI.

Ids checked: eq-q-oscillator, eq-weyl-module, eq-fock-module, eq-borel-intertwiner, thm-t-triangle (hypothesis (ii)),
prp-sg-t-unique, prp-two-borels, eq-l-operator, prp-l-operator.

Conventions: Part VI coproduct Delta(e) = e (x) 1 + k (x) e, Delta(f) = f (x) k^-1 + 1 (x) f; q = q_phys = -e^{-i pi lambda},
lambda = 1.37; y = e^{lambda theta}; kappa = q - q^-1; helpers in defect_common.py (q_oscillator_module, fock_module_b,
borel_l_operator).  The universal R-matrix of Part VI lies in U_q(b^-) (x) U_q(b^+); it is obtained from the
Khoroshkin-Tolstoy formula as written by Boos, Goehmann, Kluemper, Nirov and Razumov (arXiv:1205.1631, eqs.
(4.12)-(4.28)), whose algebra is the book's with q -> 1/q and k_i = q^{h_i}; their R is the book's R_21^-1.

Checks:
 1. eq-q-oscillator on the Weyl and Fock modules (Z = a+ a- - kappa^-2 vanishes on Weyl, Z|n> = -q^2n/kappa^2 on Fock);
    the homomorphisms U_q(b^+-) -> Osc_q of eq-weyl-module/eq-fock-module respect k e k^-1 = q^a e and the q-Serre
    relations on both modules.  A shift module of U_q(b^+) with u_n = e_1 e_0|_n satisfies q-Serre iff
    u_n = A + B q^2n + C q^-2n (random A, B, C pass; u_n = A + B q^4n fails).
 2. KT universal R, truncated at height 3 (root vectors alpha, delta - alpha, alpha + delta, 2 delta - alpha and the
    imaginary root delta), on V(y1) (x) V(y2) equals P Rcheck up to a scalar with error O(x^4), x = y2/y1: halving x
    divides the residual by 16.  Controls: without the imaginary-root factor the error is O(x^2) (ratio 4), without the
    height-3 real roots O(x^3) (ratio 8).  Same for R' = R_21^-1 (expansion in y1/y2).
 3. On V(y) (x) W with W a Weyl module: (pi (x) rho)(R) = scalar * borel_l_operator(+1) (U_q(b^+)) and
    (pi (x) rho)(R_21^-1) = scalar * borel_l_operator(-1) (U_q(b^-)), to O(x^4); the b^- scalar tends to q^-s/2
    relative to the ch. 37 matrix sg_tmatrix (the universal-R normalization).
 4. Fock modules: the U_q(b^+) intertwiner V (x) F -> F (x) V is unique on the window and equals borel_l_operator(+1)
    (BLZ's q-oscillator L-operator); likewise for U_q(b^-); the truncated KT R agrees with both to O(x^4).
 5. eq-triangle-rtt holds for the Fock and Weyl L-operators on a window.
 6. thm-t-triangle, hypothesis (ii): the U_q(b)-intertwiners V(y1) (x) V(y2) (x) W -> W (x) V(y2) (x) V(y1) form a
    one-dimensional space on the window interior (Weyl b^-, Weyl b^+, Fock b^+).
 7. prp-sg-t-unique: the U_q(b^-) Weyl intertwiner is unique, = sg_tmatrix up to a scalar, for s = 0, 0.3, 1;
    no intertwiner exists if k_0 k_1 = 1.7 on W (level zero forced).
 8. prp-two-borels: (a) Weston's L^(r0,0,0) (1006.1555 (2.9)) is a U_q(b^+) intertwiner W (x) V, Delta -> Delta', in
    Part VI conventions; its transpose at zeta = i y/y_d (his (4.15) with footnote 10) has the gauge invariants of the
    book's U_q(b^-) matrix with the same s and y/y_d.  (b) psi -> -psi in BCZ's two-label basis |a, b> (hep-th/0506169
    (3.57)-(3.59)), with the soliton on the right conjugated, maps the b^- matrix to one that, after the diagonal gauge
    q^{m(m+1)/2 - n m}, no longer depends on the second label and has the gauge invariants of the b^+ matrix with the same
    y_d and offset s_+ = -1 - s.  (c) s != 0: with the s-independent scalar factor Phi in the universal-R normalization
    q^{-s/2} * (ch. 37 matrix), unitarity and defect crossing hold for s = 0.3 and s = 1; without q^{-s/2} crossing fails.
 9. The type-I_- matrix intertwines no extension of U_q(b^-) by operators E_j on W with the Hopf coproduct
    e_j (x) 1 + k_j (x) E_j (best least-squares residual 0.88).  (The exact statement, that the generic Weyl module of b^-
    does not extend to U_q(sl2^), is proved in the text from [e_1, f_0] = 0 and [e_1, f_1] = (k_1 - k_1^-1)/kappa.)
10. Classical limit: on a Weyl module U = q^{h/2} and W_+- = q^{-+h/2} kappa a_+- q^{-+h/2} obey U W_+- = q^{+-1} W_+- U,
    W_+ W_- = 1, and L^- = diag(U^-1, U)[1 - i x (E_21 W_+ + E_12 W_-)]: the form of the type-I_- defect matrix L_d of
    ch. 38 (eq-defect-lax-matrix) with U^-1 <-> e^{b(phi - psi)/2sqrt2}, W_+- <-> e^{+-b phi/sqrt2}, -ix <-> sigma/lambda;
    RTT holds for generic q (q = e^{0.05i}, e^{0.6i}), so q -> 1 gives the classical Sklyanin bracket of L_d.
11. thm-t-triangle (ii) without truncation (exact, sympy, generic rational point): End(V(y1) (x) V(y2) (x) W) = scalars for
    the bilateral Weyl modules of both Borels (raising generators invertible; loops raise-then-lower transported to the block
    of total charge 0 leave a one-dimensional commutant) and for F^+ (its one-dimensional lowest block generates the module).
12. prp-sg-t-unique needs q not a root of unity: at q = e^{i pi/3} (q^4 != 1) the shift S^3 commutes with the W^- generators,
    so (S^3 (x) 1) T is a second intertwiner.
Windows: Weyl n = -6..6 (equations on |n| <= 5), Fock n = 0..12 (n <= 9 for check 6); 16 checks, runtime about 85 seconds.
"""
import numpy as np
from boundary_common import report, q_phys, sl2_spin, rcheck, prop_residual, cop
from defect_common import (weyl_module_b, fock_module_b, q_oscillator_module, borel_l_operator, flip_VW,
                           masked_intertwiners, transmission_window, sg_tmatrix, sg_transmission_scalar, borel_kinds)

LAM = 1.37
q = q_phys(LAM)
kap = q - 1/q
qn = lambda n: (q**n - q**-n)/(q - 1/q)
kron = np.kron
S = 0.3
NW = 6


def V(y):
    R = sl2_spin(1, q, y)
    R['h'] = np.array([1.0, -1.0])
    return R


def serre_res(X1, X0, mask):
    s = 0
    for A, B in ((X1, X0), (X0, X1)):
        M = A@A@A@B - qn(3)*A@A@B@A + qn(3)*A@B@A@A - B@A@A@A
        s = max(s, np.abs(M[np.ix_(mask, mask)]).max())
    return s


# ---------------------------------------------------------------- 1. q-oscillator, Weyl and Fock modules
res_osc, res_z, res_b = 0, 0, 0
for kind in ('weyl', 'fock'):
    O = q_oscillator_module(kind, q, 12, S)
    m2 = O['interior'] & ((np.abs(O['n']) <= 10) if kind == 'weyl' else (O['n'] <= 10))
    ap, am, H = O['ap'], O['am'], np.diag(q**O['h'])
    rel = q*ap@am - am@ap/q - np.eye(len(O['n']))/kap
    res_osc = max(res_osc, np.abs(rel[np.ix_(m2, m2)]).max(),
                  np.abs(H@ap@np.linalg.inv(H) - q**2*ap).max(), np.abs(H@am@np.linalg.inv(H) - am/q**2).max())
    Z = np.diag(ap@am) - 1/kap**2
    want = np.zeros(len(O['n'])) if kind == 'weyl' else -q**(2.0*O['n'])/kap**2
    res_z = max(res_z, np.abs((Z - want)[m2]).max())
    for sign in (1, -1):
        W = fock_module_b(sign, q, 1.3*np.exp(0.2j), 12, S, kind)
        x = 'e' if sign == 1 else 'f'
        X1, X0, K1, K0 = W[(x, 1)], W[(x, 0)], W[('k', 1)], W[('k', 0)]
        a = 2 if sign == 1 else -2
        res_b = max(res_b, np.abs(K1@X1@K0 - q**a*X1).max(), np.abs(K1@X0@K0 - q**-a*X0).max(),
                    np.abs(K0@K1 - np.eye(W['dim'])).max(), serre_res(X1, X0, m2))
report(f'eq-q-oscillator q a+ a- - q^-1 a- a+ = 1/(q - q^-1) on the Weyl and Fock modules [{res_osc:.0e}]; Z = 0 on Weyl, '
       f'Z|n> = -q^2n/kappa^2 on Fock [{res_z:.0e}]; U_q(b^+-) -> Osc_q respects k x k^-1 and q-Serre [{res_b:.0e}]',
       res_osc < 1e-12 and res_z < 1e-12 and res_b < 1e-9)

ns = np.arange(-8, 9); m4 = np.abs(ns) <= 4
up = np.diag(np.ones(len(ns) - 1), -1)


def shift_serre(u):
    e1 = up.astype(complex); e0 = (up.T*np.r_[0, u[1:]]).astype(complex)    # e_0 |n> = u_n |n-1>
    return serre_res(e1, e0, m4)


rng = np.random.default_rng(1)
A_, B_, C_ = rng.normal(size=3) + 1j*rng.normal(size=3)
ok_gen = shift_serre(A_ + B_*q**(2.0*ns) + C_*q**(-2.0*ns))
bad = shift_serre(A_ + B_*q**(4.0*ns))
report(f'shift modules: q-Serre holds for u_n = A + B q^2n + C q^-2n [{ok_gen:.0e}] and fails for A + B q^4n [{bad:.2f}] '
       '(Weyl: B = C = 0; q-oscillator Fock: C = 0, A + B = 0; Weston type II: B, C != 0)', ok_gen < 1e-9 and bad > 1e-3)

# ---------------------------------------------------------------- 2. KT universal R (Part VI) to height 3
def roots_e(M):
    e1, e0 = M[('e', 1)], M[('e', 0)]
    ed = e1@e0 - q**2*e0@e1                                     # e'_delta (BGKNR (4.13), q -> 1/q)
    return dict(a=e1, ma=e0, d=ed, ad=(e1@ed - ed@e1)/qn(2), mad=(ed@e0 - e0@ed)/qn(2))


def roots_f(M):
    f1, f0 = M[('f', 1)], M[('f', 0)]
    fd = f0@f1 - q**-2*f1@f0                                    # f'_delta (BGKNR (4.17), q -> 1/q)
    return dict(a=f1, ma=f0, d=fd, ad=(fd@f1 - f1@fd)/qn(2), mad=(f0@fd - fd@f0)/qn(2))


def expq(X, p, nmax=3):
    out = np.eye(X.shape[0], dtype=complex); P = np.eye(X.shape[0], dtype=complex); fac = 1
    for n in range(1, nmax + 1):
        fac *= (1 - p**n)/(1 - p); P = P @ X; out = out + P/fac
    return out


def R_univ(A, B, imag=True, height3=True):
    """R = q^{h(x)h/2} prod_k exp_{q^-2}(kap f_{d-a+kd} (x) e_{d-a+kd}) exp(kap sum k/[2k] f_kd (x) e_kd)
           prod_k^<- exp_{q^-2}(kap f_{a+kd} (x) e_{a+kd}),  truncated: k = 0, 1 real, k = 1 imaginary."""
    fA, eB = roots_f(A), roots_e(B)
    K = np.diag(np.array([q**(a*b/2) for a in A['h'] for b in B['h']]))
    I = np.eye(A['dim']*B['dim'])
    ks = (0, 1) if height3 else (0,)
    left = I
    for k in ks:
        key = 'ma' if k == 0 else 'mad'
        left = left @ expq(kap*kron(fA[key], eB[key]), q**-2)
    mid = I + (kap/qn(2)*kron(fA['d'], eB['d']) if imag else 0)
    right = I
    for k in reversed(ks):
        key = 'a' if k == 0 else 'ad'
        right = right @ expq(kap*kron(fA[key], eB[key]), q**-2)
    return K @ left @ mid @ right


def R_prime(A, B, imag=True, height3=True):
    """R' = R_21^-1 = prod_k exp_{q^2}(-kap e_{a+kd} (x) f_{a+kd}) exp(-kap sum k/[2k] e_kd (x) f_kd)
            prod_k^<- exp_{q^2}(-kap e_{d-a+kd} (x) f_{d-a+kd}) q^{-h(x)h/2}  (BGKNR (4.23) with q -> 1/q)."""
    eA, fB = roots_e(A), roots_f(B)
    K = np.diag(np.array([q**(-a*b/2) for a in A['h'] for b in B['h']]))
    I = np.eye(A['dim']*B['dim'])
    ks = (0, 1) if height3 else (0,)
    left = I
    for k in ks:
        key = 'a' if k == 0 else 'ad'
        left = left @ expq(-kap*kron(eA[key], fB[key]), q**2)
    mid = I - (kap/qn(2)*kron(eA['d'], fB['d']) if imag else 0)
    right = I
    for k in reversed(ks):
        key = 'ma' if k == 0 else 'mad'
        right = right @ expq(-kap*kron(eA[key], fB[key]), q**2)
    return left @ mid @ right @ K


P22 = flip_VW(2, 2)


def vv_res(x, prime=False, **kw):
    if not prime:
        A, B = V(np.exp(0.3j)/x), V(1.0)
        L = R_univ(A, B, **kw)
    else:
        A, B = V(x*np.exp(0.3j)), V(1.0)
        L = R_prime(A, B, **kw)
    return prop_residual(L, P22.T @ rcheck(A, B)[0])[0]


def ratio(f, **kw):
    r1, r2 = f(0.02, **kw), f(0.01, **kw)
    return r1, r1/r2


r_full, k_full = ratio(vv_res)
r_imag, k_imag = ratio(vv_res, imag=False)
r_h3, k_h3 = ratio(vv_res, height3=False)
r_pr, k_pr = ratio(vv_res, prime=True)
report(f'KT universal R (Part VI, height <= 3) on V(y1) (x) V(y2) = P Rcheck up to a scalar: residual {r_full:.1e} at '
       f'x = y2/y1 = 0.02, ratio {k_full:.2f} on halving x (O(x^4)); controls: no imaginary root {k_imag:.2f} (O(x^2)), '
       f'no height-3 roots {k_h3:.2f} (O(x^3)); R_21^-1 {k_pr:.2f}',
       r_full < 1e-6 and abs(k_full - 16) < 1 and abs(k_imag - 4) < 0.3 and abs(k_h3 - 8) < 0.5 and abs(k_pr - 16) < 1)


# ---------------------------------------------------------------- 3./4. Weyl and Fock modules
def l_res(x, sign, kind, s=S, yW=1.3, **kw):
    """KT R (sign +1) or R_21^-1 (sign -1) on V(y) (x) W against borel_l_operator, x = yW/y resp. y/yW."""
    W = fock_module_b(sign, q, yW, 12 if kind == 'fock' else NW, s, kind)
    y = yW/x*np.exp(0.2j) if sign == 1 else yW*x*np.exp(0.2j)
    Vy = V(y)
    L = R_univ(Vy, W, **kw) if sign == 1 else R_prime(Vy, W, **kw)
    Lc = borel_l_operator(sign, y, W)
    msk = np.repeat(((np.abs(W['n']) <= NW - 3) if kind == 'weyl' else (W['n'] <= 8))[None, :], 2, axis=0).ravel()
    r, c = prop_residual(L[np.ix_(msk, msk)], Lc[np.ix_(msk, msk)])
    return r, c


out = {}
for sign in (1, -1):
    for kind in ('weyl', 'fock'):
        (r1, c1), (r2, c2) = l_res(0.02, sign, kind), l_res(0.01, sign, kind)
        out[(sign, kind)] = (r1, r1/r2)
ok = all(v[0] < 1e-6 and abs(v[1] - 16) < 1.5 for v in out.values())
report('KT R evaluated on V (x) (Weyl or Fock module) = scalar x borel_l_operator (eq-l-operator) to O(x^4): '
       + ', '.join(f'{"b+" if k[0] == 1 else "b-"} {k[1]} {v[0]:.0e} (ratio {v[1]:.1f})' for k, v in out.items()), ok)

# b^- Weyl: the universal-R normalization is q^{-s/2} times the ch. 37 matrix
W = fock_module_b(-1, q, 1.3, NW, S, 'weyl')
y = 1.3*0.003*np.exp(0.2j)
L = R_prime(V(y), W)
T37 = sg_tmatrix(np.log(y)/LAM, LAM, q, np.log(1.3)/LAM, NW, s=S, scalar=False)[0]
msk = np.repeat((np.abs(W['n']) <= NW - 3)[None, :], 2, axis=0).ravel()
r, c = prop_residual(L[np.ix_(msk, msk)], (flip_VW(2, W['dim']).T @ T37)[np.ix_(msk, msk)])
report(f'R_21^-1 on V (x) W(b^-) = c x (ch. 37 matrix of eq-sg-transmission-matrix) with c = {c:.4f} -> q^(-s/2) = '
       f'{q**(-S/2):.4f} as y/y_d -> 0 [{r:.0e}]', r < 1e-8 and abs(c - q**(-S/2)) < 1e-4)


def intertwiner_window(Vy, W, m_eq=1, m_blk=3):
    kinds = borel_kinds(W['sign'])
    pairs = [(cop(Vy, W, kd, j), cop(W, Vy, kd, j)) for kd, j in kinds]
    nVW = np.repeat(W['n'][None, :], 2, axis=0).ravel(); nWV = np.repeat(W['n'][:, None], 2, axis=1).ravel()
    top = W['n'].max()
    if W['n'].min() < 0:
        eqr, eqc, br, bc = np.abs(nWV) <= top - m_eq, np.abs(nVW) <= top - m_eq, np.abs(nWV) <= top - m_blk, np.abs(nVW) <= top - m_blk
    else:
        eqr, eqc, br, bc = nWV <= top - m_eq, nVW <= top - m_eq, nWV <= top - m_blk, nVW <= top - m_blk
    null, _ = masked_intertwiners(pairs, eqr, eqc)
    blk = np.ix_(br, bc)
    st = np.array([X[blk].ravel() for X in null])
    sv = np.linalg.svd(st, compute_uv=False)
    _, _, Vh = np.linalg.svd(st)
    T = np.zeros_like(null[0]); T[blk] = Vh[0].reshape(br.sum(), bc.sum())
    return T, int((sv > 1e-5).sum()), blk


rk, dev = [], 0
for sign in (1, -1):
    F = fock_module_b(sign, q, 1.1*np.exp(0.3j), 12, S, 'fock')
    y = 0.8*np.exp(-0.2j)
    T, rank, blk = intertwiner_window(V(y), F)
    Lc = flip_VW(2, F['dim']) @ borel_l_operator(sign, y, F)
    rk.append(rank); dev = max(dev, prop_residual(T[blk], Lc[blk])[0])
report(f'Fock modules of U_q(b^+-): the intertwiner V(y) (x) F -> F (x) V(y) is unique on the window (ranks {rk}) and '
       f'equals the closed form eq-l-operator (BLZ L-operator for b^+) [{dev:.0e}]', rk == [1, 1] and dev < 1e-9)

# ---------------------------------------------------------------- 5. RTT (eq-triangle-rtt) for the L-operators
def rtt_res(sign, kind, y1, y2, yW=1.2):
    W = fock_module_b(sign, q, yW, 10, S, kind); dW = W['dim']
    R = P22 @ rcheck(V(y1), V(y2))[0]           # R = P Rcheck on V(y1) (x) V(y2)
    L1 = borel_l_operator(sign, y1, W); L2 = borel_l_operator(sign, y2, W)
    I2, IW = np.eye(2), np.eye(dW)
    # spaces 1 (x) 2 (x) W; L_13 = flip_23 (L1 (x) 1) flip_23
    F23 = kron(I2, flip_VW(2, dW).T)            # V (x) W (x) V -> V (x) V (x) W
    L13 = F23 @ kron(L1, I2) @ F23.T
    L23 = kron(I2, L2)
    R12 = kron(R, IW)
    lhs, rhs = R12 @ L13 @ L23, L23 @ L13 @ R12
    n = W['n']; ok = (np.abs(n) <= n.max() - 2) if n.min() < 0 else (n <= n.max() - 2)
    msk = np.tile(ok, 4)
    return np.abs((lhs - rhs)[np.ix_(msk, msk)]).max()/np.abs(lhs[np.ix_(msk, msk)]).max()


rr = max(rtt_res(s_, k_, 0.7*np.exp(0.1j), 1.9*np.exp(-0.3j)) for s_ in (1, -1) for k_ in ('weyl', 'fock'))
report(f'eq-triangle-rtt R_12 L_13 L_23 = L_23 L_13 R_12 for the Weyl and Fock L-operators of both Borels [{rr:.0e}]',
       rr < 1e-10)


# ---------------------------------------------------------------- 6. thm-t-triangle hypothesis (ii)
def charge_intertwiners(pairs, hin, hout, rows, cols, blk_rows, blk_cols):
    """Dimension of {X : X A = B X on the masked rows/columns for all (A, B) in pairs, X charge conserving}, restricted
    to the block (blk_rows, blk_cols): unknowns only where hout[i] = hin[j] (the k-equations), dense SVD."""
    idx = [(i, j) for i in range(len(hout)) for j in range(len(hin)) if abs(hout[i] - hin[j]) < 1e-9]
    R, C = np.where(rows)[0], np.where(cols)[0]
    rpos = {r: a for a, r in enumerate(R)}; cpos = {c: b for b, c in enumerate(C)}
    blocks = []
    for A, B in pairs:
        M = np.zeros((len(R)*len(C), len(idx)), complex)
        for u, (i, j) in enumerate(idx):
            if i in rpos:                                  # (E_ij A)[i, c] = A[j, c]
                a = rpos[i]
                M[a*len(C):(a + 1)*len(C), u] += A[j, C]
            if j in cpos:                                  # -(B E_ij)[r, j] = -B[r, i]
                b = cpos[j]
                M[np.arange(len(R))*len(C) + b, u] -= B[R, i]
        blocks.append(M)
    M = np.vstack(blocks)
    G = M.conj().T @ M                                  # Gram matrix: null space = eigenvectors with eigenvalue ~ 0
    ev, Vg = np.linalg.eigh(G)
    sv = np.sqrt(np.clip(ev, 0, None))
    Vh = Vg.conj().T
    null = Vh[sv/sv.max() < 1e-7].conj()
    keep = np.array([blk_rows[i] and blk_cols[j] for i, j in idx])
    st = null[:, keep]
    return int((np.linalg.svd(st, compute_uv=False) > 1e-6).sum()) if len(st) else 0


def tensor(A, B, kinds):
    g = {k: cop(A, B, k[0], k[1]) for k in kinds}
    g['dim'] = A['dim']*B['dim']
    return g


def hom_dim(W, y1, y2):
    kinds = borel_kinds(W['sign'])
    V1, V2 = V(y1), V(y2)
    full = [(k, j) for k in 'efk' for j in (0, 1)]
    V12 = tensor(V1, V2, full); V21 = tensor(V2, V1, full)
    pairs = [(cop(V12, W, kd, j), cop(W, V21, kd, j)) for kd, j in kinds if kd != 'k']
    hVV = np.array([a + b for a in (1, -1) for b in (1, -1)])
    hin = np.add.outer(hVV, W['h']).ravel(); hout = np.add.outer(W['h'], hVV).ravel()
    nin = np.tile(W['n'], 4); nout = np.repeat(W['n'], 4)
    top = W['n'].max()
    f = (lambda n, m: np.abs(n) <= top - m) if W['n'].min() < 0 else (lambda n, m: n <= top - m)
    return charge_intertwiners(pairs, hin, hout, f(nout, 1), f(nin, 1), f(nout, 3), f(nin, 3))


dims = [hom_dim(fock_module_b(-1, q, 1.3, 6, S, 'weyl'), 0.7*np.exp(0.2j), 1.6),
        hom_dim(fock_module_b(1, q, 1.3, 6, S, 'weyl'), 0.7*np.exp(0.2j), 1.6),
        hom_dim(fock_module_b(1, q, 1.3, 9, S, 'fock'), 0.7*np.exp(0.2j), 1.6)]
report(f'thm-t-triangle (ii): Hom_U_q(b)(V(y1) (x) V(y2) (x) W, W (x) V(y2) (x) V(y1)) is one-dimensional on the window '
       f'interior for Weyl b^-, Weyl b^+, Fock b^+ (dimensions {dims})', dims == [1, 1, 1])

# ---------------------------------------------------------------- 7. prp-sg-t-unique
ok_u, dev_u = True, 0
for s_, th in ((0.0, -0.4), (0.3, 0.6), (1.0, 0.2)):
    T37, Vy, _ = sg_tmatrix(th, LAM, q, 0.31, NW, s=s_, scalar=False)
    W37 = fock_module_b(-1, q, np.exp(LAM*0.31), NW, s_, 'weyl')        # = weyl_module_b with mu = i/(q - q^-1)
    T, rank = transmission_window(Vy, W37)
    n = W37['n']
    blk = np.ix_(np.repeat(np.abs(n) <= NW - 1, 2), np.repeat((np.abs(n) <= NW - 1)[None, :], 2, axis=0).ravel())
    ok_u &= rank == 1
    dev_u = max(dev_u, prop_residual(T[blk], T37[blk])[0])
Wbad = weyl_module_b(-1, q, 1.3, NW, s=0.0, mu=1j/kap, kappa=1.7)
_, rank_bad = transmission_window(V(1.1), Wbad)
report(f'prp-sg-t-unique: the U_q(b^-) intertwiner V(y) (x) W -> W (x) V(y) is unique and equals the ch. 37 matrix for '
       f's = 0, 0.3, 1 [{dev_u:.0e}]; with k_0 k_1 = 1.7 there is none (rank {rank_bad})', ok_u and dev_u < 1e-9 and rank_bad == 0)


# ---------------------------------------------------------------- 8a. Weston's L-operator
def weston_check(zeta, r0, K=7):
    """Weston 1006.1555 (2.9) with r1 = r2 = 0 on W (x) V: a|j> = |j-1>, a*|j> = |j+1>, q^D|j> = q^j|j>,
    e_0 = z a*/kappa, e_1 = z a/kappa, t_1 = r0 q^-2D, t_0 = q^2D/r0 (z = zeta_W); V = book V(zeta_V)."""
    js = np.arange(-K, K + 1); d = len(js)
    a = np.diag(np.ones(d - 1), 1).astype(complex); ast = a.T.copy()
    qD = np.diag(q**js.astype(float)).astype(complex); qDi = np.linalg.inv(qD)
    zW, zV = zeta, 1.0
    Wm = {'e0': zW*ast/kap, 'e1': zW*a/kap, 't1': r0*qDi@qDi, 't0': qD@qD/r0}
    Vr = sl2_spin(1, q, zV)
    Vm = {'e1': Vr[('e', 1)], 'e0': Vr[('e', 0)], 't1': Vr[('k', 1)], 't0': Vr[('k', 0)]}
    z = zW/zV
    I = np.eye(d)
    blocks = [[I @ qD, -z*ast @ (r0*qDi)], [-z*a @ qD, I @ (r0*qDi)]]
    E = lambda i, k: np.outer(np.eye(2)[i], np.eye(2)[k])
    L = sum(kron(blocks[i][k], E(i, k)) for i in range(2) for k in range(2))       # on W (x) V
    msk = np.repeat(np.abs(js) <= K - 2, 2)
    res = 0
    for i in (0, 1):
        e, t = 'e%d' % i, 't%d' % i
        Dl = kron(Wm[e], np.eye(2)) + kron(Wm[t], Vm[e])
        Dp = kron(Wm[e], Vm[t]) + kron(I, Vm[e])
        res = max(res, np.abs((L @ Dl - Dp @ L)[np.ix_(msk, msk)]).max())
    # transpose (Weston's footnote 10), as a map on W (x) V; book orientation V (x) W -> W (x) V is T = X (flipped input)
    X = L.T
    Tb = X @ flip_VW(2, d)                       # V (x) W -> W (x) V
    # book labels n = -j:  A_n: v+|n> -> |n> v+ etc.
    idx = {j: k for k, j in enumerate(js)}
    def el(jout, aout, jin, ain):
        return Tb[idx[jout]*2 + aout, ain*d + idx[jin]]
    inv = []
    for n in range(-2, 3):
        j = -n
        A = el(j, 0, j, 0); D = el(j, 1, j, 1)
        A1 = el(j - 1, 0, j - 1, 0)
        C = el(j - 1, 1, j, 0)                   # v+ |n> -> |n+1> v-   (n+1 <-> j-1)
        B1 = el(j, 0, j - 1, 1)                  # v- |n+1> -> |n> v+
        D1 = el(j - 1, 1, j - 1, 1)
        inv.append((A1/A, D/A, C*B1/(A*D1)))
    return res, inv


zeta = 1j*np.exp(LAM*(0.4 - 0.31))                 # Weston (4.15): zeta = i e^{gamma(theta - eta)} = i y/y_d
res_w, inv = weston_check(zeta, q**S)
yr = zeta/1j
want = [(1/q, q**(2*n + S), -(yr)**2) for n in range(-2, 3)]          # book b^-: A_{n+1}/A_n, D_n/A_n, C_n B_{n+1}/(A_n D_{n+1})
dev_w = max(abs(a - b) for u, v in zip(inv, want) for a, b in zip(u, v))
report(f'Weston (2.9) with r = (q^s, 0, 0) intertwines U_q(b^+) on W (x) V, Delta -> Delta\', in Part VI conventions '
       f'[{res_w:.0e}]; its transpose at zeta = i y/y_d has the gauge invariants of the book\'s U_q(b^-) matrix, '
       f'A_n+1/A_n = 1/q, D_n/A_n = q^(2n+s), C_n B_n+1/(A_n D_n+1) = -(y/y_d)^2 [{dev_w:.0e}]', res_w < 1e-12 and dev_w < 1e-10)

# ---------------------------------------------------------------- 8b. psi -> -psi: b^- -> b^+ in BCZ's two labels
K2 = 5
lab = [(n, m) for n in range(-K2, K2 + 1) for m in range(-K2, K2 + 1)]
pos = {l: i for i, l in enumerate(lab)}
d2 = len(lab)
y, yd, s_ = 0.7*np.exp(0.2j), 1.3, S
Am = lambda n: q**-n; Dm = lambda n: q**(n + s_); Cm = lambda n: -1j*(y/yd)*q**-n; Bm = lambda n: -1j*(y/yd)*q**(n + s_)
# b^- two-label matrix: soliton unconverted (a,b)->(a-1,b-1): n same, m-1; converted: n+1, m same; antisoliton mirrored
amp = {}                                           # (in label, in particle) -> list of (out label, out particle, amplitude)
for (n, m) in lab:
    amp[((n, m), 0)] = [((n, m - 1), 0, Am(n)), ((n + 1, m), 1, Cm(n))]
    amp[((n, m), 1)] = [((n - 1, m), 0, Bm(n)), ((n, m + 1), 1, Dm(n))]
# psi -> -psi: (a, b) -> (a, -b), i.e. (n, m) -> (n', m') = (-m, -n); conjugate the outgoing particle
amp2 = {}
for ((n, m), a), outs in amp.items():
    amp2[((-m, -n), a)] = [((-mo, -no), 1 - b, c) for ((no, mo), b, c) in outs]
g = lambda n, m: q**(m*(m + 1)/2 - n*m)            # BCZ-type diagonal gauge (cf. their (3.58))
gauged = {}
for ((n, m), a), outs in amp2.items():
    gauged[((n, m), a)] = [((no, mo), b, g(no, mo)*c/g(n, m)) for ((no, mo), b, c) in outs]
# read off: unconverted soliton (0->0), converted (0->1), etc., normalized by the unconverted soliton at the same label
mdep, invs = 0, []
for n in range(-2, 3):
    vals = []
    for m in range(-2, 3):
        o0 = dict(((b,), c) for (_, b, c) in gauged[((n, m), 0)])
        o1 = dict(((b,), c) for (_, b, c) in gauged[((n, m), 1)])
        vals.append((o0[(0,)], o0[(1,)], o1[(0,)], o1[(1,)]))
    vals = np.array(vals)
    mdep = max(mdep, np.abs(vals - vals[0]).max()/np.abs(vals[0]).max())
    invs.append(vals[0])
invs = np.array(invs)                               # rows n = -2..2: (A'_n, C'_n, B'_n, D'_n) with v+->v+, v+->v-, v-->v+, v-->v-
sp = -1 - s_
Ap = invs[:, 0]; Cp = invs[:, 1]; Bp = invs[:, 2]; Dp = invs[:, 3]
chk = max(max(abs(Ap[i + 1]/Ap[i] - q) for i in range(4)),
          max(abs(Dp[i]/Ap[i] - q**(-2*(i - 2) - sp)) for i in range(5)),
          max(abs(Cp[i]*Bp[i + 1]/(Ap[i]*Dp[i + 1]) + (yd/y)**2) for i in range(4)))
report(f'psi -> -psi (BCZ two-label basis |a,b>, outgoing soliton conjugated): after the gauge q^(m(m+1)/2 - nm) the '
       f'b^- matrix no longer depends on the second label [{mdep:.0e}] and has the b^+ gauge invariants A_n+1/A_n = q, '
       f'D_n/A_n = q^(-2n-s+), C_n B_n+1/(A_n D_n+1) = -(y_d/y)^2 with the same y_d and s+ = -1-s [{chk:.0e}]',
       mdep < 1e-12 and chk < 1e-12)

# ---------------------------------------------------------------- 8c. s != 0: crossing with the universal-R normalization
C_plain = np.array([[0, 1], [1, 0]])
TD = 0.31


def Tfull(th, s_, norm):
    T, V_, W_ = sg_tmatrix(th, LAM, q, TD, NW, s=s_, scalar=False)
    return (q**(-s_/2) if norm else 1.0)*sg_transmission_scalar(th - TD, LAM)*T, W_


def crossing_res(s_, norm, th):
    T, W_ = Tfull(th, s_, norm)
    Tb = np.linalg.inv(Tfull(th - 1j*np.pi, s_, norm)[0])
    dW = W_['dim']; ns_ = W_['n']; res, ref = [], []
    for k in np.where(np.abs(ns_) <= NW - 2)[0]:
        for a in range(2):
            for b in range(2):
                for k2 in range(dW):
                    lhs = T[k2*2 + b, a*dW + k]
                    rhs = sum(C_plain[dd, b]*Tb[c*dW + k2, k*2 + dd]*C_plain[a, c] for c in range(2) for dd in range(2))
                    res.append(lhs - rhs); ref.append(lhs)
    return np.abs(res).max()/np.abs(ref).max()


def unit_res(s_, th):
    T, W_ = Tfull(th, s_, True)
    inn = np.repeat((np.abs(W_['n']) <= NW - 2)[None, :], 2, axis=0).ravel()
    U = T.conj().T @ T
    return np.abs(U[np.ix_(inn, inn)] - np.eye(inn.sum())).max()


cr = max(crossing_res(s_, True, th) for s_ in (0.3, 1.0) for th in (0.4, -0.3))
cr_bad = min(crossing_res(s_, False, th) for s_ in (0.3, 1.0) for th in (0.4, -0.3))
ur = max(unit_res(s_, th) for s_ in (0.3, 1.0) for th in (0.4, -0.3))
report(f's != 0: Phi(theta - theta_d) q^(-s/2) x (ch. 37 matrix), i.e. the universal-R normalization with an '
       f's-independent Phi, is unitary [{ur:.0e}] and satisfies defect crossing for s = 0.3, 1 [{cr:.0e}]; without '
       f'q^(-s/2) crossing fails [{cr_bad:.2f}]', ur < 1e-9 and cr < 1e-8 and cr_bad > 1e-2)

# ---------------------------------------------------------------- 9. no extension of the Borel symmetry by e_j on W
def extension_res(j, th=0.4, s_=0.0):
    """Best E_j on W with T (e_j (x) 1 + k_j (x) E_j) = (E_j (x) 1 + k_j (x) e_j) T on the window interior (least squares)."""
    T, Vy, _ = sg_tmatrix(th, LAM, q, TD, NW, s=s_, scalar=False)
    W_ = fock_module_b(-1, q, np.exp(LAM*TD), NW, s_, 'weyl'); d = W_['dim']
    e, k, kW = Vy[('e', j)], Vy[('k', j)], W_[('k', j)]
    rhs = kron(kW, e) @ T - T @ kron(e, np.eye(d))
    cols = []
    for a in range(d):
        for b in range(d):
            E = np.zeros((d, d), complex); E[a, b] = 1
            cols.append((T @ kron(k, E) - kron(E, np.eye(2)) @ T).ravel())
    Mx = np.array(cols).T
    n = W_['n']
    msk = np.outer(np.repeat(np.abs(n) <= NW - 2, 2), np.repeat((np.abs(n) <= NW - 2)[None, :], 2, axis=0).ravel()).ravel()
    sol, *_ = np.linalg.lstsq(Mx[msk], rhs.ravel()[msk], rcond=None)
    return np.linalg.norm(Mx[msk] @ sol - rhs.ravel()[msk])/np.linalg.norm(rhs.ravel()[msk])


ext = min(extension_res(j) for j in (0, 1))
report(f'the type-I_- matrix intertwines no extension of U_q(b^-) by e_j acting on W with the Hopf coproduct '
       f'e_j (x) 1 + k_j (x) E_j (best least-squares residual {ext:.2f}); the exact no-extension argument for generic q is in sec-perturbed-defects', ext > 0.1)

# ---------------------------------------------------------------- 10. classical limit: Weyl-type structure of L^-
# On a Weyl module, U = q^{h/2}, W_+- = q^{-+h/2} kappa a_+- q^{-+h/2}:  U W_+- = q^{+-1} W_+- U,  W_+ W_- = 1, and
# L^- = diag(U^-1, U) [1 - i x (E_21 (x) W_+ + E_12 (x) W_-)], the form of the type-I_- defect matrix L_d of ch. 38 with
# U^-1 <-> e^{b(phi-psi)/2sqrt2}, W_+- <-> e^{+-b phi/sqrt2}, -ix <-> sigma/lambda.  The L-operator solves RTT for every q.
def weyl_structure(qq, y=0.7*np.exp(0.2j), yd=1.3):
    W_ = fock_module_b(-1, qq, yd, NW, S, 'weyl'); d = W_['dim']; kp = qq - 1/qq
    U = np.diag(qq**(W_['h']/2)); Ui = np.linalg.inv(U)
    Wp = Ui @ (kp*W_['ap']) @ Ui; Wm = U @ (kp*W_['am']) @ U
    inner = np.abs(W_['n']) <= NW - 2; blk = np.ix_(inner, inner)
    rel = max(np.abs((U@Wp - qq*Wp@U)[blk]).max(), np.abs((U@Wm - Wm@U/qq)[blk]).max(),
              np.abs((Wp@Wm - np.eye(d))[blk]).max(), np.abs((Wm@Wp - np.eye(d))[blk]).max())
    x = y/yd
    E = lambda a, b: np.outer(np.eye(2)[a], np.eye(2)[b])
    fact = kron(E(0, 0), Ui) + kron(E(1, 1), U)
    Lf = fact @ (np.eye(2*d) - 1j*x*(kron(E(1, 0), Wp) + kron(E(0, 1), Wm)))
    L = borel_l_operator(-1, y, W_)
    msk = np.repeat(inner[None, :], 2, axis=0).ravel()
    return rel, np.abs((Lf - L)[np.ix_(msk, msk)]).max()


def rtt_generic(qq, y1=0.7*np.exp(0.1j), y2=1.9*np.exp(-0.3j)):
    out = 0
    for sign in (1, -1):
        W_ = fock_module_b(sign, qq, 1.2, 10, S, 'weyl'); dW = W_['dim']
        R = P22 @ rcheck(sl2_spin(1, qq, y1), sl2_spin(1, qq, y2))[0]
        L1, L2 = borel_l_operator(sign, y1, W_), borel_l_operator(sign, y2, W_)
        I2, IW = np.eye(2), np.eye(dW)
        F23 = kron(I2, flip_VW(2, dW).T)
        L13 = F23 @ kron(L1, I2) @ F23.T; L23 = kron(I2, L2); R12 = kron(R, IW)
        lhs, rhs = R12 @ L13 @ L23, L23 @ L13 @ R12
        msk = np.tile(np.abs(W_['n']) <= 8, 4)
        out = max(out, np.abs((lhs - rhs)[np.ix_(msk, msk)]).max()/np.abs(lhs[np.ix_(msk, msk)]).max())
    return out


cl = [weyl_structure(qq) for qq in (q, np.exp(0.05j), np.exp(0.6j))]
rg = max(rtt_generic(qq) for qq in (np.exp(0.05j), np.exp(0.6j)))
report(f'classical limit: on a Weyl module U = q^(h/2), W_+- = q^(-+h/2) kappa a_+- q^(-+h/2) satisfy U W_+- = q^(+-1) W_+- U, '
       f'W_+ W_- = 1 [{max(c[0] for c in cl):.0e}], and L^- = diag(U^-1, U)[1 - ix(E_21 W_+ + E_12 W_-)], the form of L_d '
       f'(ch. 38) [{max(c[1] for c in cl):.0e}]; RTT holds for generic q = e^(0.05i), e^(0.6i) [{rg:.0e}], so q -> 1 gives '
       f'the classical L-operator', max(c[0] for c in cl) < 1e-12 and max(c[1] for c in cl) < 1e-12 and rg < 1e-10)


# ---------------------------------------------------------------- 11. thm-t-triangle (ii) on the infinite modules (exact)
def hom_exact():
    """Exact Hom-dimension argument without truncation (sympy, generic rational point q = 2, y1 = 3/2, y2 = 5/4, y_W = 7/6,
    s = 0).  Hom(V1 (x) V2 (x) W, W (x) V2 (x) V1) = End(V1 (x) V2 (x) W) (compose with an invertible product of local
    intertwiners).  An endomorphism X preserves the total charge l; block l of the Weyl modules has the four states (a1, a2),
    n = l - (a1 + a2)/2.  The raising generator A_l (e_1 for b^+, f_0 for b^-; iterated Part VI coproduct) is invertible
    between blocks, so X_{l+1} = A_l X_l A_l^-1, and X_l commutes with the loop B_{l+1} A_l (B: e_0, resp. f_1).  Transported
    to l = 0, three loops leave a one-dimensional commutant: X is a scalar.  Fock F^+: the one-dimensional lowest block
    (l = -1) generates the blocks l = 0, 1 (ranks 3, 4) and A_l is invertible for l >= 1, so X is again a scalar."""
    import itertools
    import sympy as sp
    qq, ys, yw = sp.Integer(2), (sp.Rational(3, 2), sp.Rational(5, 4)), sp.Rational(7, 6)
    mu, kp = sp.I/(qq - 1/qq), qq - 1/qq
    pairs = list(itertools.product((1, -1), repeat=2))

    def gen(sign, node, kind, l):
        t = 1 if node == 1 else -1                     # k_node = q^{t h}
        step = t*sign                                  # e_1 and f_0 raise l, e_0 and f_1 lower it
        ok = lambda l_, a: kind == 'weyl' or l_ - sum(a)//2 >= 0
        src = [a for a in pairs if ok(l, a)]; dst = [a for a in pairs if ok(l + step, a)]
        M = sp.zeros(len(dst), len(src))
        for c, a in enumerate(src):
            n = l - sum(a)//2
            for i in range(2):                         # particle i: e_j = y E_12 / y E_21, f_j = E_21/y / E_12/y
                if a[i] == -step:
                    b = list(a); b[i] = step
                    amp = ys[i]*qq**(t*sum(a[:i])) if sign == 1 else qq**(-t*(sum(a[i + 1:]) + 2*n))/ys[i]
                    M[dst.index(tuple(b)), c] += amp
            if sign == 1:                              # e_1|n> = mu y|n+1>, e_0|n> = mu y|n-1> (Weyl); i y a_+-/kappa (Fock)
                if kind == 'weyl':
                    amp = mu*yw*qq**(t*sum(a))
                elif step == 1:
                    amp = sp.I*yw/kp*qq**(t*sum(a))
                elif n > 0:
                    amp = sp.I*yw/kp*(1 - qq**(2*n))*qq**(t*sum(a))
                else:
                    continue
            else:                                      # f_1|n> = mu/y |n-1>, f_0|n> = mu/y |n+1>
                amp = mu/yw
            M[dst.index(tuple(a)), c] += amp
        return M

    dims, dets = [], []
    for sign in (1, -1):
        up, down = (1, 0) if sign == 1 else (0, 1)
        A = lambda l: gen(sign, up, 'weyl', l); B = lambda l: gen(sign, down, 'weyl', l)
        T = sp.eye(4); cons = []
        for l in range(3):
            C = T.inv()*B(l + 1)*A(l)*T
            cons.append(sp.kronecker_product(C.T, sp.eye(4)) - sp.kronecker_product(sp.eye(4), C))
            T = A(l)*T
        dims.append(16 - sp.Matrix.vstack(*cons).rank())
        dets += [A(l).det() for l in (-2, 0, 3)]
    A = lambda l: gen(1, 1, 'fock', l); B = lambda l: gen(1, 0, 'fock', l)
    low = A(-1); loop0 = B(1)*A(0)
    reach0 = sp.Matrix.hstack(low, loop0*low, loop0**2*low)
    seed = A(0)*reach0; loop1 = B(2)*A(1)
    ranks = (reach0.rank(), sp.Matrix.hstack(seed, loop1*seed).rank())
    dets += [A(l).det() for l in (1, 2, 3)]
    return dims, ranks, all(d != 0 for d in dets)


hdims, hranks, hdet = hom_exact()
report(f'thm-t-triangle (ii) without truncation, exactly at a generic rational point: End(V(y1) (x) V(y2) (x) W) is the scalars '
       f'for the bilateral Weyl modules of b^+ and b^- (commutant of three transported loops: dimensions {hdims}) and for the '
       f'Fock module F^+ (lowest block generates: ranks {hranks} = (3, 4)); raising generators invertible: {hdet}',
       hdims == [1, 1] and hranks == (3, 4) and hdet)

# ---------------------------------------------------------------- 12. prp-sg-t-unique fails at a root of unity
q6 = np.exp(1j*np.pi/3)                                # q^4 != 1, q^6 = 1
dev6, dev_gen = 0, 0
for qq_, store in ((q6, 'root'), (q, 'generic')):
    W6 = fock_module_b(-1, qq_, 1.3, NW, S, 'weyl')
    n = W6['n']; d6 = W6['dim']
    S3 = np.zeros((d6, d6)); S3[np.arange(3, d6), np.arange(d6 - 3)] = 1        # |n> -> |n+3>
    cols = np.abs(n) <= NW - 4
    res = max(np.abs((S3 @ W6[g] - W6[g] @ S3)[:, cols]).max() for g in (('f', 0), ('f', 1), ('k', 0), ('k', 1)))
    if store == 'root':
        dev6 = res
    else:
        dev_gen = res
report(f'prp-sg-t-unique needs q not a root of unity: at q = e^(i pi/3) (q^4 != 1, q^6 = 1) the shift S^3 commutes with f_0, f_1, '
       f'k_0, k_1 on W^- [{dev6:.0e}], so (S^3 (x) 1) T is a second, independent intertwiner; control: at q_phys [S^3, k_1] != 0 '
       f'[{dev_gen:.2f}]', dev6 < 1e-12 and dev_gen > 1e-2)
