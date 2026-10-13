"""Chapter 41 (sec-soliton-transmission): soliton transmission matrices of the a_r^(1) defects and of the a_2^(2)
(Tzitzeica) defect, from the Borel symmetry of ch. 40.

Conventions: Part VI coproduct Delta(e) = e (x) 1 + k (x) e, Delta(f) = f (x) k^-1 + 1 (x) f (boundary_common.cop);
q = q_phys = -e^{-i pi gamma} (eq-q-physical), gamma = omega = 4 pi/beta^2 - 1 = 1.37 unless stated; principal gradation
y = e^{gamma theta}, y_d = e^{gamma theta_d}; V = vector solitons (species 1, weights e_i), Vbar = species r (weights -e_i)
(boundary_common.an_vector); T: V (x) W -> W (x) V, tensor order = spatial order.  The defect space W is the Weyl module
of eq-dt-weyl (defect_common.an_weyl_module): labels N + s, N in the root lattice of a_r, offset s (class p: s = lambda_p),
   k_j |N> = q^{alpha_j.(N+s)} |N>,   f_j |N> = (mu/y_d) q^{(N+s)_{j+1}} |N - alpha_j>   (U_q(b^-), type I_-),
                                       e_j |N> = mu y_d q^{(N+s)_{j+1}} |N + alpha_j>     (U_q(b^+), type I_+),
mu = -e^{-i pi gamma/2}/(q - q^-1), so that X = -(q - q^-1) mu y/y_d = e^{gamma(theta - theta_d - i pi/2)} (CZ's xhat).
Windows: root-lattice balls |N|^2 <= R2 (defect_common.an_root_window); equations on 'interior' labels, comparisons on
'deep' labels.  Window sizes (labels/interior/deep): a_2 R2=18: 37/19/7; a_2 R2=14: 31/13/1; a_3 R2=18: 177/55/13;
a_3 R2=14: 135/43/1; a_4 R2=8: 171/21/1; a_4 R2=12: 331/51/1.

Checks (book ids; 28 PASS lines):
 (a) eq-dt-weyl: Cartan relations, q-Serre relations and level zero on the Weyl modules of U_q(b^-) and U_q(b^+) for
     a_2, a_3 (orientations +1 and -1, generic s); without the q-powers (orient = 0) q-Serre fails (control).
 (b) prp-an-transmission (1)-(2), eq-dt-tv, eq-dt-tvbar, eq-dt-tplus: the U_q(b^-) intertwiners V (x) W -> W (x) V and
     Vbar (x) W -> W (x) Vbar have a one-dimensional null space on the whole window for a_2 (R2=18) and a_3 (R2=18),
     generic s, and equal the closed forms on deep labels; the closed forms intertwine for a_2, a_3, a_4 and s = 0,
     lambda_1, generic; U_q(b^+) (a_2, s = 0 and generic): unique, Vbar linear = eq-dt-tplus, V full; eq-dt-tplus
     intertwines for a_2, a_3, a_4; orientation -1 (cz2007's D^T): V full, Vbar linear with pi^-1 (charge-conjugate
     defect); with k_0 -> 1.3 k_0 there is no intertwiner (level zero).
 (c) prp-an-transmission (3), thm-t-triangle: the closed forms satisfy eq-triangle-braid with boundary_common.rcheck
     (a_2: V-V, Vbar-Vbar, V-Vbar; a_3: V-V, Vbar-Vbar), control: the matrices built at -q; hypothesis (ii) of
     thm-t-triangle for a_2 (generic s, R2=14): Hom(V(y1) (x) V(y2) (x) W, W (x) V(y2) (x) V(y1)) has nullity 1.
 (d) prp-transmission-crossing, eq-dt-crossing, eq-dt-bootstrap, eq-dt-g: crossing (I), (II) with C = 1 and
     c_I = 1 + (-q)^r X^h, c_II = 1 - q X^h exactly (fitted a = 1), r = 2, 3, 4, s = 0, lambda_1, generic (s-independent);
     the sign of q: dual V(y) = Vbar(y e^{i pi gamma}) with C = 1 at q_phys; at -q_phys no conjugation for h odd and
     C' = diag((-1)^i) for h even; crossing at -q fails with C = 1 (all r) and with C' (h odd), holds with C' (h even);
     the scalar factor g: cz2009a (6.8), (I), (II) to 1e-9 with a real positive c_r (gamma = 1.37, 0.6); cz2009a (6.10) as
     printed: P_K(z)/P_K(0) = prod_{k<K} F_k(z)/F_k(0) Gamma(a_K - z)/Gamma(a_K) exactly, a drift like (h gamma K)^{-z},
     eq-dt-g converges; the bootstrap Vbar = Lambda^r V gives M_Vbar with ratio 1 (a_2, a_3).
 (e) prp-a2-unstable, eq-dt-pole, eq-dt-window: simple poles of g at theta_1 = Delta - i(pi r/h + pi/(h gamma)) and of gbar
     at Delta - i(pi/h + pi/(h gamma)); det of each charge block of M_V = q^r (1 + (-1/q)^r X^h), rank h - 1 at theta_1
     (a_2, a_3, a_4); eq-dt-window, open at both ends, from the signs of E at the poles (6.13) for the pole in the strip
     -pi/2 < Im theta_a < 0 (h = 3..6), Im E = 0 at the end 1/gamma = a - h/2 (included by cz2009a (6.15)), the signs recur
     in remote strips, the a_2 windows in book labels, and the failure of cz2009a's text; the kernel of the residue:
     T_23 Rc_12 on the kernel of M_V(theta_1) is the matrix of the module with offset s + lambda_1 with X^h -> q^2 X^h
     (a_2, s = 0 and lambda_1); the shift theta_d -> theta_d + 2 pi i/h + 2 pi i k/(h gamma) is fixed only up to the integer k,
     and k = 0, 1 give opposite signs of Im E at the a_5 real-absorption candidate (species 1, class 1).
 (f) prp-transmission-semiclassical: (2 pi/(h(gamma + 1))) d arg g/d theta -> ln|z_1(theta - theta_d)| of eq-dc-ar-delay
     (a_2, a_3; g and gbar), deviations at gamma = 10, 40, 160 decreasing like 1/gamma; arg z_1(theta_d) = pi/h; the
     conversion threshold of exr-dt-semiclassical for a soliton a distance eps inside the edge of its sector.
 (g) prp-a22-transmission: CZ 2011's U_q(a_2^(2)) (d = (4, 1)) in Part VI form; 3-dim representation; Borel module with
     F(N) = f1 q^{-2N-1} + f2 q^{2N+1} (quintic q-Serre included), F = 1 and (-1)^N fail; unique U_q(b^+) intertwiners for
     four (f1, f2); triangle relation with rcheck (control: perturbed R); conjugation dual V(x) = V(x x_c), x_c = -q^-6,
     antidiagonal C; crossing at the algebraic crossing point with C, at -q with C(-q); fails at -q with C(q).
 (h) prp-an-transmission (4), the selection rule: M_V has only v_i -> v_i, v_{i-1} (pi); M_Vbar has all h entries,
     |entry j| = |X|^j; the classically forbidden ones (1 <= j <= h - 2) relative to the largest vanish as gamma -> oo at
     fixed theta != theta_d and equal 1 at theta = theta_d; for r >= 2 the matrix part is not unitary for real theta.
 (i) prp-an-transmission (2)-(3), eq-dt-reversed, uniqueness for every r: exact (Fractions, q = 2) loop conditions on one
     sector leave a one-dimensional space for a_2, a_3, both types, V and Vbar, s = 0 and generic; hypothesis (ii) of
     thm-t-triangle for a_2 (81 unknowns, nullity 1); the closed forms satisfy the sector equations exactly (h = 3, 4, 5);
     the structure of Omega_j used in the proof for every r (h = 3..6); the sector determinants, independent of sector and s;
     eq-dt-reversed (exact) and its coefficient qX in the basis q^{-|n|^2/2}|n> (cz2007 (5.10)); the type-I_+ intertwiner for
     V as the inverse of eq-dt-reversed with offset -s, labels -N, X -> -X_+ (Chevalley involution; exact, h = 3..6); at
     q^5 = 1 the shift by 5 alpha_1 commutes with W, and charge-conserving intertwiners are still unique (q^5 = 1, gamma = 1.37);
     prp-an-transmission (4) for every r: the terms of order 0 and 1 of the universal R-matrix, R = q^t(1 + (q - q^-1)
     sum f_j (x) e_j + ...) for type I_+ and R_21^-1 for type I_-, on V, Vbar (x) W equal the matrix parts (h = 3, 4, 5);
     eq-triangle-braid for a_4 on complete charge sectors, both types, VV, VbarVbar, VVbar, VbarV (control: -q).
 (j) cz2011 (3.12): with z = i(theta - nu)/xi and the Gamma functions grouped at equal k, rho-hat solves (3.8) and (3.10)
     when the spectral variable of (3.2) is x = e^{2 pi theta/xi}, so that tau^2 x = e^{-2 pi i z}; with the printed
     x = q^{2 pi theta/xi} (control) (3.10) fails; its outside Gamma functions have poles at nu - i(xi/2 + pi/3),
     nu - i(xi/2 + 2 pi/3), not at the printed xi; the decaying-energy windows pi/3 < xi < 4 pi/3 and -pi/3 < xi < 2 pi/3;
     as printed the products converge only after normalization at z = 0, and then drift like K^{-2z}.
Runtime: about 200 seconds.
"""
import time
import numpy as np
import scipy.sparse as sps
from scipy.special import loggamma
from boundary_common import (report, q_phys, an_vector, rcheck, prop_residual, cop, dual, intertwiners, nodes_of,
                             _qbinom, qnum)
from defect_common import (an_root_window, an_weyl_module, an_vweights, an_solve_intertwiner, an_T_vector,
                           an_T_plus_vbar, an_block_inverse, defect_crossing_residual, an_tensor_rep,
                           an_triangle_residual, an_scalar_g, an_scalar_gbar, an_scalar_constant, richardson_sum)

T0 = time.time()
LAM = 1.37
q = q_phys(LAM)
KAP = q - 1/q
MU = -np.exp(-1j*np.pi*LAM/2)/KAP
I3 = np.eye(3)


def lam1(h):
    return np.eye(h)[0] - 1/h


def s_generic(h):
    s = np.array([0.3, -0.1, 0.25, -0.45, 0.17][:h])
    return s - s.mean()


def Xof(theta, theta_d, gam=LAM):
    return np.exp(gam*(theta - theta_d - 1j*np.pi/2))


# ---------------------------------------------------------------- (a) Weyl modules
def cartan_affine(h):
    A = 2*np.eye(h)
    for i in range(h):
        A[i, (i + 1) % h] = A[(i + 1) % h, i] = -1
    return A


def borel_relations(W, R2small):
    """Relative residuals of k_i X_j k_i^-1 = q^{+-a_ij} X_j and of the q-Serre relations (each relative to the size of
    its terms), applied to the labels with |N|^2 <= R2small (for which all products stay inside the window); and the
    level-zero residual |k_0 ... k_r - 1|."""
    h = W['h']; qq = W['q']; A = cartan_affine(h)
    X = 'f' if W['sign'] == -1 else 'e'
    cols = np.where((W['labels']**2).sum(1) <= R2small + 1e-9)[0]
    P = np.zeros((W['dim'], len(cols)), complex); P[cols, np.arange(len(cols))] = 1
    G = [sps.csr_matrix(W[(X, j)]) for j in range(h)]
    kd = [np.diag(W[('k', j)]) for j in range(h)]
    res = 0
    for i in range(h):
        for j in range(h):
            XP = G[j] @ P
            lhs = kd[i][:, None]*XP/kd[i][cols][None, :]
            res = max(res, np.abs(lhs - qq**(A[i, j]*(1 if X == 'e' else -1))*XP).max()/np.abs(XP).max())
            if i != j:
                m = int(1 - A[i, j])
                terms = []
                for t in range(m + 1):
                    v = P
                    for _ in range(t):
                        v = G[i] @ v
                    v = G[j] @ v
                    for _ in range(m - t):
                        v = G[i] @ v
                    terms.append((-1)**t*_qbinom(m, t, qq)*v)
                res = max(res, np.abs(sum(terms)).max()/max(np.abs(tt).max() for tt in terms))
    return res, np.abs(np.prod(kd, axis=0) - 1).max()


resA = 0; lev = 0
for r, R2, R2s in ((2, 50, 6), (3, 40, 2)):
    h = r + 1
    for sign in (-1, 1):
        for orient in (1, -1):
            W = an_weyl_module(r, sign, q, 1.3, R2, s=s_generic(h), mu=MU, orient=orient)
            a, b = borel_relations(W, R2s)
            resA = max(resA, a); lev = max(lev, b)
W0 = an_weyl_module(2, -1, q, 1.3, 50, s=s_generic(3), mu=MU, orient=0)
bad, _ = borel_relations(W0, 6)
report(f'eq-dt-weyl: Weyl modules of U_q(b^-) and U_q(b^+) for a_2, a_3 (orientations +-1, generic s): Cartan and '
       f'q-Serre relations [{resA:.0e}], level zero k_0...k_r = 1 [{lev:.0e}]; without the q-powers q-Serre fails '
       f'[{bad:.2f}] (all relative)', resA < 1e-12 and lev < 1e-12 and bad > 0.1)


# ---------------------------------------------------------------- (b) intertwiners
def nullity(rel, tol=1e-6):
    return int((np.asarray(rel) < tol).sum())


def deep_compare(T, Tc, W, dV):
    rows = np.repeat(W['deep'], dV); cols = np.tile(W['deep'], dV)
    return prop_residual(T[np.ix_(rows, cols)], Tc[np.ix_(rows, cols)])


def intertwining_residual(T, V, W, kind=None):
    """max |T Delta_{V,W}(b) - Delta_{W,V}(b) T| / max|T| on interior labels, b = f_j (b^-) or e_j (b^+), j = 0..r."""
    kind = kind or ('f' if W['sign'] == -1 else 'e')
    dV = V['dim']
    rows = np.repeat(W['interior'], dV); cols = np.tile(W['interior'], dV)
    Ts = sps.csr_matrix(T); res = 0
    for j in range(W['h']):
        A = sps.csr_matrix(cop(V, W, kind, j)); B = sps.csr_matrix(cop(W, V, kind, j))
        D = (Ts @ A - B @ Ts).toarray()[np.ix_(rows, cols)]
        res = max(res, np.abs(D).max())
    return res/np.abs(T).max()


y, yd = 0.7, 1.3
out = []
ok = True
for r, R2 in ((2, 18), (3, 18)):
    h = r + 1
    W = an_weyl_module(r, -1, q, yd, R2, s=s_generic(h), mu=MU)
    X = -KAP*MU*y/yd
    for conj in (False, True):
        V = an_vector(r, q, y, conj=conj); wv = an_vweights(h, conj)
        rank, T, rel = an_solve_intertwiner(V, V, W, wv, wv)
        res, _ = deep_compare(T, an_T_vector(W, X, conj), W, h)
        out.append((r, conj, nullity(rel), rank, res, rel[1]))
        ok &= nullity(rel) == 1 and rank == 1 and res < 1e-12 and rel[1] > 1e-2
report('prp-an-transmission (1): U_q(b^-) intertwiners V (x) W -> W (x) V and Vbar (x) W -> W (x) Vbar are unique on the '
       'window (nullity, deep rank, residual of eq-dt-tv/eq-dt-tvbar, gap) '
       + '; '.join(f'a_{r} {"Vbar" if c else "V"}: {n}, {k}, {e:.0e}, {g:.2f}' for r, c, n, k, e, g in out), ok)

res_cf = 0
for r, R2 in ((2, 18), (3, 14), (4, 8)):
    h = r + 1
    for s in (np.zeros(h), lam1(h), s_generic(h)):
        W = an_weyl_module(r, -1, q, yd, R2, s=s, mu=MU)
        for conj in (False, True):
            res_cf = max(res_cf, intertwining_residual(an_T_vector(W, -KAP*MU*y/yd, conj), an_vector(r, q, y, conj=conj), W))
report(f'eq-dt-tv, eq-dt-tvbar intertwine U_q(b^-) for a_2, a_3, a_4 and s = 0, lambda_1, generic [{res_cf:.0e}]',
       res_cf < 1e-12)

# U_q(b^+): uniqueness, Vbar linear (closed form), V full
okp = True; outp = []
for s in (np.zeros(3), s_generic(3)):
    W = an_weyl_module(2, 1, q, yd, 18, s=s, mu=MU)
    Xp = KAP*MU*yd/y
    for conj in (False, True):
        V = an_vector(2, q, y, conj=conj); wv = an_vweights(3, conj)
        rank, T, rel = an_solve_intertwiner(V, V, W, wv, wv)
        okp &= nullity(rel) == 1 and rank == 1
        # number of distinct output charges per input on deep labels
        k0 = W['look'][(0, 0, 0)]
        nz = sum(abs(T[kk*3 + b, 0*W['dim'] + k0]) > 1e-8*np.abs(T).max() for kk in range(W['dim']) for b in range(3))
        if conj:
            res, _ = deep_compare(T, an_T_plus_vbar(W, Xp), W, 3)
            okp &= res < 1e-12 and nz == 2
            outp.append(f'Vbar: {nz} entries, closed form {res:.0e}')
        else:
            okp &= nz == 3
            outp.append(f'V: {nz} entries')
    okp &= intertwining_residual(an_T_plus_vbar(W, Xp), an_vector(2, q, y, conj=True), W) < 1e-12
for r, R2 in ((3, 14), (4, 8)):
    W = an_weyl_module(r, 1, q, yd, R2, s=s_generic(r + 1), mu=MU)
    okp &= intertwining_residual(an_T_plus_vbar(W, KAP*MU*yd/y), an_vector(r, q, y, conj=True), W) < 1e-12
report('prp-an-transmission (2), eq-dt-tplus: U_q(b^+) intertwiners are unique (a_2, s = 0 and generic); '
       + '; '.join(outp[:2]) + '; eq-dt-tplus intertwines for a_2, a_3, a_4, generic s', okp)

# orientation -1 (cz2007's convention D^T): V full, Vbar linear with the inverse rotation; level zero
outo = []; oko = True
for conj in (False, True):
    W = an_weyl_module(2, -1, q, yd, 18, s=s_generic(3), mu=MU, orient=-1)
    V = an_vector(2, q, y, conj=conj); wv = an_vweights(3, conj)
    rank, T, rel = an_solve_intertwiner(V, V, W, wv, wv)
    look = W['look']; k0 = look[(0, 0, 0)]; sg = -1 if conj else 1
    moves = [b for b in range(3) if abs(T[look[tuple(map(int, sg*(I3[0] - I3[b])))]*3 + b, k0]) > 1e-8*np.abs(T).max()]
    outo.append(f'{"Vbar" if conj else "V"}: 1 -> {[b + 1 for b in moves]}')
    oko &= nullity(rel) == 1 and moves == ([0, 1] if conj else [0, 1, 2])
Wk = an_weyl_module(2, -1, q, yd, 18, s=s_generic(3), mu=MU)
Wk[('k', 0)] = 1.3*Wk[('k', 0)]
V = an_vector(2, q, y); wv = an_vweights(3, False)
_, _, relk = an_solve_intertwiner(V, V, Wk, wv, wv)
report(f'orientation -1 (cz2007\'s D^T): unique, {outo[0]}, {outo[1]} (vb_1 -> vb_2 is pi^-1: the charge-conjugate '
       f'defect, species a <-> h - a); level zero: with k_0 -> 1.3 k_0 there is no intertwiner (smallest relative '
       f'singular value {relk[0]:.2f})', oko and relk[0] > 1e-3)


# ---------------------------------------------------------------- (c) triangle relation, hypothesis (ii)
def Rc(r, ya, yb, ca, cb, qq=q):
    return rcheck(an_vector(r, qq, ya, conj=ca), an_vector(r, qq, yb, conj=cb))[0]


y1, y2 = 0.7, 1.9*np.exp(0.2j)
tri = []; okt = True
for r, R2 in ((2, 18), (3, 18)):
    h = r + 1
    W = an_weyl_module(r, -1, q, yd, R2, s=s_generic(h), mu=MU)
    for c1, c2 in ((False, False), (True, True), (False, True))[:3 if r == 2 else 2]:
        T1 = an_T_vector(W, -KAP*MU*y1/yd, c1); T2 = an_T_vector(W, -KAP*MU*y2/yd, c2)
        res, _ = an_triangle_residual(Rc(r, y1, y2, c1, c2), T1, T2, W, h, h)
        tri.append(res); okt &= res < 1e-12
# control: the closed forms built at -q_phys (module and particles) against the book's R-check at q_phys
qm = -q
W = an_weyl_module(2, -1, qm, yd, 18, s=s_generic(3), mu=-np.exp(-1j*np.pi*LAM/2)/(qm - 1/qm))
ctrl, _ = an_triangle_residual(Rc(2, y1, y2, False, False), an_T_vector(W, -KAP*MU*y1/yd), an_T_vector(W, -KAP*MU*y2/yd),
                               W, 3, 3)
report('thm-t-triangle for a_r: eq-dt-tv, eq-dt-tvbar satisfy eq-triangle-braid with the book R-check '
       f'(a_2: V-V, Vbar-Vbar, V-Vbar; a_3: V-V, Vbar-Vbar) [{max(tri):.0e}]; control: matrices at -q [{ctrl:.2f}]',
       okt and ctrl > 1e-2)

t1 = time.time()
W = an_weyl_module(2, -1, q, 1.3, 14, s=np.array([0.3, -0.1, -0.2]), mu=MU)
V1, V2 = an_vector(2, q, 0.7), an_vector(2, q, 1.9)
win = [I3[a] + I3[b] for a in range(3) for b in range(3)]
rank, _, rel = an_solve_intertwiner(an_tensor_rep(V1, V2), an_tensor_rep(V2, V1), W, win, win)
report(f'thm-t-triangle (ii) for a_2: Hom_U_q(b^-)(V(y1) (x) V(y2) (x) W, W (x) V(y2) (x) V(y1)) has nullity '
       f'{nullity(rel)} on the window (31 labels, gap {rel[1]:.2f}; {time.time() - t1:.0f} s)', nullity(rel) == 1 and rank == 1)
print(f'# (a)-(c) done after {time.time() - T0:.0f} s')


# ---------------------------------------------------------------- (d) crossing, scalar factor, bootstrap
def module(r, R2, s, qq=q, theta_d=0.25):
    mu = -np.exp(-1j*np.pi*LAM/2)/(qq - 1/qq)
    return an_weyl_module(r, -1, qq, np.exp(LAM*theta_d), R2, s=s, mu=mu)


def crossing(W, X, which, C):
    """(I): T = M_V(X), conjugate M_Vbar at X(theta - i pi) = -q_phys X;  (II): roles exchanged."""
    h = W['h']; Xs = -q_phys(LAM)*X
    first, second = (False, True) if which == 'I' else (True, False)
    T = an_T_vector(W, X, first)
    Tinv, done = an_block_inverse(an_T_vector(W, Xs, second), W, h, an_vweights(h, second))
    return defect_crossing_residual(T, Tinv, done, W['interior'], W['dim'], h, C)


WIN = {2: 18, 3: 14, 4: 8}
thetas = (0.4, -0.3, 0.1 + 0.2j, 0.9)
fits = {}; resmax = 0
for r in (2, 3, 4):
    h = r + 1
    for which in ('I', 'II'):
        for sname, s in (('0', np.zeros(h)), ('lambda_1', lam1(h)), ('generic', s_generic(h))):
            W = module(r, WIN[r], s)
            cs, xs = [], []
            for th in thetas:
                X = Xof(th, 0.25)
                res, c = crossing(W, X, which, np.eye(h))
                resmax = max(resmax, res); cs.append(c); xs.append(X**h)
            A = np.array([[1, x] for x in xs]); sol = np.linalg.lstsq(A, np.array(cs), rcond=None)[0]
            fits[(r, which, sname)] = (sol[0], sol[1]/sol[0], np.abs(A @ sol - np.array(cs)).max())
okI = all(abs(fits[(r, 'I', s)][1] - (-q)**r) < 1e-10 and abs(fits[(r, 'II', s)][1] + q) < 1e-10
          and fits[(r, 'I', s)][2] < 1e-10 and fits[(r, 'II', s)][2] < 1e-10
          for r in (2, 3, 4) for s in ('0', 'lambda_1', 'generic'))
a_I = [fits[(r, 'I', s)][0] for r in (2, 3, 4) for s in ('0', 'lambda_1', 'generic')]
a_II = [fits[(r, 'II', s)][0] for r in (2, 3, 4) for s in ('0', 'lambda_1', 'generic')]
okI &= max(abs(np.array(a_I) - 1)) < 1e-10 and max(abs(np.array(a_II) - 1)) < 1e-10
report(f'eq-dt-crossing: crossing (I), (II) holds with C = 1 [{resmax:.0e}] with c_I = a(1 + b X^h), b = (-q)^r, and '
       f'c_II = a(1 + b X^h), b = -q, a = 1 [{max(abs(np.array(a_I + a_II) - 1)):.0e}], for r = 2, 3, 4 and s = 0, '
       f'lambda_1, generic: the crossing scalars do not depend on s', okI and resmax < 1e-10)


def conjugations(r, qq, yy=0.8):
    """Intertwiners dual(V(y)) -> Vbar(y e^{i pi gamma}) (all generators)."""
    A = dual(an_vector(r, qq, yy)); B = an_vector(r, qq, yy*np.exp(1j*np.pi*LAM), conj=True)
    null, _ = intertwiners([(A[(kd, j)], B[(kd, j)]) for j in nodes_of(B) for kd in 'efk'])
    return null


okc = True; outc = []
for r in (2, 3, 4):
    h = r + 1; Cp = np.diag([(-1.0)**i for i in range(h)])
    nq = conjugations(r, q); nm = conjugations(r, -q)
    okc &= len(nq) == 1 and np.allclose(nq[0]/nq[0][0, 0], np.eye(h), atol=1e-10)
    if h % 2:
        okc &= len(nm) == 0
    else:
        okc &= len(nm) == 1 and np.allclose(nm[0]/nm[0][0, 0], Cp, atol=1e-10)
    W = module(r, WIN[r], np.zeros(h), qq=-q)
    rs = [max(crossing(W, Xof(th, 0.25), 'I', C)[0] for th in (0.4, -0.3)) for C in (np.eye(h), Cp)]
    outc.append(f'a_{r}: {rs[0]:.2f}, {rs[1]:.0e}' if h % 2 == 0 else f'a_{r}: {rs[0]:.2f}, {rs[1]:.2f}')
    okc &= rs[0] > 0.1 and (rs[1] < 1e-10 if h % 2 == 0 else rs[1] > 0.1)
report('the sign of q: dual V(y) = Vbar(y e^{i pi gamma}) with C = 1 at q_phys (r = 2, 3, 4); at -q_phys there is no '
       'such conjugation for h odd and C\' = diag((-1)^i) for h even; crossing (I) at -q_phys with C = 1 and with C\' '
       '(residuals) ' + '; '.join(outc), okc)

# scalar factor
okg = True; outg = []
for r in (2, 3, 4):
    h = r + 1
    for gam in (LAM, 0.6):
        qq = q_phys(gam); Q = 1/qq
        c = an_scalar_constant(r, gam, qq)
        g = lambda D: an_scalar_g(D, r, gam, c=c)
        gb = lambda D: an_scalar_gbar(D, r, gam, c=c)
        Xh = lambda D: np.exp(gam*(D - 1j*np.pi/2))**h
        d68 = max(abs(g(D + 2j*np.pi)/g(D)/((1 + (-Q)**r*Xh(D))/(1 + (-Q)**(r + 2)*Xh(D))) - 1) for D in (0.3, -0.4, 0.2 + 0.1j))
        dII = max(abs(gb(D)*g(D - 1j*np.pi)*(1 - qq*Xh(D)) - 1) for D in (0.5, -0.6, 0.1 + 0.2j))
        dI = max(abs(g(D)*gb(D - 1j*np.pi)*(1 + (-qq)**r*Xh(D)) - 1) for D in (0.5, -0.6, 0.1 + 0.2j))
        okg &= max(d68, dI, dII) < 1e-8 and abs(np.imag(c)) < 1e-12 and np.real(c) > 0
        outg.append(max(d68, dI, dII))
        if gam == LAM:
            outg.append(f'c_{r} = {np.real(c):.4f}')
report(f'eq-dt-g: g(theta + 2 pi i)/g = cz2009a (6.8), and (I) g(theta) gbar(theta - i pi) c_I = 1, (II) gbar(theta) '
       f'g(theta - i pi) c_II = 1 with gbar = prod g, for r = 2, 3, 4, gamma = 1.37, 0.6 [{max(x for x in outg if not isinstance(x, str)):.0e}]; '
       + ', '.join(x for x in outg if isinstance(x, str)) + ' (real, positive)', okg)


def lg_cz(z, r, gam, K):
    """log of cz2009a (6.10) as printed, truncated at k = K."""
    h = r + 1; k = np.arange(1, K + 1)
    return loggamma(0.5 + r*gam/2 - z) + np.sum(loggamma(0.5 + (h*k + r/2)*gam - z) + loggamma(0.5 + (h*k - 1 - r/2)*gam + z)
                                                - loggamma(0.5 + (h*k - r/2)*gam - z) - loggamma(0.5 + (h*k - r/2)*gam + z))


def lg_book(z, r, gam, K):
    h = r + 1; k = np.arange(K)
    a = 0.5 + r*gam/2 + k*h*gam; b = a + gam
    return np.sum(loggamma(a + z) + loggamma(a - z) - loggamma(b + z) - loggamma(b - z) - 2*loggamma(a) + 2*loggamma(b))


# exact finite identity: P_K(z)/P_K(0) = prod_{k<K} F_k(z)/F_k(0) * Gamma(a_K - z)/Gamma(a_K), a_K = 1/2 + (hK + r/2) gamma
idk = 0
for r in (1, 2, 3, 4):
    for gam in (LAM, 0.6):
        for z in (0.43 + 0.72j, -0.31 + 1.1j, 0.2 - 0.4j):
            for K in (3, 17, 60):
                aK = 0.5 + ((r + 1)*K + r/2)*gam
                d = (lg_cz(z, r, gam, K) - lg_cz(0, r, gam, K)) - (lg_book(z, r, gam, K) + loggamma(aK - z) - loggamma(aK))
                idk = max(idk, abs(complex(d.real, (d.imag + np.pi) % (2*np.pi) - np.pi)))       # log defined mod 2 pi i
dcz = []; dbk = []
for r in (2, 3):
    for D in (0.3, -0.7):
        z = 1j*(r + 1)*LAM*(D - 1j*np.pi/2)/(2*np.pi)
        K = 2000
        d = (lg_cz(z, r, LAM, 2*K) - lg_cz(z, r, LAM, K)) - (lg_cz(0, r, LAM, 2*K) - lg_cz(0, r, LAM, K))
        dcz.append(abs(d + z*np.log(2))); dbk.append(abs(lg_book(z, r, LAM, 2*K) - lg_book(z, r, LAM, K)))
report(f'cz2009a (6.10) as printed: P_K(z)/P_K(0) = prod_(k<K) F_k(z)/F_k(0) Gamma(a_K - z)/Gamma(a_K) exactly '
       f'(r = 1..4, two couplings, three z, K = 3, 17, 60) [{idk:.0e}]; between cutoffs K = 2000 and 2K the z-dependent part '
       f'of its logarithm changes by -z ln 2 [{max(dcz):.0e}] (|z ln 2| = 0.7 at D = 0.3), a drift like (h gamma K)^-z; paired '
       f'at equal k (eq-dt-g) the change is {max(dbk):.0e} = O(1/K)', idk < 1e-9 and max(dcz) < 1e-3 and max(dbk) < 1e-3)


def fusion_map(r, qq, yy, shifts):
    """Phi: Vbar(y) -> V(y e^{i pi gamma s_1}) (x) ... (x) V(y e^{i pi gamma s_r})."""
    reps = [an_vector(r, qq, yy*np.exp(1j*np.pi*LAM*sh)) for sh in shifts]
    P = reps[0]
    for R_ in reps[1:]:
        P = an_tensor_rep(P, R_)
    B = an_vector(r, qq, yy, conj=True)
    null, _ = intertwiners([(B[(kd, j)], P[(kd, j)]) for j in range(r + 1) for kd in 'efk'])
    return null


def fused_check(r, R2, th, theta_d=0.25):
    """T of r V's at theta + i pi (2m - r + 1)/h (rightmost crosses first) on the image of Phi, against M_Vbar(theta)."""
    h = r + 1
    W = module(r, R2, s_generic(h), theta_d=theta_d); dW = W['dim']
    shifts = [(2*m - r + 1)/h for m in range(r)]
    yy = np.exp(LAM*th)
    null = fusion_map(r, q, yy, shifts)
    Phi = null[0]
    Ts = [an_T_vector(W, Xof(th + 1j*np.pi*sh, theta_d)) for sh in shifts]
    TVb = an_T_vector(W, Xof(th, theta_d), True)
    out = []
    deep = np.where(W['deep'])[0]
    for i in range(h):
        for k in deep[:3]:
            v = np.kron(Phi[:, i], np.eye(dW)[k])
            for p in range(r - 1, -1, -1):
                v = np.einsum('ij,ajb->aib', Ts[p], v.reshape(h**p, h*dW, h**(r - 1 - p))).reshape(-1)
            col = TVb[:, i*dW + k]
            target = sum(col[k2*h + j]*np.kron(np.eye(dW)[k2], Phi[:, j]) for k2 in range(dW) for j in range(h)
                         if col[k2*h + j] != 0)
            out.append(prop_residual(v, target))
    return len(null), max(o[0] for o in out), [o[1] for o in out]


fb = []
for r, R2 in ((2, 18), (3, 14)):
    for th in (0.3, -0.5 + 0.1j):
        n, res, cs = fused_check(r, R2, th)
        fb.append((n, res, max(abs(np.array(cs) - 1))))
report(f'eq-dt-bootstrap: Vbar is the image of Phi: Vbar(y) -> V(theta - i pi (r-1)/h) (x) ... (x) V(theta + i pi (r-1)/h) '
       f'(unique), and r V\'s crossing in turn act on it as M_Vbar(theta) [{max(f[1] for f in fb):.0e}] with ratio 1 '
       f'[{max(f[2] for f in fb):.0e}] (a_2, a_3): the scalar factors multiply, gbar = prod g',
       all(f[0] == 1 and f[1] < 1e-12 and f[2] < 1e-12 for f in fb))
print(f'# (d) done after {time.time() - T0:.0f} s')


# ---------------------------------------------------------------- (e) poles, residues, windows, kernel bootstrap
def sector_blocks(T, W, conj=False, interior=True):
    """Blocks of a charge-conserving T: V (x) W -> W (x) V on the sectors whose labels are all interior (or all in the
    window): (total charge Q, input indices, output indices, block)."""
    h = W['h']; dW = W['dim']; look = W['look']; I = np.eye(h); wv = an_vweights(h, conj)
    out = []
    for k, N in enumerate(W['labels']):
        if not W['interior'][k]:
            continue
        Qc = N + wv[0]
        labs = [tuple(map(int, np.round(Qc - wv[i]))) for i in range(h)]
        if all(t in look and (W['interior'][look[t]] or not interior) for t in labs):
            ins = [i*dW + look[t] for i, t in enumerate(labs)]
            outs = [look[t]*h + i for i, t in enumerate(labs)]
            out.append((Qc, ins, outs, T[np.ix_(outs, ins)]))
    return out


okp = True; outp = []
for r in (2, 3, 4):
    h = r + 1
    D1 = 1j*np.pi/2 - 1j*(np.pi*r/h + np.pi/(h*LAM))           # theta_1 - theta_d
    Dr = 1j*np.pi/2 - 1j*(np.pi/h + np.pi/(h*LAM))             # pole of gbar
    c = an_scalar_constant(r, LAM, q)
    e1, e2 = 1e-4, 1e-5
    pg = abs(e1*an_scalar_g(D1 + e1, r, LAM, c=c)/(e2*an_scalar_g(D1 + e2, r, LAM, c=c)) - 1)
    pgb = abs(e1*an_scalar_gbar(Dr + e1, r, LAM, c=c)/(e2*an_scalar_gbar(Dr + e2, r, LAM, c=c)) - 1)
    W = module(r, WIN[r] if r < 4 else 12, s_generic(h))
    X1 = Xof(D1 + 0.25, 0.25)
    bl = sector_blocks(an_T_vector(W, X1), W)
    svr = [np.linalg.svd(B, compute_uv=False) for _, _, _, B in bl]
    small = max(s[-1]/s[0] for s in svr); second = min(s[-2]/s[0] for s in svr)
    Xg = Xof(0.37 + 0.1j, 0.25)
    dets = max(abs(np.linalg.det(B)/(q**r*(1 + (-1/q)**r*Xg**h)) - 1) for _, _, _, B in sector_blocks(an_T_vector(W, Xg), W))
    okp &= pg < 1e-3 and pgb < 1e-3 and small < 1e-12 and second > 1e-3 and dets < 1e-12 and len(bl) > 0
    outp.append(f'a_{r}: {len(bl)} sectors, sigma_min {small:.0e}, next {second:.2f}')
report('eq-dt-pole: g has a simple pole at theta_1 = Delta - i(pi r/h + pi/(h gamma)) and gbar at Delta - i(pi/h + pi/(h gamma)) '
       '(Delta = theta_d + i pi/2); each charge block of M_V has det = q^r (1 + (-1/q)^r X^h), which vanishes at theta_1, '
       'so the residue has rank h - 1 in every sector: ' + '; '.join(outp), okp)


def theta_pole(a, h, gam, eta=0.4):
    return eta + 1j*np.pi/2 - 1j*(np.pi*(h - a)/h + np.pi/(h*gam))      # cz2009a (6.13), Delta = eta + i pi/2


def unstable(a, h, gam, eta=0.4):
    """Re E > 0 and Im E < 0 (strictly), for the pole in the strip -pi/2 < Im theta_a < 0 next to the real axis."""
    th = theta_pole(a, h, gam, eta)
    E = np.cosh(th)
    return E.real > 1e-14 and E.imag < -1e-14 and -np.pi/2 < th.imag < 0


gams = np.exp(np.linspace(np.log(0.3), np.log(50), 4001))      # theta_a stays in the strip Im < 2 pi
okw = True
for h in (3, 4, 5, 6):
    for a in range(1, h):
        gg = [g_ for g_ in gams if min(abs(1/g_ - a), abs(1/g_ - a + h/2)) > 1e-9]       # the ends are tested below
        pred = np.array([max(0, a - h/2) < 1/g_ < a for g_ in gg])
        okw &= np.array_equal(pred, np.array([unstable(a, h, g_) for g_ in gg]))
ends = [np.cosh(theta_pole(a, h, 1/(a - h/2))) for h in (3, 4, 5, 6) for a in range(1, h) if a > h/2]
far = np.cosh(theta_pole(1, 3, 1/6.5, 0.7))           # 1/gamma = 6.5: decaying signs, but the pole is 2 pi i further down
s1 = [g_ for g_ in gams if unstable(1, 3, g_)]; s2 = [g_ for g_ in gams if unstable(2, 3, g_)]
text = [g_ for g_ in gams if 1/g_ < 2 and not unstable(1, 5, g_)]
report(f'eq-dt-window = cz2009a (6.13)-(6.14) re-derived: for the pole in -pi/2 < Im theta_a < 0, species a is unstable iff '
       f'max(0, a - h/2) < 1/gamma < a (h = 3..6), both ends excluded: at 1/gamma = a - h/2 the pole is real, |Im E| <= '
       f'{max(abs(e.imag) for e in ends):.0e} (cz2009a (6.15) includes this end); outside that strip the signs recur (h = 3, a = 1, '
       f'1/gamma = 6.5: E/M = {far:.3f}); a_2 in book labels: species 1 for gamma > {min(s1):.3f}, species 2 for {min(s2):.3f} < '
       f'gamma < {max(s2):.3f} (the window of cz2007, whose species 1 is book species 2); the text of cz2009a, "1/gamma < (h-1)/2", '
       f'fails for a_4, a = 1, at {len(text)} grid points with 1 <= 1/gamma < 2',
       okw and max(abs(e.imag) for e in ends) < 1e-14 and far.real > 0 and far.imag < 0 and abs(min(s1) - 1) < 3e-3
       and abs(min(s2) - 0.5) < 3e-3 and abs(max(s2) - 2) < 3e-3 and len(text) > 0)

# kernel bootstrap (a_2): T_23 Rc_12 on the kernel of M_V(theta_1) in V(theta_1) (x) W, classes 0 -> 1 and 1 -> 2
def cycle(Ad, Cd, i, Qk, h):
    num = den = 1; Qv = np.array(Qk)
    for _ in range(h):
        key = (i, tuple(map(int, Qv))); num *= Cd[key]; den *= Ad[key]
        b = (i - 1) % h; Qv = Qv + np.eye(h)[i] - np.eye(h)[b]; i = b
    return num/den


def kernel_bootstrap(s, thd=0.25, th=0.6, r=2, R2=18):
    h = r + 1
    W = module(r, R2, s, theta_d=thd); dW = W['dim']; look = W['look']
    th1 = thd + 1j*np.pi/2 - 1j*(np.pi*r/h + np.pi/(h*LAM))
    kern = {}; ksv = 0
    for Qc, ins, outs, B in sector_blocks(an_T_vector(W, Xof(th1, thd)), W, interior=False):
        u, sv, vh = np.linalg.svd(B)
        vec = np.zeros(h*dW, complex); vec[ins] = vh[-1].conj()
        kern[tuple(map(int, Qc))] = vec; ksv = max(ksv, sv[-1]/sv[0])
    R = Rc(r, np.exp(LAM*th), np.exp(LAM*th1), False, False)            # V(theta) (x) V(theta_1) -> V(theta_1) (x) V(theta)
    Op = np.kron(np.eye(h), an_T_vector(W, Xof(th, thd))) @ np.kron(R, np.eye(dW))
    keys = [(Qk, j) for Qk in kern for j in range(h)]
    Bm = np.array([np.kron(kern[Qk], np.eye(h)[j]) for Qk, j in keys]).T
    el = {}; fit = 0
    for Qk in kern:
        if not W['deep'][look[tuple(map(int, np.array(Qk) - np.eye(h)[0]))]]:
            continue
        for i in range(h):
            vout = Op @ np.kron(np.eye(h)[i], kern[Qk])
            coef = np.linalg.lstsq(Bm, vout, rcond=None)[0]
            fit = max(fit, np.linalg.norm(Bm @ coef - vout)/np.linalg.norm(vout))
            for (Q2, j), cf in zip(keys, coef):
                if abs(cf) > 1e-9*np.abs(coef).max():
                    el[(i, Qk, j, Q2)] = cf
    trans = sorted(set((j - i) % h for (i, Qk, j, Q2) in el))
    A = {(i, Qk): v for (i, Qk, j, Q2), v in el.items() if i == j and Qk == Q2}
    C = {(i, Qk): v for (i, Qk, j, Q2), v in el.items() if j == (i - 1) % h}
    Q0 = tuple([1] + [0]*(h - 1))
    s_new = s + lam1(h)                                                # offset of the excited module (labels Q - e_1)
    ratios = max(abs(A[(i, Q0)]/A[(0, Q0)]/q**(-((np.array(Q0) + s)[i] - (np.array(Q0) + s)[0])) - 1) for i in range(h))
    Am = {(i, Qk): q**(-(np.array(Qk) - np.eye(h)[0] + s_new)[i]) for (i, Qk) in A}
    inv_ratio = cycle(A, C, 0, Q0, h)/cycle(Am, {kk: Xof(th, thd) for kk in C}, 0, Q0, h)
    return ksv, fit, trans, ratios, abs(inv_ratio/q**2 - 1)


kb = [kernel_bootstrap(s) for s in (np.zeros(3), lam1(3))]
report(f'prp-a2-unstable (3): at theta_1 the kernel of M_V in V(theta_1) (x) W is one-dimensional in each sector (sigma_min '
       f'{max(k[0] for k in kb):.0e}) and invariant under T_23 Rc_12 [{max(k[1] for k in kb):.0e}]; on it the transmission '
       f'has the entries v_i -> v_i, v_(i-1) only, the diagonal ratios of the module with offset s + lambda_1 '
       f'[{max(k[3] for k in kb):.0e}], and q^2 times its cycle invariant at the same theta_d [{max(k[4] for k in kb):.0e}]: '
       f'X^h -> q^2 X^h, theta_d -> theta_d + 2 pi i/h (a_2, classes 0 -> 1 and 1 -> 2)',
       all(k[0] < 1e-12 and k[1] < 1e-10 and k[2] == [0, 2] and k[3] < 1e-10 and k[4] < 1e-10 for k in kb))
# the matrix part fixes the shift of theta_d only modulo 2 pi i/(h gamma); the representatives give different candidate poles
okk = True
for h in (3, 4, 6):
    for gam in (1.37, 0.6, 7.3):
        okk &= all(abs(np.exp(-h*gam*(2j*np.pi/h + 2j*np.pi*k/(h*gam))) - q_phys(gam)**2) < 1e-12 for k in range(-2, 3))
h6 = [theta_pole(1, 6, LAM, 0.7) + 2j*np.pi/6 + 2j*np.pi*k/(6*LAM) for k in (0, 1)]       # a = p = 1 of a_5: a + 2p = h/2
E6 = [np.cosh(t) for t in h6]
report(f'the kernel fixes theta_d -> theta_d + 2 pi i/h + 2 pi i k/(h gamma) only up to the integer k: every k gives X^h -> q^2 X^h '
       f'(h = 3, 4, 6; gamma = 1.37, 0.6, 7.3; k = -2..2), and each tends to 2 pi i/h as gamma -> oo; with g(theta - shift), '
       f'species 1 of a_5 in class 1 has its pole at theta_d - i pi/(h gamma) for k = 0 (E/M = {E6[0]:.4f}) and at '
       f'theta_d + i pi/(h gamma) for k = 1 (E/M = {E6[1]:.4f}): opposite signs of Im E',
       okk and abs(h6[0] - (0.7 - 1j*np.pi/(6*LAM))) < 1e-12 and abs(h6[1] - (0.7 + 1j*np.pi/(6*LAM))) < 1e-12
       and E6[0].imag < 0 < E6[1].imag)
print(f'# (e) done after {time.time() - T0:.0f} s')


# ---------------------------------------------------------------- (f) semiclassics
def zcl(u, a, h):
    """eq-dc-ar-delay: z_a(u) = (e^-u - i e^{i pi a/h})/(e^-u - i e^{-i pi a/h})."""
    return (np.exp(-u) - 1j*np.exp(1j*np.pi*a/h))/(np.exp(-u) - 1j*np.exp(-1j*np.pi*a/h))


devs = {}
for r in (2, 3):
    h = r + 1
    for gam in (10.0, 40.0, 160.0):
        dv = db = 0
        for D in (-0.7, 0.5, 1.1):
            e = 1e-5
            ph = np.angle(an_scalar_g(D + e, r, gam)/an_scalar_g(D - e, r, gam))/(2*e)
            phb = np.angle(an_scalar_gbar(D + e, r, gam)/an_scalar_gbar(D - e, r, gam))/(2*e)
            pred = np.log(abs(zcl(D, 1, h)))
            dv = max(dv, abs(2*np.pi/(h*(gam + 1))*ph - pred)); db = max(db, abs(2*np.pi/(h*(gam + 1))*phb - pred))
        devs[(r, gam)] = (dv, db)
okf = all(devs[(r, 160.0)][0] < 0.01 and devs[(r, 160.0)][1] < 0.01 and 3 < devs[(r, 10.0)][0]/devs[(r, 40.0)][0] < 5
          and 3 < devs[(r, 40.0)][0]/devs[(r, 160.0)][0] < 5 for r in (2, 3))
argz = max(abs(np.angle(zcl(0.0, 1, h)) - np.pi/h) for h in (3, 4, 5))
thr = 0                                                     # exr-dt-semiclassical: arg z_1(u_*) = eps at the threshold u_*
for h in (3, 4, 5, 8):
    for eps in (0.03, 0.2, np.pi/(2*h), np.pi/h):
        ustar = -0.5*np.log(np.sin(2*np.pi/h)/np.tan(eps) - np.cos(2*np.pi/h))
        thr = max(thr, abs(np.angle(zcl(ustar, 1, h)) - eps), max(ustar, 0))
    us = np.linspace(-6, 6, 1201)
    thr = max(thr, 0 if np.all(np.diff(np.angle(zcl(us, 1, h))) > 0) else 1)
report('prp-transmission-semiclassical: (2 pi/(h(gamma + 1))) d arg g/d theta -> ln|z_1(theta - theta_d)| (eq-dc-ar-delay), '
       'deviations (g, gbar) at gamma = 10, 40, 160: '
       + '; '.join(f'a_{r}: ' + ', '.join(f'{devs[(r, gm)][0]:.4f}/{devs[(r, gm)][1]:.4f}' for gm in (10.0, 40.0, 160.0))
                   for r in (2, 3)) + f'; arg z_1(theta_d) = pi/h [{argz:.0e}]; arg z_1 increases with theta, and a soliton a '
       f'distance eps inside the edge of its sector converts iff theta - theta_d > -ln(sin(2 pi/h)/tan eps - cos(2 pi/h))/2 <= 0 '
       f'(h = 3, 4, 5, 8) [{thr:.0e}]', okf and argz < 1e-14 and thr < 1e-12)
print(f'# (f) done after {time.time() - T0:.0f} s')


# ---------------------------------------------------------------- (g) a_2^(2) (Corrigan-Zambon 2011, Sec. 4)
CART22 = [[2, -1], [-4, 2]]          # CZ 2011 (4.1): a_01 = -1, a_10 = -4
D22 = [4, 1]                         # d_0 = 4, d_1 = 1


def qbin(m, k, t):
    num = np.prod([qnum(m - s_, t) for s_ in range(k)]) if k else 1
    den = np.prod([qnum(s_ + 1, t) for s_ in range(k)]) if k else 1
    return num/den


def rel22(rep, qq, kinds=('e', 'f'), mask=None):
    """Relations of U_q(a_2^(2)) in Part VI form: k_i e_j k_i^-1 = q^{d_i a_ij} e_j, [e_i, f_j] = delta_ij (k_i - k_i^-1)/(q_i - q_i^-1),
    q-Serre with q_i = q^{d_i} (orders 2 and 5).  Relative to the size of the terms."""
    sub = (lambda M: M[np.ix_(mask, mask)]) if mask is not None else (lambda M: M)
    res = 0
    for i in range(2):
        qi = qq**D22[i]; ki = rep[('k', i)]; kii = np.linalg.inv(ki)
        for j in range(2):
            for X, sg in (('e', 1), ('f', -1)):
                if X in kinds:
                    res = max(res, np.abs(sub(ki @ rep[(X, j)] @ kii - qq**(sg*D22[i]*CART22[i][j])*rep[(X, j)])).max()
                              / np.abs(sub(rep[(X, j)])).max())
            if 'e' in kinds and 'f' in kinds:
                cmm = rep[('e', i)] @ rep[('f', j)] - rep[('f', j)] @ rep[('e', i)]
                want = (ki - kii)/(qi - 1/qi) if i == j else 0*cmm
                res = max(res, np.abs(sub(cmm - want)).max())
            if i != j:
                m = 1 - CART22[i][j]
                for X in kinds:
                    Xi, Xj = rep[(X, i)], rep[(X, j)]
                    terms = [(-1)**t*qbin(m, t, qi)*np.linalg.matrix_power(Xi, m - t) @ Xj @ np.linalg.matrix_power(Xi, t)
                             for t in range(m + 1)]
                    den = max(np.abs(sub(tt)).max() for tt in terms)
                    if den > 0:                      # (on the 3-dim representation all quintic terms vanish)
                        res = max(res, np.abs(sub(sum(terms))).max()/den)
    return res


def vec3(qq, x):
    """CZ 2011 (4.5)-(4.6), homogeneous gradation (4.10) on node 0, in Part VI form e_i = K_i X_i^+, f_i = X_i^- K_i^-1, k_i = K_i^2."""
    K1 = np.diag([qq, 1, 1/qq]).astype(complex)
    X1p = np.array([[0, 1, 0], [0, 0, 1], [0, 0, 0]], complex)
    X1m = (qq + 1/qq)*np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], complex)
    K0 = np.diag([qq**-2, 1, qq**2]).astype(complex)
    X0p = x*np.array([[0, 0, 0], [0, 0, 0], [1, 0, 0]], complex)
    X0m = np.array([[0, 0, 1], [0, 0, 0], [0, 0, 0]], complex)/x
    g = {('e', 1): K1 @ X1p, ('f', 1): X1m @ np.linalg.inv(K1), ('k', 1): K1 @ K1,
         ('e', 0): K0 @ X0p, ('f', 0): X0m @ np.linalg.inv(K0), ('k', 0): K0 @ K0}
    g.update(dim=3, hw=0)
    return g


def borel22(qq, N, f1, f2, F=None):
    """CZ 2011 (4.7)-(4.9): X1+ = a^dag, X0+ = a a, K1 = q^N, K0 = q^-2N, a|j> = F(j)|j-1>, F = f1 q^{-2N-1} + f2 q^{2N+1};
    U_q(b^+) in Part VI form.  Labels j = -N..N (the topological charge)."""
    js = np.arange(-N, N + 1); d = len(js)
    Ff = (lambda j: f1*qq**(-2*j - 1) + f2*qq**(2*j + 1)) if F is None else F
    a = np.zeros((d, d), complex); ad = np.zeros((d, d), complex)
    for i in range(d - 1):
        ad[i + 1, i] = 1; a[i, i + 1] = Ff(js[i + 1])
    K1 = np.diag(qq**js.astype(float)).astype(complex); K0 = np.diag(qq**(-2.0*js)).astype(complex)
    g = {('e', 1): K1 @ ad, ('k', 1): K1 @ K1, ('e', 0): K0 @ a @ a, ('k', 0): K0 @ K0}
    g.update(dim=d, hw=None, n=js)
    return g


M3 = np.array([1, 0, -1])


def solve22(V, W, margin=3, tol=1e-7):
    """U_q(b^+) intertwiners V (x) W -> W (x) V (charge m + j conserved); returns (rank on |j| <= N - margin, T there)."""
    dW = W['dim']; js = W['n']; N = js.max(); din = 3*dW
    unk = [(int(j + M3[a] - M3[b] + N)*3 + b, a*dW + k) for a in range(3) for k, j in enumerate(js) for b in range(3)
           if abs(j + M3[a] - M3[b]) <= N]
    nu = len(unk)
    inn = np.abs(js) <= N - 2
    rowmask = np.repeat(inn, 3); colmask = np.tile(inn, 3)
    rows = []
    for i in (0, 1):
        A = cop(V, W, 'e', i); B = cop(W, V, 'e', i)
        Mf = np.zeros((dW*3, din, nu), complex)
        for u, (r0, c0) in enumerate(unk):
            Mf[r0, :, u] += A[c0, :]; Mf[:, c0, u] -= B[:, r0]
        rows.append(Mf[np.ix_(rowmask, colmask)].reshape(-1, nu))
    _, sv, Vh = np.linalg.svd(np.vstack(rows), full_matrices=False)
    null = Vh[sv/sv[0] < tol].conj()
    deep = np.abs(js) <= N - margin
    dmask = np.array([deep[r0 // 3] and deep[c0 % dW] for r0, c0 in unk])
    Uu, s2, _ = np.linalg.svd(null[:, dmask], full_matrices=False)
    vec = Uu[:, 0].conj() @ null
    T = np.zeros((dW*3, din), complex)
    for (r0, c0), val in zip(unk, vec):
        if deep[r0 // 3] and deep[c0 % dW]:
            T[r0, c0] = val
    return int((s2 > 1e-6*s2[0]).sum()), T


q22 = np.exp(-0.77j)
res3 = rel22(vec3(q22, 1.3), q22)
mods = [rel22(borel22(q22, 8, f1, f2), q22, kinds=('e',), mask=np.abs(np.arange(-8, 9)) <= 2)
        for f1, f2 in ((1.0, 0.7), (1.0, 0.0), (0.0, 1.0), (0.4 + 0.3j, -1.1))]
ctl = [rel22(borel22(q22, 8, 0, 0, F=F), q22, kinds=('e',), mask=np.abs(np.arange(-8, 9)) <= 2)
       for F in (lambda j: 1.0, lambda j: (-1.0)**j)]
ranks = [solve22(vec3(q22, 0.8), borel22(q22, 9, f1, f2))[0] for f1, f2 in ((1.0, 0.7), (1.0, 0.0), (0.0, 1.0), (0.4 + 0.3j, -1.1))]
report(f'a_2^(2): CZ 2011 algebra in Part VI form, 3-dim representation [{res3:.0e}]; Borel module (4.7)-(4.9) satisfies '
       f'the U_q(b^+) relations incl. the quintic q-Serre relation [{max(mods):.0e}], F = 1 and F = (-1)^N do not '
       f'[{min(ctl):.2f}] (no Weyl, type-I, module); intertwiners V (x) W -> W (x) V unique for four (f1, f2) (ranks {ranks})',
       res3 < 1e-12 and max(mods) < 1e-12 and min(ctl) > 0.1 and ranks == [1, 1, 1, 1])

# triangle relation with the Izergin-Korepin R-check of the 3-dim representation
N22 = 12
okt22 = True; tr22 = []
for f1, f2 in ((1.0, 0.7), (1.0, 0.0)):
    W = borel22(q22, N22, f1, f2)
    x1, x2 = 0.8, 1.9 + 0.3j
    _, T1 = solve22(vec3(q22, x1), W); _, T2 = solve22(vec3(q22, x2), W)
    R, nul = rcheck(vec3(q22, x1), vec3(q22, x2))
    dW = W['dim']; js = W['n']; IW = np.eye(dW)
    L = np.kron(IW, R) @ np.kron(T1, I3) @ np.kron(I3, T2)
    Rr = np.kron(T2, I3) @ np.kron(I3, T1) @ np.kron(R, IW)
    blk = np.ix_(np.abs(np.repeat(js, 9)) <= N22 - 6, np.abs(np.tile(js, 9)) <= N22 - 6)
    res, _ = prop_residual(L[blk], Rr[blk])
    rng = np.random.default_rng(1)
    bad, _ = prop_residual((np.kron(IW, R + 0.05*rng.standard_normal(R.shape)) @ np.kron(T1, I3) @ np.kron(I3, T2))[blk], Rr[blk])
    tr22.append((res, bad)); okt22 &= nul == 1 and res < 1e-12 and bad > 1e-3
report(f'a_2^(2): the intertwiners satisfy eq-triangle-braid with rcheck(V(x1), V(x2)) [{max(t[0] for t in tr22):.0e}]; '
       f'control with a perturbed R-check [{min(t[1] for t in tr22):.2f}]', okt22)


def conj22(qq):
    xc = -qq**-6
    A = dual(vec3(qq, 0.9)); B = vec3(qq, 0.9*xc)
    null, _ = intertwiners([(A[(kd, j)], B[(kd, j)]) for j in (0, 1) for kd in 'efk'])
    return xc, null


def block_inv22(T, W, N, margin):
    dW = W['dim']
    Tinv = np.zeros((3*dW, dW*3), complex); done = np.zeros(3*dW, bool)
    for Qc in range(-N, N + 1):
        mem = [(a, Qc - M3[a]) for a in range(3) if abs(Qc - M3[a]) <= N - margin]
        if len(mem) < 3:
            continue
        ins = [a*dW + int(j + N) for a, j in mem]; outs = [int(j + N)*3 + a for a, j in mem]
        Tinv[np.ix_(ins, outs)] = np.linalg.inv(T[np.ix_(outs, ins)]); done[ins] = True
    return Tinv, done


def crossing22(qq, f1, f2, x, xc, C, N=10, margin=3):
    W = borel22(qq, N, f1, f2)
    _, T = solve22(vec3(qq, x), W, margin); _, Ts = solve22(vec3(qq, x/xc), W, margin)
    k0 = N; T = T/T[k0*3, k0]; Ts = Ts/Ts[k0*3, k0]
    Tinv, done = block_inv22(Ts, W, N, margin)
    ok = np.abs(W['n']) <= N - margin - 2
    return defect_crossing_residual(T, Tinv, done, ok, W['dim'], 3, C)


xc, Cn = conj22(q22)
xcm, Cnm = conj22(-q22)
okx = len(Cn) == 1 and len(Cnm) == 1 and np.allclose(np.abs(Cn[0]) > 1e-9, np.fliplr(np.eye(3)) > 0)
cr = {}
for name, qq, Cm_ in (('q, C(q)', q22, Cn[0]), ('-q, C(-q)', -q22, Cnm[0]), ('-q, C(q)', -q22, Cn[0])):
    cr[name] = max(crossing22(qq, 1.0, 0.7, x, -qq**-6, Cm_)[0] for x in (0.8, 1.6, 0.5 + 0.4j))
report(f'a_2^(2): dual V(x) = V(x x_c), x_c = -q^-6, with an antidiagonal C (unique); crossing eq-defect-crossing with '
       f'x(theta - i pi) = x/x_c holds with C [{cr["q, C(q)"]:.0e}], at -q with C(-q) [{cr["-q, C(-q)"]:.0e}], and fails '
       f'at -q with C(q) [{cr["-q, C(q)"]:.2f}]', okx and cr['q, C(q)'] < 1e-10 and cr['-q, C(-q)'] < 1e-10
       and cr['-q, C(q)'] > 0.1)
print(f'# (g) done after {time.time() - T0:.0f} s')


# ---------------------------------------------------------------- (h) selection rule and tunnelling
okh = True; outh = []
for r in (2, 3, 4):
    h = r + 1
    W = module(r, WIN[r], s_generic(h)); look = W['look']; k0 = look[tuple([0]*h)]; I = np.eye(h)
    X = Xof(0.6, 0.25)
    TV, TVb = an_T_vector(W, X), an_T_vector(W, X, True)
    nV = sum(abs(TV[look[tuple(map(int, I[0] - I[b]))]*h + b, k0]) > 0 for b in range(h))
    mods_ = [abs(TVb[look[tuple(map(int, -I[0] + I[b]))]*h + b, k0]) for b in range(h)]
    okh &= nV == 2 and max(abs(mods_[j] - abs(X)**j) for j in range(h)) < 1e-12
    outh.append(f'a_{r}: V {nV}, Vbar {sum(m_ > 0 for m_ in mods_)}')
supp = []
for gam in (1.37, 10.0, 40.0):
    for D in (0.3, -0.3, 0.0):
        x_ = np.exp(gam*D); h = 4
        supp.append(max(x_**j for j in range(1, h - 1))/max(1, x_**(h - 1)))
okh &= supp[2] == 1 and supp[5] == 1 and supp[6] < 1e-4 and supp[7] < 1e-4 and supp[6] < supp[3] < supp[0]
report('prp-an-transmission (4), selection rule: entries per input ' + ', '.join(outh) + '; |Vbar entry j| = |X|^j; '
       f'the classically forbidden entries (1 <= j <= h - 2, a_3) relative to the largest: at theta - theta_d = +-0.3, '
       f'gamma = 1.37, 10, 40: {supp[0]:.2f}, {supp[3]:.0e}, {supp[6]:.0e}; at theta = theta_d: 1', okh)
nonu = []
for r in (2, 3):
    h = r + 1
    W = module(r, WIN[r], s_generic(h))
    for X in (Xof(0.6, 0.25), Xof(-0.4, 0.25)):
        for Qc, ins, outs, B in sector_blocks(an_T_vector(W, X), W)[:3]:
            G = B.conj().T @ B
            nonu.append(np.abs(G - np.diag(np.diag(G))).max()/abs(G[0, 0]))
report(f'for r >= 2 the matrix part is not unitary for real theta: M^+ M has off-diagonal entries of relative size '
       f'{min(nonu):.2f} or more in every block tested (a_2, a_3), unlike sine-Gordon (ch. 37)', min(nonu) > 0.1)
print(f'# (h) done after {time.time() - T0:.0f} s')


# ---------------------------------------------------------------- (i) uniqueness for every r: sector machinery
# A charge-conserving map P (x) W -> W (x) P is a family of d x d blocks between the sectors of total charge Q (states
# p (x) |Q - wt(p)>).  Exact version: Fractions at a rational point (q = 2), integer labels N (sum free), offsets t_c = q^{s_c}.
from fractions import Fraction as Fr
import itertools


def qp(qq, n):
    return qq**n if n >= 0 else 1/qq**(-n)


def alph(h, j):
    a = [0]*h
    if j == 0:
        a[h - 1] += 1; a[0] -= 1
    else:
        a[j - 1] += 1; a[j] -= 1
    return tuple(a)


def mmul(A, B):
    return [[sum((A[i][k]*B[k][j] for k in range(len(B))), 0*A[0][0]) for j in range(len(B[0]))] for i in range(len(A))]


def minv(A):
    n = len(A); M = [list(A[i]) + [1 if i == j else 0 for j in range(n)] for i in range(n)]
    for c in range(n):
        pv = next(i for i in range(c, n) if M[i][c] != 0); M[c], M[pv] = M[pv], M[c]
        pc = M[c][c]; M[c] = [x/pc for x in M[c]]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f_ = M[i][c]; M[i] = [a - f_*b for a, b in zip(M[i], M[c])]
    return [row[n:] for row in M]


def xrank(rows, nc):
    M = [list(r) for r in rows if any(x != 0 for x in r)]; rk = 0
    for c in range(nc):
        pv = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if pv is None:
            continue
        M[rk], M[pv] = M[pv], M[rk]; pc = M[rk][c]; M[rk] = [x/pc for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f_ = M[i][c]; M[i] = [a - f_*b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk


def prep(h, qq, yy, conj):
    """V(y) or Vbar(y) as in an_vector: node j moves v_j <-> v_{j+1} (node 0: v_h <-> v_1)."""
    zero = 0*qq
    R = {'dim': h, 'wt': [tuple((-1 if conj else 1)*(c == i) for c in range(h)) for i in range(h)]}
    for j in range(h):
        a, b = (h - 1, 0) if j == 0 else (j - 1, j)
        F = [[zero]*h for _ in range(h)]; E = [[zero]*h for _ in range(h)]
        if conj:
            F[a][b] = 1/yy; E[b][a] = yy
        else:
            F[b][a] = 1/yy; E[a][b] = yy
        al = alph(h, j)
        K_ = [[zero]*h for _ in range(h)]
        for i in range(h):
            K_[i][i] = qp(qq, sum(x*w for x, w in zip(al, R['wt'][i])))
        R[('f', j)], R[('e', j)], R[('k', j)] = F, E, K_
    return R


def ptensor(A, B):
    dA, dB = A['dim'], B['dim']
    kron = lambda X, Y: [[X[i1][j1]*Y[i2][j2] for j1 in range(dA) for j2 in range(dB)] for i1 in range(dA) for i2 in range(dB)]
    IA = [[1 if i == j else 0 for j in range(dA)] for i in range(dA)]; IB = [[1 if i == j else 0 for j in range(dB)] for i in range(dB)]
    R = {'dim': dA*dB, 'wt': [tuple(x + y for x, y in zip(wa, wb)) for wa in A['wt'] for wb in B['wt']]}
    for j in range(len(A['wt'][0])):
        add = lambda X, Y: [[x + y for x, y in zip(r1, r2)] for r1, r2 in zip(X, Y)]
        R[('e', j)] = add(kron(A[('e', j)], IB), kron(A[('k', j)], B[('e', j)]))
        R[('f', j)] = add(kron(A[('f', j)], minv(B[('k', j)])), kron(IA, B[('f', j)]))
        R[('k', j)] = kron(A[('k', j)], B[('k', j)])
    return R


class Wmod:
    """eq-dt-weyl with integer labels N: x_j|N> = m q^{orient n_{j+1}} |N + sign alpha_j>, k_j|N> = q^{alpha_j.n}, q^{n_c} = t_c q^{N_c};
    x = f (sign -1, U_q(b^-)) or e (sign +1, U_q(b^+))."""
    def __init__(self, h, qq, m, t, sign, orient=1):
        self.h, self.q, self.m, self.t, self.sign, self.orient = h, qq, m, t, sign, orient

    def qn(self, N, c):
        return self.t[c]*qp(self.q, N[c])

    def k(self, j, N):
        out = 1
        for c, x in enumerate(alph(self.h, j)):
            if x:
                out = out*self.qn(N, c)**x
        return out

    def x(self, j, N):
        al = alph(self.h, j); c2 = al.index(-1)
        return self.m*self.qn(N, c2)**self.orient, tuple(n + self.sign*x for n, x in zip(N, al))


def smap(P, W, j, Q, pfirst):
    """Block of Delta(x_j) from the sector Q of P (x) W (pfirst) or W (x) P to the sector Q + sign alpha_j."""
    al = alph(W.h, j); Qn = tuple(c + W.sign*x for c, x in zip(Q, al)); d = P['dim']
    src = [(p, tuple(c - w for c, w in zip(Q, P['wt'][p]))) for p in range(d)]
    look = {(p, tuple(c - w for c, w in zip(Qn, P['wt'][p]))): p for p in range(d)}
    M = [[0*W.q]*d for _ in range(d)]
    X = 'f' if W.sign == -1 else 'e'
    for i, (p, N) in enumerate(src):
        cW, N2 = W.x(j, N); kW = W.k(j, N); kP = P[('k', j)][p][p]
        for p2 in range(d):
            a_ = P[(X, j)][p2][p]
            if a_ != 0:
                # f: f x k^-1 + 1 x f;  e: e x 1 + k x e  (particle first), and the mirror forms for W (x) P
                M[look[(p2, N)]][i] += (a_/kW if pfirst else a_) if X == 'f' else (a_ if pfirst else kW*a_)
        M[look[(p, N2)]][i] += (cW if pfirst else cW/kP) if X == 'f' else (kP*cW if pfirst else cW)
    return M, Qn


def spath(P, W, Q, nodes, pfirst):
    d = P['dim']; M = [[1 if i == j else 0 for j in range(d)] for i in range(d)]
    for j in nodes:
        A, Q = smap(P, W, j, Q, pfirst); M = mmul(A, M)
    return M, Q


def loop_rows(Ps, Pt, W, Q):
    """Linear conditions L_B T = T L_A on the block T_Q, for the commutator loops of two nodes and the loop through all nodes."""
    h = W.h; d = Ps['dim']; rows = []
    loops = [(list(jk), list(jk[::-1])) for jk in itertools.combinations(range(h), 2)] + [(list(range(h)), None)]
    for a_, b_ in loops:
        LA, LB = spath(Ps, W, Q, a_, True)[0], spath(Pt, W, Q, a_, False)[0]
        if b_ is not None:
            LA = mmul(minv(spath(Ps, W, Q, b_, True)[0]), LA); LB = mmul(minv(spath(Pt, W, Q, b_, False)[0]), LB)
        for i in range(d):
            for k in range(d):
                r_ = [0*W.q]*(d*d)
                for c in range(d):
                    r_[c*d + k] += LB[i][c]; r_[i*d + c] -= LA[c][k]
                rows.append(r_)
    return rows


def closed_block(h, qq, X, t, Q, kind):
    """Blocks of eq-dt-tv ('V'), eq-dt-tvbar ('Vbar'), eq-dt-tplus ('Vbar+'), eq-dt-reversed ('rev') on the sector Q."""
    T = [[0*qq]*h for _ in range(h)]
    for i in range(h):
        N = tuple(Q[c] + (c == i)*(-1 if kind == 'V' else 1) for c in range(h))
        qn = lambda c: t[c]*qp(qq, N[c])
        if kind == 'V':
            T[i][i] = 1/qn(i); T[(i - 1) % h][i] = X
        elif kind == 'Vbar':
            pw = 1
            for j in range(h):
                b = (i + j) % h; pw = pw*qn(b); T[b][i] = X**j*pw
        elif kind == 'Vbar+':
            b = (i - 1) % h; T[i][i] = 1/qn(i); T[b][i] = X*qn(i)/qn(b)/qq
        else:
            b = (i + 1) % h; T[i][i] = qn(i); T[b][i] = X*qn(i)/qn(b)
    return T


t1 = time.time()
qx, mx, yx = Fr(2), Fr(3, 7), Fr(5, 4)
TT = {3: (Fr(3, 2), Fr(5, 7), Fr(14, 15)), 4: (Fr(3, 2), Fr(5, 7), Fr(2, 3), Fr(7, 5))}
nul = []
for h in (3, 4):
    for t in ((Fr(1),)*h, TT[h]):
        for sign in (-1, 1):
            for conj in (False, True):
                P = prep(h, qx, yx, conj)
                nul.append(h*h - xrank(loop_rows(P, P, Wmod(h, qx, mx, t, sign), tuple([1] + [0]*(h - 1))), h*h))
nul2 = []
for t in ((Fr(1),)*3, TT[3]):
    P12 = ptensor(prep(3, qx, yx, False), prep(3, qx, Fr(11, 3), False))
    P21 = ptensor(prep(3, qx, Fr(11, 3), False), prep(3, qx, yx, False))
    nul2.append(81 - xrank(loop_rows(P12, P21, Wmod(3, qx, mx, t, -1), (1, 1, -2)), 81))
okx = True                           # existence: the closed forms satisfy the sector equations exactly
for h in (3, 4, 5):
    t = (Fr(3, 2), Fr(5, 7)) + (Fr(1),)*(h - 3) + (Fr(14, 15),)
    for kind, sign, conj in (('V', -1, False), ('Vbar', -1, True), ('Vbar+', 1, True)):
        W = Wmod(h, qx, mx, t, sign); P = prep(h, qx, yx, conj)
        X = -(qx - 1/qx)*mx*yx if sign == -1 else (qx - 1/qx)*mx/yx
        for Q in (tuple([1] + [0]*(h - 1)), tuple([2, -1] + [0]*(h - 2))):
            for j in range(h):
                A, Qn = smap(P, W, j, Q, True); B, _ = smap(P, W, j, Q, False)
                okx &= mmul(closed_block(h, qx, X, t, Qn, kind), A) == mmul(B, closed_block(h, qx, X, t, Q, kind))
report(f'prp-an-transmission (3), exact: on one sector the loop conditions (commutator loops and the loop through all nodes) '
       f'leave a one-dimensional space for a_2 and a_3, both types, V and Vbar, s = 0 and a generic offset (nullities {set(nul)}; '
       f'q = 2, rational parameters, so by semicontinuity for generic ones); hypothesis (ii) of thm-t-triangle for a_2, '
       f'V(y1) (x) V(y2) (x) W -> W (x) V(y2) (x) V(y1): nullity {set(nul2)} (81 unknowns); eq-dt-tv, eq-dt-tvbar, eq-dt-tplus satisfy '
       f'the sector equations exactly (h = 3, 4, 5) [{time.time() - t1:.0f} s]', set(nul) == {1} and set(nul2) == {1} and okx)

# the proof for every r: Omega_j has no pure-W term, image/kernel are coordinate lines/hyperplanes (random complex parameters)
rng = np.random.default_rng(11); omok = True
for h in (3, 4, 5, 6):
    for sign in (-1, 1):
        for conj in (False, True):
            qc = np.exp(1j*rng.uniform(0.3, 1.2))*rng.uniform(0.8, 1.2)
            tc = np.exp(rng.normal(size=h) + 1j*rng.normal(size=h)); tc = tc/np.prod(tc)**(1/h)
            W = Wmod(h, qc, complex(*rng.normal(size=2)), tuple(tc), sign); P = prep(h, qc, complex(*rng.normal(size=2)), conj)
            Q = tuple([1] + [0]*(h - 1))
            for j in range(h):
                k = (j + 1) % h; tgt = (j + 1) % h                       # 0-based index of v_{j+2}
                A1 = np.array(spath(P, W, Q, [k, j], True)[0]); A2 = np.array(spath(P, W, Q, [j, k], True)[0])
                Om = A1 - (A2/qc if sign == -1 else qc*A2)
                sv = np.linalg.svd(Om, compute_uv=False)
                rows_ = [i for i in range(h) if np.abs(Om[i]).max() > 1e-12*np.abs(Om).max()]
                cols_ = [i for i in range(h) if np.abs(Om[:, i]).max() > 1e-12*np.abs(Om).max()]
                line = (sign == -1) != conj                                 # image a line for (b-, V) and (b+, Vbar)
                omok &= sv[1] < 1e-12*sv[0] and (rows_ == [tgt] if line else cols_ == [tgt])
# sector determinants (projected labels n = N + s): independent of the sector and of s
qd = q; dets = 0
for h in (3, 4, 5):
    r_ = h - 1; X_ = 0.83*np.exp(0.4j)
    for s_ in (np.zeros(h), s_generic(h)):
        for N0 in (np.zeros(h), np.eye(h)[0] - np.eye(h)[1], np.eye(h)[2] - np.eye(h)[0]):
            mu_ = lambda i: np.eye(h)[i] - 1/h
            for kind, form in (('V', qd**r_*(1 + (-1/qd)**r_*X_**h)), ('Vbar', (qd - X_**h)**r_),
                               ('Vbar+', qd**(-r_)*(1 + (-qd)**r_*X_**h)), ('rev', qd**r_*(1 + (-1)**r_*qd*X_**h))):
                B = np.zeros((h, h), complex)
                for i in range(h):
                    n = (N0 + s_ - mu_(i) + mu_(0)) if kind == 'V' else (N0 + s_ + mu_(i) - mu_(0))
                    if kind == 'V':
                        B[i, i] = qd**(-n[i]); B[(i - 1) % h, i] = X_
                    elif kind == 'Vbar':
                        pw = 1
                        for j in range(h):
                            b = (i + j) % h; pw = pw*qd**n[b]; B[b, i] = X_**j*pw
                    elif kind == 'Vbar+':
                        b = (i - 1) % h; B[i, i] = qd**(-n[i]); B[b, i] = X_*qd**(n[i] - n[b] - 1)
                    else:
                        b = (i + 1) % h; B[i, i] = qd**n[i]; B[b, i] = X_*qd**(n[i] - n[b])
                dets = max(dets, abs(np.linalg.det(B)/form - 1))
report(f'prp-an-transmission (3), the proof for every r: Omega_j = Delta(x_j)Delta(x_(j+1)) - q^-+1 Delta(x_(j+1))Delta(x_j) has '
       f'rank 1 on each sector, with image the line of v_(j+2) (b-, V; b+, Vbar) or kernel the hyperplane i != j+2 (b-, Vbar; '
       f'b+, V), h = 3..6 at random complex parameters; the sector determinants are q^r(1 + (-1/q)^r X^h) (eq-dt-tv), '
       f'(q - X^h)^r (eq-dt-tvbar), q^-r(1 + (-q)^r X_+^h) (eq-dt-tplus) and q^r(1 + (-1)^r q X^h) (eq-dt-reversed), '
       f'independent of the sector and of s (h = 3, 4, 5) [{dets:.0e}]', omok and dets < 1e-12)

# type I_+, V: the inverse of eq-dt-reversed (offset -s, labels -N, X -> -X_+), carried over by the Chevalley involution
okv = True; okr = True
for h in (3, 4, 5, 6):
    t = (Fr(3, 2), Fr(5, 7)) + (Fr(1),)*(h - 3) + (Fr(14, 15),); tinv = tuple(1/x for x in t)
    Wp = Wmod(h, qx, mx, t, +1); P = prep(h, qx, yx, False); Xp_ = -(qx - 1/qx)*mx/yx
    TQ = lambda Q: minv(closed_block(h, qx, Xp_, tinv, tuple(-c for c in Q), 'rev'))
    for Q in (tuple([1] + [0]*(h - 1)), tuple([2, -1] + [0]*(h - 2))):
        for j in range(h):
            A, Qn = smap(P, Wp, j, Q, True); B, _ = smap(P, Wp, j, Q, False)
            okv &= mmul(TQ(Qn), A) == mmul(B, TQ(Q))
        okv &= all(x != 0 for row in TQ(Q) for x in row)
    if h <= 5:                         # eq-dt-reversed itself, on the reversed-orientation module of U_q(b^-)
        Wr = Wmod(h, qx, mx, t, -1, orient=-1); Pb = prep(h, qx, yx, True); Xr = -(qx - 1/qx)*mx*yx
        for Q in (tuple([1] + [0]*(h - 1)), tuple([2, -1] + [0]*(h - 2))):
            for j in range(h):
                A, Qn = smap(Pb, Wr, j, Q, True); B, _ = smap(Pb, Wr, j, Q, False)
                okr &= mmul(closed_block(h, qx, Xr, t, Qn, 'rev'), A) == mmul(B, closed_block(h, qx, Xr, t, Q, 'rev'))
gauge = all(sum((n_ - a_)**2 for n_, a_ in zip(n, alph(h, i))) - sum(x*x for x in n) == -2*(sum(n_*a_ for n_, a_ in zip(n, alph(h, i)))) + 2
            for h in (3, 4, 5) for i in range(h) for n in itertools.product(range(-2, 3), repeat=h))
report('eq-dt-reversed intertwines U_q(b^-) on the reversed module (exact, h = 3, 4, 5); in the basis q^(-|n|^2/2)|n> its '
       'converting coefficient is q X, since |n - alpha_i|^2 - |n|^2 = -2 alpha_i.n + 2 (cz2007 (5.10) with xhat = qX); '
       'prp-an-transmission (2): the inverse of its blocks, with offset -s, labels -N and X -> -X_+, is the type-I_+ intertwiner '
       'for V, all h^2 entries nonzero (exact, h = 3..6)', okr and okv and gauge)

# root of unity: the shift by 5 alpha_1 at q^5 = 1, and uniqueness among charge-conserving maps
q5 = np.exp(2j*np.pi/5); W5 = Wmod(3, q5, 0.7 - 0.2j, (1.0, 1.0, 1.0), -1); sres = 0
for N in ((1, -1, 0), (0, 2, -2), (3, 0, -3)):
    NS = tuple(n + 5*a for n, a in zip(N, alph(3, 1)))
    for j in range(3):
        sres = max(sres, abs(W5.x(j, N)[0] - W5.x(j, NS)[0]), abs(W5.k(j, N) - W5.k(j, NS)))
sv5 = []
for qq_ in (q5, q_phys(0.6), q_phys(LAM)):
    for h in (3, 4):
        for sign in (-1, 1):
            for conj in (False, True):
                P = prep(h, qq_, 0.8 + 0.3j, conj)
                M_ = np.array(loop_rows(P, P, Wmod(h, qq_, 0.43 - 0.21j, (1.0,)*h, sign), tuple([1] + [0]*(h - 1))), complex)
                sv_ = np.linalg.svd(M_, compute_uv=False); sv5.append((sv_[-1]/sv_[0], sv_[-2]/sv_[0]))
report(f'roots of unity: at q^5 = 1 (gamma = 0.6) the shift |N> -> |N + 5 alpha_1> commutes with k_j and f_j [{sres:.0e}], so the '
       f'k_j-equations do not enforce charge conservation; among charge-conserving maps the loop conditions still leave one '
       f'dimension at q^5 = 1 and at gamma = 1.37 (q^200 = 1), h = 3, 4, both types, V and Vbar (smallest relative singular value '
       f'{max(x[0] for x in sv5):.0e}, next {min(x[1] for x in sv5):.2f})',
       sres < 1e-14 and max(x[0] for x in sv5) < 1e-12 and min(x[1] for x in sv5) > 0.05)
# prp-an-transmission (4) for every r: the universal R-matrix.  First-order terms: R = q^t (1 + (q - q^-1) sum_j f_j (x) e_j
# + ...) on Vbar/V (x) W^+ (type I_+) and R_21^{-1} = (1 - (q - q^-1) sum_j e_j (x) f_j + ...) q^{-t} on V/Vbar (x) W^- (type
# I_-) agree exactly with the matrix parts at orders 0 and 1 (so the scalar series is 1 + O(x^2)); and a direct check of
# eq-triangle-braid for a_4 on complete charge sectors, both types, VV, VbarVbar, VVbar, VbarV.
def tri_block(h, qq, m, t, yy, typ, conj, Q):
    """Sector block of the type-I_-+ matrix part for V or Vbar (basis p (x) |Q - wt_p> -> |Q - wt_p'> (x) p')."""
    kap = qq - 1/qq
    if typ == '-':
        return np.array(closed_block(h, qq, -kap*m*yy, t, Q, 'Vbar' if conj else 'V'), complex)
    if conj:
        return np.array(closed_block(h, qq, kap*m/yy, t, Q, 'Vbar+'), complex)
    return np.linalg.inv(np.array(closed_block(h, qq, -kap*m/yy, tuple(1/x for x in t), tuple(-c for c in Q), 'rev'), complex))


rng = np.random.default_rng(3); fo = 0
for h in (3, 4, 5):
    qc = np.exp(1j*rng.uniform(.3, 1.2))*rng.uniform(.8, 1.2); kc = qc - 1/qc
    sc = rng.normal(size=h); sc -= sc.mean(); tc = tuple(qc**x for x in sc)
    mc, yc = complex(*rng.normal(size=2)), complex(*rng.normal(size=2))
    for typ in '-+':
        for conj in (False, True):
            sg = -1 if conj else 1; wt = [tuple(sg*(c == i) for c in range(h)) for i in range(h)]
            P = prep(h, qc, yc, conj); Q = tuple([2 + sg, -1] + [0]*(h - 2)); Q = Q[:-1] + (Q[-1] - sum(Q) + sg,)
            qn = lambda N, c: tc[c]*qp(qc, N[c])
            cart = lambda p, N: qn(N, p)**(-sg if typ == '-' else sg)      # q^{-+wt.(N+s)}: the Cartan factor q^{-+t}
            Z = np.zeros((h, h), complex); O = np.zeros((h, h), complex)
            for p_ in range(h):
                N = tuple(c - w for c, w in zip(Q, wt[p_])); Z[p_, p_] = cart(p_, N)
                for j in range(h):
                    al = alph(h, j); c2 = al.index(-1)
                    A_ = P[('e', j)] if typ == '-' else P[('f', j)]
                    for p2 in range(h):
                        if A_[p2][p_] != 0:
                            if typ == '-':
                                O[p2, p_] += -kc*A_[p2][p_]*mc*qn(N, c2)*cart(p_, N)
                            else:
                                O[p2, p_] += kc*A_[p2][p_]*mc*qn(N, c2)*cart(p2, tuple(n + a for n, a in zip(N, al)))
            ep = 1e-6
            T1 = (tri_block(h, qc, mc*ep, tc, yc, typ, conj, Q) - tri_block(h, qc, -mc*ep, tc, yc, typ, conj, Q))/(2*ep)
            T0_ = tri_block(h, qc, mc*1e-300, tc, yc, typ, conj, Q)
            fo = max(fo, np.abs(T0_ - Z).max()/np.abs(Z).max(), np.abs(T1 - O).max()/np.abs(O).max())


def tri_sector(h, qq, m, t, ya, yb, typ, ca, cb, Q, R):
    """|LHS - RHS|/|LHS| of eq-triangle-braid on the sector Q of V_a(ya) (x) V_b(yb) (x) W; R = R-check V_a (x) V_b."""
    wa = [tuple((-1 if ca else 1)*(c == i) for c in range(h)) for i in range(h)]
    wb = [tuple((-1 if cb else 1)*(c == i) for c in range(h)) for i in range(h)]
    sub = lambda Q_, w: tuple(x - y for x, y in zip(Q_, w))
    TA = {i: tri_block(h, qq, m, t, ya, typ, ca, sub(Q, wb[i])) for i in range(h)}    # V_a crossing, V_b = v_i aside
    TB = {i: tri_block(h, qq, m, t, yb, typ, cb, sub(Q, wa[i])) for i in range(h)}
    d = h*h; L = np.zeros((d, d), complex); Rm = np.zeros((d, d), complex)
    for i in range(h):
        for k in range(h):
            st = np.zeros((h, h), complex)                       # (1 x T_b), then (T_a x 1): W (x) V_a (x) V_b
            for k1 in range(h):
                st[:, k1] += TB[i][k1, k]*TA[k1][:, i]
            L[:, i*h + k] = R @ st.reshape(d)                    # then (1 x R)
            st1 = R[:, i*h + k].reshape(h, h)                    # (R x 1): V_b (x) V_a (x) W, (k', i')
            st2 = np.array([TA[k1] @ st1[k1] for k1 in range(h)])                    # (1 x T_a)
            Rm[:, i*h + k] = np.array([TB[i2] @ st2[:, i2] for i2 in range(h)]).T.reshape(d)   # (T_b x 1)
    return np.linalg.norm(L - Rm)/np.linalg.norm(L)


t1 = time.time(); h5 = 5; ya, yb = 0.7, 1.9*np.exp(0.2j); s5 = s_generic(h5); trs = []; ctr5 = []; nul5 = set()
qm = -q; mum = -np.exp(-1j*np.pi*LAM/2)/(qm - 1/qm)
for ca, cb in ((False, False), (True, True), (False, True), (True, False)):
    R5, n5 = rcheck(an_vector(h5 - 1, q, ya, conj=ca), an_vector(h5 - 1, q, yb, conj=cb)); nul5.add(n5)
    for typ in '-+':
        m5 = MU/1.3 if typ == '-' else MU*1.3; mm = mum/1.3 if typ == '-' else mum*1.3
        for Q in ((1, 0, 0, 0, 0), (2, -1, 1, 0, -2)):
            Q = Q[:-1] + (Q[-1] + (-1 if ca else 1) + (-1 if cb else 1) - sum(Q),)
            trs.append(tri_sector(h5, q, m5, tuple(q**x for x in s5), ya, yb, typ, ca, cb, Q, R5))
        ctr5.append(tri_sector(h5, qm, mm, tuple(qm**x for x in s5), ya, yb, typ, ca, cb, Q, R5))
report(f'prp-an-transmission (4) for every r: the terms of order 0 and 1 of R (type I_+) and R_21^-1 (type I_-) on '
       f'V, Vbar (x) W are exactly those of the matrix parts (h = 3, 4, 5, random complex parameters) [{fo:.0e}]; a_4 on '
       f'complete charge sectors: eq-triangle-braid holds for both types and VV, VbarVbar, VVbar, VbarV, generic s '
       f'[{max(trs):.0e}], R-check nullity {nul5}; control: matrices at -q with R at q [{min(ctr5):.2f}] '
       f'({time.time() - t1:.0f} s)', fo < 1e-9 and max(trs) < 1e-12 and nul5 == {1} and min(ctr5) > 0.1)
print(f'# (i) done after {time.time() - T0:.0f} s')


# ---------------------------------------------------------------- (j) cz2011 (3.12): poles and functional equations
def cz11_log(z, c, K=200):
    """log of the displayed rho-hat of cz2011 (3.12), Gamma functions grouped at equal k, each group divided by its value at
    z = 0: groups Phi_k(z) Phi_k(-z), Phi_k(u) = G(1/2+u+6kc+c) G(1/2+u+6kc+2c)/(G(1/2+u+6kc+4c) G(1/2+u+6kc+5c)), c = pi/(3 xi)."""
    k = np.arange(4*K)[:, None]
    def Phi(u):
        A = 0.5 + u + 6*k*c
        return loggamma(A + c) + loggamma(A + 2*c) - loggamma(A + 4*c) - loggamma(A + 5*c)
    z = np.atleast_1d(np.asarray(z, complex))[None, :]
    terms = Phi(z) + Phi(-z) - 2*Phi(0*z)
    return np.array([richardson_sum(terms[:, i], K) for i in range(terms.shape[1])])


fe = 0; ctl = 1e9
for xi in (0.9, 1.7):
    c = np.pi/(3*xi); zs = np.array([0.11 + 0.07j, -0.23 + 0.31j, 0.05 - 0.4j])
    # (3.10): rho(theta) = (1 + tau^2 x) rho(theta + i pi/3) rho(theta - i pi/3); theta +- i pi/3 <-> z -+ c, tau^2 x = e^{-2 pi i z}
    b10 = cz11_log(zs, c) + loggamma(0.5 + zs) + loggamma(0.5 - zs) - cz11_log(zs - c, c) - cz11_log(zs + c, c)
    # (3.8): rho(theta) rho(theta - i pi) (1 + tau^2 q^-4 x)(1 - tau^2 q^-2 x) = 1; theta - i pi <-> z + 3c
    b8 = cz11_log(zs, c) + cz11_log(zs + 3*c, c) - loggamma(0.5 + zs + c) - loggamma(0.5 + zs + 2*c) \
        - loggamma(0.5 - zs - c) - loggamma(0.5 - zs - 2*c)
    w10 = cz11_log(zs, c) + loggamma(0.5 + zs) + loggamma(0.5 - zs) - cz11_log(zs - 2*c, c) - cz11_log(zs + 2*c, c)
    fe = max(fe, np.abs(np.exp(b10 - b10[0]) - 1).max(), np.abs(np.exp(b8 - b8[0]) - 1).max())
    ctl = min(ctl, np.abs(np.exp(w10 - w10[0]) - 1).max())
xi = 0.9; nu = 0.3 + 0.5j*np.pi
zp = [1j*(nu - 1j*(xi*f_ + np.pi*m_/3) - nu)/xi for f_, m_ in ((0.5, 1), (0.5, 2), (1, 1), (1, 2))]
args = [0.5 - zp[i] + (i % 2 + 1)*np.pi/(3*xi) for i in range(4)]     # Gamma arguments: 0 at xi/2, -1/2 at the printed xi
xs = np.linspace(-3, 6, 18001)
sg = lambda ph: (np.cosh(0.4 + 1j*ph).real > 0) & (np.cosh(0.4 + 1j*ph).imag < 0)
w1 = xs[sg(np.pi/6 - xs/2)]; w2 = xs[sg(-np.pi/6 - xs/2)]


def cz11_printed(z, c, K):
    k = np.arange(1, K + 1)
    return (loggamma(0.5 - z + c) + loggamma(0.5 - z + 2*c)
            + np.sum(loggamma(0.5 + z + 6*k*c - 4*c) + loggamma(0.5 + z + 6*k*c - 5*c) - loggamma(0.5 + z + 6*k*c - 2*c)
                     - loggamma(0.5 + z + 6*k*c - c) + loggamma(0.5 - z + 6*k*c + 2*c) + loggamma(0.5 - z + 6*k*c + c)
                     - loggamma(0.5 - z + 6*k*c - 2*c) - loggamma(0.5 - z + 6*k*c - c)))


# control: the printed spectral variable of (3.2), x = q^{2 pi theta/xi}, tau = e^{-pi nu/xi}, theta = nu - i xi z, with
# log q = i 4 pi^2/beta^2 = i pi (c + 1/2) (q = e^{i 4 pi^2/beta^2}, c = pi/(3 xi) = 4 pi/beta^2 - 1/2); rho = e^{i pi z} rho-hat
pr = 0
for xi in (0.9, 1.7):
    c = np.pi/(3*xi); nu_r = 0.4; zs = np.array([0.11 + 0.07j, -0.23 + 0.31j, 0.05 - 0.4j])
    Lr = lambda z_: 1j*np.pi*z_ + cz11_log(z_, c)
    th = nu_r - 1j*xi*zs
    t2x_pr = np.exp(-2*np.pi*nu_r/xi + (2*np.pi*th/xi)*1j*np.pi*(c + 0.5))
    t2x_ok = np.exp(-2j*np.pi*zs)
    bpr = Lr(zs) - np.log(1 + t2x_pr) - Lr(zs - c) - Lr(zs + c)
    bok = Lr(zs) - np.log(1 + t2x_ok) - Lr(zs - c) - Lr(zs + c)
    fe = max(fe, np.abs(np.exp(bok - bok[0]) - 1).max())
    pr = max(pr, np.abs(np.exp(bpr - bpr[0]) - 1).max())
c = np.pi/(3*0.9); z = 0.3 + 0.5j
drift = (cz11_printed(z, c, 2000) - cz11_printed(z, c, 1000)) - (cz11_printed(0, c, 2000) - cz11_printed(0, c, 1000))
report(f'cz2011 (3.12): grouped at equal k, rho-hat solves (3.10) and (3.8) up to constants with z = i(theta - nu)/xi, tau^2 x = '
       f'e^(-2 pi i z), i.e. x = e^(2 pi theta/xi) in (3.2) (xi = 0.9, 1.7) [{fe:.0e}]; controls: doubled shifts [{ctl:.2f}], the printed '
       f'x = q^(2 pi theta/xi) of (3.2) in (3.10) [{pr:.2f}]; the arguments of the Gamma functions '
       f'outside the products vanish at nu - i(xi/2 + pi/3), nu - i(xi/2 + 2 pi/3) [arguments {abs(args[0]):.0e}, {abs(args[1]):.0e}] and equal '
       f'-1/2 at the printed nu - i(xi + pi/3), nu - i(xi + 2 pi/3); Re E > 0, Im E < 0 (eta > 0) for {w1.min()/np.pi:.3f} pi < xi < '
       f'{w1.max()/np.pi:.3f} pi and {w2.min()/np.pi:.3f} pi < xi < {w2.max()/np.pi:.3f} pi; as printed, normalized at z = 0, the '
       f'products drift like K^-2z [{abs(drift + 2*z*np.log(2)):.0e}]',
       fe < 1e-8 and ctl > 1e-2 and pr > 1e-2 and abs(args[0]) < 1e-12 and abs(args[1]) < 1e-12 and abs(args[2] + 0.5) < 1e-12
       and abs(args[3] + 0.5) < 1e-12 and abs(w1.min() - np.pi/3) < 1e-3 and abs(w1.max() - 4*np.pi/3) < 1e-3
       and abs(w2.min() + np.pi/3) < 1e-3 and abs(w2.max() - 2*np.pi/3) < 1e-3 and abs(drift + 2*z*np.log(2)) < 1e-3)
print(f'# total runtime {time.time() - T0:.0f} s')
