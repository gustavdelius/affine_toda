"""Sine-Gordon soliton transmission through a jump defect: quantum conventions (Part VIII, ch. 37, 40).

Part VI conventions (boundary_common.py): Delta(e) = e (x) 1 + k (x) e, Delta(f) = f (x) k^-1 + 1 (x) f,
q = q_phys = -exp(-i pi lambda), lambda = 8 pi/beta_SG^2 - 1 = omega (sec-qg-gradation), y = exp(lambda theta);
spin-1/2 basis (v_+ soliton, v_- antisoliton).  Defect space W, Weyl modules and the matrix elements A, B, C, D of
T: V(y) (x) W -> W (x) V(y) as in defect_common.py.  BCZ = Bowcock-Corrigan-Zambon hep-th/0506169, KL = Konik-LeClair.

Checks (book ids they support):
 1. sec-qg-gradation: for sine-Gordon 2T = 2(4 pi/beta^2) - 2 gives T = lambda = omega; q_BCZ = e^{i pi gamma} with
    gamma = lambda equals -1/q_phys; the book's R-check(y1/y2) equals BCZ's S/a ((3.2)-(3.4)) entry by entry.   [tbl-defect-dictionary]
 2. Weyl modules of U_q(b^+) and U_q(b^-): Cartan relations and q-Serre on the window interior; e_0 e_1 = e_1 e_0
    (f_0 f_1 = f_1 f_0) central.                                                                          [eq-weyl-module]
 3. Borel intertwiners T: unique up to scalar on the window (s = 0 and s = 0.3), closed forms
      b^-:  A_n = q^-n,  D_n = q^(n+s),  C_n = -i (y/y_d) q^-n,      B_n = -i (y/y_d) q^(n+s)
      b^+:  A_n = q^n,   D_n = q^-(n+s), C_n = i q^-(1+s) (y_d/y) q^-n, B_n = i q^-1 (y_d/y) q^n
    (mu = i/(q - q^-1)); level zero is forced: with k_0 k_1 = kappa != 1 there is no intertwiner.   [eq-borel-intertwiner, prp-sg-t-unique]
 4. Both intertwiners solve the triangle relation with the book's R-check (braid form, defect_common).       [thm-t-triangle]
 5. BCZ (3.51) (gamma = lambda, eta = theta_d, Q^2 = -q_BCZ, alpha = 2n, any nu) solves the book's triangle relation and
    equals the b^- intertwiner up to a diagonal gauge on W: gauge invariants A_{n+1}/A_n = 1/q, D_n/A_n = nu q^{2n}
    (nu = q^s), C_n B_{n+1}/(A_n D_{n+1}) = -(y/y_d)^2.  So BCZ/KL's transmission matrix is the U_q(b^-) intertwiner,
    and the character-changing amplitude grows with theta - theta_d, as for the classical type I_- defect (z < 0 for
    theta > theta_d, defect_classical_sg.py); b^+ is the time-reversed type I_+.                       [prp-two-borels, eq-sg-transmission-matrix]
 6. Sign of q: Borel intertwiners built at q' = -q_phys, 1/q_phys, -1/q_phys do NOT satisfy the triangle relation with the
    physical R-check; each satisfies it with R-check(q') (the triangle relations alone do not fix the sign of q).
Runtime: about 100 seconds.
"""
import numpy as np
from boundary_common import report, q_phys, sl2_spin, rcheck, prop_residual, cop, check_relations, cartan_a
from defect_common import (weyl_module_b, transmission_window, tmatrix_elements, tmatrix_from_elements,
                           triangle_residual_braid, borel_kinds)

LAM = 1.37
q = q_phys(LAM)
N = 5
TH_D = 0.31
YD = np.exp(LAM*TH_D)
MU = lambda qq: 1j/(qq - 1/qq)


def Rsg(z, qq=q):
    return rcheck(sl2_spin(1, qq, z), sl2_spin(1, qq, 1.0))[0]


# ---------------------------------------------------------------- 1. dictionary
beta2 = 4*np.pi/(LAM + 1)                         # book beta^2 with lambda = 4 pi/beta^2 - 1
T_grad = (2*4*np.pi/beta2 - 2)/2                  # 2T = h(4 pi/beta^2) - h^vee, h = h^vee = 2
q_bcz = np.exp(1j*np.pi*LAM)
ok = np.isclose(T_grad, LAM) and np.isclose(q_bcz, -1/q)
th1, th2 = 0.83, -0.21
x1, x2 = np.exp(LAM*th1), np.exp(LAM*th2)
a = q_bcz*x2/x1 - x1/(q_bcz*x2); b = x1/x2 - x2/x1; c = q_bcz - 1/q_bcz
S = np.array([[a, 0, 0, 0], [0, c, b, 0], [0, b, c, 0], [0, 0, 0, a]])/a     # rows (k l), columns (m n) = S_kl^mn
R = Rsg(x1/x2)                                    # R[(m n), (k l)]: V(y1) (x) V(y2) -> V(y2) (x) V(y1)
ok &= np.allclose(R, S.T, atol=1e-12)
report('sine-Gordon: T = lambda = omega, y = e^{lambda theta}; q_BCZ = e^{i pi gamma} = -1/q_phys (gamma = lambda); '
       'R-check(y1/y2) = BCZ S/a entry by entry (transmission b/a, reflection c/a)', ok)

# ---------------------------------------------------------------- 2. modules
ok = True
for sign in (1, -1):
    W = weyl_module_b(sign, q, YD, N, s=0.3, mu=MU(q))
    I = W['interior']; blk = np.ix_(I, I)
    kinds = borel_kinds(sign)
    X = {j: W[(kinds[0][0], j)] for j in (0, 1)}
    k = {j: W[('k', j)] for j in (0, 1)}
    C = cartan_a(1)
    for i in (0, 1):
        for j in (0, 1):
            ok &= np.allclose((k[i] @ X[j] @ np.linalg.inv(k[i]))[blk], (q**(sign*C[i][j])*X[j])[blk])
    ok &= np.allclose((X[0] @ X[1] - X[1] @ X[0])[2:-2, 2:-2], 0)
    for i, j in ((0, 1), (1, 0)):
        q3 = q**2 + 1 + q**-2
        serre = (np.linalg.matrix_power(X[i], 3) @ X[j] - q3*np.linalg.matrix_power(X[i], 2) @ X[j] @ X[i]
                 + q3*X[i] @ X[j] @ np.linalg.matrix_power(X[i], 2) - X[j] @ np.linalg.matrix_power(X[i], 3))
        ok &= np.allclose(serre[4:-4, 4:-4], 0)
report('Weyl modules of U_q(b^+-): k_i x_j k_i^-1 = q^(+-a_ij) x_j, q-Serre, x_0 x_1 = x_1 x_0 (window interior)', ok)

# ---------------------------------------------------------------- 3. intertwiners and closed forms
closed = {
    -1: lambda y, s, yd: (lambda n: q**-n, lambda n: -1j*(y/yd)*q**(n + s), lambda n: -1j*(y/yd)*q**-n, lambda n: q**(n + s)),
    1: lambda y, s, yd: (lambda n: q**n, lambda n: 1j*q**-1*(yd/y)*q**n, lambda n: 1j*q**(-1 - s)*(yd/y)*q**-n, lambda n: q**-(n + s)),
}
ok_rank = True; ok_closed = True
for sign in (1, -1):
    for s in (0.0, 0.3):
        for th in (-0.4, 0.5, 1.2):
            y = np.exp(LAM*th)
            V = sl2_spin(1, q, y); W = weyl_module_b(sign, q, YD, N, s=s, mu=MU(q))
            T, rank = transmission_window(V, W)
            ok_rank &= rank == 1
            Tc = tmatrix_from_elements(*closed[sign](y, s, YD), V, W)
            mask = T != 0
            ok_closed &= prop_residual(T[mask], Tc[mask])[0] < 1e-9
report('Borel intertwiner V(y) (x) W -> W (x) V(y): unique up to scalar on the window (b^+, b^-; s = 0, 0.3; three theta)', ok_rank)
report('closed forms: b^- A = q^-n, D = q^(n+s), C = -i(y/y_d)q^-n, B = -i(y/y_d)q^(n+s); '
       'b^+ A = q^n, D = q^-(n+s), C = i q^-(1+s)(y_d/y)q^-n, B = i q^-1 (y_d/y) q^n', ok_closed)
# level zero is forced
V = sl2_spin(1, q, np.exp(LAM*0.5))
ranks = [transmission_window(V, weyl_module_b(sign, q, YD, N, s=0.0, mu=MU(q), kappa=1.7))[1] for sign in (1, -1)]
report(f'k_0 k_1 = 1.7 on W (not level zero): no intertwiner at all (interior ranks {ranks})', ranks == [0, 0])

# ---------------------------------------------------------------- 4. triangle relation
def T_borel(sign, s, qq=q, yd=YD):
    def f(y):
        V = sl2_spin(1, qq, y); W = weyl_module_b(sign, qq, yd, N, s=s, mu=MU(qq))
        return transmission_window(V, W)[0]
    return f


Wref = weyl_module_b(1, q, YD, N)
res_tri = []
for sign in (1, -1):
    for s in (0.0, 0.3):
        res_tri.append(triangle_residual_braid(Rsg, T_borel(sign, s), np.exp(LAM*0.9), np.exp(LAM*0.2), Wref)[0])
report(f'triangle relation (1 (x) R)(T(y1) (x) 1)(1 (x) T(y2)) ~ (T(y2) (x) 1)(1 (x) T(y1))(R (x) 1) for both Borels  '
       f'[max res {max(res_tri):.1e}]', max(res_tri) < 1e-9)

# ---------------------------------------------------------------- 5. BCZ (3.51)
nu = np.exp(0.4j)
Q2 = -q_bcz


def T_bcz(y):
    xt = y/YD                                     # e^{gamma (theta - eta)}
    V = sl2_spin(1, q, y)
    return tmatrix_from_elements(lambda n: nu**-0.5*Q2**n, lambda n: q_bcz**-0.5*xt, lambda n: q_bcz**-0.5*xt,
                                 lambda n: nu**0.5*Q2**-n, V, Wref)


r_bcz = triangle_residual_braid(Rsg, T_bcz, np.exp(LAM*0.9), np.exp(LAM*0.2), Wref)[0]
# gauge invariants vs the b^- intertwiner with s: q^s = nu
s_nu = np.log(nu)/np.log(q)                       # complex offset with q^s = nu
ok_inv = True; ok_dir = True
for th in (-0.3, 0.6):
    y = np.exp(LAM*th); V = sl2_spin(1, q, y)
    A1, B1, C1, D1 = tmatrix_elements(T_bcz(y), V, Wref)
    Wm = weyl_module_b(-1, q, YD, N, s=s_nu, mu=MU(q))
    A2, B2, C2, D2 = tmatrix_elements(transmission_window(V, Wm)[0], V, Wm)
    for n in (-2, 0, 1):
        for (A, B, C, D) in ((A1, B1, C1, D1), (A2, B2, C2, D2)):
            ok_inv &= np.isclose(A[n + 1]/A[n], 1/q) and np.isclose(D[n]/A[n], nu*q**(2*n))
            ok_inv &= np.isclose(C[n]*B[n + 1]/(A[n]*D[n + 1]), -(y/YD)**2)
report(f'BCZ (3.51) (gamma = lambda, eta = theta_d, Q^2 = -q_BCZ = 1/q_phys) solves the book triangle relation [res {r_bcz:.1e}]',
       r_bcz < 1e-9)
report('BCZ (3.51) = U_q(b^-) intertwiner with q^s = nu, y_d = e^{lambda theta_d}, mu^2 = -1/(q - q^-1)^2: '
       'A_{n+1}/A_n = 1/q, D_n/A_n = nu q^{2n}, C_n B_{n+1}/(A_n D_{n+1}) = -(y/y_d)^2', ok_inv)
# direction: character change dominates for theta > theta_d for b^- (type I_-, BCZ), for theta < theta_d for b^+ (type I_+)
ratio = {}
for sign in (1, -1):
    vals = []
    for th in (TH_D - 0.3, TH_D + 0.3):
        y = np.exp(LAM*th); V = sl2_spin(1, q, y); W = weyl_module_b(sign, q, YD, N, mu=MU(q))
        A_, B_, C_, D_ = tmatrix_elements(transmission_window(V, W)[0], V, W)
        vals.append(abs(C_[0]/A_[0]))
    ratio[sign] = vals
report('|C/A| (soliton -> antisoliton) grows through theta = theta_d for b^- (as classical type I_-, z_{I-} < 0 for theta > theta_d) '
       'and falls for b^+ (type I_+)', ratio[-1][1] > 1 > ratio[-1][0] and ratio[1][0] > 1 > ratio[1][1])

# ---------------------------------------------------------------- 6. the sign of q
out = []
for qq, name in ((-q, '-q_phys'), (1/q, '1/q_phys'), (-1/q, '-1/q_phys')):
    for sign in (1, -1):
        Wq = weyl_module_b(sign, qq, YD, N, mu=MU(qq))
        r_phys = triangle_residual_braid(Rsg, T_borel(sign, 0.0, qq), np.exp(LAM*0.9), np.exp(LAM*0.2), Wq)[0]
        r_own = triangle_residual_braid(lambda z: Rsg(z, qq), T_borel(sign, 0.0, qq), np.exp(LAM*0.9), np.exp(LAM*0.2), Wq)[0]
        out.append((name, sign, r_phys, r_own))
ok = all(rp > 1e-3 and ro < 1e-9 for _, _, rp, ro in out)
report('Borel intertwiners at q\' = -q_phys, 1/q_phys, -1/q_phys fail the triangle relation with the physical R-check '
       f'(min res {min(o[2] for o in out):.1e}) and pass with R-check(q\') (max res {max(o[3] for o in out):.1e})', ok)
