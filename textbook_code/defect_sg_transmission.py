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
Chapter 37 (sec-sg-transmission-matrix, sec-sg-defect-spectrum, sec-sg-breather-transmission):
 7. Scalar factor Phi (defect_common.sg_transmission_scalar = Konik-LeClair minimal solution, BCZ (3.43)-(3.44) with
    gamma = lambda): |Phi|^2 (1 + Y^2) = 1, Phi(theta) Phi(theta + i pi)(1 + Y^2) = 1, conj Phi = Phi(. + i pi);
    simple pole at theta_d - i pi/(2 lambda), none at theta_d + i pi/(2 lambda), none found in the physical strip;
    BCZ (3.42)/(3.47) as printed have a divergent phase (correction).                        [eq-sg-transmission-matrix, prp-sg-transmission]
 8. The full T (eq-sg-transmission-matrix) is unitary for real theta (also at -q: unitarity is blind to the sign);
    defect crossing eq-defect-crossing with Tbar = T(-theta)^-1 holds at q_phys with C = sigma_x, fails at -q_phys with
    sigma_x and is restored there only by C = [[0,1],[-1,0]] (the kappa = -1 conjugation).  [prp-sg-transmission]
 9. Resonance: rank-one residue in each charge sector at theta_d - i pi/(2 lambda), energy on the soliton mass shell
    with Im E < 0 (BCZ (4.1)-(4.7)); the bootstrap through the image of the residue gives the odd-sector matrix = s = 1
    Weyl intertwiner at theta_d - i pi (BCZ (4.15)), not unitary.                            [sec-sg-defect-spectrum]
10. Breathers by the bootstrap: B_1 factor = T_cl (coupling independent, label independent; BCZ (5.7)); B_2 = BCZ (5.10);
    product formula = BCZ (5.12)-(5.13) for n <= 4; unitarity, crossing, poles of T_B2 (BCZ (5.11)).
                                                                       [prp-sg-breather-transmission, eq-shg-linear-transmission]
11. Semiclassics: (pi/(lambda+1)) d arg Phi/d theta -> ln|coth((theta - theta_d)/2)| (Wigner delay of eq-sg-delay).
12. Half a soliton (CZ hep-th 1006.0939 sec. 4): S_{s B_1} by the bulk bootstrap (S_0 as a Gamma product, = eq-sg-s0) equals
    the product of two B_1 defect factors with theta_d = +-i(pi + xi/2).                        [prp-half-soliton]
13. Sinh-Gordon (Bajnok-Simon 0712.4292): Tbar K_D T = K_shG(E, F = 2i theta_d/pi) of eq-shg-reflection; their UV-IR
    relation (17) is theta_d = (1 - B/2) theta_d^cl, consistent with Corrigan-Taormina after classical fusion.
                                                                                               [cnj-defect-parameter]
Window: W truncated to n = -5..5 throughout; equations imposed on |n| <= 3 (|n| <= 2 for the odd-sector bootstrap).
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

# ======================================================================== ch. 37: unitarity, crossing, scalar factor, bootstrap
from defect_common import sg_transmission_scalar as Phi, sg_S0, sg_tmatrix, invariant_vector
kron = np.kron
TD = TH_D
C_plain = np.array([[0, 1], [1, 0]])            # s <-> sbar (the conjugation of the kappa = +1 bulk amplitude)
C_prime = np.array([[0, 1], [-1, 0]])           # the conjugation of the kappa = -1 amplitude at -q (sec-crossing-sign)


def Tcl(D):
    """eq-shg-linear-transmission (same closed form as in defect_classical_sg.py)."""
    return -1j*np.sinh(D/2 - 1j*np.pi/4)/np.sinh(D/2 + 1j*np.pi/4)


def blk(z, D):
    """BCZ hep-th/0506169 (5.9) = the book blocks (eq-blocks) with h = 2."""
    return np.sinh(D/2 + 1j*np.pi*z/4)/np.sinh(D/2 - 1j*np.pi*z/4)


# ---------------------------------------------------------------- 7. the scalar factor
devU, devC, devH = 0, 0, 0
for D in (-0.7, 0.3, 1.4):
    Y = np.exp(LAM*D)
    devU = max(devU, abs(abs(Phi(D, LAM))**2*(1 + Y**2) - 1))
    devH = max(devH, abs(np.conj(Phi(D, LAM)) - Phi(D + 1j*np.pi, LAM)))
for D in (-0.7, 0.3, 1.4, 0.2 + 0.3j, -0.4 - 0.6j):
    Y = np.exp(LAM*D)
    devC = max(devC, abs(Phi(D, LAM)*Phi(D + 1j*np.pi, LAM)*(1 + Y**2) - 1))
report(f'scalar factor (KL minimal solution, BCZ (3.43)-(3.44), gamma = lambda): |Phi|^2 (1 + Y^2) = 1 for real theta '
       f'[{devU:.0e}], Phi(theta) Phi(theta + i pi)(1 + Y^2) = 1 [{devC:.0e}], conj Phi(theta) = Phi(theta + i pi) [{devH:.0e}]',
       devU < 1e-9 and devC < 1e-8 and devH < 1e-8)
# poles: simple pole at theta_d - i pi/(2 lambda), none at theta_d + i pi/(2 lambda) (the factor 1/(1 + iY) is cancelled)
DP = -1j*np.pi/(2*LAM)
r1, r2 = 1e-5*Phi(DP + 1e-5, LAM), 1e-6*Phi(DP + 1e-6, LAM)
reg = abs(Phi(-DP + 1e-7, LAM))
grid = [abs(Phi(a + 1j*b, LAM)) for a in (-1.5, -0.4, 0, 0.6, 1.8) for b in np.linspace(0.05, np.pi/2 - 0.05, 7)]
report(f'Phi has a simple pole at theta = theta_d - i pi/(2 lambda) (residue {abs(r2):.3f}), is regular at theta_d + i pi/(2 lambda) '
       f'(|Phi| = {reg:.3f}), and is bounded on a grid of the physical strip 0 < Im theta < pi/2 (max {max(grid):.2f})',
       abs(r1 - r2) < 1e-4*abs(r2) and reg < 10 and max(grid) < 10)
# BCZ (3.42)/(3.47) as printed: the k-product has a divergent phase -(lambda D/pi) log K. The paired product r of (3.44)
# diverges on its own (like K^((lambda+1)/4)); what converges is the ratio r/rbar of (3.43) with a common cutoff, used here
from scipy.special import loggamma as LG
D = 0.6; w = LAM*D/np.pi
ref = np.log(Phi(D, LAM)*(1 + 1j*np.exp(LAM*D))) - 1j*np.pi*(1 + LAM)/4      # log r/rbar from (3.44)


def part347(K):
    k = np.arange(1, K + 1)
    return (LG(.25 - 1j*w/2) - LG(.25 + 1j*w/2) + np.sum(LG(.5 + 2*k*LAM - 1j*w) + LG(.5 + (2*k - 1)*LAM + 1j*w)
                                                       - LG(.5 + 2*k*LAM + 1j*w) - LG(.5 + (2*k - 1)*LAM - 1j*w)))


drift = [np.imag(part347(K) - ref) for K in (100, 400, 1600)]
slope = (drift[2] - drift[0])/np.log(16)
report(f'BCZ (3.47) truncated at K: phase of r/rbar drifts like -(lambda theta~/pi) log K (slope {slope:.4f} vs {-w:.4f}); '
       '(3.42)/(3.47) need an endpoint normalization; the ratio r/rbar of (3.43) with a common cutoff converges', abs(slope + w) < 2e-3)

# ---------------------------------------------------------------- 8. unitarity and crossing of the full T; the sign of q
NW = 5


def interior(W, m=2):
    return np.abs(W['n']) <= W['n'].max() - m


def unitarity_res(qq, s, th):
    T, V, W = sg_tmatrix(th, LAM, qq, TD, NW, s=s)
    dW = W['dim']
    inn = np.repeat(interior(W)[None, :], 2, axis=0).ravel()          # V (x) W ordering
    U = T.conj().T @ T
    return np.abs(U[np.ix_(inn, inn)] - np.eye(inn.sum())).max()


def crossing_res(qq, Cm, th):
    """T_{a alpha}^{beta b}(theta) - sum C_{d b} Tbar_{alpha d}^{c beta}(i pi - theta) Cinv_{a c},  Tbar(i pi - theta) =
    T(theta - i pi)^{-1} (eq-defect-unitarity, eq-defect-crossing); relative, on interior defect labels."""
    T, V, W = sg_tmatrix(th, LAM, qq, TD, NW)
    Tb = np.linalg.inv(sg_tmatrix(th - 1j*np.pi, LAM, qq, TD, NW)[0])        # W (x) V -> V (x) W
    Ci = np.linalg.inv(Cm); dW = W['dim']; ns = W['n']
    res, ref = [], []
    for k in np.where(interior(W))[0]:
        for a in range(2):
            for b in range(2):
                for k2 in range(dW):
                    lhs = T[k2*2 + b, a*dW + k]
                    rhs = sum(Cm[d, b]*Tb[c*dW + k2, k*2 + d]*Ci[a, c] for c in range(2) for d in range(2))
                    res.append(lhs - rhs); ref.append(lhs)
    return np.abs(res).max()/np.abs(ref).max()


ru = max(unitarity_res(q, s, th) for s in (0.0, 0.3) for th in (-0.4, 0.5, 1.2))
ru_m = max(unitarity_res(-q, 0.0, th) for th in (-0.4, 0.5))
report(f'T (with Phi) is unitary for real theta, theta_d (s = 0, 0.3) [{ru:.0e}]; so is its analogue at -q [{ru_m:.0e}]: '
       'unitarity does not see the sign of q', ru < 1e-9 and ru_m < 1e-9)
cq = max(crossing_res(q, C_plain, th) for th in (0.4, -0.3))
cmq = min(crossing_res(-q, C_plain, th) for th in (0.4, -0.3))
cmq2 = max(crossing_res(-q, C_prime, th) for th in (0.4, -0.3))
report(f'defect crossing (eq-defect-crossing) with Tbar = T(-theta)^-1: holds at q_phys with C = sigma_x [{cq:.0e}]; '
       f'fails at -q_phys with sigma_x [{cmq:.2f}: off-diagonal entries flip sign, no scalar factor can fix it]; '
       f'restored at -q only with C = [[0,1],[-1,0]] [{cmq2:.0e}]', cq < 1e-8 and cmq > 0.1 and cmq2 < 1e-8)

# ---------------------------------------------------------------- 9. the resonance pole and the odd sector
thP = TD + DP
M, V, W = sg_tmatrix(thP, LAM, q, TD, NW, scalar=False)
dW, ns = W['dim'], W['n']
iVW = lambda a, k: a*dW + k
Eimg = {}
rk = 0
for k, m in enumerate(ns):
    e = np.zeros(2*dW, complex); e[iVW(0, k)] = 1
    Eimg[m] = M @ e
    if k + 1 < dW:
        f = np.zeros(2*dW, complex); f[iVW(1, k + 1)] = 1
        g = M @ f
        rk = max(rk, np.linalg.norm(g - np.vdot(Eimg[m], g)/np.vdot(Eimg[m], Eimg[m])*Eimg[m])/np.linalg.norm(g))
Ms = 1.0
dE = Ms*np.cosh(thP); dP = Ms*np.sinh(thP)
report(f'at theta_P = theta_d - i pi/(2 lambda) the matrix part of T has rank one in each charge sector [{rk:.0e}] (factorizing '
       f'residue, BCZ (4.6)-(4.7)); Delta E = M cosh theta_P is on the soliton mass shell, Im Delta E = {dE.imag:.3f} M < 0 for theta_d > 0',
       rk < 1e-12 and abs(dE**2 - dP**2 - Ms**2) < 1e-12 and dE.imag < 0)
# Tbar(theta) = T(-theta)^-1 has its pole where Phi(D)(1 + Y^2) vanishes: a simple zero at D = -theta - theta_d = i pi/(2 lambda)
zr = [abs(Phi(-DP + e, LAM)*(1 + np.exp(2*LAM*(-DP + e)))) for e in (1e-5, 1e-6)]
Tbm = np.linalg.inv(sg_tmatrix(-(-TD + DP + 1e-6), LAM, q, TD, NW)[0])         # Tbar near theta = -theta_d - i pi/(2 lambda)
Tbm2 = np.linalg.inv(sg_tmatrix(-(-TD + DP + 1e-7), LAM, q, TD, NW)[0])
grow = np.abs(Tbm2).max()/np.abs(Tbm).max()
Eb = np.cosh(-TD + DP)                                                        # energy of a particle of rapidity -theta
report(f'Tbar(theta) = T(-theta)^-1 has a simple pole at theta = -theta_d - i pi/(2 lambda) (Phi (1 + Y^2) has a simple zero: '
       f'{zr[0]/zr[1]:.2f}; |Tbar| grows x{grow:.1f} per decade); a particle from the right of rapidity -theta there has '
       f'Im E = {Eb.imag:.3f} M, < 0 iff theta_d < 0', abs(zr[0]/zr[1] - 10) < 1e-3 and abs(grow - 10) < 1e-2
       and np.sign(Eb.imag) == np.sign(TD))


def excited_T(th):
    """Bootstrap (eq-defect-bootstrap, non-diagonal form): the excited defect E = image of the residue in W (x) V(theta_P);
    a soliton crosses the defect, then the absorbed soliton: R_23(theta - theta_P) T_12(theta) restricted to V (x) E."""
    T = sg_tmatrix(th, LAM, q, TD, NW, scalar=False)[0]
    X = kron(np.eye(dW), Rsg(np.exp(LAM*(th - thP)))) @ kron(T, np.eye(2))
    labels = [(m2, b) for m2 in ns for b in range(2)]
    Bm = np.array([kron(Eimg[m2], np.eye(2)[b]) for m2, b in labels]).T
    el, fit = {}, 0
    for a in range(2):
        for m in ns[np.abs(ns) <= NW - 3]:
            vout = X @ kron(np.eye(2)[a], Eimg[m])
            coef = np.linalg.lstsq(Bm, vout, rcond=None)[0]
            fit = max(fit, np.linalg.norm(Bm @ coef - vout)/np.linalg.norm(vout))
            for (m2, b), c in zip(labels, coef):
                if abs(c) > 1e-10:
                    el[(a, m, b, m2)] = c
    return el, fit


ok_odd = True; fits = []
for th in (0.6, -0.2):
    el, fit = excited_T(th)
    fits.append(fit)
    y = np.exp(LAM*th)
    A_ = {m: el[(0, m, 0, m)] for m in (-1, 0, 1)}; D_ = {m: el[(1, m, 1, m)] for m in (-1, 0, 1)}
    C_ = el[(0, 0, 1, 1)]; B_ = el[(1, 1, 0, 0)]
    shape = all(k[0] == k[2] and k[1] == k[3] or k[0] == 0 and k[2] == 1 and k[3] == k[1] + 1
                or k[0] == 1 and k[2] == 0 and k[3] == k[1] - 1 for k in el)
    ok_odd &= shape and np.isclose(A_[1]/A_[0], 1/q) and np.isclose(D_[0]/A_[0], q) and np.isclose(D_[1]/A_[1], q**3)
    ok_odd &= np.isclose(C_*B_/(A_[0]*D_[1]), -(y/YD)**2/q**2)
    # non-unitary: the two columns of a charge sector are not orthogonal for real theta
    ok_odd &= abs(np.conj(A_[0])*B_ + np.conj(C_)*D_[1]) > 1e-3*abs(A_[0]*D_[1])
report(f'odd sector by the bootstrap at theta_P [fit {max(fits):.0e}]: the excited defect transmits as the s = 1 Weyl intertwiner '
       'with theta_d -> theta_d - i pi (A_{m+1}/A_m = 1/q, D_m/A_m = q^{2m+1}, CB/AD = -(y/y_d)^2/q^2; BCZ (4.15): nu^ = nu, '
       'p^2 -> q_BCZ^2 p^2), and is not unitary for real theta', ok_odd)

# ---------------------------------------------------------------- 10. breathers by the bootstrap
def breather_T(nb, th, lam=LAM, n0=0):
    qq = q_phys(lam); u = np.pi - nb*np.pi/lam
    t1, t2 = th - 1j*u/2, th + 1j*u/2
    rel, w = invariant_vector(np.exp(lam*t1), np.exp(lam*t2), qq)
    T1, V, W = sg_tmatrix(t1, lam, qq, TD, NW); T2 = sg_tmatrix(t2, lam, qq, TD, NW)[0]
    L = kron(T1, np.eye(2)) @ kron(np.eye(2), T2)                  # V x V x W -> W x V x V
    e = np.zeros(W['dim']); e[list(W['n']).index(n0)] = 1
    vout = L @ kron(w, e); ref = kron(e, w)
    c = np.vdot(ref, vout)/np.vdot(ref, ref)
    return c, max(rel, np.linalg.norm(vout - c*ref)/np.linalg.norm(vout))


dev1, devfit = 0, 0
for lam_ in (LAM, 2.6):
    for th in (0.2, 0.9, -0.5):
        for n0 in (-1, 0, 1):
            c, fit = breather_T(1, th, lam_, n0)
            dev1 = max(dev1, abs(c - Tcl(th - TD))); devfit = max(devfit, fit)
report(f'lightest breather: the soliton bootstrap through T gives a scalar, independent of the defect label, equal to the '
       f'classical T_cl(theta - theta_d) at lambda = 1.37 and 2.6 (coupling independent; BCZ (5.7)) [{dev1:.0e}, fit {devfit:.0e}]',
       dev1 < 1e-8 and devfit < 1e-10)
dev2 = 0
for th in (0.2, 0.9):
    c, fit = breather_T(2, th)
    D = th - TD
    dev2 = max(dev2, abs(c + 1/(blk(1 + 1/LAM, D)*blk(1 - 1/LAM, D))),
               abs(c - Tcl(D + 1j*np.pi/(2*LAM))*Tcl(D - 1j*np.pi/(2*LAM))), fit)
report(f'second breather from the solitons = BCZ (5.10) -1/((1+1/lambda)(1-1/lambda)) = T_B1(theta + i pi/2lambda) '
       f'T_B1(theta - i pi/2lambda) [{dev2:.0e}]', dev2 < 1e-8)


def TBn_product(nb, D, lam):
    return np.prod([Tcl(D + 1j*np.pi*(nb + 1 - 2*k)/(2*lam)) for k in range(1, nb + 1)])


def TBn_bcz(nb, D, lam):
    if nb % 2:
        s = (nb - 1)//2
        return (-1)**(s + 1)*1j/blk(1, D)*np.prod([1/(blk(1 + 2*(s - l)/lam, D)*blk(1 - 2*(s - l)/lam, D)) for l in range(s)])
    s = nb//2
    return (-1)**s*np.prod([1/(blk(1 + (2*s - 2*l - 1)/lam, D)*blk(1 - (2*s - 2*l - 1)/lam, D)) for l in range(s)])


dev3 = max(abs(TBn_product(nb, D, 4.7) - TBn_bcz(nb, D, 4.7)) for nb in (1, 2, 3, 4) for D in (0.3, -1.1, 0.2 + 0.1j))
devuc = max(abs(TBn_product(nb, D, 4.7)*TBn_product(nb, D + 1j*np.pi, 4.7) - 1) for nb in (1, 2, 3) for D in (0.3, -1.1))
devu = max(abs(abs(TBn_product(nb, D, 4.7)) - 1) for nb in (1, 2, 3) for D in (0.3, -1.1))
# poles of T_B2 at theta~ = -i pi (1 +- 1/lambda)/2 (BCZ (5.11))
p2 = [abs(1e-7*TBn_product(2, -1j*np.pi*(1 + sg/LAM)/2 + 1e-7, LAM)) for sg in (1, -1)]
report(f'B_n = prod_k T_B1(theta + i pi (n+1-2k)/(2 lambda)) equals BCZ (5.12)-(5.13) for n = 1..4 [{dev3:.0e}]; unitary and '
       f'T_Bn(theta) T_Bn(theta + i pi) = 1 [{max(devu, devuc):.0e}]; T_B2 has poles at theta~ = -i pi(1 +- 1/lambda)/2 '
       f'(residues {p2[0]:.2f}, {p2[1]:.2f})', dev3 < 1e-12 and max(devu, devuc) < 1e-12 and min(p2) > 1e-3)

# ---------------------------------------------------------------- 11. semiclassical limit of the phase: Wigner time delay
rows = []
for lam_ in (10.0, 40.0, 160.0):
    h = 1e-5
    rows.append([(np.pi/(lam_ + 1))*np.angle(Phi(D + h, lam_, K=400)/Phi(D - h, lam_, K=400))/(2*h)
                 - np.log(abs(1/np.tanh(D/2))) for D in (-0.7, 0.5, 1.1)])
rows = np.abs(np.array(rows))
report(f'semiclassics: (pi/(lambda+1)) d arg Phi/d theta - ln|coth((theta-theta_d)/2)| = {rows.max(axis=1)[0]:.3f}, '
       f'{rows.max(axis=1)[1]:.4f}, {rows.max(axis=1)[2]:.4f} at lambda = 10, 40, 160 (-> 0 like 1/lambda): '
       'Wigner delay p dt = (8/beta^2) ln|z| of eq-sg-delay', np.all(rows[1] < rows[0]/3) and np.all(rows[2] < rows[1]/3))

# ---------------------------------------------------------------- 12. half a soliton, quantum version (CZ 2010 sec. 4)
import mpmath as mp
mp.mp.dps = 15


def S0_int(theta):
    L_ = mp.mpf(LAM); xs = mp.pi/L_
    f = lambda t: mp.sinh((mp.pi - xs)*t/2)*mp.sin(theta*t)/(t*mp.sinh(xs*t/2)*mp.cosh(mp.pi*t/2))
    return complex(-mp.exp(-1j*mp.quad(f, [0, 1, 5, 20, mp.inf])))


devS0 = max(abs(sg_S0(th, LAM) - S0_int(th)) for th in (0.23, 0.61))
xi = np.pi/LAM; u1 = np.pi - xi


def S_sB1(Th):
    t1, t2 = -1j*u1/2, 1j*u1/2
    _, w = invariant_vector(np.exp(LAM*t1), np.exp(LAM*t2), q)
    y = np.exp(LAM*Th)
    X = kron(np.eye(2), Rsg(y/np.exp(LAM*t2))) @ kron(Rsg(y/np.exp(LAM*t1)), np.eye(2))
    sc = sg_S0(Th - t1, LAM)*sg_S0(Th - t2, LAM)
    vin = kron(np.eye(2)[0], w); vout = sc*(X @ vin); ref = kron(w, np.eye(2)[0])
    return np.vdot(ref, vout)/np.vdot(ref, ref)


devh = 0
for Th in (0.4, -0.9):
    exact = (np.sinh(Th) + 1j*np.cos(xi/2))/(np.sinh(Th) - 1j*np.cos(xi/2))
    pair = Tcl(Th - 1j*(np.pi + xi/2))*Tcl(Th + 1j*(np.pi + xi/2))
    devh = max(devh, abs(S_sB1(Th) - exact), abs(pair - exact))
report(f'S_0 as a Gamma product = eq-sg-s0 [{devS0:.0e}]; soliton-B_1 amplitude by the bulk bootstrap = (sinh + i cos(xi/2))/(sinh - '
       f'i cos(xi/2)) = T_B1(theta_d = i(pi + xi/2)) T_B1(theta_d = -i(pi + xi/2)) [{devh:.0e}]', devS0 < 1e-7 and devh < 1e-7)

# ---------------------------------------------------------------- 13. sinh-Gordon: Bajnok-Simon fusion onto a Dirichlet boundary
Bb = 0.37                                        # book sinh-Gordon B (eq-sinh-gordon-smatrix); B_BS = B/2
h2 = lambda x, th: np.sinh(th/2 + 1j*np.pi*x/4)/np.sinh(th/2 - 1j*np.pi*x/4)


def K_shG(th, E, F):
    """eq-shg-reflection."""
    return h2(1, th)*h2(2 - Bb/2, th)*h2(1 + Bb/2, th)/(h2(1 - E, th)*h2(1 + E, th)*h2(1 - F, th)*h2(1 + F, th))


def K_D(th, E):
    return h2(1, th)*h2(2 - Bb/2, th)*h2(1 + Bb/2, th)/(h2(1 - E, th)*h2(1 + E, th))


devK, devTT = 0, 0
for td in (0.45, -0.8):
    for E in (0.3, 0.7j):
        for th in (0.35, 1.2, 0.2 + 0.4j):
            T_ = Tcl(th - td); Tb_ = Tcl(1j*np.pi - th - td)              # Tbar(theta) = T(i pi - theta), depends on theta + theta_d
            devK = max(devK, abs(Tb_*K_D(th, E)*T_ - K_shG(th, E, 2j*td/np.pi)))
            devTT = max(devTT, abs(Tcl(th - td)*Tcl(th + 1j*np.pi - td) - 1))       # Tbar(-theta) T(theta) = 1
# UV-IR: Bajnok-Simon (17) mu = b^2 vartheta/8pi with theta_d = -B_BS vartheta (their (16)) and theta_d^cl = -mu (sigma = e^mu)
b2 = 8*np.pi*Bb/(2 - Bb)                         # book B = 2 b^2/(8 pi + b^2), b = beta_SG of sinh-Gordon
vt = 0.83; mu = b2*vt/(8*np.pi)
ratio = (-(Bb/2)*vt)/(-mu)
# Corrigan-Taormina E, F = (a0 +- a1)(1 - B/2), C_j = cos(pi a_j), with the classical fusion C_{0,1} = cosh(theta_d^cl -+ b phi_D/2)
tcl, phD = -mu, 0.4
a0 = 1j*(tcl - np.sqrt(b2)*phD/2)/np.pi; a1 = 1j*(tcl + np.sqrt(b2)*phD/2)/np.pi
okCT = np.isclose(np.cos(np.pi*a0), np.cosh(tcl - np.sqrt(b2)*phD/2)) and np.isclose((a0 + a1)*(1 - Bb/2), 2j*(ratio*tcl)/np.pi)
report(f'Bajnok-Simon: Tbar(theta) K_D(theta; E) T(theta) = K_shG(theta; E, F = 2i theta_d/pi) (eq-shg-reflection) [{devK:.0e}]; '
       f'T T(+i pi) = 1 [{devTT:.0e}]; their (17) is theta_d = (1 - B/2) theta_d^cl [ratio {ratio:.4f} vs {1 - Bb/2:.4f}], '
       'the same factor as in Corrigan-Taormina E = (a0+a1)(1-B/2) after classical fusion (E, F enter K_shG symmetrically)',
       devK < 1e-12 and devTT < 1e-12 and np.isclose(ratio, 1 - Bb/2) and okCT)
