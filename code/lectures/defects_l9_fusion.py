"""Checks for Lecture 9 of the defect lecture notes, "Folding, fusing and the classical breather"
(lecture-notes/defects/09-folding-fusing.qmd).

Conventions of the lecture. Classical sine-Gordon theory as in Lectures 2-8: d_pm = d_t pm d_x, d_+ d_- Phi =
-m0^2 sin Phi, B = -2 m0 [sigma (cos(Sigma/2) - 1) + (cos(Delta/2) - 1)/sigma], sigma = e^{-theta_d}, Sigma = Phi + Psi,
Delta = Phi - Psi, defect conditions Phi_x = Psi_t - B_Phi, Psi_x = Phi_t + B_Psi at x = 0; z(u) = -coth(u/2),
u = theta - theta_d. Quantum: lambda = 8 pi/beta^2 - 1, q = -e^{-i pi lambda}, y = e^{lambda theta},
y_d = e^{lambda theta_d}, Y = y/y_d, the transmission matrix T of Lecture 7 (eq-l7-T) with the scalar factor Phi of
eq-l7-Phi, and the bulk Rcheck of Lecture 7. Boundaries: K(theta): V(theta) -> V(-theta) with the defect space W on the
right; the Ghoshal-Zamolodchikov matrix in the book's parametrization (eta, vartheta); Lee-Yang blocks
(x) = sinh(theta/2 + i pi x/12)/sinh(theta/2 - i pi x/12). Written from these definitions, independently of
code/defect_doubling.py, defect_axioms.py and defect_sg_transmission.py. One PASS or FAIL line per check; exit status 1
on any FAIL.

  (a) Folding the line. The doubled boundary K = [[0, Tbar], [T, 0]] on V + V' satisfies boundary unitarity exactly when
      T(theta) Tbar(-theta) = 1, and Ghoshal-Zamolodchikov crossing-unitarity k(theta) = Scheck(2 theta) k(-theta),
      k = sum K_{b abar}(i pi/2 - theta) v_a v_b, exactly when T(theta) = Tbar(i pi - theta): for the Lee-Yang defect and
      the classical factor both hold, for T e^{0.2 theta} unitarity holds and both forms of crossing fail together. For
      the sine-Gordon soliton matrix with its scalar factor, the doubled crossing-unitarity, component by component with
      the defect labels unchanged, holds at q = -e^{-i pi lambda} and fails at -q. The operator of Bajnok and George
      multiplies a particle state by Tbar(i pi/2 - theta) = T(theta + i pi/2), which for the classical factor is the
      factor of the defect turned as in Lecture 8 (theta_d -> theta_d - i pi/2).
  (b) The reflected copy and the other Borel half. P Rcheck P intertwines the opposite coproduct and not the coproduct.
      On a window, the y-independent operators Z, kappa with L (e1 (x) 1 + k1 (x) Z) = (1 (x) Z + e0 (x) kappa) L (L = P T,
      the generator e1 (x) 1 + k1 (x) e0 of (1 (x) tau) Delta(U_q(b+))) exist, with kappa |n> = q^{-2n-2-s} |n+2>
      = C B^-1 and Z |n> = a q^{-2n} |n+1>, a = i y_d q^{-s}/(q^2 - 1), up to a constant fixed by kappa Z kappa^-1 =
      q^2 Z; the e0 generator gives kappa' = B C^-1, Z' |n> = b q^{2n} |n-1>, b = i y_d q^s/(q^2 - 1); kappa kappa' = 1
      and Z Z' = Z' Z = a b q^2. The untwisted e1 (x) 1 + k1 (x) e1 has no solution.
  (c) Fusing onto a boundary. K_f = Tbar (1 (x) K) T satisfies the reflection equation on V (x) V (x) W (interior of a
      window) for the Ghoshal-Zamolodchikov matrix at generic parameters and for the Dirichlet matrix, and boundary
      unitarity; with Y multiplied by 1.1 in T only, it fails. For scalar factors (Lee-Yang), K_f is crossing-unitary.
  (d) Lee-Yang. S = f_{2/3} = -(2)(4); K_b = Tbar_b K_(1) T_b; among the frozen particles f_{2/3}(theta -+ theta_d) the
      fusion equation has only the solutions theta_d = +- theta_b (mod 2 pi i), theta_b = i pi (3 - b)/6; T(theta)/T(-theta)
      = Tbar T for any T with unitarity and crossing. K_b has a simple pole at theta = i pi (b+1)/6, the pole of Tbar_b,
      with T_b and K_(1) finite there. Fusion maps the defect bootstrap to the boundary bootstrap: Tbar^exc K_(1) T^exc =
      S(theta - iu) S(theta + iu) K_b = K_{b+2}, T^exc = T_{b+2}. The energies: eps_bdry(Phi(b)) - eps_bdry(1) =
      m sin(b pi/6) = m cosh theta_b.
  (e) Freezing is fusing onto Dirichlet (sympy). B(Phi, psi0) = -M cos((Phi - Phi0)/2) + const with
      (M/4 m0) e^{+- i Phi0/2} = cosh(theta_d +- i psi0/2), and every complex C0 arises from real (theta_d, psi0); a
      constant Psi solves the bulk equation only if sin psi0 = 0. At a Dirichlet boundary the flux of q_- + q_+ is
      -8 m0^2 sin(psi0) d_t Phi_x. The fused boundary Phi_x = -B_Phi(Phi, psi0) conserves (1/8m0^3) int_{x<0} (q_- + q_+)
      + D_{-3} + D_{+3} - (sin psi0/m0) Psi_x, with Lecture 6's D_{-+3} at Psi = psi0, Psi_t = 0, Psi_x = Phi_t + B_Psi;
      without the last term, or with q_- - q_+, it fails.
  (f) Quantum fusion onto Dirichlet (sympy). On each charge sector, Tbar (1 (x) K_D) T = (4/y_d)/(1 + y^-2 y_d^-2) times
      the Ghoshal-Zamolodchikov matrix at (eta_D, vartheta = lambda theta_d) after a diagonal change of basis; vartheta =
      -lambda theta_d fails. Numerically on a window the gauge invariants do not depend on n. The scalar factors agree:
      sigma(i lambda theta_d, -i theta) = (1 + y_d^-2) Phi(theta - theta_d)/[(1 + y^-2 y_d^-2) Phi(-theta - theta_d)], with
      sigma from its Gamma product (extrapolated) and Phi from Lecture 7's; it fails for vartheta = (lambda+1) theta_d. The
      integral representation of ln sigma used in the proof agrees with the product, and so does the identity
      int_0^inf (1 - cos ku)/(k sinh(pi k/2 lambda)) dk = ln cosh(lambda u).
  (g) Sinh-Gordon. With T = Tcl, Tbar(theta) = Tcl(i pi - theta), Tbar K_D T = K(theta; E, F = 2 i theta_d/pi) for the
      Corrigan-Delius factor in h = 2 blocks; F = -2 i theta_d/pi gives the same (even in F), F = i theta_d/pi fails.
  (h) What is left of hbar. The classical fused boundary has C_{0,1} = cosh(theta_cl -+ ... ) and the classical limit of the
      boundary UV-IR relation, cos((eta -+ i vartheta)/(lambda+1)) = C_{1,0} (each node with its own parameter; the
      book's sign of phi_0, eq-bsg-dirichlet), gives eta = -(lambda+1) psi0/2, vartheta = (lambda+1) theta_cl; with vartheta = lambda theta_d this is lambda theta_d = (lambda+1) theta_cl. At m0 = 0 the defect
      transmits left-movers unchanged and flips right-movers, so Sigma(0, t) = 2 beta phi_L. Weights: h(e^{i Sigma/2}) =
      beta^2/8 pi = 1/(lambda+1), and (8 pi/beta^2)^2 h = lambda + 1, the spin of the non-local currents; theta_cl =
      (1 - h) theta_d is lambda theta_d = (lambda+1) theta_cl; continued to sinh-Gordon it is (1 - B/2) theta_cl, whose
      expansion is Bajnok and Simon's one-loop 1 - b^2/8 pi; the Lee-Yang couplings of Bajnok, Hollo and Watts are
      proportional to e^{+-(1-h) theta_b} with 1 - h = 6/5.
  (i) The classical breather. T_{B_n} = prod_k Tcl(theta + i pi (n+1-2k)/(2 lambda)) satisfies the bulk bootstrap with
      B_1, unitarity and crossing. At fixed alpha, with n pi/(2 lambda) = pi/2 - alpha and lambda -> inf:
      (pi/lambda) ln T_{B_n} -> int_alpha^{pi-alpha} ln z(theta - theta_d + i x) dx, T_{B_{n+1}}/T_{B_n} -> z/|z| at
      theta + i alpha, and the Wigner delay (d delta_n/d theta)/(m_n sinh theta) -> ln|z(u + i alpha)|/(m0 cos alpha sinh
      theta), with errors that fall like 1/lambda^2, 1/lambda and 1/lambda. The classical breather with Y_{1,2} -> z(u -+
      ... ) Y_{1,2} on the right satisfies the defect conditions (mpmath, 25 digits); with the conjugate factors it fails;
      its envelope is displaced by -ln|z(u + i alpha)|/(m0 cos alpha cosh theta).
  (j) Folding the diagram (a_r, h = 3..6). The reflection rho of the diagram satisfies rho pi rho = pi^-1 and maps the
      mode of particle k to that of h - k, so the mirror image transmits k with T_{h-k}. For a defect at sigma1 followed
      by the mirror at sigma2, the factors of k and h - k agree for all theta exactly when sigma2 = +- sigma1; at sigma2 =
      -sigma1 the pair is transparent; at sigma2 = sigma1 the folded particle crosses with 1/f_{k/h}(theta - theta_d),
      and folded solitons with (cosh u + sin gamma_k)/(cosh u - sin gamma_k) = z_k z_{h-k}, which for h = 3 is Corrigan
      and Zambon's Tzitzeica delay coth(u/2 - i pi/12) coth(u/2 + i pi/12). For h = 2 the pair is the type-II product.

Runtime: about two minutes.
"""

import cmath
import math
import warnings

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.optimize import minimize
from scipy.special import loggamma

warnings.filterwarnings('ignore')
OK = True


def report(name, ok):
    global OK
    OK = OK and bool(ok)
    print(("PASS " if ok else "FAIL ") + name)


rng = np.random.default_rng(9)
mp.mp.dps = 25
I = 1j
PI = math.pi

# ----------------------------------------------------------------------------------------------------------------
# common functions
# ----------------------------------------------------------------------------------------------------------------


def zcl(u):
    """the classical soliton factor of Lecture 4, z(u) = -coth(u/2)"""
    return -np.cosh(u / 2) / np.sinh(u / 2)


def Tcl(u):
    """the classical transmission factor of linear waves, Tcl(u) = z(u + i pi/2)"""
    return zcl(u + I * PI / 2)


def blk(x, th):
    """Lee-Yang blocks (x), period 12"""
    return np.sinh(th / 2 + I * PI * x / 12) / np.sinh(th / 2 - I * PI * x / 12)


def f23(th):
    return (np.sinh(th) + I * math.sin(2 * PI / 3)) / (np.sinh(th) - I * math.sin(2 * PI / 3))


def fx(x, th):
    return (np.sinh(th) + I * np.sin(PI * x)) / (np.sinh(th) - I * np.sin(PI * x))


def thb(b):
    return I * PI * (3 - b) / 6


def K1(th):
    return blk(1, th) * blk(3, th) / blk(4, th)


def Kb(th, b):
    return K1(th) * blk(1 + b, th) * blk(5 - b, th) / (blk(1 - b, th) * blk(5 + b, th))


def lnPhi_gamma(u, lam, N=300, nmax=12):
    """Lecture 7's scalar factor: the ratio r_N/rbar_N at a common cutoff, with a Stirling tail"""
    a = -1j * lam * u / (2 * np.pi)
    k = np.arange(N)

    def lr(aa):
        return (loggamma(k * lam + 0.25 + aa) + loggamma((k + 1) * lam + 0.75 + aa)
                - loggamma((k + 0.5) * lam + 0.25 + aa) - loggamma((k + 0.5) * lam + 0.75 + aa))

    s = np.sum(lr(a) - lr(-a))
    am = mp.mpc(a)
    tail = mp.mpc(0)
    for n in range(2, nmax + 1):
        C = mp.mpc(0)
        for sj, cj in zip((1, 1, -1, -1), (0.25, lam + 0.75, lam / 2 + 0.25, lam / 2 + 0.75)):
            C += sj * (mp.bernpoly(n + 1, cj + am) - mp.bernpoly(n + 1, cj - am))
        tail += C * (-1) ** (n + 1) / mp.mpf(n * (n + 1)) * mp.zeta(n, N) / mp.mpf(lam) ** n
    return 1j * np.pi * (1 + lam) / 4 - np.log(1 + 1j * np.exp(lam * u)) + s + complex(tail)


def V_num(y, q):
    return dict(k1=np.diag([q, 1 / q]), k0=np.diag([1 / q, q]), e1=np.array([[0, y], [0, 0]]),
                f1=np.array([[0, 0], [1 / y, 0]]), e0=np.array([[0, 0], [y, 0]]), f0=np.array([[0, 1 / y], [0, 0]]))


def Delta_num(A, B, x):
    dA, dB = A['k1'].shape[0], B['k1'].shape[0]
    if x[0] == 'e':
        return np.kron(A[x], np.eye(dB)) + np.kron(A['k' + x[1]], B[x])
    if x[0] == 'f':
        return np.kron(A[x], np.linalg.inv(B['k' + x[1]])) + np.kron(np.eye(dA), B[x])
    return np.kron(A[x], B[x])


P4 = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], dtype=complex)


def Rcheck(q, ya, yb):
    """Lecture 7's Rcheck: V(ya) (x) V(yb) -> V(yb) (x) V(ya), basis v+v+, v+v-, v-v+, v-v-"""
    z = ya / yb
    al = q * (z ** 2 - 1) / (q ** 2 * z ** 2 - 1)
    be = (q ** 2 - 1) * z / (q ** 2 * z ** 2 - 1)
    R = np.zeros((4, 4), dtype=complex)
    R[0, 0] = R[3, 3] = 1
    R[1, 1] = be
    R[2, 1] = al
    R[1, 2] = al
    R[2, 2] = be
    return R


class Window:
    def __init__(self, K):
        self.K = K
        self.ns = np.arange(-K, K + 1)
        self.d = 2 * K + 1

    def i(self, n):
        return n + self.K


def T_win(Wn, q, s, Y):
    """T: V (x) W -> W (x) V (matrix part); V (x) W index a d + i(n), W (x) V index 2 i(n) + a; a = 0 v+, 1 v-"""
    d = Wn.d
    T = np.zeros((2 * d, 2 * d), dtype=complex)
    for n in Wn.ns:
        T[2 * Wn.i(n), Wn.i(n)] = q ** (-n)
        if n + 1 <= Wn.K:
            T[2 * Wn.i(n + 1) + 1, Wn.i(n)] = -I * Y * q ** (-n)
        if n - 1 >= -Wn.K:
            T[2 * Wn.i(n - 1), d + Wn.i(n)] = -I * Y * q ** (n + s)
        T[2 * Wn.i(n) + 1, d + Wn.i(n)] = q ** (n + s)
    return T


def Tinv_win(Wn, q, s, Y):
    """the exact inverse of T, W (x) V -> V (x) W, from the 2x2 blocks (without the truncation of the window)"""
    d = Wn.d
    M = np.zeros((2 * d, 2 * d), dtype=complex)
    c = 1 / (1 + Y ** 2)
    for n in Wn.ns:
        # |n> v+  ->  q^n v+|n> + i Y q^{-n-1-s} v-|n+1>
        M[Wn.i(n), 2 * Wn.i(n)] = c * q ** n
        if n + 1 <= Wn.K:
            M[d + Wn.i(n + 1), 2 * Wn.i(n)] = c * I * Y * q ** (-n - 1 - s)
        # |n> v-  ->  i Y q^{n-1} v+|n-1> + q^{-n-s} v-|n>
        if n - 1 >= -Wn.K:
            M[Wn.i(n - 1), 2 * Wn.i(n) + 1] = c * I * Y * q ** (n - 1)
        M[d + Wn.i(n), 2 * Wn.i(n) + 1] = c * q ** (-n - s)
    return M


def swap_VW(Wn):
    """P: W (x) V -> V (x) W"""
    d = Wn.d
    P = np.zeros((2 * d, 2 * d))
    for n in Wn.ns:
        for a in (0, 1):
            P[a * d + Wn.i(n), 2 * Wn.i(n) + a] = 1
    return P


def K_GZ(y, q, e1, e0):
    """the matrix part of the Ghoshal-Zamolodchikov matrix in the book's form (eq-bsg-k-matrix), V(y) -> V(1/y)"""
    kap = q - 1 / q
    return np.array([[kap * (e1 * y + e0 / y), y ** 2 - y ** -2], [y ** 2 - y ** -2, kap * (e0 * y + e1 / y)]])


# ----------------------------------------------------------------------------------------------------------------
# (a) Folding the line
# ----------------------------------------------------------------------------------------------------------------
ok_a = True


def doubled_checks(T, Tb, S, ths):
    """returns (max unitarity residual of K, max defect-unitarity residual, max GZ crossing residual,
    max defect-crossing residual) for the doubled boundary of a diagonal self-conjugate defect"""
    ru, rdu, rc, rdc = 0, 0, 0, 0
    for th in ths:
        K = lambda t: np.array([[0, Tb(t)], [T(t), 0]])
        ru = max(ru, np.abs(K(th) @ K(-th) - np.eye(2)).max())
        rdu = max(rdu, abs(T(th) * Tb(-th) - 1), abs(Tb(th) * T(-th) - 1))

        def kvec(t):  # basis e_a (x) e_b, a, b in {V, V'}, index 2a + b; abar = a
            Kc = K(I * PI / 2 - t)
            return np.array([Kc[b, a] for a in (0, 1) for b in (0, 1)])

        Sc = np.zeros((4, 4), dtype=complex)  # Scheck(2 theta): flip, amplitude S on same-copy pairs, 1 across
        for a in (0, 1):
            for b in (0, 1):
                Sc[2 * b + a, 2 * a + b] = S(2 * th) if a == b else 1
        rc = max(rc, np.abs(kvec(th) - Sc @ kvec(-th)).max())
        rdc = max(rdc, abs(T(th) - Tb(I * PI - th)))
    return ru, rdu, rc, rdc


ths_a = [0.31, -0.77, 0.4 + 0.15j, 1.3 - 0.2j]
b_a = 0.83
res_ly = doubled_checks(lambda t: f23(t - thb(b_a)), lambda t: f23(t + thb(b_a)), f23, ths_a)
thd_a = 0.37
res_cl = doubled_checks(lambda t: Tcl(t - thd_a), lambda t: Tcl(I * PI - t - thd_a),
                        lambda t: fx(-0.3, t), ths_a)
res_bad = doubled_checks(lambda t: f23(t - thb(b_a)) * np.exp(0.2 * t), lambda t: f23(t + thb(b_a)) * np.exp(0.2 * t),
                         f23, ths_a)
ok_a = ok_a and max(res_ly) < 1e-12 and max(res_cl) < 1e-12
ok_a = ok_a and res_bad[0] < 1e-12 and res_bad[1] < 1e-12 and res_bad[2] > 1e-2 and res_bad[3] > 1e-2
# the sine-Gordon soliton matrix: the doubled crossing-unitarity, component by component
lam_a = math.sqrt(7) + 0.2
yd_a = math.exp(lam_a * 0.23)


def T_entry(q, s, th, a, n, b, n2, lam=lam_a, yd=yd_a):
    """<n2| (v_b component of) T(theta) v_a |n> with the scalar factor, s = offset"""
    Y = cmath.exp(lam * th) / yd
    ph = cmath.exp(lnPhi_gamma(th - math.log(yd) / lam, lam)) * q ** (-s / 2)
    if a == 0 and b == 0 and n2 == n:
        return ph * q ** (-n)
    if a == 0 and b == 1 and n2 == n + 1:
        return ph * (-I * Y) * q ** (-n)
    if a == 1 and b == 0 and n2 == n - 1:
        return ph * (-I * Y) * q ** (n + s)
    if a == 1 and b == 1 and n2 == n:
        return ph * q ** (n + s)
    return 0


def Tbar_entry(q, s, th, b, n, a, n2, lam=lam_a, yd=yd_a):
    """<n2| (v_a component of) Tbar(theta) (|n> v_b), Tbar(theta) = T(-theta)^-1, from the 2x2 blocks"""
    thp = -th
    Y = cmath.exp(lam * thp) / yd
    ph = cmath.exp(lnPhi_gamma(thp - math.log(yd) / lam, lam)) * q ** (-s / 2)
    c = 1 / (ph * (1 + Y ** 2))
    if b == 0 and a == 0 and n2 == n:
        return c * q ** n
    if b == 0 and a == 1 and n2 == n + 1:
        return c * I * Y * q ** (-n - 1 - s)
    if b == 1 and a == 0 and n2 == n - 1:
        return c * I * Y * q ** (n - 1)
    if b == 1 and a == 1 and n2 == n:
        return c * q ** (-n - s)
    return 0


def sg_doubled_crossing(q, s):
    """max over components of |K_{V'<-V}(i pi/2 - t)_{b<-abar} - K_{V<-V'}(i pi/2 + t)_{a<-bbar}| / size, with the defect
    labels unchanged: T_{abar -> b}(i pi/2 - t) = Tbar_{bbar -> a}(i pi/2 + t)"""
    worst, size = 0, 0
    for t in (0.21, -0.45):
        for a in (0, 1):
            for b in (0, 1):
                for n in (-1, 0, 2):
                    for n2 in (n - 1, n, n + 1):
                        lhs = T_entry(q, s, I * PI / 2 - t, 1 - a, n, b, n2)
                        rhs = Tbar_entry(q, s, I * PI / 2 + t, 1 - b, n, a, n2)
                        worst = max(worst, abs(lhs - rhs))
                        size = max(size, abs(lhs))
    return worst / size


q_a = -cmath.exp(-I * PI * lam_a)
cr_phys = sg_doubled_crossing(q_a, 0.0)
cr_minus = sg_doubled_crossing(-q_a, 0.0)
cr_s = sg_doubled_crossing(q_a, 0.41)
ok_a = ok_a and cr_phys < 1e-9 and cr_s < 1e-9 and cr_minus > 0.5
# the defect operator of Bajnok and George on a particle: Tbar(i pi/2 - theta) = T(theta + i pi/2) by crossing; for the
# classical factor this is the factor of the turned defect, theta_d -> theta_d - i pi/2
for th in (0.4, -1.1, 0.2 + 0.3j):
    ok_a = ok_a and abs(Tcl(I * PI - (I * PI / 2 - th) - thd_a) - Tcl(th + I * PI / 2 - thd_a)) < 1e-12
    ok_a = ok_a and abs(Tcl(th - (thd_a - I * PI / 2)) - Tcl(th + I * PI / 2 - thd_a)) < 1e-12
    ok_a = ok_a and abs(f23(I * PI / 2 - th + thb(b_a)) - f23(th + I * PI / 2 - thb(b_a))) < 1e-12
report(f"(a) folding: doubled unitarity <=> defect unitarity, GZ crossing-unitarity <=> defect crossing (Lee-Yang, "
       f"classical factor; with T e^(0.2 theta) both crossings fail by {res_bad[2]:.2f}); sine-Gordon solitons: doubled "
       f"crossing with defect labels unchanged holds at q ({cr_phys:.0e}), fails at -q; the Bajnok-George operator "
       f"multiplies a particle by T(theta + i pi/2), the factor of the turned defect", ok_a)

# ----------------------------------------------------------------------------------------------------------------
# (b) The reflected copy and the other Borel half
# ----------------------------------------------------------------------------------------------------------------
ok_b = True
q_b = cmath.exp(I * 0.873)
y1, y2 = cmath.exp(0.4 + 0.1j), cmath.exp(-0.3)
A1, A2 = V_num(y1, q_b), V_num(y2, q_b)
Rc = Rcheck(q_b, y1, y2)
# Rcheck intertwines Delta (check our formula first)
for x in ('e1', 'e0', 'f1', 'f0', 'k1'):
    ok_b = ok_b and np.abs(Rc @ Delta_num(A1, A2, x) - Delta_num(A2, A1, x) @ Rc).max() < 1e-12
PRP = P4 @ Rc @ P4  # V(y2) (x) V(y1) -> V(y1) (x) V(y2)
worst_op, worst_plain = 0, 0
for x in ('e1', 'e0', 'f1', 'f0'):
    Dop_21 = P4 @ Delta_num(A1, A2, x) @ P4  # Delta^op on V(y2) (x) V(y1)
    Dop_12 = P4 @ Delta_num(A2, A1, x) @ P4  # Delta^op on V(y1) (x) V(y2)
    worst_op = max(worst_op, np.abs(PRP @ Dop_21 - Dop_12 @ PRP).max())
    worst_plain = max(worst_plain, np.abs(PRP @ Delta_num(A2, A1, x) - Delta_num(A1, A2, x) @ PRP).max())
ok_b = ok_b and worst_op < 1e-12 and worst_plain > 1e-2

# the second module on a window
Kb_w, s_b, yd_b = 9, 0.29, math.exp(0.6)
Wb = Window(Kb_w)
d = Wb.d


def L_win(Wn, q, s, Y):
    return swap_VW(Wn) @ T_win(Wn, q, s, Y)


def opV(M2, Wn):
    return np.kron(M2, np.eye(Wn.d))


def opW(MW):
    return np.kron(np.eye(2), MW)


def shift_op(vals, sh, Wn):
    """|n> -> vals[n] |n + sh>"""
    M = np.zeros((Wn.d, Wn.d), dtype=complex)
    for n in Wn.ns:
        if -Wn.K <= n + sh <= Wn.K:
            M[Wn.i(n + sh), Wn.i(n)] = vals[Wn.i(n)]
    return M


def solve_twisted(q, s, yd, gen, ys, untwisted=False):
    """least-squares solution for Z (shift sh1) and kappa (shift sh2) on interior columns; returns (residual, Z, kap)"""
    if gen == 'e1':
        sh1, sh2 = 1, (1 if untwisted else 2)
    else:
        sh1, sh2 = -1, (-1 if untwisted else -2)
    nun = 2 * d
    rows, rhs = [], []
    cols = [a * d + Wb.i(n) for a in (0, 1) for n in Wb.ns if abs(n) <= Wb.K - 4]
    for yv in ys:
        V = V_num(yv, q)
        L = L_win(Wb, q, s, yv / yd)
        if gen == 'e1':
            inhom = L @ opV(V['e1'], Wb)
            k_op, out_op = V['k1'], (V['e1'] if untwisted else V['e0'])
        else:
            inhom = L @ opV(V['e0'], Wb)
            k_op, out_op = V['k0'], (V['e0'] if untwisted else V['e1'])
        # L (k (x) Z) - (1 (x) Z) L - (out (x) kappa) L = -inhom, linear in the entries of Z and kappa
        blocks = []
        for j in range(d):
            e = np.zeros(d, dtype=complex)
            e[j] = 1
            Zj = shift_op(e, sh1, Wb)
            blocks.append(L @ np.kron(k_op, Zj) - opW(Zj) @ L)
        for j in range(d):
            e = np.zeros(d, dtype=complex)
            e[j] = 1
            Kj = shift_op(e, sh2, Wb)
            blocks.append(-np.kron(out_op, Kj) @ L)
        for c in cols:
            Amat = np.array([B[:, c] for B in blocks]).T
            rows.append(Amat)
            rhs.append(-inhom[:, c])
    A = np.vstack(rows)
    bvec = np.concatenate(rhs)
    sol, *_ = np.linalg.lstsq(A, bvec, rcond=None)
    res = np.linalg.norm(A @ sol - bvec) / np.linalg.norm(bvec)
    return res, sol[:nun // 2], sol[nun // 2:]


ys_b = [cmath.exp(0.3), cmath.exp(-0.5 + 0.2j), cmath.exp(0.9 - 0.1j)]
res1, Zs, kaps = solve_twisted(q_b, s_b, yd_b, 'e1', ys_b)
inner = [n for n in Wb.ns if abs(n) <= Wb.K - 5]
a_b = I * yd_b * q_b ** (-s_b) / (q_b ** 2 - 1)
b_b = I * yd_b * q_b ** s_b / (q_b ** 2 - 1)
kap_ok = max(abs(kaps[Wb.i(n)] - q_b ** (-2 * n - 2 - s_b)) for n in inner)
dZ_ok = max(abs((Zs[Wb.i(n + 1)] - Zs[Wb.i(n)]) - (a_b * q_b ** (-2 * n - 2) - a_b * q_b ** (-2 * n))) for n in inner[:-1])
res0, Zs0, kaps0 = solve_twisted(q_b, s_b, yd_b, 'e0', ys_b)
kap0_ok = max(abs(kaps0[Wb.i(n)] - q_b ** (2 * n - 2 + s_b)) for n in inner)
dZ0_ok = max(abs((Zs0[Wb.i(n + 1)] - Zs0[Wb.i(n)]) - (b_b * q_b ** (2 * n + 2) - b_b * q_b ** (2 * n))) for n in inner[:-1])
ok_b = ok_b and res1 < 1e-12 and res0 < 1e-12 and kap_ok < 1e-10 and dZ_ok < 1e-10 and kap0_ok < 1e-10 and dZ0_ok < 1e-10
# the module: kappa Z kappa^-1 = q^2 Z fixes the constant; kappa kappa' = 1; Z Z' = Z' Z = a b q^2
Kap = shift_op(np.array([q_b ** (-2 * n - 2 - s_b) for n in Wb.ns]), 2, Wb)
Kap0 = shift_op(np.array([q_b ** (2 * n - 2 + s_b) for n in Wb.ns]), -2, Wb)
Zm = shift_op(np.array([a_b * q_b ** (-2 * n) for n in Wb.ns]), 1, Wb)
Z0 = shift_op(np.array([b_b * q_b ** (2 * n) for n in Wb.ns]), -1, Wb)
ins = [Wb.i(n) for n in inner]
sub = np.ix_(ins, ins)
ok_b = ok_b and np.abs((Kap @ Kap0)[sub] - np.eye(len(ins))).max() < 1e-12
ok_b = ok_b and np.abs((Kap @ Zm - q_b ** 2 * Zm @ Kap)[sub]).max() < 1e-12
ok_b = ok_b and np.abs((Zm @ Z0)[sub] - a_b * b_b * q_b ** 2 * np.eye(len(ins))).max() < 1e-10
ok_b = ok_b and np.abs((Z0 @ Zm)[sub] - a_b * b_b * q_b ** 2 * np.eye(len(ins))).max() < 1e-10
# kappa = C B^-1 with C|n> = -i/y_d q^-n |n+1>, B|n> = -i/y_d q^{n+s} |n-1>
Cm = shift_op(np.array([-I / yd_b * q_b ** (-n) for n in Wb.ns]), 1, Wb)
Bm = shift_op(np.array([-I / yd_b * q_b ** (n + s_b) for n in Wb.ns]), -1, Wb)
Binv = shift_op(np.array([1 / (-I / yd_b * q_b ** (n + 1 + s_b)) for n in Wb.ns]), 1, Wb)
ok_b = ok_b and np.abs((Cm @ Binv - Kap)[sub]).max() < 1e-12
res_un, _, _ = solve_twisted(q_b, s_b, yd_b, 'e1', ys_b, untwisted=True)
ok_b = ok_b and res_un > 0.1
report(f"(b) the reflected copy: P Rcheck P intertwines Delta^op ({worst_op:.0e}) and not Delta; the doubled boundary "
       f"keeps (1 (x) tau) Delta(U_q(b+)) on a second module: kappa = C B^-1 = q^(-2n-2-s) shift 2, Z = a q^(-2n) "
       f"(residual {res1:.0e}), kappa' = B C^-1, kappa kappa' = 1, Z Z' = a b q^2; the untwisted e1 fails "
       f"(residual {res_un:.2f})", ok_b)

# ----------------------------------------------------------------------------------------------------------------
# (c) Fusing onto a boundary: the reflection equation
# ----------------------------------------------------------------------------------------------------------------
ok_c = True
lam_c = math.sqrt(3) + 0.41
q_c = -cmath.exp(-I * PI * lam_c)
Kc_w, s_c, thd_c = 11, 0.0, 0.27
Wc = Window(Kc_w)
dc = Wc.d
yd_c = math.exp(lam_c * thd_c)


def Kf_win(th, Kmat, bad=1.0):
    """K_f(theta) = Tbar (1_W (x) K) T: V (x) W -> V (x) W, matrix part"""
    y = cmath.exp(lam_c * th)
    T = T_win(Wc, q_c, s_c, bad * y / yd_c)
    Tb = Tinv_win(Wc, q_c, s_c, (1 / y) / yd_c)  # Tbar(theta) = T(-theta)^-1
    return Tb @ np.kron(np.eye(dc), Kmat(y)) @ T


def check_RE(Kmat, bad=1.0):
    t1, t2 = 0.47 + 0.05j, -0.18
    y1, y2 = cmath.exp(lam_c * t1), cmath.exp(lam_c * t2)
    I2W = np.eye(2 * dc)
    R = lambda ya, yb: np.kron(Rcheck(q_c, ya, yb), np.eye(dc))
    K2 = lambda t: np.kron(np.eye(2), Kf_win(t, Kmat, bad))
    lhs = R(y1, y2) @ K2(t1) @ R(y1, 1 / y2) @ K2(t2)  # Scheck(t1 - t2) K(t1) Scheck(t1 + t2) K(t2)
    rhs = K2(t2) @ R(y2, 1 / y1) @ K2(t1) @ R(y1, y2)
    cols = [a * 2 * dc + b * dc + Wc.i(n) for a in (0, 1) for b in (0, 1) for n in Wc.ns if abs(n) <= 2]
    return np.abs((lhs - rhs)[:, cols]).max() / np.abs(lhs[:, cols]).max()


e1g, e0g = 0.37 - 0.2j, 1.1 + 0.4j
KGZ = lambda y: K_GZ(y, q_c, e1g, e0g)
etaD = 0.61
KD = lambda y: np.diag([y * cmath.exp(I * etaD) + cmath.exp(-I * etaD) / y, y * cmath.exp(-I * etaD) + cmath.exp(I * etaD) / y])
# K itself satisfies the reflection equation with this Rcheck
t1, t2 = 0.47 + 0.05j, -0.18
y1, y2 = cmath.exp(lam_c * t1), cmath.exp(lam_c * t2)
for Km in (KGZ, KD):
    lhsK = Rcheck(q_c, y1, y2) @ np.kron(np.eye(2), Km(y1)) @ Rcheck(q_c, y1, 1 / y2) @ np.kron(np.eye(2), Km(y2))
    rhsK = np.kron(np.eye(2), Km(y2)) @ Rcheck(q_c, y2, 1 / y1) @ np.kron(np.eye(2), Km(y1)) @ Rcheck(q_c, y1, y2)
    ok_c = ok_c and np.abs(lhsK - rhsK).max() < 1e-10 * np.abs(lhsK).max()
re_gz, re_d, re_bad = check_RE(KGZ), check_RE(KD), check_RE(KGZ, bad=1.1)
ok_c = ok_c and re_gz < 1e-12 and re_d < 1e-12 and re_bad > 1e-3
# boundary unitarity of K_f, up to a scalar
th_u = 0.33 + 0.1j
U = Kf_win(-th_u, KGZ) @ Kf_win(th_u, KGZ)
colsU = [a * dc + Wc.i(n) for a in (0, 1) for n in Wc.ns if abs(n) <= 4]
Uc = U[:, colsU]
scal = U[Wc.i(0), Wc.i(0)]
ok_c = ok_c and np.abs(Uc - scal * np.eye(2 * dc)[:, colsU]).max() < 1e-10 * abs(scal)
# scalar fusion is crossing-unitary (Lee-Yang)
for b in (0.4, 1.7):
    for th in (0.3, -0.6 + 0.1j):
        Kf = lambda t: f23(t + thb(b)) * K1(t) * f23(t - thb(b))
        ok_c = ok_c and abs(Kf(I * PI / 2 - th) - f23(2 * th) * Kf(I * PI / 2 + th)) < 1e-12
        ok_c = ok_c and abs(K1(I * PI / 2 - th) - f23(2 * th) * K1(I * PI / 2 + th)) < 1e-12
report(f"(c) fusing: Tbar (1 (x) K) T satisfies the reflection equation for the generic GZ matrix ({re_gz:.0e}) and "
       f"for Dirichlet ({re_d:.0e}), and fails with Y -> 1.1 Y in T ({re_bad:.0e}); boundary unitarity; the fused "
       f"Lee-Yang factor is crossing-unitary", ok_c)

# ----------------------------------------------------------------------------------------------------------------
# (d) Lee-Yang
# ----------------------------------------------------------------------------------------------------------------
ok_d = True
th_d = np.array([0.3, -0.7, 1.1 + 0.2j, 0.05 + 0.4j])
ok_d = ok_d and np.abs(f23(th_d) + blk(2, th_d) * blk(4, th_d)).max() < 1e-13
found_all = True
for b in (0.5, 1.37, -2.2):
    tb = thb(b)
    ok_d = ok_d and np.abs(f23(th_d + tb) * K1(th_d) * f23(th_d - tb) - Kb(th_d, b)).max() < 1e-12
    ths = np.array([0.31, -0.83, 1.27, 0.11 + 0.2j, -0.4 + 0.1j])
    Fb = Kb(ths, b) / K1(ths)

    def G(v, Fb=Fb, ths=ths):
        td = v[0] + I * v[1]
        return np.sum(np.abs(f23(ths + td) * f23(ths - td) - Fb) ** 2)

    sols = []
    for re in np.linspace(-2.5, 2.5, 6):
        for im in np.linspace(-PI, PI, 9):
            r = minimize(G, [re, im], method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-26, maxiter=4000))
            if r.fun < 1e-18:
                td = complex(r.x[0], (r.x[1] + PI) % (2 * PI) - PI)
                if not any(abs(td - t) < 1e-5 for t in sols):
                    sols.append(td)
    tbm = complex(0, (tb.imag + PI) % (2 * PI) - PI)
    expect = [tbm, complex(0, (-tb.imag + PI) % (2 * PI) - PI)]
    found_all = found_all and len(sols) == 2 and all(min(abs(s - e) for s in sols) < 1e-6 for e in expect)
ok_d = ok_d and found_all
# T(theta)/T(-theta) = Tbar(theta) T(theta) for any T with unitarity and crossing: take T = T_b times an even CDD-like
# factor E with E(theta) E(theta + i pi) = 1 and E even: E = f_{1/5}(theta - 0.3) f_{1/5}(-theta - 0.3) fails the first,
# so test the algebraic identity on T_b itself and on T_b f_x(theta - c) (unitarity and crossing kept)
for c0 in (0.0, 0.6 + 0.2j):
    T = lambda t: f23(t - thb(0.9)) * fx(0.27, t - c0)
    Tb = lambda t: T(I * PI - t)  # crossing
    for th in (0.4, -0.3 + 0.2j):
        ok_d = ok_d and abs(T(th) * Tb(-th) - 1) < 1e-12  # unitarity holds for this T
        ok_d = ok_d and abs(T(th) / T(-th) - Tb(th) * T(th)) < 1e-12
# the pole: K_b at i u, u = pi (b + 1)/6, from Tbar_b
b = 0.5
u = PI * (b + 1) / 6
res_K, res_T = [], []
for eps in (1e-4, 1e-5):
    res_K.append(eps * Kb(I * u + eps, b))
    res_T.append(eps * f23(I * u + eps + thb(b)) * K1(I * u) * f23(I * u - thb(b)))
ok_d = ok_d and abs(res_K[1] - res_K[0]) < 1e-3 * abs(res_K[1]) and abs(res_K[1]) > 1e-2
ok_d = ok_d and abs(res_K[1] - res_T[1]) < 1e-3 * abs(res_K[1])
ok_d = ok_d and np.isfinite(K1(I * u)) and abs(K1(I * u)) > 1e-3 and abs(f23(I * u - thb(b))) > 1e-3
# fusion maps the defect bootstrap onto the boundary bootstrap
for th in (0.2, -0.9 + 0.1j, 0.6 + 0.3j):
    Texc = f23(th + I * u) * f23(th - thb(b))
    Tbexc = f23(th - I * u) * f23(th + thb(b))
    ok_d = ok_d and abs(Texc - f23(th - thb(b + 2))) < 1e-12 and abs(Tbexc - f23(th + thb(b + 2))) < 1e-12
    ok_d = ok_d and abs(Tbexc * K1(th) * Texc - f23(th - I * u) * f23(th + I * u) * Kb(th, b)) < 1e-12
    ok_d = ok_d and abs(f23(th - I * u) * f23(th + I * u) * Kb(th, b) - Kb(th, b + 2)) < 1e-12
# the energies (Bajnok and Simon, Dorey et al.)
for bb in (-2.1, 0.5, 1.9):
    e1 = (math.sqrt(3) - 1) / 2
    eb = (math.sqrt(3) - 1 + 2 * math.sin(bb * PI / 6)) / 2
    ok_d = ok_d and abs((eb - e1) - math.sin(bb * PI / 6)) < 1e-14 and abs(cmath.cosh(thb(bb)) - math.sin(bb * PI / 6)) < 1e-14
report("(d) Lee-Yang: S = -(2)(4); K_b = Tbar_b K_(1) T_b; the fusion equation for a frozen particle has only "
       "theta_d = +- theta_b; Tbar T = T(theta)/T(-theta); the pole of K_b at i pi (b+1)/6 is that of Tbar_b; fusion maps "
       "the defect bootstrap to the boundary bootstrap, K_b -> K_(b+2); eps_def = m sin(b pi/6) = m cosh theta_b", ok_d)

# ----------------------------------------------------------------------------------------------------------------
# (e) Freezing is fusing onto Dirichlet
# ----------------------------------------------------------------------------------------------------------------
ok_e = True
m0, sig, psi0 = sp.symbols('m0 sigma psi0', positive=True)
Fs, p1, p2, p3, thd_s = sp.symbols('Phi p1 p2 p3 theta_d', real=True)
Bdef = -2 * m0 * (sig * (sp.cos((Fs + psi0) / 2) - 1) + (sp.cos((Fs - psi0) / 2) - 1) / sig)
cC = 2 * sp.cosh(thd_s - I * psi0 / 2)
# B(Phi, psi0) = -2 m0 Re(c e^{i Phi/2}) + const with c = 2 cosh(theta_d - i psi0/2) at sigma = e^{-theta_d}
Bform = -2 * m0 * (cC * sp.exp(I * Fs / 2) + sp.conjugate(cC) * sp.exp(-I * Fs / 2)) / 2
diffB = sp.simplify((Bdef.subs(sig, sp.exp(-thd_s)) - Bform).rewrite(sp.exp).expand())
ok_e = ok_e and sp.simplify(sp.diff(diffB, Fs)) == 0
# every complex C0 arises: C0 = cos(psi0/2 - i theta_d)
for C0 in (0.3 + 2.1j, -1.7 - 0.4j, 5.0):
    w = cmath.acos(C0)
    p0, td = 2 * w.real, -w.imag
    ok_e = ok_e and abs(cmath.cosh(td + I * p0 / 2) - C0) < 1e-12
# a constant Psi = psi0 solves d+ d- Psi = -m0^2 sin Psi only if sin psi0 = 0 (trivially: lhs = 0)
# the spin-3 bulk densities and currents
t_s, x_s = sp.symbols('t x', real=True)
fF = sp.Function('f')(t_s, x_s)
dp = lambda ex: sp.diff(ex, t_s) + sp.diff(ex, x_s)
dm = lambda ex: sp.diff(ex, t_s) - sp.diff(ex, x_s)
qm_ = dp(dp(fF)) ** 2 - dp(fF) ** 4 / 4 + m0 ** 2 * sp.cos(fF) * dp(fF) ** 2
jm_ = m0 ** 2 * sp.cos(fF) * dp(fF) ** 2 - dp(dp(fF)) ** 2 + dp(fF) ** 4 / 4
qp_ = dm(dm(fF)) ** 2 - dm(fF) ** 4 / 4 + m0 ** 2 * sp.cos(fF) * dm(fF) ** 2
jp_ = -(m0 ** 2 * sp.cos(fF) * dm(fF) ** 2 - dm(dm(fF)) ** 2 + dm(fF) ** 4 / 4)


def onshell(ex):
    ex = sp.expand(ex)
    for _ in range(4):
        ex = ex.subs(sp.Derivative(fF, (t_s, 3), x_s), sp.diff(fF, t_s, x_s, x_s, x_s) - m0 ** 2 * sp.diff(sp.sin(fF), t_s, x_s))
        ex = ex.subs(sp.Derivative(fF, (t_s, 3)), sp.diff(fF, t_s, x_s, x_s) - m0 ** 2 * sp.diff(sp.sin(fF), t_s))
        ex = ex.subs(sp.Derivative(fF, (t_s, 2), x_s), sp.diff(fF, x_s, x_s, x_s) - m0 ** 2 * sp.diff(sp.sin(fF), x_s))
        ex = ex.subs(sp.Derivative(fF, (t_s, 2)), sp.diff(fF, x_s, x_s) - m0 ** 2 * sp.sin(fF))
        ex = sp.expand(ex)
    return sp.simplify(ex)


ok_e = ok_e and onshell(sp.diff(qm_, t_s) + sp.diff(jm_, x_s)) == 0 and onshell(sp.diff(qp_, t_s) + sp.diff(jp_, x_s)) == 0


def currents_at(phx, phxt, phtt):
    """j_- and j_+ at x = 0 from the jets Phi, Phi_t = p1, Phi_tt, Phi_x, Phi_xt; Phi_xx from the field equation"""
    phxx = phtt + m0 ** 2 * sp.sin(Fs)
    dpF, dmF = p1 + phx, p1 - phx
    dp2, dm2 = phtt + 2 * phxt + phxx, phtt - 2 * phxt + phxx
    jm0 = m0 ** 2 * sp.cos(Fs) * dpF ** 2 - dp2 ** 2 + dpF ** 4 / 4
    jp0 = -(m0 ** 2 * sp.cos(Fs) * dmF ** 2 - dm2 ** 2 + dmF ** 4 / 4)
    return jm0, jp0, dpF, dmF


# the Dirichlet boundary: Phi = psi0 fixed, Phi_t = Phi_tt = 0, Phi_x = w, Phi_xt = w1
w_s, w1_s = sp.symbols('w w1')
jmD, jpD, _, _ = currents_at(w_s, w1_s, 0)
jmD, jpD = jmD.subs({p1: 0, Fs: psi0}), jpD.subs({p1: 0, Fs: psi0})
ok_e = ok_e and sp.simplify(jmD + jpD + 8 * m0 ** 2 * sp.sin(psi0) * w1_s) == 0
# the fused boundary
BF = sp.diff(Bdef, Fs)
pps = sp.Symbol('pps')
BPsi = sp.diff(Bdef.subs(psi0, pps), pps).subs(pps, psi0)
gx = -BF
jm0, jp0, dpF, dmF = currents_at(gx, sp.diff(gx, Fs) * p1, p2)
Psix = p1 + BPsi
Sig, Dl = Fs + psi0, Fs - psi0
Dm3 = sig / (2 * m0 ** 2) * sp.cos(Sig / 2) * dpF * Psix + sig ** 3 * (sp.cos(Sig / 2) - sp.Rational(2, 3) * sp.cos(Sig / 2) ** 3)
Dp3 = -1 / (2 * m0 ** 2 * sig) * sp.cos(Dl / 2) * dmF * (-Psix) + sig ** -3 * (sp.cos(Dl / 2) - sp.Rational(2, 3) * sp.cos(Dl / 2) ** 3)


def ddt(ex):
    return sp.diff(ex, Fs) * p1 + sp.diff(ex, p1) * p2 + sp.diff(ex, p2) * p3


def fused_residual(c, kb):
    dQ = -(jm0 + c * jp0) / (8 * m0 ** 3) + ddt(Dm3 + c * Dp3 + kb * sp.sin(psi0) / m0 * Psix)
    return sp.simplify(sp.expand(sp.expand_trig(dQ.rewrite(sp.exp))))


ok_e = ok_e and fused_residual(1, -1) == 0
ok_e = ok_e and fused_residual(1, 0) != 0 and fused_residual(-1, -1) != 0
report("(e) freezing is fusing: B(Phi, psi0) is the GZ potential with C0, C1 = cosh(theta_d -+ i psi0/2), every "
       "boundary arises; Dirichlet flux of q_- + q_+ = -8 m0^2 sin psi0 d_t Phi_x; the fused boundary conserves "
       "(1/8m0^3) int (q_- + q_+) + D_-3 + D_+3 - (sin psi0/m0) Psi_x exactly (sympy), not without the last term, "
       "nor with q_- - q_+", ok_e)

# ----------------------------------------------------------------------------------------------------------------
# (f) Quantum fusion onto Dirichlet
# ----------------------------------------------------------------------------------------------------------------
ok_f = True
ys, yds, es, qns, qss, ss = sp.symbols('y y_d e qn q s', nonzero=True)
Ysym, Ypsym = ys / yds, 1 / (ys * yds)
Amat = lambda Z: sp.Matrix([[1, -I * Z], [-I * Z, 1]])
Dn = sp.diag(1 / qns, qns * qss ** (1 + ss))
kp_s, km_s = ys * es + 1 / (ys * es), ys / es + es / ys
Kf_s = sp.simplify(Dn.inv() * Amat(Ypsym).inv() * sp.diag(kp_s, km_s) * Amat(Ysym) * Dn)
cl, sl = (ys + 1 / ys) / 2, (ys - 1 / ys) / (2 * I)
ce, se = (es + 1 / es) / 2, (es - 1 / es) / (2 * I)


def KGZ_sym(ydv):
    ch, sh = (ydv + 1 / ydv) / 2, (ydv - 1 / ydv) / 2
    return sp.Matrix([[cl * ce * ch - sl * se * sh, -sl * cl], [-sl * cl, cl * ce * ch + sl * se * sh]])


Dp_ = sp.diag(1, -1 / es)
cand = 4 / yds / (1 + 1 / (ys ** 2 * yds ** 2)) * Dn.inv() * Dp_ * KGZ_sym(yds) * Dp_.inv() * Dn
ok_f = ok_f and sp.simplify(Kf_s - cand) == sp.zeros(2, 2)
Kw = KGZ_sym(1 / yds)  # vartheta = -lambda theta_d
ok_f = ok_f and sp.simplify(Kf_s[0, 0] / Kf_s[1, 1] - Kw[0, 0] / Kw[1, 1]) != 0
# numerically on a window: the gauge invariants do not depend on n
th_f = 0.41 + 0.07j
Kfw = Kf_win(th_f, KD)
y_f = cmath.exp(lam_c * th_f)
invs = []
for n in range(-3, 4):
    i1, i2 = Wc.i(n), dc + Wc.i(n + 1)
    k11, k22, k12, k21 = Kfw[i1, i1], Kfw[i2, i2], Kfw[i1, i2], Kfw[i2, i1]
    invs.append((k11 / k22, k12 * k21 / (k11 * k22)))
invs = np.array(invs)
ok_f = ok_f and np.abs(invs - invs[0]).max() < 1e-12
P0p = lambda y, yd, eta: ((y + 1 / y) / 2 * math.cos(eta) * (yd + 1 / yd) / 2 - (y - 1 / y) / (2 * I) * math.sin(eta) * (yd - 1 / yd) / 2)
P0m = lambda y, yd, eta: ((y + 1 / y) / 2 * math.cos(eta) * (yd + 1 / yd) / 2 + (y - 1 / y) / (2 * I) * math.sin(eta) * (yd - 1 / yd) / 2)
Q0 = lambda y: -(y - 1 / y) / (2 * I) * (y + 1 / y) / 2
ok_f = ok_f and abs(invs[0][0] - P0p(y_f, yd_c, etaD) / P0m(y_f, yd_c, etaD)) < 1e-12
ok_f = ok_f and abs(invs[0][1] - Q0(y_f) ** 2 / (P0p(y_f, yd_c, etaD) * P0m(y_f, yd_c, etaD))) < 1e-12
yd_wrong = math.exp(thd_c)  # vartheta = theta_d, without lambda
ok_f = ok_f and abs(invs[0][0] - P0p(y_f, yd_wrong, etaD) / P0m(y_f, yd_wrong, etaD)) > 1e-3
# the scalar factors


def lnsigma_prod(x, u, lam, N):
    l = np.arange(1, N + 1)

    def lp(uu):
        w = lam * uu / np.pi
        return (loggamma(0.5 + x / np.pi + (2 * l - 1) * lam - w) + loggamma(0.5 - x / np.pi + (2 * l - 1) * lam - w)
                - loggamma(0.5 - x / np.pi + (2 * l - 2) * lam - w) - loggamma(0.5 + x / np.pi + 2 * l * lam - w))

    return np.log(np.cos(x) / np.cos(x + lam * u)) + np.sum(lp(u) - lp(-u))


def lnsigma(x, u, lam):
    v = np.array([lnsigma_prod(x, u, lam, N) for N in (400, 800, 1600, 3200)])
    r1 = 2 * v[1:] - v[:-1]
    r2 = (4 * r1[1:] - r1[:-1]) / 3
    return r2[-1]


def lnsigma_int(x, u, lam):
    """ln sigma = ln[cos x/cos(x + lam u)] - int_0^inf dt/t sinh(lam u t/pi) sinh(x t/pi + lam t/2)/(sinh(t/2) cosh(lam t/2))"""
    def g(t):
        if t == 0:
            return 0j
        return np.sinh(lam * u * t / PI) * np.sinh(x * t / PI + lam * t / 2) / (t * np.sinh(t / 2) * np.cosh(lam * t / 2))
    re = quad(lambda t: g(t).real, 0, 200, limit=2000, epsabs=1e-13, epsrel=1e-12)[0]
    im = quad(lambda t: g(t).imag, 0, 200, limit=2000, epsabs=1e-13, epsrel=1e-12)[0]
    return np.log(np.cos(x) / np.cos(x + lam * u)) - (re + I * im)


worst_s, worst_bad, worst_int = 0, 1e9, 0
for lam in (math.sqrt(5), math.e + 0.3):
    for thd in (0.4, -0.7):
        for th in (0.3, 0.9 + 0.1j, -0.5):
            uu = -I * th
            y, yd = cmath.exp(lam * th), math.exp(lam * thd)
            rhs = math.log(1 + yd ** -2) - cmath.log(1 + y ** -2 * yd ** -2) + lnPhi_gamma(th - thd, lam) - lnPhi_gamma(-th - thd, lam)
            worst_s = max(worst_s, abs(cmath.exp(lnsigma(I * lam * thd, uu, lam)) - cmath.exp(rhs)))
            worst_bad = min(worst_bad, abs(cmath.exp(lnsigma(I * (lam + 1) * thd, uu, lam)) - cmath.exp(rhs)))
            if th.imag == 0 if isinstance(th, complex) else True:
                worst_int = max(worst_int, abs(lnsigma_int(I * lam * thd, uu, lam) - lnsigma(I * lam * thd, uu, lam)))
ok_f = ok_f and worst_s < 1e-7 and worst_bad > 1e-3 and worst_int < 1e-7
for lam in (1.7, 4.3):
    for uv in (0.3, 1.2):
        val = quad(lambda k: (1 - math.cos(k * uv)) / (k * math.sinh(PI * k / (2 * lam))), 0, 400 * lam, limit=4000)[0]
        ok_f = ok_f and abs(val - math.log(math.cosh(lam * uv))) < 1e-8
report(f"(f) quantum fusion onto Dirichlet: on each charge sector Tbar (1 (x) K_D) T = (4/y_d)/(1 + y^-2 y_d^-2) "
       f"x GZ(eta_D, lambda theta_d) after a diagonal change of basis (sympy), not with -lambda theta_d; gauge invariants "
       f"independent of n; scalar factors agree ({worst_s:.0e}; with (lambda+1) theta_d off by {worst_bad:.0e}); the "
       f"integral form of ln sigma ({worst_int:.0e}) and int (1 - cos ku)/(k sinh) = ln cosh(lambda u)", ok_f)

# ----------------------------------------------------------------------------------------------------------------
# (g) Sinh-Gordon
# ----------------------------------------------------------------------------------------------------------------
ok_g = True


def blk2(x, th):
    """h = 2 blocks (x) = sinh(theta/2 + i pi x/4)/sinh(theta/2 - i pi x/4)"""
    return np.sinh(th / 2 + I * PI * x / 4) / np.sinh(th / 2 - I * PI * x / 4)


def K_shG(th, B, E, F):
    return blk2(1, th) * blk2(2 - B / 2, th) * blk2(1 + B / 2, th) / (blk2(1 - E, th) * blk2(1 + E, th) * blk2(1 - F, th) * blk2(1 + F, th))


def K_D_shG(th, B, E):
    return blk2(1, th) * blk2(2 - B / 2, th) * blk2(1 + B / 2, th) / (blk2(1 - E, th) * blk2(1 + E, th))


Bv, Ev, thd_g = 0.37, 0.21 + 0.1j, 0.53
for th in (0.3, -0.8, 0.4 + 0.2j):
    fus = Tcl(I * PI - th - thd_g) * K_D_shG(th, Bv, Ev) * Tcl(th - thd_g)
    ok_g = ok_g and abs(fus - K_shG(th, Bv, Ev, 2 * I * thd_g / PI)) < 1e-12
    ok_g = ok_g and abs(fus - K_shG(th, Bv, Ev, -2 * I * thd_g / PI)) < 1e-12
    ok_g = ok_g and abs(fus - K_shG(th, Bv, Ev, I * thd_g / PI)) > 1e-3
    # the block identity (1 - F)(1 + F) = (cos(pi F/2) - i sinh)/(cos(pi F/2) + i sinh)
    F = 2 * I * thd_g / PI
    ok_g = ok_g and abs(blk2(1 - F, th) * blk2(1 + F, th) - (cmath.cos(PI * F / 2) - I * cmath.sinh(th)) / (cmath.cos(PI * F / 2) + I * cmath.sinh(th))) < 1e-12
    ok_g = ok_g and abs(Tcl(I * PI - th - thd_g) * Tcl(th - thd_g) - (math.cosh(thd_g) + I * cmath.sinh(th)) / (math.cosh(thd_g) - I * cmath.sinh(th))) < 1e-12
report("(g) sinh-Gordon: Tbar K_D T = K(theta; E, F = 2 i theta_d/pi) for the Corrigan-Delius factor (even in F; "
       "F = i theta_d/pi fails); Tbar T = (cosh theta_d + i sinh theta)/(cosh theta_d - i sinh theta)", ok_g)

# ----------------------------------------------------------------------------------------------------------------
# (h) What is left of hbar
# ----------------------------------------------------------------------------------------------------------------
ok_h = True
for lam in (3.7, 11.2):
    for thcl, p0 in ((0.6, 0.8), (-0.35, 2.3)):
        C0, C1 = cmath.cosh(thcl + I * p0 / 2), cmath.cosh(thcl - I * p0 / 2)
        # classical UV-IR: cos((eta - i vartheta)/(lambda+1)) = C1, cos((eta + i vartheta)/(lambda+1)) = C0
        eta, vth = -(lam + 1) * p0 / 2, (lam + 1) * thcl
        ok_h = ok_h and abs(cmath.cos((eta - I * vth) / (lam + 1)) - C1) < 1e-12 and abs(cmath.cos((eta + I * vth) / (lam + 1)) - C0) < 1e-12
        # with vartheta = lambda theta_d
        thd = vth / lam
        ok_h = ok_h and abs(lam * thd - (lam + 1) * thcl) < 1e-12
        # the weight argument: theta_cl = (1 - h) theta_d with h = 1/(lambda+1)
        h = 1 / (lam + 1)
        ok_h = ok_h and abs((1 - h) * thd - thcl) < 1e-12
        # the soliton parameter: y_d = sigma^{-(lambda+1)} = sigma^{-8 pi/beta^2}
        ok_h = ok_h and abs(math.exp(lam * thd) - math.exp(-thcl) ** (-(lam + 1))) < 1e-9 * math.exp(lam * thd)
# weights: e^{i Sigma/2} = e^{i beta phi_L} has h = beta^2/8 pi; the soliton exponentials (power 8 pi/beta^2) have h =
# 8 pi/beta^2 = lambda + 1, the spin of the non-local currents (charges of spin lambda)
for beta2 in (2.1, 5.3):
    lam = 8 * PI / beta2 - 1
    h = beta2 / (8 * PI)
    ok_h = ok_h and abs(h - 1 / (lam + 1)) < 1e-14 and abs((8 * PI / beta2) ** 2 * h - (lam + 1)) < 1e-12
# m0 = 0: the conditions Phi_x = Psi_t, Psi_x = Phi_t at x = 0, with Phi = phiL(t + x) + phiR(t - x) on x < 0 and
# Psi = psiL(t + x) + psiR(t - x) on x > 0, are psiL' = phiL', psiR' = -phiR' (sympy)
tt = sp.symbols('tt', real=True)
fL, fR, gL, gR = [sp.Function(nm) for nm in ('fL', 'fR', 'gL', 'gR')]
PhiS = fL(tt + x_s) + fR(tt - x_s)
PsiS = gL(tt + x_s) + gR(tt - x_s)
c1 = (sp.diff(PhiS, x_s) - sp.diff(PsiS, tt)).subs(x_s, 0).doit()
c2 = (sp.diff(PsiS, x_s) - sp.diff(PhiS, tt)).subs(x_s, 0).doit()
dfL, dfR, dgL, dgR = sp.symbols('dfL dfR dgL dgR')
c1s = c1.subs({sp.Derivative(fL(tt), tt): dfL, sp.Derivative(fR(tt), tt): dfR, sp.Derivative(gL(tt), tt): dgL,
               sp.Derivative(gR(tt), tt): dgR})
c2s = c2.subs({sp.Derivative(fL(tt), tt): dfL, sp.Derivative(fR(tt), tt): dfR, sp.Derivative(gL(tt), tt): dgL,
               sp.Derivative(gR(tt), tt): dgR})
sol = sp.solve([c1s, c2s], [dgL, dgR], dict=True)
ok_h = ok_h and len(sol) == 1 and sp.simplify(sol[0][dgL] - dfL) == 0 and sp.simplify(sol[0][dgR] + dfR) == 0
# sinh-Gordon: beta^2 -> -b^2, theta_d = theta_cl/(1 - h) = (1 - B/2) theta_cl with B = 2 b^2/(8 pi + b^2); expansion
b2 = sp.symbols('b2', positive=True)
hS = -b2 / (8 * sp.pi)
Bbook = 2 * b2 / (8 * sp.pi + b2)
ok_h = ok_h and sp.simplify(1 / (1 - hS) - (1 - Bbook / 2)) == 0
ok_h = ok_h and sp.simplify(sp.series(1 / (1 - hS), b2, 0, 2).removeO() - (1 - b2 / (8 * sp.pi))) == 0
# Lee-Yang (h = -1/5): Bajnok, Hollo and Watts' couplings mu = -(|h_c|/2) m^{6/5} e^{-i(b-2)pi/5}, mubar = its
# e^{+i(b-2) pi/5} partner are e^{+-(6/5) theta_b} times b-independent constants
hc = 0.685289
vals = []
for bb in (-3.0, -1.2, 0.5, 1.8, -3 + 2j):
    mu = -hc / 2 * cmath.exp(-I * (bb - 2) * PI / 5)
    mub = -hc / 2 * cmath.exp(I * (bb - 2) * PI / 5)
    vals.append((mu / cmath.exp(1.2 * thb(bb)), mub / cmath.exp(-1.2 * thb(bb))))
vals = np.array(vals)
ok_h = ok_h and np.abs(vals - vals[0]).max() < 1e-12
report("(h) what is left of hbar: the classical UV-IR limit maps the fused C0, C1 to eta = -(lambda+1) psi0/2, vartheta "
       "= (lambda+1) theta_cl, so lambda theta_d = (lambda+1) theta_cl and y_d = sigma^-(lambda+1); at m0 = 0 the defect "
       "passes left-movers and flips right-movers; h = beta^2/8pi = 1/(lambda+1), theta_cl = (1-h) theta_d; sinh-Gordon "
       "(1 - B/2) and Bajnok-Simon's 1 - b^2/8pi; Bajnok-Hollo-Watts couplings ~ e^{+-(6/5) theta_b}", ok_h)

# ----------------------------------------------------------------------------------------------------------------
# (i) The classical breather
# ----------------------------------------------------------------------------------------------------------------
ok_i = True


def TBn(th, n, lam, thd):
    return np.prod([Tcl(th - thd + I * PI * (n + 1 - 2 * k) / (2 * lam)) for k in range(1, n + 1)])


lam_i, thd_i = 6.3, 0.2
for n in (1, 2, 3, 5):
    for th in (0.4, -0.7 + 0.1j):
        ok_i = ok_i and abs(TBn(th, n + 1, lam_i, thd_i) - TBn(th + I * PI / (2 * lam_i), n, lam_i, thd_i)
                            * TBn(th - I * PI * n / (2 * lam_i), 1, lam_i, thd_i)) < 1e-12
        ok_i = ok_i and abs(TBn(th, n, lam_i, thd_i) * TBn(th + I * PI, n, lam_i, thd_i) - 1) < 1e-12
    ok_i = ok_i and abs(abs(TBn(0.37, n, lam_i, thd_i)) - 1) < 1e-12


def lnz_branch(th, xs, thd):
    """a continuous branch of ln z(theta - theta_d + i x) along xs"""
    vals = np.log(zcl(th - thd + I * xs))
    return np.real(vals) + I * np.unwrap(np.imag(vals))


alpha = 0.5
th_i = 0.45
errs1, errs2, errs3 = [], [], []
xs = np.linspace(alpha, PI - alpha, 20001)
lz = lnz_branch(th_i, xs, thd_i)
integral = np.trapz(lz, xs)
zt = zcl(th_i - thd_i + I * alpha)
delay_cl = math.log(abs(zt)) / (math.cos(alpha) * math.sinh(th_i))  # m0 = 1
for n in (6, 12, 24, 48):
    lam = n * PI / (PI - 2 * alpha)  # so that n pi/(2 lam) = pi/2 - alpha
    ks = np.arange(1, n + 1)
    pts = PI / 2 + PI * (n + 1 - 2 * ks) / (2 * lam)
    lnT = np.sum(lnz_branch(th_i, np.sort(pts), thd_i))
    errs1.append(abs(PI / lam * lnT - integral))
    ratio = TBn(th_i, n + 1, lam, thd_i) / TBn(th_i, n, lam, thd_i)
    errs2.append(abs(ratio - zt / abs(zt)))
    h_ = 1e-5
    dlt = (np.angle(TBn(th_i + h_, n, lam, thd_i) / TBn(th_i - h_, n, lam, thd_i))) / (2 * h_)
    M = (lam + 1) / PI
    mn = 2 * M * math.sin(n * PI / (2 * lam))
    errs3.append(abs(dlt / (mn * math.sinh(th_i)) - delay_cl))
errs1, errs2, errs3 = map(np.array, (errs1, errs2, errs3))
ok_i = ok_i and errs1[-1] < 1e-3 and np.all(errs1[1:] < errs1[:-1] / 3)
ok_i = ok_i and errs2[-1] < 0.05 and np.all(errs2[1:] < errs2[:-1])
ok_i = ok_i and errs3[-1] < 0.05 * abs(delay_cl) and np.all(errs3[1:] < errs3[:-1])
# the classical breather through the classical defect (mpmath)
mp.mp.dps = 30
m0n, thd_n, thb_n, al_n = mp.mpf(1), mp.mpf('0.3'), mp.mpf('0.55'), mp.mpf('0.7')
sign = mp.e ** (-thd_n)


def breather(x, t, c, z1=1, z2=1):
    r1, r2 = thb_n + 1j * al_n, thb_n - 1j * al_n
    Y1 = z1 * mp.e ** (m0n * (x * mp.cosh(r1) - t * mp.sinh(r1)) + c)
    Y2 = z2 * mp.e ** (m0n * (x * mp.cosh(r2) - t * mp.sinh(r2)) + mp.conj(c))
    return 4 * mp.atan((Y1 + Y2) / (1 + mp.tan(al_n) ** 2 * Y1 * Y2))


cb = mp.mpc('0.2', '0.4')
z1 = -mp.coth((thb_n + 1j * al_n - thd_n) / 2)
z2 = -mp.coth((thb_n - 1j * al_n - thd_n) / 2)


def defect_resid(zz1, zz2):
    worst = mp.mpf(0)
    for t in (mp.mpf('-1.3'), mp.mpf('0.2'), mp.mpf('1.9')):
        Ph = lambda x, tt: breather(x, tt, cb)
        Ps = lambda x, tt: breather(x, tt, cb, zz1, zz2)
        dpP = mp.diff(lambda tt: Ph(0, tt), t) + mp.diff(lambda x: Ph(x, t), 0)
        dpS = mp.diff(lambda tt: Ps(0, tt), t) + mp.diff(lambda x: Ps(x, t), 0)
        dmP = mp.diff(lambda tt: Ph(0, tt), t) - mp.diff(lambda x: Ph(x, t), 0)
        dmS = mp.diff(lambda tt: Ps(0, tt), t) - mp.diff(lambda x: Ps(x, t), 0)
        S_, D_ = Ph(0, t) + Ps(0, t), Ph(0, t) - Ps(0, t)
        r1 = (dpP - dpS) + 2 * sign * m0n * mp.sin(S_ / 2)
        r2 = (dmP + dmS) - 2 * m0n / sign * mp.sin(D_ / 2)
        worst = max(worst, abs(r1), abs(r2))
    return worst


r_good = defect_resid(z1, z2)
r_bad = defect_resid(mp.conj(z1), mp.conj(z2))
ok_i = ok_i and r_good < mp.mpf(10) ** -20 and r_bad > mp.mpf(10) ** -3
ok_i = ok_i and abs(mp.im(breather(mp.mpf('0.3'), mp.mpf('0.1'), cb, z1, z2))) < mp.mpf(10) ** -25
# the envelope: |Y1| = 1 at x cos(alpha) cosh(theta) m0 = t ... + const; on the right shifted by -ln|z1|/(m0 cos alpha cosh theta)
shift = -mp.log(abs(z1)) / (m0n * mp.cos(al_n) * mp.cosh(thb_n))
r1c = thb_n + 1j * al_n
x_left = (-mp.re(cb)) / (m0n * mp.re(mp.cosh(r1c)))
x_right = (-mp.re(cb) - mp.log(abs(z1))) / (m0n * mp.re(mp.cosh(r1c)))
ok_i = ok_i and abs((x_right - x_left) - shift) < mp.mpf(10) ** -25
mp.mp.dps = 25
report(f"(i) the classical breather: T_Bn satisfies the bootstrap with B_1, unitarity and crossing; at fixed alpha, "
       f"(pi/lambda) ln T_Bn -> int ln z (errors {errs1[0]:.0e} -> {errs1[-1]:.0e}), T_B(n+1)/T_Bn -> z/|z| "
       f"({errs2[-1]:.0e}), the Wigner delay -> ln|z(u + i alpha)|/(m0 cos alpha sinh theta) ({errs3[-1] / abs(delay_cl):.0e} "
       f"relative); the classical breather crosses with Y -> z(u -+ i alpha) Y (residual {float(r_good):.0e}; conjugate "
       f"factors fail) and its envelope is displaced by -ln|z|/(m0 cos alpha cosh theta)", ok_i)

# ----------------------------------------------------------------------------------------------------------------
# (j) Folding the diagram
# ----------------------------------------------------------------------------------------------------------------
ok_j = True


def Tk(th, sig, h, k):
    g = PI * k / h
    return (np.exp(-th) + sig * np.exp(I * g)) / (np.exp(-th) + sig * np.exp(-I * g))


def zk(u, h, k):
    g = PI * k / h
    return (np.exp(-u) - I * np.exp(I * g)) / (np.exp(-u) - I * np.exp(-I * g))


for h in (3, 4, 5, 6):
    pi_m = np.zeros((h, h))
    for kk in range(h):
        pi_m[kk, (kk + 1) % h] = 1  # (pi v)_k = v_{k+1}
    rho = np.zeros((h, h))
    for kk in range(h):
        rho[kk, (-kk) % h] = 1  # (rho v)_k = v_{-k}: rho alpha_j = alpha_{-j} up to the labelling of nodes
    ok_j = ok_j and np.abs(rho @ pi_m @ rho - np.linalg.inv(pi_m)).max() < 1e-14
    om = cmath.exp(2 * PI * I / h)
    for k in range(1, h):
        vk = np.array([om ** (k * kk) for kk in range(h)])  # pi v_k = om^k v_k
        ok_j = ok_j and np.abs(pi_m @ vk - om ** k * vk).max() < 1e-12
        rv = rho @ vk
        ok_j = ok_j and np.abs(pi_m @ rv - om ** (h - k) * rv).max() < 1e-12  # rho v_k is the mode of h - k
    # pairs
    s1 = 0.7
    ths = np.array([0.3, -0.9, 1.4, 0.2 + 0.3j])
    for k in range(1, h):
        res_plus = np.abs(Tk(ths, s1, h, k) * Tk(ths, s1, h, h - k) - Tk(ths, s1, h, h - k) * Tk(ths, s1, h, k)).max()
        res_minus = np.abs(Tk(ths, s1, h, k) * Tk(ths, -s1, h, h - k) - Tk(ths, s1, h, h - k) * Tk(ths, -s1, h, k)).max()
        trans = np.abs(Tk(ths, s1, h, k) * Tk(ths, -s1, h, h - k) - 1).max()
        ok_j = ok_j and res_plus < 1e-12 and res_minus < 1e-12 and trans < 1e-12
        if 2 * k != h:
            for s2 in (0.5, 0.9, 1.3, -0.4):
                ok_j = ok_j and np.abs(Tk(ths, s1, h, k) * Tk(ths, s2, h, h - k) - Tk(ths, s1, h, h - k) * Tk(ths, s2, h, k)).max() > 1e-4
        thd = -math.log(s1)
        ok_j = ok_j and np.abs(Tk(ths, s1, h, k) * Tk(ths, s1, h, h - k) - 1 / fx(k / h, ths - thd)).max() < 1e-12
        us = np.array([0.3, -1.1, 0.8 + 0.2j])
        g = PI * k / h
        ok_j = ok_j and np.abs(zk(us, h, k) * zk(us, h, h - k) - (np.cosh(us) + math.sin(g)) / (np.cosh(us) - math.sin(g))).max() < 1e-12
        ok_j = ok_j and np.abs(zk(us, h, k) * zk(us, h, h - k) - (1 / fx(k / h, us - I * PI / 2))).max() < 1e-12
us = np.array([0.3, -1.1, 0.8 + 0.2j])
cz = (np.cosh(us / 2 - I * PI / 12) / np.sinh(us / 2 - I * PI / 12)) * (np.cosh(us / 2 + I * PI / 12) / np.sinh(us / 2 + I * PI / 12))
ok_j = ok_j and np.abs(zk(us, 3, 1) * zk(us, 3, 2) - cz).max() < 1e-12
# h = 2: the pair (theta_d +- tau) is the type-II product; at tau = 0 it is Tcl^2 = 1/f_{1/2}
for th in (0.3, -0.4 + 0.1j):
    ok_j = ok_j and abs(Tk(th, math.exp(-0.4), 2, 1) - Tcl(th - 0.4)) < 1e-12
    ok_j = ok_j and abs(Tcl(th - 0.4) ** 2 - 1 / fx(0.5, th - 0.4)) < 1e-12
report("(j) folding the diagram: rho pi rho = pi^-1 and rho maps mode k to h - k (h = 3..6); a defect at sigma1 and "
       "the mirror at sigma2 give equal factors to k and h - k exactly for sigma2 = +- sigma1; -sigma1 is transparent; "
       "+sigma1 gives 1/f_{k/h}(theta - theta_d), and folded solitons (cosh u + sin gamma)/(cosh u - sin gamma), for h = 3 "
       "Corrigan-Zambon's Tzitzeica delay; h = 2 gives Tcl^2", ok_j)

raise SystemExit(0 if OK else 1)
