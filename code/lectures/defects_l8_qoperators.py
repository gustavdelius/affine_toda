"""Checks for Lecture 8 of the defect lecture notes, "Defects, L-operators and Q-operators"
(lecture-notes/defects/08-q-operators.qmd).

Conventions of the lecture. Classical sine-Gordon theory as in Lectures 2-6: d_pm = d_t pm d_x, d_+ d_- Phi =
-m0^2 sin Phi; the type-I_- defect has B = -2 m0 [sigma (cos(Sigma/2) - 1) + (cos(Delta/2) - 1)/sigma], sigma =
e^{-theta_d}, Sigma = Phi + Psi, Delta = Phi - Psi, and the conditions Phi_x = Psi_t - B_Phi, Psi_x = Phi_t + B_Psi at
x = 0; Baecklund's B_a: d_+(Phi - Psi) = 2 a m0 sin(Sigma/2), d_-(Phi + Psi) = -(2 m0/a) sin(Delta/2). Quantum: the
conventions of Lecture 7, lambda = 8 pi/beta^2 - 1, q = -e^{-i pi lambda}, y = e^{lambda theta}, y_d = e^{lambda theta_d},
Y = y/y_d; the L-operator L = P T of Lecture 7, L(v+|n>) = q^-n (v+|n> - iY v-|n+1>), L(v-|n>) = q^{n+s}(-iY v+|n-1>
+ v-|n>), and R = P Rcheck on (auxiliary, site), R(v+ v-) = alpha v+ v- + beta v- v+, alpha = q(z^2-1)/(q^2 z^2-1),
beta = z(q^2-1)/(q^2 z^2-1), z = y_a/y_j. On N solitons of rapidities theta_j: t_omega(theta) = Tr_a[omega^{h_a}
R_aN ... R_a1], M(y_d) = L_N^-1 ... L_1^-1, D(y_d) = <0| M |0>, and Q(y_d) = prod_j (1 + Y_j^2) D(y_d). Written from
these definitions, independently of code/defect_transfer.py and code/defect_universal_r.py. One PASS or FAIL line per
check; exit status 1 on any FAIL.

  (a) Turning the defect (sympy): the light-cone conditions are Phi_x = Psi_t - B_Phi, Psi_x = Phi_t + B_Psi; continued
      to Euclidean time (d_t = i d_tau), rotated by +90 degrees so that Phi lies before the line, and continued back,
      they become Phi_t = Psi_x - i B_Phi, Psi_t = Phi_x + i B_Psi, which is B_a with a = -i sigma (rapidity
      theta_d - i pi/2) and not with a = +i sigma; F = int (Phi Psi_x - i B) dx satisfies dF/dPhi = Phi_t and
      -dF/dPsi = Psi_t on these conditions. Turned the other way, so that Psi lies before the line, they say that
      Phi = B_(-i sigma) Psi: in both cases the field after the line is B_(-i sigma) of the field before.
  (b) Two solitons under the turned defect (mpmath, 30 digits): with Phi the two-soliton solution of Lecture 2,
      Psi = the same with Y_k -> z_k Y_k, z_k = -coth((theta_k - theta_d + i pi/2)/2), solves B_{-i sigma} at t = 0
      (residual 1e-28); the opposite sign of z_k fails; |z_k| = 1; beta^2 E and the momentum integral of Psi equal those
      of Phi, and 8 m0 (cosh theta_1 + cosh theta_2).
  (c) The monodromy is preserved: in the gauge chi = e^{i Phi h/4} chi~ of Lecture 6, with a~_x = (i/4) d_+Phi h +
      (m0/4)(lam s + s_{-2 Phi}/lam) and L~ = 1 + (i sigma/lam) s_{-Sigma}, the transports satisfy T~^Phi L~(-l) =
      L~(l) T~^Psi on [-l, l] for the pair of (b) (relative 1e-9), so tr T^Phi = tr T^Psi; it fails with sigma -> -i sigma.
  (d) Charge covariance: S^-1 L S = L q^{-h} on V (x) W (S|n> = |n+1>, h on the incoming soliton);
      <n| M |n> = q^{mn} <0| M |0> in the sector of charge m; D(y; s) = q^{-s (N - m)/2} D(y; 0).
  (e) The trace per state: the averages (2K+1)^-1 sum_{|n|<=K} <n| omega^{h_W} M |n> tend to omega^s D for omega = pm
      q^{-m/2} and to 0, like 1/K, for other unimodular omega, while the partial sums of the plain trace do not converge;
      for banded X, Y with bounded entries the averages of XY and YX differ by O(1/K).
  (f) Commutation: for N = 2..7, random rapidities, complex theta_d and s = 0.37, D commutes with t_omega(theta) in every
      sector of charge m for omega = pm q^{-m/2}; it does not for omega = 1 (m != 0), nor for the product of the L
      without inverses, nor for the inverses in the reversed order; [D(y_d), D(y_d')] = 0.
  (g) Fusion (sympy, identically in n with Q = q^n): at y = i y_d the operator R_aj(y/y_j) L_jW^-1(y_j/y_d) on
      V_j (x) V_a (x) W maps V_j (x) S into itself, S spanned by e_n = v+|n> - q^{-2n-s-1} v-|n+1>; in the basis q^n e_n
      it acts on S as [(1 + q^2 Y^2)/(1 + Y^2)] L^-1(qY) with offset s + 1, and on the quotient, spanned by the classes
      f_n of v-|n>, in the basis q^-n f_n, as q^-1 L^-1(Y/q) with offset s - 1. At y = 1.1 i y_d, S is not invariant.
  (h) The TQ relation t_omega(theta_d + i pi/(2 lam)) D(y_d) = q^{-N/2} D(q y_d) + q^{-N/2} prod_j [(1 + q^2 Y_j^2)/(1 +
      Y_j^2)] D(y_d/q), omega = q^{-m/2}, for N = 1..7 in every sector, s = 0.37 and complex theta_d; it fails at
      theta = theta_d + 1.1 i pi/(2 lam) and with q^{pm 2} in place of q^{pm 1}.
  (i) Baxter's Q (N = 5, 6): the eigenvalues of Q(y_d) on the common eigenvectors of t_omega are polynomials in 1/y_d^2 of
      degree (N - |m|)/2; their zeros satisfy the Bethe equations prod_j (1 + z_j^2)/(1 + q^-2 z_j^2) Q(q y_k) =
      -Q(y_k/q), z_j = y_j/y_k; the eigenvalues of t_omega(theta) at arbitrary theta follow from Baxter's formula with
      y_d = -i y; on the lines Im theta_d = pi/2 and pi/2 + pi/(2 lam) the left factor is a phase.
  (j) Limits: as theta_d -> +inf, D -> q^{-s(N-m)/2}, a function of the charge alone; as theta_d -> -inf,
      prod_j (Y_j/i) D tends on the neutral sector to charge conjugation times phases, and to 0 on the charged sectors.

Runtime: about a minute.
"""
import cmath
import math
import warnings

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

warnings.filterwarnings('ignore')
OK = True


def report(name, ok):
    global OK
    OK = OK and bool(ok)
    print(("PASS " if ok else "FAIL ") + name)


# ----------------------------------------------------------------------------------------------------------------
# (a) Turning the defect
# ----------------------------------------------------------------------------------------------------------------
m0s, sig = sp.symbols('m0 sigma', positive=True)
Ph, Ps, a_ = sp.symbols('Phi Psi a')
Pt, Px, St, Sx = sp.symbols('Phi_t Phi_x Psi_t Psi_x')
Bdef = -2 * m0s * (sig * (sp.cos((Ph + Ps) / 2) - 1) + (sp.cos((Ph - Ps) / 2) - 1) / sig)
BP, BS = sp.diff(Bdef, Ph), sp.diff(Bdef, Ps)
lc1 = (Pt + Px) - (St + Sx) + 2 * sig * m0s * sp.sin((Ph + Ps) / 2)
lc2 = (Pt - Px) + (St - Sx) - 2 * m0s / sig * sp.sin((Ph - Ps) / 2)
sol = sp.solve([lc1, lc2], [Px, Sx], dict=True)[0]
ok_a = sp.simplify(sol[Px] - (St - BP)) == 0 and sp.simplify(sol[Sx] - (Pt + BS)) == 0
# Euclidean continuation d_t = i d_tau; rotation (x, tau) -> (x', tau') = (-tau, x), so d_x = d_tau', d_tau = -d_x';
# back to real time tau' = i t', d_tau' = -i d_t'. Hence d_x -> -i d_t' and d_t = i d_tau -> -i d_x'.
# Phi_x = Psi_t - B_Phi  ->  -i Phi_t' = -i Psi_x' - B_Phi;  Psi_x = Phi_t + B_Psi  ->  -i Psi_t' = -i Phi_x' + B_Psi
tp, xp = sp.symbols("Phi_t' Psi_x'")
r1 = sp.solve(-sp.I * Pt - (-sp.I * Sx - BP), Pt)[0]
r2 = sp.solve(-sp.I * St - (-sp.I * Px + BS), St)[0]
ok_a = ok_a and sp.simplify(r1 - (Sx - sp.I * BP)) == 0 and sp.simplify(r2 - (Px + sp.I * BS)) == 0
bt1 = (Pt + Px) - (St + Sx) - 2 * a_ * m0s * sp.sin((Ph + Ps) / 2)
bt2 = (Pt - Px) + (St - Sx) + 2 * m0s / a_ * sp.sin((Ph - Ps) / 2)
solbt = sp.solve([bt1, bt2], [Pt, St], dict=True)[0]


def bt_mismatch(aval):
    e1 = sp.simplify(((solbt[Pt] - r1).subs(a_, aval)).rewrite(sp.exp))
    e2 = sp.simplify(((solbt[St] - r2).subs(a_, aval)).rewrite(sp.exp))
    return e1, e2


ok_a = ok_a and bt_mismatch(-sp.I * sig) == (0, 0) and bt_mismatch(sp.I * sig) != (0, 0)
Btil = sp.I * Bdef
ok_a = ok_a and sp.simplify(Sx - sp.diff(Btil, Ph) - r1) == 0 and sp.simplify(Px + sp.diff(Btil, Ps) - r2) == 0
# turned the other way, (x, tau) -> (tau, -x), Psi lies before the line: d_x -> -d_tau' -> i d_t', d_t -> i d_x'
r1o = sp.solve(sp.I * Pt - (sp.I * Sx - BP), Pt)[0]  # Phi_x = Psi_t - B_Phi
r2o = sp.solve(sp.I * St - (sp.I * Px + BS), St)[0]  # Psi_x = Phi_t + B_Psi
# B_a from Psi (before) to Phi (after): d_+(Psi - Phi) = 2 a m0 sin(Sigma/2), d_-(Psi + Phi) = -(2 m0/a) sin((Psi - Phi)/2)
bto = sp.solve([(St + Sx) - (Pt + Px) - 2 * a_ * m0s * sp.sin((Ph + Ps) / 2),
                (St - Sx) + (Pt - Px) + 2 * m0s / a_ * sp.sin((Ps - Ph) / 2)], [Pt, St], dict=True)[0]
ok_a = ok_a and all(sp.simplify(((bto[k] - r).subs(a_, -sp.I * sig)).rewrite(sp.exp)) == 0
                    for k, r in ((Pt, r1o), (St, r2o)))
report("(a) the defect conditions turned through 90 degrees (Phi before the line) are Phi_t = Psi_x - i B_Phi, "
       "Psi_t = Phi_x + i B_Psi: Baecklund's B_a with a = -i sigma, rapidity theta_d - i pi/2 (not a = +i sigma); "
       "F = int (Phi Psi_x - i B) generates them; turned the other way, Phi = B_(-i sigma) Psi", ok_a)

# ----------------------------------------------------------------------------------------------------------------
# (b) Two solitons under the turned defect
# ----------------------------------------------------------------------------------------------------------------
mp.mp.dps = 30
m0 = mp.mpf(1)
th1, th2, c1, c2 = mp.mpf('0.7'), mp.mpf('-0.4'), mp.mpf('0.3'), mp.mpf('-1.1')
thd = mp.mpf('0.25')
sg = mp.e ** (-thd)
ap = -1j * sg
thp = thd - 1j * mp.pi / 2


def W2(x, t, z1=1, z2=1):
    Y1 = z1 * mp.e ** (m0 * (x * mp.cosh(th1) - t * mp.sinh(th1)) + c1)
    Y2 = z2 * mp.e ** (m0 * (x * mp.cosh(th2) - t * mp.sinh(th2)) + c2)
    return -mp.coth((th1 - th2) / 2) * (Y1 - Y2) / (1 + Y1 * Y2)


def bt_resid(z1, z2):
    Phi = lambda x, t: 4 * mp.atan(W2(x, t))  # noqa: E731
    Psi = lambda x, t: 4 * mp.atan(W2(x, t, z1, z2))  # noqa: E731
    res = 0
    for x0 in (-1.3, 0.2, 0.9):
        x0 = mp.mpf(x0)
        d = lambda f, nx, nt: mp.diff(f, (x0, mp.mpf(0)), (nx, nt))  # noqa: E731
        Pp, Pm = d(Phi, 0, 1) + d(Phi, 1, 0), d(Phi, 0, 1) - d(Phi, 1, 0)
        Sp, Sm = d(Psi, 0, 1) + d(Psi, 1, 0), d(Psi, 0, 1) - d(Psi, 1, 0)
        P0, S0 = Phi(x0, 0), Psi(x0, 0)
        res = max(res, abs((Pp - Sp) - 2 * ap * m0 * mp.sin((P0 + S0) / 2)),
                  abs((Pm + Sm) + 2 * m0 / ap * mp.sin((P0 - S0) / 2)))
    return res


zk = [-mp.coth((th - thp) / 2) for th in (th1, th2)]
res_good = bt_resid(*zk)
res_bad = bt_resid(-zk[0], -zk[1])


def charges(z1, z2):
    def dens(x):
        W = W2(x, 0, z1, z2)
        Wx = mp.diff(lambda xx: W2(xx, 0, z1, z2), x)
        Wt = mp.diff(lambda tt: W2(x, tt, z1, z2), 0)
        Fx, Ft = 4 * Wx / (1 + W ** 2), 4 * Wt / (1 + W ** 2)
        return Ft ** 2 / 2 + Fx ** 2 / 2 + m0 ** 2 * 8 * W ** 2 / (1 + W ** 2) ** 2, -Ft * Fx
    pts = [-mp.inf, -3, 0, 3, mp.inf]
    return mp.quad(lambda x: dens(x)[0], pts), mp.quad(lambda x: dens(x)[1], pts)


E0, P0 = charges(1, 1)
E1, P1 = charges(*zk)
ok_b = (res_good < 1e-25 and res_bad > 1e-2 and all(abs(abs(z) - 1) < 1e-25 for z in zk)
        and abs(E1 - E0) < 1e-20 and abs(P1 - P0) < 1e-20 and abs(E0 - 8 * m0 * (mp.cosh(th1) + mp.cosh(th2))) < 1e-20
        and abs(P0 - 8 * m0 * (mp.sinh(th1) + mp.sinh(th2))) < 1e-20)
report(f"(b) two solitons under the turned defect: Y_k -> z_k Y_k with z_k = -coth((theta_k - theta_d + i pi/2)/2), "
       f"|z_k| = 1, solves B_(-i sigma) at t = 0 (residual {mp.nstr(res_good, 2)}; other sign {mp.nstr(res_bad, 2)}); "
       f"E and P unchanged ({mp.nstr(abs(E1 - E0), 2)}, {mp.nstr(abs(P1 - P0), 2)})", ok_b)

# ----------------------------------------------------------------------------------------------------------------
# (c) The monodromy is preserved (gauge e^{i Phi h/4}; all phases rational in W = tan(Phi/4))
# ----------------------------------------------------------------------------------------------------------------
hM = np.diag([1, -1]).astype(complex)
sM = np.array([[0, 1], [1, 0]], dtype=complex)
eM, fM = np.array([[0, 1], [0, 0]], dtype=complex), np.array([[0, 0], [1, 0]], dtype=complex)
th1f, th2f, c1f, c2f, thdf = 0.7, -0.4, 0.3, -1.1, 0.25
z1f, z2f = [complex(z) for z in zk]


def fields(x, z1, z2):
    """e^{i Phi/2} and d_+ Phi at (x, 0), with tan(Phi/4) = W, derivatives exact."""
    Cc = -1 / math.tanh((th1f - th2f) / 2)
    Y1 = z1 * cmath.exp(x * math.cosh(th1f) + c1f)
    Y2 = z2 * cmath.exp(x * math.cosh(th2f) + c2f)
    W = Cc * (Y1 - Y2) / (1 + Y1 * Y2)

    def dW(d1, d2):
        return Cc * ((d1 - d2) * (1 + Y1 * Y2) - (Y1 - Y2) * (d1 * Y2 + Y1 * d2)) / (1 + Y1 * Y2) ** 2
    Wx = dW(Y1 * math.cosh(th1f), Y2 * math.cosh(th2f))
    Wt = dW(-Y1 * math.sinh(th1f), -Y2 * math.sinh(th2f))
    e_half = (1 + 1j * W) / (1 - 1j * W)  # e^{i Phi/2}
    dplus = 4 * (Wt + Wx) / (1 + W ** 2)
    return e_half, dplus


def atil(x, lam_, z1, z2):
    e_half, dplus = fields(x, z1, z2)
    s_m2 = e_half ** -2 * eM + e_half ** 2 * fM  # s_{-2 Phi}
    return 0.25j * dplus * hM + 0.25 * (lam_ * sM + s_m2 / lam_)


def transport(lam_, z1, z2, ell):
    def rhs(x, v):
        T = v.reshape(2, 2)
        return (-atil(x, lam_, z1, z2) @ T).reshape(-1)
    sol_ = solve_ivp(rhs, (-ell, ell), np.eye(2, dtype=complex).reshape(-1), rtol=1e-12, atol=1e-13, method='DOP853')
    return sol_.y[:, -1].reshape(2, 2)


def Ltil(x, lam_, sig_):
    eP, _ = fields(x, 1, 1)
    eS, _ = fields(x, z1f, z2f)
    eSig = eP * eS  # e^{i Sigma/2}
    return np.eye(2) + sig_ / lam_ * (eSig ** -1 * eM + eSig * fM)  # s_{-Sigma}


ell, lam_spec = 6.0, 0.8 + 0.5j
TP_ = transport(lam_spec, 1, 1, ell)
TS_ = transport(lam_spec, z1f, z2f, ell)
sig_f = 1j * math.exp(-thdf)
rel_good = np.abs(TP_ @ Ltil(-ell, lam_spec, sig_f) - Ltil(ell, lam_spec, sig_f) @ TS_).max() / np.abs(TP_).max()
rel_bad = np.abs(TP_ @ Ltil(-ell, lam_spec, -sig_f) - Ltil(ell, lam_spec, -sig_f) @ TS_).max() / np.abs(TP_).max()
ok_c = rel_good < 1e-8 and rel_bad > 1e-3
report(f"(c) the turned defect preserves the monodromy: T~^Phi L~(-l) = L~(l) T~^Psi with sigma' = i sigma "
       f"(relative {rel_good:.0e}); with -i sigma {rel_bad:.1e}", ok_c)

# ----------------------------------------------------------------------------------------------------------------
# quantum chain
# ----------------------------------------------------------------------------------------------------------------
lam = math.sqrt(5) + 0.3
q = -cmath.exp(-1j * math.pi * lam)


def site_Linv(q_, s, Y, K):
    """A[b, a] (d x d) for L^-1(Y) on V (x) W, window n = -K..K; b, a: 0 = v+, 1 = v-."""
    ns = np.arange(-K, K + 1)
    d = 2 * K + 1
    A = np.zeros((2, 2, d, d), dtype=complex)
    den = 1 + Y ** 2
    A[0, 0] = np.diag(q_ ** ns.astype(complex)) / den
    A[1, 1] = np.diag(q_ ** (-ns - s).astype(complex)) / den
    for i, n in enumerate(ns):
        if i + 1 < d:  # v+|n> -> v-|n+1>
            A[1, 0][i + 1, i] = 1j * Y * q_ ** (-n - 1 - s) / den
        if i >= 1:  # v-|n> -> v+|n-1>
            A[0, 1][i - 1, i] = 1j * Y * q_ ** (n - 1) / den
    return A


def site_L(q_, s, Y, K):
    ns = np.arange(-K, K + 1)
    d = 2 * K + 1
    A = np.zeros((2, 2, d, d), dtype=complex)
    A[0, 0] = np.diag(q_ ** (-ns).astype(complex))
    A[1, 1] = np.diag(q_ ** (ns + s).astype(complex))
    for i, n in enumerate(ns):
        if i + 1 < d:
            A[1, 0][i + 1, i] = -1j * Y * q_ ** (-n)
        if i >= 1:
            A[0, 1][i - 1, i] = -1j * Y * q_ ** (n + s)
    return A


def mpo(sites, d):
    """sites[j] = A_j[b, a] (d x d), applied in the order j = 1..N. Returns X[b, a] = (A_N ... A_1)[b, a], slot 1 the
    most significant bit."""
    X = np.eye(d, dtype=complex)[None, None]
    for A in sites:
        nb = X.shape[0]
        Xn = np.zeros((2 * nb, 2 * nb, d, d), dtype=complex)
        for b2 in (0, 1):
            for a2 in (0, 1):
                if A[b2, a2].any():
                    Xn[b2::2, a2::2] = np.einsum('mn,bank->bamk', A[b2, a2], X)
        X = Xn
    return X


def M_full(q_, s, ys, yd, K, inverse=True, reverse=False):
    f = site_Linv if inverse else site_L
    sites = [f(q_, s, yj / yd, K) for yj in ys]
    if reverse:  # L_1^-1 ... L_N^-1: apply site N first; build with reversed slots, then permute back
        N = len(ys)
        X = mpo(sites[::-1], 2 * K + 1)
        perm = [int(format(v, f'0{N}b')[::-1], 2) for v in range(2 ** N)]
        return X[np.ix_(perm, perm)]
    return mpo(sites, 2 * K + 1)


def D_op(q_, s, ys, yd, **kw):
    K = len(ys) + 1
    return M_full(q_, s, ys, yd, K, **kw)[:, :, K, K]


def R_site(q_, y, yj):
    z = y / yj
    al = q_ * (z ** 2 - 1) / (q_ ** 2 * z ** 2 - 1)
    be = z * (q_ ** 2 - 1) / (q_ ** 2 * z ** 2 - 1)
    B = np.zeros((2, 2, 2, 2), dtype=complex)  # B[b_site, a_site][aux_out, aux_in]
    B[0, 0][0, 0] = 1
    B[1, 1][1, 1] = 1
    B[1, 1][0, 0] = al
    B[0, 0][1, 1] = al
    B[0, 1][1, 0] = be
    B[1, 0][0, 1] = be
    return B


def t_op(q_, ys, y, om):
    X = mpo([R_site(q_, y, yj) for yj in ys], 2)
    return om * X[:, :, 0, 0] + X[:, :, 1, 1] / om


def charge_list(N):
    return np.array([sum(1 if ((v >> (N - 1 - k)) & 1) == 0 else -1 for k in range(N)) for v in range(2 ** N)])


rng = np.random.default_rng(8)

# ----------------------------------------------------------------------------------------------------------------
# (d) Charge covariance
# ----------------------------------------------------------------------------------------------------------------
ok_d = True
K = 6
d = 2 * K + 1
s0, Yv = 0.37, 0.6 + 0.3j
A = site_L(q, s0, Yv, K)
Sop = np.diag(np.ones(d - 1), -1)  # |n> -> |n+1>
inner = slice(2, d - 2)
for b in (0, 1):
    for a in (0, 1):
        lhs = np.linalg.pinv(Sop) @ A[b, a] @ Sop  # S^-1 L S on the column of incoming v_a
        rhs = A[b, a] * q ** (-(1 if a == 0 else -1))
        ok_d = ok_d and np.abs((lhs - rhs)[inner, inner]).max() < 1e-12
N = 4
ys4 = np.exp(lam * rng.normal(size=N) * 0.6)
yd4 = cmath.exp(lam * (0.21 + 0.05j))
Kc = 2 * N + 2
Mf = M_full(q, s0, ys4, yd4, Kc)
ch4 = charge_list(N)
for n in (-2, -1, 1, 3):
    for v in range(2 ** N):
        for w in range(2 ** N):
            if ch4[v] == ch4[w]:
                lhs = Mf[v, w][Kc + n, Kc + n]
                rhs = q ** (ch4[v] * n) * Mf[v, w][Kc, Kc]
                ok_d = ok_d and abs(lhs - rhs) < 1e-12
D_s = D_op(q, s0, ys4, yd4)
D_0 = D_op(q, 0.0, ys4, yd4)
ok_d = ok_d and np.abs(D_s - np.diag(q ** (-s0 * (N - ch4) / 2)) @ D_0).max() < 1e-12
report("(d) charge covariance S^-1 L S = L q^-h; <n|M|n> = q^{mn} <0|M|0> in each sector; "
       "D(y; s) = q^{-s(N-m)/2} D(y; 0)", ok_d)

# ----------------------------------------------------------------------------------------------------------------
# (e) The trace per state
# ----------------------------------------------------------------------------------------------------------------
ok_e = True
sec = np.where(ch4 == 2)[0]
mS = 2
Kw = 46
Mw = M_full(q, s0, ys4, yd4, Kw)  # window n = -46..46; <n|M|n> is exact for |n| <= 42


def avg_trace(om, Kmax):
    tot = np.zeros((len(sec), len(sec)), dtype=complex)
    for n in range(-Kmax, Kmax + 1):
        blk = np.array([[Mw[v, w][Kw + n, Kw + n] for w in sec] for v in sec])
        tot += om ** (2 * n + s0) * blk
    return tot, tot / (2 * Kmax + 1)


Dsec = D_s[np.ix_(sec, sec)]
for om in (q ** (-mS / 2), -q ** (-mS / 2)):
    _, av = avg_trace(om, 40)
    ok_e = ok_e and np.abs(av - om ** s0 * Dsec).max() < 1e-12
om_bad = cmath.exp(0.77j)
phi_bad = cmath.phase(om_bad ** 2 * q ** mS)
bound = (2 / abs(1 - cmath.exp(1j * phi_bad)) + 1) * np.abs(Dsec).max()
part = [avg_trace(om_bad, Kk) for Kk in (10, 20, 40)]
devs = [np.abs(av).max() for _, av in part]
ok_e = ok_e and all(np.abs(tot).max() <= bound for tot, _ in part) and devs[2] < 0.05 * np.abs(Dsec).max()
# cyclicity of the average for banded matrices with bounded entries
Kb = 60
db = 2 * Kb + 1
Xb = sum(np.diag(np.exp(1j * rng.uniform(0, 6, db - abs(k))), k) for k in (-2, -1, 0, 1, 2))
Yb = sum(np.diag(np.exp(1j * rng.uniform(0, 6, db - abs(k))), k) for k in (-1, 0, 3))
cyc = [abs(np.trace((Xb @ Yb - Yb @ Xb)[Kb - Kk:Kb + Kk + 1, Kb - Kk:Kb + Kk + 1])) / (2 * Kk + 1) for Kk in (10, 40)]
ok_e = ok_e and cyc[1] < 0.5 * cyc[0] + 1e-12
report(f"(e) trace per state: averages -> omega^s D for omega = pm q^(-m/2), -> 0 otherwise, partial sums bounded "
       f"(|avg| at K = 10, 20, 40: {devs[0]:.3f}, {devs[1]:.3f}, {devs[2]:.3f}); the plain trace grows like K; "
       f"cyclic up to O(1/K) for banded matrices", ok_e)

# ----------------------------------------------------------------------------------------------------------------
# (f) Commutation
# ----------------------------------------------------------------------------------------------------------------
ok_f = True
worst, worst_bad, worst_noinv, worst_rev = 0.0, np.inf, np.inf, np.inf
for N in range(2, 8):
    ysN = np.exp(lam * rng.normal(size=N) * 0.6)
    ydN = cmath.exp(lam * (0.21 + 0.05j))
    DN = D_op(q, s0, ysN, ydN)
    chN = charge_list(N)
    yth = cmath.exp(lam * (0.4 + 0.3j))
    if N in (3, 4):
        Dni = D_op(q, s0, ysN, ydN, inverse=False)
        Drv = D_op(q, s0, ysN, ydN, reverse=True)
    for m in range(-N, N + 1, 2):
        idx = np.where(chN == m)[0]
        Ds = DN[np.ix_(idx, idx)]
        for om in (q ** (-m / 2), -q ** (-m / 2)):
            ts = t_op(q, ysN, yth, om)[np.ix_(idx, idx)]
            worst = max(worst, np.abs(Ds @ ts - ts @ Ds).max() / np.abs(Ds).max() / np.abs(ts).max())
        if len(idx) > 1:
            ts = t_op(q, ysN, yth, q ** (-m / 2))[np.ix_(idx, idx)]
            sc = np.abs(Ds).max() * np.abs(ts).max()
            if m != 0:
                tb = t_op(q, ysN, yth, 1.0)[np.ix_(idx, idx)]
                worst_bad = min(worst_bad, np.abs(Ds @ tb - tb @ Ds).max() / sc)
            if N in (3, 4):
                for Dx, which in ((Dni, 'ni'), (Drv, 'rv')):
                    Dxs = Dx[np.ix_(idx, idx)]
                    r = np.abs(Dxs @ ts - ts @ Dxs).max() / np.abs(Dxs).max() / np.abs(ts).max()
                    if which == 'ni':
                        worst_noinv = min(worst_noinv, r)
                    else:
                        worst_rev = min(worst_rev, r)
    if N == 5:
        DN2 = D_op(q, s0, ysN, cmath.exp(lam * (-0.5 + 0.2j)))
        ok_f = ok_f and np.abs(DN @ DN2 - DN2 @ DN).max() < 1e-12 * np.abs(DN).max() * np.abs(DN2).max()
ok_f = ok_f and worst < 1e-13 and worst_bad > 1e-3 and worst_noinv > 1e-3 and worst_rev > 1e-3
report(f"(f) [D, t_omega] = 0 for omega = pm q^(-m/2), N = 2..7, every sector (relative {worst:.0e}); "
       f"fails for omega = 1 ({worst_bad:.1e}), without inverses ({worst_noinv:.1e}), in reversed order "
       f"({worst_rev:.1e}); [D(y_d), D(y_d')] = 0", ok_f)

# ----------------------------------------------------------------------------------------------------------------
# (g) Fusion, exactly in n
# ----------------------------------------------------------------------------------------------------------------
qs, Ys, ss, Qn = sp.symbols('q Y s Q', nonzero=True)  # Q stands for q^n


def Linv_apply(vec, Y_):
    out = {}

    def add(key, val):
        out[key] = out.get(key, 0) + val
    for (bj, a, k), co in vec.items():
        nk = Qn * qs ** k
        if bj == 0:
            add((0, a, k), co * nk / (1 + Y_ ** 2))
            add((1, a, k + 1), co * sp.I * Y_ / (nk * qs * qs ** ss) / (1 + Y_ ** 2))
        else:
            add((0, a, k - 1), co * sp.I * Y_ * nk / qs / (1 + Y_ ** 2))
            add((1, a, k), co / (nk * qs ** ss) / (1 + Y_ ** 2))
    return out


def R_apply(vec, zz):
    al = qs * (zz ** 2 - 1) / (qs ** 2 * zz ** 2 - 1)
    be = zz * (qs ** 2 - 1) / (qs ** 2 * zz ** 2 - 1)
    out = {}

    def add(key, val):
        out[key] = out.get(key, 0) + val
    for (bj, a, k), co in vec.items():
        if a == bj:
            add((bj, a, k), co)
        else:
            add((bj, a, k), co * al)
            add((a, bj, k), co * be)
    return out


def Cn(k):
    return -qs ** (-ss - 1) / (Qn ** 2 * qs ** (2 * k))


def reduce_mod_S(vec):
    out = {}
    for (bj, a, k), co in vec.items():
        key = (bj, 1, k + 1) if a == 0 else (bj, 1, k)
        out[key] = out.get(key, 0) + (-Cn(k) * co if a == 0 else co)
    return {kk: sp.simplify(v) for kk, v in out.items() if sp.simplify(v) != 0}


def fused(vec, zfac=1):
    """R_aj(y/y_j) L_jW^-1(Y) with y/y_j = zfac * i/Y (zfac = 1: y = i y_d)."""
    return R_apply(Linv_apply(vec, Ys), zfac * sp.I / Ys)


ok_g = True
for bj in (0, 1):  # invariance of S
    ok_g = ok_g and reduce_mod_S(fused({(bj, 0, 0): 1, (bj, 1, 1): Cn(0)})) == {}
# block on S in the basis q^n e_n: compare with ((1 + q^2 Y^2)/(1 + Y^2)) L^-1(qY) at s + 1
for bj in (0, 1):
    out = fused({(bj, 0, 0): 1, (bj, 1, 1): Cn(0)})
    # coefficient of v_b (x) e_{n+k}: read off the v+ component |n+k>
    coeff = {}
    for (b2, a2, k2), co in out.items():
        if a2 == 0:
            coeff[(b2, k2)] = sp.simplify(co * qs ** (-k2))  # coefficient of q^{n+k} e_{n+k}
    pref = (1 + qs ** 2 * Ys ** 2) / (1 + Ys ** 2)
    Yp, sp1 = qs * Ys, ss + 1
    if bj == 0:
        tgt = {(0, 0): Qn / (1 + Yp ** 2), (1, 1): sp.I * Yp / (Qn * qs * qs ** sp1) / (1 + Yp ** 2)}
    else:
        tgt = {(0, -1): sp.I * Yp * Qn / qs / (1 + Yp ** 2), (1, 0): 1 / (Qn * qs ** sp1) / (1 + Yp ** 2)}
    ok_g = ok_g and set(coeff) == set(tgt) and all(sp.simplify(coeff[k_] - pref * tgt[k_]) == 0 for k_ in tgt)
# block on the quotient, basis q^-n f_n: compare with q^-1 L^-1(Y/q) at s - 1
for bj in (0, 1):
    red = reduce_mod_S(fused({(bj, 1, 0): 1}))
    coeff = {(b2, k2): sp.simplify(co * qs ** k2) for (b2, a2, k2), co in red.items()}  # of q^-(n+k) f_(n+k)
    Ypp, sm1 = Ys / qs, ss - 1
    if bj == 0:
        tgt = {(0, 0): Qn / (1 + Ypp ** 2), (1, 1): sp.I * Ypp / (Qn * qs * qs ** sm1) / (1 + Ypp ** 2)}
    else:
        tgt = {(0, -1): sp.I * Ypp * Qn / qs / (1 + Ypp ** 2), (1, 0): 1 / (Qn * qs ** sm1) / (1 + Ypp ** 2)}
    ok_g = ok_g and set(coeff) == set(tgt) and all(sp.simplify(coeff[k_] - tgt[k_] / qs) == 0 for k_ in tgt)
# away from the special point S is not invariant
leak = reduce_mod_S(fused({(0, 0, 0): 1, (0, 1, 1): Cn(0)}, zfac=sp.Rational(11, 10)))
ok_g = ok_g and leak != {}
report("(g) fusion at y = i y_d: S = span{v+|n> - q^(-2n-s-1) v-|n+1>} invariant; on S the block is "
       "(1+q^2 Y^2)/(1+Y^2) L^-1(qY) with s+1, on the quotient q^-1 L^-1(Y/q) with s-1 (sympy, all n); "
       "S not invariant at y = 1.1 i y_d", ok_g)

# ----------------------------------------------------------------------------------------------------------------
# (h) The TQ relation
# ----------------------------------------------------------------------------------------------------------------
ok_h = True
worst_tq, best_bad1, best_bad2 = 0.0, np.inf, np.inf
for N in range(1, 8):
    ysN = np.exp(lam * rng.normal(size=N) * 0.6)
    ydN = cmath.exp(lam * (0.13 - 0.04j))
    chN = charge_list(N)
    D0, Dp, Dm = (D_op(q, s0, ysN, y_) for y_ in (ydN, q * ydN, ydN / q))
    Dp2, Dm2 = D_op(q, s0, ysN, q ** 2 * ydN), D_op(q, s0, ysN, ydN / q ** 2)
    Y2 = (ysN / ydN) ** 2
    c2 = q ** (-N / 2) * np.prod((1 + q ** 2 * Y2) / (1 + Y2))
    for m in range(-N, N + 1, 2):
        idx = np.where(chN == m)[0]
        om = q ** (-m / 2)
        sub = lambda X: X[np.ix_(idx, idx)]  # noqa: E731
        tt = sub(t_op(q, ysN, 1j * ydN, om))
        lhs = tt @ sub(D0)
        rhs = q ** (-N / 2) * sub(Dp) + c2 * sub(Dm)
        worst_tq = max(worst_tq, np.abs(lhs - rhs).max() / np.abs(lhs).max())
        if N >= 3 and len(idx) > 1:
            tt_bad = sub(t_op(q, ysN, 1j * ydN * cmath.exp(0.1j * math.pi / 2), om))
            best_bad1 = min(best_bad1, np.abs(tt_bad @ sub(D0) - rhs).max() / np.abs(lhs).max())
            rhs2 = q ** (-N / 2) * sub(Dp2) + c2 * sub(Dm2)
            best_bad2 = min(best_bad2, np.abs(lhs - rhs2).max() / np.abs(lhs).max())
ok_h = worst_tq < 1e-12 and best_bad1 > 1e-3 and best_bad2 > 1e-3
report(f"(h) TQ relation for N = 1..7 in every sector (relative {worst_tq:.0e}); fails at theta_d + 1.1 i pi/(2 lam) "
       f"({best_bad1:.1e}) and with q^(pm2) shifts ({best_bad2:.1e})", ok_h)

# ----------------------------------------------------------------------------------------------------------------
# (i) Baxter's Q: zeros are Bethe roots
# ----------------------------------------------------------------------------------------------------------------
ok_i = True
bethe_worst, baxter_worst, phase_worst, deg_ok = 0.0, 0.0, 0.0, True
for N, ths in ((5, [-0.9, -0.4, 0.05, 0.5, 0.95]), (6, [-1.0, -0.55, -0.15, 0.2, 0.6, 1.05])):
    ysN = np.exp(lam * np.array(ths))
    chN = charge_list(N)
    for m in range(N % 2, N + 1, 2):
        idx = np.where(chN == m)[0]
        om = q ** (-m / 2)
        tg = t_op(q, ysN, cmath.exp(lam * (0.13 + 0.07j)), om)[np.ix_(idx, idx)]
        ev, Vv = np.linalg.eig(tg)
        Vi = np.linalg.inv(Vv)
        Mg = (N - abs(m)) // 2
        npts = Mg + 4
        w2s = 0.6 * np.exp(1j * (np.arange(npts) * 2 * np.pi / npts + 0.3))
        Qs = []
        for w2 in w2s:
            yd_ = 1 / np.sqrt(w2)
            Qm = np.prod(1 + (ysN / yd_) ** 2) * D_op(q, 0.0, ysN, yd_)[np.ix_(idx, idx)]
            dd = Vi @ Qm @ Vv
            ok_i = ok_i and np.abs(dd - np.diag(np.diag(dd))).max() < 1e-9 * np.abs(dd).max()
            Qs.append(np.diag(dd))
        Qs = np.array(Qs)
        for st in range(len(idx)):
            cf = np.linalg.solve(np.vander(w2s, npts, increasing=True), Qs[:, st])
            deg_ok = deg_ok and np.abs(cf[Mg + 1:]).max() < 1e-9 * np.abs(cf).max() and abs(cf[Mg]) > 1e-8 * np.abs(cf).max()
            cM = cf[:Mg + 1]
            Qf = lambda y_: np.polyval(cM[::-1], 1 / y_ ** 2)  # noqa: E731
            for w2k in (np.roots(cM[::-1]) if Mg > 0 else []):
                yk = 1 / np.sqrt(w2k)
                z2 = (ysN / yk) ** 2
                lhs = np.prod((1 + z2) / (1 + z2 / q ** 2)) * Qf(q * yk)
                bethe_worst = max(bethe_worst, abs(lhs + Qf(yk / q)) / abs(lhs))
            # Baxter's formula for the eigenvalue of t at theta, with y_d = -i y
            for thv in (0.31 - 0.2j, -0.7 + 0.1j):
                yv = cmath.exp(lam * thv)
                ydb = -1j * yv
                tval = (Vi @ t_op(q, ysN, yv, om)[np.ix_(idx, idx)] @ Vv)[st, st]
                z2 = (ysN / ydb) ** 2
                bax = q ** (-N / 2) * (np.prod((1 + z2) / (1 + z2 / q ** 2)) * Qf(q * ydb) + Qf(ydb / q)) / Qf(ydb)
                baxter_worst = max(baxter_worst, abs(tval - bax) / abs(tval))
    # on the line Im theta_d = pi/2 the factor (1 + z^2)/(1 + q^-2 z^2) is a phase
    for xr in (-0.8, 0.1, 0.9):
        for im in (0.5 * math.pi, 0.5 * math.pi + 0.5 * math.pi / lam):  # the line and the line half a period away
            yk = cmath.exp(lam * (xr + 1j * im))
            z2 = (ysN / yk) ** 2
            phase_worst = max(phase_worst, abs(abs(np.prod((1 + z2) / (1 + z2 / q ** 2))) - 1))
ok_i = ok_i and deg_ok and bethe_worst < 1e-10 and baxter_worst < 1e-10 and phase_worst < 1e-12
report(f"(i) Baxter's Q: eigenvalues of Q(y_d) are polynomials in 1/y_d^2 of degree (N-|m|)/2; their zeros solve the "
       f"Bethe equations (relative {bethe_worst:.0e}); Baxter's formula gives t at any theta ({baxter_worst:.0e}); "
       f"on Im theta_d = pi/2 and pi/2 + pi/(2 lam) the Bethe factor is a phase", ok_i)

# ----------------------------------------------------------------------------------------------------------------
# (j) Limits
# ----------------------------------------------------------------------------------------------------------------
ok_j = True
N, sj = 4, 0.3
ysj = np.exp(lam * np.array([-0.6, -0.1, 0.35, 0.8]))
chj = charge_list(N)
errs = []
for thv in (3.0, 4.5):
    Dj = D_op(q, sj, ysj, math.exp(lam * thv))
    errs.append(np.abs(Dj - np.diag(q ** (-sj * (N - chj) / 2))).max())
ok_j = ok_j and errs[1] < 1e-7 and errs[1] < 1e-2 * errs[0]
sec0 = np.where(chj == 0)[0]
oth = np.where(chj != 0)[0]
lim = []
for thv in (-3.0, -4.5):
    yd_ = math.exp(lam * thv)
    Dn = np.prod(ysj / yd_ / 1j) * D_op(q, sj, ysj, yd_)
    blk = Dn[np.ix_(sec0, sec0)]
    flip = [2 ** N - 1 - v for v in sec0]  # charge conjugation: flip every soliton
    Cperm = np.array([[1.0 if sec0[i] == flip[j] else 0.0 for j in range(len(sec0))] for i in range(len(sec0))])
    phases = np.diag(blk @ Cperm.T)
    lim.append((np.abs(blk - np.diag(phases) @ Cperm).max(), np.abs(np.abs(phases) - 1).max(),
                np.abs(Dn[np.ix_(oth, oth)]).max()))
ok_j = ok_j and lim[1][0] < 1e-7 and lim[1][1] < 1e-7 and lim[1][2] < 1e-7 and lim[1][2] < 1e-2 * lim[0][2]
report(f"(j) limits: theta_d -> +inf, D -> q^(-s(N-m)/2) (error {errs[1]:.0e}); theta_d -> -inf, prod(Y/i) D -> "
       f"charge conjugation times phases on m = 0, -> 0 on m != 0 ({lim[1][2]:.0e})", ok_j)

raise SystemExit(0 if OK else 1)
