"""Exact S-matrices of the Delius-Grisaru-Penati-Zanon superalgebra theories
(sec-super-scattering): the boson-fermion sector of A^(2)(0,2n-1) and the
fermion amplitudes of C^(2)(n+1).

Notation of the chapter (eq-dgpz-amplitudes): h = 2n - 1, gamma = -sin(pi Omega),
B = 2 Omega the block parameter of eq-blocks, and

    R = i gamma T / sinh(h theta/2),        D = (-1)^(n+1) i gamma T / cosh(h theta/2),
    B = (1 + 2 i gamma / sinh(h theta)) T,  F = (-1 + 2 i gamma / sinh(h theta)) T,
    T = Pi(theta) Pi(i pi - theta)          (Delius et al. 1991, (2.34)-(2.41)).

The two-particle operator Sc(theta) on V (x) V, V = span(phi_n even, psi odd),
maps the in-state basis vector a (x) b to the out-state basis, with the
particles' spatial order exchanged (the braid form): Sc(phi phi) = B phi phi +
D psi psi, Sc(psi psi) = -D phi phi + F psi psi, Sc(psi phi) = T phi psi + R psi phi,
Sc(phi psi) = T psi phi + R phi psi.  At Omega = 0 it is the graded permutation.

Checks:
 1. T has no poles in the physical strip, satisfies crossing T(i pi - theta) = T(theta)
    and the unitarity relation T(theta)T(-theta) = sinh^2(h theta/2)/(sinh^2(h theta/2) + gamma^2).
 2. The eight Yang-Baxter equations (2.7) of Delius et al. (and their B <-> F partners).
 3. The braid Yang-Baxter equation for Sc as an even operator on the graded space,
    and its failure for the ungraded assignment (free limit = ordinary permutation).
 4. Unitarity Sc(theta)Sc(-theta) = 1 and, for real theta, Sc Sc^dagger = 1.
 5. Crossing: B, F, T symmetric under theta -> i pi - theta, D(i pi - theta) = i R(theta);
    the factor i is the crossing phase of the Majorana spinors u(theta + i pi) = i v(theta)
    of Delius et al. (A.2).
 6. B - F = 2T and the free-fermion condition B F + T^2 - D^2 - R^2 = 0 (their (B.1)).
 7. Poles: B, F, D have poles at theta = i pi (2(n-a)-1)/h, at s = m_a^2 = 8 sin^2(pi a/h),
    with residues of B and F equal to i times a positive number and a rank-one residue
    matrix, as for a bound state in both channels; R, B, F have poles at theta = 2 pi i a/h.
 8. Tree level: the order-Omega terms of T, R, D, B agree with the tree amplitudes (2.5) of
    Delius et al. with Omega = beta^2/4pi, and the phi_n phi_n -> phi_n phi_n amplitude B
    agrees with the tree amplitude computed here from the bosonic potential of the
    A^(2)(0,2n-1) Lagrangian (n = 2, 3).  Delius et al. fix only gamma^2 = sin^2(pi Omega);
    the tree amplitudes require gamma = -sin(pi Omega).
 9. The bootstrap of the phi_a phi_n amplitude: (2.45) equals the block product (2.47) = (16).
10. C^(2)(n+1): the fermion amplitudes (3.10) (n even) and (3.13) (n odd) are unitary,
    satisfy the stated crossing relations, and bootstrap to (3.15)/(3.16).
"""
import itertools
import os
import sys
import numpy as np
import sympy as sp
import mpmath as mp

mp.mp.dps = 30

FAILS = []


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok:
        FAILS.append(name)


# ------------------------------------------------------------------ blocks
def paren(x, th, h):
    return np.sinh(th / 2 + 1j * np.pi * x / (2 * h)) / np.sinh(th / 2 - 1j * np.pi * x / (2 * h))


def brace(x, th, h, Bp):
    """Block {x} of eq-blocks with parameter B = 2 Omega."""
    return paren(x - 1, th, h) * paren(x + 1, th, h) / (paren(x - 1 + Bp, th, h) * paren(x + 1 - Bp, th, h))


# ------------------------------------------------------- the function T(theta)
def log_Pi(th, h, Om, N=40, K=14):
    """log Pi(theta) of Delius et al. (2.41), eta = -i h theta/(2 pi): the first N
    factors exactly, and the tail l > N from the Stirling series
    log Gamma(z + a) ~ (z + a - 1/2) log z - z + log(2 pi)/2
                       + sum_k (-1)^k B_k(a)/(k(k-1) z^(k-1)),   z = h l,
    summed over l with Hurwitz zeta functions.  The (z + a - 1/2) log z terms and the
    k = 2 terms cancel between numerator and denominator."""
    eta = -1j * h * mp.mpc(th) / (2 * mp.pi)
    Om = mp.mpf(Om)
    A = [eta - Om + 1 - h, eta + Om - h, eta + 1 - mp.mpf(h) / 2, eta - mp.mpf(h) / 2]
    Bv = [eta - Om + 1 - mp.mpf(h) / 2, eta + Om - mp.mpf(h) / 2, eta + 1 - h, eta - h]
    s = mp.fsum(mp.loggamma(h * l + x) for l in range(1, N + 1) for x in A) - \
        mp.fsum(mp.loggamma(h * l + x) for l in range(1, N + 1) for x in Bv)
    tail = 0
    for k in range(3, K + 1):
        c = mp.fsum(mp.bernpoly(k, x) for x in A) - mp.fsum(mp.bernpoly(k, x) for x in Bv)
        tail += (-1) ** k * c / (k * (k - 1)) * mp.mpf(h) ** (1 - k) * mp.zeta(k - 1, N + 1)
    return complex(s + tail)


def Tfun(th, h, Om):
    return np.exp(log_Pi(th, h, Om) + log_Pi(1j * np.pi - th, h, Om))


def amps(th, n, Om, T=None):
    h = 2 * n - 1
    g = -np.sin(np.pi * Om)          # sign fixed by the tree-level amplitudes (check 8)
    if T is None:
        T = Tfun(th, h, Om)
    R = 1j * g * T / np.sinh(h * th / 2)
    D = (-1) ** (n + 1) * 1j * g * T / np.cosh(h * th / 2)
    B = (1 + 2j * g / np.sinh(h * th)) * T
    F = (-1 + 2j * g / np.sinh(h * th)) * T
    return B, F, D, R, T


def Sc(th, n, Om, graded=True, T=None):
    """Matrix of the braid-form operator in the basis (phiphi, phipsi, psiphi, psipsi);
    column k is the image of the k-th in-basis vector."""
    B, F, D, R, T = amps(th, n, Om, T)
    if not graded:
        F = -F
    M = np.zeros((4, 4), dtype=complex)
    M[0, 0], M[3, 0] = B, D              # phi phi -> B phi phi + D psi psi
    M[0, 3], M[3, 3] = -D, F             # psi psi -> -D phi phi + F psi psi
    M[1, 2], M[2, 2] = T, R              # psi phi -> T phi psi + R psi phi
    M[2, 1], M[1, 1] = T, R              # phi psi -> T psi phi + R phi psi
    return M


def kron(A, Bm):
    return np.kron(A, Bm)


def mnorm(M):
    return np.max(np.abs(M))


I2 = np.eye(2)
Om0 = (np.sqrt(5) - 1) / 4          # generic, irrational
pts = [0.37 + 0.21j, -0.52 + 1.13j, 0.81 + 2.4j]

# ======================================================== 1. the function T
for n in (2, 3):
    h = 2 * n - 1
    ok_cross = all(abs(Tfun(1j * np.pi - th, h, Om0) / Tfun(th, h, Om0) - 1) < 1e-9 for th in pts)
    ok_unit = all(abs(Tfun(th, h, Om0) * Tfun(-th, h, Om0)
                      - np.sinh(h * th / 2) ** 2 / (np.sinh(h * th / 2) ** 2 + np.sin(np.pi * Om0) ** 2)) < 1e-9
                  for th in pts)
    # no poles in the strip: the Gamma arguments of the numerators of Pi(theta) and
    # Pi(i pi - theta) have positive real part for 0 < Re eta < h/2 (structural), and |T| is
    # bounded on a grid of the strip
    grid = [abs(Tfun(x + 1j * np.pi * y, h, Om0)) for x in (-1.5, 0.3) for y in (0.05, 0.31, 0.5, 0.77, 0.95)]
    report(f'n = {n}: T(theta) = Pi(theta)Pi(i pi - theta) obeys crossing and the unitarity relation (2.37), '
           f'and is bounded in the physical strip (max |T| on the grid {float(max(grid)):.3f})',
           ok_cross and ok_unit and max(grid) < 50)

# ================================================ 2.-6. algebraic identities
rng = np.random.default_rng(7)
for n in (2, 3):
    ok_yb, ok_braid, ok_ungr, ok_unit, ok_herm, ok_cross, ok_bf, ok_free, ok_b1 = [True] * 9
    for _ in range(4):
        t1 = complex(rng.uniform(-1, 1), rng.uniform(0.1, 1))
        t2 = complex(rng.uniform(-1, 1), rng.uniform(0.1, 1))
        # the common factor T cancels from the Yang-Baxter equations; use T = 1 there
        b1, f1, d1, r1, s1 = amps(t1, n, Om0, T=1)
        b2, f2, d2, r2, s2 = amps(t2, n, Om0, T=1)
        b3, f3, d3, r3, s3 = amps(t1 + t2, n, Om0, T=1)
        for (B1, F1, B2, F2, B3, F3) in ((b1, f1, b2, f2, b3, f3), (f1, b1, f2, b2, f3, b3)):
            eqs = [B1 * r3 * B2 - d1 * F3 * d2 - (r1 * B3 * r2 + s1 * r3 * s2),
                   B1 * s3 * r2 - d1 * d3 * s2 - (s1 * B3 * r2 + r1 * r3 * s2),
                   B1 * B3 * d2 + d1 * r3 * F2 - (r1 * d3 * B2 + s1 * s3 * d2),
                   B1 * d3 * s2 + d1 * s3 * r2 - (s1 * d3 * B2 + r1 * s3 * d2)]
            ok_yb &= max(abs(e) for e in eqs) < 1e-12
        L = kron(Sc(t1, n, Om0, T=1), I2) @ kron(I2, Sc(t1 + t2, n, Om0, T=1)) @ kron(Sc(t2, n, Om0, T=1), I2)
        Rr = kron(I2, Sc(t2, n, Om0, T=1)) @ kron(Sc(t1 + t2, n, Om0, T=1), I2) @ kron(I2, Sc(t1, n, Om0, T=1))
        ok_braid &= mnorm(L - Rr) < 1e-12
        L = kron(Sc(t1, n, Om0, False, 1), I2) @ kron(I2, Sc(t1 + t2, n, Om0, False, 1)) @ kron(Sc(t2, n, Om0, False, 1), I2)
        Rr = kron(I2, Sc(t2, n, Om0, False, 1)) @ kron(Sc(t1 + t2, n, Om0, False, 1), I2) @ kron(I2, Sc(t1, n, Om0, False, 1))
        ok_ungr &= mnorm(L - Rr) > 1e-3
    for th in pts[:2] + [0.63, -1.7]:
        ok_unit &= mnorm(Sc(th, n, Om0) @ Sc(-th, n, Om0) - np.eye(4)) < 1e-9
        B, F, D, R, T = amps(th, n, Om0)
        Bc, Fc, Dc, Rc, Tc = amps(1j * np.pi - th, n, Om0)
        ok_cross &= max(abs(Bc - B), abs(Fc - F), abs(Tc - T), abs(Dc - 1j * R)) < 1e-9
        ok_bf &= abs(B - F - 2 * T) < 1e-12
        ok_free &= abs(B * F + T ** 2 + D ** 2 - R ** 2) < 1e-10      # Fan-Wu: d1 d2 = D (-D)
        ok_b1 &= abs(B * F + T ** 2 - D ** 2 - R ** 2) > 1e-3       # (B.1) as printed
    for th in (0.63, -1.7, 2.2):
        S = Sc(th, n, Om0)
        ok_herm &= mnorm(S @ S.conj().T - np.eye(4)) < 1e-9
    report(f'n = {n}: the Yang-Baxter equations (2.7) of Delius et al. and their B <-> F partners', ok_yb)
    report(f'n = {n}: braid Yang-Baxter equation for the graded operator Sc (eq-dgpz-ybe)', ok_braid)
    report(f'n = {n}: with the ungraded assignment (free limit = ordinary permutation) the braid '
           f'Yang-Baxter equation fails', ok_ungr)
    report(f'n = {n}: unitarity Sc(theta)Sc(-theta) = 1, and Sc Sc^dagger = 1 for real theta', ok_unit and ok_herm)
    report(f'n = {n}: crossing B, F, T symmetric and D(i pi - theta) = i R(theta)', ok_cross)
    report(f'n = {n}: B - F = 2T and the free-fermion condition a1 a2 + b1 b2 = c1 c2 + d1 d2 of the '
           f'eight-vertex weights of Sc, i.e. BF + T^2 + D^2 - R^2 = 0 since d1 d2 = D(-D); the sign '
           f'printed in (B.1) of Delius et al., BF + T^2 - D^2 - R^2 = 0, fails', ok_bf and ok_free and ok_b1)

# Majorana spinors of Delius et al. (A.2): u = (e^{theta/2}, e^{-theta/2})/sqrt2, v = (e^{theta/2}, -e^{-theta/2})/sqrt2
u = lambda th: np.array([np.exp(th / 2), np.exp(-th / 2)]) / np.sqrt(2)
v = lambda th: np.array([np.exp(th / 2), -np.exp(-th / 2)]) / np.sqrt(2)
report('crossing phase of the Majorana spinors: u(theta + i pi) = i v(theta)',
       all(mnorm(u(th + 1j * np.pi) - 1j * v(th)) < 1e-14 for th in (0.3, -1.1)))

# ================================================================= 7. poles
for n in (2, 3):
    h = 2 * n - 1
    ok = True
    for Om in (0.1, Om0, 0.6):
        for a in range(1, n):
            ts = 1j * np.pi * (2 * (n - a) - 1) / h
            eps = 1e-7
            Bv, Fv, Dv, Rv, Tv = amps(ts + eps, n, Om)
            resB, resF, resD = Bv * eps, Fv * eps, Dv * eps
            # s-channel bound state: residues i x positive in the phi phi and psi psi channels, and a
            # residue matrix [[B, -D], [D, F]] of rank one (one state coupling to both channels)
            ok &= (resB / 1j).real > 1e-3 and abs((resB / 1j).imag) < 1e-6
            ok &= (resF / 1j).real > 1e-3 and abs((resF / 1j).imag) < 1e-6
            ok &= abs(resB * resF + resD ** 2) < 1e-6 and abs(Rv) < 1e3 and abs(Tv) > 1e-3
            s = 2 * 2 * (1 + np.cos((2 * (n - a) - 1) * np.pi / h))      # m_n^2 = 2
            ok &= abs(s - 8 * np.sin(np.pi * a / h) ** 2) < 1e-12
            tu = 2j * np.pi * a / h
            Bv, Fv, Dv, Rv, Tv = amps(tu + eps, n, Om)
            ok &= abs(Rv * eps) > 1e-3 and abs(Bv * eps) > 1e-3 and abs(Dv) < 1e3
    report(f'n = {n}: simple poles of B, F, D at theta = i pi(2(n-a)-1)/h, i.e. s = m_a^2 = 8 sin^2(pi a/h), '
           f'with residues i x positive in B and F and a residue matrix of rank one (Omega = 0.1, 0.31, 0.6); '
           f'poles of R, B, F at theta = 2 pi i a/h, a = 1..n-1', ok)

# ============================================================ 8. tree level
for n in (2, 3):
    h = 2 * n - 1
    th = 0.71
    s = 4 * (1 + np.cosh(th))
    uu = 4 * (1 - np.cosh(th))
    ma2 = [8 * np.sin(np.pi * a / h) ** 2 for a in range(1, n)]

    def first_order(fun, d=1e-5):
        # central difference about Omega = 0 (all amplitudes are smooth in Omega there)
        return (fun(d) - fun(-d)) / (2 * d)
    dT = first_order(lambda O: amps(th, n, O)[4])
    dR = first_order(lambda O: amps(th, n, O)[3])
    dD = first_order(lambda O: amps(th, n, O)[2])
    dB = first_order(lambda O: amps(th, n, O)[0])
    b2 = 4 * np.pi        # beta^2 = 4 pi Omega
    tT = -1j * b2 / (2 * np.sinh(th) * h)
    tR = -(1j * b2 / np.sinh(th)) * (np.cosh(th / 2) / h) * sum(
        (-1) ** a * ma2[a - 1] * np.cos(np.pi * a / h) / (uu - ma2[a - 1]) for a in range(1, n))
    tD = -(1j * b2 / np.sinh(th)) * (np.sinh(th / 2) / h) * sum(
        (-1) ** a * ma2[a - 1] * np.cos(np.pi * a / h) / (s - ma2[a - 1]) for a in range(1, n))
    tB = (1j * b2 / np.sinh(th)) * (-1 / h + 2 * np.sinh(th) ** 2 / h * sum(
        1 / (s - m2) + 1 / (uu - m2) for m2 in ma2))
    report(f'n = {n}: order-Omega terms of T, R, D, B equal the tree amplitudes (2.5) of Delius et al. '
           f'with beta^2 = 4 pi Omega', max(abs(dT - tT), abs(dR - tR), abs(dD - tD), abs(dB - tB)) < 1e-6)

# independent tree amplitude for phi_n phi_n -> phi_n phi_n from the bosonic potential
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def boson_tree_B(n, th):
    """Tree amplitude S = 1 + i M/(4 m^2 sinh theta) for the boson degenerate with the
    fermion, from V = sum_j q_j e^{alpha_j.phi} + e^{2 alpha.phi} (Delius et al. (1.5),
    beta = 1 after phi -> beta phi; the coupling enters as beta^2 per loop order)."""
    E = sp.eye(n)
    e = [E[:, i] for i in range(n)]
    roots = [-e[0] - e[1]] + [e[i] - e[i + 1] for i in range(n - 1)] + [2 * e[n - 1]]
    q = [1, 1] + [2] * (n - 2) + [1]
    x = sp.symbols(f'x0:{n}')
    X = sp.Matrix(x)
    V = sum(qq * sp.exp((r.T * X)[0, 0]) for qq, r in zip(q, roots))
    z = {xx: 0 for xx in x}
    H = np.array(sp.Matrix(n, n, lambda i, j: sp.diff(V, x[i], x[j])).subs(z), dtype=float)
    w, U = np.linalg.eigh(H)
    # the boson phi_n of mass^2 2, degenerate with the fermion (the eigenvalue 2 is simple)
    k_n = [k for k in range(n) if abs(w[k] - 2) < 1e-9]
    assert len(k_n) == 1
    k_n = k_n[0]
    d3 = {}
    d4 = 0.0
    vn = U[:, k_n]
    for k in range(n):
        vk = U[:, k]
        c3 = 0.0
        for qq, r in zip(q, roots):
            rv = np.array(r, dtype=float).ravel()
            c3 += qq * (rv @ vn) ** 2 * (rv @ vk)
        d3[k] = c3
    for qq, r in zip(q, roots):
        rv = np.array(r, dtype=float).ravel()
        d4 += qq * (rv @ vn) ** 4
    m2 = 2.0
    s_ = 2 * m2 * (1 + np.cosh(th))
    u_ = 2 * m2 * (1 - np.cosh(th))
    M = -d4 - sum(d3[k] ** 2 * (1 / (s_ - w[k]) + 1 / (0 - w[k]) + 1 / (u_ - w[k])) for k in range(n))
    return 1j * M / (4 * m2 * np.sinh(th))      # coefficient of beta^2


for n in (2, 3):
    th = 0.71
    dB = (amps(th, n, 1e-5)[0] - amps(th, n, -1e-5)[0]) / 2e-5
    lag = boson_tree_B(n, th) * 4 * np.pi       # d/dOmega with beta^2 = 4 pi Omega
    report(f'n = {n}: the order-Omega term of B equals the tree amplitude computed here from the bosonic '
           f'potential sum q_j e^(alpha_j.phi) + e^(2 alpha.phi) (contact, s-, t- and u-channel exchange)',
           abs(dB - lag) < 1e-6)

# ============================================================== 9. bootstrap
for n in (2, 3):
    h = 2 * n - 1
    Bp = 2 * Om0
    g = -np.sin(np.pi * Om0)
    ok = True
    for a in range(1, n):
        tb = np.pi * (2 * (n - a) - 1) / h / 2              # bar theta^n_an = theta^a_nn / 2
        for th in (0.4, 0.2 + 0.3j):
            lhs = (2 * g ** 2 / (1 + 1j * (-1) ** (n - a) * np.sinh(h * th)) + 1
                   - 4 * g ** 2 / np.cosh(h * th) ** 2) * Tfun(th + 1j * tb, h, Om0) * Tfun(th - 1j * tb, h, Om0)
            rhs = np.prod([brace(h / 2 - a + 2 * p - 1, th, h, Bp) for p in range(1, a + 1)])
            ok &= abs(lhs - rhs) < 1e-8
    report(f'n = {n}: the bootstrap amplitude (2.45) equals the block product prod_p {{h/2 - a + 2p - 1}} '
           f'(2.47) with B = 2 Omega', ok)

# ============================================== 10. C^(2)(n+1) fermion amplitudes
for n in (2, 3, 4, 5):
    h = 2 * n
    Bp = 2 * Om0
    ok = True
    xs1 = range(1, 2 * n - 2, 4) if n % 2 == 0 else range(3, 2 * n - 2, 4)
    xs2 = range(3, 2 * n, 4) if n % 2 == 0 else range(1, 2 * n, 4)
    for th in (0.37, -0.4 + 0.9j):
        A1 = -np.prod([brace(x, th, h, Bp) for x in xs1])          # S_xx^xx (n even), S_++^++ (n odd)
        A2 = -np.prod([brace(x, th, h, Bp) for x in xs2])          # S_xxbar (n even), S_+-^-+ (n odd)
        A1c = -np.prod([brace(x, 1j * np.pi - th, h, Bp) for x in xs1])
        A2c = -np.prod([brace(x, 1j * np.pi - th, h, Bp) for x in xs2])
        A1m = A1 * (-np.prod([brace(x, -th, h, Bp) for x in xs1]))
        A2m = A2 * (-np.prod([brace(x, -th, h, Bp) for x in xs2]))
        if n % 2 == 0:        # complex fermions chi, chibar: crossing exchanges the two amplitudes
            cross = abs(A1c - A2) < 1e-12
        else:                 # real chi_+-: each amplitude is crossing symmetric
            cross = abs(A1c - A1) < 1e-12 and abs(A2c - A2) < 1e-12
        ok &= cross and abs(A1m - 1) < 1e-12 and abs(A2m - 1) < 1e-12
    report(f'C^(2)({n}+1): the fermion amplitudes ({"3.10" if n % 2 == 0 else "3.13"}) are unitary; crossing '
           f'{"exchanges S_chichi and S_chichibar" if n % 2 == 0 else "maps S_++ and S_+- to themselves"} (h = {h})', ok)

sys.exit(1 if FAILS else 0)
