#!/usr/bin/env python3
"""Integration agent B: independent checks for Lecture 3 (half-line lectures).

Written from the lecture's and the book's definitions (no referee scripts imported).
Conventions: m = beta = 1; a_n^(1) with h = n+1, roots alpha_j = e_j - e_{j+1 mod h}
in R^h (so sum_j alpha_j = 0, alpha_j . alpha_{j+-1} = -1 for h >= 3).
Real coupling: phi = -sum_j alpha_j ln tau_j; boundary condition
  phi'(0) = -sum_j C_j alpha_j exp(alpha_j . phi(0)/2).
Imaginary coupling: phi = i sum_j alpha_j ln tau_j; boundary condition (relative to v_a)
  phi'(0) = i sum_j C_j alpha_j exp(i alpha_j . (phi(0) - v_a)/2).
"""
import itertools
import numpy as np
from scipy.integrate import quad

OK = True


def rep(name, ok, detail=""):
    global OK
    OK = OK and ok
    print(("PASS" if ok else "FAIL"), name, detail)


def roots(h):
    R = []
    for j in range(h):
        v = np.zeros(h)
        v[j] += 1.0
        v[(j + 1) % h] -= 1.0
        R.append(v)
    return R


# ---------------------------------------------------------------------------
# 1. Static stationary pair: which sign pattern does the real field satisfy?
# ---------------------------------------------------------------------------
def vac_tau(x, chi, h, a):
    s = np.sin(np.pi * a / h)
    ma = 2.0 * s
    u = np.exp(ma * x)
    T, dT = [], []
    for j in range(h):
        c = 2.0 * np.cos(chi + 2.0 * np.pi * j * a / h) / (1.0 + s)
        d = (1.0 - s) / (1.0 + s)
        T.append(1.0 + c * u + d * u * u)
        dT.append(ma * (c * u + 2.0 * d * u * u))
    return np.array(T), np.array(dT)


def bc_residuals(chi, h, a):
    R = roots(h)
    T, dT = vac_tau(0.0, chi, h, a)
    phi0 = -sum(R[j] * np.log(T[j]) for j in range(h))
    dphi0 = -sum(R[j] * dT[j] / T[j] for j in range(h))
    out = {}
    for C in itertools.product((1, -1), repeat=h):
        r = dphi0 + sum(C[j] * R[j] * np.exp(0.5 * R[j] @ phi0) for j in range(h))
        out[C] = np.linalg.norm(r)
    return out


def pattern_formula(chi, h, a):
    # C_j = sgn(c_{j+1} c_{j-1}), c_k = cos(chi/2 + pi k a/h), integer (unreduced) k
    c = lambda k: np.cos(chi / 2.0 + np.pi * k * a / h)
    return tuple(int(np.sign(c(j + 1) * c(j - 1))) for j in range(h))


def pattern_formula_reduced(chi, h, a):
    # the lecture's (wrong) reading: reduce indices mod h
    c = lambda k: np.cos(chi / 2.0 + np.pi * (k % h) * a / h)
    return tuple(int(np.sign(c(j + 1) * c(j - 1))) for j in range(h))


def topological_charge(chi, h, a):
    """Topological charge (units 2 pi / beta) of the imaginary-coupling soliton
    tau_j = 1 + zeta^{ja} e^{r + i chi}, r from -inf to +inf, continuous logs.
    phi(+inf) - phi(-inf) = i sum_j alpha_j * i*arg_j  (the real parts cancel)."""
    R = roots(h)
    lam = np.zeros(h)
    for j in range(h):
        ang = chi + 2.0 * np.pi * j * a / h
        ang_red = (ang + np.pi) % (2.0 * np.pi) - np.pi  # continuous arg reached
        lam += R[j] * (-ang_red / (2.0 * np.pi))
    return lam


def book_pattern(chi, h, a):
    R = roots(h)
    lam = topological_charge(chi, h, a)
    out = []
    for j in range(h):
        p = R[j] @ lam
        assert abs(p - round(p)) < 1e-9, p
        out.append(1 if round(p) % 2 == 0 else -1)
    return tuple(out)


def check_patterns():
    rng = np.random.default_rng(7)
    max_res_formula = 0.0
    min_res_other = 1e9
    agree_book = True
    never_uniform = True
    prod_one = True
    reduced_differs = 0
    n_samples = 0
    for h in range(3, 9):
        for a in range(1, h):
            if 2 * a == h:
                continue
            # endpoints: chi = pi - 2 pi j a/h mod 2 pi
            ends = sorted(((np.pi - 2 * np.pi * j * a / h + np.pi) % (2 * np.pi)) - np.pi
                          for j in range(h))
            ends = sorted(set(np.round(ends, 12)))
            ends = list(ends) + [ends[0] + 2 * np.pi]
            for lo, hi in zip(ends[:-1], ends[1:]):
                for frac in (0.2, 0.5, 0.8):
                    chi = lo + frac * (hi - lo)
                    chi = (chi + np.pi) % (2 * np.pi) - np.pi
                    res = bc_residuals(chi, h, a)
                    Cf = pattern_formula(chi, h, a)
                    max_res_formula = max(max_res_formula, res[Cf])
                    min_res_other = min(min_res_other, min(v for k, v in res.items() if k != Cf))
                    agree_book &= (book_pattern(chi, h, a) == Cf)
                    never_uniform &= (Cf != tuple([1] * h))
                    prod_one &= (np.prod(Cf) == 1)
                    reduced_differs += (pattern_formula_reduced(chi, h, a) != Cf)
                    n_samples += 1
    rep("static pair satisfies BC with C_j = sgn(c_{j+1}c_{j-1}) (integer indices), h=3..8",
        max_res_formula < 1e-12, f"max res {max_res_formula:.2e}, {n_samples} samples")
    rep("every other sign pattern fails", min_res_other > 1e-3, f"min res {min_res_other:.3f}")
    rep("pattern equals book prp-solitonic-vacua C_j=(-1)^{alpha_j.lambda}", agree_book)
    rep("pattern never uniform, product +1", never_uniform and prod_one)
    rep("reducing indices mod h gives a different (wrong) pattern in some samples",
        reduced_differs > 0, f"{reduced_differs} of {n_samples}")


def check_a2_table():
    table = {(-np.pi, -np.pi / 3): (-1, 1, -1),
             (-np.pi / 3, np.pi / 3): (1, -1, -1),
             (np.pi / 3, np.pi): (-1, -1, 1)}
    ok = True
    worst = 0.0
    for (lo, hi), C in table.items():
        for frac in np.linspace(0.05, 0.95, 7):
            chi = lo + frac * (hi - lo)
            res = bc_residuals(chi, 3, 1)
            best = min(res, key=res.get)
            ok &= (best == C) and res[C] < 1e-12
            worst = max(worst, res[C])
    rep("a_2 interval table (-pi,-pi/3):(-,+,-); (-pi/3,pi/3):(+,-,-); (pi/3,pi):(-,-,+)",
        ok, f"max res {worst:.2e}")
    res = bc_residuals(0.25, 3, 1)
    rep("chi=0.25 gives (+1,-1,-1); lecture's (-1,+1,-1) fails",
        res[(1, -1, -1)] < 1e-12 and res[(-1, 1, -1)] > 1.0,
        f"{res[(1, -1, -1)]:.1e} vs {res[(-1, 1, -1)]:.3f}")
    # epsilon lists at chi = 0 (k = -1..3)
    eps = [int(np.sign(np.cos(0.0 / 2 + np.pi * k / 3))) for k in range(-1, 4)]
    rep("chi=0: (eps_-1..eps_3) = (+,+,+,-,-)", eps == [1, 1, 1, -1, -1], str(eps))


def check_energy():
    R = roots(3)
    worst = 0.0
    for chi, C in ((0.25, (1, -1, -1)), (2.0, (-1, -1, 1)), (-1.3, (-1, 1, -1)), (0.9, (1, -1, -1))):
        def dens(x):
            T, dT = vac_tau(x, chi, 3, 1)
            ph = -sum(R[j] * np.log(T[j]) for j in range(3))
            dph = -sum(R[j] * dT[j] / T[j] for j in range(3))
            return 0.5 * dph @ dph + sum(np.exp(R[j] @ ph) - 1.0 for j in range(3))
        bulk = quad(dens, -np.inf, 0.0, epsabs=1e-13, epsrel=1e-13, limit=200)[0]
        T, _ = vac_tau(0.0, chi, 3, 1)
        ph0 = -sum(R[j] * np.log(T[j]) for j in range(3))
        bdy = 2.0 * sum(C[j] * (np.exp(0.5 * R[j] @ ph0) - 1.0) for j in range(3))
        worst = max(worst, abs(bulk + bdy - 2.0 * (4.0 - 3.0 * np.sqrt(3.0))))
    rep("a_2 E_vac = 2(4-3 sqrt3) with the table patterns", worst < 1e-9, f"err {worst:.1e}")
    # a_3, a=1: general formula E = 2(2k - h m_a)
    R4 = roots(4)
    worst = 0.0
    for chi in (0.3, 1.9, -2.2, -0.8):
        C = pattern_formula(chi, 4, 1)
        k = sum(1 for c in C if c == -1)
        def dens(x):
            T, dT = vac_tau(x, chi, 4, 1)
            ph = -sum(R4[j] * np.log(T[j]) for j in range(4))
            dph = -sum(R4[j] * dT[j] / T[j] for j in range(4))
            return 0.5 * dph @ dph + sum(np.exp(R4[j] @ ph) - 1.0 for j in range(4))
        bulk = quad(dens, -np.inf, 0.0, epsabs=1e-13, epsrel=1e-13, limit=200)[0]
        T, _ = vac_tau(0.0, chi, 4, 1)
        ph0 = -sum(R4[j] * np.log(T[j]) for j in range(4))
        bdy = 2.0 * sum(C[j] * (np.exp(0.5 * R4[j] @ ph0) - 1.0) for j in range(4))
        ma = 2.0 * np.sin(np.pi / 4)
        worst = max(worst, abs(bulk + bdy - 2.0 * (2 * k - 4 * ma)))
    rep("a_3 (a=1) E_vac = (2m/b^2)(2k - h m_a/m)", worst < 1e-9, f"err {worst:.1e}")


def check_inertia():
    # I = int_{-inf}^0 |d phi/d chi|^2 dx for a_2, a=1, m=beta=1; claim beta^2 m_1 I = 6
    R = roots(3)
    m1 = np.sqrt(3.0)
    vals = []
    for chi in (0.0, 0.4, 0.9):
        def integrand(x):
            s = np.sqrt(3.0) / 2
            u = np.exp(m1 * x)
            dphi = np.zeros(3)
            T, _ = vac_tau(x, chi, 3, 1)
            for j in range(3):
                dT = -2.0 * np.sin(chi + 2 * np.pi * j / 3) / (1 + s) * u
                dphi += -R[j] * dT / T[j]
            return dphi @ dphi
        I = quad(integrand, -np.inf, 0.0, epsabs=1e-13, epsrel=1e-13, limit=400)[0]
        vals.append(m1 * I)
    rep("a_2 modulus inertia beta^2 m_1 I = 6", max(abs(v - 6.0) for v in vals) < 1e-8,
        str(np.round(vals, 10)))


# ---------------------------------------------------------------------------
# 2. Positivity lower bound for the vacuum pair
# ---------------------------------------------------------------------------
def check_positivity():
    ok = True
    for s in (0.3, np.sqrt(3) / 2, 0.97):
        u = np.linspace(1e-6, 1 - 1e-6, 2001)
        bound = ((1 - u) ** 2 + s * (1 - u * u)) / (1 + s)
        worst_tau = 1 - 2 * u / (1 + s) + (1 - s) / (1 + s) * u * u
        ok &= np.max(np.abs(worst_tau - bound)) < 1e-14 and np.min(bound) > 0
    rep("tau_j >= ((1-u)^2+s(1-u^2))/(1+s) > 0 on 0<u<1 (equality at cos=-1)", ok)
    # refute the old claim that both roots lie at u >= 1: chi + 2 pi j a/h = 0
    s = np.sqrt(3) / 2
    r = np.roots([(1 - s) / (1 + s), 2 / (1 + s), 1.0])
    rep("roots at cos=+1 are -1 and -(1+s)/(1-s) (negative, not u>=1)",
        np.allclose(sorted(r.real), sorted([-1.0, -(1 + s) / (1 - s)])), str(r))


# ---------------------------------------------------------------------------
# 3. Moving mirror pair: corrected derivative, k s = 2 m sin^2, constant ratio,
#    and the branch picked by phi - v_a equals the integer-index product.
# ---------------------------------------------------------------------------
def check_mirror_square():
    worst_deriv = 0.0
    worst_ratio = 0.0
    worst_branch = 0.0
    worst_old_deriv = 0.0
    for h in (3, 4, 5):
        R = roots(h)
        zeta = np.exp(2j * np.pi / h)
        for a in range(1, h):
            for C in (1, -1):
                theta = 0.6
                xi1 = 0.31 - 0.42j
                s = np.sin(np.pi * a / h) / np.cosh(theta)
                ma = 2 * np.sin(np.pi * a / h)
                k = ma * np.cosh(theta)
                w = ma * np.sinh(theta)
                A = 1 - s * s
                S = -2 * np.log(1 + C * s)
                xi2 = S - xi1
                D = 1 + C * s
                assert abs(k * s - 2 * np.sin(np.pi * a / h) ** 2) < 1e-14

                def taus(x, t):
                    E1 = np.exp(k * x - w * t + xi1)
                    E2 = np.exp(k * x + w * t + xi2)
                    T = np.array([1 + zeta ** (j * a) * E1 + zeta ** (-j * a) * E2 + A * E1 * E2
                                  for j in range(h)])
                    dT = np.array([k * (zeta ** (j * a) * E1 + zeta ** (-j * a) * E2 + 2 * A * E1 * E2)
                                   for j in range(h)])
                    return T, dT

                # initial boundary vacuum v_a from the unwrapped log at early time
                xs = np.linspace(-30.0, 0.0, 3001)
                for t in (-1.2, 0.1, 0.9):
                    T0, _ = taus(xs, -25.0)
                    L0 = np.log(np.abs(T0)) + 1j * np.unwrap(np.angle(T0), axis=1)
                    # imaginary coupling: alpha_j.phi = i (2L_j - L_{j+1} - L_{j-1}); v_a at x=0,t->-inf
                    T, dT = taus(xs, t)
                    L = np.log(np.abs(T)) + 1j * np.unwrap(np.angle(T), axis=1)
                    # continue in t at x = 0 as well: unwrap along a t-path
                    tt = np.linspace(-25.0, t, 4001)
                    Tt = np.array([taus(0.0, tv)[0] for tv in tt]).T
                    Lt = np.log(np.abs(Tt)) + 1j * np.unwrap(np.angle(Tt), axis=1)
                    # value of L at (x=0, t=-25) from x-unwrapping, then follow t
                    Lt = Lt - Lt[:, :1] + L0[:, -1:]
                    Lx0 = Lt[:, -1]
                    # vacuum part at t->-inf: alpha_j.v_a / (2 pi) must be integer-ish (coweight)
                    # e^{i alpha_j.(phi - v_a)/2} = exp(-(2L_j - L_{j+1} - L_{j-1})/2) * e^{-i alpha_j.v_a/2}
                    a0 = np.array([(2 * L0[j, -1] - L0[(j + 1) % h, -1] - L0[(j - 1) % h, -1]) for j in range(h)])
                    # alpha_j . v_a = i * a0_j (imag-coupling phi = i sum alpha ln tau, beta=1)
                    av = 1j * a0
                    # alpha_j . v_a lies in 2 pi Z (coweight lattice); snap to it
                    assert np.max(np.abs(av / (2 * np.pi) - np.round(av.real / (2 * np.pi)))) < 1e-6
                    av = 2 * np.pi * np.round(av.real / (2 * np.pi))
                    branch = np.array([np.exp(-(2 * Lx0[j] - Lx0[(j + 1) % h] - Lx0[(j - 1) % h]) / 2)
                                       * np.exp(-0.5j * av[j]) for j in range(h)])
                    Tb, dTb = taus(0.0, t)
                    U = np.exp(-w * t + (xi1 - xi2) / 2)
                    psi = np.log(U) / 2j
                    c = lambda kk: np.cos(psi + np.pi * kk * a / h)
                    alg = np.array([4 * c(j + 1) * c(j - 1) / D for j in range(h)])
                    worst_branch = max(worst_branch, np.max(np.abs(branch * Tb - alg) / np.abs(alg)))
                    # corrected derivative and the printed (old) one
                    for j in range(h):
                        z = zeta ** (j * a) * U
                        new = k / D * (z + 1 / z + 2 * (1 - C * s))
                        old = k / D * (z + 1 / z + 2 * (1 - C * s) / (1 + C * s))
                        worst_deriv = max(worst_deriv, abs(dTb[j] - new))
                        worst_old_deriv = max(worst_old_deriv, abs(dTb[j] - old))
                    ratio = (dTb - C * alg) / Tb
                    worst_ratio = max(worst_ratio, np.max(np.abs(ratio - (k - C))))
    rep("d_x tau_j(0) = k/(1+Cs) (zU + 1/(zU) + 2(1-Cs))", worst_deriv < 1e-12, f"{worst_deriv:.1e}")
    rep("printed derivative with 2(1-Cs)/(1+Cs) is wrong", worst_old_deriv > 1e-2, f"{worst_old_deriv:.3f}")
    rep("(d_x tau_j - m C sqrt)/tau_j = k - m C, all j", worst_ratio < 1e-11, f"{worst_ratio:.1e}")
    rep("branch fixed by phi - v_a equals 4 c_{j+1} c_{j-1}/(1+Cs) (integer indices)",
        worst_branch < 1e-9, f"{worst_branch:.1e}")


# ---------------------------------------------------------------------------
# 4. Classical reflection factor in blocks
# ---------------------------------------------------------------------------
def block(x, th, h):
    return np.sinh(th / 2 + 1j * np.pi * x / (2 * h)) / np.sinh(th / 2 - 1j * np.pi * x / (2 * h))


def check_linear_reflection():
    worst_new = 0.0
    diff_old_C0 = 0.0
    for h in (3, 4, 5, 6):
        for a in range(1, h):
            for th in (0.37, 0.8 + 0.3j):
                for C in (-1, 0, 1):
                    sig = np.sin(np.pi * a / h)
                    rat = (np.sinh(th) - 1j * C * sig) / (np.sinh(th) + 1j * C * sig)
                    new = (-block(a, th, h) * block(h - a, th, h)) ** (-C)
                    old = -(block(a, th, h) * block(h - a, th, h)) ** (-C)
                    worst_new = max(worst_new, abs(rat - new))
                    if C == 0:
                        diff_old_C0 = max(diff_old_C0, abs(rat - old))
    rep("K_cl = [-(a)(h-a)]^{-C} for C = -1, 0, 1", worst_new < 1e-13, f"{worst_new:.1e}")
    rep("printed -[(a)(h-a)]^{-C} fails at C=0 by 2", abs(diff_old_C0 - 2) < 1e-14)
    # Robin derivation: ip(1-K) = -C m_a^2/(2m) (1+K)
    worst = 0.0
    for h, a, C, th in ((3, 1, 1, 0.52), (5, 2, -1, 0.3 + 0.2j)):
        ma = 2 * np.sin(np.pi * a / h)
        p = ma * np.sinh(th)
        K = (np.sinh(th) - 1j * C * np.sin(np.pi * a / h)) / (np.sinh(th) + 1j * C * np.sin(np.pi * a / h))
        worst = max(worst, abs(1j * p * (1 - K) + C * ma ** 2 / 2 * (1 + K)))
    rep("Robin condition ip(1-K) = -C m_a^2/(2m)(1+K)", worst < 1e-14)
    # Exercise 3.3: a = h/2, C = -1: pole at theta = i pi/2 (strip edge)
    K = lambda th: (np.sinh(th) + 1j) / (np.sinh(th) - 1j)
    eps = 1e-6
    e1, e2 = 1e-3, 1e-4
    r1, r2 = abs(K(0.5j * np.pi + e1)) * e1 ** 2, abs(K(0.5j * np.pi + e2)) * e2 ** 2
    rep("a=h/2, C=-1: pole at theta=i pi/2 (strip edge), K*eps^2 -> 4 (double in theta)",
        abs(r1 - 4) < 1e-5 and abs(r2 - 4) < 1e-7, f"{r1:.8f} {r2:.8f}")


# ---------------------------------------------------------------------------
# 5. Breathers (Exercise 3.2)
# ---------------------------------------------------------------------------
def check_breathers():
    y = np.linspace(-6, 0, 6001)
    ok = True
    # C=+1, rho=0: lower bound and wall zero at psi=pi
    for s in (0.4, 0.9, 1.3):
        tau_pi = 1 + 2 * np.exp(y) * (-1) / (1 + s) + (1 - s) / (1 + s) * np.exp(2 * y)
        bound = ((1 - np.exp(y)) ** 2 + s * (1 - np.exp(2 * y))) / (1 + s)
        ok &= np.max(np.abs(tau_pi - bound)) < 1e-14 and np.min(bound[:-1]) > 0 and abs(bound[-1]) < 1e-14
    rep("C=+1, rho=0: tau >= bound > 0 on y<0, tau=0 at y=0, psi=pi", ok)
    ok = True
    for s in (0.3, 0.7):
        for rho in (0.0, 0.5):
            tau = 1 - 2 * np.exp(y) * np.cosh(rho) / (1 - s) + (1 + s) / (1 - s) * np.exp(2 * y)
            if rho == 0.0:
                fac = (1 - np.exp(y)) * (1 - (1 + s) / (1 - s) * np.exp(y))
                ok &= np.max(np.abs(tau - fac)) < 1e-13
            ok &= np.min(tau[:-1]) < 0  # interior zero
    for s in (1.2, 2.0):
        for rho in (0.0, 0.5):
            tau = 1 + 2 * np.exp(y) * np.cosh(rho) / (1 - s) + (1 + s) / (1 - s) * np.exp(2 * y)
            ok &= tau[-1] < 0 and abs(tau[-1] - (2 + 2 * np.cosh(rho)) / (1 - s)) < 1e-13
    rep("C=-1: zero on y<0 (s<1 psi=pi factorization; s>1 psi=0 tau(0)<0)", ok)


if __name__ == "__main__":
    check_patterns()
    check_a2_table()
    check_energy()
    check_inertia()
    check_positivity()
    check_mirror_square()
    check_linear_reflection()
    check_breathers()
    print("ALL PASS" if OK else "SOME FAIL")
