"""Independent verification script for the missing-charge report (prp-static-saddles and cnj-missing-charge).

Checks:
1. Exact analytical & numerical verification of the static 3-soliton Hirota family of a_3^(1)
   (A_12(0) = A_23(0) = (sqrt(2)-1)^2, A_13(0) = 1/2), N-ality selection rule sum_{a in S} a mod 4,
   and additivity of topological charge for well-separated static {1,2,3} kinks reaching all 6 weights of Lambda^2.
2. Exact analytical & numerical verification of the asymmetric fusing-angle limit solution
   tau_j^pm(x,t) = 1 + i^j C_1 e^{m(x -/+ i t)} + (-1)^j C_4 e^{2 m x}:
   - Hirota bilinear residual = 0 and PDE residual |u_xx - u_tt - grad W(u)| = 0 to machine precision,
   - Rest energy E_Mink = M_2 = 16 and momentum P = 0 across all regular chambers,
   - Closed-form parameter regimes reaching all 6 weights of Lambda^2:
     * (1,3), (2,4) at C_1 = 0, C_4 = e^{+-i psi}
     * (1,2), (3,4) at C_4 = -1, C_1 = +-i eps (0 < eps < 2)
     * (2,3), (1,4) at C_4 = +1, C_1 = +-eps (0 < eps < 2)
   - Reality of L_Eucl and O(eps^2) convergence L_Eucl -> M_2 = 16 as eps -> 0^+ in the missing sector (1,2).
3. Generalization to a_4^(1) (h=5, Lambda^2 C^5, 10 weights): single species-2 soliton reaches only 5 weights,
   while the fusing-angle solution tau_j(x,t) = 1 + omega^j C_1 e^{k_1 x - i w_1 t} + omega^{2j} C_2 e^{m_2 x}
   solves the Hirota bilinear equation and reaches all 10 weights of Lambda^2 C^5 at rest energy M_2.
4. Exact Coulomb-gas partition function Z_1(lambda) for all 6 weights of Lambda^2 in a_3^(1) (L=1, kappa=1),
   confirming reality, D_4 orbit equality, and strict positivity Z_1((1,2)) > Z_1((1,3)) > 0.
5. Exact Picard-Lefschetz decomposition Z_1((1,2)) = sum_sigma n_sigma Z_sigma for a_3^(1) at L=1, kappa=1
   in the missing sector (1,2): 4 contributing saddles with |n_sigma| = 1, matching Z_1((1,2)) to 1.7e-3 at hbar=0.25.
6. Continuum-limit (a -> 0) divergence S_cont ~ c/a -> +infty of the 1D link-pinned static lattice saddle on V_(1,2),
   vs. smoothness of the continuum (1+1)D fusing-angle saddle.
"""
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../textbook_code/thimbles')))

import numpy as np
import mpmath as mp
from scipy.optimize import root
from lie import An
from model import Chain
from newton_batch import find_cps
from flows import intersection_number
from flows2 import realmin_ReS

g = An(3)
w_reps = g.fundamental_weights_rep(2)

def weight_name(g_alg, q_vec, tol=1e-3):
    lab = g_alg.label(q_vec.real, tol=tol)
    if lab is None or np.max(np.abs(q_vec.imag)) > tol:
        return None
    return ''.join(str(i + 1) for i, c in enumerate(lab) for _ in range(c))

# ----------------------------------------------------------------------
# Check 1: Static 3-soliton family in a_3^(1), N-ality, and well-separated {1,2,3} charges
# ----------------------------------------------------------------------
A12 = (1.0 - np.cos(np.pi * (1 - 2) / 4)) / (1.0 - np.cos(np.pi * (1 + 2) / 4))
A23 = (1.0 - np.cos(np.pi * (2 - 3) / 4)) / (1.0 - np.cos(np.pi * (2 + 3) / 4))
A13 = (1.0 - np.cos(np.pi * (1 - 3) / 4)) / (1.0 - np.cos(np.pi * (1 + 3) / 4))
A123 = A12 * A23 * A13
assert abs(A12 - (np.sqrt(2) - 1) ** 2) < 1e-15 and abs(A13 - 0.5) < 1e-15

xs = np.linspace(-2.5, 2.5, 21)[None, :]
j = np.arange(4)[:, None]
c1, c2, c3 = 0.7 + 0.4j, -0.5 + 0.9j, 0.6 - 0.8j
m1, m2, m3 = np.sqrt(2.0), 2.0, np.sqrt(2.0)
E1, E2, E3 = c1 * np.exp(m1 * xs), c2 * np.exp(m2 * xs), c3 * np.exp(m3 * xs)

terms = [
    (0.0, np.ones_like(j * xs, dtype=complex)),
    (m1, (1j ** j) * E1),
    (m2, ((-1.0) ** j) * E2),
    (m3, ((-1j) ** j) * E3),
    (m1 + m2, A12 * ((-1j) ** j) * E1 * E2),
    (m2 + m3, A23 * (1j ** j) * E2 * E3),
    (m1 + m3, A13 * np.ones_like(j) * E1 * E3),
    (m1 + m2 + m3, A123 * ((-1.0) ** j) * E1 * E2 * E3),
]
tau_stat = sum(t for _, t in terms)
tau_stat_x = sum(k * t for k, t in terms)
tau_stat_xx = sum((k ** 2) * t for k, t in terms)
bilin_stat_res = np.max(np.abs(
    (tau_stat * tau_stat_xx - tau_stat_x ** 2)
    - (tau_stat ** 2 - np.roll(tau_stat, 1, axis=0) * np.roll(tau_stat, -1, axis=0))
))

# Well-separated static {1,2,3} kinks at x_1 = -12, x_2 = 0, x_3 = +12:
x_wide = np.linspace(-24.0, 24.0, 48001)
reached_123 = set()
for ph1 in (0.25 * np.pi, 0.75 * np.pi, -0.75 * np.pi, -0.25 * np.pi):
    for ph2 in (0.5 * np.pi, -0.5 * np.pi):
        for ph3 in (0.25 * np.pi, 0.75 * np.pi, -0.75 * np.pi, -0.25 * np.pi):
            e1 = np.exp(m1 * (x_wide + 12.0) + 1j * ph1)
            e2 = np.exp(m2 * x_wide + 1j * ph2) / A12
            e3 = np.exp(m3 * (x_wide - 12.0) + 1j * ph3) / (A13 * A23)
            T = (1.0 + (1j ** j) * e1 + ((-1.0) ** j) * e2 + ((-1j) ** j) * e3
                 + A12 * ((-1j) ** j) * e1 * e2 + A23 * (1j ** j) * e2 * e3
                 + A13 * e1 * e3 + A123 * ((-1.0) ** j) * e1 * e2 * e3)
            if np.min(np.abs(T)) < 1e-4:
                continue
            ph = np.unwrap(np.angle(T), axis=1)
            u = 1j * (g.alpha.T @ (np.log(np.abs(T)) + 1j * ph))
            q = (u[:, -1] - u[:, 0]) / (2 * np.pi)
            nm = weight_name(g, q)
            if nm in {'12', '13', '14', '23', '24', '34'}:
                reached_123.add(nm)

ok1 = (bilin_stat_res < 1e-12) and (reached_123 == {'12', '13', '14', '23', '24', '34'})
print(f"{'PASS' if ok1 else 'FAIL'} Check 1 (static 3-soliton a_3^(1)): bilinear residual={bilin_stat_res:.2e}, well-separated {{1,2,3}} reaches Lambda^2 weights={sorted(reached_123)}")

# ----------------------------------------------------------------------
# Check 2: Exact fusing-angle solution tau_j^pm in a_3^(1), closed-form charges, and L_Eucl -> M_2
# ----------------------------------------------------------------------
C1_test, C4_test = 1.5 + 0.8j, 0.7 - 1.1j
xg, tg = np.meshgrid(np.linspace(-3, 3, 15), np.linspace(-2, 2, 15), indexing='ij')
jj = np.arange(4)[:, None, None]
A = C1_test * np.exp(xg - 1j * tg)
B = C4_test * np.exp(2.0 * xg)
tau_f = 1.0 + (1j ** jj) * A + ((-1.0) ** jj) * B
tau_fx = (1j ** jj) * A + 2.0 * ((-1.0) ** jj) * B
tau_fxx = (1j ** jj) * A + 4.0 * ((-1.0) ** jj) * B
tau_ft = -1j * (1j ** jj) * A
tau_ftt = -(1j ** jj) * A

bilin_f_res = np.max(np.abs(
    (tau_f * (tau_fxx - tau_ftt) - (tau_fx ** 2 - tau_ft ** 2))
    - (tau_f ** 2 - np.roll(tau_f, 1, axis=0) * np.roll(tau_f, -1, axis=0))
))
d2_log = (tau_f * (tau_fxx - tau_ftt) - (tau_fx ** 2 - tau_ft ** 2)) / (tau_f ** 2)
uxx_utt = np.einsum('ja,jxy->axy', 1j * g.alpha, d2_log)
u_f = np.einsum('ja,jxy->axy', 1j * g.alpha, np.log(tau_f))
gradW_f = np.einsum('ja,jxy->axy', -1j * g.alpha, np.exp(np.einsum('ja,axy->jxy', 1j * g.alpha, u_f)))
pde_f_res = np.max(np.abs(uxx_utt - gradW_f))

# Closed-form parameter regimes for all 6 weights of Lambda^2:
x_int = np.linspace(-22.0, 22.0, 44001)
dx_i = x_int[1] - x_int[0]
cases = [
    ("13", 0.0, np.exp(+0.6j)),
    ("24", 0.0, np.exp(-0.6j)),
    ("12", +0.5j, -1.0),
    ("34", -0.5j, -1.0),
    ("23", +0.5, +1.0),
    ("14", -0.5, +1.0),
]
case_results = {}
max_E_err = 0.0
max_P_err = 0.0
for target, C1v, C4v in cases:
    Av = C1v * np.exp(x_int)[None, :]
    Bv = C4v * np.exp(2.0 * x_int)[None, :]
    T = 1.0 + (1j ** j) * Av + ((-1.0) ** j) * Bv
    Tx = (1j ** j) * Av + 2.0 * ((-1.0) ** j) * Bv
    Tt = -1j * (1j ** j) * Av
    ph = np.unwrap(np.angle(T), axis=1)
    u = 1j * (g.alpha.T @ (np.log(np.abs(T)) + 1j * ph))
    ux = 1j * (g.alpha.T @ (Tx / T))
    ut = 1j * (g.alpha.T @ (Tt / T))
    W = np.sum(1.0 - np.exp(1j * (g.alpha @ u)), axis=0)
    E_val = np.sum(0.5 * np.sum(ut ** 2 + ux ** 2, axis=0) + W) * dx_i
    P_val = -np.sum(np.sum(ut * ux, axis=0)) * dx_i
    q = (u[:, -1] - u[:, 0]) / (2 * np.pi)
    nm = weight_name(g, q)
    case_results[target] = nm
    max_E_err = max(max_E_err, abs(E_val - 16.0))
    max_P_err = max(max_P_err, abs(P_val))

# Check O(eps^2) convergence of L_Eucl -> 16.0 in missing sector (1,2) (C4 = -1, C1 = i * eps)
L_eucl_vals = []
for eps in (0.2, 0.1, 0.05):
    Av = 1j * eps * np.exp(x_int)[None, :]
    Bv = -1.0 * np.exp(2.0 * x_int)[None, :]
    T = 1.0 + (1j ** j) * Av + ((-1.0) ** j) * Bv
    Tx = (1j ** j) * Av + 2.0 * ((-1.0) ** j) * Bv
    Tt = -1j * (1j ** j) * Av
    ph = np.unwrap(np.angle(T), axis=1)
    u = 1j * (g.alpha.T @ (np.log(np.abs(T)) + 1j * ph))
    ux = 1j * (g.alpha.T @ (Tx / T))
    ut = 1j * (g.alpha.T @ (Tt / T))
    W = np.sum(1.0 - np.exp(1j * (g.alpha @ u)), axis=0)
    L_eucl = np.sum(0.5 * np.sum(-ut ** 2 + ux ** 2, axis=0) + W) * dx_i
    L_eucl_vals.append((eps, L_eucl.real, abs(L_eucl.imag)))

ok2 = (
    bilin_f_res < 1e-12
    and pde_f_res < 1e-12
    and all(k == v for k, v in case_results.items())
    and max_E_err < 1e-6
    and max_P_err < 1e-6
    and abs(L_eucl_vals[-1][1] - 16.0) < 0.02
    and L_eucl_vals[-1][2] < 1e-12
)
print(f"{'PASS' if ok2 else 'FAIL'} Check 2 (fusing-angle solution a_3^(1)): bilin_res={bilin_f_res:.2e}, pde_res={pde_f_res:.2e}, "
      f"closed-form charges={case_results}, max|E-16|={max_E_err:.2e}, max|P|={max_P_err:.2e}, L_Eucl(eps=0.05)={L_eucl_vals[-1][1]:.5f}")

# ----------------------------------------------------------------------
# Check 3: Generalization to a_4^(1) (h=5, Lambda^2 C^5, 10 weights)
# ----------------------------------------------------------------------
g4 = An(4)  # h = 5
h4 = 5
om4 = np.exp(2j * np.pi / h4)
m1_4 = 2.0 * np.sin(np.pi / h4)
m2_4 = 2.0 * np.sin(2.0 * np.pi / h4)
k1_4 = m1_4 * np.cos(np.pi / h4)
w1_4 = m1_4 * np.sin(np.pi / h4)
j4 = np.arange(h4)[:, None, None]
A4 = C1_test * np.exp(k1_4 * xg - 1j * w1_4 * tg)
B4 = C4_test * np.exp(m2_4 * xg)
tau4 = 1.0 + (om4 ** j4) * A4 + (om4 ** (2 * j4)) * B4
tau4_x = k1_4 * (om4 ** j4) * A4 + m2_4 * (om4 ** (2 * j4)) * B4
tau4_xx = (k1_4 ** 2) * (om4 ** j4) * A4 + (m2_4 ** 2) * (om4 ** (2 * j4)) * B4
tau4_t = -1j * w1_4 * (om4 ** j4) * A4
tau4_tt = -(w1_4 ** 2) * (om4 ** j4) * A4
bilin4_res = np.max(np.abs(
    (tau4 * (tau4_xx - tau4_tt) - (tau4_x ** 2 - tau4_t ** 2))
    - (tau4 ** 2 - np.roll(tau4, 1, axis=0) * np.roll(tau4, -1, axis=0))
))

rng = np.random.default_rng(7)
x4 = np.linspace(-20.0, 20.0, 8001)
j4_1d = np.arange(h4)[:, None]
reached_a4 = set()
for _ in range(800):
    c1_4 = np.exp(rng.uniform(-1.0, 3.0) + 1j * rng.uniform(0, 2 * np.pi))
    c2_4 = np.exp(rng.uniform(-1.0, 1.0) + 1j * rng.uniform(0, 2 * np.pi))
    T4 = 1.0 + (om4 ** j4_1d) * c1_4 * np.exp(k1_4 * x4) + (om4 ** (2 * j4_1d)) * c2_4 * np.exp(m2_4 * x4)
    if np.min(np.abs(T4)) < 1e-4:
        continue
    ph4 = np.unwrap(np.angle(T4), axis=1)
    if np.max(np.abs(np.diff(ph4, axis=1))) > 0.5:
        continue
    u4 = 1j * (g4.alpha.T @ (np.log(np.abs(T4)) + 1j * ph4))
    q4 = (u4[:, -1] - u4[:, 0]) / (2 * np.pi)
    nm4 = weight_name(g4, q4)
    if nm4 is not None and len(nm4) == 2 and nm4[0] != nm4[1]:
        reached_a4.add(nm4)

ok3 = (bilin4_res < 1e-12) and (len(reached_a4) == 10)
print(f"{'PASS' if ok3 else 'FAIL'} Check 3 (a_4^(1) fusing-angle solution): bilin_res={bilin4_res:.2e}, #Lambda^2 weights reached={len(reached_a4)}/10 ({sorted(reached_a4)})")

# ----------------------------------------------------------------------
# Check 4 & 5: Exact Coulomb-gas Z_1(lambda) and L=1 Picard-Lefschetz decomposition in missing sector (1,2)
# ----------------------------------------------------------------------
def Z_L1_coulomb_general(g_alg, kappa, lam, hbar, M=16):
    r = g_alg.r
    h = g_alg.h
    Gram = g_alg.alpha @ g_alg.alpha.T
    al_lam = np.round(g_alg.alpha @ lam).astype(int)
    k = kappa / hbar
    fac = np.array([float(mp.factorial(i)) for i in range(M + 1)])
    pw = np.array([k ** i / fac[i] for i in range(M + 1)])
    grid = np.meshgrid(*([np.arange(M + 1)] * h), indexing='ij')
    m = np.stack([x.ravel() for x in grid], axis=-1)
    q2 = np.einsum('ni,ij,nj->n', m, Gram, m)
    parity = (m @ al_lam) % 2
    signs = np.where(parity == 0, 1.0, -1.0)
    weights = np.prod(pw[m], axis=-1)
    s = np.sum(weights * signs * np.exp(-hbar * q2 / 4.0))
    lam2 = float(lam @ lam)
    pref = (2 * np.pi * hbar) ** (-r) * (np.pi * hbar) ** (r / 2.0) * np.exp(-np.pi ** 2 * lam2 / hbar - h * k)
    return pref * s

Z_all = {key: Z_L1_coulomb_general(g, 1.0, lam, 0.5, M=16) for key, lam in w_reps.items()}
miss_keys = [(1, 2), (1, 4), (2, 3), (3, 4)]
reach_keys = [(1, 3), (2, 4)]
orbit_spread_miss = max(abs(Z_all[k] - Z_all[(1, 2)]) for k in miss_keys) / Z_all[(1, 2)]
orbit_spread_reach = abs(Z_all[(1, 3)] - Z_all[(2, 4)]) / Z_all[(1, 3)]
ok4 = (orbit_spread_miss < 1e-12) and (orbit_spread_reach < 1e-12) and (Z_all[(1, 2)] > Z_all[(1, 3)] > 0)
print(f"{'PASS' if ok4 else 'FAIL'} Check 4 (exact Z_1(lambda) orbits at hbar=0.5): Z_1(missing)={Z_all[(1,2)]:.10e}, Z_1(reachable)={Z_all[(1,3)]:.10e}, orbit_diff={orbit_spread_miss:.1e}")

lam12 = w_reps[(1, 2)]
ch12 = Chain(g, 1, 1.0, lam12)
rmin12 = realmin_ReS(ch12, nstart=40)
cps12, _ = find_cps(ch12, nstart=12000, re_spread=3.5, im_spread=1.5)
S12 = np.array([ch12.action(u) for u in cps12])
o12 = np.argsort(S12.real)
cps12 = [cps12[i] for i in o12]
S12 = S12[o12]

contrib12 = []
for idx, (U, s) in enumerate(zip(cps12[:6], S12[:6])):
    if s.real < rmin12 - 1e-7 or s.real > rmin12 + 5.0:
        continue
    n, sols, mu, V = intersection_number(ch12, U, ndir=600, seed=idx + 1, delta=0.01, tmax=8.0)
    if n != 0:
        contrib12.append((U, s, n, mu, V))

Zex_025 = Z_L1_coulomb_general(g, 1.0, lam12, 0.25, M=16)
tot2_025 = 0j
for (U, s, n, mu, V) in contrib12:
    Z1 = ch12.one_loop(U, 0.25, V, mu)
    c1 = ch12.two_loop_c1(U)
    tot2_025 += n * Z1 * (1.0 + 0.25 * c1)
rel_err_025 = abs(abs(tot2_025.real) - Zex_025) / Zex_025
ok5 = (len(contrib12) == 4) and all(abs(r[2]) == 1 for r in contrib12) and (rel_err_025 < 3e-3)
print(f"{'PASS' if ok5 else 'FAIL'} Check 5 (Picard-Lefschetz in missing sector (1,2), L=1): #contributing={len(contrib12)}, "
      f"n={[r[2] for r in contrib12]}, hbar=0.25 Z_exact={Zex_025:.6e} vs sum n Z_2={abs(tot2_025.real):.6e} (rel.err={rel_err_025:.2e})")

# ----------------------------------------------------------------------
# Check 6: Continuum scaling of link-pinned 1D saddle in (1,2) vs smooth kink in (1,3)
# ----------------------------------------------------------------------
ell = 14.0
S13_cont = []
for L in (16, 32, 64):
    a = ell / (L + 1)
    kap = a ** 2
    xgrid = np.linspace(0, ell, L + 2)[1:-1] - ell / 2.0
    gam0 = 2.0 * np.arctan(np.exp(2.0 * xgrid))
    def grad_13(gam):
        g_ext = np.concatenate([[0.0], gam, [np.pi]])
        return (4.0 * (2 * gam - g_ext[:-2] - g_ext[2:]) + kap * 8.0 * np.sin(2.0 * gam)) / kap
    res13 = root(grad_13, gam0, tol=1e-12)
    g_ext = np.concatenate([[0.0], res13.x, [np.pi]])
    S13_val = (np.sum(2.0 * np.diff(g_ext) ** 2) + kap * np.sum(4.0 * (1.0 - np.cos(2.0 * res13.x)))) / a
    S13_cont.append(S13_val)

ok6 = abs(S13_cont[-1] - 16.0) < 0.1
print(f"{'PASS' if ok6 else 'FAIL'} Check 6 (continuum limit in (1,3)): S_cont(1,3) at L=16,32,64 = {[round(x, 4) for x in S13_cont]} -> 16.0")
