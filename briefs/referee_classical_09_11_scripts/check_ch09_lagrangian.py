#!/usr/bin/env python3
"""
Independent referee verification script for Chapter 9:
Ranks 8, 10, 11, 12, 14, 16.
Does not import any code from textbook_code/.
"""

import numpy as np
import sympy as sp

print("=" * 72)
print("CHECK 1 (Rank 8): Spin-3 current of sinh-Gordon and sine-Gordon, and kink evaluation")
print("=" * 72)

phi, p_p, p_pp, p_ppp, m, beta = sp.symbols("phi p_p p_pp p_ppp m beta", real=True)

# 1. Sinh-Gordon: d_+ d_- phi = -(m^2/beta) sinh(beta*phi)
F_shg = -(m**2 / beta) * sp.sinh(beta * phi)
F_shg_p = sp.diff(F_shg, phi) * p_p

T4_shg = p_pp**2 + (beta**2 / 4) * p_p**4
X2_shg = m**2 * sp.cosh(beta * phi) * p_p**2

d_minus_T4_shg = sp.diff(T4_shg, p_p) * F_shg + sp.diff(T4_shg, p_pp) * F_shg_p
d_plus_X2_shg = sp.diff(X2_shg, phi) * p_p + sp.diff(X2_shg, p_p) * p_pp
res_shg = sp.simplify(d_minus_T4_shg + d_plus_X2_shg)
print("1. Sinh-Gordon conservation residual d_- T_4 + d_+ X_2 =", res_shg)

# 2. Sine-Gordon: d_+ d_- phi = -(m^2/beta) sin(beta*phi)
F_sg = -(m**2 / beta) * sp.sin(beta * phi)
F_sg_p = sp.diff(F_sg, phi) * p_p

# Test what happens if we literally only "replace beta^2 by -beta^2" (keeping cosh(beta*phi)):
T4_sg = p_pp**2 - (beta**2 / 4) * p_p**4
X2_sg_wrong = m**2 * sp.cosh(beta * phi) * p_p**2
X2_sg_right = m**2 * sp.cos(beta * phi) * p_p**2

d_minus_T4_sg = sp.diff(T4_sg, p_p) * F_sg + sp.diff(T4_sg, p_pp) * F_sg_p
d_plus_X2_wrong = sp.diff(X2_sg_wrong, phi) * p_p + sp.diff(X2_sg_wrong, p_p) * p_pp
d_plus_X2_right = sp.diff(X2_sg_right, phi) * p_p + sp.diff(X2_sg_right, p_p) * p_pp

print(
    "2. Sine-Gordon with ONLY beta^2 -> -beta^2 (keeping cosh): d_- T_4 + d_+ X_2 =",
    sp.simplify(d_minus_T4_sg + d_plus_X2_wrong),
)
print(
    "3. Sine-Gordon with beta -> i*beta (beta^2 -> -beta^2 AND cosh -> cos): d_- T_4 + d_+ X_2 =",
    sp.simplify(d_minus_T4_sg + d_plus_X2_right),
)

# 3. Evaluate Q_3 = int_{-infty}^{+infty} (T_4 + X_2) dx on the sine-Gordon kink of @sec-sg-classical:
# phi(x,t) = (4/beta) arctan(exp(m (x cosh theta - t sinh theta)))
x_s, t_s, theta_s = sp.symbols("x t theta", real=True)
u_s = m * (x_s * sp.cosh(theta_s) - t_s * sp.sinh(theta_s))
# Notice: d_t u_s = -m sinh(theta), d_x u_s = m cosh(theta)
# So d_+ u_s = (d_t + d_x) u_s = m (cosh(theta) - sinh(theta)) = m e^{-theta}!
# And d_- u_s = (d_t - d_x) u_s = -m (cosh(theta) + sinh(theta)) = -m e^{+theta}!
# With phi = (4/beta) arctan(e^u):
# dphi/du = (2/beta) sech(u), d^2phi/du^2 = -(2/beta) sech(u) tanh(u),
# cos(beta*phi) = cos(4 arctan(e^u)) = 1 - 2 sech^2(u) (Wait: let's check cos(4 arctan(e^u))!).
u_sym = sp.Symbol("u", real=True)
phi_u = (4 / beta) * sp.atan(sp.exp(u_sym))
dphi_du = sp.simplify(sp.diff(phi_u, u_sym))
d2phi_du2 = sp.simplify(sp.diff(dphi_du, u_sym))
cos_bphi = sp.simplify(sp.cos(beta * phi_u))
print("   On SG kink: dphi/du =", dphi_du, ", d^2phi/du^2 =", d2phi_du2, ", cos(beta*phi) =", cos_bphi)

# Since d_+ = (m e^{-theta}) d/du:
p_p_kink = (m * sp.exp(-theta_s)) * dphi_du
p_pp_kink = (m * sp.exp(-theta_s)) ** 2 * d2phi_du2
T4_kink = sp.simplify(p_pp_kink**2 - (beta**2 / 4) * p_p_kink**4)
X2_kink = sp.simplify(m**2 * cos_bphi * p_p_kink**2)
print("   On SG kink: T_4 =", T4_kink)
print("   On SG kink: X_2 =", X2_kink)
# Change integration variable from x to u = m (x cosh theta - t sinh theta) at fixed t: dx = du / (m cosh theta)
int_T4 = sp.integrate(T4_kink, (u_sym, -sp.oo, sp.oo)) / (m * sp.cosh(theta_s))
int_X2 = sp.integrate(X2_kink, (u_sym, -sp.oo, sp.oo)) / (m * sp.cosh(theta_s))
Q3_kink = sp.simplify(int_T4 + int_X2)
print("   int T_4 dx =", sp.simplify(int_T4))
print("   int X_2 dx =", sp.simplify(int_X2))
print("   Q_3 = int (T_4 + X_2) dx =", Q3_kink)

print("\n" + "=" * 72)
print("CHECK 2 (Rank 16): Lax pair zero-curvature, boost law, and Exercise 9.1")
print("=" * 72)

# Check Exercise 9.1 (a_1^(1) 2x2 matrices):
# alpha_1 = sqrt(2), alpha_0 = -sqrt(2), n_0 = n_1 = 1 => c_0 = c_1 = sqrt(1 * 2 / 2) = 1.
lam = sp.Symbol("lam")
E_a1 = sp.Matrix([[0, 1], [0, 0]])
E_a0 = sp.Matrix([[0, 0], [1, 0]])
F_a1 = E_a1.T
F_a0 = E_a0.T
H_mat = (1 / sp.sqrt(2)) * sp.Matrix([[1, 0], [0, -1]])

phi_p, phi_m, phi_pm = sp.symbols("phi_p phi_m phi_pm")
A_plus = (beta / 2) * phi_p * H_mat + m * lam * (
    sp.exp(beta * sp.sqrt(2) * phi / 2) * E_a1
    + sp.exp(-beta * sp.sqrt(2) * phi / 2) * E_a0
)
A_minus = -(beta / 2) * phi_m * H_mat - (m / lam) * (
    sp.exp(beta * sp.sqrt(2) * phi / 2) * F_a1
    + sp.exp(-beta * sp.sqrt(2) * phi / 2) * F_a0
)

# [d_+ + A_+, d_- + A_-] = d_+ A_- - d_- A_+ + [A_+, A_-]
d_plus_A_minus = sp.diff(A_minus, phi) * phi_p + sp.diff(A_minus, phi_m) * phi_pm
d_minus_A_plus = sp.diff(A_plus, phi) * phi_m + sp.diff(A_plus, phi_p) * phi_pm
comm = A_plus * A_minus - A_minus * A_plus
F_curv = sp.simplify(d_plus_A_minus - d_minus_A_plus + comm)
print("Zero-curvature matrix [d_+ + A_+, d_- + A_-] for a_1^(1):")
print("  (0,0) entry =", F_curv[0, 0])
print("  (0,1) entry =", F_curv[0, 1])
print("  (1,0) entry =", F_curv[1, 0])
print("  (1,1) entry =", F_curv[1, 1])
print(
    "  Vanishes iff phi_pm =",
    sp.solve(F_curv[0, 0], phi_pm)[0],
    "= -(m^2/beta) sqrt(2) sinh(sqrt(2) beta phi) (matches eq-toda-eom!).",
)

print("\n" + "=" * 72)
print("CHECK 3 (Rank 14): Particle masses for a_n^(1), d_n^(1), e_8^(1) and sum rule")
print("=" * 72)

# Check d_n^(1) for n = 4, 5, 6:
for n_val in [4, 5, 6]:
    h_val = 2 * n_val - 2
    # Simple roots of d_n^(1) in R^n:
    # alpha_1 = e_1 - e_2, ..., alpha_{n-1} = e_{n-1} - e_n, alpha_n = e_{n-1} + e_n, alpha_0 = -(e_1 + e_2)
    e = np.eye(n_val)
    alphas = [-(e[0] + e[1])]
    for j in range(n_val - 1):
        alphas.append(e[j] - e[j + 1])
    alphas.append(e[n_val - 2] + e[n_val - 1])
    alphas = np.array(alphas)
    n_kac = np.array([1, 1] + [2] * (n_val - 3) + [1, 1], dtype=float)
    assert np.allclose(n_kac @ alphas, 0)
    M2 = sum(n_kac[j] * np.outer(alphas[j], alphas[j]) for j in range(n_val + 1))
    evals = np.sort(np.linalg.eigvalsh(M2))
    masses_num = np.sqrt(evals)
    masses_tab = np.sort(
        [2 * np.sqrt(2) * np.sin(np.pi * a / h_val) for a in range(1, n_val - 1)]
        + [np.sqrt(2), np.sqrt(2)]
    )
    assert np.allclose(masses_num, masses_tab)
    assert np.isclose(np.sum(evals), 2 * h_val)
    print(f"d_{n_val}^(1) (h={h_val}): masses and sum rule 2h={2*h_val} match to {np.max(np.abs(masses_num - masses_tab)):.2e}")

# Check e_8^(1) (h = 30):
# Affine E_8 Cartan matrix: chain 0-1-2-3-4-5-6-7 with node 8 attached to node 5 (where 0 is attached to 1)
# Wait: standard E_8 has nodes 1..8 with branch at 5 (or 3) and Kac labels 2,3,4,5,6,4,2,3 and n_0=1.
A_e8 = np.array(
    [
        [2, -1, 0, 0, 0, 0, 0, 0, 0],  # 0 (n=1) connected to 1 (n=2)
        [-1, 2, -1, 0, 0, 0, 0, 0, 0],  # 1 (n=2)
        [0, -1, 2, -1, 0, 0, 0, 0, 0],  # 2 (n=3)
        [0, 0, -1, 2, -1, 0, 0, 0, 0],  # 3 (n=4)
        [0, 0, 0, -1, 2, -1, 0, 0, 0],  # 4 (n=5)
        [0, 0, 0, 0, -1, 2, -1, 0, -1],  # 5 (n=6) connected to 4, 6, 8
        [0, 0, 0, 0, 0, -1, 2, -1, 0],  # 6 (n=4)
        [0, 0, 0, 0, 0, 0, -1, 2, 0],  # 7 (n=2)
        [0, 0, 0, 0, 0, -1, 0, 0, 2],  # 8 (n=3)
    ],
    dtype=float,
)
n_e8 = np.array([1, 2, 3, 4, 5, 6, 4, 2, 3], dtype=float)
assert np.allclose(A_e8 @ n_e8, 0) and np.sum(n_e8) == 30
K_e8 = np.outer(np.sqrt(n_e8), np.sqrt(n_e8)) * A_e8
evals_e8 = np.sort(np.linalg.eigvalsh(K_e8))[1:]
masses_e8 = np.sqrt(evals_e8)
m1_e8 = masses_e8[0]
m2_e8 = 2 * m1_e8 * np.cos(np.pi / 5)
m3_e8 = 2 * m1_e8 * np.cos(np.pi / 30)
m4_e8 = 2 * m2_e8 * np.cos(7 * np.pi / 30)
m5_e8 = 2 * m2_e8 * np.cos(2 * np.pi / 15)
m6_e8 = 2 * m2_e8 * np.cos(np.pi / 30)
m7_e8 = 4 * m2_e8 * np.cos(np.pi / 5) * np.cos(7 * np.pi / 30)
m8_e8 = 4 * m2_e8 * np.cos(np.pi / 5) * np.cos(2 * np.pi / 15)
tab_e8 = np.array([m1_e8, m2_e8, m3_e8, m4_e8, m5_e8, m6_e8, m7_e8, m8_e8])
print(
    f"e_8^(1) (h=30): max diff between eigenvalues and table formulas = {np.max(np.abs(masses_e8 - tab_e8)):.2e}, "
    f"sum m_a^2 = {np.sum(evals_e8):.6f} (2h = 60)"
)
assert np.allclose(masses_e8, tab_e8)

print("\n" + "=" * 72)
print("CHECK 4 (Rank 11 & Rank 12): Couplings C_{abc}, Area Rule, and Dorey's Rule")
print("=" * 72)


def triangle_area(a, b, c):
    s = 0.5 * (a + b + c)
    val = s * (s - a) * (s - b) * (s - c)
    if val <= 1e-12:
        return 0.0
    return np.sqrt(val)


# 1) Check a_2^(1) (Exercise 9.5 & 9.6):
# Simple roots of a_2^(1):
alphas_a2 = np.array(
    [
        [-1.0 / np.sqrt(2.0), -np.sqrt(1.5)],
        [np.sqrt(2.0), 0.0],
        [-1.0 / np.sqrt(2.0), np.sqrt(1.5)],
    ]
)
omega3 = np.exp(2j * np.pi / 3.0)
# Z_3 charge eigenvectors e_1, e_2 with e_1 . e_2 = 1, e_2 = bar(e_1):
e1_a2 = np.array([1.0, -1j]) / np.sqrt(2.0)
e2_a2 = np.array([1.0, +1j]) / np.sqrt(2.0)
es_a2 = [e1_a2, e2_a2]
print("a_2^(1) cubic couplings in Z_3 charge basis (m=1, beta=1):")
for a in range(2):
    for b in range(a, 2):
        for c in range(b, 2):
            C_val = sum(
                np.dot(alphas_a2[j], es_a2[a])
                * np.dot(alphas_a2[j], es_a2[b])
                * np.dot(alphas_a2[j], es_a2[c])
                for j in range(3)
            )
            area_rhs = (4.0 / np.sqrt(3.0)) * triangle_area(
                np.sqrt(3.0), np.sqrt(3.0), np.sqrt(3.0)
            )
            print(
                f"  |C_{{{a+1}{b+1}{c+1}}}| = {abs(C_val):.6f}  (4/sqrt(3) * Delta = {area_rhs:.6f})"
            )

# Check what happens in a REAL orthonormal basis of a_2^(1):
r1_a2 = np.array([1.0, 0.0])
r2_a2 = np.array([0.0, 1.0])
rs_a2 = [r1_a2, r2_a2]
print("a_2^(1) cubic couplings in REAL orthonormal basis (m=1, beta=1):")
for a in range(2):
    for b in range(a, 2):
        for c in range(b, 2):
            C_val = sum(
                np.dot(alphas_a2[j], rs_a2[a])
                * np.dot(alphas_a2[j], rs_a2[b])
                * np.dot(alphas_a2[j], rs_a2[c])
                for j in range(3)
            )
            print(f"  |C^real_{{{a+1}{b+1}{c+1}}}| = {abs(C_val):.6f}")

# 2) Check d_4^(1) (even n = 4, where all particles are self-conjugate and the 3 outer particles 1, 3, 4 are degenerate with m = sqrt(2)):
# Look at 09-lagrangian.qmd:130:
# "For degenerate masses (a_n, d_n, e_6) the e_a must be charge eigenvectors, the complex eigenvectors of the centre of the affine diagram (of Z_h for a_n), normalized by e_a . e_{bar a} = 1 so that phi_{bar a} = bar(phi_a). Particle a is the eigenvector with the centre character of node a. In a real orthonormal basis of a degenerate pair neither (-@eq-area-rule) nor @prp-dorey-rule holds."
# WAIT A SECOND!
# In d_{2m}^(1) (such as d_4^(1) or d_6^(1)), the centre of D_{2m} is Z_2 x Z_2, whose characters are ALL REAL (+1 and -1)!
# Because every character of Z_2 x Z_2 is REAL (+1 or -1), the charge eigenvectors e_a in d_{2m}^(1) ARE REAL ORTHONORMAL VECTORS (with bar(a) = a and e_a . e_a = 1)!
# Let's check d_4^(1) right now!
e_d4 = np.eye(4)
alphas_d4 = np.array(
    [
        -(e_d4[0] + e_d4[1]),  # alpha_0 (n=1)
        e_d4[0] - e_d4[1],     # alpha_1 (n=1, vector)
        e_d4[1] - e_d4[2],     # alpha_2 (n=2, central)
        e_d4[2] - e_d4[3],     # alpha_3 (n=1, spinor s)
        e_d4[2] + e_d4[3],     # alpha_4 (n=1, spinor s')
    ]
)
n_d4 = np.array([1.0, 1.0, 2.0, 1.0, 1.0])
M2_d4 = sum(n_d4[j] * np.outer(alphas_d4[j], alphas_d4[j]) for j in range(5))
# Look at the eigenvectors of M2_d4 Separated by the Z_2 x Z_2 centre of the affine diagram:
# In d_4^(1):
# - alpha_0 + alpha_1 = -2 e_2, alpha_0 - alpha_1 = -2 e_1
# - alpha_3 + alpha_4 = 2 e_3,  alpha_3 - alpha_4 = -2 e_4
# So the eigenvectors of the Z_2 x Z_2 centre of the affine diagram are:
#   u_1 = e_1 = (1, 0, 0, 0)         (mass sqrt(2), character (-1, +1, +1))
#   u_s = e_4 = (0, 0, 0, 1)         (mass sqrt(2), character (+1, -1, +1))
#   u_c = (e_2 + e_3)/sqrt(2)        (mass sqrt(2), character (+1, +1, -1))
#   u_2 = (e_2 - e_3)/sqrt(2)        (mass sqrt(6), central particle, character (+1, +1, +1))
u1 = np.array([1.0, 0.0, 0.0, 0.0])
us = np.array([0.0, 0.0, 0.0, 1.0])
us_prime = np.array([0.0, 1.0, 1.0, 0.0]) / np.sqrt(2.0)
u2 = np.array([0.0, 1.0, -1.0, 0.0]) / np.sqrt(2.0)
basis_d4 = [("1", u1, np.sqrt(2)), ("s", us, np.sqrt(2)), ("s'", us_prime, np.sqrt(2)), ("2(central)", u2, np.sqrt(6))]

print("\nd_4^(1) (h=6) in the REAL orthonormal Z_2 x Z_2 character eigenbasis:")
for i, (na, va, ma) in enumerate(basis_d4):
    for j in range(i, 4):
        nb, vb, mb = basis_d4[j]
        for k in range(j, 4):
            nc, vc, mc = basis_d4[k]
            C_val = sum(
                n_d4[l]
                * np.dot(alphas_d4[l], va)
                * np.dot(alphas_d4[l], vb)
                * np.dot(alphas_d4[l], vc)
                for l in range(5)
            )
            area_rhs = (4.0 / np.sqrt(6.0)) * triangle_area(ma, mb, mc)
            if abs(C_val) > 1e-10:
                print(
                    f"  |C_{{{na},{nb},{nc}}}| = {abs(C_val):.6f},  (4/sqrt(6))*Delta = {area_rhs:.6f},  diff = {abs(abs(C_val)-area_rhs):.2e}"
                )

# Also check Dorey's rule and Figure 9.1 (ch9_f13_dorey_d4.svg) for D_4:
# In D_4, bicolour: central node 2 is white, outer nodes 1, 3, 4 are black (or vice versa).
# If outer nodes 1, 3, 4 are white and central node 2 is black:
# gamma_1 = alpha_1, gamma_3 = alpha_3, gamma_4 = alpha_4, gamma_2 = s_white(alpha_2) = alpha_2 + alpha_1 + alpha_3 + alpha_4.
# Wait: for ANY choice of bicolouring of D_4:
# Does gamma_1 + gamma_1' + gamma_c = 0 have two roots in the SAME orbit Omega_1, or in TWO DIFFERENT outer orbits (or the same)?
# Let's compute the 4 Coxeter orbits Omega_1, Omega_2, Omega_3, Omega_4 of D_4!
alphas_D4_fin = alphas_d4[1:]  # alpha_1, alpha_2, alpha_3, alpha_4
C_D4 = alphas_D4_fin @ alphas_D4_fin.T
fund_weights_D4 = np.linalg.inv(C_D4) @ alphas_D4_fin


def refl(alpha):
    return np.eye(4) - np.outer(alpha, alpha)


# White = {1, 3, 4} (indices 0, 2, 3), Black = {2} (index 1)
s_white = refl(alphas_D4_fin[0]) @ refl(alphas_D4_fin[2]) @ refl(alphas_D4_fin[3])
s_black = refl(alphas_D4_fin[1])
w_cox = s_black @ s_white
w_inv = np.linalg.inv(w_cox)

orbits_D4 = {}
names_D4 = ["1", "2(c)", "3(s)", "4(s')"]
for idx, name in enumerate(names_D4):
    gam = (np.eye(4) - w_inv) @ fund_weights_D4[idx]
    orb = []
    cur = gam.copy()
    for p in range(6):
        orb.append(np.round(cur, 8))
        cur = w_cox @ cur
    orbits_D4[name] = orb

print("\nD_4 Coxeter orbits (h=6): checking which triples (a, b, c) have alpha + beta + gamma = 0:")
for i, na in enumerate(names_D4):
    for j in range(i, 4):
        nb = names_D4[j]
        for k in range(j, 4):
            nc = names_D4[k]
            found = False
            for r1 in orbits_D4[na]:
                for r2 in orbits_D4[nb]:
                    for r3 in orbits_D4[nc]:
                        if np.linalg.norm(r1 + r2 + r3) < 1e-6:
                            found = True
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                print(f"  Dorey fusing allowed for ({na}, {nb}, {nc})")
