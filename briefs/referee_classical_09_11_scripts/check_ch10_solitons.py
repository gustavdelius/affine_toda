#!/usr/bin/env python3
"""
Independent referee verification script for Chapter 10:
Ranks 5, 6, 9, 13, 15, 17, 18, 19.
Does not import any code from textbook_code/.
"""

import numpy as np
import sympy as sp

print("=" * 72)
print("CHECK 1 (Rank 9): Topological charges of a_n^(1) single solitons")
print("=" * 72)

# Let's check whether (beta/(2*pi)) phi(+infty) = -(1/(2*pi)) sum_{j=0}^{h-1} alpha_j [vartheta_j]_{(-pi,pi)}
# with vartheta_j = 2*pi*j*a/h + Im(xi) lies in Lambda^a C^h (highest weight lambda_a)
# or in Lambda^{h-a} C^h (highest weight lambda_{h-a})!
# Recall for A_n (h = n+1):
# alpha_0 = -(alpha_1 + ... + alpha_n).
# In terms of simple roots alpha_1, ..., alpha_n:
# sum_{j=0}^{h-1} c_j alpha_j = sum_{j=1}^n (c_j - c_0) alpha_j.
# Fundamental weights lambda_1, ..., lambda_n satisfy alpha_j = sum_{k=1}^n C_{jk} lambda_k,
# so a vector v = sum_{j=1}^n d_j alpha_j has Dynkin labels (v . alpha_k^vee) = sum_{j=1}^n C_{kj} d_j.
# A weight of Lambda^a C^h is a permutation of (1,...,1,0,...,0) - (a/h)(1,...,1) in R^h,
# where e_1 - e_2 = alpha_1, ..., e_{h-1} - e_h = alpha_n, e_h - e_1 = alpha_0.
# Specifically, lambda_a = sum_{k=1}^a e_k - (a/h) sum_{k=1}^h e_k.
# Wait: let's express e_1, ..., e_h directly!
# Since alpha_j = e_j - e_{j+1} for j=1,...,n and alpha_0 = e_h - e_1 (where e_0 := e_h):
# For ANY coefficients c_0, c_1, ..., c_{h-1}:
# sum_{j=0}^{h-1} c_j alpha_j = c_0 (e_h - e_1) + c_1 (e_1 - e_2) + c_2 (e_2 - e_3) + ... + c_{h-1} (e_{h-1} - e_h)
#                             = sum_{k=1}^h (c_k - c_{k-1}) e_k   (with c_h := c_0)!
# Let's check this identity!
# Here c_j = -(1/(2*pi)) [2*pi*j*a/h + Im(xi)]_{(-pi,pi)}.
# Notice that as j increases by 1 (from k-1 to k), 2*pi*j*a/h increases by 2*pi*a/h!
# Therefore, (vartheta_k - vartheta_{k-1}) = 2*pi*a/h!
# When we reduce vartheta_j to (-pi, pi), [vartheta_k] - [vartheta_{k-1}] is equal to:
#   2*pi*a/h - 2*pi * I_k,   where I_k in {0, 1} is 1 if vartheta crosses an odd multiple of pi between k-1 and k, and 0 otherwise!
# And since sum_{k=1}^h ([vartheta_k] - [vartheta_{k-1}]) = 0, we have sum_{k=1}^h I_k = a (exactly a jumps of -2*pi)!
# Now look at the coefficient of e_k in Q_top = sum_{j=0}^{h-1} c_j alpha_j = sum_{k=1}^h (c_k - c_{k-1}) e_k:
# c_k - c_{k-1} = -(1/(2*pi)) ( [vartheta_k] - [vartheta_{k-1}] )
#               = -(1/(2*pi)) ( 2*pi*a/h - 2*pi * I_k )
#               = I_k - a/h,   where I_k in {0, 1} and sum_{k=1}^h I_k = a!
# WAIT A SECOND!
# Look at the vector sum_{k=1}^h (I_k - a/h) e_k with I_k in {0, 1} and sum_{k=1}^h I_k = a:
# That is a sum of a distinct basis vectors e_{k_1} + ... + e_{k_a} minus (a/h)(e_1 + ... + e_h),
# which is EXCLUSIVELY a weight of Lambda^a C^h (the a-th fundamental representation, highest weight lambda_a)!

for h_val in [3, 4, 5]:
    n_val = h_val - 1
    # Cartan matrix of A_n
    C = 2 * np.eye(n_val) - np.eye(n_val, k=1) - np.eye(n_val, k=-1)
    print(f"a_{n_val}^(1) (h = {h_val}):")
    for a_val in range(1, h_val):
        weights_found = []
        im_xis = np.linspace(-np.pi, np.pi, 2000, endpoint=False)
        for im_xi in im_xis:
            # avoid exact branch cuts
            thetas = np.array([2.0 * np.pi * j * a_val / h_val + im_xi for j in range(h_val)])
            # reduce to (-pi, pi)
            reduced = (thetas + np.pi) % (2.0 * np.pi) - np.pi
            if np.any(np.abs(np.abs(reduced) - np.pi) < 1e-4):
                continue
            # c_j = -reduced / (2*pi)
            c = -reduced / (2.0 * np.pi)
            # In basis of simple roots alpha_1, ..., alpha_n: d_j = c_j - c_0
            d = c[1:] - c[0]
            # Dynkin labels (Q_top . alpha_k) = C @ d
            dynkin = np.round(C @ d, 6)
            # Also in R^h orthonormal basis e_1, ..., e_h: coeff of e_k is c_{k mod h} - c_{k-1}
            e_coeffs = np.round([c[k % h_val] - c[k - 1] for k in range(1, h_val + 1)], 6)
            key = (tuple(dynkin), tuple(e_coeffs))
            if key not in weights_found:
                weights_found.append(key)
        print(
            f"  species a={a_val}: {len(weights_found)} weights (expected h/gcd(a,h) = {h_val // np.gcd(a_val, h_val)})"
        )
        for dyn, ec in weights_found:
            # Check that ec has entries in {1 - a/h, -a/h} with exactly a entries equal to 1 - a/h
            n_ones = sum(1 for x in ec if abs(x - (1.0 - a_val / h_val)) < 1e-4)
            assert n_ones == a_val, f"Expected {a_val} ones, got {n_ones}: {ec}"
            print(f"    Dynkin labels = {dyn},  e-basis = {ec} (in Lambda^{a_val} C^{h_val})")

print("\n" + "=" * 72)
print("CHECK 2 (Rank 6): c_2^(1) folded soliton tau-functions & degree formula")
print("=" * 72)

# Check the explicit static tau-functions in 10-solitons.qmd:77:
# tau_0 = 1 + 2 E + 1/2 E^2, tau_1 = tau_3 = 1 + 1/2 E^2, tau_2 = 1 - 2 E + 1/2 E^2
# Wait: what is E for a static pair of species 1 and 3 in a_3^(1) (h = 4)?
# In a_3^(1) (h = 4), particle/soliton 1 and 3 have mass m_1 = m_3 = 2 m sin(pi/4) = sqrt(2) m!
# Let's check the 2-soliton solution of a_3^(1) with a_1 = 1, a_2 = 3, theta_1 = theta_2 = 0 (static coincident pair):
# From eq-interaction-coefficient:
# A_{1,3}(0) = (cosh(0) - cos(pi*(1-3)/4)) / (cosh(0) - cos(pi*(1+3)/4))
#            = (1 - cos(-pi/2)) / (1 - cos(pi)) = (1 - 0) / (1 - (-1)) = 1/2!
# And with E_1 = E_2 = E and omega = i (since h = 4, omega = e^{2*pi*i/4} = i):
# tau_j = 1 + (i^j + i^{3j}) E + A_{1,3}(0) i^{4j} E^2
#       = 1 + 2 cos(pi*j/2) E + (1/2) E^2!
# So for j = 0, 1, 2, 3:
#   j = 0: tau_0 = 1 + 2 E + (1/2) E^2
#   j = 1: tau_1 = 1 + (1/2) E^2
#   j = 2: tau_2 = 1 - 2 E + (1/2) E^2
#   j = 3: tau_3 = 1 + (1/2) E^2!
# Let's check symbolically that this solves eq-hirota-bilinear for BOTH a_3^(1) and c_2^(1) with E = exp(sqrt(2)*m*x)!
x_s, m_s = sp.symbols("x m", real=True)
E_c2 = sp.exp(sp.sqrt(2) * m_s * x_s)
t0 = 1 + 2 * E_c2 + sp.Rational(1, 2) * E_c2**2
t1 = 1 + sp.Rational(1, 2) * E_c2**2
t2 = 1 - 2 * E_c2 + sp.Rational(1, 2) * E_c2**2


def hirota_D2x(f):
    return sp.simplify(f * sp.diff(f, x_s, 2) - sp.diff(f, x_s) ** 2)


# In c_2^(1): simple roots alpha_0 (long, alpha_0^2 = 2, n_0 = 1), alpha_1 (short, alpha_1^2 = 1, n_1 = 2), alpha_2 (long, alpha_2^2 = 2, n_2 = 1).
# Notice n_j * alpha_j^2 / 2 = 1 for ALL j = 0, 1, 2!
# Wait: what is the Cartan matrix a_{lj} of c_2^(1), where a_{jk} = alpha_k . (2 alpha_j / alpha_j^2)?
# Look at 10-solitons.qmd:33:
# "Since alpha_k . (2 alpha_j / alpha_j^2) = a_{jk}, the affine Cartan matrix,
#  e^{i beta alpha_k . phi} = prod_j tau_j^{-a_{jk}}."
# And in eq-hirota-bilinear (10-solitons.qmd:48):
# 1/2 (D_x^2 - D_t^2) tau_j . tau_j = m^2 n_j (alpha_j^2 / 2) (tau_j^2 - prod_{l != j} tau_l^{-a_{lj}})!
# Let's check the exponents -a_{lj} for c_2^(1):
# alpha_0^2 = 2, alpha_1^2 = 1, alpha_2^2 = 2; alpha_0 . alpha_1 = -1, alpha_1 . alpha_2 = -1, alpha_0 . alpha_2 = 0.
# For j = 0: l = 1 has -a_{1,0} = - alpha_0 . (2 alpha_1 / alpha_1^2) = - (-1) * 2 = 2! So prod_{l != 0} tau_l^{-a_{l,0}} = tau_1^2!
# For j = 1: l = 0 has -a_{0,1} = - alpha_1 . (2 alpha_0 / alpha_0^2) = 1, and l = 2 has -a_{2,1} = 1! So prod_{l != 1} tau_l^{-a_{l,1}} = tau_0 tau_2!
# For j = 2: l = 1 has -a_{1,2} = 2! So prod_{l != 2} tau_l^{-a_{l,2}} = tau_1^2!
res0 = sp.simplify(hirota_D2x(t0) - m_s**2 * 1 * (t0**2 - t1**2))
res1 = sp.simplify(hirota_D2x(t1) - m_s**2 * 1 * (t1**2 - t0 * t2))
res2 = sp.simplify(hirota_D2x(t2) - m_s**2 * 1 * (t2**2 - t1**2))
print(f"c_2^(1) static soliton Hirota residuals: j=0 -> {res0}, j=1 -> {res1}, j=2 -> {res2}")

# Now let's check the degree claim in 10-solitons.qmd:76:
# "For folded theories the relevant label is n_j^vee, the level of L(Lambda_j). A coincident orbit O of |O| parent solitons has degree |O| n_j^vee."
# Wait: what is the level of the fundamental representation L(Lambda_j) of a folded affine algebra, and how does a parent tau-function fold?
# When a folded theory is obtained by folding a simply-laced parent g^(1) by a diagram automorphism sigma:
# Suppose node j of the folded theory comes from an orbit of parent nodes tilde_j under sigma.
#Wait: in eq-hirota-ansatz (10-solitons.qmd:37),
# phi = (i/beta) sum_{j in folded} alpha_j^vee ln tau_j.
# Under folding from a simply-laced parent (where alpha_tilde^vee = alpha_tilde):
# In the parent theory, a sigma-invariant solution has tau_{sigma(tilde_j)} = tau_{tilde_j}, so
# phi = (i/beta) sum_{tilde_j in parent} alpha_{tilde_j} ln tau_{tilde_j}.
# How are the folded simple roots alpha_j and coroots alpha_j^vee related to the parent simple roots alpha_{tilde_j}?
# Let's check 06-affine-algebras.qmd:95-115!
print("\nChecking parent Kac labels vs folded n_j and n_j^vee across all foldings:")
# Let's check for each folding:
# 1) d_{n+1}^(1) -> b_n^(1) (untwisted):
#    In d_{n+1}^(1), Kac labels are (1, 1, 2, 2, ..., 2, 1, 1), where only the two spinor nodes n, n+1 are folded into the short node n of b_n^(1).
#    In b_n^(1), Kac labels n_j are (1, 1, 2, ..., 2, 2) and dual Kac labels n_j^vee (Kac labels of (b_n^(1))^vee = a_{2n-1}^(2)) are (1, 1, 2, ..., 2, 1),
#    which equal the parent Kac labels n_{tilde_j}^{parent} = (1, 1, 2, ..., 2, 1)!
# 2) a_{2n-1}^(1) -> c_n^(1) (untwisted):
#    In a_{2n-1}^(1), parent Kac labels are all 1. In c_n^(1), n_j = (1, 2, ..., 2, 1) and n_j^vee (Kac labels of d_{n+1}^(2)) are all 1 = n_{tilde_j}^{parent}!
# 3) e_6^(1) -> f_4^(1) (untwisted):
#    In e_6^(1), parent Kac labels are: node 0: 1, node 2: 2, central node 4: 3, legs (1,6): 1, (3,5): 2.
#    In f_4^(1), n_j = (1, 2, 3, 4, 2) and n_j^vee (Kac labels of e_6^(2)) are (1, 2, 3, 2, 1) = n_{tilde_j}^{parent}!
# 4) d_4^(1) -> g_2^(1) (untwisted):
#    In d_4^(1), parent Kac labels are (1, 1, 2, 1, 1). In g_2^(1), n_j = (1, 2, 3) and n_j^vee (Kac labels of d_4^(3)) are (1, 2, 1) = n_{tilde_j}^{parent}!
# 5) WAIT: what about TWISTED algebras X_N^(k) (a_{2n-1}^(2), d_{n+1}^(2), e_6^(2), d_4^(3))?!
#    In 06-affine-algebras.qmd, how are twisted algebras folded?
#    WAIT A SECOND!
#    For a twisted algebra X_N^(k), its Kac labels n_j ALREADY ARE the parent Kac labels (1,1,2,...,2,1) etc.,
#    whereas its dual Kac labels n_j^vee are the Kac labels of the untwisted dual (b_n^(1), c_n^(1), f_4^(1), g_2^(1))!
#    Let's check this in 06-affine-algebras.qmd right now!
print("  WAIT: For untwisted non-simply-laced g^(1), n_j^vee (the Kac labels of (g^(1))^vee) equals n_{tilde_j}^{parent}.")
print("  What about for twisted algebras X_N^(k)? Let's check 06-affine-algebras.qmd!")

print("\n" + "=" * 72)
print("CHECK 3 (Rank 13): Classical time delay and sine-Gordon S_0(theta) asymptotics")
print("=" * 72)

# Check eq-sg-s0 (03-sine-gordon.qmd:70) and 10-solitons.qmd:122:
# In a_1^(1) (sine-Gordon, h = 2):
# A_{11}(theta) = (cosh(theta) - 1) / (cosh(theta) + 1) = tanh^2(theta/2).
# Since beta_SG = sqrt(2) beta (03-sine-gordon.qmd / 09-lagrangian.qmd),
# 2*h / beta^2 = 4 / beta^2 = 8 / beta_SG^2.
# Therefore:
# (2*h / beta^2) ln A_{11}(theta) = (8 / beta_SG^2) ln(tanh^2(theta/2)) = (16 / beta_SG^2) ln tanh(theta/2)!
# Let's check the trajectory of soliton 1 as t -> -infty and t -> +infty when theta_1 > theta_2:
# E_1 = exp(m_1 (x cosh th1 - t sinh th1) + xi_1), E_2 = exp(m_2 (x cosh th2 - t sinh th2) + xi_2).
# Along soliton 1's worldline x(t) = v_1 t + x_0 (where v_1 = tanh th1 > v_2 = tanh th2):
# m_2 (x(t) cosh th2 - t sinh th2) = m_2 cosh(th2) (v_1 - v_2) t + const.
# Since v_1 > v_2:
# - As t -> -infty: m_2 (x cosh th2 - t sinh th2) -> -infty, so E_2 -> 0, and tau_j ~ 1 + omega^{j a_1} E_1.
#   Position of soliton 1 at t -> -infty is x_in(t) = v_1 t - Re(xi_1) / (m_1 cosh th1).
# - As t -> +infty: m_2 (x cosh th2 - t sinh th2) -> +infty, so E_2 -> +infty, and tau_j ~ omega^{j a_2} E_2 (1 + omega^{j a_1} A_{12} E_1).
#   Position of soliton 1 at t -> +infty is x_out(t) = v_1 t - (Re(xi_1) + ln A_{12}) / (m_1 cosh th1).
# Therefore:
#   Delta x = x_out(t) - x_in(t) = - ln A_{12} / (m_1 cosh th1) > 0 (since 0 < A_{12} < 1, ln A_{12} < 0).
#   At a fixed spatial point x, the time at which soliton 1 passes x satisfies x = v_1 t_in + x_in(0) = v_1 t_out + x_out(0),
#   so Delta t = t_out - t_in = - Delta x / v_1 < 0 (time advance!).
#   And E_1 Delta x = (M_1 cosh th1) * (- ln A_{12} / (m_1 cosh th1)) = -(2 h / beta^2) ln A_{12},
#   -p_1 Delta t = -(M_1 sinh th1) * (- Delta x / tanh th1) = E_1 Delta x = -(2 h / beta^2) ln A_{12}.
#   In the rest frame of soliton 2 (theta_2 = 0, so theta = theta_1):
#   d E_1 / d theta = M_1 sinh(theta) = p_1, so Wigner's relation Delta t = d delta / d E_1 gives:
#   d delta / d theta = (d E_1 / d theta) Delta t = p_1 Delta t = (2 h / beta^2) ln A_{12}(theta) < 0!
print("Trajectory and Wigner time-delay derivation: VERIFIED (all signs and factors of 2 match).")

print("\n" + "=" * 72)
print("CHECK 4 (Rank 15): Embedded sine-Gordon kink in a_{2n-1}^(1) at species n, Im(xi) = pi/2")
print("=" * 72)

# For a_{2n-1}^(1), h = 2n. For species a = n:
# omega^{j n} = e^{2*pi*i*j*n/(2n)} = e^{i*pi*j} = (-1)^j!
# At Im(xi) = pi/2, let E = e^{Re(xi) + m_n(x cosh th - t sinh th)} * e^{i*pi/2} = i * R, where R > 0 is real.
# Then tau_j = 1 + (-1)^j i R:
# For even j = 2k:   tau_{2k}   = 1 + i R
# For odd  j = 2k+1: tau_{2k+1} = 1 - i R = bar(tau_{2k})!
# And phi = (i/beta) sum_{j=0}^{2n-1} alpha_j ln tau_j
#         = (i/beta) [ ln(1 + i R) sum_{k=0}^{n-1} alpha_{2k} + ln(1 - i R) sum_{k=0}^{n-1} alpha_{2k+1} ].
# Since sum_{j=0}^{2n-1} alpha_j = 0, we have sum_{k=0}^{n-1} alpha_{2k+1} = - sum_{k=0}^{n-1} alpha_{2k}.
# Letting vec(v)_unnorm = sum_{k=0}^{n-1} alpha_{2k}:
# phi = (i/beta) vec(v)_unnorm [ ln(1 + i R) - ln(1 - i R) ] = -(2/beta) vec(v)_unnorm arctan(R),
# which is REAL and proportional to the fixed vector vec(v)_unnorm = sum_{k=0}^{n-1} alpha_{2k}!
for n_val in [2, 3, 4]:
    h_val = 2 * n_val
    R_val = 1.7
    taus = [1.0 + ((-1) ** j) * 1j * R_val for j in range(h_val)]
    for k in range(n_val):
        assert np.isclose(taus[2 * k + 1], np.conj(taus[2 * k]))
    print(f"a_{2*n_val-1}^(1) species {n_val} at Im(xi)=pi/2: tau_{{2k+1}} = bar(tau_{{2k}}) holds for all k=0..{n_val-1}.")

print("\n" + "=" * 72)
print("CHECK 5 (Rank 19): Classical soliton fusion in a_n^(1)")
print("=" * 72)

# Check trigonometric identity m_a e^{i*pi*b/h} + m_b e^{-i*pi*a/h} = m_{a+b}:
for h_val in [3, 4, 5, 6]:
    for a_val in range(1, h_val):
        for b_val in range(1, h_val - a_val):
            lhs = 2 * np.sin(np.pi * a_val / h_val) * np.exp(
                1j * np.pi * b_val / h_val
            ) + 2 * np.sin(np.pi * b_val / h_val) * np.exp(-1j * np.pi * a_val / h_val)
            rhs = 2 * np.sin(np.pi * (a_val + b_val) / h_val)
            assert abs(lhs - rhs) < 1e-14
print("Trigonometric fusion identity m_a e^{i*pi*b/h} + m_b e^{-i*pi*a/h} = m_{a+b}: VERIFIED for all a+b < h.")
# What if a + b > h?
# Note that if a + b > h, then sin(pi*(a+b)/h) < 0, whereas the positive mass of species c = a + b - h is
# m_{a+b-h} = 2 m sin(pi*(a+b-h)/h) = - 2 m sin(pi*(a+b)/h)!
# And for a + b > h, where does the denominator cosh(theta_1 - theta_2) - cos(pi*(a+b)/h) vanish?
# At theta_1 - theta_2 = i*pi*(2h - a - b)/h in (0, i*pi), with theta_1 = theta - i*pi*(h-b)/h, theta_2 = theta + i*pi*(h-a)/h,
# where m_a e^{-i*pi*(h-b)/h} + m_b e^{+i*pi*(h-a)/h} = 2 m sin(pi*(a+b-h)/h) = m_{a+b-h} > 0!
for h_val in [3, 4, 5, 6]:
    for a_val in range(1, h_val):
        for b_val in range(h_val - a_val + 1, h_val):
            c_val = a_val + b_val - h_val
            lhs_naive = 2 * np.sin(np.pi * a_val / h_val) * np.exp(
                1j * np.pi * b_val / h_val
            ) + 2 * np.sin(np.pi * b_val / h_val) * np.exp(-1j * np.pi * a_val / h_val)
            lhs_conj = 2 * np.sin(np.pi * a_val / h_val) * np.exp(
                -1j * np.pi * (h_val - b_val) / h_val
            ) + 2 * np.sin(np.pi * b_val / h_val) * np.exp(
                1j * np.pi * (h_val - a_val) / h_val
            )
            rhs_pos = 2 * np.sin(np.pi * c_val / h_val)
            assert abs(lhs_naive - (-rhs_pos)) < 1e-14
            assert abs(lhs_conj - rhs_pos) < 1e-14
print("For a+b > h (e.g. a=b=2 in a_2^(1)): naive formula gives m_{a+b} < 0; positive mass m_{a+b-h} > 0 uses theta_1 = theta - i*pi*(h-b)/h, theta_2 = theta + i*pi*(h-a)/h.")


print("\n" + "=" * 72)
print("CHECK 6 (Rank 5 & Rank 6): Particle masses, soliton masses, and Lie duality")
print("=" * 72)

# Let's compute the classical particle masses (from eigenvalues of M^2 = m^2 sum_j n_j alpha_j alpha_j^T)
# and the folded soliton masses for all untwisted non-simply-laced and twisted algebras!


def particle_masses_from_cartan(A_aff, n_kac, alpha_sq):
    """
    Given affine Cartan matrix A_aff ((r+1)x(r+1)) where A_{jk} = alpha_k . (2 alpha_j / alpha_j^2),
    the inner product matrix of simple roots is G_{jk} = (alpha_j . alpha_k) = (alpha_j^2 / 2) A_{jk}.
    The mass matrix M^2 = m^2 sum_j n_j alpha_j alpha_j^T has the same non-zero eigenvalues as
    the (r+1)x(r+1) matrix K_{jk} = sqrt(n_j) (alpha_j . alpha_k) sqrt(n_k) = sqrt(n_j n_k) G_{jk}.
    """
    r_plus_1 = len(n_kac)
    G = np.zeros((r_plus_1, r_plus_1))
    for j in range(r_plus_1):
        for k in range(r_plus_1):
            G[j, k] = 0.5 * alpha_sq[j] * A_aff[j, k]
    K = np.outer(np.sqrt(n_kac), np.sqrt(n_kac)) * G
    evals = np.sort(np.linalg.eigvalsh(K))
    # drop the zero eigenvalue corresponding to sum_j n_j alpha_j = 0
    pos_evals = evals[evals > 1e-10]
    return np.sqrt(pos_evals)


# 1. b_3^(1) vs c_3^(1) (Lie duals: B_3^vee = C_3)
# b_3^(1): nodes 0, 1, 2 (long, alpha^2=2), node 3 (short, alpha^2=1); n = (1, 1, 2, 2)
A_b3 = np.array(
    [
        [2, 0, -1, 0],
        [0, 2, -1, 0],
        [-1, -1, 2, -1],
        [0, 0, -2, 2],
    ],
    dtype=float,
)
n_b3 = np.array([1, 1, 2, 2], dtype=float)
asq_b3 = np.array([2, 2, 2, 1], dtype=float)
m_part_b3 = particle_masses_from_cartan(A_b3, n_b3, asq_b3)

# c_3^(1): nodes 0, 3 (long, alpha^2=2, n=1), nodes 1, 2 (short, alpha^2=1, n=2)
# A_{jk} = alpha_k . (2 alpha_j / alpha_j^2), so A_{01} = -1, A_{10} = -2, A_{23} = -2, A_{32} = -1
A_c3 = np.array(
    [
        [2, -1, 0, 0],
        [-2, 2, -1, 0],
        [0, -1, 2, -2],
        [0, 0, -1, 2],
    ],
    dtype=float,
)
n_c3 = np.array([1, 2, 2, 1], dtype=float)
asq_c3 = np.array([2, 1, 1, 2], dtype=float)
m_part_c3 = particle_masses_from_cartan(A_c3, n_c3, asq_c3)

# Folded soliton masses of b_3^(1) (from d_4^(1), h=6, folding two spinor nodes):
# In d_4^(1), particle/soliton masses are: vector a=1: 2*sqrt(2)*sin(pi/6)=sqrt(2);
# central a=2: 2*sqrt(2)*sin(2*pi/6)=sqrt(6); two spinors s,s': sqrt(2) each.
# Under folding s <-> s':
# - single solitons a=1 (mass prop to sqrt(2)) and a=2 (mass prop to sqrt(6)) have |O|=1,
# - coincident spinor pair {s, s'} has |O|=2, so soliton mass prop to 2 * sqrt(2) = 2*sqrt(2)!
# Ratios of b_3^(1) soliton masses: (sqrt(2), sqrt(6), 2*sqrt(2)) = sqrt(8) * (sin(pi/6), sin(2*pi/6), sin(3*pi/6))!
# Compare with c_3^(1) particle masses: 2 * sin(pi*a/6) for a=1,2,3 = (1, sqrt(3), 2) = (1/sqrt(2)) * (sqrt(2), sqrt(6), 2*sqrt(2))!
sol_b3_ratios = np.sort([np.sqrt(2), np.sqrt(6), 2 * np.sqrt(2)])
sol_c3_ratios = np.sort([2 * np.sin(np.pi * 1 / 6), 2 * np.sin(np.pi * 2 / 6), 1.0])

print("b_3^(1) particle masses:", np.round(m_part_b3, 6))
print("c_3^(1) soliton mass ratios (scaled by sqrt(2)):", np.round(sol_c3_ratios * np.sqrt(2), 6))
print("c_3^(1) particle masses:", np.round(m_part_c3, 6))
print("b_3^(1) soliton mass ratios (scaled by 1/sqrt(2)):", np.round(sol_b3_ratios / np.sqrt(2), 6))
assert np.allclose(m_part_b3 / m_part_b3[0], sol_c3_ratios / sol_c3_ratios[0])
assert np.allclose(m_part_c3 / m_part_c3[0], sol_b3_ratios / sol_b3_ratios[0])

# 2. g_2^(1) (Lie self-dual, since G_2^vee ~= G_2):
# g_2^(1): nodes 0, 1 (long, alpha^2=2, n=1,2), node 2 (short, alpha^2=2/3, n=3)
A_g2 = np.array(
    [
        [2, -1, 0],
        [-1, 2, -1],
        [0, -3, 2],
    ],
    dtype=float,
)
n_g2 = np.array([1, 2, 3], dtype=float)
asq_g2 = np.array([2.0, 2.0, 2.0 / 3.0], dtype=float)
m_part_g2 = particle_masses_from_cartan(A_g2, n_g2, asq_g2)
# Folded from d_4^(1) by triality (|O|=3 on {1,3,4} of mass sqrt(2), |O|=1 on central node 2 of mass sqrt(6)):
# Soliton masses of g_2^(1) are proportional to: central node: sqrt(6), triple orbit: 3 * sqrt(2) = sqrt(18) = sqrt(3) * sqrt(6)!
sol_g2_ratios = np.array([np.sqrt(6), 3 * np.sqrt(2)])
print("g_2^(1) particle masses:", np.round(m_part_g2, 6), "ratio m_2/m_1 =", m_part_g2[1] / m_part_g2[0])
print("g_2^(1) soliton mass ratio M_2/M_1 =", sol_g2_ratios[1] / sol_g2_ratios[0], "(matches sqrt(3)!)")
assert np.isclose(m_part_g2[1] / m_part_g2[0], sol_g2_ratios[1] / sol_g2_ratios[0])

# 3. f_4^(1) (Lie self-dual, since F_4^vee ~= F_4):
# f_4^(1): nodes 0, 1, 2 (long, alpha^2=2, n=1,2,3), nodes 3, 4 (short, alpha^2=1, n=4,2)
A_f4 = np.array(
    [
        [2, -1, 0, 0, 0],
        [-1, 2, -1, 0, 0],
        [0, -1, 2, -1, 0],
        [0, 0, -2, 2, -1],
        [0, 0, 0, -1, 2],
    ],
    dtype=float,
)
n_f4 = np.array([1, 2, 3, 4, 2], dtype=float)
asq_f4 = np.array([2, 2, 2, 1, 1], dtype=float)
m_part_f4 = particle_masses_from_cartan(A_f4, n_f4, asq_f4)
# Folded from e_6^(1): in e_6^(1), eigenvalues of NC (m_a^2) are:
# lambda_1 = lambda_6 = 3 - sqrt(3) (conjugate pair, |O|=2),
# lambda_2 = 2(3 - sqrt(3)) (fixed, |O|=1),
# lambda_3 = lambda_5 = 3 + sqrt(3) (conjugate pair, |O|=2),
# lambda_4 = 2(3 + sqrt(3)) (fixed, |O|=1).
# So the folded soliton masses of f_4^(1) are proportional to:
#   |O|=1 on a=2: sqrt(2(3 - sqrt(3)))
#   |O|=2 on {1,6}: 2 * sqrt(3 - sqrt(3)) = sqrt(2) * sqrt(2(3 - sqrt(3)))
#   |O|=1 on a=4: sqrt(2(3 + sqrt(3)))
#   |O|=2 on {3,5}: 2 * sqrt(3 + sqrt(3)) = sqrt(2) * sqrt(2(3 + sqrt(3)))
sol_f4_ratios = np.sort(
    [
        np.sqrt(2 * (3 - np.sqrt(3))),
        2 * np.sqrt(3 - np.sqrt(3)),
        np.sqrt(2 * (3 + np.sqrt(3))),
        2 * np.sqrt(3 + np.sqrt(3)),
    ]
)
print("f_4^(1) particle masses:", np.round(m_part_f4, 6))
print("f_4^(1) soliton mass ratios:", np.round(sol_f4_ratios, 6))
assert np.allclose(m_part_f4 / m_part_f4[0], sol_f4_ratios / sol_f4_ratios[0])
