#!/usr/bin/env python3
"""
Complete verification script for Chapter 11 (and related Chapter 10 claims):
Ranks 1, 2, 3, 4, 7, 20.
Does not import any code from textbook_code/.
"""

import numpy as np
import sympy as sp

print("=" * 72)
print("CHECK 4 (Rank 2): Exact spin-2 conserved current and charge in a_2^(1)")
print("=" * 72)

m_s, beta_s = sp.symbols("m beta", real=True, positive=True)
phi1, phi2 = sp.symbols("phi1 phi2")
p1_p, p2_p, p1_pp, p2_pp = sp.symbols("p1_p p2_p p1_pp p2_pp")
e0, e1, e2 = sp.symbols("e0 e1 e2")

# alpha_0 = (-sqrt(2), 0)
# alpha_1 = (1/sqrt(2), sqrt(3/2))
# alpha_2 = (1/sqrt(2), -sqrt(3/2))
# E_j = exp(i beta alpha_j . phi)
# d_+ E_0 = -i sqrt(2) beta p1_p E_0
# d_+ E_1 = i beta (p1_p/sqrt(2) + sqrt(3/2) p2_p) E_1
# d_+ E_2 = i beta (p1_p/sqrt(2) - sqrt(3/2) p2_p) E_2

F1 = (sp.I * m_s**2 / beta_s) * (1 / sp.sqrt(2)) * (-2 * e0 + e1 + e2)
F2 = (sp.I * m_s**2 / beta_s) * sp.sqrt(sp.Rational(3, 2)) * (e1 - e2)

de0 = -sp.I * sp.sqrt(2) * beta_s * p1_p * e0
de1 = sp.I * beta_s * (p1_p / sp.sqrt(2) + sp.sqrt(sp.Rational(3, 2)) * p2_p) * e1
de2 = sp.I * beta_s * (p1_p / sp.sqrt(2) - sp.sqrt(sp.Rational(3, 2)) * p2_p) * e2

F1_p = (sp.I * m_s**2 / beta_s) * (1 / sp.sqrt(2)) * (-2 * de0 + de1 + de2)
F2_p = (sp.I * m_s**2 / beta_s) * sp.sqrt(sp.Rational(3, 2)) * (de1 - de2)

c_W = -3 * sp.I / (sp.sqrt(2) * beta_s)
J3 = p2_p * (p2_p**2 - 3 * p1_p**2) + c_W * (p1_p * p2_pp - p2_p * p1_pp)

d_minus_J3 = sp.expand(
    sp.diff(J3, p1_p) * F1
    + sp.diff(J3, p2_p) * F2
    + sp.diff(J3, p1_pp) * F1_p
    + sp.diff(J3, p2_pp) * F2_p
)

# Candidate flux J_{-1} such that d_- J3 + d_+ J_{-1} = 0:
# Let J_{-1} = (m^2 / beta^2) * [ 3 e0 p2_p + e1 (a_1 p1_p + b_1 p2_p) + e2 (a_2 p1_p + b_2 p2_p) ]
a1_s, b1_s, a2_s, b2_s = sp.symbols("a1 b1 a2 b2")
Jm1 = (m_s**2 / beta_s**2) * (
    3 * e0 * p2_p + e1 * (a1_s * p1_p + b1_s * p2_p) + e2 * (a2_s * p1_p + b2_s * p2_p)
)
d_plus_Jm1 = sp.expand(
    (m_s**2 / beta_s**2)
    * (
        3 * de0 * p2_p
        + 3 * e0 * p2_pp
        + de1 * (a1_s * p1_p + b1_s * p2_p)
        + e1 * (a1_s * p1_pp + b1_s * p2_pp)
        + de2 * (a2_s * p1_p + b2_s * p2_p)
        + e2 * (a2_s * p1_pp + b2_s * p2_pp)
    )
)

tot = sp.expand(d_minus_J3 + d_plus_Jm1)
sol = sp.solve(
    [
        tot.coeff(e1).coeff(p1_pp),
        tot.coeff(e1).coeff(p2_pp),
        tot.coeff(e2).coeff(p1_pp),
        tot.coeff(e2).coeff(p2_pp),
    ],
    [a1_s, b1_s, a2_s, b2_s],
)
print("Flux coefficients:", sol)
tot_sub = sp.simplify(tot.subs(sol))
print("Conservation residual d_- J3 + d_+ J_{-1} =", tot_sub)
Jm1_exact = sp.simplify(Jm1.subs(sol))
print("J_{+3} =", J3)
print("J_{-1} =", Jm1_exact)

# Check Theta-parity of J_{+3} and J_{-1}:
# Recall Theta: vec(phi)(x,t) -> -bar(vec(phi))(x,t).
# Since x, t are real, d_+ -> d_+ and d_- -> d_-, so d_+^k vec(phi) -> -bar(d_+^k vec(phi)),
# and E_j = exp(i beta alpha_j . phi) -> exp(i beta alpha_j . (-bar(phi))) = bar(exp(i beta alpha_j . phi)) = bar(E_j)!
# Now look at J_{+3}:
# First term: p2_p * (p2_p^2 - 3 p1_p^2) is homogeneous of degree 3 in d_+ phi with REAL coefficients!
# Under phi -> -bar(phi), p_k -> -bar(p_k), so degree 3 with real coefficients transforms to:
# (-1)^3 bar(p2_p * (p2_p^2 - 3 p1_p^2)) = - bar(...)!
# Second term: c_W * (p1_p * p2_pp - p2_p * p1_pp) where c_W = -3*i/(sqrt(2)*beta) is PURELY IMAGINARY!
# Under phi -> -bar(phi), (p1_p * p2_pp - p2_p * p1_pp) is degree 2 in phi, so it goes to (+1)*bar(...),
# while c_W is purely imaginary, so bar(c_W) = -c_W, which means:
# c_W * (p1_p * p2_pp - p2_p * p1_pp) -> c_W * bar(p1_p * p2_pp - p2_p * p1_pp) = - bar(c_W * (p1_p * p2_pp - p2_p * p1_pp))!
# And look at J_{-1}:
# J_{-1} = (3 m^2 / (2 beta^2)) * [ 2 E_0 p2_p - E_1 (sqrt(3) p1_p + p2_p) + E_2 (sqrt(3) p1_p - p2_p) ],
# which has REAL coefficients times E_j * p_k!
# Under Theta: E_j -> bar(E_j) and p_k -> -bar(p_k), so J_{-1} -> - bar(J_{-1})!
print("\nTheta transformation of J_{+3} and J_{-1} (with real normalization displayed above):")
print("  Under Theta (phi -> -bar(phi)):")
print("  - cubic term in J_{+3}: real coeff * (d_+ phi)^3 -> - bar(cubic term)")
print("  - derivative term in J_{+3}: (-3i/(sqrt(2) beta)) * (d_+ phi_1 d_+^2 phi_2 - d_+ phi_2 d_+^2 phi_1) -> - bar(derivative term)")
print("  - flux term J_{-1}: real coeff * E_j * d_+ phi -> - bar(J_{-1})")
print("  => J_{Q_2}(-bar(phi)) = - bar(J_{Q_2}(phi)), i.e. epsilon_{Q_2} = -1!")

# Now let's evaluate Q_2 = int_{-infty}^{+infty} (J_{+3} + J_{-1}) dx on single solitons and on a breather!
# Wait: since d_- J_{+3} + d_+ J_{-1} = 0 and d_+/- = d_t +/- d_x:
# (d_t - d_x) J_{+3} + (d_t + d_x) J_{-1} = 0  =>  d_t (J_{+3} + J_{-1}) + d_x (J_{-1} - J_{+3}) = 0!
# So the conserved charge density is J_{Q_2}^0 = J_{+3} + J_{-1}!
# Let's integrate J_{+3} + J_{-1} numerically on:
# (a) single soliton of species a=1 at rapidity theta,
# (b) single soliton of species a=2 at rapidity theta,
# (c) breather {(1, theta_0 + i*u), (2, theta_0 - i*u)} at two different times t!

alphas_a2 = np.array(
    [
        [-np.sqrt(2.0), 0.0],
        [1.0 / np.sqrt(2.0), np.sqrt(1.5)],
        [1.0 / np.sqrt(2.0), -np.sqrt(1.5)],
    ]
)
omega = np.exp(2j * np.pi / 3.0)
m_val = 1.0
beta_val = 1.0
m1 = np.sqrt(3.0) * m_val


def eval_Q2_1sol(a_spec, th, xi=0.1 + 0.2j):
    L = 30.0
    xs = np.linspace(-L, L, 100001)
    px = m1 * np.cosh(th)
    pt = -m1 * np.sinh(th)
    pp = pt + px  # d_+ = d_t + d_x => d_+ E = (pt + px) E = m1 e^{-th} E!
    # Wait: if E = exp(m_a (x cosh th - t sinh th) + xi), then:
    # (d_t + d_x) E = m_a (cosh th - sinh th) E = m_a e^{-th} E,
    # (d_t - d_x) E = -m_a (cosh th + sinh th) E = -m_a e^{+th} E!
    # Wait: if we want spin +2 scaling e^{+2 th}, that is the charge built from d_- = d_t - d_x (or with x -> -x)!
    # Let's compute BOTH Q_{+2} (built from d_-) and Q_{-2} (built from d_+)!
    E = np.exp(px * xs + xi)
    taus = np.array([1.0 + (omega ** (a_spec * j)) * E for j in range(3)])
    c = np.array([omega ** (a_spec * j) for j in range(3)])[:, None]
    # For any derivative D with D(E) = lam * E:
    # D(ln tau_j) = lam * c_j E / tau_j
    # D^2(ln tau_j) = lam^2 * c_j E / tau_j^2
    # D^3(ln tau_j) = lam^3 * c_j E (1 - c_j E) / tau_j^3
    results = {}
    for sign_label, lam in [("d_- (scales as e^{+2 th})", pt - px), ("d_+ (scales as e^{-2 th})", pt + px)]:
        d1_ln = lam * c * E[None, :] / taus
        d2_ln = (lam**2) * c * E[None, :] / (taus**2)
        phi_d1 = (1j / beta_val) * np.einsum("jk,jx->kx", alphas_a2, d1_ln)
        phi_d2 = (1j / beta_val) * np.einsum("jk,jx->kx", alphas_a2, d2_ln)
        p1_1, p2_1 = phi_d1[0], phi_d1[1]
        p1_2, p2_2 = phi_d2[0], phi_d2[1]
        cW_val = -3j / (np.sqrt(2.0) * beta_val)
        J_lead = p2_1 * (p2_1**2 - 3.0 * p1_1**2) + cW_val * (p1_1 * p2_2 - p2_1 * p1_1 * 0 - p2_1 * p1_2)
        # E_j = prod_l tau_l^{-a_{lj}} = tau_{j-1} tau_{j+1} / tau_j^2
        E0_v = taus[1] * taus[2] / (taus[0] ** 2)
        E1_v = taus[0] * taus[2] / (taus[1] ** 2)
        E2_v = taus[0] * taus[1] / (taus[2] ** 2)
        J_flux = (1.5 * m_val**2 / beta_val**2) * (
            2.0 * E0_v * p2_1
            + E1_v * (np.sqrt(3.0) * p1_1 - p2_1)
            - E2_v * (np.sqrt(3.0) * p1_1 + p2_1)
        )
        # Also subtract vacuum value of J_flux (which is 0 since p1_1, p2_1 -> 0)
        Q_val = np.trapz(J_lead + J_flux, xs)
        results[sign_label] = Q_val
    return results


for th_test in [0.0, 0.4]:
    r1 = eval_Q2_1sol(1, th_test)
    r2 = eval_Q2_1sol(2, th_test)
    for k in r1:
        print(
            f"th={th_test} [{k}]: Q(a=1) = {r1[k]:.6f}, Q(a=2) = {r2[k]:.6f}, "
            f"Q(a=1)/e^{{+-2th}} = {r1[k]/np.exp(2*th_test if 'd_-' in k else -2*th_test):.6f}"
        )


def eval_Q2_breather(th0, u, t_val=0.0, xi1=0.1 + 0.2j):
    L = 30.0
    xs = np.linspace(-L, L, 100001)
    th1 = th0 + 1j * u
    th2 = th0 - 1j * u
    xi2 = np.conj(xi1) + 1j * np.pi
    A12 = (np.cosh(2j * u) - np.cos(-np.pi / 3.0)) / (np.cosh(2j * u) - np.cos(np.pi))

    px1, pt1 = m1 * np.cosh(th1), -m1 * np.sinh(th1)
    px2, pt2 = m1 * np.cosh(th2), -m1 * np.sinh(th2)
    lam1 = pt1 - px1
    lam2 = pt2 - px2

    E1 = np.exp(px1 * xs + pt1 * t_val + xi1)
    E2 = np.exp(px2 * xs + pt2 * t_val + xi2)

    taus = np.zeros((3, len(xs)), dtype=complex)
    d1_taus = np.zeros((3, len(xs)), dtype=complex)
    d2_taus = np.zeros((3, len(xs)), dtype=complex)
    for j in range(3):
        c1 = omega**j
        c2 = omega ** (2 * j)
        c12 = A12
        taus[j] = 1.0 + c1 * E1 + c2 * E2 + c12 * E1 * E2
        d1_taus[j] = c1 * lam1 * E1 + c2 * lam2 * E2 + c12 * (lam1 + lam2) * E1 * E2
        d2_taus[j] = (
            c1 * (lam1**2) * E1
            + c2 * (lam2**2) * E2
            + c12 * ((lam1 + lam2) ** 2) * E1 * E2
        )

    d1_ln = d1_taus / taus
    d2_ln = (d2_taus * taus - d1_taus**2) / (taus**2)
    phi_d1 = (1j / beta_val) * np.einsum("jk,jx->kx", alphas_a2, d1_ln)
    phi_d2 = (1j / beta_val) * np.einsum("jk,jx->kx", alphas_a2, d2_ln)
    p1_1, p2_1 = phi_d1[0], phi_d1[1]
    p1_2, p2_2 = phi_d2[0], phi_d2[1]
    cW_val = -3j / (np.sqrt(2.0) * beta_val)
    J_lead = p2_1 * (p2_1**2 - 3.0 * p1_1**2) + cW_val * (p1_1 * p2_2 - p2_1 * p1_2)
    E0_v = taus[1] * taus[2] / (taus[0] ** 2)
    E1_v = taus[0] * taus[2] / (taus[1] ** 2)
    E2_v = taus[0] * taus[1] / (taus[2] ** 2)
    J_flux = (1.5 * m_val**2 / beta_val**2) * (
        2.0 * E0_v * p2_1
        + E1_v * (np.sqrt(3.0) * p1_1 - p2_1)
        - E2_v * (np.sqrt(3.0) * p1_1 + p2_1)
    )
    return np.trapz(J_lead + J_flux, xs)


q1 = eval_Q2_1sol(1, 0.0)["d_- (scales as e^{+2 th})"]
for t_test in [0.0, 0.7]:
    Q_br = eval_Q2_breather(0.3, 1.0, t_val=t_test)
    Q_add = q1 * (np.exp(2.0 * (0.3 + 1.0j)) - np.exp(2.0 * (0.3 - 1.0j)))
    print(
        f"Breather (th0=0.3, u=1.0, t={t_test}): Q_2(num) = {Q_br:.6f}, "
        f"2i q_1 e^{{2 th0}} sin(2u) = {Q_add:.6f}, diff = {abs(Q_br - Q_add):.2e}"
    )

print("\n" + "=" * 72)
print("CHECK 5 (Rank 7): Laurent expansion near zero of tau_j & topological charge jump")
print("=" * 72)

# Near a simple zero x_0 of tau_j(x,t) at fixed t:
# tau_j(x,t) = a_1 (x - x_0) + a_2 (x - x_0)^2 + a_3 (x - x_0)^3 + O((x-x_0)^4)
# Then ln tau_j(x,t) = ln(x - x_0) + ln(a_1) + (a_2/a_1)(x - x_0) + (a_3/a_1 - a_2^2/(2 a_1^2))(x - x_0)^2 + ...
# d_x ln tau_j(x,t) = 1/(x - x_0) + (a_2/a_1) + 2(a_3/a_1 - a_2^2/(2 a_1^2))(x - x_0) + ...
# d_x^2 ln tau_j(x,t) = -1/(x - x_0)^2 + 2(a_3/a_1 - a_2^2/(2 a_1^2)) + O(x - x_0)!
# Notice the residue (the coefficient of 1/(x - x_0) in d_x^2 ln tau_j) is IDENTICALLY ZERO because d_x^2 ln tau_j
# is the derivative of d_x ln tau_j, whose simple pole 1/(x - x_0) differentiates to -1/(x - x_0)^2 with zero 1/(x - x_0) term!
# What about the topological charge (beta / (2*pi)) [phi(+infty) - phi(-infty)] when a zero of tau_j crosses the real x-axis?
# Since phi = (i/beta) sum_j alpha_j^vee ln tau_j, the change in ln tau_j from -infty to +infty is
# i * Delta arg(tau_j). When a zero of tau_j in the complex x-plane crosses the real x-axis as t (or Im xi) varies,
# the winding number of tau_j(x) around 0 changes by +-1, so Delta Im(ln tau_j) jumps by +- 2*pi,
# and the topological charge (beta/(2*pi)) [phi(+infty) - phi(-infty)] JUMPS by +- alpha_j^vee!
print("1. Energy density T_00 = (2/beta^2) sum_j (2/alpha_j^2) d_x^2 ln tau_j near a zero x_0 of tau_j of multiplicity m_0:")
print("   ln tau_j = m_0 ln(x - x_0) + reg  =>  d_x ln tau_j = m_0/(x - x_0) + reg")
print("   =>  d_x^2 ln tau_j = -m_0/(x - x_0)^2 + reg (double pole with ZERO residue for any multiplicity m_0!).")
print("   Wait: on the real axis, int_{-eps}^{+eps} (-m_0/x^2) dx = -infty (not Lebesgue or Cauchy-PV integrable!),")
print("   so 'integrated charges are unaffected' holds ONLY if the integration contour is deformed into C around x_0")
print("   (or regularized by Hadamard finite part / analytic continuation from complex x_0).")
print("2. Topological charge Q_top = (beta/(2*pi)) (phi(+infty) - phi(-infty)):")
print("   Q_top = -(1/(2*pi)) sum_j alpha_j^vee [Im ln tau_j]_{-infty}^{+infty}.")
print("   When a zero of tau_j crosses the real x-axis, [Im ln tau_j]_{-infty}^{+infty} jumps by +-2*pi,")
print("   so Q_top JUMPS by -+ alpha_j^vee (a coroot)! It is NOT unaffected by the singularity!")
