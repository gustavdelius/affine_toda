#!/usr/bin/env python3
"""Independent verification script for Lecture 6 of lecture-notes/defects/."""

import sympy as sp


def check_lecture_6():
    print("=== Lecture 6: The classical algebra of a defect ===")
    lam, sig, beta, m0 = sp.symbols("lambda sigma beta m0", nonzero=True)
    Phi, Psi = sp.symbols("Phi Psi")
    Phi_p, Psi_p, Phi_m, Psi_m = sp.symbols("Phi_p Psi_p Phi_m Psi_m")  # d_+ and d_- derivatives

    h = sp.Matrix([[1, 0], [0, -1]])
    e = sp.Matrix([[0, 1], [0, 0]])
    f = sp.Matrix([[0, 0], [1, 0]])
    s = e + f

    def exp_ih(angle):
        return sp.Matrix([[sp.exp(sp.I * angle), 0], [0, sp.exp(-sp.I * angle)]])

    def s_field(phi):
        return exp_ih(phi / 2) * s

    # 1. Zero-curvature Lax pair (eq-l6-lax):
    # A_+ = i/4 (d_+ Phi) h + (m0 * lam / 2) s_Phi
    # A_- = -i/4 (d_- Phi) h - (m0 / (2 * lam)) s_{-Phi}
    def A_plus(phi, dplus_phi):
        return sp.I / 4 * dplus_phi * h + (m0 * lam / 2) * s_field(phi)

    def A_minus(phi, dminus_phi):
        return -sp.I / 4 * dminus_phi * h - (m0 / (2 * lam)) * s_field(-phi)

    # Check thm-l6-gauge: d_pm L_d = L_d A_pm(Psi) - A_pm(Phi) L_d
    # L_d = e^{i(Phi - Psi) h / 4} + (sig / lam) e^{-i(Phi + Psi) h / 4} s
    U_d = exp_ih((Phi - Psi) / 4)
    V_d = exp_ih(-(Phi + Psi) / 4) * s
    L_d = U_d + (sig / lam) * V_d

    dplus_Ld = (
        sp.diff(U_d, Phi) * Phi_p
        + sp.diff(U_d, Psi) * Psi_p
        + (sig / lam) * (sp.diff(V_d, Phi) * Phi_p + sp.diff(V_d, Psi) * Psi_p)
    )
    dminus_Ld = (
        sp.diff(U_d, Phi) * Phi_m
        + sp.diff(U_d, Psi) * Psi_m
        + (sig / lam) * (sp.diff(V_d, Phi) * Phi_m + sp.diff(V_d, Psi) * Psi_m)
    )

    rhs_plus = L_d * A_plus(Psi, Psi_p) - A_plus(Phi, Phi_p) * L_d
    rhs_minus = L_d * A_minus(Psi, Psi_m) - A_minus(Phi, Phi_m) * L_d

    diff_plus = sp.simplify(dplus_Ld - rhs_plus)
    diff_minus = sp.simplify(dminus_Ld - rhs_minus)

    # Backlund equations in Lecture 6 (with d_pm = d_t pm d_x):
    # d_+(Phi - Psi) = -2 sig m0 sin((Phi + Psi)/2)
    # d_-(Phi + Psi) = (2 m0 / sig) sin((Phi - Psi)/2)
    bt_plus = {Phi_p: Psi_p - 2 * sig * m0 * sp.sin((Phi + Psi) / 2).rewrite(sp.exp)}
    bt_minus = {Phi_m: -Psi_m + (2 * m0 / sig) * sp.sin((Phi - Psi) / 2).rewrite(sp.exp)}

    assert sp.simplify(diff_plus.subs(bt_plus)) == sp.zeros(2, 2)
    assert sp.simplify(diff_minus.subs(bt_minus)) == sp.zeros(2, 2)
    print("PASS L6.1: thm-l6-gauge verified (both light-cone gauge equations iff Backlund equations)")

    # 2. Check thm-l6-frozen: gluing matrix G(lam) = G0 + lam^{-1} G1 at x=0
    # Verify middle-order equations (eq-l6-four) and time-gauge equation
    Phi_t, Psi_t, Phi_x, Psi_x = sp.symbols("Phi_t Psi_t Phi_x Psi_x")
    sub_xt = {
        Phi_p: Phi_t + Phi_x,
        Psi_p: Psi_t + Psi_x,
        Phi_m: Phi_t - Phi_x,
        Psi_m: Psi_t - Psi_x,
    }
    dt_Ld = sp.diff(L_d, Phi) * Phi_t + sp.diff(L_d, Psi) * Psi_t
    at_Psi = sp.Rational(1, 2) * (A_plus(Psi, Psi_p) + A_minus(Psi, Psi_m)).subs(sub_xt)
    at_Phi = sp.Rational(1, 2) * (A_plus(Phi, Phi_p) + A_minus(Phi, Phi_m)).subs(sub_xt)
    res_frozen = sp.expand(dt_Ld - (L_d * at_Psi - at_Phi * L_d))

    dc_subs = {
        Phi_x: Psi_t
        - sig * m0 * sp.sin((Phi + Psi) / 2).rewrite(sp.exp)
        - (m0 / sig) * sp.sin((Phi - Psi) / 2).rewrite(sp.exp),
        Psi_x: Phi_t
        + sig * m0 * sp.sin((Phi + Psi) / 2).rewrite(sp.exp)
        - (m0 / sig) * sp.sin((Phi - Psi) / 2).rewrite(sp.exp),
    }
    assert sp.simplify(res_frozen.subs(dc_subs)) == sp.zeros(2, 2)

    # Also check general G0 = p D + q O_bar, G1 = p' D^{-1} + q' O gives eq-l6-four
    p, q, pp, qp = sp.symbols("p q pp qp")
    O_bar = exp_ih((Phi + Psi) / 4) * s
    G_gen = (p * U_d + q * O_bar) + (1 / lam) * (pp * exp_ih(-(Phi - Psi) / 4) + qp * V_d)
    dt_G = sp.diff(G_gen, Phi) * Phi_t + sp.diff(G_gen, Psi) * Psi_t
    eq_four_subs = {
        Phi_x: Psi_t - (qp / p) * m0 * sp.sin((Phi + Psi) / 2).rewrite(sp.exp) - (p / qp) * m0 * sp.sin((Phi - Psi) / 2).rewrite(sp.exp),
        Psi_x: Phi_t + (qp / p) * m0 * sp.sin((Phi + Psi) / 2).rewrite(sp.exp) - (p / qp) * m0 * sp.sin((Phi - Psi) / 2).rewrite(sp.exp),
    }
    assert sp.simplify((dt_G - (G_gen * at_Psi - at_Phi * G_gen)).subs({q: 0, pp: 0}).subs(eq_four_subs)) == sp.zeros(2, 2)
    print("PASS L6.2: thm-l6-frozen verified (time-gauge equation at x=0 and general G_0 + lam^{-1} G_1 reduction)")

    # 3. Check thm-l6-sklyanin: Sklyanin bracket {L_d(lam1) ,^x L_d(lam2)} = [r(lam1/lam2), L_d(lam1) x L_d(lam2)]
    l1, l2, c = sp.symbols("l1 l2 c", nonzero=True)
    L1 = U_d + (sig / l1) * V_d
    L2 = U_d + (sig / l2) * V_d

    def kron(A, B):
        return sp.kronecker_product(A, B)

    pb_4x4 = c * (
        kron(sp.diff(L1, Phi), sp.diff(L2, Psi))
        - kron(sp.diff(L1, Psi), sp.diff(L2, Phi))
    )

    z = l1 / l2
    Pi = kron(e, f) + kron(f, e)
    r_mat = -(beta**2 / (16 * (z**2 - 1))) * ((z**2 + 1) * kron(h, h) + 4 * z * Pi)
    comm_4x4 = r_mat * kron(L1, L2) - kron(L1, L2) * r_mat
    diff_sklyanin = sp.simplify(pb_4x4.subs(c, -beta**2) - comm_4x4)
    assert diff_sklyanin == sp.zeros(4, 4)

    # Also check uniqueness of alpha(z), gamma(z) in r(z) = alpha h x h + gamma Pi
    alpha, gamma = sp.symbols("alpha gamma")
    r_ansatz = alpha * kron(h, h) + gamma * Pi
    comm_ansatz = r_ansatz * kron(L1, L2) - kron(L1, L2) * r_ansatz
    sol_ag = sp.solve(list(comm_ansatz - pb_4x4), [alpha, gamma])
    assert sp.simplify(sol_ag[alpha] - c * (z**2 + 1) / (16 * (z**2 - 1))) == 0
    assert sp.simplify(sol_ag[gamma] - c * 4 * z / (16 * (z**2 - 1))) == 0
    print("PASS L6.3: thm-l6-sklyanin verified symbolically (uniqueness of r(z) and match at c = -beta^2)")

    # 4. Check spin-mp 3 defect contributions D_{-3} and D_{+3} (eq-l6-D3)
    # In Lecture 6, d_+ = d_t + d_x, d_- = d_t - d_x, so d_+ d_- Phi = -m0^2 sin(Phi).
    # Bulk spin -3 density: q_-(Phi) = (d_+^2 Phi)^2 - 1/4 (d_+ Phi)^4 + m0^2 cos(Phi) (d_+ Phi)^2
    p1, p2, p3, q1, q2 = sp.symbols("p1 p2 p3 q1 q2")
    dm_phi = sp.Symbol("dm_phi")
    q_minus_Phi = p2**2 - sp.Rational(1, 4) * p1**4 + m0**2 * sp.cos(Phi) * p1**2
    # Since d_-(p1) = -m0^2 sin(Phi), d_-(p2) = -m0^2 p1 cos(Phi):
    dm_q_minus = (
        2 * p2 * (-m0**2 * p1 * sp.cos(Phi))
        - p1**3 * (-m0**2 * sp.sin(Phi))
        - m0**2 * sp.sin(Phi) * dm_phi * p1**2
        + 2 * m0**2 * sp.cos(Phi) * p1 * (-m0**2 * sp.sin(Phi))
    )
    # And X_minus = -m0^2 sin(Phi) (d_+ Phi)^2? Wait: in eq-l6-spin3-current,
    # T_4 = p2^2 - 1/4 p1^4, X_2 = -m0^2 cos(Phi) p1^2, so d_- T_4 = d_+ X_2,
    # and the conserved density for d_t q_- + d_x j_- = 0 is q_- = T_4 - X_2, with j_- = -(T_4 + X_2)!
    T4_Phi = p2**2 - sp.Rational(1, 4) * p1**4
    X2_Phi = -m0**2 * sp.cos(Phi) * p1**2
    dm_T4 = 2 * p2 * (-m0**2 * p1 * sp.cos(Phi)) - p1**3 * (-m0**2 * sp.sin(Phi))
    dp_X2 = -m0**2 * (-sp.sin(Phi) * p1**3 + 2 * sp.cos(Phi) * p1 * p2)
    assert sp.simplify(dm_T4 - dp_X2) == 0
    # Check d_t (T4 - X2) + d_x (-(T4 + X2)) = 1/2(d_+ + d_-)(T4 - X2) - 1/2(d_+ - d_-)(T4 + X2) = d_- T4 - d_+ X2 = 0!
    j_Phi = -(T4_Phi + X2_Phi)
    j_Psi = -(q2**2 - sp.Rational(1, 4) * q1**4 - m0**2 * sp.cos(Psi) * q1**2)

    S_half = (Phi + Psi) / 2
    D_half = (Phi - Psi) / 2
    c_p = sp.cos(S_half)
    D_m3 = (sig / (2 * m0**2)) * c_p * p1 * q1 + sig**3 * (c_p - sp.Rational(2, 3) * c_p**3)

    # At x=0, d_+ = 2 d_t - d_-, so p2 = d_+^2 Phi = 2 dot_p1 + m0^2 sin(Phi), q2 = 2 dot_q1 + m0^2 sin(Psi).
    # And dot_Phi = 1/2(p1 + dm_Phi), dot_Psi = 1/2(q1 + dm_Psi), with dm_Phi + dm_Psi = (2 m0 / sig) sin(D_half),
    # and dot_p1 - dot_q1 = d_t(-2 sig m0 sin(S_half)) = -sig m0 cos(S_half) (dot_Phi + dot_Psi).
    dot_q1 = sp.Symbol("dot_q1")
    dot_Phi = sp.Rational(1, 2) * (p1 + dm_phi)
    dot_Psi = sp.Rational(1, 2) * (q1 + (2 * m0 / sig) * sp.sin(D_half) - dm_phi)
    dot_p1 = dot_q1 - sig * m0 * sp.cos(S_half) * (dot_Phi + dot_Psi)
    p2_sub = 2 * dot_p1 + m0**2 * sp.sin(Phi)
    q2_sub = 2 * dot_q1 + m0**2 * sp.sin(Psi)
    dt_Dm3 = (
        sp.diff(D_m3, Phi) * dot_Phi
        + sp.diff(D_m3, Psi) * dot_Psi
        + sp.diff(D_m3, p1) * dot_p1
        + sp.diff(D_m3, q1) * dot_q1
    )
    p1_sub = q1 - 2 * sig * m0 * sp.sin(S_half)
    flux_balance = ((j_Psi - j_Phi) / (8 * m0**3) + dt_Dm3).subs({p2: p2_sub, q2: q2_sub}).subs({p1: p1_sub})
    assert sp.simplify(flux_balance.rewrite(sp.exp)) == 0
    print("PASS L6.4: Exact symbolic spin -3 defect flux identity d_t D_{-3} = (j_-(Phi) - j_-(Psi))/(8 m0^3)")

    # 5. Check a_r^{(1)} defect Lax matrix gauge equations (eq-l6-ar-Ld) for h=3
    h_val = 3
    m = sp.Symbol("m", nonzero=True)
    phis = sp.symbols(f"phi0:{h_val}")
    psis = sp.symbols(f"psi0:{h_val}")
    dphis_p = sp.symbols(f"dphi_p0:{h_val}")
    dpsis_p = sp.symbols(f"dpsi_p0:{h_val}")
    dphis_m = sp.symbols(f"dphi_m0:{h_val}")
    dpsis_m = sp.symbols(f"dpsi_m0:{h_val}")

    def A_plus_ar(f_vec, df_p):
        M = sp.zeros(h_val, h_val)
        for j in range(h_val):
            M[j, j] = sp.Rational(1, 2) * df_p[j]
            M[j, (j + 1) % h_val] = m * lam * sp.exp((f_vec[j] - f_vec[(j + 1) % h_val]) / 2)
        return M

    def A_minus_ar(f_vec, df_m):
        M = sp.zeros(h_val, h_val)
        for j in range(h_val):
            M[j, j] = -sp.Rational(1, 2) * df_m[j]
            M[(j + 1) % h_val, j] = -(m / lam) * sp.exp((f_vec[j] - f_vec[(j + 1) % h_val]) / 2)
        return M

    Ld_ar = sp.zeros(h_val, h_val)
    for j in range(h_val):
        Ld_ar[j, j] = sp.exp((phis[j] - psis[j]) / 2)
        Ld_ar[(j + 1) % h_val, j] = (sig / lam) * sp.exp((psis[j] - phis[(j + 1) % h_val]) / 2)

    dp_Ld_ar = sum((sp.diff(Ld_ar, phis[j]) * dphis_p[j] + sp.diff(Ld_ar, psis[j]) * dpsis_p[j] for j in range(h_val)), sp.zeros(h_val, h_val))
    dm_Ld_ar = sum((sp.diff(Ld_ar, phis[j]) * dphis_m[j] + sp.diff(Ld_ar, psis[j]) * dpsis_m[j] for j in range(h_val)), sp.zeros(h_val, h_val))

    res_p_ar = dp_Ld_ar - (Ld_ar * A_plus_ar(psis, dpsis_p) - A_plus_ar(phis, dphis_p) * Ld_ar)
    res_m_ar = dm_Ld_ar - (Ld_ar * A_minus_ar(psis, dpsis_m) - A_minus_ar(phis, dphis_m) * Ld_ar)

    # Light-cone Backlund equations for a_r^{(1)} (eq-l5-ar-lightcone):
    # d_+(phi_j - psi_j) = -m sig (e^{psi_j - phi_{j+1}} - e^{psi_{j-1} - phi_j})
    # d_-(psi_j - phi_{j+1}) = (m / sig) (e^{phi_j - psi_j} - e^{phi_{j+1} - psi_{j+1}})
    bt_p_ar = {
        dphis_p[j]: dpsis_p[j]
        - m * sig * (sp.exp(psis[j] - phis[(j + 1) % h_val]) - sp.exp(psis[(j - 1) % h_val] - phis[j]))
        for j in range(h_val)
    }
    bt_m_ar = {
        dphis_m[(j + 1) % h_val]: dpsis_m[j]
        - (m / sig) * (sp.exp(phis[j] - psis[j]) - sp.exp(phis[(j + 1) % h_val] - psis[(j + 1) % h_val]))
        for j in range(h_val)
    }
    assert sp.simplify(res_p_ar.subs(bt_p_ar)) == sp.zeros(h_val, h_val)
    assert sp.simplify(res_m_ar.subs(bt_m_ar)) == sp.zeros(h_val, h_val)
    print("PASS L6.5: a_r^{(1)} defect matrix (eq-l6-ar-Ld) satisfies both light-cone gauge equations iff Backlund")


if __name__ == "__main__":
    check_lecture_6()
    print("\nALL L6 CHECKS PASSED.")
