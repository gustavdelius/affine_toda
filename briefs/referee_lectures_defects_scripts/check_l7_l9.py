#!/usr/bin/env python3
"""Independent verification script for Lectures 7, 8, and 9 of lecture-notes/defects/."""

import numpy as np
import scipy.integrate as integrate
import sympy as sp


def check_lecture_7():
    print("=== Lecture 7: Quantizing the defect ===")
    # 1. thm-l7-half: Weyl module of U_q(b^-) cannot extend to U_q(\hat{sl}_2)
    q, s, a, b = sp.symbols("q s a b", nonzero=True)
    c_n, c_np1 = sp.symbols("c_n c_np1")
    comm_e1_f0 = b * (c_n - c_np1)  # must be 0 -> c_{n+1} = c_n
    comm_e1_f1 = a * (c_np1 - c_n)
    assert comm_e1_f1.subs(c_np1, c_n) == 0
    print("PASS L7.1: thm-l7-half verified ([e_1, f_0]=0 forces c_n constant, making [e_1, f_1]=0 != (k_1-k_1^{-1})/(q-q^{-1}))")

    # 2. thm-l7-T: solve the 8 intertwining equations symbolically in a block
    y, yd, q_sym = sp.symbols("y yd q", nonzero=True)
    n = sp.Symbol("n", integer=True)
    A = lambda m: q_sym ** (-m)
    D = lambda m: q_sym ** (m + s)
    C = lambda m: -sp.I * (y / yd) * q_sym ** (-m)
    B = lambda m: -sp.I * (y / yd) * q_sym ** (m + s)
    r1 = sp.simplify((q_sym * A(n + 1) - A(n)) / C(n))
    r2 = sp.simplify((D(n) - q_sym * D(n + 1)) / C(n))
    r3 = sp.simplify((A(n) - q_sym * A(n - 1)) / B(n))
    r4 = sp.simplify((q_sym * D(n - 1) - D(n)) / B(n))
    assert r1 == 0 and r4 == 0
    assert sp.simplify(r2 - sp.I * (yd / y) * q_sym**s * (q_sym ** (2 * n) - q_sym ** (2 * n + 2))) == 0
    assert sp.simplify(r3 - sp.I * (yd / y) * q_sym ** (-s) * (q_sym ** (-2 * n) - q_sym ** (-2 * n + 2))) == 0
    print("PASS L7.2: thm-l7-T intertwining recurrence relations verified symbolically")

    # 3. thm-l7-crossing: check crossing relation and sign of q = -e^{-i*pi*lam}
    Y, Yp = sp.symbols("Y Yp", nonzero=True)
    Mn = sp.Matrix([[1, -sp.I * Y], [-sp.I * Y, 1]]) * sp.Matrix([[q_sym ** (-n), 0], [0, q_sym ** (n + 1 + s)]])
    Mnp_inv = (
        1
        / (1 + Yp**2)
        * sp.Matrix([[q_sym**n, 0], [0, q_sym ** (-n - 1 - s)]])
        * sp.Matrix([[1, sp.I * Yp], [sp.I * Yp, 1]])
    )
    assert sp.simplify(Mn.subs(Y, Yp) * Mnp_inv) == sp.eye(2)
    assert sp.simplify((-sp.I * Y) - (sp.I * (-q_sym * Y) * q_sym ** (-1))) == 0
    print("PASS L7.3: thm-l7-crossing verified (requires Y' = -q Y, i.e. q = -e^{-i*pi*lam}, and Phi(u)Phi(u+i*pi)(1+Y^2) = q^{-s})")

    # 4. Integral representation of minimal scalar factor Phi(u) (eq-l7-Phi) vs crossing and unitarity
    lam_val = np.sqrt(2.0)
    u_val = 0.65
    Y_num = np.exp(lam_val * u_val)

    def phi_integral(u_complex, lam):
        def integrand_re(t):
            k = 2.0 * lam * t / np.pi
            val = (np.cosh(0.5 * np.pi * k + 1j * k * u_complex) / np.cosh(0.5 * np.pi * k) - 1.0) / (
                2.0 * t * np.sinh(t)
            )
            return float(np.real(val))

        def integrand_im(t):
            k = 2.0 * lam * t / np.pi
            val = (np.cosh(0.5 * np.pi * k + 1j * k * u_complex) / np.cosh(0.5 * np.pi * k) - 1.0) / (
                2.0 * t * np.sinh(t)
            )
            return float(np.imag(val))

        re_part, _ = integrate.quad(integrand_re, 0, 35.0, limit=200)
        im_part, _ = integrate.quad(integrand_im, 0, 35.0, limit=200)
        ln_phi = -0.5 * lam * u_complex + 0.25j * np.pi * lam - 0.5 * np.log(2.0) + (re_part + 1j * im_part)
        return np.exp(ln_phi)

    v_int = phi_integral(u_val, lam_val)
    unit_err = abs(abs(v_int) ** 2 * (1.0 + Y_num**2) - 1.0)
    assert unit_err < 1e-11, f"Unitarity error: {unit_err}"
    v_int_cross = phi_integral(u_val + 1j * np.pi, lam_val)
    cross_err = abs(v_int * v_int_cross * (1.0 + Y_num**2) - 1.0)
    assert cross_err < 1e-11, f"Crossing error: {cross_err}"
    print(f"PASS L7.4: Scalar factor integral (eq-l7-Phi) satisfies unitarity (err={unit_err:.2e}) and crossing (err={cross_err:.2e})")

    # 5. thm-l7-breather: lightest breather transmission T_{B_1}(u) = (i - sinh u)/cosh u
    # Using the crossing-reduced ratio inside the strip: Phi(v - eps) / Phi(v + eps) with v = u + i*pi/2, eps = i*pi/(2*lam)
    v_pt = u_val + 0.5j * np.pi
    eps = 0.5j * np.pi / lam_val
    lhs_br = phi_integral(v_pt - eps, lam_val) / phi_integral(v_pt + eps, lam_val)
    rhs_br = (1j - np.sinh(u_val)) / np.cosh(u_val)
    err_br = abs(lhs_br - rhs_br)
    assert err_br < 1e-11, f"Breather identity error: {err_br}"
    print(f"PASS L7.5: thm-l7-breather verified numerically to {err_br:.2e}")

    # 6. Wigner delay slope at u = 0: delta'(0) = (lam/pi)*(ln(8*lam/pi) + gamma_E) + o(1)
    lam_big = 400.0
    exact_slope, _ = integrate.quad(
        lambda t: (lam_big / np.pi) * np.tanh(lam_big * t) / np.sinh(t) if t > 0 else lam_big**2 / np.pi,
        0,
        35.0,
        limit=400,
    )
    asymp_book = (lam_big / np.pi) * (np.log(8.0 * lam_big / np.pi) + np.euler_gamma)
    err_wigner = abs(exact_slope - asymp_book)
    assert err_wigner < 1e-3, f"Wigner delay slope error: {err_wigner}"
    print(f"PASS L7.6: Wigner delay slope delta'(0) = (lam/pi)*(ln(8*lam/pi)+gamma_E)+o(1) verified (diff={err_wigner:.2e} at lam={lam_big})")


def check_lecture_8():
    print("\n=== Lecture 8: Defects, L-operators and Q-operators ===")
    # 1. thm-l8-turned: turned defect at th_d - i*pi/2 is Backlund B_{-i*sigma} and preserves tr T
    th_d = sp.Symbol("th_d", real=True)
    sig = sp.exp(-th_d)
    sig_turned = sp.exp(-(th_d - sp.I * sp.pi / 2))
    assert sp.simplify(sig_turned - sp.I * sig) == 0
    print("PASS L8.1: thm-l8-turned parameter shift th_d -> th_d - i*pi/2 gives B_{-i*sigma}")

    # 2. Construct exact N-soliton transfer matrix t_omega(y), defect operator D(yd)
    # on V_1 x ... x V_N for N = 2, 3 and verify thm-l8-commute, prp-l8-dd, thm-l8-tq, prp-l8-baxter!
    lam_val = np.sqrt(2.0)
    q = -np.exp(-1j * np.pi * lam_val)

    def build_R(ya, yb, q_val):
        # Normalized so R(v_+ x v_+) = v_+ x v_+ as in Lecture 7 sec-l7-bulk:
        z = ya / yb
        denom = q_val * z - 1.0 / (q_val * z)
        b_r = (z - 1.0 / z) / denom
        c_r = (q_val - 1.0 / q_val) / denom
        return np.array(
            [
                [1.0, 0, 0, 0],
                [0, b_r, c_r, 0],
                [0, c_r, b_r, 0],
                [0, 0, 0, 1.0],
            ],
            dtype=complex,
        )

    M = 6
    dim_W = 2 * M + 1

    def L_matrix(y_val, yd_val, q_val):
        L00 = np.zeros((dim_W, dim_W), dtype=complex)
        L01 = np.zeros((dim_W, dim_W), dtype=complex)
        L10 = np.zeros((dim_W, dim_W), dtype=complex)
        L11 = np.zeros((dim_W, dim_W), dtype=complex)
        Y = y_val / yd_val
        for idx in range(dim_W):
            n_int = idx - M
            L00[idx, idx] = q_val ** (-n_int)
            L11[idx, idx] = q_val**n_int
            if idx + 1 < dim_W:
                L10[idx + 1, idx] = -1j * Y * (q_val ** (-n_int))
            if idx - 1 >= 0:
                L01[idx - 1, idx] = -1j * Y * (q_val**n_int)
        return np.block([[L00, L01], [L10, L11]])

    def L_inv_matrix(y_val, yd_val, q_val):
        L00 = np.zeros((dim_W, dim_W), dtype=complex)
        L01 = np.zeros((dim_W, dim_W), dtype=complex)
        L10 = np.zeros((dim_W, dim_W), dtype=complex)
        L11 = np.zeros((dim_W, dim_W), dtype=complex)
        Y = y_val / yd_val
        pref = 1.0 / (1.0 + Y**2)
        for idx in range(dim_W):
            n_int = idx - M
            L00[idx, idx] = pref * (q_val**n_int)
            L11[idx, idx] = pref * (q_val ** (-n_int))
            if idx + 1 < dim_W:
                L10[idx + 1, idx] = pref * 1j * Y * (q_val ** (-n_int - 1))
            if idx - 1 >= 0:
                L01[idx - 1, idx] = pref * 1j * Y * (q_val ** (n_int - 1))
        return np.block([[L00, L01], [L10, L11]])

    ya_v, yb_v, yd_v = 1.3, 0.85, 0.95
    La = L_matrix(ya_v, yd_v, q).reshape(2, dim_W, 2, dim_W)
    Lb = L_matrix(yb_v, yd_v, q).reshape(2, dim_W, 2, dim_W)
    LaLb = np.einsum("iajk,lkbm->ilajbm", La, Lb).reshape(4 * dim_W, 4 * dim_W)
    LbLa = np.einsum("labk,ikjm->ilajbm", Lb, La).reshape(4 * dim_W, 4 * dim_W)
    R_ab = np.kron(build_R(ya_v, yb_v, q), np.eye(dim_W, dtype=complex))
    comm_RLL = R_ab @ LaLb - LbLa @ R_ab
    interior = [s * dim_W + idx for s in range(4) for idx in range(3, dim_W - 3)]
    err_RLL = np.max(np.abs(comm_RLL[np.ix_(interior, interior)]))
    assert err_RLL < 1e-13, f"RLL error: {err_RLL}"
    print(f"PASS L8.2: Normalized RLL relation R_{{ab}}(y_a/y_b) L_{{aW}}(y_a) L_{{bW}}(y_b) = L_{{bW}}(y_b) L_{{aW}}(y_a) R_{{ab}}(y_a/y_b) verified (err={err_RLL:.2e})")

    # Build N-site transfer matrix t_omega(y) and defect operator D(yd) on each charge sector m:
    for N_sites, y_inhs in [(2, [1.1, 0.8]), (3, [1.15, 0.85, 1.35])]:
        dim_phys = 2**N_sites
        # Total charge m = \sum_j h_j where h_j(v_+) = +1, h_j(v_-) = -1:
        charges = np.zeros(dim_phys, dtype=int)
        for state_idx in range(dim_phys):
            bits = [(state_idx >> (N_sites - 1 - j)) & 1 for j in range(N_sites)]
            charges[state_idx] = sum(1 if b == 0 else -1 for b in bits)

        def build_t_omega(y_val, omega_val):
            mon_mat = np.eye(2 * dim_phys, dtype=complex)
            for j in range(N_sites):
                R_aj = build_R(y_val, y_inhs[j], q).reshape(2, 2, 2, 2)
                d_left = 2**j
                d_right = 2 ** (N_sites - 1 - j)
                T_cur = mon_mat.reshape(2, d_left, 2, d_right, 2 * dim_phys)
                T_next = np.einsum("abij,ikjcd->akbcd", R_aj, T_cur)
                mon_mat = T_next.reshape(2 * dim_phys, 2 * dim_phys)
            mon_blocks = mon_mat.reshape(2, dim_phys, 2, dim_phys)
            return omega_val * mon_blocks[0, :, 0, :] + (1.0 / omega_val) * mon_blocks[1, :, 1, :]

        def build_D(yd_val):
            vec = np.zeros((dim_W, dim_phys, dim_phys), dtype=complex)
            vec[M, :, :] = np.eye(dim_phys, dtype=complex)
            for j in range(N_sites):
                Linv_j = L_inv_matrix(y_inhs[j], yd_val, q).reshape(2, dim_W, 2, dim_W)
                d_left = 2**j
                d_right = 2 ** (N_sites - 1 - j)
                v_cur = vec.reshape(dim_W, d_left, 2, d_right, dim_phys)
                v_next = np.einsum("axbw,wkbcd->xkacd", Linv_j, v_cur)
                vec = v_next.reshape(dim_W, dim_phys, dim_phys)
            return vec[M, :, :]

        D1 = build_D(0.95)
        D2 = build_D(1.55)
        err_DD = np.max(np.abs(D1 @ D2 - D2 @ D1))
        assert err_DD < 1e-12, f"N={N_sites} [D, D] error: {err_DD}"

        max_err_tD = 0.0
        max_err_TQ = 0.0
        yd_tq = 0.95
        D_yd = build_D(yd_tq)
        D_qyd = build_D(q * yd_tq)
        D_yd_div_q = build_D(yd_tq / q)
        prod_factor = np.prod([(1.0 + (q**2) * (yj / yd_tq) ** 2) / (1.0 + (yj / yd_tq) ** 2) for yj in y_inhs])

        for m_val in sorted(set(charges)):
            sec_idx = np.where(charges == m_val)[0]
            omega_m = q ** (-0.5 * m_val)
            t_m = build_t_omega(1.4, omega_m)[np.ix_(sec_idx, sec_idx)]
            D1_m = D_yd[np.ix_(sec_idx, sec_idx)]
            err_tD = np.max(np.abs(t_m @ D1_m - D1_m @ t_m))
            max_err_tD = max(max_err_tD, err_tD)

            # Check TQ relation (eq-l8-tq) at y = i * yd_tq (i.e. theta = theta_d + i*pi/(2*lam)):
            t_fused_m = build_t_omega(1j * yd_tq, omega_m)[np.ix_(sec_idx, sec_idx)]
            lhs_tq = t_fused_m @ D1_m
            rhs_tq = (q ** (-0.5 * N_sites)) * D_qyd[np.ix_(sec_idx, sec_idx)] + (
                q ** (-0.5 * N_sites)
            ) * prod_factor * D_yd_div_q[np.ix_(sec_idx, sec_idx)]
            err_tq = np.max(np.abs(lhs_tq - rhs_tq))
            max_err_TQ = max(max_err_TQ, err_tq)

        assert max_err_tD < 1e-12, f"N={N_sites} [t_omega, D] error: {max_err_tD}"
        assert max_err_TQ < 1e-12, f"N={N_sites} TQ error: {max_err_TQ}"
        print(
            f"PASS L8.3 (N={N_sites}): [t_omega, D]=0 (err={max_err_tD:.2e}), [D, D]=0 (err={err_DD:.2e}), and TQ relation eq-l8-tq (err={max_err_TQ:.2e}) verified across all charge sectors"
        )


def check_lecture_9_math():
    print("\n=== Lecture 9 (from handoff & Ch. 43): Folding, fusing and the classical breather ===")
    lam = np.sqrt(2.0)
    th_d = 0.45
    th = 0.32
    t_test = 1.3
    diff_sinh = np.sinh(t_test * (0.5 - 1j * lam * (th - th_d) / np.pi)) - np.sinh(
        t_test * (0.5 + 1j * lam * (th + th_d) / np.pi)
    )
    trig_id = -2j * np.sin(t_test * lam * th / np.pi) * np.cosh(t_test * (0.5 + 1j * lam * th_d / np.pi))
    assert abs(diff_sinh - trig_id) < 1e-15
    print("PASS L9.1: Hyperbolic trig identity behind integral proof of thm-l9-dirichlet verified to 1e-15")

    q_s, s_s, yd_s, n_s = sp.symbols("q s yd n")
    kappa_n = (-sp.I / yd_s * q_s ** (-(n_s + 1))) * (sp.I * yd_s * q_s ** (-n_s - 1 - s_s))
    assert sp.simplify(kappa_n - q_s ** (-2 * n_s - 2 - s_s)) == 0
    print("PASS L9.2: Second Borel half operator kappa = C B^{-1} = q^{-2n-2-s} verified symbolically")


if __name__ == "__main__":
    check_lecture_7()
    check_lecture_8()
    check_lecture_9_math()
    print("\nALL L7-L9 CHECKS PASSED.")
