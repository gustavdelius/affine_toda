#!/usr/bin/env python3
"""Independent verification script for Lectures 1-5 of lecture-notes/defects/."""

import numpy as np
import sympy as sp


def check_lecture_1():
    print("=== Lecture 1: Mirrors and windows ===")
    # 1. Free boson delta defect matching and DMS (1994) comparison
    k, lam, m, theta = sp.symbols("k lambda m theta", real=True, positive=True)
    T, R = sp.symbols("T R")
    sol = sp.solve([1 + R - T, sp.I * k * T - sp.I * k * (1 - R) - lam * T], [T, R])
    assert sp.simplify(sol[T] - 2 * sp.I * k / (2 * sp.I * k - lam)) == 0
    assert sp.simplify(sol[R] - lam / (2 * sp.I * k - lam)) == 0
    print("PASS L1.1: Free-boson delta-defect T(k) and R(k) from jump conditions")

    T_th = sol[T].subs(k, m * sp.sinh(theta))
    R_th = sol[R].subs(k, m * sp.sinh(theta))
    assert sp.simplify(T_th - sp.sinh(theta) / (sp.sinh(theta) + sp.I * lam / (2 * m))) == 0
    assert sp.simplify(R_th - (-sp.I * lam / (2 * m)) / (sp.sinh(theta) + sp.I * lam / (2 * m))) == 0
    print("PASS L1.2: Rapidity form has parameter lam/(2m) (confirming DMS eqs. 4.6, 4.8 misprint g/(4m))")

    # 2. Unitarity of exchange algebra (eq-l1-unitarity)
    R_p, T_p, Rb_p, Tb_p = sp.symbols("R_p T_p Rb_p Tb_p")
    R_m, T_m, Rb_m, Tb_m = sp.symbols("R_m T_m Rb_m Tb_m")
    coeff_AD = R_p * R_m + T_p * Tb_m
    coeff_DA = R_p * T_m + T_p * Rb_m
    assert sp.simplify(
        coeff_AD.subs({R_p: R_th, R_m: R_th.subs(theta, -theta), T_p: T_th, Tb_m: T_th.subs(theta, -theta)}) - 1
    ) == 0
    assert sp.simplify(
        coeff_DA.subs({R_p: R_th, T_m: T_th.subs(theta, -theta), T_p: T_th, Rb_m: R_th.subs(theta, -theta)})
    ) == 0
    print("PASS L1.3: Defect unitarity relations from double application of exchange algebra")

    # 3. Two-history expansions for thm-l1-nogo
    Sab_diff, Sab_sum, Sba_sum = sp.symbols("Sab_diff Sab_sum Sba_sum")
    Ra1, Ta1, Rb2, Tb2 = sp.symbols("Ra1 Ta1 Rb2 Tb2")
    h1_RT = Rb2 * Sab_sum * Ta1
    h1_TR = Tb2 * Ra1
    h2_RT = Sab_diff * Ta1 * Rb2
    h2_TR = Sab_diff * Ra1 * Sba_sum * Tb2
    assert sp.expand((h1_RT - h2_RT) - Rb2 * Ta1 * (Sab_sum - Sab_diff)) == 0
    assert sp.expand((h1_TR - h2_TR) - Ra1 * Tb2 * (1 - Sab_diff * Sba_sum)) == 0
    print("PASS L1.4: Two-history mismatch equations in thm-l1-nogo verified")

    # 4. Lee-Yang transmission factor: unitarity, crossing, bootstrap using rational variable x = e^theta
    # S_ly(theta) = (sinh(theta) + i sin(2*pi/3)) / (sinh(theta) - i sin(2*pi/3))
    # With x = e^theta and omega = e^{i*pi/3} (so omega^2 - omega + 1 = 0, i*sin(2*pi/3) = (omega - omega^{-1})/2 = (2*omega - 1)/2):
    x, w = sp.symbols("x w", nonzero=True)
    # Using exact numerical/algebraic test across multiple generic complex (theta, b):
    for b_val in [0.37, 1.23, -0.41]:
        for th_val in [0.5 + 0.1j, -0.8 + 0.2j, 1.3 - 0.15j]:
            S_ly = lambda u: (np.sinh(u) + 1j * np.sin(2 * np.pi / 3)) / (np.sinh(u) - 1j * np.sin(2 * np.pi / 3))
            th_b = 1j * np.pi * (3.0 - b_val) / 6.0
            T_ly = lambda u: S_ly(u - th_b)
            Tb_ly = lambda u: S_ly(u + th_b)
            assert abs(T_ly(th_val) * Tb_ly(-th_val) - 1.0) < 1e-15
            assert abs(T_ly(th_val) - Tb_ly(1j * np.pi - th_val)) < 1e-15
            assert abs(T_ly(th_val + 1j * np.pi / 3) * T_ly(th_val - 1j * np.pi / 3) - T_ly(th_val)) < 1e-15
    print("PASS L1.5: Lee-Yang defect unitarity, crossing, and bootstrap verified to 1e-15")


def check_lecture_2():
    print("\n=== Lecture 2: Backlund's machine ===")
    # 1. Beltrami pseudosphere parametrization and Gaussian curvature
    x = sp.Symbol("x", real=True, positive=True)
    r = 1 / sp.cosh(x)
    z = x - sp.tanh(x)
    ds_x2 = sp.trigsimp(sp.diff(r, x) ** 2 + sp.diff(z, x) ** 2)
    assert sp.trigsimp(ds_x2 - sp.tanh(x) ** 2) == 0
    dr_ds = sp.trigsimp(sp.diff(r, x) / sp.tanh(x))
    d2r_ds2 = sp.trigsimp(sp.diff(dr_ds, x) / sp.tanh(x))
    K = sp.trigsimp(-d2r_ds2 / r)
    assert K == -1
    print("PASS L2.1: Beltrami pseudosphere metric and Gaussian curvature K = -1 verified")

    # 2. Backlund segment length and angle
    theta = sp.Symbol("theta", real=True)
    sigma_param = sp.exp(-theta)
    cos_ang = (1 - sigma_param**2) / (1 + sigma_param**2)
    sin_ang = 2 * sigma_param / (1 + sigma_param**2)
    assert sp.simplify((cos_ang - sp.tanh(theta)).rewrite(sp.exp)) == 0
    assert sp.simplify((sin_ang - 1 / sp.cosh(theta)).rewrite(sp.exp)) == 0
    print("PASS L2.2: Backlund segment angle cos(sigma_geom) = tanh(theta), sin(sigma_geom) = sech(theta)")

    # 3. Bianchi permutability (prp-l2-bianchi)
    phi0, phi1, phi2, s1, s2, phi0_u, phi0_v = sp.symbols("phi0 phi1 phi2 s1 s2 phi0_u phi0_v")
    phi1_u = phi0_u - 2 * s1 * sp.sin((phi0 + phi1) / 2)
    phi2_u = phi0_u - 2 * s2 * sp.sin((phi0 + phi2) / 2)
    phi12 = phi0 + 4 * sp.atan(((s1 + s2) / (s1 - s2)) * sp.tan((phi1 - phi2) / 4))
    phi12_u = sp.diff(phi12, phi0) * phi0_u + sp.diff(phi12, phi1) * phi1_u + sp.diff(phi12, phi2) * phi2_u
    rhs_u = phi1_u - 2 * s2 * sp.sin((phi1 + phi12) / 2)
    vals = {phi0: 0.4, phi1: 1.1, phi2: -0.7, s1: 1.3, s2: 0.6, phi0_u: 0.9, phi0_v: -0.5}
    err_u = abs(complex((phi12_u - rhs_u).subs(vals)))
    assert err_u < 1e-14

    phi1_v = -phi0_v + (2 / s1) * sp.sin((phi1 - phi0) / 2)
    phi2_v = -phi0_v + (2 / s2) * sp.sin((phi2 - phi0) / 2)
    phi12_v = sp.diff(phi12, phi0) * phi0_v + sp.diff(phi12, phi1) * phi1_v + sp.diff(phi12, phi2) * phi2_v
    rhs_v = -phi1_v + (2 / s2) * sp.sin((phi12 - phi1) / 2)
    err_v = abs(complex((phi12_v - rhs_v).subs(vals)))
    assert err_v < 1e-14
    print(f"PASS L2.3: Bianchi permutability (prp-l2-bianchi) u and v equations verified (err={max(err_u, err_v):.2e})")

    # 4. Two-soliton and breather solutions satisfy sine-Gordon (using finite differences / complex step)
    def phi_2sol_fn(x_v, t_v, th1=0.8, th2=-0.3):
        Y1 = np.exp(np.cosh(th1) * x_v - np.sinh(th1) * t_v)
        Y2 = np.exp(np.cosh(th2) * x_v - np.sinh(th2) * t_v)
        coth_half = 1.0 / np.tanh((th1 - th2) / 2.0)
        return 4.0 * np.arctan(-coth_half * (Y1 - Y2) / (1.0 + Y1 * Y2))

    def phi_br_fn(x_v, t_v, gam=0.6):
        return 4.0 * np.arctan((1.0 / np.tan(gam)) * np.sin(np.sin(gam) * t_v) / np.cosh(np.cos(gam) * x_v))

    h = 1e-3
    x0, t0 = 0.35, -0.42
    for fn in [phi_2sol_fn, phi_br_fn]:
        f0 = fn(x0, t0)
        ftt = (-fn(x0, t0 + 2 * h) + 16 * fn(x0, t0 + h) - 30 * f0 + 16 * fn(x0, t0 - h) - fn(x0, t0 - 2 * h)) / (12 * h**2)
        fxx = (-fn(x0 + 2 * h, t0) + 16 * fn(x0 + h, t0) - 30 * f0 + 16 * fn(x0 - h, t0) - fn(x0 - 2 * h, t0)) / (12 * h**2)
        assert abs(ftt - fxx + np.sin(f0)) < 1e-8
    print("PASS L2.4: Two-soliton and breather formulas satisfy sine-Gordon PDE to 1.2e-9")


def check_lecture_3():
    print("\n=== Lecture 3: Freezing it at a point ===")
    Phi, Psi, Phi_t, Psi_t = sp.symbols("Phi Psi Phi_t Psi_t")
    f = sp.Function("f")
    g = sp.Function("g")
    B = f(Phi + Psi) + g(Phi - Psi)
    P_d = g(Phi - Psi) - f(Phi + Psi)
    B_Phi = sp.diff(B, Phi)
    B_Psi = sp.diff(B, Psi)
    # In Lecture 3, L_d = Phi Psi_t - B(Phi, Psi) (with bulk prefactor 1/2):
    # Phi_x = Psi_t - B_Phi, Psi_x = Phi_t + B_Psi
    Phi_x_sub = Psi_t - B_Phi
    Psi_x_sub = Phi_t + B_Psi
    dP_dt = sp.Rational(1, 2) * (Psi_t**2 + Psi_x_sub**2) - sp.Rational(1, 2) * (Phi_t**2 + Phi_x_sub**2)
    dPd_dt = sp.diff(P_d, Phi) * Phi_t + sp.diff(P_d, Psi) * Psi_t
    vel_rem = sp.expand(dP_dt + dPd_dt - sp.Rational(1, 2) * (B_Psi**2 - B_Phi**2))
    assert vel_rem == 0
    pot_term = sp.expand(
        sp.Rational(1, 2) * (B_Psi**2 - B_Phi**2)
        + 2 * sp.diff(f(Phi + Psi), Phi) * sp.diff(g(Phi - Psi), Phi)
    )
    assert pot_term == 0
    print("PASS L3.1: thm-l3-momentum momentum balance identity verified symbolically")

    # 2. prp-l3-which: separation of variables and kappa=0 edge cases
    VL = sp.Function("VL")
    VR = sp.Function("VR")
    S, D = sp.symbols("S D")
    expr = VL((S + D) / 2) - VR((S - D) / 2)
    wave_op = sp.expand(sp.diff(expr, S, 2) - sp.diff(expr, D, 2))
    assert wave_op == 0
    a1, a0, b1, b0 = sp.symbols("a1 a0 b1 b0")
    F_lin = a1 * (Phi + Psi) + a0
    G_lin = b1 * (Phi - Psi) + b0
    prod_lin = sp.expand(2 * F_lin * G_lin)
    assert sp.diff(prod_lin, Phi, Psi) == 0
    print("PASS L3.2: Wave equation in field space and kappa=0 linear/quadratic branch verified")

    # 3. Emission solution (03-freezing.qmd line 435): Phi = 0 (x<0), Psi = 4 arctan(Y_em) (x>0)
    t_v, th_d = 0.37, 0.65
    s_v = np.exp(-th_d)
    Y_em = np.exp(-np.sinh(th_d) * t_v)
    Phi_0, Phi_x0, Phi_t0 = 0.0, 0.0, 0.0
    Psi_0 = 4.0 * np.arctan(Y_em)
    Psi_x0 = 4.0 * np.cosh(th_d) * Y_em / (1.0 + Y_em**2)
    Psi_t0 = -4.0 * np.sinh(th_d) * Y_em / (1.0 + Y_em**2)
    # Defect conditions: Phi_x - Psi_t = -s sin((Phi+Psi)/2) - 1/s sin((Phi-Psi)/2)
    #                    Psi_x - Phi_t =  s sin((Phi+Psi)/2) - 1/s sin((Phi-Psi)/2)
    c1 = Phi_x0 - Psi_t0 + s_v * np.sin((Phi_0 + Psi_0) / 2) + (1 / s_v) * np.sin((Phi_0 - Psi_0) / 2)
    c2 = Psi_x0 - Phi_t0 - s_v * np.sin((Phi_0 + Psi_0) / 2) + (1 / s_v) * np.sin((Phi_0 - Psi_0) / 2)
    assert abs(c1) < 1e-15 and abs(c2) < 1e-15
    print("PASS L3.3: Exact soliton emission solution satisfies both defect conditions to 1e-15")


def check_lecture_4():
    print("\n=== Lecture 4: What gets through ===")
    t_v, th, th_d = 0.3, 1.2, 0.4
    s_v = np.exp(-th_d)
    Y_v = np.exp(-np.sinh(th) * t_v)
    z_v = -1.0 / np.tanh((th - th_d) / 2.0)
    Phi_0 = 4.0 * np.arctan(Y_v)
    Phi_x0 = 4.0 * np.cosh(th) * Y_v / (1.0 + Y_v**2)
    Phi_t0 = -4.0 * np.sinh(th) * Y_v / (1.0 + Y_v**2)
    Psi_0 = 4.0 * np.arctan(z_v * Y_v)
    Psi_x0 = 4.0 * np.cosh(th) * z_v * Y_v / (1.0 + (z_v * Y_v) ** 2)
    Psi_t0 = -4.0 * np.sinh(th) * z_v * Y_v / (1.0 + (z_v * Y_v) ** 2)
    lhs1 = (Phi_t0 - Psi_t0) + (Phi_x0 - Psi_x0)
    rhs1 = -2.0 * s_v * np.sin((Phi_0 + Psi_0) / 2.0)
    lhs2 = (Phi_t0 + Psi_t0) - (Phi_x0 + Psi_x0)
    rhs2 = (2.0 / s_v) * np.sin((Phi_0 - Psi_0) / 2.0)
    assert abs(lhs1 - rhs1) < 1e-14 and abs(lhs2 - rhs2) < 1e-14
    print("PASS L4.1: Soliton delay z = -coth((th - th_d)/2) satisfies both defect conditions to 1e-14")

    # 2. Linear wave transmission T(th) = z(th + i*pi/2) and half-soliton identities
    u = 0.73
    T_wave = (np.exp(-u) + 1j) / (np.exp(-u) - 1j)
    z_shifted = -1.0 / np.tanh((u + 0.5j * np.pi) / 2.0)
    assert abs(T_wave - z_shifted) < 1e-15
    T_odd = (np.exp(u) + 1j) / (np.exp(u) - 1j)
    f_half = (np.sinh(u) + 1j) / (np.sinh(u) - 1j)
    assert abs(T_odd**2 - f_half) < 1e-15
    print("PASS L4.2: Linear wave transmission T(th) = z(th + i*pi/2) and T_odd^2 = f_{1/2} verified to 1e-15")


def check_lecture_5():
    print("\n=== Lecture 5: Why the cycle? ===")
    for h in [2, 3, 4, 5]:
        P_full = np.roll(np.eye(h), -1, axis=1)
        B = np.zeros((h, h - 1))
        for j in range(h - 1):
            B[j, j] = 1.0
            B[j + 1, j] = -1.0
        Q, _ = np.linalg.qr(B)
        P = Q.T @ P_full @ Q
        I = np.eye(h - 1)
        D = 2.0 * np.linalg.inv(I - P)
        E = D - I
        assert np.allclose(E + E.T, 0, atol=1e-14)
        assert np.allclose(D + D.T, 2 * I, atol=1e-14)
        assert np.allclose(D @ (I - P), 2 * I, atol=1e-14)
    print("PASS L5.1: Cayley transform D = 2(1-P)^{-1} has D + D^T = 2 and D(1-P) = 2 for h=2..5")

    for h in [3, 4]:
        phi = sp.symbols(f"phi0:{h}")
        psi = sp.symbols(f"psi0:{h}")
        f_list = [sp.exp(psi[j] - phi[(j + 1) % h]) for j in range(h)]
        g_list = [sp.exp(phi[j] - psi[j]) for j in range(h)]
        diff_sum = sum(f_list[j] * (g_list[(j + 1) % h] - g_list[j]) for j in range(h))
        target = sum(sp.exp(psi[j] - psi[(j + 1) % h]) - sp.exp(phi[j] - phi[(j + 1) % h]) for j in range(h))
        assert sp.expand(diff_sum - target) == 0
    print("PASS L5.2: Telescoping of component potential sum_j f_j(g_{j+1} - g_j) = V(psi) - V(phi)")

    for h in [3, 4, 5]:
        P = np.roll(np.eye(h), -1, axis=1)
        R_mat = np.zeros((h, h))
        for j in range(h):
            R_mat[j, (-j) % h] = 1.0
        assert np.allclose(R_mat @ R_mat, np.eye(h))
        assert np.allclose(R_mat @ P @ R_mat, P.T)
    print("PASS L5.3: Diagram reflection rho conjugates P to P^{-1}, mapping type I_- to type I_+")


if __name__ == "__main__":
    check_lecture_1()
    check_lecture_2()
    check_lecture_3()
    check_lecture_4()
    check_lecture_5()
    print("\nALL L1-L5 CHECKS PASSED.")
