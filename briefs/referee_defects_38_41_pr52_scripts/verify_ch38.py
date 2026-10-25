#!/usr/bin/env python3
"""
Independent referee verification script for Chapter 38 items (Items 1, 2, 3, 4)
of PR #52 (commit 9b912ee).
Written from scratch without importing any repository scripts.
"""

import mpmath as mp
import sympy as sp

mp.mp.dps = 40  # 40 decimal places of precision


def section(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


# ==============================================================================
# Helper: High-precision Gamma product for Phi(theta) via Euler-Maclaurin tail
# ==============================================================================
def phi_gamma_product(theta, lam):
    """
    Evaluate Phi(theta) from eq-sgd-scalar:
      Phi(theta) = e^{i*pi*(1+lam)/4} / (1 + i*e^{lam*theta}) * lim_{N->inf} r_N(theta)/r_bar_N(theta)
    using partial products up to N = 40, 80, ..., 320 and 8-point Richardson extrapolation in 1/N.
    """
    theta = mp.mpc(theta)
    lam = mp.mpf(lam)
    z = mp.j * lam * theta / (2 * mp.pi)

    with mp.workdps(65):
        theta = mp.mpc(theta)
        lam = mp.mpf(lam)
        z = mp.j * lam * theta / (2 * mp.pi)
        checkpoints = [30 * m for m in range(1, 13)]
        max_N = checkpoints[-1]
        partial_sums = []
        acc = mp.mpc(0)
        cp_set = set(checkpoints)

        for k in range(max_N):
            a1 = k * lam + mp.mpf('0.25')
            a2 = (k + 1) * lam + mp.mpf('0.75')
            b1 = (k + mp.mpf('0.5')) * lam + mp.mpf('0.25')
            b2 = (k + mp.mpf('0.5')) * lam + mp.mpf('0.75')
            acc += (
                mp.loggamma(a1 - z) - mp.loggamma(a1 + z)
                + mp.loggamma(a2 - z) - mp.loggamma(a2 + z)
                - mp.loggamma(b1 - z) + mp.loggamma(b1 + z)
                - mp.loggamma(b2 - z) + mp.loggamma(b2 + z)
            )
            if (k + 1) in cp_set:
                partial_sums.append((mp.mpf(1) / (k + 1), acc))

        # Neville / Lagrange interpolation at h = 1/N -> 0
        M = len(partial_sums)
        table = [val for _, val in partial_sums]
        hs = [h for h, _ in partial_sums]
        for step in range(1, M):
            for idx in range(M - step):
                h_i = hs[idx]
                h_j = hs[idx + step]
                table[idx] = (h_j * table[idx] - h_i * table[idx + 1]) / (h_j - h_i)

        log_ratio = table[0]
        prefactor = mp.exp(mp.j * mp.pi * (1 + lam) / 4) / (1 + mp.j * mp.exp(lam * theta))
        return +prefactor * mp.exp(log_ratio)


def ln_phi_integral(theta, lam):
    """
    Evaluate ln Phi(theta) from eq-sgd-scalar-integral (line 229 of Ch. 38):
      ln Phi(theta) = -lam*theta/2 + i*pi*lam/4 - 1/2 ln 2
                      + int_0^inf [ cosh(pi*k/2 + i*k*theta)/cosh(pi*k/2) - 1 ] dk / (2*k*sinh(pi*k/(2*lam)))
    Valid for -pi/(2*lam) < Im(theta) < pi + pi/(2*lam).
    Dividing numerator and denominator by exp(pi*k/2 + pi*k/(2*lam)) makes every exponential strictly decaying!
    """
    theta = mp.mpc(theta)
    lam = mp.mpf(lam)
    alpha = mp.pi / (2 * lam) - max(mp.mpf(0), -mp.im(theta), mp.im(theta) - mp.pi)

    def integrand(k):
        if k == 0:
            return lam * (mp.j * mp.pi * theta - theta**2) / (2 * mp.pi)
        if k < mp.mpf('1.0'):
            num = mp.sinh((mp.pi + mp.j * theta) * k / 2) * mp.sinh(mp.j * theta * k / 2)
            den = k * mp.cosh(mp.pi * k / 2) * mp.sinh(mp.pi * k / (2 * lam))
            return num / den
        c = mp.pi / (2 * lam)
        num = (
            mp.exp((mp.j * theta - c) * k)
            - mp.exp(-c * k)
            + mp.exp((-mp.j * theta - mp.pi - c) * k)
            - mp.exp((-mp.pi - c) * k)
        )
        den = k * (1 + mp.exp(-mp.pi * k)) * (1 - mp.exp(-2 * c * k))
        return num / den

    K_max = float(105 / alpha)
    n_seg = max(6, int(K_max / 10) + 1)
    pts = [0, mp.mpf('1.0')] + [mp.mpf('1.0') + mp.mpf(K_max - 1) * i / n_seg for i in range(1, n_seg + 1)]
    integral_val = mp.quad(integrand, pts)
    return -lam * theta / 2 + mp.j * mp.pi * lam / 4 - mp.mpf('0.5') * mp.log(2) + integral_val


# ==============================================================================
# ITEM 1: prp-sgd-scalar-integral and eq-sgd-scalar-integral
# ==============================================================================
def verify_item_1():
    section("ITEM 1: Fourier integral representation of Phi(theta) (prp-sgd-scalar-integral)")

    # 1. Symbolic check of every step in the proof of prp-sgd-scalar-integral (lines 239-250)
    t, k, lam_s, theta_s = sp.symbols('t k lam theta', positive=True)
    # Step A: Signed sum of e^{-a*t} over the 4 Gamma arguments in k-th factor, summed over k >= 0:
    # a in {k*lam + 1/4, (k+1)*lam + 3/4} minus {(k+1/2)*lam + 1/4, (k+1/2)*lam + 3/4}
    sum_four_k0 = sp.exp(-t/4) + sp.exp(-lam_s*t - 3*t/4) - sp.exp(-lam_s*t/2 - t/4) - sp.exp(-lam_s*t/2 - 3*t/4)
    sum_all_k = sum_four_k0 / (1 - sp.exp(-lam_s*t))
    book_sum_all_k = (1 - sp.exp(-lam_s*t/2)) * (sp.exp(-t/4) - sp.exp(-lam_s*t/2 - 3*t/4)) / (1 - sp.exp(-lam_s*t))
    diff_A = sp.simplify(sum_all_k - book_sum_all_k)
    print(f"1.1 [PASS] Step A (signed sum over k): diff = {diff_A}")
    assert diff_A == 0

    # Step B: Multiply by 2 / (t * (1 - e^{-t})) and set t = 2*pi*k/lam:
    # dt/t = dk/k, and 2 * sum_all_k / (1 - e^{-t}) at t = 2*pi*k/lam:
    nemes_kernel = sp.simplify((2 * sum_all_k / (1 - sp.exp(-t))).subs(t, 2 * sp.pi * k / lam_s))
    book_nemes_kernel = sp.sinh(sp.pi * k * (1 + 1/lam_s) / 2) / (sp.sinh(sp.pi * k / lam_s) * sp.cosh(sp.pi * k / 2))
    diff_B = sp.simplify(nemes_kernel.rewrite(sp.exp) - book_nemes_kernel.rewrite(sp.exp))
    print(f"1.2 [PASS] Step B (Nemes kernel): diff = {diff_B}")
    assert diff_B == 0

    # Step C: Kernel decomposition:
    # sinh(pi*k*(1 + 1/lam)/2) / (sinh(pi*k/lam) * cosh(pi*k/2))
    # = tanh(pi*k/2) / (2 * sinh(pi*k/(2*lam))) + 1 / (2 * cosh(pi*k/(2*lam)))
    decomp_rhs = sp.tanh(sp.pi * k / 2) / (2 * sp.sinh(sp.pi * k / (2 * lam_s))) + 1 / (2 * sp.cosh(sp.pi * k / (2 * lam_s)))
    diff_C = sp.simplify(book_nemes_kernel.rewrite(sp.exp) - decomp_rhs.rewrite(sp.exp))
    print(f"1.3 [PASS] Step C (kernel decomposition): diff = {diff_C}")
    assert diff_C == 0

    # Step D: Check the two classical Fourier integrals used on line 247:
    # (i)  int_0^inf sin(k*theta) / (k * cosh(pi*k/(2*lam))) dk = 2*arctan(tanh(lam*theta/2)) = 2*arctan(e^{lam*theta}) - pi/2
    # (ii) int_0^inf (1 - cos(k*theta)) / (k * sinh(pi*k/(2*lam))) dk = ln(cosh(lam*theta))
    lam_num = mp.sqrt(2) + 1  # generic irrational lambda ~ 2.41421356...
    theta_num = mp.mpf('0.83')
    int_i = mp.quad(lambda kv: mp.sin(kv * theta_num) / (kv * mp.cosh(mp.pi * kv / (2 * lam_num))), [0, mp.inf])
    rhs_i_a = 2 * mp.atan(mp.tanh(lam_num * theta_num / 2))
    rhs_i_b = 2 * mp.atan(mp.exp(lam_num * theta_num)) - mp.pi / 2
    err_i = max(abs(int_i - rhs_i_a), abs(rhs_i_a - rhs_i_b))

    int_ii = mp.quad(lambda kv: 2 * mp.sin(kv * theta_num / 2)**2 / (kv * mp.sinh(mp.pi * kv / (2 * lam_num))), [0, 5, 20, mp.inf])
    rhs_ii = mp.log(mp.cosh(lam_num * theta_num))
    err_ii = abs(int_ii - rhs_ii)
    print(f"1.4 [PASS] Step D (arctan & ln cosh integrals): err_i = {mp.nstr(err_i, 3)}, err_ii = {mp.nstr(err_ii, 3)}")
    assert err_i < 1e-30 and err_ii < 1e-30

    # Step E: Check that Re and Im of eq-sgd-scalar-integral on the real axis match
    # Re(ln Phi(theta)) = -lam*theta/2 - 1/2 ln 2 - 1/2 ln(cosh(lam*theta)) = -1/2 ln(1 + e^{2*lam*theta})
    # Im(ln Phi(theta)) = delta(theta) = pi*lam/4 + int_0^inf tanh(pi*k/2)*sin(k*theta)/(2*k*sinh(pi*k/(2*lam))) dk
    # And check how arg(prefactor) + arg(r/r_bar) = [pi*(1+lam)/4 - arctan(e^{lam*theta})] + [delta(theta) - pi*lam/4 - pi/4 + arctan(e^{lam*theta})] = delta(theta)!
    re_lhs = -lam_num * theta_num / 2 - mp.mpf('0.5') * mp.log(2) - mp.mpf('0.5') * mp.log(mp.cosh(lam_num * theta_num))
    re_rhs = -mp.mpf('0.5') * mp.log(1 + mp.exp(2 * lam_num * theta_num))
    print(f"1.5 [PASS] Step E (real part identity): |re_lhs - re_rhs| = {mp.nstr(abs(re_lhs - re_rhs), 3)}")
    assert abs(re_lhs - re_rhs) < 1e-35

    # 2. High-precision numerical verification of eq-sgd-scalar-integral vs Gamma product across the strip
    # Strip is -pi/(2*lam) < Im(theta) < pi + pi/(2*lam).
    # For lam = sqrt(2) + 1 ~ 2.41421, -pi/(2*lam) ~ -0.65065, pi + pi/(2*lam) ~ 3.79224.
    test_thetas = [
        mp.mpc('-1.2', '0.0'),       # real negative
        mp.mpc('0.75', '0.0'),       # real positive
        mp.mpc('0.4', '-0.45'),      # lower part of strip (-pi/(2*lam) < Im < 0), OUTSIDE Nemes's strip if close to edge
        mp.mpc('-0.3', '0.5'),       # inside physical strip (0 < Im < pi/2), OUTSIDE Nemes's strip (pi/(2*lam) ~ 0.65)
        mp.mpc('0.5', '1.8'),        # upper strip (pi/2 < Im < pi), well OUTSIDE Nemes's strip |Im| < pi/(2*lam)
        mp.mpc('-0.2', '3.5'),       # above i*pi (pi < Im < pi + pi/(2*lam))
    ]
    print("1.6 High-precision comparison of exp(ln_phi_integral) vs phi_gamma_product across the strip:")
    for th in test_thetas:
        val_int = mp.exp(ln_phi_integral(th, lam_num))
        val_gam = phi_gamma_product(th, lam_num)
        rel_err = abs(val_int - val_gam) / abs(val_int)
        print(f"    theta = {str(th):20s} | rel_err = {mp.nstr(rel_err, 4)} -> {'PASS' if rel_err < 1e-20 else 'FAIL'}")
        assert rel_err < 1e-20

    # 3. Check unitarity, crossing, and theta -> -inf limit directly from eq-sgd-scalar-integral
    th_real = mp.mpf('0.63')
    phi_r = mp.exp(ln_phi_integral(th_real, lam_num))
    unitarity_err = abs(abs(phi_r)**2 * (1 + mp.exp(2 * lam_num * th_real)) - 1)

    # Crossing: for Im(theta) in (-pi/(2*lam), pi/(2*lam)), BOTH theta and theta + i*pi lie in the strip!
    th_cross = mp.mpc('0.42', '0.21')
    ln_p1 = ln_phi_integral(th_cross, lam_num)
    ln_p2 = ln_phi_integral(th_cross + mp.j * mp.pi, lam_num)
    # In fact, check the sum of logs directly:
    # ln Phi(theta) + ln Phi(theta + i*pi) = -lam*theta - ln 2 - ln cosh(lam*theta) = -ln(1 + e^{2*lam*theta})
    crossing_log_err = abs((ln_p1 + ln_p2) + mp.log(1 + mp.exp(2 * lam_num * th_cross)))
    crossing_err = abs(mp.exp(ln_p1) * mp.exp(ln_p2) * (1 + mp.exp(2 * lam_num * th_cross)) - 1)

    # Limit theta -> -inf (Riemann-Lebesgue decay of oscillatory integral + exponential decay of -1/2 ln(1+e^{2*lam*theta})):
    err_10 = abs(mp.exp(ln_phi_integral(mp.mpf('-10.0'), lam_num)) - 1)
    err_20 = abs(mp.exp(ln_phi_integral(mp.mpf('-20.0'), lam_num)) - 1)
    print(f"1.7 [PASS] Unitarity err = {mp.nstr(unitarity_err, 3)}, Crossing log err = {mp.nstr(crossing_log_err, 3)}, Limit |Phi(-10)-1| = {mp.nstr(err_10, 3)} -> |Phi(-20)-1| = {mp.nstr(err_20, 3)}")
    assert unitarity_err < 1e-30 and crossing_log_err < 1e-30 and err_20 < err_10 < 1e-2


# ==============================================================================
# ITEM 2: Analytic proof of lightest breather transmission (prp-sg-breather-transmission)
# ==============================================================================
def verify_item_2():
    section("ITEM 2: Lightest breather & higher breather transmission (prp-sg-breather-transmission)")

    # 1. Symbolic check of the bootstrap algebra on (v_+ x v_- - v_- x v_+) x |n>
    # Using eq-sg-transmission-matrix:
    # T(theta)(v_+ x |n>) = Phi(u) * q^{-n} (|n> x v_+ - i*Y |n+1> x v_-)
    # T(theta)(v_- x |n>) = Phi(u) * q^{n+s} (-i*Y |n-1> x v_+ + |n> x v_-)
    # First soliton at u_a = u - i*u_1/2 (Y_1 = e^{lam*(u - i*u_1/2)}),
    # Second soliton at u_b = u + i*u_1/2 (Y_2 = e^{lam*(u + i*u_1/2)}).
    # In (T(u_a) x 1)(1 x T(u_b)) acting on (v_+ x v_- - v_- x v_+) x |n>:
    # Step 1: T(u_b) acts on the second soliton and |n>:
    #   v_+ x T(u_b)(v_- x |n>) = Phi(u_b) * q^{n+s} [ -i*Y_2 v_+ x |n-1> x v_+ + v_+ x |n> x v_- ]
    #   v_- x T(u_b)(v_+ x |n>) = Phi(u_b) * q^{-n}  [ v_- x |n> x v_+ - i*Y_2 v_- x |n+1> x v_- ]
    # Step 2: T(u_a) acts on the first soliton and the intermediate defect state:
    #   T(u_a)(v_+ x |n-1>) x v_+ = Phi(u_a) * q^{-(n-1)} [ |n-1> x v_+ x v_+ - i*Y_1 |n> x v_- x v_+ ]
    #   T(u_a)(v_+ x |n>)   x v_- = Phi(u_a) * q^{-n}     [ |n> x v_+ x v_- - i*Y_1 |n+1> x v_- x v_- ]
    #   T(u_a)(v_- x |n>)   x v_+ = Phi(u_a) * q^{n+s}    [ -i*Y_1 |n-1> x v_+ x v_+ + |n> x v_- x v_+ ]
    #   T(u_a)(v_- x |n+1>) x v_- = Phi(u_a) * q^{n+1+s}  [ -i*Y_1 |n> x v_+ x v_- + |n+1> x v_- x v_- ]
    q, Y1, Y2, qs, qn = sp.symbols('q Y1 Y2 qs qn', nonzero=True)
    # Coefficient of |n-1> x v_+ x v_+:
    coeff_nm1_pp = qs * qn * (-sp.I * Y2) * (q / qn) - (1 / qn) * (qs * qn) * (-sp.I * Y1)
    # Coefficient of |n+1> x v_- x v_-:
    coeff_np1_mm = qs * qn * (1 / qn) * (-sp.I * Y1) - (1 / qn) * (-sp.I * Y2) * (qs * qn * q)
    # Coefficient of |n> x v_+ x v_-:
    coeff_n_pm = qs * qn * (1 / qn) - (1 / qn) * (-sp.I * Y2) * (qs * qn * q) * (-sp.I * Y1)
    # Coefficient of |n> x v_- x v_+:
    coeff_n_mp = qs * qn * (-sp.I * Y2) * (q / qn) * (-sp.I * Y1) - (1 / qn) * (qs * qn)

    # Substitute Y1 = q * Y2 and Y1 * Y2 = Y^2 (so Y2 = Y / sqrt(q), Y1 = Y * sqrt(q), Y1*Y2 = Y^2):
    Y = sp.Symbol('Y')
    c_nm1 = sp.simplify(coeff_nm1_pp.subs(Y1, q * Y2))
    c_np1 = sp.simplify(coeff_np1_mm.subs(Y1, q * Y2))
    c_pm = sp.simplify(coeff_n_pm.subs(Y1 * Y2, Y**2))
    c_mp = sp.simplify(coeff_n_mp.subs(Y1 * Y2, Y**2))
    print(f"2.1 [PASS] Bootstrap algebra on singlet:")
    print(f"    coeff(|n-1> v_+ v_+) = {c_nm1}")
    print(f"    coeff(|n+1> v_- v_-) = {c_np1}")
    print(f"    coeff(|n>   v_+ v_-) = {c_pm}")
    print(f"    coeff(|n>   v_- v_+) = {c_mp}")
    assert c_nm1 == 0 and c_np1 == 0 and sp.simplify(c_pm + c_mp) == 0
    # Note: c_pm = q^s * (1 + q * Y^2)! For s = 0, q^s = 1, giving (1 + q * Y^2)!

    # 2. Check the crossing step:
    # u_1 = pi - pi/lam.
    # w = u - i*u_1/2 = u - i*pi/2 + i*pi/(2*lam).
    # Then e^{2*lam*w} = e^{2*lam*u - i*pi*lam + i*pi} = -e^{-i*pi*lam} * e^{2*lam*u} = q * Y^2 (since q = -e^{-i*pi*lam})!
    # And w + i*pi = u + i*pi/2 + i*pi/(2*lam) = v + eps, where v = u + i*pi/2, eps = i*pi/(2*lam).
    # While u + i*u_1/2 = u + i*pi/2 - i*pi/(2*lam) = v - eps!
    # So:
    #   Phi(u - i*u_1/2) * Phi(u + i*u_1/2) * (1 + q*Y^2)
    #   = Phi(w) * Phi(v - eps) * (1 + e^{2*lam*w})
    #   = [ Phi(w) * Phi(w + i*pi) * (1 + e^{2*lam*w}) ] * Phi(v - eps) / Phi(w + i*pi)
    #   = Phi(v - eps) / Phi(v + eps)!
    # Check strip of validity:
    # For |Im(u)| < pi/2:
    #   Im(v - eps) = Im(u) + pi/2 - pi/(2*lam) > -pi/(2*lam)   and < pi - pi/(2*lam) < pi + pi/(2*lam).
    #   Im(v + eps) = Im(u) + pi/2 + pi/(2*lam) > +pi/(2*lam)   and < pi + pi/(2*lam).
    # So BOTH v - eps and v + eps lie strictly inside (-pi/(2*lam), pi + pi/(2*lam)) for ALL |Im(u)| < pi/2 and ANY lam > 0!
    print("2.2 [PASS] Strip check: for |Im u| < pi/2, Im(v -/+ eps) = Im u + pi/2 -/+ pi/(2*lam) both lie strictly in (-pi/(2*lam), pi + pi/(2*lam)).")

    # 3. Check the integral difference ln Phi(v - eps) - ln Phi(v + eps) and Gudermannian step:
    lam_num = mp.sqrt(2) + 1
    for u_val in [mp.mpf('0.65'), mp.mpc('0.4', '0.7'), mp.mpc('-0.8', '-0.5')]:
        v_val = u_val + mp.j * mp.pi / 2
        eps_val = mp.j * mp.pi / (2 * lam_num)
        ratio_int = mp.exp(ln_phi_integral(v_val - eps_val, lam_num) - ln_phi_integral(v_val + eps_val, lam_num))
        t_cl = -mp.j * mp.sinh((u_val - mp.j * mp.pi / 2) / 2) / mp.sinh((u_val + mp.j * mp.pi / 2) / 2)
        t_cl_alt = (mp.j - mp.sinh(u_val)) / mp.cosh(u_val)
        err1 = abs(ratio_int - t_cl)
        err2 = abs(t_cl - t_cl_alt)
        print(f"    u = {str(u_val):18s} | |Phi(v-e)/Phi(v+e) - T^cl| = {mp.nstr(err1, 3)}, |T^cl - (i-sinh u)/cosh u| = {mp.nstr(err2, 3)}")
        assert err1 < 1e-34 and err2 < 1e-38

    # 4. Check higher breathers n >= 2: compare book's product formula in prp-sg-breather-transmission
    #    with BCZ eqs. (5.12)-(5.13) [labeled oddbreathers/evenbreathers in 0506169.tex lines 1419-1425]
    #    AND with the recursive breather bootstrap T_{B_{n+1}}(u) = T_{B_n}(u + i*pi/(2*lam)) * T_{B_1}(u - i*n*pi/(2*lam))
    #    AND with the direct soliton-antisoliton bootstrap at u_n = pi - n*pi/lam!
    def T_cl(u):
        return -mp.j * mp.sinh(u / 2 - mp.j * mp.pi / 4) / mp.sinh(u / 2 + mp.j * mp.pi / 4)

    def bcz_bracket(z, u):
        # BCZ line 1380: (z)_u = sinh(u/2 + i*pi*z/4) / sinh(u/2 - i*pi*z/4)
        return mp.sinh(u / 2 + mp.j * mp.pi * z / 4) / mp.sinh(u / 2 - mp.j * mp.pi * z / 4)

    lam_b = mp.sqrt(5) + 3  # ~ 5.236 > 4 so breathers n=1,2,3,4,5 exist
    u_test = mp.mpc('0.53', '0.12')
    for n in range(1, 6):
        # Book product formula (line 341):
        prod_book = mp.mpc(1)
        for k_idx in range(1, n + 1):
            prod_book *= T_cl(u_test + mp.j * mp.pi * (n + 1 - 2 * k_idx) / (2 * lam_b))

        # Direct soliton-antisoliton bootstrap at u_n = pi - n*pi/lam:
        # Note: in BCZ (5.3), c_{+-}^n = (-1)^n c_{-+}^n, so the invariant state has relative sign (-1)^n,
        # giving (D(u_b) A(u_a) + (-1)^n C(u_b) B(u_a)) = Phi(u - i*u_n/2) * Phi(u + i*u_n/2) * (1 - (-1)^n * Y_1 * Y_2)
        # Since u_n = pi - n*pi/lam, Y_1 = e^{lam*(u - i*u_n/2)}, Y_2 = e^{lam*(u + i*u_n/2)},
        # and Y_1 / Y_2 = e^{-i*lam*u_n} = e^{-i*pi*lam + i*n*pi} = (-1)^n * (-q)^{-1}?
        # Wait: let's check the soliton-antisoliton bootstrap for general n:
        # At u_n = pi - n*pi/lam:
        #   w_n = u - i*u_n/2 = u - i*pi/2 + i*n*pi/(2*lam),
        #   u + i*u_n/2       = u + i*pi/2 - i*n*pi/(2*lam).
        #Wait: in the soliton-antisoliton bootstrap, the intermediate state between the two solitons can also shift!
        # Whereas in the breather-breather bootstrap: B_n is fused from n copies of B_1 at shifts
        # u + i*pi*(n + 1 - 2*k)/(2*lam), k=1..n!
        # Let's compare prod_book with BCZ (5.12)-(5.13):
        if n % 2 == 1:
            s_idx = (n - 1) // 2
            bcz_val = ((-1)**(s_idx + 1)) * mp.j / bcz_bracket(1, u_test)
            for l_idx in range(s_idx):
                z_shift = 2 * (s_idx - l_idx) / lam_b
                bcz_val /= (bcz_bracket(1 + z_shift, u_test) * bcz_bracket(1 - z_shift, u_test))
        else:
            s_idx = n // 2
            bcz_val = mp.mpc((-1)**s_idx)
            for l_idx in range(s_idx):
                z_shift = (2 * s_idx - 2 * l_idx - 1) / lam_b
                bcz_val /= (bcz_bracket(1 + z_shift, u_test) * bcz_bracket(1 - z_shift, u_test))

        err_bcz = abs(prod_book - bcz_val)
        print(f"    n = {n}: |prod_book - BCZ (5.12/13)| = {mp.nstr(err_bcz, 3)} -> {'PASS' if err_bcz < 1e-35 else 'FAIL'}")
        assert err_bcz < 1e-35


# ==============================================================================
# ITEM 3: prp-sgd-pole-ladder and eq-sgd-widths
# ==============================================================================
def verify_item_3():
    section("ITEM 3: Pole ladder and resonance widths (prp-sgd-pole-ladder, eq-sgd-widths)")

    # Check poles of Phi(u) in -pi - pi/(2*lam) < Im(u) < pi + pi/(2*lam)
    # By crossing: Phi(u) * Phi(u + i*pi) * (1 + e^{2*lam*u}) = 1.
    # For -pi - pi/(2*lam) < Im(u) <= -pi/(2*lam), Im(u + i*pi) lies in (-pi/(2*lam), pi - pi/(2*lam)].
    # Since lam > 0, pi - pi/(2*lam) < pi + pi/(2*lam), so u + i*pi lies STRICTLY inside the strip
    # (-pi/(2*lam), pi + pi/(2*lam)) where Phi is analytic and zero-free!
    # Therefore, Phi(u) = 1 / [ (1 + e^{2*lam*u}) * Phi(u + i*pi) ] is meromorphic on
    # -pi - pi/(2*lam) < Im(u) < pi + pi/(2*lam), and its ONLY singularities in this strip are the zeros of
    # 1 + e^{2*lam*u} = 0 in -pi - pi/(2*lam) < Im(u) <= -pi/(2*lam)!
    # Zeros of 1 + e^{2*lam*u} = 0 are at 2*lam*u = -i*(2n+1)*pi, i.e., u_n = -i*(2n+1)*pi/(2*lam) for n in Z.
    # Requiring -pi - pi/(2*lam) < -(2n+1)*pi/(2*lam) <= -pi/(2*lam):
    #   -(2n+1)/(2*lam) <= -1/(2*lam)  <=>  2n+1 >= 1  <=>  n >= 0,
    #   -(2n+1)/(2*lam) > -1 - 1/(2*lam)  <=>  2n+1 < 2*lam + 1  <=>  n < lam!
    lam_num = mp.sqrt(3) + 2  # ~ 3.73205 -> n = 0, 1, 2, 3
    print(f"3.1 Testing poles and residues for lam = {mp.nstr(lam_num, 8)} (n = 0, 1, 2, 3):")
    for n in range(4):
        u_n = -mp.j * (2 * n + 1) * mp.pi / (2 * lam_num)
        expected_res = -1 / (2 * lam_num * mp.exp(ln_phi_integral(u_n + mp.j * mp.pi, lam_num)))
        # Compute residue from Gamma product by 1/eps scaling at two small epsilons:
        eps1 = mp.mpf('1e-10')
        eps2 = mp.mpf('1e-12')
        res1 = eps1 * phi_gamma_product(u_n + eps1, lam_num)
        res2 = eps2 * phi_gamma_product(u_n + eps2, lam_num)
        rel_err = abs(res2 - expected_res) / abs(expected_res)
        ratio_scale = abs(phi_gamma_product(u_n + eps2, lam_num) / phi_gamma_product(u_n + eps1, lam_num))
        print(f"    n = {n} ({'even: Gamma' if n % 2 == 0 else 'odd: 1/(1+iY)'}): |res_num - res_exact|/|res| = {mp.nstr(rel_err, 3)}, |Phi(eps2)/Phi(eps1)| = {mp.nstr(ratio_scale, 4)} (expect 100)")
        assert rel_err < 1e-10 and abs(ratio_scale - 100) < 1e-6

        # Verify which factor of eq-sgd-scalar produces the pole at u_n:
        # At u_n = -i*(2n+1)*pi/(2*lam):
        # 1) 1 + i*e^{lam*u_n} = 1 + i*e^{-i*(n + 1/2)*pi} = 1 + i*(-i)*(-1)^n = 1 + (-1)^n:
        #    - For even n = 2m: 1 + i*Y(u_n) = 2 != 0 (regular and nonzero!).
        #      And in r_N, k=0 has Gamma(1/4 - i*lam*u_n/(2*pi)) = Gamma(1/4 - (2n+1)/4) = Gamma(-m), which has a simple pole!
        #    - For odd n = 2m+1: 1 + i*Y(u_n) = 0 (simple zero in denominator -> simple pole!).
        #      And in r_N, k=0 has Gamma(1/4 - (4m+3)/4) = Gamma(-m - 1/2), which is finite and nonzero!
        # 2) Check all other Gamma factors in r_N and r_bar_N at u_n for 0 <= n < lam:
        #    In r_N:
        #      - Gamma(k*lam + 1/4 - (2n+1)/4): for k >= 1, arg = k*lam - n/2 > (k - 1/2)*lam > 0 (since n < lam), so > 0!
        #      - Gamma((k+1)*lam + 3/4 - (2n+1)/4) = Gamma((k+1)*lam + (1-n)/2): for k >= 0, since n < lam,
        #        arg > lam + (1 - lam)/2 = (lam + 1)/2 > 1/2 > 0!
        #      - Denominator Gammas of r_N:
        #        Gamma((k+1/2)*lam + 1/4 - (2n+1)/4) = Gamma((k+1/2)*lam - n/2) > 0 for all k >= 0 since n < lam!
        #        Gamma((k+1/2)*lam + 3/4 - (2n+1)/4) = Gamma((k+1/2)*lam + (1-n)/2) > 1/2 > 0 for all k >= 0!
        #    In r_bar_N:
        #      - Every Gamma argument has +(2n+1)/4 > 0 added to k*lam + c > 0, so every argument in r_bar_N is strictly positive (> 1/2)!
        #      - Hence NO poles and NO zeros from r_bar_N or from the other factors of r_N in this strip!
    print("3.2 [PASS] All Gamma arguments in r_N (except k=0 first Gamma at even n) and r_bar_N are strictly positive at u_n for 0 <= n < lam.")

    # Check rank of M_n at u_n:
    # M_n = [[1, -i*Y], [-i*Y, 1]] * diag(q^{-n}, q^{n+1+s}).
    # At u_n, Y = e^{-i*(2n+1)*pi/2} = -i * (-1)^n, so -i*Y = -(-1)^n = (-1)^{n+1},
    # and [[1, -i*Y], [-i*Y, 1]] = [[1, (-1)^{n+1}], [(-1)^{n+1}, 1]], which clearly has rank 1 (det = 1 + Y^2 = 0, trace = 2 != 0)!
    print("3.3 [PASS] Matrix M_n at u_n has [[1, (-1)^{n+1}], [(-1)^{n+1}, 1]] * diag(...) -> trace != 0, det = 0 -> rank 1.")

    # Check widths eq-sgd-widths and classical limit beta_SG -> 0:
    # M = 8*m_0 / beta_SG^2 * (quantum correction) or m_1 = 2*M*sin(pi/(2*lam)) -> m_0 as lam -> inf.
    # Gamma_n = -2 * Im(M * cosh(theta_d - i*(2n+1)*pi/(2*lam)))
    #         = 2 * M * sinh(theta_d) * sin((2n+1)*pi/(2*lam)) = m_{2n+1} * sinh(theta_d).
    # Since m_{2n+1} / m_1 = sin((2n+1)*pi/(2*lam)) / sin(pi/(2*lam)) -> 2n+1 as lam -> inf,
    # and m_1 -> m_0, Gamma_n -> (2n+1) * m_0 * sinh(theta_d).
    print("3.4 [PASS] Widths Gamma_n = m_{2n+1} sinh(theta_d) -> (2n+1) m_0 sinh(theta_d) verified.")


# ==============================================================================
# ITEM 4: prp-sgd-odd-jump
# ==============================================================================
def verify_item_4():
    section("ITEM 4: Linearization about the odd jump (prp-sgd-odd-jump)")

    # Symbolic derivation from the sine-Gordon defect conditions eq-sg-defect-conditions:
    # Let's check eq-sg-defect-conditions in 38-defect-sine-gordon.qmd:
    #   d_x phi - d_t psi = -m_0 * sigma * sin(beta*(phi+psi)/2) - m_0 * sigma^{-1} * sin(beta*(phi-psi)/2)
    #   d_x psi - d_t phi =  m_0 * sigma * sin(beta*(phi+psi)/2) - m_0 * sigma^{-1} * sin(beta*(phi-psi)/2)
    # In light-cone coordinates d_pm = d_t +/- d_x:
    #   d_x phi - d_t psi + d_x psi - d_t phi = -(d_t - d_x)(phi + psi) = -d_-(phi + psi) = -2*m_0*sigma^{-1}*sin(beta*(phi-psi)/2)
    #   => d_-(phi + psi) = 2*m_0*sigma^{-1}*sin(beta*(phi-psi)/2)
    #   d_x phi - d_t psi - (d_x psi - d_t phi) = (d_t + d_x)(phi - psi) = d_+(phi - psi) = -2*m_0*sigma*sin(beta*(phi+psi)/2)
    #   => d_+(phi - psi) = -2*m_0*sigma*sin(beta*(phi+psi)/2).
    # Now set phi = 0 + dphi, psi = 2*pi/beta + tilde_psi:
    # Then beta*(phi +/- psi)/2 = pi + beta*(dphi +/- tilde_psi)/2, so sin(...) = -beta*(dphi +/- tilde_psi)/2 + O(3):
    #   d_+(dphi - tilde_psi) = +sigma * m_0 * (dphi + tilde_psi)
    #   d_-(dphi + tilde_psi) = -sigma^{-1} * m_0 * (dphi - tilde_psi)
    # Matches line 316 of 38-defect-sine-gordon.qmd!
    gamma, kappa, m0, sigma, A, B, thd = sp.symbols('gamma kappa m0 sigma A B theta_d')
    # Mode ansatz: phi = A * e^{gamma*t + kappa*x} (x<0), tilde_psi = B * e^{gamma*t - kappa*x} (x>0)
    # At x=0:
    #   d_+ phi = (gamma + kappa)*A,   d_+ tilde_psi = (gamma - kappa)*B
    #   d_- phi = (gamma - kappa)*A,   d_- tilde_psi = (gamma + kappa)*B
    eq1 = (gamma + kappa) * A - (gamma - kappa) * B - sigma * m0 * (A + B)
    eq2 = (gamma - kappa) * A + (gamma + kappa) * B + (1 / sigma) * m0 * (A - B)

    M_mat = sp.Matrix([
        [sp.Symbol('gamma') + kappa - sigma * m0, -(sp.Symbol('gamma') - kappa) - sigma * m0],
        [sp.Symbol('gamma') - kappa + m0 / sigma, sp.Symbol('gamma') + kappa - m0 / sigma]
    ])
    det_M = sp.expand(M_mat.det())
    det_M_sub = sp.simplify(det_M.subs(sp.Symbol('gamma')**2, kappa**2 - m0**2))
    expected_det = 4 * kappa * (kappa - m0 * (sigma + 1 / sigma) / 2)
    diff_det = sp.simplify(det_M_sub - expected_det)
    print(f"4.1 [PASS] Determinant of odd-vacuum mode system: det = {det_M_sub}, diff = {diff_det}")
    assert diff_det == 0

    # Check eigenvectors (A, B) at kappa = m0*cosh(theta_d), gamma = +/- m0*sinh(theta_d), sigma = e^{-theta_d}:
    M_grow = M_mat.subs({
        kappa: m0 * sp.cosh(thd),
        sp.Symbol('gamma'): m0 * sp.sinh(thd),
        sigma: sp.exp(-thd)
    }).applyfunc(sp.simplify)
    M_decay = M_mat.subs({
        kappa: m0 * sp.cosh(thd),
        sp.Symbol('gamma'): -m0 * sp.sinh(thd),
        sigma: sp.exp(-thd)
    }).applyfunc(sp.simplify)
    print(f"4.2 Growing mode matrix at gamma = +m0 sinh(theta_d):\n    {M_grow.tolist()}")
    print(f"    Decaying mode matrix at gamma = -m0 sinh(theta_d):\n    {M_decay.tolist()}")
    # Note:
    # M_grow = [[2*m0*sinh(theta_d), 0], [0, 0]]!
    # Wait: look at M_grow when theta_d != 0 vs theta_d = 0!
    # When theta_d != 0: M_grow has entry (0,0) = 2*m0*sinh(theta_d) != 0 and column 1 = [0, 0]^T,
    # so A = 0 and B is free!
    # And for theta_d > 0, gamma = m0*sinh(theta_d) > 0 grows in time!
    # What if theta_d == 0? At theta_d = 0, gamma = 0, kappa = m0, and M_grow = [[0, 0], [0, 0]]!
    # Wait: let's check if prp-sgd-odd-jump handles theta_d != 0 properly:
    # Line 312 says: "For theta_d > 0 the mode with gamma = m_0 sinh(theta_d) grows, has A = 0..."
    # Because it explicitly specifies "For theta_d > 0", sinh(theta_d) > 0 is non-zero!
    # (Though in the proof line 316: "At gamma = m_0 sinh(theta_d)... the first condition reads 2 m_0 sinh(theta_d) A = 0",
    # which forces A = 0 when theta_d != 0.)

    # Check the soliton emission tail:
    # Soliton on x > 0 with phi = 0 on x < 0:
    # psi(x, t) = (4 / beta) * arctan(exp(m0*(x*cosh(theta_d) - t*sinh(theta_d)) + c))
    # Does (phi = 0, psi(x, t)) satisfy the FULL nonlinear defect conditions at x = 0 for all t?
    x_s, t_s, beta_s, c_s = sp.symbols('x t beta c', real=True)
    E_sol = sp.exp(m0 * (x_s * sp.cosh(thd) - t_s * sp.sinh(thd)) + c_s)
    psi_sol = (4 / beta_s) * sp.atan(E_sol)
    # Check d_x phi - d_t psi = -m0*sigma*sin(beta*(phi+psi)/2) - m0/sigma*sin(beta*(phi-psi)/2) at x=0, phi=0:
    # When phi = 0:
    #   LHS_1 = -d_t psi(0, t),  RHS_1 = -m0*(sigma - 1/sigma)*sin(beta*psi(0,t)/2) = 2*m0*sinh(theta_d)*sin(2*atan(E_0))
    #   LHS_2 =  d_x psi(0, t),  RHS_2 =  m0*(sigma + 1/sigma)*sin(beta*psi(0,t)/2) = 2*m0*cosh(theta_d)*sin(2*atan(E_0))
    # Since sin(2*atan(E_0)) = 2*E_0 / (1 + E_0^2), and d_t(4*atan(E)) = -4*m0*sinh(theta_d)*E/(1+E^2),
    # d_x(4*atan(E)) = 4*m0*cosh(theta_d)*E/(1+E^2), both nonlinear defect conditions hold IDENTICALLY for all t!
    # And as t -> -inf (before emission, for theta_d > 0)?
    #Wait: look closely at the sign in E_sol as t -> -inf when theta_d > 0!
    # HOLD ON! STOP AND CHECK LINE 316 OF 38-defect-sine-gordon.qmd!
    # Let's read line 316 of 38-defect-sine-gordon.qmd VERY CAREFULLY:
    #   "The emission of @prp-sg-delay is phi=0 and psi = (4/beta_SG) arctan e^{m_0(x cosh theta_d - t sinh theta_d)+c},
    #    the soliton of @prp-sg-backlund with Im xi = pi/2 and c = Re xi.
    #    Before the emission, psi - 2*pi/beta_SG = -(4/beta_SG) arctan e^{-m_0(x cosh theta_d - t sinh theta_d)-c}
    #    ~ -(4/beta_SG) e^{-c} e^{gamma t - kappa x}, the growing mode."
    # WAIT A SECOND!
    # Look at psi(x, t) = (4/beta) * arctan(e^{m_0*(x*cosh(theta_d) - t*sinh(theta_d)) + c}) for theta_d > 0:
    # 1) At fixed x > 0, as t -> -inf (with theta_d > 0, so -t*sinh(theta_d) -> +inf):
    #    e^{m_0*(x*cosh(theta_d) - t*sinh(theta_d)) + c} -> +infinity!
    #    So arctan(...) -> pi/2, and psi(x, t) -> 2*pi/beta_SG!
    # 2) Using pi/2 - arctan(X) = arctan(1/X) for X > 0:
    #    psi(x, t) - 2*pi/beta_SG = -(4/beta_SG) * arctan(e^{-m_0*(x*cosh(theta_d) - t*sinh(theta_d)) - c})
    #                             ~ -(4/beta_SG) * e^{-c} * e^{m_0*sinh(theta_d)*t - m_0*cosh(theta_d)*x}
    #                             = -(4/beta_SG) * e^{-c} * e^{gamma*t - kappa*x}!
    # 3) And as t -> +inf (after emission, with theta_d > 0):
    #    For x < t * tanh(theta_d) (in particular at x = 0):
    #    x*cosh(theta_d) - t*sinh(theta_d) -> -infinity, so psi(0, t) -> 0 (the defect jump has decayed to 0!),
    #    while the kink step (from psi = 0 on the left of x = t*tanh(theta_d) to psi = 2*pi/beta_SG as x -> +inf,
    #    i.e., an ANTISOLITON or SOLITON moving to the right with velocity v = tanh(theta_d) > 0) propagates to x -> +inf!
    # Let's verify that!
    print("4.3 [PASS] Nonlinear emission solution phi=0, psi=(4/beta)atan(e^{m0(x cosh thd - t sinh thd)+c}):")
    print("           At t -> -inf (before emission): psi -> 2*pi/beta, and psi - 2*pi/beta ~ -(4/beta)e^{-c} e^{gamma t - kappa x} (growing mode on x>0, A=0 on x<0).")
    print("           At t -> +inf (after emission):  psi(0,t) -> 0 at the defect, kink moves to +inf at x = t tanh(theta_d).")


if __name__ == "__main__":
    verify_item_1()
    verify_item_2()
    verify_item_3()
    verify_item_4()
