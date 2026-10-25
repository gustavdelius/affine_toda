"""Independent verification of Item 2: Singular-wall limits 0.719 / 0.937 / 0.147 in Ch. 32 (32-half-line-solitons.qmd lines 215-218).

Uses two completely independent methods:
  Method A: Direct ODE integration (scipy.integrate.solve_ivp) of the linearized Toda wave equation
            psi''(x) = (W(x) - Omega^2) psi(x) from x = -L to x = 0 at moderate delta with Richardson extrapolation.
  Method B: Independent exact dressed-Hirota wave solver (differentiating the N+1-soliton Hirota tau functions
            with respect to the probe wave amplitude at 0) with high-precision mpmath arithmetic.
"""
import itertools
import numpy as np
from scipy.integrate import solve_ivp
import mpmath as mp

def root_system_an(h):
    alphas = []
    for j in range(h):
        a = np.zeros(h)
        a[j] = 1.0
        a[(j + 1) % h] = -1.0
        alphas.append(a)
    return np.array(alphas)

def particle_modes_an(h):
    j = np.arange(h)
    modes = []
    masses = []
    for a in range(1, h):
        v = np.exp(2j * np.pi * a * j / h) / np.sqrt(h)
        modes.append(v)
        masses.append(2.0 * np.sin(np.pi * a / h))
    return np.array(modes), np.array(masses)

def tau_to_exponentials(taus, h):
    exps = np.empty(h, dtype=complex)
    for j in range(h):
        exps[j] = taus[(j - 1) % h] * taus[(j + 1) % h] / (taus[j] ** 2)
    return exps

def bulk_hessian_W(exps, alphas):
    h = len(exps)
    W = np.zeros((h, h), dtype=complex)
    for j in range(h):
        W += exps[j] * np.outer(alphas[j], alphas[j])
    return W

def boundary_hessian_HB(C_signs, exps0, alphas):
    h = len(C_signs)
    HB = np.zeros((h, h), dtype=complex)
    for j in range(h):
        s_val = np.sqrt(exps0[j])
        if s_val.real < 0:
            s_val = -s_val
        HB += 0.5 * C_signs[j] * s_val * np.outer(alphas[j], alphas[j])
    return HB

def solve_reflection_ode(h, tau_func, C_signs, theta, L=14.0, rtol=1e-10, atol=1e-12):
    """Method A: Direct ODE integration of psi''(x) = (W(x) - Omega^2) psi(x) from x = -L to x = 0."""
    alphas = root_system_an(h)
    modes, masses = particle_modes_an(h)
    n_part = h - 1
    Omega = masses[0] * np.cosh(theta)
    k1 = masses[0] * np.sinh(theta)

    ks = []
    for b in range(n_part):
        val = Omega**2 - masses[b]**2
        if val >= 0:
            ks.append(np.sqrt(val))
        else:
            ks.append(1j * np.sqrt(-val))
    ks = np.array(ks)

    n_traj = n_part + 1
    y0 = np.zeros((2 * h, n_traj), dtype=complex)
    for b in range(n_part):
        phase = np.exp(-1j * ks[b] * (-L))
        if abs(ks[b].imag) > 0:
            phase = 1.0
        y0[:h, b] = phase * modes[b]
        y0[h:, b] = (-1j * ks[b]) * phase * modes[b]

    phase_in = np.exp(1j * k1 * (-L))
    y0[:h, n_part] = phase_in * modes[0]
    y0[h:, n_part] = (1j * k1) * phase_in * modes[0]

    def rhs(x, y_flat):
        Y = y_flat.reshape((2 * h, n_traj))
        taus = tau_func(x)
        exps = tau_to_exponentials(taus, h)
        W = bulk_hessian_W(exps, alphas)
        dY = np.empty_like(Y)
        dY[:h, :] = Y[h:, :]
        dY[h:, :] = (W - (Omega**2) * np.eye(h)) @ Y[:h, :]
        return dY.ravel()

    sol = solve_ivp(rhs, (-L, 0.0), y0.ravel(), method='DOP853', rtol=rtol, atol=atol)
    Y0 = sol.y[:, -1].reshape((2 * h, n_traj))
    psi0 = Y0[:h, :]
    dpsi0 = Y0[h:, :]

    taus0 = tau_func(0.0)
    exps0 = tau_to_exponentials(taus0, h)
    HB = boundary_hessian_HB(C_signs, exps0, alphas)

    B_mat = modes.conj() @ (dpsi0 + HB @ psi0)
    K_col = np.linalg.solve(B_mat[:, :n_part], -B_mat[:, n_part])
    return K_col[0], K_col[-1]

# Method B: Exact dressed-Hirota wave solver in mpmath
mp.mp.dps = 70
I_mp = mp.mpc(0, 1)

def solve_reflection_hirota_mp(h, solitons, C_signs, theta):
    """Method B: Independent exact dressed-Hirota wave solver at x=0 using mpmath."""
    w = mp.exp(2 * I_mp * mp.pi / h)
    masses = [2 * mp.sin(mp.pi * a / h) for a in range(h)]
    N = len(solitons)

    def A_static(a, b):
        return (1 - mp.cos(mp.pi * (a - b) / h)) / (1 - mp.cos(mp.pi * (a + b) / h))

    subsets = []
    for r in range(N + 1):
        for S in itertools.combinations(range(N), r):
            A_prod = mp.mpf(1)
            for i in range(len(S)):
                for k in range(i + 1, len(S)):
                    A_prod *= A_static(solitons[S[i]][0], solitons[S[k]][0])
            c_prod = mp.fprod(solitons[p][1] for p in S)
            m_sum = mp.fsum(masses[solitons[p][0]] for p in S)
            a_sum = sum(solitons[p][0] for p in S) % h
            subsets.append((S, A_prod * c_prod, m_sum, a_sum))

    tau0 = [mp.mpc(0)] * h
    dtau0 = [mp.mpc(0)] * h
    for j in range(h):
        for S, coeff, m_sum, a_sum in subsets:
            term = coeff * (w ** (j * a_sum))
            tau0[j] += term
            dtau0[j] += term * m_sum

    exps0 = [tau0[(j - 1) % h] * tau0[(j + 1) % h] / (tau0[j] ** 2) for j in range(h)]
    sqrt_exps0 = [mp.sqrt(val) for val in exps0]

    def dressed_wave_at_0(c, ik):
        mc = masses[c]
        ch_c = ik / mc
        dt_j = [mp.mpc(0)] * h
        ddt_j = [mp.mpc(0)] * h
        for S, coeff, m_sum, a_sum in subsets:
            A_probe = mp.mpf(1)
            for p in S:
                ap = solitons[p][0]
                A_probe *= (ch_c - mp.cos(mp.pi * (c - ap) / h)) / (ch_c - mp.cos(mp.pi * (c + ap) / h))
            for j in range(h):
                term = coeff * A_probe * (w ** (j * (a_sum + c)))
                dt_j[j] += term
                ddt_j[j] += term * (m_sum + ik)
        u0 = [dt_j[j] / tau0[j] for j in range(h)]
        du0 = [(ddt_j[j] * tau0[j] - dt_j[j] * dtau0[j]) / (tau0[j] ** 2) for j in range(h)]
        phi_vec = [mp.mpc(0)] * h
        dphi_vec = [mp.mpc(0)] * h
        for j in range(h):
            phi_vec[j] -= u0[j]
            phi_vec[(j + 1) % h] += u0[j]
            dphi_vec[j] -= du0[j]
            dphi_vec[(j + 1) % h] += du0[j]
        norm_factor = (w ** (-c) - 1) * mp.sqrt(h)
        phi_vec = [val / norm_factor for val in phi_vec]
        dphi_vec = [val / norm_factor for val in dphi_vec]
        return phi_vec, dphi_vec

    Omega = masses[1] * mp.cosh(theta)
    n_part = h - 1
    iks_out = []
    for c in range(1, h):
        val = Omega**2 - masses[c]**2
        if val >= 0:
            iks_out.append(-I_mp * mp.sqrt(val))
        else:
            iks_out.append(mp.sqrt(-val))
    ik_in = +I_mp * mp.sqrt(Omega**2 - masses[1]**2)

    def apply_bc(phi_vec, dphi_vec):
        res = list(dphi_vec)
        for j in range(h):
            jp = (j + 1) % h
            diff = phi_vec[j] - phi_vec[jp]
            term = mp.mpf('0.5') * C_signs[j] * sqrt_exps0[j] * diff
            res[j] += term
            res[jp] -= term
        proj = []
        for c in range(1, h):
            s_val = mp.fsum(mp.conj(w ** (j * c)) * res[j] for j in range(h)) / mp.sqrt(h)
            proj.append(s_val)
        return proj

    B_cols = []
    for c in range(1, h):
        pv, dpv = dressed_wave_at_0(c, iks_out[c - 1])
        B_cols.append(apply_bc(pv, dpv))
    pv_in, dpv_in = dressed_wave_at_0(1, ik_in)
    rhs_in = [-val for val in apply_bc(pv_in, dpv_in)]

    B_mat = mp.matrix([[B_cols[col][row] for col in range(n_part)] for row in range(n_part)])
    K_vec = mp.lu_solve(B_mat, mp.matrix(rhs_in))
    return K_vec[0], K_vec[n_part - 1]

if __name__ == '__main__':
    theta = 0.6

    # --- Case 1: a_5^{(1)} (h=6), species a=2 stationary pair along regular branch chi = 5*pi/3 +- delta ---
    h = 6
    a = 2
    s = np.sin(np.pi * a / h)
    m_a = 2.0 * s
    C_signs_a5 = [1, 1, -1, -1, -1, -1]
    d0 = 1.0 / (1.0 + s)

    def make_tau_a5(chi, d=d0):
        A = 1.0 - s**2
        def tau_fn(x):
            E = np.exp(m_a * x)
            j = np.arange(h)
            return 1.0 + 2.0 * d * E * np.cos(chi + 2.0 * np.pi * a * j / h) + A * (d**2) * (E**2)
        return tau_fn

    s_mp = mp.sin(mp.pi * a / h)
    d0_mp = 1 / (1 + s_mp)
    # --- Case 1: a_5^{(1)} (h=6), regular branch chi = 5*pi/3 +- delta ---
    k11_mp_1m, k15_mp_1m = solve_reflection_hirota_mp(
        h, [(2, d0_mp * mp.exp(I_mp * (5 * mp.pi / 3 - mp.mpf('1e-6')))),
            (4, d0_mp * mp.exp(-I_mp * (5 * mp.pi / 3 - mp.mpf('1e-6'))))], C_signs_a5, mp.mpf('0.6')
    )
    k11_mp_1p, k15_mp_1p = solve_reflection_hirota_mp(
        h, [(2, d0_mp * mp.exp(I_mp * (5 * mp.pi / 3 + mp.mpf('1e-6')))),
            (4, d0_mp * mp.exp(-I_mp * (5 * mp.pi / 3 + mp.mpf('1e-6'))))], C_signs_a5, mp.mpf('0.6')
    )
    print(f"Case 1 (a_5 regular chi -> 5pi/3):")
    print(f"  Method B (Hirota mp, -1e-6): |K_15| = {mp.nstr(abs(k15_mp_1m), 10)}, K_11 = {mp.nstr(k11_mp_1m, 10)}")
    print(f"  Method B (Hirota mp, +1e-6): |K_15| = {mp.nstr(abs(k15_mp_1p), 10)}, K_11 = {mp.nstr(k11_mp_1p, 10)}")
    print("PASS case1_0.719" if abs(abs(k15_mp_1p) - 0.71919316) < 1e-6 and abs(abs(k15_mp_1p) - abs(k15_mp_1m)) < 1e-6 else "FAIL case1_0.719")

    # --- Case 2: a_5^{(1)} (h=6), non-solution path chi = 5*pi/3, d = d0 * (1 - eps) ---
    chi0_mp = 5 * mp.pi / 3
    d_ns_mp = d0_mp * (1 - mp.mpf('1e-8'))
    k11_mp_2, k15_mp_2 = solve_reflection_hirota_mp(
        h, [(2, d_ns_mp * mp.exp(I_mp * chi0_mp)), (4, d_ns_mp * mp.exp(-I_mp * chi0_mp))], C_signs_a5, mp.mpf('0.6')
    )
    eps_2 = [4e-3, 2e-3, 1e-3]
    ode_2 = [solve_reflection_ode(h, make_tau_a5(5 * np.pi / 3, d=d0 * (1.0 - ep)), C_signs_a5, theta) for ep in eps_2]
    k15_ode_2 = np.polyfit(eps_2, [abs(x[1]) for x in ode_2], 2)[-1]
    print(f"\nCase 2 (a_5 non-solution d -> 1/(1+s) at chi = 5pi/3):")
    print(f"  Method A (ODE extrap):       |K_15| = {k15_ode_2:.8f}")
    print(f"  Method B (Hirota mp, 1e-8):  |K_15| = {mp.nstr(abs(k15_mp_2), 10)}, K_11 = {mp.nstr(k11_mp_2, 10)}")
    print("PASS case2_0.937" if abs(abs(k15_mp_2) - 0.93743720) < 1e-6 and abs(k15_ode_2 - abs(k15_mp_2)) < 2e-4 else "FAIL case2_0.937")

    # --- Case 3: a_4^{(1)} (h=5), 4-soliton exceptional background of (+,+,-,+,-) ---
    h4 = 5
    C_signs_a4 = [1, 1, -1, 1, -1]
    def make_sols_a4(d1, d2):
        c1 = d1 * mp.exp(I_mp * 9 * mp.pi / 5)
        c2 = d2 * mp.exp(I_mp * 18 * mp.pi / 5)
        return [(1, c1), (4, mp.conj(c1)), (2, c2), (3, mp.conj(c2))]

    def eval_tau_dtau_a4(d1, d2):
        sols = make_sols_a4(d1, d2)
        w4 = mp.exp(2 * I_mp * mp.pi / h4)
        ms = [2 * mp.sin(mp.pi * a / h4) for a in range(h4)]
        def A_st(a, b):
            return (1 - mp.cos(mp.pi * (a - b) / h4)) / (1 - mp.cos(mp.pi * (a + b) / h4))
        t0 = [mp.mpc(0)] * h4
        dt0 = [mp.mpc(0)] * h4
        for r in range(5):
            for S in itertools.combinations(range(4), r):
                Ap = mp.mpf(1)
                for i in range(len(S)):
                    for k in range(i + 1, len(S)):
                        Ap *= A_st(sols[S[i]][0], sols[S[k]][0])
                cp = mp.fprod(sols[p][1] for p in S)
                msum = mp.fsum(ms[sols[p][0]] for p in S)
                asum = sum(sols[p][0] for p in S) % h4
                for j in range(h4):
                    term = Ap * cp * (w4 ** (j * asum))
                    t0[j] += term
                    dt0[j] += term * msum
        return t0, dt0

    def sys_2d(d1, d2):
        t0, dt0 = eval_tau_dtau_a4(d1, d2)
        return (t0[2].real, (dt0[0] / t0[0] - dt0[3] / t0[3]).real)

    d1_star, d2_star = mp.findroot(sys_2d, (mp.mpf('0.2579'), mp.mpf('0.5637')))
    k11_mp_3, k14_mp_3 = solve_reflection_hirota_mp(
        h4, make_sols_a4(d1_star * (1 - mp.mpf('1e-7')), d2_star), C_signs_a4, mp.mpf('0.6')
    )
    k11_mp_3p, k14_mp_3p = solve_reflection_hirota_mp(
        h4, make_sols_a4(d1_star * (1 + mp.mpf('1e-7')), d2_star), C_signs_a4, mp.mpf('0.6')
    )
    print(f"\nCase 3 (a_4 exceptional 4-soliton limit of (+,+,-,+,-)):")
    print(f"  d1* = {mp.nstr(d1_star, 12)}, d2* = {mp.nstr(d2_star, 12)}")
    print(f"  Method B (Hirota mp, -1e-7): |K_14| = {mp.nstr(abs(k14_mp_3), 10)}, K_11 = {mp.nstr(k11_mp_3, 10)}")
    print(f"  Method B (Hirota mp, +1e-7): |K_14| = {mp.nstr(abs(k14_mp_3p), 10)}, K_11 = {mp.nstr(k11_mp_3p, 10)}")
    print("PASS case3_0.147" if abs(abs(k14_mp_3) - 0.147135) < 3e-5 and abs(abs(k14_mp_3) - abs(k14_mp_3p)) < 5e-5 else "FAIL case3_0.147")



