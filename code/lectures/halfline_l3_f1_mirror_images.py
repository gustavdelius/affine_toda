#!/usr/bin/env python3
"""
Half-Line Lecture 3 figure & verification script:
  lecture-notes/half-line/figures/l3_f1_mirror_images.svg

Verifies (m = beta = 1, a_2^(1), species a = 1, sampled parameters):
  1. The mirror pair with shift xi_1 + xi_2 = -ln(1-s^2) (Neumann) or -2 ln(1+Cs) (C = +-1)
     satisfies the imaginary-coupling boundary condition at x = 0, relative to the initial
     boundary vacuum v_a computed from the solution.
  2. The blackboard steps for C = +-1: d_x tau_j(0,t) = k/(1+Cs) (zeta^{ja}U + zeta^{-ja}/U + 2(1-Cs)),
     k s = 2 m sin^2(pi a/h), and (d_x tau_j - m C sqrt(tau_{j+1} tau_{j-1}))/tau_j = k - m C with
     the square root 4 c_{j+1} c_{j-1}/(1+Cs) at integer indices j +- 1.
  3. Time delay and virtual collision point: Delta t_N = (Delta t_+1 + Delta t_-1)/2 and
     d = -(xi_1 + xi_2 + ln A)/(2 m_a cosh theta) for C = +1.
  4. Real-coupling stationary pair (eq-lec3-vacuum-pair): for several moduli chi in each of the
     three intervals of the lecture's table, the field is positive-tau (regular) on x <= 0, solves
     the static field equation, satisfies the boundary condition with exactly the sign pattern
     that the table assigns to that interval (and with none of the other seven patterns), the
     cosine-product formula eq-lec3-cj-from-chi at integer indices gives the same pattern, and
     the energy is E_vac = 2 (4 - 3 sqrt 3).
  5. Linearized classical reflection factor: the rational form satisfies the linearized
     boundary condition, and equals [-(a)(h-a)]^{-C} for C = -1, 0, +1.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib.patches as patches

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "figures"))
import style as st  # noqa: E402
from style import (  # noqa: E402
    plt,
    GREY,
    LIGHTGREY,
    BLUE,
    RED,
    GREEN,
    ORANGE,
    PURPLE,
    setup,
    save,
    report,
    all_ok,
)

OUT_DIR = st.ROOT / "lecture-notes" / "half-line" / "figures"

ZETA = np.exp(2.0j * np.pi / 3.0)


# ---------------------------------------------------------------------------
# Hirota mirror-pair helpers for a_2^(1) (h = 3, m_1 = sqrt(3)*m)
# ---------------------------------------------------------------------------
def tau_mirror_a2(x, t, theta, xi1, C):
    """
    Return tau_j(x,t) and d_x tau_j(x,t) (j=0,1,2) for a species-1 soliton
    (theta_1 = +theta) and its species-2 mirror antisoliton (theta_2 = -theta)
    in a_2^(1) (with m = 1, m_1 = sqrt(3)).
    In the Hirota ansatz, E_r = exp(m_a(x cosh theta_r - t sinh theta_r) + xi_r),
    so BOTH E_1 and E_2 have +m_1*cosh(theta)*x!
    """
    m1 = np.sqrt(3.0)
    s = np.sin(np.pi / 3.0) / np.cosh(theta)
    A = 1.0 - s**2
    if C == 0:
        shift = -np.log(1.0 - s**2)
    else:
        shift = -2.0 * np.log(1.0 + C * s)
    xi2 = shift - xi1

    k = m1 * np.cosh(theta)
    w = m1 * np.sinh(theta)
    E1 = np.exp(k * x - w * t + xi1)
    E2 = np.exp(k * x + w * t + xi2)

    taus = np.zeros((3, np.size(x)), dtype=complex)
    dx_taus = np.zeros((3, np.size(x)), dtype=complex)
    for j in range(3):
        z = ZETA**j
        zb = ZETA ** (-j)
        taus[j] = 1.0 + z * E1 + zb * E2 + A * E1 * E2
        dx_taus[j] = k * (z * E1 + zb * E2 + 2.0 * A * E1 * E2)
    return taus, dx_taus


def tau_vac_a2(x, chi):
    """
    Real-coupling continued stationary mirror pair for a_2^(1) (eq-vacuum-pair):
      tau_j(x) = 1 + 2*cos(chi + 2*pi*j/3)/(1+s)*e^{m_1 x} + (1-s)/(1+s)*e^{2*m_1 x}
    with m = 1, m_1 = sqrt(3), s = sqrt(3)/2.
    """
    m1 = np.sqrt(3.0)
    s = 0.5 * np.sqrt(3.0)
    ex = np.exp(m1 * x)
    ex2 = np.exp(2.0 * m1 * x)
    taus = []
    dx_taus = []
    d2x_taus = []
    for j in range(3):
        c = 2.0 * np.cos(chi + 2.0 * np.pi * j / 3.0) / (1.0 + s)
        d = (1.0 - s) / (1.0 + s)
        tj = 1.0 + c * ex + d * ex2
        dxtj = m1 * (c * ex + 2.0 * d * ex2)
        d2xtj = (m1**2) * (c * ex + 4.0 * d * ex2)
        taus.append(tj)
        dx_taus.append(dxtj)
        d2x_taus.append(d2xtj)
    return np.array(taus), np.array(dx_taus), np.array(d2x_taus)


# ---------------------------------------------------------------------------
# Verification checks
# ---------------------------------------------------------------------------
def run_checks():
    roots = [
        np.array([-1.0, 0.0, 1.0]),
        np.array([1.0, -1.0, 0.0]),
        np.array([0.0, 1.0, -1.0]),
    ]
    theta = 0.65
    xi1 = 0.31 - 0.42j

    # Check 1: Mirror pair satisfies imaginary-coupling boundary condition at x=0
    x_unwrap = np.linspace(-25.0, 0.0, 1001)
    max_err_bc = 0.0
    for C in (0, 1, -1):
        T0, _ = tau_mirror_a2(x_unwrap, -20.0, theta, xi1, C)
        L0 = np.log(np.abs(T0)) + 1.0j * np.unwrap(np.angle(T0), axis=1)
        phi_t0_raw = 1.0j * sum(roots[j] * L0[j, -1] for j in range(3))
        # Round to the exact (2*pi) * coweight lattice vacuum
        coords = np.array(
            [np.real(np.dot(roots[k], phi_t0_raw) / (2.0 * np.pi)) for k in (1, 2)]
        )
        A_mat = np.array([roots[1], roots[2], np.ones(3)])
        v_a = np.linalg.solve(
            A_mat, np.concatenate([2.0 * np.pi * np.round(coords), [0.0]])
        )
        for t in (-1.5, -0.3, 0.4, 1.8):
            T, dT = tau_mirror_a2(x_unwrap, t, theta, xi1, C)
            L = np.log(np.abs(T)) + 1.0j * np.unwrap(np.angle(T), axis=1)
            phi_0 = 1.0j * sum(roots[j] * L[j, -1] for j in range(3)) - v_a
            dx_phi_0 = 1.0j * sum(
                roots[j] * (dT[j, -1] / T[j, -1]) for j in range(3)
            )
            rhs = 1.0j * sum(
                C * roots[j] * np.exp(0.5j * np.dot(roots[j], phi_0))
                for j in range(3)
            )
            max_err_bc = max(max_err_bc, np.max(np.abs(dx_phi_0 - rhs)))
    report(
        f"a_2^(1) mirror pair satisfies boundary condition at x=0 for C=0,+1,-1 (max err={max_err_bc:.2e})",
        max_err_bc < 1e-13,
    )

    # Check 2: blackboard steps of the C = +-1 perfect-square proof
    m1 = np.sqrt(3.0)
    sa = np.sin(np.pi / 3.0)
    max_err_steps = 0.0
    for C in (1, -1):
        s_ = sa / np.cosh(theta)
        k_ = m1 * np.cosh(theta)
        D = 1.0 + C * s_
        max_err_steps = max(max_err_steps, abs(k_ * s_ - 2.0 * sa**2))
        for t in (-0.7, 0.2, 1.1):
            T, dT = tau_mirror_a2(np.array([0.0]), t, theta, xi1, C)
            xi2 = (-2.0 * np.log(1.0 + C * s_)) - xi1
            U = np.exp(-m1 * np.sinh(theta) * t + 0.5 * (xi1 - xi2))
            psi = np.log(U) / 2.0j
            c = lambda kk: np.cos(psi + np.pi * kk / 3.0)
            for j in range(3):
                z = ZETA**j * U
                d_pred = k_ / D * (z + 1.0 / z + 2.0 * (1.0 - C * s_))
                root = 4.0 * c(j + 1) * c(j - 1) / D  # integer indices j +- 1
                ratio = (dT[j, 0] - C * root) / T[j, 0]
                max_err_steps = max(
                    max_err_steps,
                    abs(dT[j, 0] - d_pred),
                    abs(root**2 - T[(j + 1) % 3, 0] * T[(j - 1) % 3, 0]),
                    abs(ratio - (k_ - C)),
                )
    report(
        f"C=+-1 steps: d_x tau_j(0) with 2(1-Cs), k s = 2m sin^2, ratio = k - mC (max err={max_err_steps:.2e})",
        max_err_steps < 1e-10,
    )

    # Check 3: Time delay and virtual collision point relations
    m1 = np.sqrt(3.0)
    s = np.sin(np.pi / 3.0) / np.cosh(theta)
    dt_p1 = 2.0 * np.log(1.0 + s) / (m1 * np.sinh(theta))
    dt_m1 = 2.0 * np.log(1.0 - s) / (m1 * np.sinh(theta))
    dt_N = np.log(1.0 - s**2) / (m1 * np.sinh(theta))
    err_dt_avg = abs(dt_N - 0.5 * (dt_p1 + dt_m1))
    d_p1 = np.log((1.0 + s) / (1.0 - s)) / (2.0 * m1 * np.cosh(theta))
    A = 1.0 - s**2
    shift_p1 = -2.0 * np.log(1.0 + s)
    d_from_shift = -(shift_p1 + np.log(A)) / (2.0 * m1 * np.cosh(theta))
    err_d = abs(d_p1 - d_from_shift)
    report(
        f"Time delay average Delta t_N == (Delta t_+ + Delta t_-)/2 and virtual collision d "
        f"(err_dt={err_dt_avg:.2e}, err_d={err_d:.2e})",
        max(err_dt_avg, err_d) < 1e-15,
    )

    # Check 4: Real-coupling stationary-soliton vacuum against the lecture's interval table
    from scipy.integrate import quad

    # the table printed in the lecture: interval of chi -> (C_0, C_1, C_2)
    table = [
        ((-np.pi, -np.pi / 3.0), (-1.0, 1.0, -1.0)),
        ((-np.pi / 3.0, np.pi / 3.0), (1.0, -1.0, -1.0)),
        ((np.pi / 3.0, np.pi), (-1.0, -1.0, 1.0)),
    ]
    all_patterns = [
        tuple(float(c) for c in pat)
        for pat in np.array(np.meshgrid([1, -1], [1, -1], [1, -1])).T.reshape(-1, 3)
    ]
    E_exact = 2.0 * (4.0 - 3.0 * np.sqrt(3.0))
    max_err_vac_eom = 0.0
    max_err_vac_bc = 0.0
    min_err_other = np.inf
    max_err_vac_energy = 0.0
    formula_ok = True
    min_tau = np.inf
    x_grid = np.linspace(-12.0, 0.0, 601)
    for (lo, hi), C_pat in table:
        for frac in (0.15, 0.5, 0.85):
            chi_val = lo + frac * (hi - lo)
            t_x, dxt_x, d2xt_x = tau_vac_a2(-0.7, chi_val)
            phi_x = -sum(roots[j] * np.log(t_x[j]) for j in range(3))
            d2x_phi = -sum(
                roots[j] * (d2xt_x[j] / t_x[j] - (dxt_x[j] / t_x[j]) ** 2)
                for j in range(3)
            )
            grad_V = sum(
                roots[j] * np.exp(np.dot(roots[j], phi_x)) for j in range(3)
            )
            max_err_vac_eom = max(max_err_vac_eom, np.max(np.abs(d2x_phi - grad_V)))
            tg, _, _ = tau_vac_a2(x_grid, chi_val)
            min_tau = min(min_tau, np.min(tg))

            t0, dxt0, _ = tau_vac_a2(0.0, chi_val)
            phi_0 = -sum(roots[j] * np.log(t0[j]) for j in range(3))
            dx_phi_0 = -sum(roots[j] * (dxt0[j] / t0[j]) for j in range(3))

            def bc_err(pat):
                rhs = -sum(
                    pat[j] * roots[j] * np.exp(0.5 * np.dot(roots[j], phi_0))
                    for j in range(3)
                )
                return np.max(np.abs(dx_phi_0 - rhs))

            max_err_vac_bc = max(max_err_vac_bc, bc_err(C_pat))
            for pat in all_patterns:
                if pat != C_pat:
                    min_err_other = min(min_err_other, bc_err(pat))

            # eq-lec3-cj-from-chi with integer (unreduced) indices j +- 1
            cc = lambda kk: np.cos(0.5 * chi_val + np.pi * kk / 3.0)
            C_formula = tuple(float(np.sign(cc(j + 1) * cc(j - 1))) for j in range(3))
            formula_ok = formula_ok and (C_formula == C_pat)

            def edens(xv, cv=chi_val):
                tv, dxtv, _ = tau_vac_a2(xv, cv)
                phv = -sum(roots[j] * np.log(tv[j]) for j in range(3))
                dxphv = -sum(roots[j] * (dxtv[j] / tv[j]) for j in range(3))
                return 0.5 * np.sum(dxphv**2) + sum(
                    np.exp(np.dot(roots[j], phv)) - 1.0 for j in range(3)
                )

            bulk_E, _ = quad(edens, -np.inf, 0.0, epsabs=1e-12, epsrel=1e-12, limit=200)
            bdy_E = 2.0 * sum(
                C_pat[j] * (np.exp(0.5 * np.dot(roots[j], phi_0)) - 1.0)
                for j in range(3)
            )
            max_err_vac_energy = max(max_err_vac_energy, abs(bulk_E + bdy_E - E_exact))

    report(
        f"a_2^(1) stationary vacuum: tau_j > 0 on x<=0 (min={min_tau:.3f}), static EOM (err={max_err_vac_eom:.2e})",
        min_tau > 0.0 and max_err_vac_eom < 1e-13,
    )
    report(
        f"a_2^(1) stationary vacuum satisfies the BC with the table's pattern in each chi interval "
        f"(err={max_err_vac_bc:.2e}), all other patterns fail (min err={min_err_other:.3f}), "
        f"and eq-lec3-cj-from-chi at integer indices gives the table ({formula_ok})",
        max_err_vac_bc < 1e-12 and min_err_other > 1e-2 and formula_ok,
    )
    report(
        f"a_2^(1) stationary vacuum energy E_vac = 2(4-3*sqrt(3)) in all three intervals (err={max_err_vac_energy:.2e})",
        max_err_vac_energy < 1e-10,
    )

    # Check 5: Linearized classical reflection factor K_a^cl(theta) around phi=0
    def block(xb, th, h=3):
        return np.sinh(th / 2 + 1j * np.pi * xb / (2 * h)) / np.sinh(
            th / 2 - 1j * np.pi * xb / (2 * h)
        )

    err_lin = 0.0
    err_blk = 0.0
    for th_test in (0.52, 0.31 + 0.27j):
        p = m1 * np.sinh(th_test)
        for C in (1.0, 0.0, -1.0):
            K_cl = (np.sinh(th_test) - 1.0j * C * np.sin(np.pi / 3.0)) / (
                np.sinh(th_test) + 1.0j * C * np.sin(np.pi / 3.0)
            )
            lhs_lin = 1.0j * p * (1.0 - K_cl)
            rhs_lin = -C * (m1**2 / 2.0) * (1.0 + K_cl)
            err_lin = max(err_lin, abs(lhs_lin - rhs_lin))
            K_blk = (-block(1, th_test) * block(2, th_test)) ** (-C)
            err_blk = max(err_blk, abs(K_cl - K_blk))
    report(
        f"Linearized K_a^cl solves d_x eps|_0 = -(C m_a^2/2m) eps|_0 (err={err_lin:.2e}) and equals "
        f"[-(a)(h-a)]^(-C) for C=+1,0,-1 (err={err_blk:.2e})",
        err_lin < 1e-14 and err_blk < 1e-14,
    )


def arrow_mid(ax, x, y, frac, color):
    idx = int(len(x) * frac)
    ax.annotate(
        "",
        xy=(x[idx + 1], y[idx + 1]),
        xytext=(x[idx - 1], y[idx - 1]),
        arrowprops=dict(
            arrowstyle="->,head_width=0.18,head_length=0.22",
            color=color,
            lw=1.4,
        ),
    )


# ---------------------------------------------------------------------------
# Figure generation
# ---------------------------------------------------------------------------
def make_figure():
    setup()
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(7.4, 3.2), gridspec_kw={"width_ratios": [1.05, 1.15]}
    )

    # -----------------------------------------------------------------------
    # Panel (a): Mirror-image collision and time delay for C = +1 (a_2^(1))
    # -----------------------------------------------------------------------
    theta = 0.6
    m1 = np.sqrt(3.0)
    v = np.tanh(theta)
    s = np.sin(np.pi / 3.0) / np.cosh(theta)
    dt_p1 = 2.0 * np.log(1.0 + s) / (m1 * np.sinh(theta))
    d_p1 = np.log((1.0 + s) / (1.0 - s)) / (2.0 * m1 * np.cosh(theta))

    ax1.add_patch(
        patches.Rectangle(
            (0.0, -2.2),
            1.1,
            4.4,
            facecolor=LIGHTGREY,
            edgecolor="none",
            alpha=0.45,
        )
    )
    ax1.plot([0.0, 0.0], [-2.2, 2.2], color=PURPLE, lw=2.2, zorder=4)
    ax1.text(
        0.52,
        1.88,
        r"virtual region $x>0$",
        fontsize=7.2,
        color=GREY,
        ha="center",
    )
    ax1.text(
        -0.06,
        1.92,
        r"wall $x=0$",
        fontsize=7.5,
        color=PURPLE,
        ha="right",
    )

    t_coll = d_p1 / v
    t_in_phys = np.linspace(-2.0, 0.0, 60)
    x_in_phys = v * t_in_phys
    t_in_virt = np.linspace(0.0, t_coll, 40)
    x_in_virt = v * t_in_virt
    ax1.plot(
        x_in_phys, t_in_phys, color=BLUE, lw=1.9, label=r"Incoming soliton $a$"
    )
    arrow_mid(ax1, x_in_phys, t_in_phys, 0.5, BLUE)
    ax1.plot(x_in_virt, t_in_virt, color=BLUE, lw=1.2, ls="--")

    t_mir = np.linspace(-2.0, t_coll, 60)
    x_mir = d_p1 - v * (t_mir - t_coll)
    ax1.plot(
        x_mir,
        t_mir,
        color=GREY,
        lw=1.3,
        ls="--",
        label=r"Mirror image $\bar a$ ($x>0$)",
    )

    ax1.plot(d_p1, t_coll, "o", color=GREEN, ms=5.5, zorder=6)
    ax1.annotate(
        r"virtual collision $x=d>0$",
        xy=(d_p1, t_coll),
        xytext=(0.18, -0.55),
        fontsize=7.0,
        color=GREEN,
        arrowprops=dict(arrowstyle="->", color=GREEN, lw=0.8),
    )

    t_out_phys = np.linspace(dt_p1, 2.1, 60)
    x_out_phys = -v * (t_out_phys - dt_p1)
    ax1.plot(
        x_out_phys,
        t_out_phys,
        color=RED,
        lw=1.9,
        label=r"Reflected antisoliton $\bar a$",
    )
    arrow_mid(ax1, x_out_phys, t_out_phys, 0.55, RED)

    t_ref = np.linspace(0.0, 1.6, 40)
    ax1.plot(-v * t_ref, t_ref, color=GREY, lw=1.0, ls=":")

    ax1.annotate(
        "",
        xy=(-0.08, dt_p1),
        xytext=(-0.08, 0.0),
        arrowprops=dict(arrowstyle="<->", color=ORANGE, lw=1.2),
    )
    ax1.text(
        -0.16,
        0.5 * dt_p1,
        r"$\Delta t > 0$",
        fontsize=7.8,
        color=ORANGE,
        ha="right",
        va="center",
    )

    ax1.set_xlim(-1.35, 1.05)
    ax1.set_ylim(-2.1, 2.15)
    ax1.set_xlabel(r"Position $m x$")
    ax1.set_ylabel(r"Time $m t$")
    ax1.set_title(
        r"(a) Mirror collision & time delay ($C=+1$, $a_2^{(1)}$)",
        fontsize=8.8,
        pad=6,
    )
    ax1.legend(loc="lower left", fontsize=6.7)

    # -----------------------------------------------------------------------
    # Panel (b): Real-coupling stationary-soliton vacuum profiles tau_j(x)
    # -----------------------------------------------------------------------
    x_plot = np.linspace(-3.2, 0.0, 300)
    chi_demo = 0.25
    taus_x, _, _ = tau_vac_a2(x_plot, chi_demo)
    colors_tau = [BLUE, GREEN, ORANGE]
    labels_tau = [
        r"$\tau_0(x)$ ($C_0=+1$)",
        r"$\tau_1(x)$ ($C_1=-1$)",
        r"$\tau_2(x)$ ($C_2=-1$)",
    ]
    for j in range(3):
        ax2.plot(
            x_plot, taus_x[j], color=colors_tau[j], lw=1.8, label=labels_tau[j]
        )

    ax2.axvline(0.0, color=PURPLE, lw=2.0)
    ax2.axhline(1.0, color=GREY, lw=0.8, ls=":")
    ax2.text(
        -2.95,
        1.07,
        r"bulk vacuum $\tau_j\to 1$ ($\vec\phi\to 0$)",
        fontsize=7.2,
        color=GREY,
    )
    ax2.text(
        -1.65,
        0.22,
        r"Modulus $\chi$: $E_{\rm vac} = \frac{2m}{\beta^2}(4-3\sqrt{3})$"
        "\n"
        r"independent of $\chi$ on $|\chi|<\frac{\pi}{3}$",
        fontsize=7.3,
        color=GREY,
        ha="center",
        bbox=dict(
            boxstyle="round,pad=0.35",
            facecolor="none",
            edgecolor=GREY,
            lw=0.8,
        ),
    )

    ax2.set_xlim(-3.2, 0.12)
    ax2.set_ylim(0.0, 2.15)
    ax2.set_xlabel(r"Position $m x \leq 0$")
    ax2.set_ylabel(r"Hirota vacuum factors $\tau_j(x)$")
    ax2.set_title(
        r"(b) Real-coupling solitonic vacuum ($\chi=0.25$, $a_2^{(1)}$)",
        fontsize=8.8,
        pad=6,
    )
    ax2.legend(loc="upper left", fontsize=7.0)

    fig.tight_layout()
    save(fig, "l3_f1_mirror_images", OUT_DIR)


if __name__ == "__main__":
    run_checks()
    make_figure()
    sys.exit(0 if all_ok() else 1)
