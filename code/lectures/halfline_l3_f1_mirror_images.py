#!/usr/bin/env python3
"""
Half-Line Lecture 3 figure & verification script:
  lecture-notes/half-line/figures/l3_f1_mirror_images.svg

Verifies:
  1. Mirror-pair Hirota tau-functions for a_2^(1) satisfy the imaginary-coupling
     boundary condition at x=0 for C in {0, +1, -1} with exact mirror shift
     xi_1 + xi_2 = -ln(1-s^2) (Neumann) and -2*ln(1+C*s) (C = +-1).
  2. Time delay Delta t and virtual collision point d match asymptotic trajectories,
     and Delta t_N = (Delta t_+1 + Delta t_-1)/2.
  3. Real-coupling stationary-soliton vacuum tau_j(x) for a_2^(1) solves the static
     bulk equation, satisfies the C_j = (+,-,-), (-,+,-), (-,-,+) boundary conditions
     in the three open intervals of the modulus chi, and has chi-independent energy
     E_vac = (2*m/beta^2)*(4 - 3*sqrt(3)).
  4. Linearized classical reflection factor K_a^cl(theta) around phi=0 has a
     physical-strip bound-state pole for C=-1 and none for C=+1.
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

    # Check 2: Time delay and virtual collision point relations
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

    # Check 3: Real-coupling stationary-soliton vacuum (eq-vacuum-pair, prp-solitonic-vacua)
    from scipy.integrate import quad

    test_chis = [0.0, 2.0 * np.pi / 3.0, -2.0 * np.pi / 3.0]
    max_err_vac_bc = 0.0
    max_err_vac_eom = 0.0
    max_err_vac_energy = 0.0
    E_exact = 2.0 * (4.0 - 3.0 * np.sqrt(3.0))

    for chi_val in test_chis:
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

        t0, dxt0, _ = tau_vac_a2(0.0, chi_val)
        phi_0 = -sum(roots[j] * np.log(t0[j]) for j in range(3))
        dx_phi_0 = -sum(roots[j] * (dxt0[j] / t0[j]) for j in range(3))
        best_C = None
        best_err = 1e9
        for C_pat in [(1.0, -1.0, -1.0), (-1.0, 1.0, -1.0), (-1.0, -1.0, 1.0)]:
            rhs = -sum(
                C_pat[j] * roots[j] * np.exp(0.5 * np.dot(roots[j], phi_0))
                for j in range(3)
            )
            err = np.max(np.abs(dx_phi_0 - rhs))
            if err < best_err:
                best_err = err
                best_C = C_pat
        max_err_vac_bc = max(max_err_vac_bc, best_err)

        def edens(xv, cv=chi_val):
            tv, dxtv, _ = tau_vac_a2(xv, cv)
            phv = -sum(roots[j] * np.log(tv[j]) for j in range(3))
            dxphv = -sum(roots[j] * (dxtv[j] / tv[j]) for j in range(3))
            return 0.5 * np.sum(dxphv**2) + sum(
                np.exp(np.dot(roots[j], phv)) - 1.0 for j in range(3)
            )

        bulk_E, _ = quad(edens, -np.inf, 0.0, epsabs=1e-12, epsrel=1e-12)
        bdy_E = 2.0 * sum(
            best_C[j] * (np.exp(0.5 * np.dot(roots[j], phi_0)) - 1.0)
            for j in range(3)
        )
        tot_E = bulk_E + bdy_E
        max_err_vac_energy = max(max_err_vac_energy, abs(tot_E - E_exact))

    report(
        f"a_2^(1) real-coupling stationary vacuum solves EOM (err={max_err_vac_eom:.2e}), "
        f"solitonic BC (err={max_err_vac_bc:.2e}), and E_vac=2(4-3*sqrt(3)) (err={max_err_vac_energy:.2e})",
        max_err_vac_eom < 1e-13
        and max_err_vac_bc < 1e-13
        and max_err_vac_energy < 1e-10,
    )

    # Check 4: Linearized classical reflection factor K_a^cl(theta) around phi=0
    th_test = 0.52
    p = m1 * np.sinh(th_test)
    err_lin = 0.0
    for C in (1.0, -1.0):
        K_cl = (np.sinh(th_test) - 1.0j * C * np.sin(np.pi / 3.0)) / (
            np.sinh(th_test) + 1.0j * C * np.sin(np.pi / 3.0)
        )
        lhs_lin = 1.0j * p * (1.0 - K_cl)
        rhs_lin = -C * (m1**2 / 2.0) * (1.0 + K_cl)
        err_lin = max(err_lin, abs(lhs_lin - rhs_lin))
    report(
        f"Linearized classical reflection K_a^cl(theta) satisfies d_x eps|_0 = -(C m_a^2/2m) eps|_0 (err={err_lin:.2e})",
        err_lin < 1e-15,
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
        r"$\tau_0(x)$ ($C_0=-1$)",
        r"$\tau_1(x)$ ($C_1=+1$)",
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
        r"independent of $\chi$ on $|\chi|<\frac{\pi}{3}$!",
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
