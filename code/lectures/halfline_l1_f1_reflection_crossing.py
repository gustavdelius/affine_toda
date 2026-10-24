#!/usr/bin/env python3
"""
Half-Line Lecture 1 figure & verification script:
  lecture-notes/half-line/figures/l1_f1_reflection_crossing.svg

Verifies:
  1. Equivalence of S/K, braid (R_check), and twisted reflection equations.
  2. Equivalence of the three diagonal boundary crossing-unitarity equations
     (book form, Corrigan-Delius, and Delius-Gandenberger).
  3. Boundary scaling Lee-Yang reflection factor K_b(theta):
     - boundary unitarity K_b(theta) K_b(-theta) = 1
     - boundary crossing-unitarity K_b(i*pi/2 - theta) = S_LY(2*theta) K_b(i*pi/2 + theta)
     - quadratic boundary bootstrap K_b(theta) = S_LY(2*theta) K_b(theta+i*pi/3) K_b(theta-i*pi/3)
       (and why the overall + sign is fixed)
     - residues of K_{(1)}(theta) at theta = i*pi/2 ((i/2)*g^2 with g = 2*i*sqrt(2*sqrt(3)-3))
       and at theta = i*pi/6 (-(i/2)*Gamma*g)
     - boundary bound-state shift K_b^(1)(theta) = K_{b+2}(theta)
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


# ---------------------------------------------------------------------------
# Scaling Lee-Yang building blocks (h = 6 => period 2h = 12)
# ---------------------------------------------------------------------------
def blk6(theta, x):
    """Block (x) = sinh(theta/2 + i*pi*x/12) / sinh(theta/2 - i*pi*x/12)."""
    return np.sinh(0.5 * theta + 1.0j * np.pi * x / 12.0) / np.sinh(
        0.5 * theta - 1.0j * np.pi * x / 12.0
    )


def S_LY(theta):
    """Bulk scaling Lee-Yang S-matrix: f_{2/3}(theta) = -(2)(4)."""
    return -blk6(theta, 2.0) * blk6(theta, 4.0)


def K_LY(theta, b):
    """Dorey-Pocklington-Tateo-Watts boundary Lee-Yang reflection factor K_b(theta)."""
    base = blk6(theta, 1.0) * blk6(theta, 3.0) / blk6(theta, 4.0)
    b_dep = (
        blk6(theta, 1.0 + b)
        * blk6(theta, 5.0 - b)
        / (blk6(theta, 1.0 - b) * blk6(theta, 5.0 + b))
    )
    return base * b_dep


# ---------------------------------------------------------------------------
# Verification checks
# ---------------------------------------------------------------------------
def run_checks():
    # Check 1: Equivalence of S/K and braid reflection equations on C^2 x C^2
    eta = 0.73 + 0.19j
    P = np.array(
        [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]],
        dtype=complex,
    )

    def R_6v(u):
        a = np.sinh(u + eta)
        b = np.sinh(u)
        c = np.sinh(eta)
        return np.array(
            [[a, 0, 0, 0], [0, b, c, 0], [0, c, b, 0], [0, 0, 0, a]],
            dtype=complex,
        )

    def R_check(u):
        return P @ R_6v(u)

    xi_B = 0.41 - 0.12j

    def K_diag(u):
        return np.diag([np.sinh(xi_B + u), np.sinh(xi_B - u)])

    u1, u2 = 0.37 + 0.11j, -0.29 + 0.23j
    I2 = np.eye(2, dtype=complex)

    R12_diff = R_6v(u1 - u2)
    R21_sum = P @ R_6v(u1 + u2) @ P
    R12_sum = R_6v(u1 + u2)
    R21_diff = P @ R_6v(u1 - u2) @ P
    K1 = np.kron(K_diag(u1), I2)
    K2 = np.kron(I2, K_diag(u2))

    lhs_sk = R12_diff @ K1 @ R21_sum @ K2
    rhs_sk = K2 @ R12_sum @ K1 @ R21_diff
    err_sk = np.max(np.abs(lhs_sk - rhs_sk))

    K1_u2 = np.kron(K_diag(u2), I2)
    lhs_braid = R_check(u1 - u2) @ K1 @ R_check(u1 + u2) @ K1_u2
    rhs_braid = K1_u2 @ R_check(u1 + u2) @ K1 @ R_check(u1 - u2)
    err_braid = np.max(np.abs(lhs_braid - rhs_braid))

    err_rel = np.max(np.abs(P @ (lhs_sk - rhs_sk) - (lhs_braid - rhs_braid)))
    report(
        f"S/K and braid reflection equations match (max err={max(err_sk, err_braid, err_rel):.2e})",
        max(err_sk, err_braid, err_rel) < 1e-13,
    )

    # Check 2: Equivalence of the three diagonal boundary crossing-unitarity equations
    th = 0.43 - 0.17j
    b_val = -0.37
    ip2 = 0.5j * np.pi
    err_book = abs(
        K_LY(ip2 - th, b_val) - S_LY(2.0 * th) * K_LY(ip2 + th, b_val)
    )
    err_cd = abs(
        K_LY(th + ip2, b_val) * K_LY(th - ip2, b_val) * S_LY(2.0 * th) - 1.0
    )
    err_dg = abs(
        K_LY(th, b_val) * K_LY(th + 1.0j * np.pi, b_val) - S_LY(2.0 * th)
    )
    report(
        f"Three diagonal boundary crossing-unitarity forms (max err={max(err_book, err_cd, err_dg):.2e})",
        max(err_book, err_cd, err_dg) < 1e-13,
    )

    # Check 3: Boundary Lee-Yang unitarity and quadratic bootstrap (shift i*pi/3)
    ip3 = 1.0j * np.pi / 3.0
    err_unit = abs(K_LY(th, b_val) * K_LY(-th, b_val) - 1.0)
    err_boot_plus = abs(
        K_LY(th, b_val)
        - K_LY(th + ip3, b_val) * K_LY(th - ip3, b_val) * S_LY(2.0 * th)
    )
    K_minus = -K_LY(th, b_val)
    err_boot_minus = abs(
        K_minus
        - (-K_LY(th + ip3, b_val)) * (-K_LY(th - ip3, b_val)) * S_LY(2.0 * th)
    )
    report(
        f"Boundary Lee-Yang unitarity & quadratic bootstrap sign fixing (err+={err_boot_plus:.2e})",
        err_unit < 1e-13
        and err_boot_plus < 1e-13
        and abs(err_boot_minus - 2.0 * abs(K_LY(th, b_val))) < 1e-12,
    )

    # Check 4: Residues of K_{(1)}(theta) = K_0(theta) at i*pi/2 and i*pi/6
    eps = 1e-8
    res_pi2 = 0.5 * (K_LY(ip2 + eps, 0.0) - K_LY(ip2 - eps, 0.0)) * eps
    res_pi6 = (
        0.5
        * (
            K_LY(1.0j * np.pi / 6.0 + eps, 0.0)
            - K_LY(1.0j * np.pi / 6.0 - eps, 0.0)
        )
        * eps
    )
    res_bulk = (
        0.5
        * (
            S_LY(2.0j * np.pi / 3.0 + eps)
            - S_LY(2.0j * np.pi / 3.0 - eps)
        )
        * eps
    )
    Gamma_sq = res_bulk / 1.0j  # S ~ i * Gamma^2 / (theta - 2*i*pi/3)
    g_sq_exact = -4.0 * (2.0 * np.sqrt(3.0) - 3.0)
    err_pi2 = abs(res_pi2 - 0.5j * g_sq_exact)
    err_pi6 = abs(res_pi6**2 + 0.25 * Gamma_sq * g_sq_exact)
    report(
        f"Lee-Yang residues at i*pi/2 and i*pi/6 (err_pi2={err_pi2:.2e}, err_pi6={err_pi6:.2e})",
        max(err_pi2, err_pi6) < 1e-9,
    )

    # Check 5: Excited boundary reflection factor K_b^(1)(theta) = K_{b+2}(theta)
    u_pole = 1.0j * np.pi * (b_val + 1.0) / 6.0
    K_excited = S_LY(th + u_pole) * S_LY(th - u_pole) * K_LY(th, b_val)
    err_exc = abs(K_excited - K_LY(th, b_val + 2.0))
    report(
        f"Boundary bound-state bootstrap K_b^(1)(theta) == K_{{b+2}}(theta) (err={err_exc:.2e})",
        err_exc < 1e-13,
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
    # Panel (a): 90-degree Euclidean trick & Boundary Crossing-Unitarity
    # -----------------------------------------------------------------------
    ax1.set_xlim(-2.6, 2.6)
    ax1.set_ylim(-1.9, 2.3)
    ax1.axis("off")
    ax1.set_title(
        r"(a) $90^\circ$ Euclidean rotation: wall $\to$ boundary state $|B\rangle$",
        fontsize=9,
        pad=6,
    )

    wall_x = -0.65
    ax1.add_patch(
        patches.Rectangle(
            (wall_x, 0.25),
            0.22,
            1.65,
            facecolor=LIGHTGREY,
            edgecolor="none",
            hatch="///",
        )
    )
    ax1.plot([wall_x, wall_x], [0.25, 1.9], color=PURPLE, lw=2.2, zorder=4)
    ax1.text(
        wall_x + 0.08,
        1.98,
        r"wall $x=0$",
        fontsize=7.5,
        color=PURPLE,
        ha="center",
    )

    x_in = np.linspace(-2.2, wall_x, 60)
    t_in = np.linspace(0.4, 1.075, 60)
    x_out = np.linspace(wall_x, -2.2, 60)
    t_out = np.linspace(1.075, 1.75, 60)
    ax1.plot(x_in, t_in, color=BLUE, lw=1.8)
    arrow_mid(ax1, x_in, t_in, 0.55, BLUE)
    ax1.plot(x_out, t_out, color=RED, lw=1.8)
    arrow_mid(ax1, x_out, t_out, 0.55, RED)
    ax1.plot(wall_x, 1.075, "o", color=PURPLE, ms=4.5, zorder=5)
    ax1.text(-2.15, 0.28, r"$a(\theta)$", fontsize=8, color=BLUE)
    ax1.text(-2.15, 1.78, r"$b(-\theta)$", fontsize=8, color=RED)
    ax1.text(
        wall_x - 0.62,
        1.075,
        r"$K_a^b(\theta)$",
        fontsize=7.5,
        color=GREY,
        va="center",
    )

    ax1.annotate(
        "",
        xy=(0.45, 1.08),
        xytext=(-0.2, 1.08),
        arrowprops=dict(
            arrowstyle="->,head_width=0.18,head_length=0.22",
            color=GREY,
            lw=1.2,
        ),
    )
    ax1.text(0.125, 1.25, r"$90^\circ$", fontsize=7.5, color=GREY, ha="center")
    ax1.text(
        0.125,
        0.82,
        r"$\theta\to i\frac{\pi}{2}-\theta$",
        fontsize=7,
        color=GREY,
        ha="center",
    )

    base_y = 0.45
    ax1.add_patch(
        patches.Rectangle(
            (0.65, base_y - 0.2),
            1.85,
            0.2,
            facecolor=LIGHTGREY,
            edgecolor="none",
            hatch="///",
        )
    )
    ax1.plot([0.65, 2.5], [base_y, base_y], color=PURPLE, lw=2.2, zorder=4)
    ax1.text(
        1.58,
        base_y - 0.38,
        r"boundary state $|B\rangle$ at $t_E=0$",
        fontsize=7.5,
        color=PURPLE,
        ha="center",
    )

    v_x, v_y = 1.58, base_y
    ax1.plot(v_x, v_y, "o", color=PURPLE, ms=4.5, zorder=5)
    s_t = np.linspace(0, 1, 80)
    c1_x = v_x + 0.55 * np.sin(np.pi * s_t) - 0.65 * s_t**1.4
    c1_y = v_y + 1.38 * s_t
    c2_x = v_x - 0.55 * np.sin(np.pi * s_t) + 0.65 * s_t**1.4
    c2_y = v_y + 1.38 * s_t
    ax1.plot(c1_x, c1_y, color=BLUE, lw=1.8)
    arrow_mid(ax1, c1_x, c1_y, 0.82, BLUE)
    ax1.plot(c2_x, c2_y, color=RED, lw=1.8)
    arrow_mid(ax1, c2_x, c2_y, 0.82, RED)
    cross_y = v_y + 1.38 * 0.69
    ax1.plot(v_x, cross_y, "o", color=GREEN, ms=4.5, zorder=6)
    ax1.text(
        v_x + 0.45,
        cross_y,
        r"$S(2\theta)$",
        fontsize=8,
        color=GREEN,
        va="center",
    )
    ax1.text(
        v_x + 0.42,
        v_y + 0.25,
        r"$K(i\frac{\pi}{2}+\theta)$",
        fontsize=7.5,
        color=GREY,
    )

    ax1.text(
        0.0,
        -1.15,
        r"Boundary crossing-unitarity:"
        "\n"
        r"$K_a\!\left(i\frac{\pi}{2}-\theta\right)"
        r" = S_{ba}(2\theta)\,K_b\!\left(i\frac{\pi}{2}+\theta\right)$",
        fontsize=8,
        color=GREY,
        ha="center",
        va="center",
        bbox=dict(
            boxstyle="round,pad=0.4",
            facecolor="none",
            edgecolor=GREY,
            lw=0.8,
        ),
    )

    # -----------------------------------------------------------------------
    # Panel (b): Scaling Lee-Yang pole trajectories in the physical strip
    # -----------------------------------------------------------------------
    b_grid = np.linspace(-3.0, 2.0, 400)
    ax2.axhspan(0.0, 3.0, color=LIGHTGREY, alpha=0.35, label="Physical strip")
    ax2.axhline(3.0, color=GREY, lw=0.8, ls="--")
    ax2.axhline(0.0, color=GREY, lw=0.8)

    ax2.plot(
        b_grid,
        np.ones_like(b_grid),
        color=GREY,
        lw=1.2,
        ls=":",
        label=r"Fixed one-point pole $\theta = i\frac{\pi}{6}$",
    )

    mask1 = (b_grid >= -1.0) & (b_grid <= 2.0)
    ax2.plot(
        b_grid[mask1],
        b_grid[mask1] + 1.0,
        color=GREEN,
        lw=1.8,
        label=r"Bound-state pole $i\frac{\pi}{6}(b+1)$",
    )

    mask2 = (b_grid >= 1.0) & (b_grid <= 2.0)
    ax2.plot(
        b_grid[mask2],
        b_grid[mask2] - 1.0,
        color=ORANGE,
        lw=1.8,
        ls="-.",
        label=r"Second pole $i\frac{\pi}{6}(b-1)$ ($b>1$)",
    )

    ax2.plot(0.0, 1.0, "o", color=PURPLE, ms=5.5, zorder=6)
    ax2.annotate(
        r"$b=0$: identity wall" "\n" r"($(1-b)^{-1}$ zero cancels $p_1$)",
        xy=(0.0, 1.0),
        xytext=(-2.7, 1.85),
        fontsize=7.0,
        color=PURPLE,
        arrowprops=dict(arrowstyle="->", color=PURPLE, lw=0.8),
    )

    ax2.plot(-1.0, 0.0, "s", color=GREEN, ms=4.5, zorder=6)
    ax2.annotate(
        r"threshold $b=-1$",
        xy=(-1.0, 0.05),
        xytext=(-2.65, 0.55),
        fontsize=7.5,
        color=GREEN,
        arrowprops=dict(arrowstyle="->", color=GREEN, lw=0.8),
    )

    ax2.set_xlim(-3.0, 2.0)
    ax2.set_ylim(-0.25, 3.35)
    ax2.set_xlabel(r"Boundary parameter $b$")
    ax2.set_ylabel(r"Pole position $\operatorname{Im}\theta \,/\, (\pi/6)$")
    ax2.set_yticks([0, 1, 2, 3])
    ax2.set_yticklabels([r"$0$", r"$\pi/6$", r"$\pi/3$", r"$\pi/2$"])
    ax2.set_title(r"(b) Scaling Lee–Yang poles of $K_b(\theta)$", fontsize=9, pad=6)
    ax2.legend(loc="upper left", fontsize=6.8)

    fig.tight_layout()
    save(fig, "l1_f1_reflection_crossing", OUT_DIR)


if __name__ == "__main__":
    run_checks()
    make_figure()
    sys.exit(0 if all_ok() else 1)
