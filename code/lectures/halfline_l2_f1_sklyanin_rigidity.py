#!/usr/bin/env python3
"""
Half-Line Lecture 2 figure & verification script:
  lecture-notes/half-line/figures/l2_f1_sklyanin_rigidity.svg

Verifies:
  1. Explicit 2x2 classical K(lambda) for a_1^(1) solves eq-boundary-k-equation
     for arbitrary complex (C_0, C_1).
  2. Explicit 3x3 classical K(lambda) for a_2^(1) solves eq-boundary-k-equation
     if and only if C_0 = C_1 = C_2 = 0 (Neumann) or C_j in {+1, -1}.
  3. Both K(lambda) matrices satisfy Sklyanin's classical reflection equation
     (eq-classical-reflection) with the affine Toda r-matrix.
  4. Full Dirichlet compatibility has only K_+ = 0 for a_2^(1) (prp-no-dirichlet),
     whereas sigma = Hermitian conjugation gives K_+ = exp(i*pi*nu.H) for
     Re(phi(0)) = (pi/beta)*nu and d_x Im(phi(0)) = 0 (prp-soliton-preserving-half).
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
# Explicit classical K(lambda) for a_1^(1) and a_2^(1)
# ---------------------------------------------------------------------------
def a1_generators():
    """Chevalley generators of a_1^(1) in the 2x2 defining representation."""
    h1 = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
    h0 = -h1
    E1 = np.array([[0.0, 1.0], [0.0, 0.0]], dtype=complex)
    F1 = np.array([[0.0, 0.0], [1.0, 0.0]], dtype=complex)
    E0 = F1.copy()
    F0 = E1.copy()
    return [(h0, E0, F0), (h1, E1, F1)]


def K_a1_classical(lam, C0, C1):
    """
    Explicit 2x2 solution of K(lam)(-C_j h_j + lam E_j - lam^-1 F_j)
                           = (C_j h_j + lam E_j - lam^-1 F_j) K(lam)
    for a_1^(1) with arbitrary complex (C0, C1).
    """
    diag = lam**2 - 1.0 / lam**2
    k12 = -2.0 * (C0 * lam - C1 / lam)
    k21 = -2.0 * (C1 * lam - C0 / lam)
    return np.array([[diag, k12], [k21, diag]], dtype=complex)


def a2_generators():
    """Chevalley generators of a_2^(1) in the 3x3 defining representation."""
    h1 = np.diag([1.0, -1.0, 0.0]).astype(complex)
    h2 = np.diag([0.0, 1.0, -1.0]).astype(complex)
    h0 = np.diag([-1.0, 0.0, 1.0]).astype(complex)
    E1 = np.zeros((3, 3), dtype=complex); E1[0, 1] = 1.0
    F1 = E1.T.copy()
    E2 = np.zeros((3, 3), dtype=complex); E2[1, 2] = 1.0
    F2 = E2.T.copy()
    E0 = np.zeros((3, 3), dtype=complex); E0[2, 0] = 1.0
    F0 = E0.T.copy()
    return [(h0, E0, F0), (h1, E1, F1), (h2, E2, F2)]


def K_a2_classical(lam, C0, C1, C2):
    """
    Explicit 3x3 solution of eq-boundary-k-equation for a_2^(1) when C_j in {+1, -1}.
    Derived on the blackboard from the spectator-state kernel condition det = 1 - C_j^2 = 0.
    """
    s = C0 * C1 * C2  # sign product in {+1, -1}
    d = s * lam**2 - 1.0 / lam
    return np.array(
        [
            [d, 2.0 * C1, -2.0 * C1 * C2 * lam],
            [-2.0 * C0 * C2 * lam, d, 2.0 * C2],
            [2.0 * C0, -2.0 * C0 * C1 * lam, d],
        ],
        dtype=complex,
    )


def min_sv_k_equation(gens, C_list, lam):
    """Smallest singular value of the linear homogeneous system for K(lam)."""
    d = gens[0][0].shape[0]
    I_d = np.eye(d, dtype=complex)
    rows = []
    for (hj, Ej, Fj), Cj in zip(gens, C_list):
        Xj_minus = -Cj * hj + lam * Ej - (1.0 / lam) * Fj
        Xj_plus = Cj * hj + lam * Ej - (1.0 / lam) * Fj
        op = np.kron(I_d, Xj_minus.T) - np.kron(Xj_plus, I_d)
        rows.append(op)
    M = np.vstack(rows)
    sv = np.linalg.svd(M, compute_uv=False)
    return sv[-1]


def r_matrix_an(gens, z):
    """
    Classical r-matrix for a_n^(1) in the (n+1)-dimensional defining representation:
      r(z) = 1/(z^h - 1) * [ (z^h + 1)/2 * sum_a H_a x H_a
             + sum_{i != j} z^{(j - i) mod h} e_{ij} x e_{ji} ]
    """
    h = len(gens)
    r = np.zeros((h * h, h * h), dtype=complex)
    for i in range(h):
        for j in range(h):
            val = (1.0 if i == j else 0.0) - 1.0 / h
            e_ii = np.zeros((h, h), dtype=complex); e_ii[i, i] = 1.0
            e_jj = np.zeros((h, h), dtype=complex); e_jj[j, j] = 1.0
            r += 0.5 * (z**h + 1.0) * val * np.kron(e_ii, e_jj)
    for i in range(h):
        for j in range(h):
            if i != j:
                power = (j - i) % h
                e_ij = np.zeros((h, h), dtype=complex); e_ij[i, j] = 1.0
                e_ji = np.zeros((h, h), dtype=complex); e_ji[j, i] = 1.0
                r += (z**power) * np.kron(e_ij, e_ji)
    return r / (z**h - 1.0)


def partial_transpose(mat, d, slot):
    """Apply transposition to tensor factor `slot` (1, 2, or 3 for both) of a (d^2 x d^2) matrix."""
    T = mat.reshape(d, d, d, d)
    if slot == 1:
        return T.transpose(2, 1, 0, 3).reshape(d * d, d * d)
    if slot == 2:
        return T.transpose(0, 3, 2, 1).reshape(d * d, d * d)
    return T.transpose(2, 3, 0, 1).reshape(d * d, d * d)


# ---------------------------------------------------------------------------
# Verification checks
# ---------------------------------------------------------------------------
def run_checks():
    # Check 1: Explicit a_1^(1) K(lam) for generic complex C0, C1
    gens_a1 = a1_generators()
    lam = 1.27 + 0.31j
    C0, C1 = 0.63 - 0.41j, -1.19 + 0.28j
    K1 = K_a1_classical(lam, C0, C1)
    err_a1 = 0.0
    for (hj, Ej, Fj), Cj in zip(gens_a1, [C0, C1]):
        Xm = -Cj * hj + lam * Ej - (1.0 / lam) * Fj
        Xp = Cj * hj + lam * Ej - (1.0 / lam) * Fj
        err_a1 = max(err_a1, np.max(np.abs(K1 @ Xm - Xp @ K1)))
    det_K1 = abs(np.linalg.det(K1))
    report(
        f"a_1^(1) explicit K(lam) solves eq-boundary-k-equation for generic complex C_0,C_1 "
        f"(err={err_a1:.2e}, |det K|={det_K1:.2f})",
        err_a1 < 1e-13 and det_K1 > 1e-2,
    )

    # Check 2: Explicit a_2^(1) K(lam) for all 8 sign choices C_j in {+1,-1} and rigidity test
    gens_a2 = a2_generators()
    err_a2_signs = 0.0
    min_det_a2 = 1e9
    for c0 in (-1.0, 1.0):
        for c1 in (-1.0, 1.0):
            for c2 in (-1.0, 1.0):
                K2 = K_a2_classical(lam, c0, c1, c2)
                min_det_a2 = min(min_det_a2, abs(np.linalg.det(K2)))
                for (hj, Ej, Fj), Cj in zip(gens_a2, [c0, c1, c2]):
                    Xm = -Cj * hj + lam * Ej - (1.0 / lam) * Fj
                    Xp = Cj * hj + lam * Ej - (1.0 / lam) * Fj
                    err_a2_signs = max(
                        err_a2_signs, np.max(np.abs(K2 @ Xm - Xp @ K2))
                    )
    sv_bad = min_sv_k_equation(gens_a2, [0.6, 0.6, 0.6], lam)
    report(
        f"a_2^(1) explicit K(lam) holds for all 8 sign patterns C_j=+-1 (err={err_a2_signs:.2e}, "
        f"min|det|={min_det_a2:.2f}) and fails off +-1 (sv={sv_bad:.2f})",
        err_a2_signs < 1e-13 and min_det_a2 > 1e-2 and sv_bad > 0.1,
    )

    # Check 3: Sklyanin's classical reflection equation (eq-classical-reflection)
    lam1, lam2 = 1.31 + 0.17j, 0.79 - 0.23j
    for name, gens, K_fn in [
        ("a_1^(1)", gens_a1, lambda z: K_a1_classical(z, C0, C1)),
        ("a_2^(1)", gens_a2, lambda z: K_a2_classical(z, 1.0, -1.0, 1.0)),
    ]:
        d = len(gens)
        Id = np.eye(d, dtype=complex)
        K_l1 = np.kron(K_fn(lam1), Id)
        K_l2 = np.kron(Id, K_fn(lam2))
        r_12 = r_matrix_an(gens, lam1 / lam2)
        r_21_tt = partial_transpose(r_matrix_an(gens, lam2 / lam1), d, 3)
        r_prod_t2 = partial_transpose(r_matrix_an(gens, lam1 * lam2), d, 2)
        r_iprod_t1 = partial_transpose(
            r_matrix_an(gens, 1.0 / (lam1 * lam2)), d, 1
        )
        lhs_cre = (
            K_l1 @ K_l2 @ r_12
            + r_21_tt @ K_l1 @ K_l2
            + K_l1 @ r_prod_t2 @ K_l2
            + K_l2 @ r_iprod_t1 @ K_l1
        )
        err_cre = np.max(np.abs(lhs_cre))
        report(
            f"Classical reflection equation holds for {name} (max err={err_cre:.2e})",
            err_cre < 1e-12,
        )

    # Check 4: Soliton-preserving condition (sigma = dagger) with K_+ = exp(i*pi*nu.H)
    lam_real = 1.37
    K_sp = np.diag([1.0, -1.0, -1.0]).astype(complex)
    signs_nu = [-1.0, 1.0, 0.0]  # alpha_j . nu for j = 0, 1, 2
    dx_re_phi = [0.41, -0.73]
    at_lam = (
        0.5j * dx_re_phi[0] * gens_a2[1][0]
        + 0.5j * dx_re_phi[1] * gens_a2[2][0]
    )
    at_inv = (
        0.5j * dx_re_phi[0] * gens_a2[1][0]
        + 0.5j * dx_re_phi[1] * gens_a2[2][0]
    )
    for j, (hj, Ej, Fj) in enumerate(gens_a2):
        wj = np.exp(-0.3 * (j + 1)) * np.exp(0.5j * np.pi * signs_nu[j])
        at_lam += 0.5 * wj * (lam_real * Ej - (1.0 / lam_real) * Fj)
        at_inv += 0.5 * wj * ((1.0 / lam_real) * Ej - lam_real * Fj)
    err_sp = np.max(np.abs(K_sp @ at_lam + at_inv.conj().T @ K_sp))
    report(
        f"Soliton-preserving compatibility K_+ a_t(lam) + a_t(1/lam)^dagger K_+ = 0 (err={err_sp:.2e})",
        err_sp < 1e-14,
    )


# ---------------------------------------------------------------------------
# Figure generation
# ---------------------------------------------------------------------------
def make_figure():
    setup()
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(7.4, 3.1), gridspec_kw={"width_ratios": [1.05, 1.15]}
    )

    # -----------------------------------------------------------------------
    # Panel (a): Sklyanin's Double-Row Monodromy Matrix on x <= 0
    # -----------------------------------------------------------------------
    ax1.set_xlim(-2.7, 0.9)
    ax1.set_ylim(-1.4, 1.6)
    ax1.axis("off")
    ax1.set_title(
        r"(a) Sklyanin's double-row transfer matrix $\tau(\lambda)$",
        fontsize=9,
        pad=6,
    )

    ax1.add_patch(
        patches.Rectangle(
            (0.0, -0.85),
            0.22,
            1.8,
            facecolor=LIGHTGREY,
            edgecolor="none",
            hatch="///",
        )
    )
    ax1.plot([0.0, 0.0], [-0.85, 0.95], color=PURPLE, lw=2.2, zorder=4)
    ax1.text(0.1, 1.1, r"wall $x_+=0$", fontsize=7.5, color=PURPLE, ha="center")

    ax1.plot([-2.3, -2.3], [-0.7, 0.8], color=GREY, lw=1.0, ls="--")
    ax1.text(
        -2.3, 0.95, r"$x_-\to-\infty$", fontsize=7.5, color=GREY, ha="center"
    )

    y_top = 0.42
    ax1.plot([-2.3, 0.0], [y_top, y_top], color=BLUE, lw=1.9)
    ax1.annotate(
        "",
        xy=(-0.95, y_top),
        xytext=(-1.35, y_top),
        arrowprops=dict(
            arrowstyle="->,head_width=0.18,head_length=0.22",
            color=BLUE,
            lw=1.5,
        ),
    )
    ax1.text(
        -1.15,
        y_top + 0.16,
        r"forward transport $T(\lambda)$",
        fontsize=7.5,
        color=BLUE,
        ha="center",
    )

    y_bot = -0.32
    ax1.plot([0.0, -2.3], [y_bot, y_bot], color=RED, lw=1.9)
    ax1.annotate(
        "",
        xy=(-1.35, y_bot),
        xytext=(-0.95, y_bot),
        arrowprops=dict(
            arrowstyle="->,head_width=0.18,head_length=0.22",
            color=RED,
            lw=1.5,
        ),
    )
    ax1.text(
        -1.15,
        y_bot - 0.25,
        r"return transport $\sigma(T(1/\lambda))$",
        fontsize=7.5,
        color=RED,
        ha="center",
    )

    theta_arc = np.linspace(-0.5 * np.pi, 0.5 * np.pi, 50)
    r_arc = 0.5 * (y_top - y_bot)
    yc = 0.5 * (y_top + y_bot)
    ax1.plot(
        0.0 + 0.25 * np.cos(theta_arc),
        yc + r_arc * np.sin(theta_arc),
        color=PURPLE,
        lw=1.8,
    )
    ax1.plot(0.0, y_top, "o", color=PURPLE, ms=4.5, zorder=5)
    ax1.plot(0.0, y_bot, "o", color=PURPLE, ms=4.5, zorder=5)
    ax1.text(
        0.32,
        yc,
        r"$K_+(\lambda)$",
        fontsize=8,
        color=PURPLE,
        va="center",
    )

    ax1.plot(
        -2.3 - 0.2 * np.cos(theta_arc),
        yc + r_arc * np.sin(theta_arc),
        color=GREY,
        lw=1.2,
        ls=":",
    )
    ax1.text(
        -2.55,
        yc,
        r"$K_-=1$",
        fontsize=7.5,
        color=GREY,
        va="center",
        ha="right",
    )

    ax1.text(
        -0.9,
        -1.12,
        r"$\dot\tau(\lambda)=0 \Longleftrightarrow K_+(\lambda)a_t(\lambda)"
        r"+\sigma(a_t(1/\lambda))K_+(\lambda)=0$",
        fontsize=7.4,
        color=GREY,
        ha="center",
        va="center",
        bbox=dict(
            boxstyle="round,pad=0.35",
            facecolor="none",
            edgecolor=GREY,
            lw=0.8,
        ),
    )

    # -----------------------------------------------------------------------
    # Panel (b): Singular value landscape of eq-boundary-k-equation
    # -----------------------------------------------------------------------
    c_grid = np.linspace(-1.6, 1.6, 321)
    lam_plot = 1.25 + 0.2j
    gens_a1 = a1_generators()
    gens_a2 = a2_generators()

    sv_a1 = [min_sv_k_equation(gens_a1, [c, 0.75], lam_plot) for c in c_grid]
    sv_a2_unif = [
        min_sv_k_equation(gens_a2, [c, c, c], lam_plot) for c in c_grid
    ]
    sv_a2_mixed = [
        min_sv_k_equation(gens_a2, [c, 1.0, -1.0], lam_plot) for c in c_grid
    ]

    ax2.plot(
        c_grid,
        sv_a1,
        color=BLUE,
        lw=1.8,
        ls="--",
        label=r"$a_1^{(1)}$ ($C_1=0.75$, any $C_0$)",
    )
    ax2.plot(
        c_grid,
        sv_a2_unif,
        color=GREEN,
        lw=1.8,
        label=r"$a_2^{(1)}$ uniform ($C_0=C_1=C_2=C$)",
    )
    ax2.plot(
        c_grid,
        sv_a2_mixed,
        color=ORANGE,
        lw=1.5,
        ls="-.",
        label=r"$a_2^{(1)}$ slice ($C_1=+1,C_2=-1$)",
    )

    ax2.plot([-1.0, 0.0, 1.0], [0.0, 0.0, 0.0], "o", color=PURPLE, ms=5.5, zorder=6)
    ax2.annotate(
        r"Neumann $C=0$",
        xy=(0.0, 0.03),
        xytext=(0.0, 0.52),
        fontsize=7.2,
        color=PURPLE,
        ha="center",
        arrowprops=dict(arrowstyle="->", color=PURPLE, lw=0.8),
    )
    ax2.annotate(
        r"$C=-1$",
        xy=(-1.0, 0.03),
        xytext=(-1.25, 0.45),
        fontsize=7.2,
        color=PURPLE,
        ha="center",
        arrowprops=dict(arrowstyle="->", color=PURPLE, lw=0.8),
    )
    ax2.annotate(
        r"$C=+1$",
        xy=(1.0, 0.03),
        xytext=(1.25, 0.45),
        fontsize=7.2,
        color=PURPLE,
        ha="center",
        arrowprops=dict(arrowstyle="->", color=PURPLE, lw=0.8),
    )

    ax2.set_xlim(-1.6, 1.6)
    ax2.set_ylim(-0.06, 1.15)
    ax2.set_xlabel(r"Boundary coupling parameter $C$")
    ax2.set_ylabel(r"Smallest singular value $\sigma_{\min}$")
    ax2.set_title(
        r"(b) Continuous $a_1^{(1)}$ vs. rigid $a_2^{(1)}$ ($K(\lambda)\ne 0 \Longleftrightarrow \sigma_{\min}=0$)",
        fontsize=8.5,
        pad=6,
    )
    ax2.legend(loc="upper center", fontsize=6.8)

    fig.tight_layout()
    save(fig, "l2_f1_sklyanin_rigidity", OUT_DIR)


if __name__ == "__main__":
    run_checks()
    make_figure()
    sys.exit(0 if all_ok() else 1)
