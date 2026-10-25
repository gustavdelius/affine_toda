#!/usr/bin/env python3
"""
Half-Line Lecture 4 figure & verification script:
  lecture-notes/half-line/figures/l4_f1_coideal_sg.svg

Verifies, at sampled generic parameters (lambda = sqrt 2, so q is not a root of unity):
  1. Left coideal coproduct Delta(b_j) = (b_j - eps_j k_j) x 1 + k_j x b_j
     in V(y_1) x V(y_2) (prp-coideal-left).
  2. Sine-Gordon intertwiner K(y) pi_y(b_j) = pi_{1/y}(b_j) K(y) and the braid reflection
     equation (prp-sg-k-unique); with general coefficients c_j and eps = 0, a nonzero K exists
     for generic y iff c_0 c_1 = q^-2 (smallest singular value of the 8x4 system), and
     c_0 = c_1 = -1/q is also allowed (Exercise 4.1).
  3. Sign of q: the six-vertex braid matrix obeys R(-q) = (Pi x 1) R(q) (Pi x 1), and
     K(eps; -q) = K(-eps; q). The matrix part of boundary crossing-unitarity,
     k(theta) ~ R(2 theta) k(-theta) with k = sum K_{abar}^b(i pi/2 - theta) v_a x v_b,
     holds at q = -exp(-i pi lambda) with plain C, fails at +exp(-i pi lambda) with plain C,
     and holds there with the twisted C = ((0,1),(-1,0)). Proportionality only: the scalar
     factor is not tested.
  4. Bound-state submodules: at y_n = i(-1)^n e^{-i eta}, v_n = (1, (-1)^n e^{-vartheta}) is a
     common eigenvector of pi(b_1), pi(b_0) with beta_{1,0} = -cos(bar_eta +- i vartheta)/sin(pi lambda),
     bar_eta = pi(lambda+1) - eta; common eigenvectors exist exactly at
     y^2 in {-e^{+-2 i eta}, -e^{+-2 vartheta}} (roots of the resultant).
  5. Parametrization and Neumann point: K/(-4i) = [[P0+, Q0], [Q0, P0-]] with (eta, vartheta);
     at eta = pi(lambda+1)/2 the coideal value is 1/(2 cos(pi lambda/2)), nu_n = pi/2 - n pi/lambda,
     E_n - E_0 = M sin(n pi/lambda) = m_{2n}/2, interior poles exactly for 0 < n < lambda/2.
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
# Sine-Gordon quantum group & coideal helpers
# ---------------------------------------------------------------------------
def sg_rep(y, q):
    """
    2x2 evaluation representation V(y) of U_q(a_1^(1)) in the principal gradation:
      e_1 = y E, f_1 = F/y, k_1 = diag(q, q^-1)
      e_0 = y F, f_0 = E/y, k_0 = diag(q^-1, q)
    """
    E = np.array([[0.0, 1.0], [0.0, 0.0]], dtype=complex)
    F = np.array([[0.0, 0.0], [1.0, 0.0]], dtype=complex)
    k1 = np.diag([q, 1.0 / q]).astype(complex)
    k0 = np.diag([1.0 / q, q]).astype(complex)
    return {
        "e1": y * E,
        "f1": F / y,
        "k1": k1,
        "e0": y * F,
        "f0": E / y,
        "k0": k0,
    }


def sg_b_gens(y, q, eps0, eps1, c0=None, c1=None):
    """Coideal generators b_j = e_j + c_j f_j k_j + eps_j k_j on V(y)."""
    if c0 is None:
        c0 = 1.0 / q
    if c1 is None:
        c1 = 1.0 / q
    r = sg_rep(y, q)
    b0 = r["e0"] + c0 * (r["f0"] @ r["k0"]) + eps0 * r["k0"]
    b1 = r["e1"] + c1 * (r["f1"] @ r["k1"]) + eps1 * r["k1"]
    return b0, b1


def K_sg_coideal(y, q, eps0, eps1):
    """Ghoshal-Zamolodchikov K-matrix from eq-sg-k-coideal."""
    k11 = (q - 1.0 / q) * (eps1 * y + eps0 / y)
    k22 = (q - 1.0 / q) * (eps0 * y + eps1 / y)
    off = y**2 - 1.0 / (y**2)
    return np.array([[k11, off], [off, k22]], dtype=complex)


def R_check_sg(z, q):
    """
    Bulk 6-vertex braid R-matrix R_check(z) on V(y_1) x V(y_2) where z = y_1 / y_2,
    normalized to 1 on highest weights |s> x |s>.
    Intertwines the Part VI coproduct Delta(e_j) = e_j x 1 + k_j x e_j.
    """
    denom = q * z - 1.0 / (q * z)
    b_val = (z - 1.0 / z) / denom
    c12 = (q - 1.0 / q) / denom
    return np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, c12, b_val, 0.0],
            [0.0, b_val, c12, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=complex,
    )


# ---------------------------------------------------------------------------
# Verification checks
# ---------------------------------------------------------------------------
def run_checks():
    lam = np.sqrt(2.0)  # irrational lambda to avoid roots of unity
    q_phys = -np.exp(-1.0j * np.pi * lam)
    y1 = np.exp(lam * (0.37 + 0.14j))
    y2 = np.exp(lam * (-0.29 + 0.21j))
    eps0, eps1 = 0.53 - 0.19j, -0.41 + 0.31j
    I2 = np.eye(2, dtype=complex)

    # Check 1: Left coideal coproduct identity (prp-coideal-left)
    r1 = sg_rep(y1, q_phys)
    r2 = sg_rep(y2, q_phys)
    b0_1, b1_1 = sg_b_gens(y1, q_phys, eps0, eps1)
    b0_2, b1_2 = sg_b_gens(y2, q_phys, eps0, eps1)

    De1 = np.kron(r1["e1"], I2) + np.kron(r1["k1"], r2["e1"])
    Df1 = np.kron(r1["f1"], np.linalg.inv(r2["k1"])) + np.kron(I2, r2["f1"])
    Dk1 = np.kron(r1["k1"], r2["k1"])
    Db1_direct = De1 + (1.0 / q_phys) * (Df1 @ Dk1) + eps1 * Dk1
    Db1_coideal = np.kron(b1_1 - eps1 * r1["k1"], I2) + np.kron(r1["k1"], b1_2)
    err_coid = np.max(np.abs(Db1_direct - Db1_coideal))
    report(
        f"Left coideal coproduct Delta(b_j) == (b_j - eps_j k_j) x 1 + k_j x b_j (err={err_coid:.2e})",
        err_coid < 1e-14,
    )

    # Check 2: Intertwining equation, c0*c1 = q^-2 condition, and reflection equation
    b0_y1, b1_y1 = sg_b_gens(y1, q_phys, eps0, eps1)
    b0_iy1, b1_iy1 = sg_b_gens(1.0 / y1, q_phys, eps0, eps1)
    K1 = K_sg_coideal(y1, q_phys, eps0, eps1)
    K2 = K_sg_coideal(y2, q_phys, eps0, eps1)
    err_int = max(
        np.max(np.abs(K1 @ b0_y1 - b0_iy1 @ K1)),
        np.max(np.abs(K1 @ b1_y1 - b1_iy1 @ K1)),
    )
    R12_diff = R_check_sg(y1 / y2, q_phys)
    R12_sum = R_check_sg(y1 * y2, q_phys)
    lhs_re = R12_diff @ np.kron(I2, K1) @ R12_sum @ np.kron(I2, K2)
    rhs_re = np.kron(I2, K2) @ R12_sum @ np.kron(I2, K1) @ R12_diff
    err_re = np.max(np.abs(lhs_re - rhs_re))
    report(
        f"Sine-Gordon K(y) solves linear intertwining (err={err_int:.2e}) and braid reflection eq (err={err_re:.2e})",
        max(err_int, err_re) < 1e-13,
    )

    def smin_system(y, q, c0, c1):
        """Smallest singular value of K -> (K b_j(y) - b_j(1/y) K)_{j=0,1} at eps = 0."""
        b_y = sg_b_gens(y, q, 0.0, 0.0, c0, c1)
        b_iy = sg_b_gens(1.0 / y, q, 0.0, 0.0, c0, c1)
        M = np.zeros((8, 4), dtype=complex)
        for idx in range(4):
            Kb = np.zeros((2, 2), dtype=complex)
            Kb.flat[idx] = 1.0
            M[:, idx] = np.concatenate([(Kb @ b_y[j] - b_iy[j] @ Kb).ravel() for j in range(2)])
        return np.linalg.svd(M, compute_uv=False)[-1]

    s_on = smin_system(y1, q_phys, 2.3 / q_phys, 1.0 / (2.3 * q_phys))
    s_off = smin_system(y1, q_phys, 2.3 / q_phys, 1.1 / (2.3 * q_phys))
    s_neg = smin_system(y1, q_phys, -1.0 / q_phys, -1.0 / q_phys)
    report(
        f"eps=0, general c_j: nonzero K for generic y iff c0*c1 = q^-2 (sigma_min on {s_on:.1e}, "
        f"off {s_off:.3f}); c0=c1=-1/q also allowed ({s_neg:.1e})",
        s_on < 1e-13 and s_off > 1e-3 and s_neg < 1e-13,
    )

    # Check 3: sign of q. Braid matrix and K under q -> -q, and the matrix part of
    # boundary crossing-unitarity k(theta) ~ R(2 theta) k(-theta).
    Pi1 = np.kron(np.diag([1.0, -1.0]), I2).astype(complex)
    z_t = y1 / y2
    err_Rq = np.max(np.abs(R_check_sg(z_t, -q_phys) - Pi1 @ R_check_sg(z_t, q_phys) @ Pi1))
    err_Kq = np.max(np.abs(K_sg_coideal(y1, -q_phys, eps0, eps1) - K_sg_coideal(y1, q_phys, -eps0, -eps1)))
    report(
        f"q -> -q: R(-q) = (Pi x 1) R(q) (Pi x 1) (err={err_Rq:.1e}) and K(eps;-q) = K(-eps;q) (err={err_Kq:.1e})",
        max(err_Rq, err_Kq) < 1e-13,
    )

    th = 0.37 - 0.12j
    C_plain = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    C_twist = np.array([[0.0, 1.0], [-1.0, 0.0]], dtype=complex)

    def crossing_err(q_choice, C_mat):
        y_minus = np.exp(lam * (0.5j * np.pi - th))
        y_plus = np.exp(lam * (0.5j * np.pi + th))
        Km = K_sg_coideal(y_minus, q_phys, eps0, eps1)
        Kp = K_sg_coideal(y_plus, q_phys, eps0, eps1)
        # k = sum_{a,b} K_{abar}^b v_a x v_b, K_{abar}^b = (K C)[b, a]; index (a,b) -> 2a+b
        km = (Km @ C_mat).T.reshape(4)
        kp = (Kp @ C_mat).T.reshape(4)
        z_2th = np.exp(2.0 * lam * th)
        lhs_vec = R_check_sg(z_2th, q_choice) @ kp
        ratio = np.vdot(km, lhs_vec) / np.vdot(km, km)
        return np.linalg.norm(lhs_vec - ratio * km) / np.linalg.norm(lhs_vec)

    err_q_phys = crossing_err(q_phys, C_plain)
    err_q_wrong = crossing_err(-q_phys, C_plain)
    err_q_twist = crossing_err(-q_phys, C_twist)
    report(
        f"Matrix part of boundary crossing-unitarity: plain C at q=-e^{{-i*pi*lam}} (err={err_q_phys:.1e}), "
        f"fails at +e^{{-i*pi*lam}} (err={err_q_wrong:.2f}), twisted C restores it there (err={err_q_twist:.1e})",
        err_q_phys < 1e-13 and err_q_wrong > 0.1 and err_q_twist < 1e-13,
    )

    # Check 4: one-dimensional B-submodules of V(y) x C_eps (prp-boundary-module)
    eta = 1.85
    vartheta = 0.42
    e1_eta = -np.cos(eta - 1.0j * vartheta) / np.sin(np.pi * lam)
    e0_eta = -np.cos(eta + 1.0j * vartheta) / np.sin(np.pi * lam)
    bar_eta = np.pi * (lam + 1.0) - eta
    beta1_exact = -np.cos(bar_eta + 1.0j * vartheta) / np.sin(np.pi * lam)
    beta0_exact = -np.cos(bar_eta - 1.0j * vartheta) / np.sin(np.pi * lam)

    max_err_sub = 0.0
    for n_idx in (0, 1, 2, 3):
        nu_n = eta / lam - (2 * n_idx + 1) * np.pi / (2.0 * lam)
        y_refl = np.exp(-1.0j * lam * nu_n)
        max_err_sub = max(max_err_sub, abs(y_refl - 1.0j * (-1) ** n_idx * np.exp(-1.0j * eta)))
        b0_sub, b1_sub = sg_b_gens(y_refl, q_phys, e0_eta, e1_eta)
        v = np.array([1.0, (-1) ** n_idx * np.exp(-vartheta)], dtype=complex)
        max_err_sub = max(
            max_err_sub,
            np.linalg.norm(b1_sub @ v - beta1_exact * v),
            np.linalg.norm(b0_sub @ v - beta0_exact * v),
        )
    # all branches: resultant in alpha of y^2 a^2 + y E1 a - 1 and a^2 - y E0 a - y^2,
    # with E_j = eps_j (q - 1/q); as a polynomial in Y = y^2 it is (1-Y^2)^2 - Y (E0 Y + E1)(E1 Y + E0)
    E1c = e1_eta * (q_phys - 1.0 / q_phys)
    E0c = e0_eta * (q_phys - 1.0 / q_phys)
    res_poly = np.polysub(
        np.polymul([-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0]),
        np.polymul([1.0, 0.0], np.polymul([E0c, E1c], [E1c, E0c])),
    )
    Y_roots = np.sort_complex(np.roots(res_poly))
    Y_expected = np.sort_complex(
        np.array([-np.exp(2j * eta), -np.exp(-2j * eta), -np.exp(2 * vartheta), -np.exp(-2 * vartheta)])
    )
    err_branch = np.max(np.abs(Y_roots - Y_expected))
    report(
        f"Common eigenvector v_n=(1,(-1)^n e^(-vartheta)) at y_n=i(-1)^n e^(-i eta) with beta_(1,0)=-cos(bar_eta +- i vartheta)/sin(pi lam) "
        f"(max err={max_err_sub:.1e}); all branches y^2 in {{-e^(+-2i eta), -e^(+-2 vartheta)}} (err={err_branch:.1e})",
        max_err_sub < 1e-13 and err_branch < 1e-10,
    )

    # Check 5: (eta, vartheta) form of K and the Neumann point
    u_t = 0.27 + 0.05j
    K_u = K_sg_coideal(np.exp(1j * lam * u_t), q_phys, e0_eta, e1_eta)
    Pp = np.cos(lam * u_t) * np.cos(eta) * np.cosh(vartheta) - np.sin(lam * u_t) * np.sin(eta) * np.sinh(vartheta)
    Pm = np.cos(lam * u_t) * np.cos(eta) * np.cosh(vartheta) + np.sin(lam * u_t) * np.sin(eta) * np.sinh(vartheta)
    Q0 = -np.sin(lam * u_t) * np.cos(lam * u_t)
    err_param = np.max(np.abs(K_u / (-4j) - np.array([[Pp, Q0], [Q0, Pm]])))
    eps_N = -np.cos(0.5 * np.pi * (lam + 1.0)) / np.sin(np.pi * lam)
    err_N = max(
        abs(eps_N - 1.0 / (2.0 * np.cos(0.5 * np.pi * lam))),
        abs((1 - q_phys) * (1 - 1 / q_phys) - 4.0 * np.cos(0.5 * np.pi * lam) ** 2),
    )
    lam_b = 7.3  # pole kinematics only (no representation theory involved)
    eta_N = 0.5 * np.pi * (lam_b + 1.0)
    inside = []
    for n_idx in range(10):
        nu_n = eta_N / lam_b - (2 * n_idx + 1) * np.pi / (2.0 * lam_b)
        m_k = lambda k: 2.0 * np.sin(k * np.pi / (2.0 * lam_b))
        err_N = max(
            err_N,
            abs(nu_n - (0.5 * np.pi - n_idx * np.pi / lam_b)),
            abs(np.cos(nu_n) - 0.5 * m_k(2 * n_idx)),
            abs(np.cos(nu_n) - m_k(n_idx) * np.cos(n_idx * np.pi / (2.0 * lam_b))),
        )
        if 0.0 < nu_n < 0.5 * np.pi:
            inside.append(n_idx)
    ok_range = inside == [n for n in range(1, 10) if n < lam_b / 2.0]
    report(
        f"K/(-4i) = [[P0+,Q0],[Q0,P0-]] (err={err_param:.1e}); Neumann eps_*=1/(2cos(pi lam/2)), "
        f"nu_n=pi/2-n pi/lam, E_n=m_(2n)/2 (err={err_N:.1e}), interior poles n={inside}",
        err_param < 1e-13 and err_N < 1e-13 and ok_range,
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
    # Panel (a): Left coideal action on N solitons + boundary module V_[lambda]
    # -----------------------------------------------------------------------
    ax1.set_xlim(-2.9, 0.8)
    ax1.set_ylim(-1.5, 1.7)
    ax1.axis("off")
    ax1.set_title(
        r"(a) Why the wall at $x=0$ forces a left coideal",
        fontsize=8.8,
        pad=6,
    )

    ax1.add_patch(
        patches.Rectangle(
            (0.0, -0.65),
            0.25,
            1.65,
            facecolor=LIGHTGREY,
            edgecolor="none",
            hatch="///",
        )
    )
    ax1.plot([0.0, 0.0], [-0.65, 1.0], color=PURPLE, lw=2.2, zorder=4)
    ax1.text(
        0.12,
        1.15,
        r"Wall $V_{[\lambda]}$" "\n" r"($B_{\hat\epsilon}$-module)",
        fontsize=7.3,
        color=PURPLE,
        ha="center",
    )

    ax1.plot([-2.7, 0.0], [0.15, 0.15], color=GREY, lw=1.0)
    xs = [-2.25, -1.4, -0.55]
    cols = [BLUE, GREEN, ORANGE]
    labs = [
        r"$V_1(\theta_1)$",
        r"$V_2(\theta_2)$",
        r"$V_N(\theta_N)$",
    ]
    sub_labs = [
        r"$U_q(\hat{g})$",
        r"$U_q(\hat{g})$",
        r"$U_q(\hat{g})$",
    ]
    for x_pos, c_val, l_top, l_bot in zip(xs, cols, labs, sub_labs):
        ax1.plot(x_pos, 0.15, "o", color=c_val, ms=6.5, zorder=5)
        ax1.text(
            x_pos, 0.42, l_top, fontsize=7.6, color=c_val, ha="center"
        )
        ax1.text(
            x_pos, -0.22, l_bot, fontsize=7.2, color=GREY, ha="center"
        )

    ax1.text(-0.97, 0.15, r"$\cdots$", fontsize=9, color=GREY, ha="center", va="center")
    ax1.text(
        0.12,
        -0.22,
        r"$B_{\hat\epsilon}$",
        fontsize=7.8,
        color=PURPLE,
        ha="center",
    )

    ax1.annotate(
        "",
        xy=(-0.1, 0.82),
        xytext=(-2.45, 0.82),
        arrowprops=dict(
            arrowstyle="->,head_width=0.16,head_length=0.2",
            color=RED,
            lw=1.2,
            ls="--",
        ),
    )
    ax1.text(
        -1.28,
        0.95,
        r"topological tail $k_j = q^{T_j}$ accumulates to the left",
        fontsize=7.0,
        color=RED,
        ha="center",
    )

    ax1.text(
        -1.05,
        -1.05,
        r"$\Delta(b_j) = (b_j - \hat\epsilon_j k_j)\otimes 1 + k_j\otimes b_j"
        r"\;\in\; U_q(\hat{g})\otimes B_{\hat\epsilon}$"
        "\n"
        r"On $V(\theta)\otimes\mathbb{C}_{\beta}$: $\;\hat\epsilon_j \longmapsto \beta_j$ (bound-state shift!)",
        fontsize=7.3,
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
    # Panel (b): Sine-Gordon UV-IR map and bound-state poles nu_n(eta)
    # -----------------------------------------------------------------------
    lam_plot = 3.5
    eta_max = 0.5 * np.pi * (lam_plot + 1.0)
    eta_grid = np.linspace(0.0, eta_max, 300)

    ax2.axhspan(0.0, 0.5 * np.pi, color=LIGHTGREY, alpha=0.35, label="Physical strip")
    ax2.axhline(0.5 * np.pi, color=GREY, lw=0.8, ls="--")
    ax2.axhline(0.0, color=GREY, lw=0.8)

    pal = [BLUE, GREEN, ORANGE, RED]
    for n_idx in range(4):
        nu_curve = eta_grid / lam_plot - (2 * n_idx + 1) * np.pi / (2.0 * lam_plot)
        mask = (nu_curve >= 0.0) & (nu_curve <= 0.5 * np.pi)
        ax2.plot(
            eta_grid[mask] / np.pi,
            nu_curve[mask],
            color=pal[n_idx],
            lw=1.8,
            label=rf"Bound state $|{n_idx}\rangle$: $\nu_{n_idx}(\eta)$",
        )

    ax2.axvline(0.5, color=PURPLE, lw=1.1, ls=":")
    ax2.text(
        0.53,
        0.43 * np.pi,
        r"$\eta=\frac{\pi}{2}$ ($\hat\epsilon_j=0$," "\n" r"classical $C\to +1$)",
        fontsize=6.8,
        color=PURPLE,
        va="top",
    )

    ax2.plot(
        0.5 * (lam_plot + 1.0),
        0.5 * np.pi,
        "o",
        color=PURPLE,
        ms=5.5,
        zorder=6,
    )
    ax2.annotate(
        r"Neumann $\eta=\frac{\pi}{2}(\lambda+1)$" "\n" r"($\hat\epsilon_j=\epsilon_*$)",
        xy=(0.5 * (lam_plot + 1.0), 0.5 * np.pi),
        xytext=(1.05, 0.34 * np.pi),
        fontsize=6.8,
        color=PURPLE,
        arrowprops=dict(arrowstyle="->", color=PURPLE, lw=0.8),
    )

    ax2.set_xlim(0.0, 0.5 * (lam_plot + 1.0) + 0.15)
    ax2.set_ylim(-0.08, 0.56 * np.pi)
    ax2.set_xlabel(r"Boundary parameter $\eta / \pi$ ($\lambda = 3.5$)")
    ax2.set_ylabel(r"Pole position $\nu_n$ in $u = -i\theta$")
    ax2.set_yticks([0.0, 0.25 * np.pi, 0.5 * np.pi])
    ax2.set_yticklabels([r"$0$", r"$\pi/4$", r"$\pi/2$"])
    ax2.set_title(
        r"(b) Sine-Gordon boundary bound states $\nu_n(\eta)$",
        fontsize=8.8,
        pad=6,
    )
    ax2.legend(loc="lower right", fontsize=6.6)

    fig.tight_layout()
    save(fig, "l4_f1_coideal_sg", OUT_DIR)


if __name__ == "__main__":
    run_checks()
    make_figure()
    sys.exit(0 if all_ok() else 1)
