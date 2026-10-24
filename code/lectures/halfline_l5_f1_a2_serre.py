"""Lecture 5 verification and figure: a_2^(1) intertwiner rigidity, Serre obstruction, and q-Onsager.

Verifies:
 1. For a_2^(1), the 27x9 linear intertwining system K(y) pi_{V,y}(b_j) = pi_{Vbar,1/y}(b_j) K(y)
    has a non-trivial nullspace if and only if eps_j = 0 (giving K = I) or eps_j^2 = eps_*^2 =
    1/((1-q)(1-q^-1)) = 1/(4 cos^2(pi*omega/2)) (giving the closed-form Toeplitz matrices K_{s=+-1}).
 2. The closed-form matrices K_{s=+-1}(y) and K_{eps=0}(y) = I solve the soliton-conjugating
    reflection equation for all four pairs (mu, nu) in {V, Vbar}^2.
 3. The deformed quantum Serre relation F_{ij}(b_i, b_j) = q c_i b_j + eps_j((1-q)(1-q^-1) eps_i^2 - q c_i) k_i^2 k_j
    holds on V(y_1) x Vbar(y_2) x V(y_3), and its Cartan obstruction vanishes iff eps_j = 0 or eps_i^2 = eps_*^2.
 4. Boundary unitarity K(-theta)^T K(theta) = U_s I and matrix boundary crossing-unitarity
    k(theta) = Psi_s(mu) R_{Vbar Vbar}(2 theta) k(-theta) hold to machine precision.
 5. Generates lecture-notes/half-line/figures/l5_f1_a2_serre.svg.
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code" / "figures"))
from style import setup, save, report, all_ok, plt, BLUE, RED, GREEN, ORANGE, PURPLE, GREY, LIGHTGREY

OUT_DIR = ROOT / "lecture-notes" / "half-line" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def a2_rep(q, y, conj=False):
    """3-dimensional evaluation representation V(y) (conj=False) or Vbar(y) (conj=True) of U_q(a_2^(1))."""
    N = 3
    rep = {"dim": N}
    for j in range(N):
        a, b = (2, 0) if j == 0 else (j - 1, j)
        E = np.zeros((N, N), dtype=complex)
        E[a, b] = 1.0
        F = E.T
        k_diag = np.ones(N, dtype=complex)
        k_diag[a] = q
        k_diag[b] = 1.0 / q
        if not conj:
            rep[("e", j)] = y * E
            rep[("f", j)] = F / y
            rep[("k", j)] = np.diag(k_diag)
        else:
            rep[("e", j)] = y * F
            rep[("f", j)] = E / y
            rep[("k", j)] = np.diag(1.0 / k_diag)
    return rep


def cop(A, B, gen, j):
    """Part VI coproduct Delta(e) = e x 1 + k x e, Delta(f) = f x k^-1 + 1 x f, Delta(k) = k x k."""
    IA = np.eye(A["dim"], dtype=complex)
    IB = np.eye(B["dim"], dtype=complex)
    if gen == "e":
        return np.kron(A[("e", j)], IB) + np.kron(A[("k", j)], B[("e", j)])
    if gen == "f":
        kB_inv = np.diag(1.0 / np.diag(B[("k", j)]))
        return np.kron(A[("f", j)], kB_inv) + np.kron(IA, B[("f", j)])
    if gen == "k":
        return np.kron(A[("k", j)], B[("k", j)])
    raise ValueError(gen)


def tensor_rep(A, B):
    """Build tensor product representation A x B."""
    AB = {"dim": A["dim"] * B["dim"]}
    for j in range(3):
        for g in ("e", "f", "k"):
            AB[(g, j)] = cop(A, B, g, j)
    return AB


def coideal_b(rep, j, c_j, eps_j):
    """Coideal generator b_j = e_j + c_j f_j k_j + eps_j k_j on rep."""
    return rep[("e", j)] + c_j * (rep[("f", j)] @ rep[("k", j)]) + eps_j * rep[("k", j)]


def k_intertwiner_svd(q, y, eps_vec, conj_in=False):
    """Return smallest singular value and best nullvector reshaped to 3x3 for K: V_mu(y) -> V_mubar(1/y)."""
    Vin = a2_rep(q, y, conj=conj_in)
    Vout = a2_rep(q, 1.0 / y, conj=(not conj_in))
    I3 = np.eye(3, dtype=complex)
    rows = []
    for j in range(3):
        b_in = coideal_b(Vin, j, 1.0 / q, eps_vec[j])
        b_out = coideal_b(Vout, j, 1.0 / q, eps_vec[j])
        # K b_in - b_out K = 0 in vec(K) row-major: (I x b_in^T - b_out x I) vec(K) = 0
        M_j = np.kron(I3, b_in.T) - np.kron(b_out, I3)
        rows.append(M_j)
    M = np.vstack(rows)
    u, s, vh = np.linalg.svd(M)
    K_best = vh[-1].conj().reshape(3, 3)
    return s[-1], K_best


def k_a2_closed(q, y, s_sign):
    """Closed-form Toeplitz reflection matrix (eq-an-k-matrix) for n=2 (N=3), s = +-1."""
    N = 3
    w = np.sqrt(-q + 0j)
    es = 1.0 / np.sqrt((1.0 - q) * (1.0 - 1.0 / q))
    X = w * y
    K = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            if i == j:
                K[i, j] = (w * X ** (-N / 2.0) - (s_sign / w) * X ** (N / 2.0)) / (es * (1.0 / q - q))
            elif j > i:
                K[i, j] = X ** (j - i - N / 2.0)
            else:
                K[i, j] = s_sign * X ** (j - i + N / 2.0)
    return K


def rcheck_intertwiner(A, B):
    """Solve R_check Delta_{A,B}(x) = Delta_{B,A}(x) R_check normalized to 1 on highest weight (0,0)."""
    d = A["dim"] * B["dim"]
    Id = np.eye(d, dtype=complex)
    rows = []
    for j in range(3):
        for g in ("e", "f", "k"):
            X_in = cop(A, B, g, j)
            X_out = cop(B, A, g, j)
            rows.append(np.kron(Id, X_in.T) - np.kron(X_out, Id))
    M = np.vstack(rows)
    _, s, vh = np.linalg.svd(M)
    R = vh[-1].conj().reshape(d, d)
    return R / R[0, 0]


def main():
    setup()
    omega = np.sqrt(2.0) - 1.0  # generic irrational omega ~ 0.41421356
    q = -np.exp(-1j * np.pi * omega)
    eps_star = 1.0 / (2.0 * np.cos(0.5 * np.pi * omega))
    y0 = np.exp(omega * (0.43 + 0.17j))

    # 1. Check existence & closed form for K: V(y) -> Vbar(1/y)
    s_zero, K_zero = k_intertwiner_svd(q, y0, [0.0, 0.0, 0.0])
    K_zero = K_zero / K_zero[0, 0]
    err_zero = np.linalg.norm(K_zero - np.eye(3))

    s_plus, K_plus_num = k_intertwiner_svd(q, y0, [eps_star, eps_star, eps_star])
    K_plus_ex = k_a2_closed(q, y0, +1)
    c_p = np.vdot(K_plus_ex.ravel(), K_plus_num.ravel()) / np.vdot(K_plus_ex.ravel(), K_plus_ex.ravel())
    err_plus = np.linalg.norm(K_plus_num - c_p * K_plus_ex) / np.linalg.norm(K_plus_num)

    s_minus, K_minus_num = k_intertwiner_svd(q, y0, [-eps_star, eps_star, eps_star])
    K_minus_ex = k_a2_closed(q, y0, -1)
    c_m = np.vdot(K_minus_ex.ravel(), K_minus_num.ravel()) / np.vdot(K_minus_ex.ravel(), K_minus_ex.ravel())
    err_minus = np.linalg.norm(K_minus_num - c_m * K_minus_ex) / np.linalg.norm(K_minus_num)

    s_generic, _ = k_intertwiner_svd(q, y0, [0.6 * eps_star, 0.6 * eps_star, 0.6 * eps_star])
    s_mixed, _ = k_intertwiner_svd(q, y0, [0.0, eps_star, eps_star])

    ok1 = (
        max(s_zero, s_plus, s_minus, err_zero, err_plus, err_minus) < 1e-14
        and min(s_generic, s_mixed) > 0.1
    )
    report(
        f"a_2^(1) intertwiner K exists iff eps=0 or eps_j^2=eps_*^2 and matches Toeplitz closed form "
        f"(err={max(err_zero, err_plus, err_minus):.2e}, generic s_min={s_generic:.2f})",
        ok1,
    )

    # 2. Check soliton-conjugating reflection equation for all 4 pairs (mu, nu) in {V, Vbar}^2
    y1 = np.exp(omega * (0.51 - 0.12j))
    y2 = np.exp(omega * (0.23 + 0.09j))
    I3 = np.eye(3, dtype=complex)
    re_max_err = 0.0
    for s_val, eps_vec in [
        (0, [0.0, 0.0, 0.0]),
        (+1, [eps_star, eps_star, eps_star]),
        (-1, [-eps_star, eps_star, eps_star]),
    ]:
        for mu_conj in (False, True):
            for nu_conj in (False, True):
                _, K_mu_1 = k_intertwiner_svd(q, y1, eps_vec, conj_in=mu_conj)
                _, K_nu_2 = k_intertwiner_svd(q, y2, eps_vec, conj_in=nu_conj)
                R_mu_nu = rcheck_intertwiner(
                    a2_rep(q, y1 / y2, conj=mu_conj), a2_rep(q, 1.0, conj=nu_conj)
                )
                R_nu_mubar = rcheck_intertwiner(
                    a2_rep(q, y1 * y2, conj=nu_conj), a2_rep(q, 1.0, conj=(not mu_conj))
                )
                R_mu_nubar = rcheck_intertwiner(
                    a2_rep(q, y1 * y2, conj=mu_conj), a2_rep(q, 1.0, conj=(not nu_conj))
                )
                R_nubar_mubar = rcheck_intertwiner(
                    a2_rep(q, y1 / y2, conj=(not nu_conj)), a2_rep(q, 1.0, conj=(not mu_conj))
                )
                K1_2 = np.kron(I3, K_mu_1)
                K2_2 = np.kron(I3, K_nu_2)
                lhs = R_nubar_mubar @ K1_2 @ R_mu_nubar @ K2_2
                rhs = K2_2 @ R_nu_mubar @ K1_2 @ R_mu_nu
                # Check proportionality (and equality when mu == nu)
                if mu_conj == nu_conj:
                    err = np.linalg.norm(lhs - rhs) / np.linalg.norm(lhs)
                else:
                    c_re = np.vdot(rhs.ravel(), lhs.ravel()) / np.vdot(rhs.ravel(), rhs.ravel())
                    err = np.linalg.norm(lhs - c_re * rhs) / np.linalg.norm(lhs)
                re_max_err = max(re_max_err, err)

    report(
        f"Soliton-conjugating reflection equation holds for all (mu,nu) in {{V,Vbar}}^2 and all 3 K-matrices "
        f"(max err={re_max_err:.2e})",
        re_max_err < 1e-12,
    )

    # 3. Check deformed quantum Serre relation F_{ij}(b_i, b_j) on V(y1) x Vbar(y2) x V(y3)
    W = tensor_rep(tensor_rep(a2_rep(q, 0.8 + 0.3j, False), a2_rep(q, 1.2 - 0.2j, True)), a2_rep(q, 0.6 - 0.4j, False))
    two_q = q + 1.0 / q
    serre_err = 0.0
    for i, j in [(0, 1), (1, 2), (2, 0), (1, 0)]:
        c_i, c_j = 1.0 / q, 1.0 / q
        for ei, ej in [(0.73 - 0.41j, -0.52 + 0.88j), (eps_star, -eps_star)]:
            bi = coideal_b(W, i, c_i, ei)
            bj = coideal_b(W, j, c_j, ej)
            ki = W[("k", i)]
            kj = W[("k", j)]
            F_ij = bi @ bi @ bj - two_q * (bi @ bj @ bi) + bj @ bi @ bi
            rhs = q * c_i * bj + ej * ((1.0 - q) * (1.0 - 1.0 / q) * (ei ** 2) - q * c_i) * (ki @ ki @ kj)
            serre_err = max(serre_err, np.linalg.norm(F_ij - rhs))
    # Verify obstruction vanishes at ei = +-eps_star when c_i = 1/q
    bi_star = coideal_b(W, 0, 1.0 / q, eps_star)
    bj_star = coideal_b(W, 1, 1.0 / q, -eps_star)
    F_star_err = np.linalg.norm(
        bi_star @ bi_star @ bj_star - two_q * (bi_star @ bj_star @ bi_star) + bj_star @ bi_star @ bi_star - bj_star
    )
    report(
        f"Deformed Serre relation F_ij(b_i,b_j) == q c_i b_j + Obstruction * k_i^2 k_j "
        f"(err={serre_err:.2e}, at eps_* err={F_star_err:.2e})",
        max(serre_err, F_star_err) < 1e-13,
    )

    # 4. Check boundary unitarity K(-theta)^T K(theta) = U_s I and matrix crossing-unitarity
    theta_test = 0.37 - 0.11j
    T_param = 1.5 * omega  # h*omega/2 with h=3
    mu_test = -1j * T_param * theta_test / np.pi
    unit_cross_err = 0.0
    R_WbWb = rcheck_intertwiner(
        a2_rep(q, np.exp(2.0 * omega * theta_test), conj=True),
        a2_rep(q, 1.0, conj=True),
    )
    for s_sign in (+1, -1):
        K_th = k_a2_closed(q, np.exp(omega * theta_test), s_sign)
        K_mth = k_a2_closed(q, np.exp(-omega * theta_test), s_sign)
        U_s = (s_sign * np.cosh(2.0 * T_param * theta_test) - np.cos(np.pi * T_param)) / (
            2.0 * (np.sin(0.5 * np.pi * omega) ** 2)
        )
        unit_cross_err = max(unit_cross_err, np.linalg.norm(K_mth.T @ K_th - U_s * I3))

        # Matrix crossing-unitarity: k(theta) = Psi_s(mu) R_{Vbar Vbar}(2 theta) k(-theta)
        ym = np.exp(omega * (0.5j * np.pi - theta_test))
        yp = np.exp(omega * (0.5j * np.pi + theta_test))
        k_th = k_a2_closed(q, ym, s_sign).T.reshape(9)
        k_mth = k_a2_closed(q, yp, s_sign).T.reshape(9)
        if s_sign == +1:
            Psi_s = np.sin(np.pi * (0.5 * omega - mu_test)) / np.sin(np.pi * (0.5 * omega + mu_test))
        else:
            Psi_s = np.cos(np.pi * (mu_test - 0.5 * omega)) / np.cos(np.pi * (mu_test + 0.5 * omega))
        unit_cross_err = max(unit_cross_err, np.linalg.norm(k_th - Psi_s * (R_WbWb @ k_mth)))

    report(
        f"Boundary unitarity K(-th)^T K(th) = U_s I and crossing k(th) = Psi_s R(2th) k(-th) "
        f"(max err={unit_cross_err:.2e})",
        unit_cross_err < 1e-13,
    )

    # 5. Generate two-panel figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.4))

    # Left panel: Smallest singular value of 27x9 intertwiner system & Serre obstruction vs eps / eps_*
    r_grid = np.linspace(-1.6, 1.6, 121)
    s_min_vals = []
    obst_vals = []
    for r in r_grid:
        ev = r * eps_star
        sm, _ = k_intertwiner_svd(q, y0, [ev, ev, ev])
        s_min_vals.append(sm)
        # Magnitude of Cartan obstruction coefficient |eps_j ((1-q)(1-q^-1) eps_i^2 - 1)| / eps_*
        obst_vals.append(abs(r * (r ** 2 - 1.0)))

    ax1.plot(r_grid, s_min_vals, color=BLUE, lw=2.2, label=r"Smallest singular value $\sigma_{\min}(\mathcal{M}(\hat\epsilon))$")
    ax1.plot(r_grid, np.array(obst_vals) * 0.45, color=RED, lw=1.8, ls="--",
             label=r"Scaled Serre obstruction $\propto |\hat\epsilon(\hat\epsilon^2/\epsilon_*^2 - 1)|$")

    for r_pt, lbl, col in [
        (-1.0, r"$-\epsilon_*$ ($s=-1$)", PURPLE),
        (0.0, r"$0$ ($K=\mathbf{1}$)", GREEN),
        (1.0, r"$+\epsilon_*$ ($s=+1$)", ORANGE),
    ]:
        ax1.plot(r_pt, 0.0, "o", color=col, markersize=8, zorder=5)
        ax1.annotate(
            lbl,
            xy=(r_pt, 0.03),
            xytext=(r_pt, 0.52 if r_pt != 0 else 0.65),
            ha="center",
            fontsize=8.5,
            color=col,
            fontweight="bold",
            arrowprops=dict(arrowstyle="->", color=col, lw=1.1),
        )

    ax1.set_xlabel(r"Normalized coideal parameter $\hat\epsilon / \epsilon_*$")
    ax1.set_ylabel(r"Singular value / Serre obstruction")
    ax1.set_title(r"$a_2^{(1)}$ Quantum Rigidity: Intertwiner vs. Serre Obstruction")
    ax1.set_xlim(-1.6, 1.6)
    ax1.set_ylim(-0.05, 1.15)
    ax1.legend(loc="upper center", fontsize=8, framealpha=0.92)
    ax1.grid(True, color=LIGHTGREY, lw=0.6)

    # Right panel: Why first-order perturbation theory around Neumann missed eps_*
    lam_sg = 3.4
    eps_star_sg = 1.0 / (2.0 * np.cos(0.5 * np.pi * lam_sg))
    mu_grid = np.linspace(-0.25, 1.0, 250)
    # Exact UV-IR curve: eps_hat(mu) / eps_* = -cos((lam+1) arccos(mu)) / (sin(pi*lam) * eps_*)
    eps_exact_ratio = -np.cos((lam_sg + 1.0) * np.arccos(np.clip(mu_grid, -1.0, 1.0))) / (
        np.sin(np.pi * lam_sg) * eps_star_sg
    )
    # First-order tangent around Neumann (mu = 0): ratio = 1 - (lam + 1) * cot(pi*lam/2) * mu
    slope_neumann = -(lam_sg + 1.0) / np.tan(0.5 * np.pi * lam_sg)
    eps_tangent_true = 1.0 + slope_neumann * mu_grid
    eps_tangent_naive = slope_neumann * mu_grid

    ax2.plot(mu_grid, eps_exact_ratio, color=BLUE, lw=2.3,
             label=r"Exact UV–IR curve $\hat\epsilon(\mu)/\epsilon_*$ ($\lambda=3.4$)")
    ax2.plot(mu_grid, eps_tangent_true, color=GREEN, lw=1.7, ls="--",
             label=r"True $O(\mu)$ tangent at Neumann: $\epsilon_* + \kappa\mu$")
    ax2.plot(mu_grid, eps_tangent_naive, color=RED, lw=1.5, ls=":",
             label=r"Naive $g=0$ reading: $0 + \kappa\mu$ (misses $\epsilon_*$)")

    mu_zero = np.cos(np.pi / (2.0 * (lam_sg + 1.0)))
    ax2.plot(0.0, 1.0, "o", color=ORANGE, markersize=8, zorder=5)
    ax2.annotate(
        r"Neumann ($\mu=0$)" + "\n" + r"$\hat\epsilon = +\epsilon_*$ ($s=+1$)",
        xy=(0.0, 1.0),
        xytext=(0.18, 1.25),
        fontsize=8.2,
        color=ORANGE,
        fontweight="bold",
        arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.1),
    )

    ax2.plot(mu_zero, 0.0, "s", color=GREEN, markersize=7.5, zorder=5)
    ax2.annotate(
        r"Uniform limit ($\mu \to 1$)" + "\n" + r"$\hat\epsilon = 0$ ($K=\mathbf{1}$)",
        xy=(mu_zero, 0.0),
        xytext=(0.42, -0.68),
        fontsize=8.2,
        color=GREEN,
        fontweight="bold",
        arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.1),
    )

    ax2.axhline(0.0, color=GREY, lw=0.8)
    ax2.axvline(0.0, color=GREY, lw=0.8)
    ax2.set_xlabel(r"Boundary coupling $\mu = \tilde\mu / \tilde\mu_{\mathrm{crit}} \propto C_j$")
    ax2.set_ylabel(r"Normalized coideal parameter $\hat\epsilon / \epsilon_*$")
    ax2.set_title(r"Neumann ($\hat\epsilon=\epsilon_*$) vs. Uniform ($\hat\epsilon=0$) Boundary")
    ax2.set_xlim(-0.22, 1.02)
    ax2.set_ylim(-1.1, 1.65)
    ax2.legend(loc="lower left", fontsize=7.8, framealpha=0.92)
    ax2.grid(True, color=LIGHTGREY, lw=0.6)

    fig.tight_layout()
    save(fig, "l5_f1_a2_serre", OUT_DIR)

    if not all_ok():
        sys.exit(1)


if __name__ == "__main__":
    main()
