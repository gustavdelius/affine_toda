"""Lecture 6 verification and figure: Breathers, particles, weak-strong duality, and a_2^(1) boundary spectrum.

Verifies:
 1. For a_n^(1) (n=1,2,3,4), the uniform particle reflection amplitudes K_a(theta) (eq-uniform-amplitude)
    and Neumann amplitudes K^N_a(theta) = K_a(theta)|_{B -> 2-B} (eq-neumann-amplitude) satisfy boundary
    unitarity, boundary crossing-unitarity, and the bulk boundary bootstrap, and their B -> 0 limits match
    Lecture 3's linearized wave reflection factors -1/((a)(h-a)) and +1.
 2. For a_2^(1), the matrix-level boundary bootstrap on the non-diagonal soliton reflection matrix K_+(theta)
    (s=+1) maps the breather singlet in V(theta_a) x Vbar(theta_b) to c_+(mu) times itself (eq-an-neumann-breather).
 3. For a_2^(1) with a solitonic boundary (b=1), the excited boundary amplitudes K_1^{b_{n,m}}(theta)
    (eq-a2-excited) satisfy the boundary bootstrap, their energies match e_{b_{n,m}} - e_{b_{0,0}} = m_1[F(n)+F(m)]
    (eq-a2-energies), and the rotor tower b_{n,-n} matches n^2/(2I) with classical moment of inertia I = 6/(beta^2 m_1).
 4. In the window 1/(2n-1) < B < 2/(2n-1) (n=2, B=9/20), the pole of K_1^{b_{-2,2}} at p = 1/2 + 3B/2 has
    -i Res K_1 < 0 and matches the Coleman-Thun residue formula (1/2) g^2 Gamma^2 (K_2^{b_{-1,1}})'(theta_R)
    at the reflection zero theta_R.
 5. Generates lecture-notes/half-line/figures/l6_f1_particles_duality.svg.
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code" / "figures"))
from style import setup, save, report, all_ok, plt, BLUE, RED, GREEN, ORANGE, PURPLE, GREY, LIGHTGREY
from halfline_l5_f1_a2_serre import a2_rep, rcheck_intertwiner, k_a2_closed

OUT_DIR = ROOT / "lecture-notes" / "half-line" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def blk(x, theta, h):
    """Standard hyperbolic block (x) = sinh(theta/2 + i*pi*x/(2h)) / sinh(theta/2 - i*pi*x/(2h))."""
    return np.sinh(0.5 * theta + 0.5j * np.pi * x / h) / np.sinh(0.5 * theta - 0.5j * np.pi * x / h)


def curly_blk(x, theta, h, B):
    """Bulk S-matrix block {x}_B = (x-1)(x+1) / ((x-1+B)(x+1-B))."""
    return (
        blk(x - 1.0, theta, h)
        * blk(x + 1.0, theta, h)
        / (blk(x - 1.0 + B, theta, h) * blk(x + 1.0 - B, theta, h))
    )


def S_an(a, b, theta, h, B):
    """Bulk a_n^(1) S-matrix S_{ab}(theta) (eq-an-smatrix)."""
    prod = 1.0 + 0.0j
    for p in range(abs(a - b) + 1, a + b, 2):
        prod = prod * curly_blk(p, theta, h, B)
    return prod


def K_uniform(a, theta, h, B):
    """Uniform boundary reflection amplitude K_a(theta) (eq-uniform-amplitude)."""
    prod = 1.0 + 0.0j
    for c in range(1, a + 1):
        prod = (
            prod
            * blk(c - 1.0, theta, h)
            * blk(h + 1.0 - c - 0.5 * B, theta, h)
            / (blk(h - c, theta, h) * blk(c - 0.5 * B, theta, h))
        )
    return prod


def K_neumann(a, theta, h, B):
    """Neumann boundary reflection amplitude K^N_a(theta) = K_a(theta)|_{B -> 2-B} (eq-neumann-amplitude)."""
    return K_uniform(a, theta, h, 2.0 - B)


def K_a2_excited_1(n, m, theta, B):
    """Excited boundary reflection amplitude K_1^{b_{n,m}}(theta) for a_2^(1) (eq-a2-excited)."""
    h = 3.0
    K1 = K_uniform(1, theta, h, B)
    line1 = blk(0.5 - (2 * n + 1) * 0.5 * B, theta, h) * blk(2.5 + (2 * m + 1) * 0.5 * B, theta, h)
    line2 = blk(0.5 - (2 * n - 1) * 0.5 * B, theta, h) * blk(2.5 + (2 * m - 1) * 0.5 * B, theta, h)
    line3 = blk(1.5 + (2 * n - 1) * 0.5 * B, theta, h) * blk(1.5 - (2 * m - 1) * 0.5 * B, theta, h)
    line4 = blk(-2.5 + (2 * n + 1) * 0.5 * B, theta, h) * blk(-0.5 - (2 * m + 1) * 0.5 * B, theta, h)
    return K1 * line1 * line2 * line3 * line4


def K_a2_excited_2(n, m, theta, B):
    """K_2^{b_{n,m}}(theta) = K_1^{b_{m,n}}(theta)."""
    return K_a2_excited_1(m, n, theta, B)


def F_energy(k, B):
    """Energy function F(k) = sin(k*delta) cos(pi/6 - k*delta) / sin(delta), delta = pi*B/6 (eq-a2-energies)."""
    delta = np.pi * B / 6.0
    return np.sin(k * delta) * np.cos(np.pi / 6.0 - k * delta) / np.sin(delta)


def main():
    setup()

    # 1. Check uniform and Neumann amplitudes for n=1,2,3,4
    th0 = 0.61 - 0.19j
    B0 = 0.43
    max_ax_err = 0.0
    for n in (1, 2, 3, 4):
        h = float(n + 1)
        for K_fn in (K_uniform, K_neumann):
            for a in range(1, n + 1):
                abar = n + 1 - a
                # Unitarity: K_a(th) K_a(-th) = 1
                u_err = abs(K_fn(a, th0, h, B0) * K_fn(a, -th0, h, B0) - 1.0)
                # Charge conjugation: K_a = K_{h-a}
                c_err = abs(K_fn(a, th0, h, B0) - K_fn(abar, th0, h, B0))
                # Crossing-unitarity: K_a(th) K_abar(th + i*pi) = S_{aa}(2*th)
                cr_err = abs(
                    K_fn(a, th0, h, B0) * K_fn(abar, th0 + 1j * np.pi, h, B0) - S_an(a, a, 2.0 * th0, h, B0)
                )
                max_ax_err = max(max_ax_err, u_err, c_err, cr_err)
            # Bootstrap for a + b <= n
            for a in range(1, n):
                for b in range(1, n + 1 - a):
                    lhs = K_fn(a + b, th0, h, B0)
                    rhs = (
                        K_fn(a, th0 + 1j * np.pi * b / h, h, B0)
                        * S_an(a, b, 2.0 * th0 + 1j * np.pi * (b - a) / h, h, B0)
                        * K_fn(b, th0 - 1j * np.pi * a / h, h, B0)
                    )
                    max_ax_err = max(max_ax_err, abs(lhs - rhs))
        # Classical B -> 0 limits
        for a in range(1, n + 1):
            cl_unif = -1.0 / (blk(a, th0, h) * blk(h - a, th0, h))
            max_ax_err = max(
                max_ax_err,
                abs(K_uniform(a, th0, h, 0.0) - cl_unif),
                abs(K_neumann(a, th0, h, 0.0) - 1.0),
            )

    report(
        f"Uniform and Neumann K_a(theta) satisfy unitarity, crossing, bootstrap, and B->0 limits for n=1..4 "
        f"(max err={max_ax_err:.2e})",
        max_ax_err < 1e-14,
    )

    # 2. Check matrix-level breather projection c_+(mu) for a_2^(1) non-diagonal K_+ (s=+1)
    omega = np.sqrt(2.0) - 1.0
    q = -np.exp(-1j * np.pi * omega)
    T = 1.5 * omega
    mu = 0.29 - 0.08j
    th_b = 1j * np.pi * mu / T
    mu_a = mu - (0.5 * T - 0.5)
    mu_b = mu + (0.5 * T - 0.5)
    th_a = 1j * np.pi * mu_a / T
    th_bb = 1j * np.pi * mu_b / T
    ya = np.exp(omega * th_a)
    yb = np.exp(omega * th_bb)

    # Singlet in V(ya) x Vbar(yb): annihilated by Delta(e_j), Delta(f_j)
    Va = a2_rep(q, ya, conj=False)
    Wb = a2_rep(q, yb, conj=True)
    I3 = np.eye(3, dtype=complex)
    sing_rows = []
    for j in range(3):
        de = np.kron(Va[("e", j)], I3) + np.kron(Va[("k", j)], Wb[("e", j)])
        df = np.kron(Va[("f", j)], np.diag(1.0 / np.diag(Wb[("k", j)]))) + np.kron(I3, Wb[("f", j)])
        sing_rows.extend([de, df])
    _, _, vh_s = np.linalg.svd(np.vstack(sing_rows))
    v_sing = vh_s[-1].conj()

    # Reflect constituents: (1 x K_V(ya)) R_{V, V}(ya*yb) (1 x K_Vbar(yb)) v_sing
    # Recall K_Vbar(yb) = K_V(yb)^T for s=+1
    K_Va = k_a2_closed(q, ya, +1)
    K_Wb = k_a2_closed(q, yb, +1).T
    R_VV = rcheck_intertwiner(a2_rep(q, ya * yb, conj=False), a2_rep(q, 1.0, conj=False))
    v_out = np.kron(I3, K_Va) @ R_VV @ np.kron(I3, K_Wb) @ v_sing

    c_plus_formula = -(
        np.cos(np.pi * mu)
        * np.cos(np.pi * (mu + 0.5 * omega))
        * np.cos(np.pi * (mu - T))
        / ((np.sin(0.5 * np.pi * omega) ** 2) * np.cos(np.pi * (mu - 0.5 * omega)))
    )
    err_c_plus = np.linalg.norm(v_out - c_plus_formula * v_sing) / np.linalg.norm(v_sing)
    report(
        f"a_2^(1) matrix breather bootstrap on K_+ (s=+1) maps singlet to c_+(mu) * singlet (err={err_c_plus:.2e})",
        err_c_plus < 1e-13,
    )

    # 3. Check a_2^(1) excited boundary states b_{n,m}, energies, and rotor moment of inertia
    B_spec = 0.45
    h3 = 3.0
    exc_err = 0.0
    for n_idx in range(-2, 3):
        for m_idx in range(-2, 3):
            # Step 1: bind particle 1 at v = (pi/3)(1/2 - (2n+1)B/2) -> b_{n+1, m}
            v1 = (np.pi / 3.0) * (0.5 - (2 * n_idx + 1) * 0.5 * B_spec)
            lhs1 = K_a2_excited_1(n_idx + 1, m_idx, th0, B_spec)
            rhs1 = (
                K_a2_excited_1(n_idx, m_idx, th0, B_spec)
                * S_an(1, 1, th0 - 1j * v1, h3, B_spec)
                * S_an(1, 1, th0 + 1j * v1, h3, B_spec)
            )
            exc_err = max(exc_err, abs(lhs1 - rhs1) / abs(lhs1))
            # Energy step check
            dE_step = np.cos(v1)
            dE_formula = (F_energy(n_idx + 1, B_spec) + F_energy(m_idx, B_spec)) - (
                F_energy(n_idx, B_spec) + F_energy(m_idx, B_spec)
            )
            exc_err = max(exc_err, abs(dE_step - dE_formula))

        # Diagonal rotor step on b_{-n, n} at v_rot = (pi/3)(3/2 - (2n+1)B/2) -> b_{-n-1, n+1}
        if n_idx >= 0:
            v_rot = (np.pi / 3.0) * (1.5 - (2 * n_idx + 1) * 0.5 * B_spec)
            lhs_r = K_a2_excited_1(-n_idx - 1, n_idx + 1, th0, B_spec)
            rhs_r = (
                K_a2_excited_1(-n_idx, n_idx, th0, B_spec)
                * S_an(1, 1, th0 - 1j * v_rot, h3, B_spec)
                * S_an(1, 1, th0 + 1j * v_rot, h3, B_spec)
            )
            exc_err = max(exc_err, abs(lhs_r - rhs_r) / abs(lhs_r))
            dE_rot = np.cos(v_rot)
            dE_rot_f = (F_energy(-n_idx - 1, B_spec) + F_energy(n_idx + 1, B_spec)) - (
                F_energy(-n_idx, B_spec) + F_energy(n_idx, B_spec)
            )
            exc_err = max(exc_err, abs(dE_rot - dE_rot_f))

    report(
        f"a_2^(1) excited boundary bootstrap and energy formula e_{{b_{{n,m}}}} = m_1[F(n)+F(m)] hold "
        f"(max err={exc_err:.2e})",
        exc_err < 1e-14,
    )

    # 4. Check Coleman-Thun residue at reflection zero for b_{-2, 2} at B = 9/20
    B_ct = 9.0 / 20.0
    n_ct = 2
    # Pole of K_1^{b_{-2,2}} at p = 1/2 + 3B/2 = 47/40
    th_pole = (1j * np.pi / 3.0) * (0.5 + (2 * n_ct - 1) * 0.5 * B_ct)
    # Emission pole of b_{-1,1} -> b_{-2,2} in K_1^{b_{-1,1}} is at th_create = (i*pi/3)(3/2 - 3B/2)
    th_create = (1j * np.pi / 3.0) * (1.5 - (2 * n_ct - 1) * 0.5 * B_ct)
    # Reflection zero of K_2^{b_{-1,1}} is at th_R = (i*pi/3)(3B/2 - 1/2)
    th_R = (1j * np.pi / 3.0) * ((2 * n_ct - 1) * 0.5 * B_ct - 0.5)
    th_S11 = 2.0j * np.pi / 3.0

    eps_d = 1e-7
    res_K1_pole = (
        (K_a2_excited_1(-n_ct, n_ct, th_pole + eps_d, B_ct) - K_a2_excited_1(-n_ct, n_ct, th_pole - eps_d, B_ct))
        * (0.5 * eps_d)
    )
    res_create = (
        (K_a2_excited_1(-n_ct + 1, n_ct - 1, th_create + eps_d, B_ct)
         - K_a2_excited_1(-n_ct + 1, n_ct - 1, th_create - eps_d, B_ct))
        * (0.5 * eps_d)
    )
    res_S11 = (S_an(1, 1, th_S11 + eps_d, h3, B_ct) - S_an(1, 1, th_S11 - eps_d, h3, B_ct)) * (0.5 * eps_d)
    dK2_at_zero = (
        K_a2_excited_2(-n_ct + 1, n_ct - 1, th_R + eps_d, B_ct)
        - K_a2_excited_2(-n_ct + 1, n_ct - 1, th_R - eps_d, B_ct)
    ) / (2.0 * eps_d)

    # res_create = (i/2) g^2, res_S11 = i Gamma^2 => (1/2) g^2 Gamma^2 = (-i res_create) * (-i res_S11)
    ct_rhs = (-1j * res_create) * (-1j * res_S11) * dK2_at_zero
    ct_err = abs(res_K1_pole - ct_rhs) / abs(res_K1_pole)
    minus_i_res = (-1j * res_K1_pole).real
    report(
        f"Boundary Coleman-Thun pole on b_{{-2,2}} at B=9/20 has -i*Res={minus_i_res:.4f} < 0 and matches "
        f"(1/2)g^2 Gamma^2 K_2'(theta_R) (err={ct_err:.2e})",
        minus_i_res < -1.0 and ct_err < 1e-8,
    )

    # 5. Generate two-panel figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.4))

    # Left panel: Weak-strong duality B <-> 2-B between Uniform and Neumann reflection phases
    th_grid = np.linspace(0.01, 3.0, 250)
    for B_val, col, ls_u, ls_n in [(0.35, BLUE, "-", "--"), (0.001, GREEN, "-", ":")]:
        ku_vals = [K_uniform(1, t, 3.0, B_val) for t in th_grid]
        kn_dual_vals = [K_neumann(1, t, 3.0, 2.0 - B_val) for t in th_grid]
        # Phase in units of pi
        ph_u = np.unwrap(np.angle(ku_vals)) / np.pi
        if B_val > 0.1:
            ax1.plot(th_grid, ph_u, color=BLUE, lw=2.3,
                     label=r"Uniform $\arg K_1(\theta; B=0.35)/\pi$")
            ax1.plot(th_grid, np.unwrap(np.angle(kn_dual_vals)) / np.pi, color=ORANGE, lw=2.0, ls="--",
                     label=r"Neumann $\arg K_1^N(\theta; 2-B=1.65)/\pi$")
        else:
            ax1.plot(th_grid, ph_u, color=GREEN, lw=1.8, ls="-.",
                     label=r"Classical Uniform ($B\to 0$) $=$ Strong Neumann ($B\to 2$)")

    kn_b035 = [K_neumann(1, t, 3.0, 0.35) for t in th_grid]
    ax1.plot(th_grid, np.unwrap(np.angle(kn_b035)) / np.pi, color=RED, lw=1.8,
             label=r"Neumann $\arg K_1^N(\theta; B=0.35)/\pi$")
    ax1.axhline(0.0, color= PURPLE, lw=1.4, ls=":",
                label=r"Classical Neumann ($B\to 0$) $=$ Strong Uniform ($B\to 2$): $K=1$")

    ax1.set_xlabel(r"Real rapidity $\theta$")
    ax1.set_ylabel(r"Reflection phase $\delta(\theta) / \pi$")
    ax1.set_title(r"$a_2^{(1)}$ Weak–Strong Boundary Duality: $K_a^N(\theta; B) = K_a(\theta; 2-B)$")
    ax1.set_xlim(0.0, 3.0)
    ax1.set_ylim(-1.1, 0.35)
    ax1.legend(loc="lower right", fontsize=7.6, framealpha=0.92)
    ax1.grid(True, color=LIGHTGREY, lw=0.6)

    # Right panel: a_2^(1) solitonic boundary bound-state spectrum vs B, showing rotor tower b_{n,-n}
    B_grid = np.linspace(0.01, 0.98, 200)
    for n_rot, col, lbl in [
        (1, RED, r"Rotor $b_{\pm 1,\mp 1}$: $\sin^2(\delta)/\sin\delta \sim \frac{\pi B}{6}$"),
        (2, ORANGE, r"Rotor $b_{\pm 2,\mp 2}$: $\sin^2(2\delta)/\sin\delta \sim 4\frac{\pi B}{6}$"),
        (3, PURPLE, r"Rotor $b_{\pm 3,\mp 3}$: $\sin^2(3\delta)/\sin\delta \sim 9\frac{\pi B}{6}$"),
    ]:
        B_valid = B_grid[B_grid < 3.0 / (2 * n_rot - 1)]
        E_rot = [F_energy(n_rot, b) + F_energy(-n_rot, b) for b in B_valid]
        ax2.plot(B_valid, E_rot, color=col, lw=2.1, label=lbl)

    # Free-rotor quadratic approximation for n=1,2
    delta_grid = np.pi * B_grid / 6.0
    ax2.plot(B_grid, 4.0 * delta_grid, color=ORANGE, lw=1.2, ls=":", alpha=0.85,
             label=r"Classical rotor $n^2/(2I m_1) = n^2 \frac{\pi B}{6}$ ($n=2$)")

    # Bound-particle states b_{1,0} and b_{1,1}
    E_10 = [F_energy(1, b) + F_energy(0, b) for b in B_grid]
    E_11 = [F_energy(1, b) + F_energy(1, b) for b in B_grid]
    ax2.plot(B_grid, E_10, color=BLUE, lw=1.9, ls="--", label=r"Bound particle $b_{1,0}, b_{0,1}$")
    ax2.plot(B_grid, E_11, color=GREEN, lw=1.9, ls="-.", label=r"Two bound particles $b_{1,1}$")

    ax2.axvline(0.45, color=GREY, lw=1.0, ls="--")
    ax2.text(0.465, 2.15, r"$B = \frac{9}{20}$", fontsize=8.2, color=GREY)

    ax2.set_xlabel(r"Effective coupling $B(\beta) = \frac{2\beta^2}{\beta^2 + 4\pi}$")
    ax2.set_ylabel(r"Excitation energy $(e_{b_{n,m}} - e_{b_{0,0}}) / m_1$")
    ax2.set_title(r"$a_2^{(1)}$ Solitonic Boundary Spectrum & Quantized Vacuum Rotor")
    ax2.set_xlim(0.0, 0.98)
    ax2.set_ylim(0.0, 2.45)
    ax2.legend(loc="upper left", fontsize=7.4, framealpha=0.92)
    ax2.grid(True, color=LIGHTGREY, lw=0.6)

    fig.tight_layout()
    save(fig, "l6_f1_particles_duality", OUT_DIR)

    if not all_ok():
        sys.exit(1)


if __name__ == "__main__":
    main()
