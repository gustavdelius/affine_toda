"""Lecture 6 verification and figure: breathers, particles, weak-strong duality, and the a_2^(1) boundary spectrum.

Checks (one PASS/FAIL line each):
 1. For a_n^(1) (n=1,2,3,4), the uniform amplitudes K_a(theta) (eq-uniform-amplitude) and the Neumann amplitudes
    K^N_a(theta) satisfy boundary unitarity, crossing-unitarity, charge conjugation and the bulk boundary bootstrap,
    and their B -> 0 limits are Lecture 3's linearized factors -1/((a)(h-a)) and +1.
 2. The printed Neumann product equals K_a(theta; 2-B) exactly (n=1..4).
 3. Continuation omega -> -2/B: the three imaginary-coupling blocks of K_B1 and of K_B1^(+) become K_1 and K^N_1;
    the middle label of K_B1 is -h-B/2, equal to h-B/2 by the period 2h of the blocks.
 4. For a_2^(1), the normalized R_{Vbar Vbar} fixes sum_i vbar_i x vbar_i, the crossed boundary vector of K = 1.
 5. For a_2^(1), the matrix-level breather bootstrap on K_+ (s=+1) maps the singlet in V(theta_a) x Vbar(theta_b)
    to c_+(mu) times itself.
 6. CDD factors: P_x = (T/2-x)_I (T/2+x)_I is unitary and crossing-symmetric for x = omega/2, (omega+1)/2 (n=2,3,4);
    P_x(mu_a) P_x(mu_b) continues to D(chi), chi = 2x/omega; D = 1 at B = 0; replacing P_{omega/2} by P_{(omega+1)/2}
    changes K^N_1 by a nonzero O(B^2) factor and breaks exact duality.
 7. Exercise 6.1: fusion angle u_11^2 = 2 pi/3 from the mass triangle; shift and doubling identities at h = 3.
 8. Exercise 6.2: O(beta^2) term of K^N_1 for a_2; Kim's [1/2][3/2] is self-dual, has the same O(B) term as K^N_1
    and a different O(B^2) term; the deleted second expression for K^N_1 is not equal to the first.
 9. For a_2^(1) with a solitonic boundary (b=1), the excited amplitudes K_1^{b_{n,m}} (eq-a2-excited) satisfy the boundary
    bootstrap and the energy steps of e_{b_{n,m}} - e_{b_{0,0}} = m_1[F(n)+F(m)] (eq-a2-energies), including the diagonal
    steps; F telescopes, F(-1)+F(1) = sin(delta), b_{-2,2} and b_{0,1} are degenerate at B = 1/2, and the creation poles
    of b_{n+1,m} and b_{-n-1,n+1} lie in the physical strip exactly for n+1 < 1/(2B)+1/2 and n+1 < 3/(2B)+1/2.
10. The classical moment of inertia of Lecture 3's a_2 solitonic vacuum, I = int |d phi/d chi|^2 dx, equals
    6/(beta^2 m_1) at several chi, and m_1 sin^2(n delta)/sin(delta) -> n^2/(2I) as beta -> 0.
11. In the window 1/(2n-1) < B < 2/(2n-1) (n=2, B=9/20), the pole of K_1^{b_{-2,2}} at p = 1/2 + 3B/2 has
    -i Res K_1 < 0 and matches the Coleman-Thun residue formula (1/2) g^2 Gamma^2 (K_2^{b_{-1,1}})'(theta_R)
    at the reflection zero theta_R, with relative tolerance 1e-8 (finite differences).
12. Generates lecture-notes/half-line/figures/l6_f1_particles_duality.svg.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad

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
    """Proposed Neumann amplitude K^N_a(theta), the printed product (eq-neumann-amplitude)."""
    prod = 1.0 + 0.0j
    for c in range(1, a + 1):
        prod = (
            prod
            * blk(c - 1.0, theta, h)
            * blk(h - c + 0.5 * B, theta, h)
            / (blk(h - c, theta, h) * blk(c - 1.0 + 0.5 * B, theta, h))
        )
    return prod


def blk_I(a, mu, h, omega):
    """Imaginary-coupling block (a)_I = sin(pi(mu+a)/(h omega)) / sin(pi(mu-a)/(h omega))."""
    return np.sin(np.pi * (mu + a) / (h * omega)) / np.sin(np.pi * (mu - a) / (h * omega))


def P_cdd(mu, x, h, omega):
    """CDD factor P_x = (T/2 - x)_I (T/2 + x)_I, T = h omega / 2."""
    T = 0.5 * h * omega
    return blk_I(0.5 * T - x, mu, h, omega) * blk_I(0.5 * T + x, mu, h, omega)


def D_cdd(chi, theta, h, B):
    """Continued breather image of P_x: D(chi) = (h+chi+B/2)(chi-B/2)(h-chi+B/2)(-chi-B/2)."""
    return (blk(h + chi + 0.5 * B, theta, h) * blk(chi - 0.5 * B, theta, h)
            * blk(h - chi + 0.5 * B, theta, h) * blk(-chi - 0.5 * B, theta, h))


def kim_block(x, theta, h, B):
    """Kim's block [x] = (x-1/2)(x+1/2) / ((x-1/2+B/2)(x+1/2-B/2))."""
    return (blk(x - 0.5, theta, h) * blk(x + 0.5, theta, h)
            / (blk(x - 0.5 + 0.5 * B, theta, h) * blk(x + 0.5 - 0.5 * B, theta, h)))


def a2_inertia(chi):
    """beta^2 m_1 I for Lecture 3's a_2 solitonic vacuum (units m_1 = beta = 1), I = int_{-inf}^0 |d phi/d chi|^2 dx,
    phi = -sum_j alpha_j log tau_j, tau_j = 1 + 2 cos(chi + 2 pi j/3)/(1+s) e^x + (1-s)/(1+s) e^{2x}, s = sin(pi/3)."""
    s = np.sqrt(3.0) / 2.0
    alpha = np.array([[-1.0, 0.0, 1.0], [1.0, -1.0, 0.0], [0.0, 1.0, -1.0]])  # alpha_0, alpha_1, alpha_2

    def integrand(x):
        u = np.exp(x)
        r = np.empty(3)
        for j in range(3):
            ph = chi + 2.0 * np.pi * j / 3.0
            tau = 1.0 + 2.0 * np.cos(ph) / (1.0 + s) * u + (1.0 - s) / (1.0 + s) * u * u
            r[j] = -2.0 * np.sin(ph) / (1.0 + s) * u / tau
        v = -(alpha.T @ r)
        return float(v @ v)

    pts = [-40.0, -10.0, -3.0, -1.0, -0.3, -0.05, 0.0]
    return sum(quad(integrand, pts[i], pts[i + 1], limit=200, epsabs=1e-13, epsrel=1e-12)[0]
               for i in range(len(pts) - 1))


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

    # 2. The printed Neumann product is the uniform amplitude at 2 - B
    dual_err = max(abs(K_neumann(a, th0, float(n + 1), B0) - K_uniform(a, th0, float(n + 1), 2.0 - B0))
                   for n in (1, 2, 3, 4) for a in range(1, n + 1))
    report(f"Printed Neumann product K^N_a(theta; B) equals K_a(theta; 2-B) for n=1..4 (max err={dual_err:.2e})",
           dual_err < 1e-14)

    # 3. Continuation omega -> -2/B of the breather amplitudes (a_I = (2a/omega))
    cont_err = 0.0
    for n in (2, 3, 4):
        h = float(n + 1)
        B1 = 1.0 / np.sqrt(5.0)
        om = -2.0 / B1
        Tc = 0.5 * h * om
        mu_c = -1j * Tc * th0 / np.pi
        KB1 = (blk_I(-0.5 * n * om, mu_c, h, om) * blk_I(0.5 - 0.5 * h * om, mu_c, h, om)
               * blk_I(-0.5 * om - 0.5, mu_c, h, om))
        KB1p = (blk_I(-0.5 * n * om, mu_c, h, om) * blk_I(0.5 * n * om - 0.5, mu_c, h, om)
                * blk_I(0.5, mu_c, h, om))
        cont_err = max(cont_err, abs(KB1 - K_uniform(1, th0, h, B1)), abs(KB1p - K_neumann(1, th0, h, B1)),
                       abs(blk_I(0.5 - 0.5 * h * om, mu_c, h, om) - blk(-h - 0.5 * B1, th0, h)),
                       abs(blk(-h - 0.5 * B1, th0, h) - blk(h - 0.5 * B1, th0, h)))
    report(f"Continuation omega -> -2/B: K_B1 -> K_1, K_B1^(+) -> K^N_1, middle label -h-B/2 = h-B/2 (period 2h) "
           f"(max err={cont_err:.2e})", cont_err < 1e-13)

    # 4. Crossed boundary vector of K = 1 is sum_i vbar_i x vbar_i, fixed by the normalized R_{Vbar Vbar}
    omega = np.sqrt(2.0) - 1.0
    q = -np.exp(-1j * np.pi * omega)
    R_bb = rcheck_intertwiner(a2_rep(q, np.exp(2.0 * omega * (0.37 - 0.11j)), conj=True), a2_rep(q, 1.0, conj=True))
    k_one = np.eye(3, dtype=complex).reshape(9)
    kv_err = np.linalg.norm(R_bb @ k_one - k_one)
    report(f"a_2^(1): normalized R_(Vbar Vbar) fixes the crossed boundary vector sum_i vbar_i x vbar_i of K = 1 "
           f"(err={kv_err:.2e})", kv_err < 1e-12)

    # 5. Check matrix-level breather projection c_+(mu) for a_2^(1) non-diagonal K_+ (s=+1)
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

    # 6. CDD factors P_x: unitarity, crossing symmetry, continued image D(chi), O(B^2) effect, loss of duality
    cdd_err = 0.0
    for n in (2, 3, 4):
        h = float(n + 1)
        omg = np.sqrt(2.0) + 1.0 / 7.0
        Tg = 0.5 * h * omg
        m1 = 0.23 - 0.17j
        for x in (0.5 * omg, 0.5 * (omg + 1.0)):
            cdd_err = max(cdd_err, abs(P_cdd(m1, x, h, omg) * P_cdd(-m1, x, h, omg) - 1.0),
                          abs(P_cdd(0.5 * Tg - m1, x, h, omg) - P_cdd(0.5 * Tg + m1, x, h, omg)))
        B2 = 1.0 / np.sqrt(7.0)
        om = -2.0 / B2
        Tb = 0.5 * h * om
        mu_c = -1j * Tb * (0.71 + 0.13j) / np.pi
        shift = 0.5 * Tb - 0.5
        for x in (0.5 * om, 0.5 * (om + 1.0)):
            image = P_cdd(mu_c - shift, x, h, om) * P_cdd(mu_c + shift, x, h, om)
            cdd_err = max(cdd_err, abs(image - D_cdd(2.0 * x / om, 0.71 + 0.13j, h, B2)))
        cdd_err = max(cdd_err, abs(D_cdd(1.0, th0, h, 0.0) - 1.0), abs(D_cdd(np.e / 5.0, th0, h, 0.0) - 1.0))
    th_r = 0.71
    coeffs = [abs(D_cdd(1.0 - 0.5 * Bv, th_r, 3.0, Bv) / D_cdd(1.0, th_r, 3.0, Bv) - 1.0) / Bv ** 2
              for Bv in (1e-2, 1e-3)]
    B3 = 3.0 / np.sqrt(17.0)
    K_alt = K_neumann(1, th0, 3.0, B3) * D_cdd(1.0 - 0.5 * B3, th0, 3.0, B3) / D_cdd(1.0, th0, 3.0, B3)
    alt_gap = abs(K_alt - K_uniform(1, th0, 3.0, 2.0 - B3))
    report(f"CDD: P_(omega/2), P_((omega+1)/2) unitary and crossing-symmetric, continue to D(chi), D=1 at B=0; "
           f"P_((omega+1)/2) changes K^N_1 at O(B^2) (coefficient {coeffs[1]:.4f}) and is not dual (gap {alt_gap:.2f})",
           cdd_err < 1e-12 and abs(coeffs[0] - coeffs[1]) < 0.01 * coeffs[1] and coeffs[1] > 0.1 and alt_gap > 1e-3)

    # 7. Exercise 6.1: fusion angle, shift and doubling identities at h = 3
    u11 = np.arccos(1.0 / 2.0 - 1.0)  # m_2^2 = 2 m_1^2 (1 + cos u) with m_2 = m_1 gives cos u = 1/2 - 1
    th1 = 0.71 + 0.05j
    B4 = np.sqrt(2.0) / 3.0
    ex1_err = abs(u11 - 2.0 * np.pi / 3.0) + abs(np.sqrt(2.0 * (1.0 + np.cos(np.pi / 3.0))) - np.sqrt(3.0))
    for x in (B4, 2.0 - B4, 2.0, np.e / 4.0):
        ex1_err = max(ex1_err, abs(blk(x, 2.0 * th1, 3.0) + blk(0.5 * x, th1, 3.0) * blk(0.5 * x + 3.0, th1, 3.0)),
                      abs(blk(x, th1 + 1j * np.pi / 3.0, 3.0) * blk(x, th1 - 1j * np.pi / 3.0, 3.0)
                          - blk(x + 1.0, th1, 3.0) * blk(x - 1.0, th1, 3.0)))
    K2_boot = (K_uniform(1, th1 + 1j * np.pi / 3.0, 3.0, B4) * S_an(1, 1, 2.0 * th1, 3.0, B4)
               * K_uniform(1, th1 - 1j * np.pi / 3.0, 3.0, B4))
    ex1_err = max(ex1_err, abs(K2_boot - K_uniform(1, th1, 3.0, B4)))
    report(f"Exercise 6.1: u_11^2 = 2pi/3 (pi/3 would give mass sqrt3 m_1), shift and doubling identities, "
           f"K_2 = K_1 (max err={ex1_err:.2e})", ex1_err < 1e-13)

    # 8. Exercise 6.2: one loop of K^N_1 and Kim's self-dual [1/2][3/2] (a_2, h = 3)
    th2 = 0.71
    KN1 = lambda Bv: K_neumann(1, th2, 3.0, Bv)
    Kkim = lambda Bv: kim_block(0.5, th2, 3.0, Bv) * kim_block(1.5, th2, 3.0, Bv)
    eps1, eps2 = 1e-5, 1e-3
    d1N = (KN1(eps1) - KN1(-eps1)) / (2 * eps1)
    d1K = (Kkim(eps1) - Kkim(-eps1)) / (2 * eps1)
    d2N = (KN1(eps2) - 2 * KN1(0.0) + KN1(-eps2)) / eps2 ** 2
    d2K = (Kkim(eps2) - 2 * Kkim(0.0) + Kkim(-eps2)) / eps2 ** 2
    one_loop = -0.25 * np.cos(th2 / 2j) / np.sin(3 * th2 / 2j)  # coefficient of beta^2, B = beta^2/(2 pi) + ...
    second_form = -blk(2.0 + 0.5 * 0.43, th2, 3.0) / (blk(2.0, th2, 3.0) * blk(3.0 - 0.5 * 0.43, th2, 3.0))
    ok8 = (abs(d1N / (2 * np.pi) - one_loop) < 1e-8 and abs(d1N - d1K) < 1e-8 and abs(d2N - d2K) > 0.1
           and abs(Kkim(0.43) - Kkim(2.0 - 0.43)) < 1e-14 and abs(second_form - KN1(0.43)) > 0.1)
    report(f"Exercise 6.2: O(beta^2) term of K^N_1; Kim's [1/2][3/2] self-dual, same O(B), different O(B^2) "
           f"(|second-derivative difference|={abs(d2N - d2K):.3f}); deleted second form differs", ok8)

    # 9. Check a_2^(1) excited boundary states b_{n,m} and their energies
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

    # Energy identities of Exercise 6.3 and the physical-strip range of the creation poles
    d_g = np.pi * 0.37 / 6.0
    for nk in range(1, 8):
        exc_err = max(exc_err, abs(sum(np.cos(np.pi / 6.0 - (2 * k + 1) * d_g) for k in range(nk))
                                   - F_energy(nk, 0.37)))
    exc_err = max(exc_err, abs(F_energy(-1, 0.37) + F_energy(1, 0.37) - np.sin(d_g)),
                  abs(F_energy(-2, 0.5) + F_energy(2, 0.5) - F_energy(1, 0.5) - F_energy(0.0, 0.5)))
    range_ok = all(((0.5 - (2 * nk + 1) * Bv / 2 > 0) == (nk + 1 < 1 / (2 * Bv) + 0.5))
                   and ((1.5 - (2 * nk + 1) * Bv / 2 > 0) == (nk + 1 < 3 / (2 * Bv) + 0.5))
                   for Bv in (0.23, 0.45, 0.61) for nk in range(6))
    report(
        f"a_2^(1) excited boundary bootstrap and energy formula e_{{b_{{n,m}}}} = m_1[F(n)+F(m)] hold; F telescopes, "
        f"F(-1)+F(1) = sin(delta), b_(-2,2) ~ b_(0,1) at B=1/2, strip ranges of creation poles (max err={exc_err:.2e})",
        exc_err < 1e-14 and range_ok,
    )

    # 10. Classical moment of inertia of the a_2 solitonic vacuum and the rotor-like limit of the diagonal tower
    inertia = [a2_inertia(c) for c in (0.0, 0.31, -0.77, 1.7, 2.6)]
    bet2 = 1e-4
    B_small = 2 * bet2 / (bet2 + 4 * np.pi)
    d_s = np.pi * B_small / 6.0
    rot_err = max(abs(np.sin(nq * d_s) ** 2 / np.sin(d_s) / (nq ** 2 * bet2 / 12.0) - 1.0) for nq in (1, 2, 3))
    report(f"a_2 solitonic vacuum: beta^2 m_1 I = {min(inertia):.10f}..{max(inertia):.10f} (= 6) at five chi; "
           f"m_1 sin^2(n delta)/sin(delta) -> n^2/(2I) as beta -> 0 (rel err {rot_err:.1e})",
           max(abs(v - 6.0) for v in inertia) < 1e-8 and rot_err < 1e-3)

    # 11. Check Coleman-Thun residue at reflection zero for b_{-2, 2} at B = 9/20
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
        f"(1/2)g^2 Gamma^2 K_2'(theta_R) within relative tolerance 1e-8 (err={ct_err:.2e})",
        minus_i_res < -1.0 and ct_err < 1e-8,
    )

    # 12. Generate two-panel figure
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
    ax1.set_title(r"$a_2^{(1)}$ proposed amplitudes: $K_a^N(\theta; B) = K_a(\theta; 2-B)$")
    ax1.set_xlim(0.0, 3.0)
    ax1.set_ylim(-1.1, 0.35)
    ax1.legend(loc="lower right", fontsize=7.6, framealpha=0.92)
    ax1.grid(True, color=LIGHTGREY, lw=0.6)

    # Right panel: a_2^(1) solitonic boundary bound-state spectrum vs B, with the rotor-like diagonal tower b_{n,-n}
    B_grid = np.linspace(0.01, 0.98, 200)
    for n_rot, col, lbl in [
        (1, RED, r"Diagonal $b_{\pm 1,\mp 1}$: $\sin^2\delta/\sin\delta \approx \frac{\pi B}{6}$"),
        (2, ORANGE, r"Diagonal $b_{\pm 2,\mp 2}$: $\sin^2 2\delta/\sin\delta \approx 4\frac{\pi B}{6}$"),
        (3, PURPLE, r"Diagonal $b_{\pm 3,\mp 3}$: $\sin^2 3\delta/\sin\delta \approx 9\frac{\pi B}{6}$"),
    ]:
        B_valid = B_grid[B_grid < 3.0 / (2 * n_rot - 1)]
        E_rot = [F_energy(n_rot, b) + F_energy(-n_rot, b) for b in B_valid]
        ax2.plot(B_valid, E_rot, color=col, lw=2.1, label=lbl)

    # Weak-coupling (rotor-like) quadratic limit for n=2
    delta_grid = np.pi * B_grid / 6.0
    ax2.plot(B_grid, 4.0 * delta_grid, color=ORANGE, lw=1.2, ls=":", alpha=0.85,
             label=r"Rotor-like limit $n^2 \frac{\pi B}{6} = n^2/(2I m_1)$ ($n=2$)")

    # Bound-particle states b_{1,0} and b_{1,1}
    E_10 = [F_energy(1, b) + F_energy(0, b) for b in B_grid]
    E_11 = [F_energy(1, b) + F_energy(1, b) for b in B_grid]
    ax2.plot(B_grid, E_10, color=BLUE, lw=1.9, ls="--", label=r"Bound particle $b_{1,0}, b_{0,1}$")
    ax2.plot(B_grid, E_11, color=GREEN, lw=1.9, ls="-.", label=r"Two bound particles $b_{1,1}$")

    ax2.axvline(0.45, color=GREY, lw=1.0, ls="--")
    ax2.text(0.465, 0.06, r"$B = \frac{9}{20}$", fontsize=8.2, color=GREY)

    ax2.set_xlabel(r"Effective coupling $B(\beta) = \frac{2\beta^2}{\beta^2 + 4\pi}$")
    ax2.set_ylabel(r"Excitation energy $(e_{b_{n,m}} - e_{b_{0,0}}) / m_1$")
    ax2.set_title(r"$a_2^{(1)}$ solitonic boundary spectrum, rotor-like diagonal tower")
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
