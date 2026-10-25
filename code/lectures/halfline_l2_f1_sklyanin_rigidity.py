#!/usr/bin/env python3
"""
Half-Line Lecture 2 figure & verification script:
  lecture-notes/half-line/figures/l2_f1_sklyanin_rigidity.svg

Checks (one PASS/FAIL line each):
  1. a_1^(1): the explicit 2x2 K(lambda) solves the node equations eq-boundary-k-equation for sample
     complex (C_0, C_1); lambda^4 det K = lambda^8 - 4 C0 C1 lambda^6 + (4C0^2+4C1^2-2) lambda^4
     - 4 C0 C1 lambda^2 + 1 (generically invertible); at C_0 = C_1 = 0, lambda = 1 the matrix vanishes.
  2. a_2^(1): the explicit 3x3 K(lambda) solves the node equations for all 8 sign patterns C_j = +-1.
  3. a_2^(1) rigidity on a grid of 125 patterns with C_j in {0, 1, -1, generic, 1.7}, including mixed
     patterns: a nonzero K exists exactly for C = 0 and for all C_j^2 = 1 (sampled, one lambda); the
     figure's orange slice (C, +1, -1) has no K at C = 0.
  4. Both explicit K(lambda) satisfy Sklyanin's classical reflection equation (eq-classical-reflection)
     with the a_n^(1) r-matrix; K_- = 1 satisfies the equation at the left end, and at a vacuum
     a_t(lambda) + a_t(1/lambda)^T = 0 (so K_- = 1 satisfies the compatibility condition there).
  5. Dirichlet: {K, H} = 0 for all Cartan H forces K = 0 in the defining representation of a_2^(1) and
     a_3^(1), and leaves the antidiagonal K for a_1^(1).
  6. d_4^(1) vector representation (8x8): with the normalization sqrt(2 n_j/alpha_j^2) a K exists for all
     32 sign patterns; without n_j (node equation with C_j/sqrt(n_j)) none exists for |C_j| = 1, and one
     exists when C_2^2 = 2.
  7. a_2^(1) coweight flips D = diag(+-1): D K D^-1 solves the node equations with C_j -> d_a d_b C_j at the
     same lambda; the flips preserve s = C_0 C_1 C_2 (two classes of four); K^-1 solves them with -C_j.
  8. Soliton-preserving wall (sigma = Hermitian conjugation), one coweight of a_2^(1): with
     Re phi(0) in (pi/beta) nu + (2 pi/beta) (coweight lattice) and Neumann on Im phi, K_+ = exp(i pi nu.H)
     satisfies K_+ a_t(lambda) + a_t(1/lambda)^dagger K_+ = 0, while K_+ = 1 does not.
Runtime: a few seconds.
"""
from pathlib import Path
import sys
import itertools
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
def rel_min_sv(gens, C_list, lam):
    """Smallest singular value of the node system relative to the largest."""
    d = gens[0][0].shape[0]
    I_d = np.eye(d, dtype=complex)
    M = np.vstack([np.kron(I_d, (-Cj * hj + lam * Ej - Fj / lam).T) - np.kron(Cj * hj + lam * Ej - Fj / lam, I_d)
                   for (hj, Ej, Fj), Cj in zip(gens, C_list)])
    sv = np.linalg.svd(M, compute_uv=False)
    return sv[-1] / sv[0]


def node_residual(K, gens, C_list, lam):
    return max(np.max(np.abs(K @ (-Cj * hj + lam * Ej - Fj / lam) - (Cj * hj + lam * Ej - Fj / lam) @ K))
               for (hj, Ej, Fj), Cj in zip(gens, C_list))


def d4_vector_generators():
    """so(8) vector representation, basis e_1..e_4, e_-1..e_-4 (indices 0..3, 4..7).
    E_{e_i-e_j} = e_{i,j} - e_{-j,-i}, E_{e_i+e_j} = e_{i,-j} - e_{j,-i}; simple roots
    alpha_0 = -e1-e2, alpha_1 = e1-e2, alpha_2 = e2-e3, alpha_3 = e3-e4, alpha_4 = e3+e4."""
    def unit(i, j):
        m = np.zeros((8, 8), dtype=complex); m[i, j] = 1.0; return m
    p = lambda i: i - 1
    m = lambda i: 3 + i
    Emin = lambda i, j: unit(p(i), p(j)) - unit(m(j), m(i))
    Eplus = lambda i, j: unit(p(i), m(j)) - unit(p(j), m(i))
    Es = [Eplus(1, 2).T, Emin(1, 2), Emin(2, 3), Emin(3, 4), Eplus(3, 4)]
    return [(E @ E.T - E.T @ E, E, E.T.copy()) for E in Es]


def run_checks():
    gens_a1 = a1_generators()
    gens_a2 = a2_generators()
    lam = 1.27 + 0.31j

    # Check 1: a_1^(1)
    C0, C1 = 0.63 - 0.41j, -1.19 + 0.28j
    K1 = K_a1_classical(lam, C0, C1)
    err_a1 = node_residual(K1, gens_a1, [C0, C1], lam)
    det_poly = (lam**8 - 4 * C0 * C1 * lam**6 + (4 * C0**2 + 4 * C1**2 - 2) * lam**4
                - 4 * C0 * C1 * lam**2 + 1)
    err_det = abs(lam**4 * np.linalg.det(K1) - det_poly)
    zero_at_1 = np.max(np.abs(K_a1_classical(1.0, 0.0, 0.0)))
    report(
        f"a_1^(1): explicit K(lam) solves the node equations for sample complex C_0,C_1 (err={err_a1:.1e}); "
        f"lam^4 det K is the stated polynomial (err={err_det:.1e}); K = 0 at C=0, lam=1",
        err_a1 < 1e-13 and err_det < 1e-11 and zero_at_1 == 0.0,
    )

    # Check 2: a_2^(1), all 8 sign patterns
    err_a2_signs, min_det_a2 = 0.0, 1e9
    for c in itertools.product((1.0, -1.0), repeat=3):
        K2 = K_a2_classical(lam, *c)
        min_det_a2 = min(min_det_a2, abs(np.linalg.det(K2)))
        err_a2_signs = max(err_a2_signs, node_residual(K2, gens_a2, list(c), lam))
    report(
        f"a_2^(1): explicit K(lam) solves the node equations for all 8 sign patterns C_j=+-1 "
        f"(err={err_a2_signs:.1e}, min|det|={min_det_a2:.2f})",
        err_a2_signs < 1e-13 and min_det_a2 > 1e-2,
    )

    # Check 3: rigidity scan including mixed patterns, and the orange slice at C = 0
    vals = [0.0, 1.0, -1.0, 0.43 - 0.21j, 1.7]
    allowed, forbidden = [], []
    for c in itertools.product(vals, repeat=3):
        sv = rel_min_sv(gens_a2, list(c), lam)
        good = all(x == 0 for x in c) or all(abs(x**2 - 1) < 1e-12 for x in c)
        (allowed if good else forbidden).append(sv)
    sv_orange0 = rel_min_sv(gens_a2, [0.0, 1.0, -1.0], 1.25 + 0.2j)
    report(
        f"a_2^(1): on 125 sampled patterns (incl. mixed) a K exists only for C=0 or all C_j^2=1 "
        f"(max allowed sv {max(allowed):.1e}, min forbidden {min(forbidden):.3f}); slice (C,+1,-1) has none "
        f"at C=0 (sv {sv_orange0:.3f})",
        max(allowed) < 1e-12 and min(forbidden) > 1e-3 and sv_orange0 > 1e-3,
    )

    # Check 4: classical reflection equations, K_- = 1, vacuum condition
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
        r_iprod_t1 = partial_transpose(r_matrix_an(gens, 1.0 / (lam1 * lam2)), d, 1)
        lhs_cre = (
            K_l1 @ K_l2 @ r_12
            + r_21_tt @ K_l1 @ K_l2
            + K_l1 @ r_prod_t2 @ K_l2
            + K_l2 @ r_iprod_t1 @ K_l1
        )
        err_cre = np.max(np.abs(lhs_cre))
        # K_- = 1 in the equation at the left end: r^{t1 t2}(l1/l2) + r(l2/l1) + r^{t1}(l1 l2) + r^{t2}(1/(l1 l2))
        lhs_minus = (partial_transpose(r_matrix_an(gens, lam1 / lam2), d, 3) + r_matrix_an(gens, lam2 / lam1)
                     + partial_transpose(r_matrix_an(gens, lam1 * lam2), d, 1)
                     + partial_transpose(r_matrix_an(gens, 1.0 / (lam1 * lam2)), d, 2))
        err_minus = np.max(np.abs(lhs_minus))
        at_vac = lambda z, gens=gens: 0.5 * sum((z * Ej - Fj / z) for hj, Ej, Fj in gens)
        err_vac = np.max(np.abs(at_vac(lam).T + at_vac(1 / lam)))
        report(
            f"{name}: classical reflection equation for K (err={err_cre:.1e}); K_- = 1 solves the left-end "
            f"equation (err={err_minus:.1e}); vacuum a_t(lam)+a_t(1/lam)^T = 0 (err={err_vac:.1e})",
            err_cre < 1e-12 and err_minus < 1e-12 and err_vac < 1e-14,
        )

    # Check 5: Dirichlet anticommutator {K, H} = 0
    dims = {}
    for N in (2, 3, 4):
        Hs = [np.diag([1.0 if a == i else (-1.0 if a == i + 1 else 0.0) for a in range(N)]) for i in range(N - 1)]
        Id = np.eye(N)
        M = np.vstack([np.kron(Id, H.T) + np.kron(H, Id) for H in Hs])
        dims[N] = M.shape[1] - np.linalg.matrix_rank(M)
    report(
        f"Dirichlet: {{K,H}} = 0 for all Cartan H leaves dim {dims[3]} (a_2), {dims[4]} (a_3), "
        f"{dims[2]} (a_1, antidiagonal K)",
        dims[3] == 0 and dims[4] == 0 and dims[2] == 2,
    )

    # Check 6: d_4^(1) normalization
    g8 = d4_vector_generators()
    nj = [1, 1, 2, 1, 1]
    lam8 = 1.317 + 0.07j
    signs = list(itertools.product((1.0, -1.0), repeat=5))
    sv_norm = max(rel_min_sv(g8, list(sg), lam8) for sg in signs)
    sv_alt = min(rel_min_sv(g8, [c / np.sqrt(n) for c, n in zip(sg, nj)], lam8) for sg in signs)
    sv_alt2 = max(rel_min_sv(g8, [c / np.sqrt(n) for c, n in zip((sg[0], sg[1], sg[2] * np.sqrt(2), sg[3], sg[4]), nj)],
                             lam8) for sg in signs)
    report(
        f"d_4^(1) vector: K exists for all 32 patterns C_j=+-1 (sv {sv_norm:.1e}); without n_j none for "
        f"|C_j|=1 (min sv {sv_alt:.3f}); with C_2^2 = 2 it exists (sv {sv_alt2:.1e})",
        sv_norm < 1e-10 and sv_alt > 1e-3 and sv_alt2 < 1e-10,
    )

    # Check 7: coweight flips and inversion for a_2^(1)
    links = {0: (2, 0), 1: (0, 1), 2: (1, 2)}
    pats = list(itertools.product((1.0, -1.0), repeat=3))
    err_flip = 0.0
    tsets = set()
    for dvec in itertools.product((1.0, -1.0), repeat=3):
        D = np.diag(dvec)
        t = [dvec[links[j][0]] * dvec[links[j][1]] for j in range(3)]
        tsets.add(tuple(t))
        for c in pats:
            err_flip = max(err_flip, node_residual(D @ K_a2_classical(lam, *c) @ D, gens_a2,
                                                   [t[j] * c[j] for j in range(3)], lam))
    err_inv = max(node_residual(np.linalg.inv(K_a2_classical(lam, *c)), gens_a2, [-x for x in c], lam)
                  for c in pats)
    classes = {frozenset(tuple(tt[j] * c[j] for j in range(3)) for tt in tsets) for c in pats}
    report(
        f"a_2^(1): D K D^-1 solves the node equations with C_j -> d_a d_b C_j at the same lambda "
        f"(err={err_flip:.1e}); flips preserve s, {len(classes)} classes of 4; K^-1 solves them with -C_j "
        f"(err={err_inv:.1e})",
        err_flip < 1e-12 and err_inv < 1e-12 and len(classes) == 2
        and all(np.prod(tt) == 1 for tt in tsets) and all(len(cl) == 4 for cl in classes),
    )

    # Check 8: soliton-preserving wall with K_+ = exp(i pi nu.H), one coweight of a_2^(1)
    beta, mm = 0.9, 1.3
    roots = [np.array([-1.0, 0.0, 1.0]), np.array([1.0, -1.0, 0.0]), np.array([0.0, 1.0, -1.0])]
    nu = np.array([2.0, -1.0, -1.0]) / 3.0          # fundamental coweight: alpha_j . nu integer
    rng = np.random.default_rng(31)
    err_sp, err_one = 0.0, np.inf
    for trial in range(3):
        im_phi = rng.normal(size=3)
        im_phi -= im_phi.mean()
        re_phi = (np.pi / beta) * nu + (2 * np.pi / beta) * np.array([1.0, 0.0, -1.0]) * (trial - 1)
        phi = re_phi + 1j * im_phi
        lr = 1.37 + 0.1 * trial                      # real lambda
        w = [np.exp(0.5j * beta * np.dot(a, phi)) for a in roots]
        steps = lambda z: 0.5 * mm * sum(wj * (z * Ej - Fj / z) for wj, (hj, Ej, Fj) in zip(w, gens_a2))
        Kp = np.diag(np.exp(1j * np.pi * nu))         # exp(i pi nu.H)
        err_sp = max(err_sp, np.max(np.abs(Kp @ steps(lr) + steps(1 / lr).conj().T @ Kp)))
        err_one = min(err_one, np.max(np.abs(steps(lr) + steps(1 / lr).conj().T)))
    report(
        f"soliton-preserving wall: K_+ = exp(i pi nu.H) cancels the step terms at Re phi(0) = (pi/beta) nu "
        f"+ (2 pi/beta) coweights (err={err_sp:.1e}); K_+ = 1 does not ({err_one:.2f})",
        err_sp < 1e-12 and err_one > 1e-3,
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
        r"$\dot\tau(\lambda)=0 \Longleftarrow K_+(\lambda)a_t(\lambda)"
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
        label=r"$a_2^{(1)}$ equal couplings ($C_0=C_1=C_2=C$)",
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
        r"Neumann $C=0$ (green slice only)",
        xy=(0.0, 0.03),
        xytext=(0.0, 0.70),
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
