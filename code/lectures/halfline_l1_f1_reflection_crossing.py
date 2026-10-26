#!/usr/bin/env python3
"""
Half-Line Lecture 1 figure & verification script:
  lecture-notes/half-line/figures/l1_f1_reflection_crossing.svg

Checks (one PASS/FAIL line each):
  1. Braid/fixed-slot conversion: with Scheck = S_12 P = P S_21 and both K's in the second slot,
     the two sides of the braid reflection equation equal the two sides of the S/K equation term by
     term, for arbitrary 4x4 S and 2x2 K matrices; with Scheck = P S_12 they do not.
  2. A six-vertex example (diagonal K) satisfies both the S/K and the braid form (K in slot 2).
  3. Boundary crossing-unitarity in index form, K_abar^b(i pi/2 - t) = S_cd^ba(2t) K_cbar^d(i pi/2 + t):
     for the sine-Gordon coideal K (built here from the U_q(a_1^(1)) coproduct, physical q at
     lambda = sqrt 2) it holds up to a scalar with the lower indices cd, and fails with dc.
  4. Diagonal bootstrap channels for a_2^(1) particles (uniform boundary, B generic): the bulk-fusion
     bootstrap K_2 = K_1(+i pi/3) S_11(2 theta) K_1(-i pi/3) holds and fails with S_12; the excited-
     boundary dressing S_ba(th - iv) S_ba(th + iv) is unitary and crossing-compatible, while
     S_ba(th - iv) S_{b abar}(th + iv) is not.
  5. Three diagonal crossing-unitarity forms (centred, Corrigan-Delius, Delius-Gandenberger) for the
     Lee-Yang K_b at a sample rapidity.
  6. Lee-Yang block algebra of Blackboard Calculation 3: doubling identity (x)(2 theta) = -(x/2)/(6-x/2),
     S_LY(2 theta) = -(1)(2)/((4)(5)), shifted product -eps^2 (3)(5)/(2).
  7. Lee-Yang unitarity and bootstrap: holds with eps = +1, fails with eps = -1 (GZ's printed sign).
  8. Residues (1/eps scaling): Gamma^2 = -2 sqrt3 (Res S = i Gamma^2); Res_{i pi/2} K_(1) = (i/2) g^2 with
     g^2 = -4(2 sqrt3 - 3); Res_{i pi/6} K_(1) = i(3 - sqrt3), whose square is -Gamma^2 g^2/4.
  9. CDD factor of K_b: unitary, crossing-invariant and bootstrap-invariant; f_x = -(6x)(6-6x).
 10. Excited boundary: S(th - iv) S(th + iv) K_b = K_{b+2} = K_{4-b}; K_{-2} has no pole at i pi/6.
Runtime: a few seconds.
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
# Blocks
# ---------------------------------------------------------------------------
def blk(theta, x, h):
    """Block (x) = sinh(theta/2 + i pi x/(2h)) / sinh(theta/2 - i pi x/(2h))."""
    return np.sinh(0.5 * theta + 0.5j * np.pi * x / h) / np.sinh(
        0.5 * theta - 0.5j * np.pi * x / h
    )


def blk6(theta, x):
    """Lee-Yang blocks, h = 6 (period 12)."""
    return blk(theta, x, 6.0)


def S_LY(theta):
    """Bulk scaling Lee-Yang S-matrix: f_{2/3}(theta) = -(2)(4)."""
    return -blk6(theta, 2.0) * blk6(theta, 4.0)


def f_x(theta, x):
    return (np.sinh(theta) + 1j * np.sin(np.pi * x)) / (np.sinh(theta) - 1j * np.sin(np.pi * x))


def K_LY(theta, b, eps=1.0):
    """Dorey-Pocklington-Tateo-Watts boundary Lee-Yang reflection factor eps * K_b(theta)."""
    base = blk6(theta, 1.0) * blk6(theta, 3.0) / blk6(theta, 4.0)
    b_dep = (
        blk6(theta, 1.0 + b)
        * blk6(theta, 5.0 - b)
        / (blk6(theta, 1.0 - b) * blk6(theta, 5.0 + b))
    )
    return eps * base * b_dep


def residue(f, z0):
    """Simple-pole residue by symmetric limits; returns the value and its change under eps -> eps/10."""
    vals = [0.5 * e * (f(z0 + e) - f(z0 - e)) for e in (1e-5, 1e-6)]
    return vals[1], abs(vals[0] - vals[1])


# ---------------------------------------------------------------------------
# Sine-Gordon soliton multiplet: U_q(a_1^(1)) spin 1/2, principal gradation
# coproduct Delta(e) = e x 1 + k x e, Delta(f) = f x k^-1 + 1 x f, Delta(k) = k x k
# ---------------------------------------------------------------------------
def nullspace(M, tol=1e-9):
    _, sv, Vh = np.linalg.svd(M)
    rel = sv / sv[0]
    return [Vh[i].conj() for i in range(len(sv)) if rel[i] < tol]


def sg_rep(q, y):
    E = np.array([[0, 1], [0, 0]], complex)
    F = E.T.copy()
    K = np.diag([q, 1 / q]).astype(complex)
    return {("e", 1): y * E, ("f", 1): F / y, ("k", 1): K,
            ("e", 0): y * F, ("f", 0): E / y, ("k", 0): np.linalg.inv(K)}


def cop(A, B, kind, j):
    I2 = np.eye(2)
    if kind == "e":
        return np.kron(A[("e", j)], I2) + np.kron(A[("k", j)], B[("e", j)])
    if kind == "f":
        return np.kron(A[("f", j)], np.linalg.inv(B[("k", j)])) + np.kron(I2, B[("f", j)])
    return np.kron(A[("k", j)], B[("k", j)])


def intertwiner(pairs, d_in, d_out):
    """X (d_out x d_in) with X A = B X for all (A, B); row-major vec."""
    I_in, I_out = np.eye(d_in), np.eye(d_out)
    M = np.vstack([np.kron(I_out, A.T) - np.kron(B, I_in) for A, B in pairs])
    null = nullspace(M)
    assert len(null) == 1, len(null)
    return null[0].reshape(d_out, d_in)


def R_check_sg(q, z):
    """Rcheck: V(z) x V(1) -> V(1) x V(z), normalized to 1 on the highest-weight vectors."""
    A, B = sg_rep(q, z), sg_rep(q, 1.0)
    pairs = [(cop(A, B, kd, j), cop(B, A, kd, j)) for j in (0, 1) for kd in "efk"]
    R = intertwiner(pairs, 4, 4)
    return R / R[0, 0]


def K_sg(q, y, e0, e1):
    """Coideal K: V(y) -> V(1/y) with K b_j = b_j K, b_j = e_j + q^-1 f_j k_j + eps_j k_j."""
    Vi, Vo = sg_rep(q, y), sg_rep(q, 1 / y)
    b = lambda V, j, ep: V[("e", j)] + V[("f", j)] @ V[("k", j)] / q + ep * V[("k", j)]
    pairs = [(b(Vi, 0, e0), b(Vo, 0, e0)), (b(Vi, 1, e1), b(Vo, 1, e1))]
    return intertwiner(pairs, 2, 2)


def prop_res(u, v):
    c = np.vdot(v, u) / np.vdot(v, v)
    return np.linalg.norm(u - c * v) / np.linalg.norm(u)


# ---------------------------------------------------------------------------
# Verification checks
# ---------------------------------------------------------------------------
def run_checks():
    P = np.array(
        [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]],
        dtype=complex,
    )
    I2 = np.eye(2, dtype=complex)

    # Check 1: braid <-> S/K conversion on arbitrary matrices
    rng = np.random.default_rng(29)
    S_d = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    S_s = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    Ka = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    Kb = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    K1a, K2a, K2b = np.kron(Ka, I2), np.kron(I2, Ka), np.kron(I2, Kb)
    lhs_sk = S_d @ K1a @ (P @ S_s @ P) @ K2b
    rhs_sk = K2b @ S_s @ K1a @ (P @ S_d @ P)
    Sc = lambda S: S @ P
    err_ok = max(np.abs(Sc(S_d) @ K2a @ Sc(S_s) @ K2b - lhs_sk).max(),
                 np.abs(K2b @ Sc(S_s) @ K2a @ Sc(S_d) - rhs_sk).max())
    Sw = lambda S: P @ S
    err_bad = max(np.abs(Sw(S_d) @ K2a @ Sw(S_s) @ K2b - lhs_sk).max(),
                  np.abs(K2b @ Sw(S_s) @ K2a @ Sw(S_d) - rhs_sk).max())
    report(
        f"Scheck = S_12 P, K in slot 2: braid sides = S/K sides term by term for arbitrary matrices "
        f"(err={err_ok:.1e}); with Scheck = P S_12 they differ ({err_bad:.1f})",
        err_ok < 1e-12 and err_bad > 1e-3,
    )

    # Check 2: six-vertex example, both forms
    eta = 0.73 + 0.19j

    def R_6v(u):
        a = np.sinh(u + eta)
        b = np.sinh(u)
        c = np.sinh(eta)
        return np.array(
            [[a, 0, 0, 0], [0, b, c, 0], [0, c, b, 0], [0, 0, 0, a]],
            dtype=complex,
        )

    xi_B = 0.41 - 0.12j

    def K_diag(u):
        return np.diag([np.sinh(xi_B + u), np.sinh(xi_B - u)])

    u1, u2 = 0.37 + 0.11j, -0.29 + 0.23j
    K1 = np.kron(K_diag(u1), I2)
    K2_1, K2_2 = np.kron(I2, K_diag(u1)), np.kron(I2, K_diag(u2))
    S21 = lambda S: P @ S @ P
    err_sk = np.abs(R_6v(u1 - u2) @ K1 @ S21(R_6v(u1 + u2)) @ K2_2
                    - K2_2 @ R_6v(u1 + u2) @ K1 @ S21(R_6v(u1 - u2))).max()
    Rc = lambda u: R_6v(u) @ P
    err_braid = np.abs(Rc(u1 - u2) @ K2_1 @ Rc(u1 + u2) @ K2_2
                       - K2_2 @ Rc(u1 + u2) @ K2_1 @ Rc(u1 - u2)).max()
    report(
        f"six-vertex example: S/K form and braid form (K in slot 2) hold (max err={max(err_sk, err_braid):.1e})",
        max(err_sk, err_braid) < 1e-13,
    )

    # Check 3: crossing-unitarity index order, sine-Gordon coideal K
    lam = np.sqrt(2.0)
    q = -np.exp(-1j * np.pi * lam)          # physical q = -exp(-i pi lambda)
    res_ok, res_bad = [], []
    for e0, e1 in [(0.37 + 0.21j, -0.52 + 0.11j), (1.3 - 0.4j, 0.2 + 0.9j)]:
        def kvec(t):
            Km = K_sg(q, np.exp(lam * (0.5j * np.pi - t)), e0, e1)   # K[b, a] = K_a^b
            return np.array([Km[b, 1 - a] for a in range(2) for b in range(2)])  # K_abar^b
        for t in (0.29, 0.83):
            Rm = R_check_sg(q, np.exp(2 * lam * t))   # Rm[(a,b),(c,d)] = S_cd^{ba}
            res_ok.append(prop_res(kvec(t), Rm @ kvec(-t)))
            res_bad.append(prop_res(kvec(t), Rm @ P @ kvec(-t)))
    report(
        f"crossing-unitarity K_abar^b(i pi/2 - t) ~ S_cd^ba(2t) K_cbar^d(i pi/2 + t) for the sine-Gordon "
        f"coideal K (residual {max(res_ok):.1e}); with S_dc^ba it fails ({min(res_bad):.2f})",
        max(res_ok) < 1e-10 and min(res_bad) > 1e-2,
    )

    # Check 4: a_2 diagonal bootstrap channels (real coupling, h = 3)
    B = np.sqrt(2.0) / 4

    def curly(p, th):
        return (blk(th, p - 1, 3) * blk(th, p + 1, 3)
                / (blk(th, p - 1 + B, 3) * blk(th, p + 1 - B, 3)))

    Sab = {(1, 1): lambda t: curly(1, t), (2, 2): lambda t: curly(1, t),
           (1, 2): lambda t: curly(2, t), (2, 1): lambda t: curly(2, t)}

    def K1u(th):     # uniform boundary, a = 1 (= a = 2)
        return 1.0 / (blk(th, 2, 3) * blk(th, 1 - B / 2, 3) * blk(th, 3 + B / 2, 3))

    th = 0.71 + 0.13j
    bb_ok = abs(K1u(th) - K1u(th + 1j * np.pi / 3) * Sab[(1, 1)](2 * th) * K1u(th - 1j * np.pi / 3))
    bb_bad = abs(K1u(th) - K1u(th + 1j * np.pi / 3) * Sab[(1, 2)](2 * th) * K1u(th - 1j * np.pi / 3))
    v = 0.77 * np.pi / 3
    ok_d, bad_u, bad_c = 0.0, 0.0, 0.0
    for b in (1, 2):
        a = 1
        D = lambda t, b=b: Sab[(b, a)](t - 1j * v) * Sab[(b, a)](t + 1j * v)
        Db = lambda t, b=b: Sab[(3 - b, a)](t - 1j * v) * Sab[(3 - b, a)](t + 1j * v)
        W = lambda t, b=b: Sab[(b, a)](t - 1j * v) * Sab[(b, 3 - a)](t + 1j * v)
        Wb = lambda t, b=b: Sab[(3 - b, a)](t - 1j * v) * Sab[(3 - b, 3 - a)](t + 1j * v)
        ok_d = max(ok_d, abs(D(th) * D(-th) - 1), abs(D(0.5j * np.pi - th) - Db(0.5j * np.pi + th)))
        bad_u = max(bad_u, abs(W(th) * W(-th) - 1))
        bad_c = max(bad_c, abs(W(0.5j * np.pi - th) - Wb(0.5j * np.pi + th)))
    report(
        f"a_2 bootstrap channels: K_2 = K_1 S_11 K_1 holds ({bb_ok:.1e}), fails with S_12 ({bb_bad:.2f}); "
        f"S_ba S_ba dressing unitary and crossing-compatible ({ok_d:.1e}), S_ba S_(b abar) not "
        f"({bad_u:.2f}, {bad_c:.2f})",
        bb_ok < 1e-12 and bb_bad > 1e-2 and ok_d < 1e-12 and bad_u > 1e-2 and bad_c > 1e-2,
    )

    # Check 5: three diagonal crossing-unitarity forms for the Lee-Yang K_b
    th = 0.43 - 0.17j
    b_val = -0.37
    ip2 = 0.5j * np.pi
    err_book = abs(K_LY(ip2 - th, b_val) - S_LY(2.0 * th) * K_LY(ip2 + th, b_val))
    err_cd = abs(K_LY(th + ip2, b_val) * K_LY(th - ip2, b_val) * S_LY(2.0 * th) - 1.0)
    err_dg = abs(K_LY(th, b_val) * K_LY(th + 1.0j * np.pi, b_val) - S_LY(2.0 * th))
    report(
        f"Lee-Yang K_b: centred, Corrigan-Delius and Delius-Gandenberger crossing forms "
        f"(max err={max(err_book, err_cd, err_dg):.1e})",
        max(err_book, err_cd, err_dg) < 1e-13,
    )

    # Check 6: block algebra of Blackboard Calculation 3
    th = 0.71 + 0.17j
    xg = 1.37 + np.sqrt(2.0) / 10
    e_dbl = abs(blk6(2 * th, xg) + blk6(th, xg / 2) / blk6(th, 6 - xg / 2))
    e_s2 = abs(S_LY(2 * th) + blk6(th, 1) * blk6(th, 2) / (blk6(th, 4) * blk6(th, 5)))
    e_s2_printed = abs(S_LY(2 * th) + blk6(th, 2) * blk6(th, 4))
    e_shift = max(abs(K_LY(th + 1j * np.pi / 3, 0, eps) * K_LY(th - 1j * np.pi / 3, 0, eps)
                      + eps**2 * blk6(th, 3) * blk6(th, 5) / blk6(th, 2)) for eps in (1.0, -1.0))
    e_six = abs(blk6(th, 6) + 1)
    report(
        f"Lee-Yang blocks: (x)(2 theta) = -(x/2)/(6-x/2) ({e_dbl:.1e}), S(2 theta) = -(1)(2)/((4)(5)) "
        f"({e_s2:.1e}; -(2)(4) would be off by {e_s2_printed:.2f}), shifted product -eps^2(3)(5)/(2) "
        f"({e_shift:.1e}), (6) = -1",
        max(e_dbl, e_s2, e_shift, e_six) < 1e-13 and e_s2_printed > 1e-2,
    )

    # Check 7: unitarity and bootstrap sign
    th = 0.43 - 0.17j
    ip3 = 1.0j * np.pi / 3.0
    err_unit = abs(K_LY(th, b_val) * K_LY(-th, b_val) - 1.0)
    boot = lambda eps, bb: abs(K_LY(th, bb, eps)
                               - K_LY(th + ip3, bb, eps) * K_LY(th - ip3, bb, eps) * S_LY(2.0 * th))
    report(
        f"Lee-Yang unitarity ({err_unit:.1e}); bootstrap holds with eps=+1 "
        f"({max(boot(1.0, 0.0), boot(1.0, b_val)):.1e}) and fails with eps=-1 ({boot(-1.0, 0.0):.2f})",
        err_unit < 1e-13 and max(boot(1.0, 0.0), boot(1.0, b_val)) < 1e-13
        and abs(boot(-1.0, 0.0) - 2.0 * abs(K_LY(th, 0.0))) < 1e-12,
    )

    # Check 8: residues
    Rbulk, d0 = residue(S_LY, 2.0j * np.pi / 3.0)
    Gamma_sq = Rbulk / 1.0j
    R2, d1 = residue(lambda t: K_LY(t, 0.0), ip2)
    g_sq = R2 / 0.5j
    R6, d2 = residue(lambda t: K_LY(t, 0.0), 1.0j * np.pi / 6.0)
    errs = [abs(Gamma_sq + 2 * np.sqrt(3.0)),
            abs(g_sq + 4 * (2 * np.sqrt(3.0) - 3)),
            abs(R6 - 1j * (3 - np.sqrt(3.0))),
            abs(R6**2 + 0.25 * Gamma_sq * g_sq)]
    report(
        f"residues: Gamma^2 = {Gamma_sq.real:.10f} = -2 sqrt3, g^2 = {g_sq.real:.10f} = -4(2 sqrt3-3), "
        f"Res_(i pi/6) K_(1) = {R6.imag:.10f} i = i(3-sqrt3), square = -Gamma^2 g^2/4 "
        f"(max err={max(errs):.1e}; eps-change {max(d0, d1, d2):.1e})",
        max(errs) < 1e-8 and max(d0, d1, d2) < 1e-8,
    )

    # Check 9: CDD factor and f_x in blocks
    th = 0.71 + 0.17j
    bv = 0.37 + np.sqrt(2.0) / 100
    Phi = lambda t: blk6(t, 1 + bv) * blk6(t, 5 - bv) / (blk6(t, 1 - bv) * blk6(t, 5 + bv))
    e_cdd = max(abs(Phi(th) * Phi(-th) - 1), abs(Phi(ip2 - th) - Phi(ip2 + th)),
                abs(Phi(th) - Phi(th + ip3) * Phi(th - ip3)))
    xx = 1.23
    e_cross = abs(blk6(ip2 - th, xx) + blk6(ip2 + th, 6 - xx))
    e_fx = abs(f_x(th, xx / 6) + blk6(th, xx) * blk6(th, 6 - xx))
    e_ratio = abs(K_LY(th, bv) / K_LY(th, 0.0) - f_x(th, (1 + bv) / 6) / f_x(th, (1 - bv) / 6))
    report(
        f"Phi_b: unitary, CDD-crossing and bootstrap invariant ({e_cdd:.1e}); (x)(i pi/2-th) = -(6-x)(i pi/2+th) "
        f"({e_cross:.1e}); f_x = -(6x)(6-6x) ({e_fx:.1e}); K_b/K_(1) = f_((1+b)/6)/f_((1-b)/6) ({e_ratio:.1e})",
        max(e_cdd, e_cross, e_fx, e_ratio) < 1e-12,
    )

    # Check 10: excited boundary state and cancellations in the figure
    th = 0.43 - 0.17j
    u_pole = 1.0j * np.pi * (b_val + 1.0) / 6.0
    K_excited = S_LY(th + u_pole) * S_LY(th - u_pole) * K_LY(th, b_val)
    err_exc = max(abs(K_excited - K_LY(th, b_val + 2.0)), abs(K_LY(th, b_val + 2.0) - K_LY(th, 4.0 - b_val)))
    err_m2 = abs(K_LY(th, -2.0) - 1.0 / (blk6(th, 3) * blk6(th, 4) * blk6(th, 5)))
    err_0 = abs(K_LY(th, 0.0) - blk6(th, 1) * blk6(th, 3) / blk6(th, 4))
    report(
        f"K_b^(1) = K_(b+2) = K_(4-b) ({err_exc:.1e}); K_0 = K_(1) and K_(-2) = (3)^-1(4)^-1(5)^-1 "
        f"(no pole at i pi/6) ({max(err_0, err_m2):.1e})",
        err_exc < 1e-13 and max(err_0, err_m2) < 1e-13,
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
        r"$\theta=i\frac{\pi}{2}-\tilde\theta$",
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
        r"boundary state $|B\rangle$ at $\tilde t_E=0$",
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
        r"$S(2\tilde\theta)$",
        fontsize=8,
        color=GREEN,
        va="center",
    )
    ax1.text(
        v_x + 0.42,
        v_y + 0.25,
        r"$K(i\frac{\pi}{2}+\tilde\theta)$",
        fontsize=7.5,
        color=GREY,
    )

    ax1.text(
        0.0,
        -1.15,
        r"Boundary crossing-unitarity:"
        "\n"
        r"$K_a\!\left(i\frac{\pi}{2}-\tilde\theta\right)"
        r" = S_{a\bar a}(2\tilde\theta)\,K_{\bar a}\!\left(i\frac{\pi}{2}+\tilde\theta\right)$",
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
        r"$b=0$: $K_b=K_{(\mathbb{1})}$ (equal amplitudes);" "\n"
        r"zero of $(1-b)^{-1}$ cancels the moving pole",
        xy=(0.0, 1.0),
        xytext=(-2.7, 1.85),
        fontsize=7.0,
        color=PURPLE,
        arrowprops=dict(arrowstyle="->", color=PURPLE, lw=0.8),
    )

    ax2.plot(-2.0, 1.0, "x", color=GREY, ms=6, mew=1.5, zorder=6)
    ax2.annotate(
        r"$b=-2$: zero of $(1+b)$" "\n" r"cancels the pole at $i\pi/6$",
        xy=(-2.0, 0.95),
        xytext=(-2.9, 0.12),
        fontsize=6.5,
        color=GREY,
        arrowprops=dict(arrowstyle="->", color=GREY, lw=0.8),
    )

    ax2.plot(-1.0, 0.0, "s", color=GREEN, ms=4.5, zorder=6)
    ax2.annotate(
        r"threshold $b=-1$",
        xy=(-1.0, 0.05),
        xytext=(-0.75, 0.45),
        fontsize=7.5,
        color=GREEN,
        arrowprops=dict(arrowstyle="->", color=GREEN, lw=0.8),
    )

    ax2.set_xlim(-3.0, 2.0)
    ax2.set_ylim(-0.25, 3.35)
    ax2.set_xlabel(r"Boundary parameter $b$")
    ax2.set_ylabel(r"Pole position $\operatorname{Im}\theta$")
    ax2.set_yticks([0, 1, 2, 3])
    ax2.set_yticklabels([r"$0$", r"$\pi/6$", r"$\pi/3$", r"$\pi/2$"])
    ax2.set_title(r"(b) Candidate poles of the Lee–Yang $K_b(\theta)$", fontsize=9, pad=6)
    ax2.legend(loc="upper left", fontsize=6.8)

    fig.tight_layout()
    save(fig, "l1_f1_reflection_crossing", OUT_DIR)


if __name__ == "__main__":
    run_checks()
    make_figure()
    sys.exit(0 if all_ok() else 1)
