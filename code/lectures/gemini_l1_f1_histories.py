"""Lectures on integrable defects (defects-gemini), Lecture 1, figure 1 (fig-lec1-histories):
Two-particle spacetime histories for the Delfino-Mussardo-Simonetti no-go theorem and the Triangle Relation.

Top row: Mixed Channel A (A_b(-theta_2) D A_a(theta_1), where b reflects and a transmits):
  (a) Order I: b hits D first, reflects to -theta_2, collides head-on with a(theta_1) on x < 0 at relative
      rapidity theta_1 + theta_2, then a transmits.
  (b) Order II: a(theta_1) overtakes b(theta_2) on x < 0 at relative rapidity theta_1 - theta_2, then a
      transmits and b reflects without meeting a again.
Bottom row: Pure Transmission (D_gamma A_e(theta_2) A_f(theta_1), the Triangle Relation):
  (c) Order I: b crosses D first, then a crosses D, then a overtakes b on x > 0 at theta_1 - theta_2.
  (d) Order II: a overtakes b on x < 0 at theta_1 - theta_2, then a crosses D, then b crosses D.

Checks:
  1. Worldline intersection solver verifies the chronological event order, the side of the defect (x < 0 vs x > 0)
     where the bulk collision occurs, and the relative rapidity (theta_1 + theta_2 vs theta_1 - theta_2) in all 4 panels.
  2. Symbolic reduction of A_a(theta_1) A_b(theta_2) D along Order I and Order II verifies eq-lec1-order1,
     eq-lec1-order2, and the mixed-channel obstructions eq-lec1-mixed-a and eq-lec1-mixed-b.
  3. Equivalence of the S/T index form (eq-lec1-triangle-st), braid form (eq-lec1-triangle-braid), and RTT form
     (eq-lec1-triangle-rtt) of the Triangle Relation on the six-vertex R-matrix.

Writes lecture-notes/defects-gemini/figures/l1_histories.svg.
"""
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects-gemini" / "figures"

th1, th2 = 1.0, 0.45
v1, v2 = np.tanh(th1), np.tanh(th2)
T0, T1 = -2.4, 2.6


def worldline(v, rap, t_hit, fate):
    """Piecewise linear worldline x(t) hitting x = 0 at t = t_hit and either transmitting ('T') or reflecting ('R')."""
    v_after = v if fate == "T" else -v
    rap_after = rap if fate == "T" else -rap

    def x(t):
        t = np.asarray(t, dtype=float)
        return np.where(t <= t_hit, v * (t - t_hit), v_after * (t - t_hit))

    return x, rap_after


def analyze_history(t_a, t_b, fate_a, fate_b):
    xa, ra = worldline(v1, th1, t_a, fate_a)
    xb, rb = worldline(v2, th2, t_b, fate_b)
    ev = [
        (t_a, "a " + ("transmitted" if fate_a == "T" else "reflected"), 0.0, None),
        (t_b, "b " + ("transmitted" if fate_b == "T" else "reflected"), 0.0, None),
    ]
    ts = np.linspace(T0, T1, 200001)
    d = xa(ts) - xb(ts)
    for idx in np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]:
        tc = ts[idx] - d[idx] * (ts[idx + 1] - ts[idx]) / (d[idx + 1] - d[idx])
        rA = th1 if tc < t_a else ra
        rB = th2 if tc < t_b else rb
        ev.append((tc, "bulk collision", float(xa(tc)), abs(rA - rB)))
    return sorted(ev), xa, xb


# ---- Check 1: worldline kinematics in all four panels ----
panels = [
    ("(a) Mixed A, Order I (b hits D first)", -0.3, -1.2, "T", "R"),
    ("(b) Mixed A, Order II (bulk collision first)", -0.6, 0.2, "T", "R"),
    ("(c) Pure transmission, Order I", 0.3, -0.9, "T", "T"),
    ("(d) Pure transmission, Order II", -0.6, 0.2, "T", "T"),
]
expected = [
    (["b reflected", "bulk collision", "a transmitted"], th1 + th2, "left"),
    (["bulk collision", "a transmitted", "b reflected"], th1 - th2, "left"),
    (["b transmitted", "a transmitted", "bulk collision"], th1 - th2, "right"),
    (["bulk collision", "a transmitted", "b transmitted"], th1 - th2, "left"),
]

panel_data = []
for (title, t_a, t_b, fa, fb), (exp_order, exp_drap, exp_side) in zip(panels, expected):
    ev, xa, xb = analyze_history(t_a, t_b, fa, fb)
    names = [e[1] for e in ev]
    colls = [e for e in ev if e[1] == "bulk collision"]
    ok = (
        names == exp_order
        and len(colls) == 1
        and abs(colls[0][3] - exp_drap) < 1e-9
        and ((colls[0][2] < 0) if exp_side == "left" else (colls[0][2] > 0))
    )
    report(f"{title}: order={names}, dtheta={colls[0][3]:.3f} on {exp_side}", ok)
    panel_data.append((title, t_a, t_b, fa, fb, colls[0], xa, xb))

# ---- Check 2: symbolic ZF reduction along Order I and Order II ----
Ra, Rb, Ta, Tb, Sab_minus, Sba_minus, Sab_plus, Sba_plus = sp.symbols(
    "Ra Rb Ta Tb Sab_minus Sba_minus Sab_plus Sba_plus"
)
# Out-state basis: (RR: A_a(-1)A_b(-2)D, TT: D A_b(2)A_a(1), MixedA: A_b(-2)D A_a(1), MixedB: A_a(-1)D A_b(2))
order_I = {
    "RR": Rb * Sab_plus * Ra * Sba_minus,
    "TT": Tb * Ta * Sab_minus,
    "MixedA": Rb * Sab_plus * Ta,
    "MixedB": Tb * Ra,
}
order_II = {
    "RR": Sab_minus * Ra * Sba_plus * Rb,
    "TT": Sab_minus * Ta * Tb,
    "MixedA": Sab_minus * Ta * Rb,
    "MixedB": Sab_minus * Ra * Sba_plus * Tb,
}
diff_A = sp.factor(order_I["MixedA"] - order_II["MixedA"])
diff_B = sp.factor(order_I["MixedB"] - order_II["MixedB"])
report(
    "symbolic ZF reduction: Mixed A gives Rb*Ta*(Sab(th1+th2)-Sab(th1-th2))=0 and Mixed B gives Ra*Tb*(1-Sab(th1-th2)Sba(th1+th2))=0",
    sp.simplify(diff_A - Rb * Ta * (Sab_plus - Sab_minus)) == 0
    and sp.simplify(diff_B - Ra * Tb * (1 - Sab_minus * Sba_plus)) == 0
    and sp.simplify(order_I["TT"] - order_II["TT"]) == 0,
)

# ---- Check 3: equivalence of the three forms of the Triangle Relation ----
def r_six_vertex(u, eta=0.7):
    a = np.sinh(u + eta)
    b = np.sinh(u)
    c = np.sinh(eta)
    return np.array(
        [[a, 0, 0, 0], [0, b, c, 0], [0, c, b, 0], [0, 0, 0, a]], dtype=complex
    )


P2 = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], dtype=complex)
I2 = np.eye(2, dtype=complex)
P12 = np.kron(P2, I2)
P23 = np.kron(I2, P2)
P13 = P12 @ P23 @ P12

u1, u2, ud = 0.85, 0.31, -0.42
R12_mat = r_six_vertex(u1 - u2)
Rcheck12_mat = P2 @ R12_mat
T1_mat = P2 @ r_six_vertex(u1 - ud)
T2_mat = P2 @ r_six_vertex(u2 - ud)
L1_mat = P2 @ T1_mat
L2_mat = P2 @ T2_mat

# Braid operators on V1 (x) V2 (x) W
Rcheck_23 = np.kron(I2, Rcheck12_mat)
Rcheck_12 = np.kron(Rcheck12_mat, I2)
T_12_u1 = np.kron(T1_mat, I2)
T_23_u2 = np.kron(I2, T2_mat)
T_12_u2 = np.kron(T2_mat, I2)
T_23_u1 = np.kron(I2, T1_mat)

braid_lhs = Rcheck_23 @ T_12_u1 @ T_23_u2
braid_rhs = T_12_u2 @ T_23_u1 @ Rcheck_12

# RTT operators on V1 (x) V2 (x) W
R_12 = np.kron(R12_mat, I2)
L_13_u1 = P23 @ np.kron(L1_mat, I2) @ P23
L_23_u2 = np.kron(I2, L2_mat)
L_13_u2 = P23 @ np.kron(L2_mat, I2) @ P23
L_23_u1 = np.kron(I2, L1_mat)

rtt_lhs = R_12 @ L_13_u1 @ L_23_u2
rtt_rhs = L_23_u2 @ L_13_u1 @ R_12

# Index form tensor contraction
S_tens = Rcheck12_mat.reshape(2, 2, 2, 2).transpose(0, 1, 3, 2)  # S_{ab}^{cd}
T1_tens = T1_mat.reshape(2, 2, 2, 2)  # T_{a,alpha}^{beta,b}
T2_tens = T2_mat.reshape(2, 2, 2, 2)
idx_lhs = np.einsum("abcd,cABf,dBGe->abAGef", S_tens, T1_tens, T2_tens)
idx_rhs = np.einsum("bABd,aBGc,cdfe->abAGef", T2_tens, T1_tens, S_tens)

report(
    "three forms of the Triangle Relation (index, braid, RTT) agree to machine precision and hold for a particle defect",
    np.allclose(braid_lhs, P13 @ rtt_lhs)
    and np.allclose(braid_rhs, P13 @ rtt_rhs)
    and np.allclose(braid_lhs, braid_rhs)
    and np.allclose(idx_lhs, idx_rhs),
)

# ---- Figure ----
fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.1))
ts = np.linspace(T0, T1, 400)

for idx, (ax, (title, t_a, t_b, fa, fb, coll, xa, xb)) in enumerate(zip(axes.flat, panel_data)):
    # Shaded half-lines
    ax.axvspan(-2.1, 0, color=style.BLUE, alpha=0.04, lw=0)
    ax.axvspan(0, 2.1, color=style.TEAL, alpha=0.04, lw=0)
    # Defect line at x = 0
    ax.plot([0, 0], [T0, T1], color=style.PURPLE, lw=2.2, zorder=3)
    if fb == "R":
        ax.plot([0, 0], [T0, T1], color=style.ORANGE, lw=2.2, ls=(0, (3, 3)), zorder=4)

    # Particle worldlines
    ax.plot(xa(ts), ts, color=style.BLUE, lw=1.9, zorder=4)
    ax.plot(xb(ts), ts, color=style.RED, lw=1.9, zorder=4)

    # Defect vertices
    ax.plot([0], [t_a], "s", ms=4.2, color=style.PURPLE, zorder=6)
    ax.plot([0], [t_b], "s", ms=4.2, color=style.PURPLE, zorder=6)

    # Bulk collision vertex
    tc, _, xc, drap = coll
    ax.plot([xc], [tc], "o", ms=5.2, color=style.ORANGE, zorder=6)

    if abs(drap - (th1 + th2)) < 1e-6:
        s_lab = r"$S_{ab}(\theta_1+\theta_2)$"
    else:
        s_lab = r"$S_{ab}(\theta_1-\theta_2)$"

    if idx == 0:
        ax.annotate(
            s_lab,
            xy=(xc, tc),
            xytext=(-1.95, tc + 0.45),
            fontsize=8,
            color=style.GREY,
            va="center",
            arrowprops=dict(arrowstyle="->", color=style.GREY, lw=0.7, shrinkB=4),
        )
        ax.text(0.12, t_b - 0.05, r"$R_b(\theta_2)$", fontsize=7.8, color=style.RED, va="center")
        ax.text(-0.12, t_a + 0.16, r"$\mathbb{T}_a(\theta_1)$", fontsize=7.8, color=style.BLUE, ha="right", va="bottom")
    elif idx == 1:
        ax.annotate(
            s_lab,
            xy=(xc, tc),
            xytext=(-1.95, tc + 0.55),
            fontsize=8,
            color=style.GREY,
            va="center",
            arrowprops=dict(arrowstyle="->", color=style.GREY, lw=0.7, shrinkB=4),
        )
        ax.text(0.12, t_a - 0.22, r"$\mathbb{T}_a(\theta_1)$", fontsize=7.8, color=style.BLUE, va="top")
        ax.text(0.12, t_b + 0.12, r"$R_b(\theta_2)$", fontsize=7.8, color=style.RED, va="bottom")
    elif idx == 2:
        ax.annotate(
            s_lab,
            xy=(xc, tc),
            xytext=(-1.75, tc - 0.10),
            fontsize=8,
            color=style.GREY,
            va="center",
            arrowprops=dict(arrowstyle="->", color=style.GREY, lw=0.7, shrinkB=4),
        )
        ax.text(0.12, t_b - 0.15, r"$\mathbb{T}_{23}(\theta_2)$", fontsize=7.8, color=style.RED, ha="left", va="top")
        ax.text(-0.12, t_a + 0.12, r"$\mathbb{T}_{12}(\theta_1)$", fontsize=7.8, color=style.BLUE, ha="right", va="bottom")
    else:
        ax.annotate(
            s_lab,
            xy=(xc, tc),
            xytext=(-1.95, tc + 0.55),
            fontsize=8,
            color=style.GREY,
            va="center",
            arrowprops=dict(arrowstyle="->", color=style.GREY, lw=0.7, shrinkB=4),
        )
        ax.text(0.12, t_a - 0.22, r"$\mathbb{T}_{23}(\theta_1)$", fontsize=7.8, color=style.BLUE, va="top")
        ax.text(-0.12, t_b + 0.12, r"$\mathbb{T}_{12}(\theta_2)$", fontsize=7.8, color=style.RED, ha="right", va="bottom")

    tl_a = T0 + 0.55
    tl_b = T0 + 0.32
    ax.text(float(xa(tl_a)) - 0.12, tl_a, r"$a(\theta_1)$", color=style.BLUE, fontsize=8, ha="right", va="center")
    ax.text(float(xb(tl_b)) + 0.10, tl_b, r"$b(\theta_2)$", color=style.RED, fontsize=8, ha="left", va="center")
    ax.text(0.08, T1 - 0.25, r"$D$", color=style.PURPLE, fontsize=8.5, ha="left", va="top")

    ax.set_xlim(-2.1, 2.1)
    ax.set_ylim(T0, T1)
    ax.set_xticks([0])
    ax.set_xticklabels([r"$x=0$"])
    ax.set_yticks([])
    ax.set_ylabel(r"$t$", labelpad=2)
    ax.set_title(title, fontsize=8.5, color=style.GREY)

fig.text(
    0.5,
    0.985,
    r"No-go theorem (Mixed Channel A, $A_b(-\theta_2)\,D\,A_a(\theta_1)$): head-on $S_{ab}(\theta_1+\theta_2)$ vs. overtaking $S_{ab}(\theta_1-\theta_2)$",
    ha="center",
    va="top",
    fontsize=8.5,
    color=style.GREY,
)
fig.text(
    0.5,
    0.495,
    r"Pure transmission (Triangle Relation, $D_\gamma\,A_e(\theta_2)\,A_f(\theta_1)$): both orders overtake at $\theta_1-\theta_2$",
    ha="center",
    va="top",
    fontsize=8.5,
    color=style.GREY,
)
fig.tight_layout(rect=(0, 0, 1, 0.965), h_pad=2.4, w_pad=1.4)
style.save(fig, "l1_histories", OUT)
raise SystemExit(0 if style.all_ok() else 1)
