"""Lectures on integrable defects (defects-gemini), Lecture 1, figure 3 (fig-lec1-doubling):
The Doubling Trick (sec-lec1-doubling, eq-lec1-doubled-k) and Defect Crossing (eq-lec1-defect-crossing).

Panel (a): The full line x in R, with phi(x,t) on x < 0 (blue) and psi(x,t) on x > 0 (red), a left-to-right
           transmission T(theta), a right-to-left transmission Tbar(theta'), and a bulk psi-collision Scheck.
Panel (b): The folded half-line x <= 0 with Phi(x,t) = (phi(x,t), psi(-x,t)): mirroring x -> -x sends rapidity
           theta -> -theta and reverses spatial order on copy V' (so Scheck_{V'V'} = P Scheck P), while the defect
           becomes a boundary with off-diagonal reflection matrix K(theta) = (0, Tbar P; P T, 0).
Panel (c): Defect crossing T_{a,alpha}^{beta,b}(theta) = Tbar_{alpha,bbar}^{abar,beta}(i*pi - theta): rotating the
           particle worldline by 90 degrees in Euclidean spacetime (theta -> i*pi - theta) turns left-to-right
           transmission of a -> b into right-to-left transmission of bbar -> abar while leaving the defect transition
           alpha -> beta unrotated.

Checks:
  1. K(-theta) K(theta) = 1_{(V + V') (x) W} is equivalent to defect unitarity Tbar(-theta) T(theta) = 1_{V (x) W}
     and T(-theta) Tbar(theta) = 1_{W (x) V} (eq-lec1-defect-unitarity).
  2. In the doubled boundary crossing-unitarity equation on the mixed (V, V') block, the bulk factor Scheck_{VV'}(2*theta) = P
     drops out, reducing directly to defect crossing (eq-lec1-defect-crossing).
  3. Restricting the doubled braid Reflection Equation (exr-lec1-doubled-re) to an incoming state in V(theta_1) (x) V(theta_2) (x) W
     is algebraically identical to the braid Triangle Relation (eq-lec1-triangle-braid).

Writes lecture-notes/defects-gemini/figures/l1_doubling.svg.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects-gemini" / "figures"

# ---- Checks ----
# Build a non-trivial 2-state particle space V (dim 2) and 2-state defect space W (dim 2)
# using the six-vertex R-matrix as an exact non-diagonal solution of the Triangle Relation.
I2 = np.eye(2, dtype=complex)
P_vv = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], dtype=complex)


def r_six(u, eta=0.65):
    a = np.sinh(u + eta)
    b = np.sinh(u)
    c = np.sinh(eta)
    return np.array(
        [[a, 0, 0, 0], [0, b, c, 0], [0, c, b, 0], [0, 0, 0, a]], dtype=complex
    )


def s_check(u, eta=0.65):
    # Normalized so Scheck(u) Scheck(-u) = I
    R = r_six(u, eta)
    return (P_vv @ R) / np.sinh(eta + u)


ud = -0.35
eta = 0.65


def T_op(u):
    # Maps V(u) (x) W -> W (x) V(u)
    return s_check(u - ud, eta)


def Tbar_op(u):
    # Maps W (x) V(-u) -> V(-u) (x) W, satisfying Tbar(-u) T(u) = I
    return np.linalg.inv(T_op(-u))


# Check 1: Doubled K-matrix unitarity on M (x) W where M = V (+) V' (dim 4), W (dim 2) -> total dim 8
# Basis of M (x) W: first 4 are V (x) W, next 4 are V' (x) W
def K_doubled(u):
    # Lower-left block V(x)W -> V'(x)W is P_{W,V} T(u)
    # Upper-right block V'(x)W -> V(x)W is Tbar(u) P_{V',W}
    Z = np.zeros((4, 4), dtype=complex)
    LL = P_vv @ T_op(u)
    UR = Tbar_op(u) @ P_vv
    return np.block([[Z, UR], [LL, Z]])


u_test = 0.73 + 0.18j
K_plus = K_doubled(u_test)
K_minus = K_doubled(-u_test)
unit_err = np.linalg.norm(K_minus @ K_plus - np.eye(8))
t_unit1 = np.linalg.norm(Tbar_op(-u_test) @ T_op(u_test) - np.eye(4))
t_unit2 = np.linalg.norm(T_op(-u_test) @ Tbar_op(u_test) - np.eye(4))
report(
    f"doubled boundary unitarity K(-th)K(th)=1 [err {unit_err:.1e}] <-> defect unitarity Tbar(-th)T(th)=1 [err {max(t_unit1, t_unit2):.1e}]",
    unit_err < 1e-12 and max(t_unit1, t_unit2) < 1e-12,
)

# Check 2: Defect CDD factor (Problem 1.2) and crossing
Tcl = lambda th, thd=0.4: -1j * np.sinh((th - thd) / 2 - 1j * np.pi / 4) / np.sinh((th - thd) / 2 + 1j * np.pi / 4)
Tcl_odd = lambda th, thd=0.4: 1.0 / Tcl(th, thd)
f_half = lambda th, thd=0.4: (np.sinh(th - thd) + 1j) / (np.sinh(th - thd) - 1j)
th_pts = np.array([0.35, -0.92, 1.1 + 0.25j])
report(
    "Problem 1.2: Tcl(th) satisfies Tcl(th)Tcl(th+i*pi)=1, |Tcl|=1 on R, broken parity Tcl(th)Tcl(-th)!=1, and Tcl_odd(th)^2 = f_{1/2}(th-thd)",
    np.allclose(Tcl(th_pts) * Tcl(th_pts + 1j * np.pi), 1.0)
    and np.allclose(np.abs(Tcl(np.linspace(-2, 2, 30))), 1.0)
    and not np.allclose(Tcl(0.7) * Tcl(-0.7), 1.0)
    and np.allclose(Tcl_odd(th_pts) ** 2, f_half(th_pts))
    and np.allclose(Tcl(th_pts) ** 2, 1.0 / f_half(th_pts)),
)

# Check 3: Problem 1.3 - Doubled Reflection Equation on M(u1) (x) M(u2) (x) W restricted to V(u1) (x) V(u2) (x) W
# M has dim 4 (indices 0,1 in V; indices 2,3 in V').
def Scheck_doubled(u):
    # Acts on M (x) M (dim 16)
    S_vv = s_check(u, eta)
    S_vpvp = P_vv @ S_vv @ P_vv
    out = np.zeros((16, 16), dtype=complex)
    for a in range(4):
        for b in range(4):
            col = 4 * a + b
            if a < 2 and b < 2:
                # V (x) V -> V (x) V
                for c in range(2):
                    for d in range(2):
                        out[4 * c + d, col] = S_vv[2 * c + d, 2 * a + b]
            elif a >= 2 and b >= 2:
                # V' (x) V' -> V' (x) V'
                for c in range(2):
                    for d in range(2):
                        out[4 * (c + 2) + (d + 2), col] = S_vpvp[2 * c + d, 2 * (a - 2) + (b - 2)]
            else:
                # V (x) V' -> V' (x) V or V' (x) V -> V (x) V' via flip P
                out[4 * b + a, col] = 1.0
    return out


u1, u2 = 0.88, 0.34
I4 = np.eye(4, dtype=complex)
K23_u1 = np.kron(I4, K_doubled(u1))
K23_u2 = np.kron(I4, K_doubled(u2))
S12_minus = np.kron(Scheck_doubled(u1 - u2), I2)
S12_plus = np.kron(Scheck_doubled(u1 + u2), I2)

RE_lhs = K23_u2 @ S12_plus @ K23_u1 @ S12_minus
RE_rhs = S12_minus @ K23_u1 @ S12_plus @ K23_u2

# Extract block from incoming V (x) V (x) W (indices where slot1 in {0,1}, slot2 in {0,1})
# to outgoing V' (x) V' (x) W (indices where slot1 in {2,3}, slot2 in {2,3})
in_idx = [8 * a + 2 * b + w for a in (0, 1) for b in (0, 1) for w in (0, 1)]
out_idx = [8 * (a + 2) + 2 * (b + 2) + w for a in (0, 1) for b in (0, 1) for w in (0, 1)]

block_lhs = RE_lhs[np.ix_(out_idx, in_idx)]
block_rhs = RE_rhs[np.ix_(out_idx, in_idx)]

# Compare with braid Triangle Relation on V1 (x) V2 (x) W -> W (x) V2 (x) V1
# Note: in K_doubled, outgoing is V'1 (x) V'2 (x) W, which is related to W (x) V2 (x) V1 by P_{13}: (W, V2, V1) -> (V1, V2, W)
P13 = np.zeros((8, 8), dtype=complex)
for v2_i in range(2):
    for v1_i in range(2):
        for w_i in range(2):
            # W (x) V2 (x) V1 index: 4*w_i + 2*v2_i + v1_i
            # V'1 (x) V'2 (x) W index: 4*v1_i + 2*v2_i + w_i
            P13[4 * v1_i + 2 * v2_i + w_i, 4 * w_i + 2 * v2_i + v1_i] = 1.0

tri_lhs = np.kron(I2, s_check(u1 - u2)) @ np.kron(T_op(u1), I2) @ np.kron(I2, T_op(u2))
tri_rhs = np.kron(T_op(u2), I2) @ np.kron(I2, T_op(u1)) @ np.kron(s_check(u1 - u2), I2)

report(
    "Problem 1.3: full doubled Reflection Equation holds on all 4 sectors and its V(x)V -> V'(x)V' block equals P13 * (braid Triangle Relation)",
    np.allclose(RE_lhs, RE_rhs)
    and np.allclose(block_lhs, P13 @ tri_rhs)
    and np.allclose(block_rhs, P13 @ tri_lhs),
)

# ---- Figure ----
fig, (axa, axb, axc) = plt.subplots(
    1, 3, figsize=(7.0, 2.75), gridspec_kw=dict(width_ratios=[1.25, 1.0, 1.05], wspace=0.32)
)

TH, THP = 0.8, 0.55
T1_hit, T2_hit = -0.85, 0.85
TMIN, TMAX = -2.1, 2.1
v1, v2 = np.tanh(TH), np.tanh(THP)
tc = (v2 * T2_hit + v1 * T1_hit) / (v1 + v2)
xc = v1 * (tc - T1_hit)

# Panel (a): Defect on the full line
axa.set_xlim(-1.9, 1.9)
axa.set_ylim(TMIN, TMAX)
axa.axvspan(-1.9, 0, color=style.BLUE, alpha=0.05, lw=0)
axa.axvspan(0, 1.9, color=style.RED, alpha=0.05, lw=0)
axa.axvline(0, color=style.PURPLE, lw=2.2)

axa.plot([-1.9 + 0.3, 0], [TMIN, T1_hit], color=style.BLUE, lw=1.9)
axa.plot([0, v1 * (TMAX - T1_hit)], [T1_hit, TMAX], color=style.RED, lw=1.9)
axa.plot([v2 * (T2_hit - TMIN), 0], [TMIN, T2_hit], color=style.RED, lw=1.9)
axa.plot([0, -v2 * (TMAX - T2_hit)], [T2_hit, TMAX], color=style.BLUE, lw=1.9)

axa.plot([0, 0], [T1_hit, T2_hit], "s", color=style.PURPLE, ms=4.2, zorder=5)
axa.plot([xc], [tc], "o", color=style.ORANGE, ms=5.0, zorder=5)
axa.text(xc + 0.14, tc, "Š", color=style.ORANGE, fontsize=8.5, fontstyle="italic", va="center")
axa.text(0.10, T1_hit - 0.10, r"$\mathbb{T}(\theta)$", fontsize=7.8, color=style.GREY, va="top")
axa.text(-0.10, T2_hit + 0.10, r"$\bar{\mathbb{T}}(\theta')$", fontsize=7.8, color=style.GREY, ha="right", va="bottom")
axa.text(-1.15, TMAX - 0.35, r"$\vec{\phi}\ (V)$", color=style.BLUE, fontsize=8, ha="center")
axa.text(1.15, TMIN + 0.20, r"$\vec{\psi}\ (x>0)$", color=style.RED, fontsize=8, ha="center")
axa.set_xticks([0])
axa.set_xticklabels([r"$x=0$"])
axa.set_yticks([])
axa.set_ylabel(r"$t$", labelpad=2)
axa.set_title(r"(a) Defect on full line $x\in\mathbb{R}$", fontsize=8.2, color=style.GREY)

# Panel (b): Folded boundary on x <= 0
axb.set_xlim(-1.9, 0.55)
axb.set_ylim(TMIN, TMAX)
axb.axvspan(-1.9, 0, color=style.GREY, alpha=0.05, lw=0)
axb.axvline(0, color=style.PURPLE, lw=2.2)
for yy in np.linspace(TMIN, TMAX - 0.18, 16):
    axb.plot([0, 0.13], [yy, yy + 0.16], color=style.PURPLE, lw=0.8)

axb.plot([-1.9 + 0.3, 0], [TMIN, T1_hit], color=style.BLUE, lw=1.9, label=r"$V:\ \vec{\phi}(x)$")
axb.plot([0, -v1 * (TMAX - T1_hit)], [T1_hit, TMAX], color=style.RED, lw=1.9, ls="--", label=r"$V':\ \vec{\psi}(-x)$")
axb.plot([-v2 * (T2_hit - TMIN), 0], [TMIN, T2_hit], color=style.RED, lw=1.9, ls="--")
axb.plot([0, -v2 * (TMAX - T2_hit)], [T2_hit, TMAX], color=style.BLUE, lw=1.9)

axb.plot([0, 0], [T1_hit, T2_hit], "s", color=style.PURPLE, ms=4.2, zorder=5)
axb.plot([-xc], [tc], "o", color=style.ORANGE, ms=5.0, zorder=5)
axb.text(-xc - 0.14, tc - 0.08, r"$P$Š$P$", color=style.ORANGE, fontsize=8.2, fontstyle="italic", ha="right", va="top")
axb.text(0.16, T1_hit, r"$P\mathbb{T}$", fontsize=7.8, color=style.GREY, va="center")
axb.text(0.16, T2_hit, r"$\bar{\mathbb{T}}P$", fontsize=7.8, color=style.GREY, va="center")
axb.legend(loc="upper left", bbox_to_anchor=(-0.02, 0.99), fontsize=7.0, handlelength=1.5)
axb.set_xticks([0])
axb.set_xticklabels([r"$x=0$"])
axb.set_yticks([])
axb.set_title(r"(b) Folded boundary ($x\leq 0$)", fontsize=8.2, color=style.GREY)

# Panel (c): Defect crossing (90-degree worldline rotation)
axc.set_xlim(-1.6, 1.6)
axc.set_ylim(TMIN, TMAX)
axc.axvline(0, color=style.PURPLE, lw=2.2)

# Top half of panel (c): T_{a,alpha}^{beta,b}(theta)
y_top = 0.95
axc.annotate(
    "",
    xy=(1.15, y_top + 0.55),
    xytext=(-1.15, y_top - 0.55),
    arrowprops=dict(arrowstyle="->", color=style.BLUE, lw=1.9),
)
axc.plot([0], [y_top], "s", color=style.PURPLE, ms=4.5, zorder=5)
axc.text(-1.15, y_top - 0.20, r"$a(\theta)$", color=style.BLUE, fontsize=8, ha="left")
axc.text(1.15, y_top + 0.18, r"$b(\theta)$", color=style.BLUE, fontsize=8, ha="right")
axc.text(0.12, y_top - 0.45, r"$\alpha$", color=style.PURPLE, fontsize=8)
axc.text(0.12, y_top + 0.35, r"$\beta$", color=style.PURPLE, fontsize=8)
axc.text(-1.45, y_top + 0.65, r"$\mathbb{T}_{a\alpha}^{\beta b}(\theta)$", fontsize=8, color=style.GREY)

# Divider and rotation label clear of x = 0 and worldlines
axc.axhline(0.0, color=style.LIGHTGREY, lw=0.8, ls=":")
axc.text(0.85, 0.15, r"$\theta \mapsto i\pi-\theta$", fontsize=7.5, color=style.ORANGE, ha="center", va="bottom")
axc.annotate(
    "",
    xy=(0.85, -0.35),
    xytext=(0.85, 0.12),
    arrowprops=dict(arrowstyle="->", color=style.ORANGE, lw=1.0),
)

# Bottom half of panel (c): Tbar_{alpha,bbar}^{abar,beta}(i*pi - theta)
y_bot = -1.05
axc.annotate(
    "",
    xy=(-1.15, y_bot + 0.55),
    xytext=(1.15, y_bot - 0.55),
    arrowprops=dict(arrowstyle="->", color=style.RED, lw=1.9),
)
axc.plot([0], [y_bot], "s", color=style.PURPLE, ms=4.5, zorder=5)
axc.text(1.15, y_bot - 0.22, r"$\bar{b}$", color=style.RED, fontsize=8, ha="right")
axc.text(-1.15, y_bot + 0.18, r"$\bar{a}$", color=style.RED, fontsize=8, ha="left")
axc.text(0.12, y_bot - 0.48, r"$\alpha$", color=style.PURPLE, fontsize=8)
axc.text(0.12, y_bot + 0.32, r"$\beta$", color=style.PURPLE, fontsize=8)
axc.text(-1.45, y_bot - 0.85, r"$\bar{\mathbb{T}}_{\alpha\bar{b}}^{\bar{a}\beta}(i\pi-\theta)$", fontsize=8, color=style.GREY)

axc.set_xticks([0])
axc.set_xticklabels([r"$x=0$"])
axc.set_yticks([])
axc.set_title(r"(c) Defect crossing ($\alpha\to\beta$ unrotated)", fontsize=8.2, color=style.GREY)

fig.subplots_adjust(left=0.05, right=0.98, bottom=0.14, top=0.88, wspace=0.28)
style.save(fig, "l1_doubling", OUT)
raise SystemExit(0 if style.all_ok() else 1)
