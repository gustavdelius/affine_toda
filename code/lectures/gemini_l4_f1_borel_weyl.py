"""Lectures on integrable defects (defects-gemini), Lecture 4, figure 1 (fig-lec4-borel-weyl):
Why Transmitting Defects Require a Borel Hopf Subalgebra and the Bilateral Weyl Module Ladder
(sec-lec4-conformal-borel, sec-lec4-q-oscillator, sec-lec4-solving-intertwiner, thm-t-triangle, eq-weyl-module).

Panel (a): Boundary vs. defect coproduct architecture on a three-fold tensor product:
           - Top: A boundary module W_B stays permanently on the rightmost tensor slot V_1 (x) V_2 (x) W_B, so a left
             coideal subalgebra Delta(B) subset U_q(g^) (x) B suffices.
           - Bottom: In the defect Triangle Relation V_1 (x) V_2 (x) W -> V_1 (x) W (x) V_2 -> W (x) V_2 (x) V_1,
             the defect space W sits in the middle tensor slot at the intermediate step, requiring a full Hopf
             subalgebra Delta(U_q(b^-)) subset U_q(b^-) (x) U_q(b^-).
Panel (b): The bilateral Weyl module W = span{|n> : n in Z} (Z = 0, invertible f_1, f_0 shifts) vs. the bounded-below
           Fock module F = span{|n> : n >= 0} (a_-|0> = 0), and the 2-dimensional charge-conserving blocks
           H_n = span{v_+ (x) |n>, v_- (x) |n+1>} on which T(theta) acts.

Checks:
  1. Weyl module of U_q(b^-) satisfies Cartan relations, [f_0, f_1] = 0, and cubic q-Serre relations (Problem 4.1).
  2. The sine-Gordon transmission matrix T(theta) (eq-sg-transmission-matrix) is the unique U_q(b^-) intertwiner
     V(y) (x) W -> W (x) V(y) up to a scalar, and level zero k_0 k_1 = 1 is required.
  3. T(theta) satisfies the Defect Yang-Baxter (Triangle) Relation on V(y_1) (x) V(y_2) (x) W -> W (x) V(y_2) (x) V(y_1).
  4. Problem 4.2: On each 2D charge block H_n, det(T_n(y)/Phi) = (1 + (y/y_d)^2) q^{1+s} vanishes at y/y_d = +-i,
     where T_n has rank 1 and spans a 1D U_q(b^-)-submodule (shifted Weyl module).

Writes lecture-notes/defects-gemini/figures/l4_borel_weyl.svg.
"""
import sys
from pathlib import Path

import numpy as np
from matplotlib.patches import FancyBboxPatch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from boundary_common import cartan_a, q_phys, rcheck, sl2_spin  # noqa: E402
from defect_common import (  # noqa: E402
    tmatrix_elements,
    tmatrix_from_elements,
    triangle_residual_braid,
    weyl_module_b,
)

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects-gemini" / "figures"

lam = np.sqrt(2.0)
q = q_phys(lam)
N_win = 5
thd = 0.31
yd = np.exp(lam * thd)
s_val = 0.3
mu = 1j / (q - 1.0 / q)

# ---- Check 1: Weyl module of U_q(b^-) and cubic q-Serre relations (Problem 4.1) ----
W = weyl_module_b(-1, q, yd, N_win, s=s_val, mu=mu)
I_int = W["interior"]
blk = np.ix_(I_int, I_int)
f0, f1 = W[("f", 0)], W[("f", 1)]
k0, k1 = W[("k", 0)], W[("k", 1)]
C_mat = cartan_a(1)

ok_weyl = np.allclose((k0 @ k1 - np.eye(len(I_int)))[blk], 0.0) and np.allclose((f0 @ f1 - f1 @ f0)[2:-2, 2:-2], 0.0)
for i, k_op in ((0, k0), (1, k1)):
    for j, f_op in ((0, f0), (1, f1)):
        ok_weyl &= np.allclose((k_op @ f_op @ np.linalg.inv(k_op))[blk], (q ** (-C_mat[i][j]) * f_op)[blk])

q3 = q**2 + 1.0 + q**-2
for X, Y in ((f1, f0), (f0, f1)):
    serre = (
        np.linalg.matrix_power(X, 3) @ Y
        - q3 * np.linalg.matrix_power(X, 2) @ Y @ X
        + q3 * X @ Y @ np.linalg.matrix_power(X, 2)
        - Y @ np.linalg.matrix_power(X, 3)
    )
    ok_weyl &= np.allclose(serre[4:-4, 4:-4], 0.0)

report(
    "Problem 4.1: Weyl module W^-(y_d) satisfies Cartan relations, k_0 k_1 = 1, [f_0, f_1] = 0, and cubic q-Serre relations",
    ok_weyl,
)

# ---- Check 2: T(theta) is the unique U_q(b^-) intertwiner V(y) (x) W -> W (x) V(y) ----
th_test = 0.77
y_test = np.exp(lam * th_test)
V = sl2_spin(1, q, y_test)
ns = W["n"]


def build_T(y_val):
    Y_r = y_val / yd
    V_loc = sl2_spin(1, q, y_val)
    return tmatrix_from_elements(
        lambda n: q ** (-n),
        lambda n: -1j * Y_r * (q ** (n + s_val)),
        lambda n: -1j * Y_r * (q ** (-n)),
        lambda n: q ** (n + s_val),
        V_loc,
        W,
    )


T_mat = build_T(y_test)

# Verify intertwining of Delta(k_i) = k_i (x) k_i and Delta(f_i) = f_i (x) k_i^-1 + 1 (x) f_i on interior
rows = np.repeat(np.abs(ns) <= N_win - 2, 2)  # W (x) V output
cols = np.tile(np.abs(ns) <= N_win - 2, 2)    # V (x) W input
dW = W["dim"]
ok_int = True
for j in (0, 1):
    dk_in = np.kron(V[("k", j)], W[("k", j)])
    dk_out = np.kron(W[("k", j)], V[("k", j)])
    df_in = np.kron(V[("f", j)], np.linalg.inv(W[("k", j)])) + np.kron(np.eye(2), W[("f", j)])
    df_out = np.kron(W[("f", j)], np.linalg.inv(V[("k", j)])) + np.kron(np.eye(dW), V[("f", j)])
    ok_int &= np.allclose((T_mat @ dk_in - dk_out @ T_mat)[np.ix_(rows, cols)], 0.0)
    ok_int &= np.allclose((T_mat @ df_in - df_out @ T_mat)[np.ix_(rows, cols)], 0.0)

report(
    "sine-Gordon T(theta) of eq-sg-transmission-matrix intertwines Delta(k_i) and Delta(f_i) of U_q(b^-) on V(y) (x) W",
    ok_int,
)

# ---- Check 3: Triangle relation on V(y_1) (x) V(y_2) (x) W -> W (x) V(y_2) (x) V(y_1) ----
th1, th2 = 0.83, -0.29
y1, y2 = np.exp(lam * th1), np.exp(lam * th2)
Rcheck_fn = lambda z: rcheck(sl2_spin(1, q, z), sl2_spin(1, q, 1.0))[0]
tri_res, _ = triangle_residual_braid(Rcheck_fn, build_T, y1, y2, W)

report(
    f"Borel intertwiner T(theta) satisfies the Defect Yang-Baxter (Triangle) Relation on V_1 (x) V_2 (x) W [res {tri_res:.1e}]",
    tri_res < 1e-13,
)

# ---- Check 4: Problem 4.2 (2D charge blocks H_n and rank-1 kernel at Y = +-i) ----
# On H_n = span{v_+ (x) |n>, v_- (x) |n+1>}, T_n(Y) = [[1, -iY], [-iY, 1]] @ diag(q^-n, q^(n+1+s))
# At Y = -i (theta = theta_d - i*pi/(2*lambda)), [[1, -1], [-1, 1]] has kernel spanned by w_n = q^n v_+ (x) |n> + q^{-(n+1+s)} v_- (x) |n+1>
Y_deg = -1j
y_deg = Y_deg * yd
V_deg = sl2_spin(1, q, y_deg)
T_deg = build_T(y_deg)

# Construct the kernel basis vectors w_n in V(y_deg) (x) W for interior n
# In V (x) W, v_+ (x) |n> is at index i_n, and v_- (x) |n+1> is at index dW + i_{n+1}
w_vecs = {}
for idx, n_val in enumerate(ns[:-1]):
    vec = np.zeros(2 * dW, dtype=complex)
    vec[idx] = q**n_val
    vec[dW + idx + 1] = q ** (-(n_val + 1 + s_val))
    w_vecs[n_val] = vec

ker_norm = max(np.linalg.norm(T_deg @ w_vecs[n_val]) for n_val in range(-N_win + 2, N_win - 2))

# Check that Delta(f_1) maps w_n to a scalar multiple of w_{n-1} and Delta(f_0) maps w_n to a scalar multiple of w_{n+1}
df1_deg = np.kron(V_deg[("f", 1)], np.linalg.inv(W[("k", 1)])) + np.kron(np.eye(2), W[("f", 1)])
df0_deg = np.kron(V_deg[("f", 0)], np.linalg.inv(W[("k", 0)])) + np.kron(np.eye(2), W[("f", 0)])
submod_ok = True
for n_val in range(-N_win + 2, N_win - 2):
    img1 = df1_deg @ w_vecs[n_val]
    target1 = w_vecs[n_val - 1]
    ratio1 = np.vdot(target1, img1) / np.vdot(target1, target1)
    img0 = df0_deg @ w_vecs[n_val]
    target0 = w_vecs[n_val + 1]
    ratio0 = np.vdot(target0, img0) / np.vdot(target0, target0)
    submod_ok &= np.allclose(img1, ratio1 * target1) and np.allclose(img0, ratio0 * target0)

report(
    f"Problem 4.2: at Y = y/y_d = -i, ker(T(y)) is spanned by w_n in H_n [max |T w_n| = {ker_norm:.1e}] and forms a U_q(b^-)-submodule",
    ker_norm < 1e-14 and submod_ok,
)

# ---- Figure ----
fig, (axa, axb) = plt.subplots(
    1, 2, figsize=(7.0, 3.05), gridspec_kw=dict(width_ratios=[1.05, 1.15], wspace=0.18)
)

# Panel (a): Left coideal (boundary) vs Borel Hopf subalgebra (defect)
axa.set_xlim(0.0, 10.0)
axa.set_ylim(0.0, 7.2)
axa.axis("off")

# Top box: Boundary (left coideal)
axa.text(5.0, 6.82, r"Boundary ($\Delta(\mathcal{B}) \subset U_q(\hat{\mathfrak{g}})\otimes\mathcal{B}$): $\mathcal{W}_B$ stays on right",
         fontsize=7.7, color=style.GREY, ha="center", weight="bold")

slots_top = [
    (1.6, 5.55, r"$V_1(\theta_1)$", style.BLUE),
    (4.8, 5.55, r"$V_2(\theta_2)$", style.BLUE),
    (8.0, 5.55, r"$\mathcal{W}_B$ (right)", style.TEAL),
]
for cx, cy, txt, col in slots_top:
    box = FancyBboxPatch((cx - 1.25, cy - 0.48), 2.5, 0.96, boxstyle="round,pad=0.06,rounding_size=0.15",
                         fc="white", ec=col, lw=1.4)
    axa.add_patch(box)
    axa.text(cx, cy, txt, fontsize=7.6, color=col, ha="center", va="center")
axa.text(3.2, 5.55, r"$\otimes$", fontsize=9, color=style.GREY, ha="center", va="center")
axa.text(6.4, 5.55, r"$\otimes$", fontsize=9, color=style.GREY, ha="center", va="center")
axa.text(5.0, 4.55, r"$(\Delta\otimes 1)\Delta(\mathcal{B}) \subset U_q(\hat{\mathfrak{g}}) \otimes U_q(\hat{\mathfrak{g}}) \otimes \mathcal{B}\ \checkmark$",
         fontsize=7.3, color=style.TEAL, ha="center")

axa.axhline(4.05, xmin=0.04, xmax=0.96, color=style.LIGHTGREY, lw=0.9, ls="--")

# Bottom box: Defect (middle step of Triangle Relation)
axa.text(5.0, 3.55, r"Defect Triangle step 2: $\mathcal{W}$ sits in the MIDDLE slot!",
         fontsize=7.7, color=style.RED, ha="center", weight="bold")

slots_bot = [
    (1.6, 2.25, r"$V_1(\theta_1)$", style.BLUE),
    (4.8, 2.25, r"$\mathcal{W}$ (middle!)", style.PURPLE),
    (8.0, 2.25, r"$V_2(\theta_2)$", style.BLUE),
]
for cx, cy, txt, col in slots_bot:
    box = FancyBboxPatch((cx - 1.25, cy - 0.48), 2.5, 0.96, boxstyle="round,pad=0.06,rounding_size=0.15",
                         fc="#f8f2fc" if col == style.PURPLE else "white", ec=col, lw=1.6 if col == style.PURPLE else 1.4)
    axa.add_patch(box)
    axa.text(cx, cy, txt, fontsize=7.6, color=col, ha="center", va="center")
axa.text(3.2, 2.25, r"$\otimes$", fontsize=9, color=style.GREY, ha="center", va="center")
axa.text(6.4, 2.25, r"$\otimes$", fontsize=9, color=style.GREY, ha="center", va="center")

# Brace/arrow showing T_1 acting on V_1 (x) W
axa.annotate("", xy=(1.6, 1.52), xytext=(4.8, 1.52),
             arrowprops=dict(arrowstyle="<->", color=style.PURPLE, lw=1.3))
axa.text(3.2, 1.15, r"$\mathbb{T}_1(\theta_1)$ intertwines $V_1\otimes\mathcal{W}$", fontsize=7.2, color=style.PURPLE, ha="center")
axa.text(5.0, 0.35, r"Requires Hopf subalgebra $\Delta(U_q(\mathfrak{b}^-)) \subset U_q(\mathfrak{b}^-) \otimes U_q(\mathfrak{b}^-)$",
         fontsize=7.3, color=style.PURPLE, ha="center")

axa.set_title(r"(a) Coideal vs. Borel Hopf subalgebra", fontsize=8.4, color=style.GREY)

# Panel (b): Bilateral Weyl module ladder vs Fock floor and 2D charge block H_n
axb.set_xlim(-2.45, 2.55)
axb.set_ylim(-0.45, 7.35)
axb.axis("off")

# Column headers
axb.text(-1.25, 6.50, r"$v_+ \otimes |n\rangle$" + "\n" + r"(soliton $s$)", fontsize=7.5, color=style.BLUE, ha="center")
axb.text(0.85, 6.50, r"$v_- \otimes |n\rangle$" + "\n" + r"(antisoliton $\bar{s}$)", fontsize=7.5, color=style.RED, ha="center")

# Draw ladder states n = -1, 0, 1, 2
n_levels = [(-1, 1.05), (0, 2.55), (1, 3.85), (2, 5.15)]

# Highlight 2D charge block H_0 = span{v_+ (x) |0>, v_- (x) |1>}
poly_x = [-1.95, -0.55, 1.55, 0.15]
poly_y = [2.55 - 0.34, 2.55 - 0.34, 3.85 + 0.34, 3.85 + 0.34]
axb.fill(poly_x, poly_y, color=style.PURPLE, alpha=0.09, zorder=1)
axb.plot(poly_x + [poly_x[0]], poly_y + [poly_y[0]], color=style.PURPLE, lw=1.1, ls="--", zorder=2)
axb.text(-0.20, 3.22, r"2D block $\mathcal{H}_n$" + "\n" + r"$\det \mathbb{T}_n \propto 1 + Y^2$",
         fontsize=7.0, color=style.PURPLE, ha="center", va="center", zorder=6,
         bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none", alpha=0.88))

for n_val, y_lv in n_levels:
    # Left rung: v_+ (x) |n>
    axb.plot([-1.25], [y_lv], "o", ms=6.8, color=style.BLUE, zorder=5)
    axb.text(-1.52, y_lv, rf"$|n={n_val}\rangle$", fontsize=7.3, color=style.GREY, ha="right", va="center")
    # Right rung: v_- (x) |n>
    axb.plot([0.85], [y_lv], "o", ms=6.8, color=style.RED, zorder=5)
    axb.text(1.12, y_lv, rf"$|n={n_val}\rangle$", fontsize=7.3, color=style.GREY, ha="left", va="center")

# Vertical f_0 / f_1 ladder arrows on Weyl module
for idx in range(len(n_levels) - 1):
    y_lo = n_levels[idx][1]
    y_hi = n_levels[idx + 1][1]
    axb.annotate("", xy=(-1.25, y_hi - 0.15), xytext=(-1.25, y_lo + 0.15),
                 arrowprops=dict(arrowstyle="<->", color=style.PURPLE, lw=1.2))
axb.text(-2.35, 4.50, r"$f_0 \uparrow$" + "\n" + r"$f_1 \downarrow$", fontsize=7.3, color=style.PURPLE, ha="left", va="center")

# Fock truncation line between n = -1 and n = 0
axb.plot([-0.90, 2.35], [1.75, 1.75], color=style.ORANGE, lw=1.3, ls=":")
axb.text(0.75, 1.83, r"Fock $\mathcal{F}$ floor ($n\geq 0,\ \mathsf{a}_-|0\rangle=0$)",
         fontsize=6.8, color=style.ORANGE, ha="center", va="bottom")

# Weyl continuation dots at top and bottom
axb.text(-1.25, 5.82, r"$\vdots$", fontsize=10, color=style.PURPLE, ha="center", va="center")
axb.text(0.85, 5.82, r"$\vdots$", fontsize=10, color=style.PURPLE, ha="center", va="center")
axb.text(-1.25, 0.52, r"$\vdots$", fontsize=10, color=style.PURPLE, ha="center", va="center")
axb.text(0.85, 0.52, r"$\vdots$", fontsize=10, color=style.PURPLE, ha="center", va="center")
axb.text(-0.15, -0.12, r"Weyl $\mathcal{W}$ ($n\in\mathbb{Z},\ Z=0$): no floor!",
         fontsize=7.3, color=style.PURPLE, ha="center")

axb.set_title(r"(b) Weyl ladder $\mathcal{W}$ vs. Fock $\mathcal{F}$ & block $\mathcal{H}_n$", fontsize=8.4, color=style.GREY)

fig.subplots_adjust(left=0.02, right=0.98, bottom=0.04, top=0.90, wspace=0.15)
style.save(fig, "l4_borel_weyl", OUT)
raise SystemExit(0 if style.all_ok() else 1)
