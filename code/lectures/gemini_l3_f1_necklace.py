"""Lectures on integrable defects (defects-gemini), Lecture 3, figure 1 (fig-lec3-necklace):
The Necklace Selection Rule (thm-lec3-necklace, eq-lec3-necklace-d) and the Vanishing of the Cubic Spin-2 Obstruction
(eq-lec3-dk-eigenvalues, eq-lec3-cotangent-id).

Panel (a): Extended Dynkin diagram of d_4^(1): any diagram automorphism P fixes the central node alpha_2 (red), so
           1 - P is singular and D = 2(1 - P)^{-1} does not exist.
Panel (b): Extended Dynkin diagram of a_4^(1) (h = 5): the two-step permutation pi^{2} (orange dashed star) has
           invertible 1 - P, but 10 forbidden non-neighbour inner products alpha_i . D alpha_j != 0 spoil momentum
           conservation; only the single-step necklace cycle P = pi (green solid cycle) gives
           alpha_i . D alpha_j = 2 delta_{ij} - 2 delta_{i,j+1}.
Panel (c): Eigenvalues d_k = 1 + i cot(pi*k/h) of D = 2(1 - pi)^{-1} in the complex plane for h = 3, 4, 5, 6:
           all lie on the vertical line Re(d_k) = 1 (D + D^T = 2), and for every resonant mode triple
           k_1 + k_2 + k_3 = 0 (mod h), Re(d_{k_1} d_{k_2} d_{k_3}) = 0 so the cubic spin-2 obstruction vanishes.

Checks:
  1. Every diagram automorphism of d_4^(1) fixes the central node alpha_2, making det(1 - P) = 0.
  2. For a_4^(1), P = pi^{2} produces 10 forbidden non-zero entries in alpha_i . D alpha_j, whereas P = pi
     produces the exact matrix alpha_i . D alpha_j = 2 delta_{ij} - 2 delta_{i,j+1}.
  3. Explicit matrices E and D for a_2^(1) (Problem 3.1) satisfy D + D^T = 2 I_2, D D^T = (4/3) I_2, and have
     eigenvalues 1 +- i/sqrt(3).
  4. For h = 3, 4, 5, 6, every triple (k_1, k_2, k_3) in {1,...,h-1}^3 with k_1 + k_2 + k_3 = 0 (mod h) satisfies
     the cotangent identity cot(a)cot(b) + cot(b)cot(c) + cot(c)cot(a) = 1 and O_{k_1 k_2 k_3} = 0.

Writes lecture-notes/defects-gemini/figures/l3_necklace.svg.
"""
import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects-gemini" / "figures"

# ---- Check 1: d_4^(1) automorphisms fix central node ----
# Extended Cartan matrix of d_4^(1): node 2 is central, connected to 0, 1, 3, 4
C_d4 = np.array(
    [
        [2, 0, -1, 0, 0],
        [0, 2, -1, 0, 0],
        [-1, -1, 2, -1, -1],
        [0, 0, -1, 2, 0],
        [0, 0, -1, 0, 2],
    ],
    dtype=float,
)
aut_d4 = []
for perm in itertools.permutations(range(5)):
    P_mat = np.eye(5)[list(perm)]
    if np.allclose(P_mat @ C_d4 @ P_mat.T, C_d4):
        aut_d4.append(perm)

report(
    f"d_4^(1) has {len(aut_d4)} diagram automorphisms (= |S_4|) and every single one fixes the central node alpha_2",
    len(aut_d4) == 24 and all(p[2] == 2 for p in aut_d4),
)


# ---- Check 2: a_r^(1) inner product matrix for P = pi vs P = pi^{2} ----
def ar_roots(r):
    """Simple and lowest roots of a_r^(1) in R^r (h = r + 1)."""
    h = r + 1
    # Standard orthonormal basis e_0, ..., e_r of R^h projected to sum = 0
    E_std = np.eye(h)
    roots_h = np.array([E_std[j] - E_std[(j + 1) % h] for j in range(h)])
    # Orthonormal basis of the (h-1)-dimensional hyperplane sum x_i = 0
    Q, _ = np.linalg.qr(roots_h[:r].T)
    # Express all h roots in R^r
    return roots_h @ Q  # shape (h, r)


roots_a4 = ar_roots(4)  # shape (5, 4), simple roots are rows 0..3 (or 1..4)
# Clock shift pi: alpha_j -> alpha_{j-1}, so P = pi: alpha_j -> alpha_{j-1}
# In basis of first 4 roots (rows 0, 1, 2, 3), let's find the 4x4 matrix of P_step(s): alpha_j -> alpha_{(j+s) mod 5}
B_a4 = roots_a4[:4]  # 4x4 matrix whose rows are alpha_0, alpha_1, alpha_2, alpha_3


def d_matrix_a4(step):
    target = np.array([roots_a4[(j + step) % 5] for j in range(4)])
    # Row action: target = P_row @ B_a4, on column vectors in R^4: P_vec = B_a4.T @ np.linalg.inv(B_a4.T) ...
    # Solve (P_vec @ B_a4[j]) = target[j] for j=0..3 -> P_vec @ B_a4.T = target.T
    P_vec = target.T @ np.linalg.inv(B_a4.T)
    D_mat = 2.0 * np.linalg.inv(np.eye(4) - P_vec)
    # Inner product matrix M_{ij} = alpha_i . D alpha_j for i,j in 0..4
    M = roots_a4 @ D_mat @ roots_a4.T
    return D_mat, M


D_step1, M_step1 = d_matrix_a4(-1)
D_step2, M_step2 = d_matrix_a4(-2)

expected_M1 = np.zeros((5, 5))
for i in range(5):
    expected_M1[i, i] = 2.0
    expected_M1[i, (i - 1) % 5] = -2.0

# For step 2, allowed entries for partner p(i) = (i-2)%5 would only be (i,i)=2 and (i,(i-2)%5)=-2; count other non-zero entries
forbidden_step2 = sum(
    abs(M_step2[i, j]) > 1e-8
    for i in range(5)
    for j in range(5)
    if j != i and j != (i - 2) % 5
)

report(
    f"a_4^(1): single-step necklace cycle P=pi gives alpha_i.D alpha_j = 2 delta_ij - 2 delta_{{i,j+1}}; "
    f"two-step star P=pi^2 has {forbidden_step2} forbidden non-zero cross-terms",
    np.allclose(M_step1, expected_M1) and forbidden_step2 == 10,
)

# ---- Check 3: Problem 3.1 (Explicit E and D for a_2^(1)) ----
roots_a2 = ar_roots(2)  # shape (3, 2)
B_a2 = roots_a2[:2]
target_a2 = np.array([roots_a2[2], roots_a2[0]])
P_a2 = target_a2.T @ np.linalg.inv(B_a2.T)
D_a2 = 2.0 * np.linalg.inv(np.eye(2) - P_a2)
E_a2 = np.eye(2) - D_a2
M_a2 = roots_a2 @ D_a2 @ roots_a2.T
exp_M_a2 = np.array([[2.0, 0.0, -2.0], [-2.0, 2.0, 0.0], [0.0, -2.0, 2.0]])
eigs_a2 = np.linalg.eigvals(D_a2)

report(
    "Problem 3.1 (a_2^(1)): D + D^T = 2 I_2, D D^T = (4/3) I_2, E^T = -E, eig(D) = 1 +- i/sqrt(3), and alpha_i.D alpha_j = 2 delta_ij - 2 delta_{i,j+1}",
    np.allclose(D_a2 + D_a2.T, 2.0 * np.eye(2))
    and np.allclose(D_a2 @ D_a2.T, (4.0 / 3.0) * np.eye(2))
    and np.allclose(E_a2 + E_a2.T, 0.0)
    and np.allclose(M_a2, exp_M_a2)
    and min(abs(eigs_a2[0] - (1 + 1j / np.sqrt(3))), abs(eigs_a2[1] - (1 + 1j / np.sqrt(3)))) < 1e-12,
)

# ---- Check 4: Spin-2 cubic obstruction O_{k1 k2 k3} = 0 for h = 3, 4, 5, 6 ----
max_obs = 0.0
for h in (3, 4, 5, 6):
    dk = {k: 1.0 + 1j / np.tan(np.pi * k / h) for k in range(1, h)}
    for k1 in range(1, h):
        for k2 in range(1, h):
            for k3 in range(1, h):
                if (k1 + k2 + k3) % h == 0:
                    obs = dk[k1] * dk[k2] * dk[k3] + np.conj(dk[k1] * dk[k2] * dk[k3])
                    c1, c2, c3 = 1.0 / np.tan(np.pi * k1 / h), 1.0 / np.tan(np.pi * k2 / h), 1.0 / np.tan(np.pi * k3 / h)
                    cot_err = abs(c1 * c2 + c2 * c3 + c3 * c1 - 1.0)
                    max_obs = max(max_obs, abs(obs), cot_err)

report(
    f"spin-2 cubic obstruction O_{{k1 k2 k3}} = d_k1 d_k2 d_k3 + dbar_k1 dbar_k2 dbar_k3 vanishes for all h=3..6 [max err {max_obs:.1e}]",
    max_obs < 1e-12,
)

# ---- Figure ----
fig, (axa, axb, axc) = plt.subplots(
    1, 3, figsize=(7.0, 2.75), gridspec_kw=dict(width_ratios=[0.95, 1.05, 1.15], wspace=0.30)
)

# Panel (a): d_4^(1) tree diagram
axa.set_xlim(-1.45, 1.45)
axa.set_ylim(-1.45, 1.45)
axa.axis("off")

outer_angles = [np.pi / 4, 3 * np.pi / 4, -3 * np.pi / 4, -np.pi / 4]
outer_pts = [(0.92 * np.cos(a), 0.92 * np.sin(a)) for a in outer_angles]
outer_labs = [r"$\alpha_1$", r"$\alpha_0$", r"$\alpha_3$", r"$\alpha_4$"]

for (px, py), lab in zip(outer_pts, outer_labs):
    axa.plot([0, px], [0, py], color=style.GREY, lw=1.6, zorder=2)
    axa.plot([px], [py], "o", ms=7.5, color=style.BLUE, zorder=4)
    axa.text(1.26 * px, 1.26 * py, lab, fontsize=8, color=style.GREY, ha="center", va="center")

# Central node alpha_2 (fixed by all P)
axa.plot([0], [0], "o", ms=8.5, color=style.RED, zorder=5)
axa.text(0.0, -0.42, r"$\alpha_2$", fontsize=8.2, color=style.RED, ha="center", va="top")

# Curved permutation arrows among 3 outer nodes
for i in range(3):
    x1, y1 = outer_pts[i]
    x2, y2 = outer_pts[(i + 1) % 3]
    axa.annotate(
        "",
        xy=(x2, y2),
        xytext=(x1, y1),
        arrowprops=dict(
            arrowstyle="->",
            color=style.BLUE,
            lw=1.1,
            ls="--",
             shrinkA=7,
            shrinkB=7,
            connectionstyle="arc3,rad=0.25",
        ),
    )

axa.text(0.0, -1.35, r"$P\vec{\alpha}_2 = \vec{\alpha}_2 \Rightarrow \det(1-P)=0$", fontsize=7.6, color=style.RED, ha="center")
axa.set_title(r"(a) Tree $d_4^{(1)}$: centre fixed", fontsize=8.3, color=style.GREY)

# Panel (b): a_4^(1) pentagon (two-step star vs one-step necklace)
axb.set_xlim(-1.45, 1.45)
axb.set_ylim(-1.45, 1.45)
axb.axis("off")

pent_angles = [np.pi / 2 - 2 * np.pi * j / 5 for j in range(5)]
pent_pts = [(0.92 * np.cos(a), 0.92 * np.sin(a)) for j, a in enumerate(pent_angles)]

# Draw pentagon Dynkin bonds and green 1-step necklace arrows
for j in range(5):
    x1, y1 = pent_pts[j]
    x2, y2 = pent_pts[(j - 1) % 5]
    axb.annotate(
        "",
        xy=(x2, y2),
        xytext=(x1, y1),
        arrowprops=dict(
            arrowstyle="->",
            color=style.GREEN,
            lw=1.9,
            shrinkA=7,
            shrinkB=7,
            connectionstyle="arc3,rad=-0.12",
        ),
    )
    # Inner 2-step star (forbidden cross-terms)
    x3, y3 = pent_pts[(j - 2) % 5]
    axb.annotate(
        "",
        xy=(x3, y3),
        xytext=(x1, y1),
        arrowprops=dict(
            arrowstyle="->",
            color=style.ORANGE,
            lw=1.0,
            ls=":",
            shrinkA=8,
            shrinkB=8,
        ),
    )
    axb.plot([x1], [y1], "o", ms=7.5, color=style.GREEN, zorder=5)
    axb.text(1.25 * x1, 1.25 * y1, rf"$\alpha_{j}$", fontsize=8, color=style.GREY, ha="center", va="center")

axb.text(0.0, -1.24, r"green: $P=\pi$ (necklace $\checkmark$)", fontsize=7.4, color=style.GREEN, ha="center")
axb.text(0.0, -1.44, r"orange: $P=\pi^{2}$ (cross-terms $\times$)", fontsize=7.4, color=style.ORANGE, ha="center")
axb.set_title(r"(b) Necklace $a_4^{(1)}$ ($h=5$)", fontsize=8.3, color=style.GREY)

# Panel (c): Eigenvalues d_k = 1 + i cot(pi k / h) in the complex plane
axc.axhline(0.0, color=style.LIGHTGREY, lw=0.8)
axc.axvline(1.0, color=style.PURPLE, lw=1.4, ls="--", label=r"$\operatorname{Re} d_k = 1$")

series = [
    (3, "o", style.BLUE, r"$a_2^{(1)}\ (h=3)$"),
    (4, "s", style.TEAL, r"$a_3^{(1)}\ (h=4)$"),
    (5, "^", style.ORANGE, r"$a_4^{(1)}\ (h=5)$"),
    (6, "D", style.RED, r"$a_5^{(1)}\ (h=6)$"),
]
for h_val, mk, col, lab in series:
    ks = np.arange(1, h_val)
    im_vals = 1.0 / np.tan(np.pi * ks / h_val)
    axc.plot(np.ones_like(im_vals), im_vals, mk, color=col, ms=4.8, label=lab, zorder=5)

axc.text(1.10, 1.48, r"$d_k = 1 + i\cot\frac{\pi k}{h}$", fontsize=7.6, color=style.PURPLE, ha="left")
axc.text(
    1.10,
    -1.45,
    r"$k_1+k_2+k_3 \equiv 0$" + "\n" + r"$\Rightarrow \operatorname{Re}(d_{k_1}d_{k_2}d_{k_3})=0$",
    fontsize=7.0,
    color=style.GREY,
    ha="left",
)

axc.set_xlim(0.0, 2.38)
axc.set_ylim(-1.95, 1.95)
axc.set_xlabel(r"$\operatorname{Re} d_k$")
axc.set_ylabel(r"$\operatorname{Im} d_k = \cot(\pi k / h)$")
axc.set_title(r"(c) Spectrum of $\mathsf{D}$ & spin-2 survival", fontsize=8.3, color=style.GREY)
axc.legend(loc="upper left", fontsize=6.9, handlelength=1.2, borderpad=0.2)

fig.subplots_adjust(left=0.03, right=0.98, bottom=0.17, top=0.88, wspace=0.28)
style.save(fig, "l3_necklace", OUT)
raise SystemExit(0 if style.all_ok() else 1)
