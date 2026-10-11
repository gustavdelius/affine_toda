"""Figure F11 (chapter 6, sec-affine-folding): folding affine Dynkin diagrams.

Top row: a3^(1) -> c2^(1) (reflection fixing nodes 0 and 2).  Bottom row: d4^(1) -> g2^(1) (triality, permuting
three end nodes).  Nodes are labelled by their Kac labels; the roots of one orbit share a colour.

Checks (from explicit roots, not typed in): sum_j n_j alpha_j = 0 in the parent; the folded roots alpha_O
(orbit averages) have squared lengths 1 and 2/3 for the short orbits and 2 otherwise; the folded labels |O| n_O
are 1,1,2 and 1,2,3 and sum to h = 4 and 6; the Gram matrix of the folded roots gives the affine Cartan matrices of
c2^(1) and g2^(1) and the labels are its null vector.
Writes part2-algebra/figures/ch6_f11_folding.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()


def affine_roots(edges, marks, r):
    """Simply-laced finite diagram on nodes 1..r (edges), highest root with marks; returns the r+1 affine roots
    as rows of a matrix in R^r (Cholesky of the Cartan matrix); alpha_0 = -theta."""
    C = 2 * np.eye(r)
    for i, j in edges:
        C[i - 1, j - 1] = C[j - 1, i - 1] = -1
    S = np.linalg.cholesky(C)
    a0 = -np.array(marks) @ S
    return np.vstack([a0, S])


def fold(roots, orbits):
    return np.array([roots[list(o)].mean(0) for o in orbits])


# a3^(1): finite nodes 1-2-3 (node 2 in the middle), theta = a1+a2+a3; the affine node 0 closes the square.
A3 = affine_roots([(1, 2), (2, 3)], [1, 1, 1], 3)
marks3 = np.array([1, 1, 1, 1])
report("a3^(1): sum_j n_j alpha_j = 0, all alpha_j^2 = 2", np.allclose(marks3 @ A3, 0) and np.allclose((A3**2).sum(1), 2))
orb3 = [(0,), (2,), (1, 3)]                                   # fixes 0 and 2, exchanges 1 and 3
F3 = fold(A3, orb3)
lab3 = np.array([len(o) * marks3[list(o)[0]] for o in orb3])
sq3 = (F3**2).sum(1)
report("a3 -> c2: folded roots have squared lengths 2, 2, 1 and labels 1, 1, 2 (sum h = 4)",
       np.allclose(sq3, [2, 2, 1]) and list(lab3) == [1, 1, 2] and lab3.sum() == 4)
G3 = F3 @ F3.T
Cf3 = 2 * G3 / np.diag(G3)[:, None]                              # a_ij = 2 alpha_i . alpha_j / alpha_i^2
report("c2^(1): Cartan matrix of the diagram 0 => {1,3} <= 2 (nodes ordered 0, 2, {1,3}); labels are its null vector",
       np.allclose(Cf3, [[2, 0, -1], [0, 2, -1], [-2, -2, 2]]) and np.allclose(Cf3 @ lab3, 0))
report("c2^(1): sum_O |O| n_O alpha_O = 0", np.allclose(lab3 @ F3, 0))

# d4^(1): finite nodes 1,3,4 attached to the central node 2, theta = a1 + 2 a2 + a3 + a4; node 0 attached to 2.
D4 = affine_roots([(1, 2), (2, 3), (2, 4)], [1, 2, 1, 1], 4)
marks4 = np.array([1, 1, 2, 1, 1])
report("d4^(1): sum_j n_j alpha_j = 0, all alpha_j^2 = 2, h = 6",
       np.allclose(marks4 @ D4, 0) and np.allclose((D4**2).sum(1), 2) and marks4.sum() == 6)
orb4 = [(0,), (2,), (1, 3, 4)]                                # triality fixing 0 and the centre
F4 = fold(D4, orb4)
lab4 = np.array([len(o) * marks4[list(o)[0]] for o in orb4])
sq4 = (F4**2).sum(1)
report("d4 -> g2: folded squared lengths 2, 2, 2/3 and labels 1, 2, 3 (sum h = 6)",
       np.allclose(sq4, [2, 2, 2 / 3]) and list(lab4) == [1, 2, 3] and lab4.sum() == 6)
G4 = F4 @ F4.T
Cf4 = 2 * G4 / np.diag(G4)[:, None]
report("g2^(1): Cartan matrix of the diagram 0 - 2 <== {1,3,4} (nodes ordered 0, 2, {1,3,4}); labels are its null vector",
       np.allclose(Cf4, [[2, -1, 0], [-1, 2, -1], [0, -3, 2]]) and np.allclose(Cf4 @ lab4, 0))
report("g2^(1): sum_O |O| n_O alpha_O = 0", np.allclose(lab4 @ F4, 0))

# ---- figure ---------------------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.0, 6.0))
ax.set_xlim(0, 10)
ax.set_ylim(0.2, 8.8)
ax.axis("off")
ax.set_aspect("equal")
RAD = 0.3


def node(x, y, kac, col, idx=None, above=False):
    ax.add_patch(plt.Circle((x, y), RAD, fc=col, ec=col, lw=1.2, zorder=3))
    ax.text(x, y, kac, ha="center", va="center", fontsize=9.5, color="white", zorder=4, fontweight="bold")
    if idx is not None:
        ax.text(x, y + (RAD + 0.12 if above else -RAD - 0.12), idx, ha="center", va="bottom" if above else "top",
                fontsize=8, color=style.GREY)


def bond(p, q, mult=1, arrow=False):
    """Bond from p to q; an arrowhead points towards q (the shorter root) when arrow is set."""
    p, q = np.array(p, float), np.array(q, float)
    d = (q - p) / np.linalg.norm(q - p)
    n = np.array([-d[1], d[0]])
    a, b = p + d * RAD, q - d * RAD
    for o in {1: [0], 2: [-0.07, 0.07], 3: [-0.12, 0, 0.12]}[mult]:
        ax.plot(*zip(a + n * o, b + n * o), color=style.GREY, lw=1.3, zorder=1)
    if arrow:
        mid = (a + b) / 2
        ax.annotate("", xy=mid + d * 0.15, xytext=mid - d * 0.15,
                    arrowprops=dict(arrowstyle="-|>", color=style.GREY, lw=1.3, mutation_scale=15,
                                    shrinkA=0, shrinkB=0))


def fold_arrow(y):
    ax.annotate("", xy=(5.4, y), xytext=(3.9, y), arrowprops=dict(arrowstyle="-|>", color=style.GREY, lw=1.4))
    ax.text(4.65, y + 0.15, "fold", ha="center", va="bottom", fontsize=8.5, color=style.GREY)


def title(x, y, t):
    ax.text(x, y, t, ha="center", va="center", fontsize=10.5, color=style.GREY)


FIX, ORB, MID = style.BLUE, style.ORANGE, style.PURPLE

# --- a3^(1) -> c2^(1)
sq = {0: (1.2, 7.5), 1: (2.8, 7.5), 2: (2.8, 5.9), 3: (1.2, 5.9)}
for i, j in [(0, 1), (1, 2), (2, 3), (3, 0)]:
    bond(sq[i], sq[j])
for j, col in ((0, FIX), (1, ORB), (2, FIX), (3, ORB)):
    node(*sq[j], "1", col, r"$j=%d$" % j, above=j in (0, 1))
title(2.0, 8.6, r"$a_3^{(1)}$")
fold_arrow(6.7)
y, xs = 6.7, [6.3, 7.8, 9.3]
bond((xs[0], y), (xs[1], y), 2, arrow=True)
bond((xs[2], y), (xs[1], y), 2, arrow=True)
node(xs[0], y, "1", FIX, r"$\{0\}$")
node(xs[1], y, "2", ORB, r"$\{1,3\}$")
node(xs[2], y, "1", FIX, r"$\{2\}$")
title(7.8, 8.6, r"$c_2^{(1)}$")
ax.text(7.8, 5.25, r"$\vec\alpha_{\{1,3\}}^2=1$, others $2$", ha="center", va="center", fontsize=8, color=style.GREY)

# --- d4^(1) -> g2^(1)
cx, cy = 2.0, 2.3
legs = {0: (0.9, 3.2), 1: (3.1, 3.2), 3: (0.9, 1.4), 4: (3.1, 1.4)}
for p in legs.values():
    bond((cx, cy), p)
node(cx, cy, "2", MID, r"$j=2$")
node(*legs[0], "1", FIX, r"$j=0$", above=True)
for k in (1, 3, 4):
    node(*legs[k], "1", ORB, r"$j=%d$" % k, above=k == 1)
title(2.0, 4.3, r"$d_4^{(1)}$")
fold_arrow(2.3)
y = 2.3
bond((xs[0], y), (xs[1], y))
bond((xs[1], y), (xs[2], y), 3, arrow=True)
node(xs[0], y, "1", FIX, r"$\{0\}$")
node(xs[1], y, "2", MID, r"$\{2\}$")
node(xs[2], y, "3", ORB, r"$\{1,3,4\}$")
title(7.8, 4.3, r"$g_2^{(1)}$")
ax.text(7.8, 0.9, r"$\vec\alpha_{\{1,3,4\}}^2=\frac{2}{3}$, others $2$", ha="center", va="center", fontsize=8,
        color=style.GREY)
ax.text(0.2, 0.45, r"number in node: Kac label $n_j$ (parent), $|O|\,n_O$ (folded); colour = orbit", ha="left",
        va="center", fontsize=7.5, color=style.GREY)
style.save(fig, "ch6_f11_folding", style.OUT_PART2)
raise SystemExit(0 if style.all_ok() else 1)
