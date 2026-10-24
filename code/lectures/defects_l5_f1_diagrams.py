"""Lecture notes on defects, Lecture 5, figure 1 (fig-l5-diagrams): why only the cycle has type-I defects.

Four affine Dynkin diagrams with a symmetry p of their nodes, drawn as arrows from i to p(i):
(a) d_4^(1) with a 3-cycle of its outer nodes: like every tree, it has a centre that every symmetry fixes;
(b) a_3^(1) with p = pi^2: two orbits, so alpha_0 + alpha_2 is fixed and 1 - P is singular;
(c) a_4^(1) with p = pi^2: D = 2 (1 - P)^{-1} exists, but the partners are not neighbours, and forbidden
    cross terms alpha_i . D alpha_j != 0, j not in {i, p(i)}, remain;
(d) a_4^(1) with p = pi, the rotation by one step: partners are neighbours, and alpha_i . D alpha_j =
    2 delta_ij - 2 delta_{i,j+1}.
Here pi alpha_j = alpha_{j-1}, so the arrows of pi point from j to j - 1.

Checks: every automorphism of d_4^(1) fixes the centre; det(1 - P) = 0 and P fixes alpha_0 + alpha_2 for
a_3^(1), pi^2; the number of forbidden cross terms for a_4^(1), pi^2 (ten); the entries of alpha_i . D alpha_j for
a_4^(1), pi. Writes lecture-notes/defects/figures/l5_diagrams.svg.
"""
import itertools
import math
import sys
from pathlib import Path

import numpy as np
from matplotlib.patches import FancyArrowPatch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects" / "figures"


def cycle_roots(h):
    """alpha_j = e_j - e_{j+1} in an orthonormal basis of the hyperplane sum = 0 of R^h."""
    U = np.linalg.svd(np.eye(h) - np.ones((h, h)) / h)[0][:, :h - 1]
    e = np.eye(h)
    return [U.T @ (e[j] - e[(j + 1) % h]) for j in range(h)]


def perm_matrix(al, p):
    A = np.array(al[1:]).T
    Ap = np.array([al[p[i]] for i in range(1, len(al))]).T
    return Ap @ np.linalg.inv(A)


# (a) d_4^(1): Gram matrix with centre 2
Gd = 2 * np.eye(5)
for i in (0, 1, 3, 4):
    Gd[i, 2] = Gd[2, i] = -1
auts = [p for p in itertools.permutations(range(5)) if all(Gd[p[i], p[j]] == Gd[i, j] for i in range(5)
                                                               for j in range(5))]
report(f"d_4^(1): all {len(auts)} automorphisms fix the centre node", len(auts) == 24 and all(p[2] == 2 for p in auts))

# (b) a_3^(1), pi^2
al3 = cycle_roots(4)
p2 = tuple((j - 2) % 4 for j in range(4))
P = perm_matrix(al3, p2)
report("a_3^(1), pi^2: det(1 - P) = 0 and P fixes alpha_0 + alpha_2",
       abs(np.linalg.det(np.eye(3) - P)) < 1e-12 and np.allclose(P @ (al3[0] + al3[2]), al3[0] + al3[2]))

# (c), (d) a_4^(1)
al4 = cycle_roots(5)
forbidden = {}
for k in (1, 2):
    p = tuple((j - k) % 5 for j in range(5))
    P = perm_matrix(al4, p)
    D = 2 * np.linalg.inv(np.eye(4) - P)
    G = np.array([[al4[i] @ D @ al4[j] for j in range(5)] for i in range(5)])
    forbidden[k] = sum(1 for i in range(5) for j in range(5) if j not in (i, p[i]) and abs(G[i, j]) > 1e-9)
    if k == 1:
        ok_d = all(abs(G[i, j] - (2 * (i == j) - 2 * (i == (j + 1) % 5))) < 1e-12 for i in range(5) for j in range(5))
report(f"a_4^(1), pi^2: {forbidden[2]} forbidden cross terms alpha_i.D alpha_j, j not in (i, p(i))", forbidden[2] == 10)
report("a_4^(1), pi: no forbidden cross terms, alpha_i.D alpha_j = 2 delta_ij - 2 delta_{i,j+1}",
       forbidden[1] == 0 and ok_d)

# ---------------------------------------------------------------------------------------------------------
fig, axes = plt.subplots(1, 4, figsize=(7.0, 2.35))
NODE = dict(s=150, facecolor="white", edgecolor=style.GREY, linewidth=1.2, zorder=3)


def draw(ax, pos, edges, arrows, fixed=(), arrow_col=style.BLUE, rad=0.25, title=""):
    for i, j in edges:
        ax.plot(*zip(pos[i], pos[j]), color=style.GREY, lw=1.1, zorder=1)
    for i, j in arrows:
        a = FancyArrowPatch(pos[i], pos[j], connectionstyle=f"arc3,rad={rad}", arrowstyle="-|>",
                            mutation_scale=8, color=arrow_col, lw=1.0, shrinkA=7, shrinkB=7, zorder=2)
        ax.add_patch(a)
    xs, ys = zip(*pos.values())
    ax.scatter(xs, ys, **NODE)
    for i in fixed:
        ax.scatter(*pos[i], s=150, facecolor=style.RED, edgecolor=style.RED, zorder=4)
    for i, (x, y) in pos.items():
        ax.text(x, y, str(i), ha="center", va="center", fontsize=7, zorder=5,
                color="white" if i in fixed else style.GREY)
    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.45, 1.45)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=7.5, color=style.GREY, pad=2)


def polygon(h):
    return {j: (math.sin(2 * math.pi * j / h), math.cos(2 * math.pi * j / h)) for j in range(h)}


# (a)
posd = {2: (0.0, 0.0), 0: (-0.85, 0.85), 1: (0.85, 0.85), 3: (0.85, -0.85), 4: (-0.85, -0.85)}
draw(axes[0], posd, [(0, 2), (1, 2), (3, 2), (4, 2)], [(0, 1), (1, 3), (3, 0)], fixed=(2,), rad=-0.3,
     title="(a) $d_4^{(1)}$: the centre is fixed")
# (b)
pos4 = polygon(4)
draw(axes[1], pos4, [(j, (j + 1) % 4) for j in range(4)], [(j, (j - 2) % 4) for j in range(4)], rad=0.2,
     title=r"(b) $a_3^{(1)}$, $\pi^2$: two orbits")
# (c)
pos5 = polygon(5)
draw(axes[2], pos5, [(j, (j + 1) % 5) for j in range(5)], [(j, (j - 2) % 5) for j in range(5)], rad=0.12,
     arrow_col=style.ORANGE, title=r"(c) $a_4^{(1)}$, $\pi^2$: cross terms")
# (d)
draw(axes[3], pos5, [(j, (j + 1) % 5) for j in range(5)], [(j, (j - 1) % 5) for j in range(5)], rad=-0.45,
     arrow_col=style.GREEN, title=r"(d) $a_4^{(1)}$, $\pi$: a defect")
fig.tight_layout(w_pad=0.3)
style.save(fig, "l5_diagrams", OUT)
raise SystemExit(0 if style.all_ok() else 1)
