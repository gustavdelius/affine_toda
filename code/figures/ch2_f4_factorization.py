"""Figure F4 (chapter 2, sec-factorization and sec-yang-baxter): factorization.

Three particles with velocities v1 > v2 > v3 (rapidity ordering theta_1 > theta_2 > theta_3)
move on straight world lines; an elastic factorized collision does not change them, it
only multiplies by a two-particle S-matrix. Shifting the middle trajectory (which the
momentum-dependent translations e^{i a Q_s} do) changes the order of the pairwise
collisions: 12,13,23 versus 23,13,12. Independence of this order is the Yang-Baxter
equation S12 S13 S23 = S23 S13 S12 (eq-yang-baxter).

Writes part1-background/figures/ch2_f4_factorization.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
v = np.array([1.0, 0.0, -1.0])             # v1 > v2 > v3
colors = [style.BLUE, style.ORANGE, style.GREEN]


def collisions(x0):
    out = {}
    for i in range(3):
        for j in range(i + 1, 3):
            tij = (x0[j] - x0[i]) / (v[i] - v[j])     # x_i(t) = x_j(t)
            out[(i + 1, j + 1)] = (tij, x0[i] + v[i] * tij)
    return out


cfgA = np.array([-3.0, -1.0, 3.0])
cfgB = np.array([-3.0, +1.0, 3.0])
cA, cB = collisions(cfgA), collisions(cfgB)
orderA = [k for k, _ in sorted(cA.items(), key=lambda kv: kv[1][0])]
orderB = [k for k, _ in sorted(cB.items(), key=lambda kv: kv[1][0])]
report("configuration (a): collision order 12, 13, 23", orderA == [(1, 2), (1, 3), (2, 3)])
report("configuration (b): collision order 23, 13, 12", orderB == [(2, 3), (1, 3), (1, 2)])
report("same three pairs in both, all collisions at positive time",
       set(orderA) == set(orderB) and all(t > 0 for t, _ in list(cA.values()) + list(cB.values())))
report("the 1-3 collision is unchanged by shifting particle 2", cA[(1, 3)] == cB[(1, 3)])

fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.0), sharey=True)
tt = np.array([0.0, 6.0])
for ax, x0, c, order, title in (
        (axes[0], cfgA, cA, orderA, r"(a) collisions in the order $12,\,13,\,23$"),
        (axes[1], cfgB, cB, orderB, r"(b) the order $23,\,13,\,12$")):
    for i in range(3):
        ax.plot(x0[i] + v[i] * tt, tt, color=colors[i])
        ax.text(x0[i] + v[i] * 0.0, -0.45, rf"${i + 1}$", color=colors[i], ha="center", fontsize=10)
    for (i, j), (t_, x_) in c.items():
        ax.plot([x_], [t_], "o", color="white", mec=style.GREY, ms=7, zorder=5)
        ax.text(x_ + 0.28, t_ + 0.05, rf"$S_{{{i}{j}}}$", fontsize=8.5)
    ax.set_title(title, fontsize=9, color=style.GREY)
    ax.set_xlim(-4.5, 4.5)
    ax.set_ylim(-0.9, 6.2)
    ax.set_xlabel(r"$x$")
    ax.set_xticks([])
axes[0].set_ylabel(r"$t$")
axes[0].set_yticks([])
fig.text(0.5, -0.04,
         r"$S_{12}S_{13}S_{23}=S_{23}S_{13}S_{12}$:  the amplitude cannot depend on the order",
         ha="center", fontsize=9, color=style.GREY)
fig.tight_layout()
style.save(fig, "ch2_f4_factorization", style.OUT_PART1)
raise SystemExit(0 if style.all_ok() else 1)
