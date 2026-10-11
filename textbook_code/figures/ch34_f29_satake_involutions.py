"""Figure F29 (chapter 34, sec-symmetric-pairs, sec-preserving-conjugating): the two Satake diagrams of the physical
reflection matrices of a_5^(1), (X, tau) = (empty, id) and (empty, opposition).

(empty, id) is the Chevalley involution: Theta = -1 on weights, the multiplet V_a goes to its conjugate V_{h-a}
(soliton-conjugating boundary).  (empty, opposition) is the affine extension of j -> n+1-j, fixing node 0: on weights
Theta = -tau = reversal of the coordinates eps_i -> eps_{n+2-i}, V_a goes to V_a (soliton-preserving boundary).

Checks: the opposition permutation of the nodes of a_5^(1) preserves the Cartan matrix, squares to 1 and fixes node 0 and node 3 (for X = empty, Kolb's admissibility reduces to tau^2 = 1);
the weight multiset of Lambda^a C^6 is mapped to that of Lambda^{6-a} by Theta = -1
(and to itself only for a = 3) and to itself by the coordinate reversal, for all a.

Writes part7-half-line/figures/ch34_f29_satake_involutions.svg.
"""
import itertools

import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch

import style
from style import plt, report

style.setup()
n = 5
h = n + 1
# ---- checks ----------------------------------------------------------------
C = 2 * np.eye(h, dtype=int)
for j in range(h):
    C[j, (j + 1) % h] -= 1
    C[(j + 1) % h, j] -= 1
tau = [(-j) % h for j in range(h)]                    # j -> n+1-j for j >= 1, node 0 fixed
P = np.eye(h, dtype=int)[tau]
report("the opposition permutation tau(j) = -j mod 6 preserves the affine Cartan matrix of a_5^(1), tau^2 = 1", np.array_equal(P @ C @ P.T, C) and all(tau[tau[j]] == j for j in range(h)))
report("tau fixes node 0 (and node 3); pairs (1,5), (2,4)", tau[0] == 0 and tau[3] == 3 and tau[1] == 5 and tau[2] == 4)


def weights(a):
    """weights of Lambda^a C^6 as 0/1 vectors in the eps-basis, projected to the sum-zero hyperplane (times 6 to stay integer)"""
    return {tuple(h * np.array(w) - a) for w in (np.eye(h, dtype=int)[list(S)].sum(0) for S in itertools.combinations(range(h), a))}


ok_conj = ok_pres = True
for a in range(1, h):
    W = weights(a)
    neg = {tuple(-x for x in w) for w in W}
    rev = {tuple(w[::-1]) for w in W}                 # Theta = -tau: eps_i -> eps_{n+2-i}
    ok_conj &= (neg == weights(h - a)) and ((neg == W) == (a == h // 2))
    ok_pres &= rev == W
report("Theta = -1 (tau = id) maps the weights of Lambda^a to those of Lambda^{6-a}, equal to Lambda^a only for a = 3", ok_conj)
report("Theta = -tau (tau = opposition: eps_i -> eps_{7-i}) maps the weights of Lambda^a to themselves, a = 1..5", ok_pres)

# ---- figure ----------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.3))
ang = [np.pi / 2 - 2 * np.pi * j / h for j in range(h)]
pos = [(1.0 * np.cos(t), 1.0 * np.sin(t)) for t in ang]
panels = (("$\\tau=\\mathrm{id}$:  soliton-conjugating", False, r"$\Theta=-1$:   $V_a\mapsto V_{h-a}$"),
          (r"$\tau=$ opposition:  soliton-preserving", True, r"$\Theta=-\tau$:   $V_a\mapsto V_a$"))
for ax, (title, arrows, foot) in zip(axs, panels):
    ax.set_aspect("equal")
    ax.axis("off")
    for j in range(h):
        (x0, y0), (x1, y1) = pos[j], pos[(j + 1) % h]
        dx, dy = x1 - x0, y1 - y0
        L_ = np.hypot(dx, dy)
        ax.plot([x0 + 0.15 * dx / L_, x1 - 0.15 * dx / L_], [y0 + 0.15 * dy / L_, y1 - 0.15 * dy / L_], color=style.GREY, lw=1.4, zorder=1)
    for j in range(h):
        x, y = pos[j]
        ax.add_patch(Circle((x, y), 0.15, fc="none" if j else style.LIGHTGREY, ec=style.GREY, lw=1.4, zorder=3))
        ax.text(1.0 * x * 1.38, y * 1.38, str(j), ha="center", va="center", fontsize=8.5)
    if arrows:
        for j, k, rad in ((1, 5, 0.0), (2, 4, 0.0)):
            (x0, y0), (x1, y1) = pos[j], pos[k]
            ax.add_patch(FancyArrowPatch((x0 * 0.82, y0), (x1 * 0.82, y1), arrowstyle="<->", mutation_scale=9, color=style.RED, lw=1.5, zorder=4,
                                          connectionstyle="arc3,rad=0"))
        ax.text(0, 0.0, "nodes 0 and 3 fixed", fontsize=7, ha="center", va="center", color=style.RED)
        ax.text(0, -0.32, "", fontsize=6)
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.95, 1.7)
    ax.set_title(title, fontsize=8.5, pad=2)
    ax.text(0, -1.8, foot, ha="center", va="center", fontsize=8)
fig.tight_layout()
style.save(fig, "ch34_f29_satake_involutions", style.OUT_PART7)
raise SystemExit(0 if style.all_ok() else 1)
