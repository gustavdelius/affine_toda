"""Figure F14 (chapter 10, sec-topological-charges): topological charges of single a_n^(1) solitons.

Charges q = beta phi(+inf)/(2 pi) = -(1/2 pi) sum_j alpha_j [theta_j]_(-pi,pi), theta_j = 2 pi j a/h + Im xi, are
computed by sweeping Im xi over a full period.  Panels: the four weights of Lambda^1 C^4 (a tetrahedron), all
reached by the species-1 soliton, and the six weights of Lambda^2 C^4 (an octahedron), of which the species-2
soliton reaches only two.

Checks: for a_3, a_4, a_5 and every species a, the number of distinct charges equals h/gcd(a,h); every charge is a
weight e_S - (a/h)(1,...,1) of Lambda^a C^h (|S| = a), up to an overall sign convention printed below; the
number of weights is binom(h,a); a_2 reaches all 3 weights, a_3 species 2 reaches 2 of 6.
Writes part3-classical/figures/ch10_f14_topological_charges.svg.
"""
import itertools
from math import comb, gcd

import numpy as np

import style
from style import plt, report

style.setup()


def charges(h, a, n=2000):
    """Distinct charges (as rounded tuples) and the order in which they appear as Im xi runs over (0, 2 pi)."""
    out, order = {}, []
    for im in (np.arange(n) + 0.5) * 2 * np.pi / n + 1e-4:
        th = 2 * np.pi * np.arange(h) * a / h + im
        c = (th + np.pi) % (2 * np.pi) - np.pi                  # reduced to (-pi, pi)
        # alpha_j = e_j - e_{j+1} (indices mod h), j = 0..h-1  ->  coefficient of e_j is c_j - c_{j-1}
        q = -(c - np.roll(c, 1)) / (2 * np.pi)
        k = tuple(np.round(q, 6) + 0.0)
        if k not in out:
            out[k] = np.array(q)
            order.append(k)
    return out, order


def weights(h, a, sign=1):
    ws = []
    for S in itertools.combinations(range(h), a):
        w = np.full(h, -a / h)
        w[list(S)] += 1
        ws.append(tuple(np.round(sign * w, 6) + 0.0))
    return ws


ok_count = ok_wt = True
sign_seen = set()
for h in (3, 4, 5, 6):
    for a in range(1, h):
        ch, order = charges(h, a)
        ok_count &= len(ch) == h // gcd(a, h)
        wp, wm = set(weights(h, a, 1)), set(weights(h, a, -1))
        if set(ch) <= wp:
            sign_seen.add("+")
        elif set(ch) <= wm:
            sign_seen.add("-")
        else:
            ok_wt = False
report("a_2..a_5: number of distinct charges of species a is h/gcd(a,h)", ok_count)
report("every charge is a weight of Lambda^a C^h (sign %s)" % "".join(sorted(sign_seen)), ok_wt and len(sign_seen) == 1)
sgn = 1 if "+" in sign_seen else -1
report("a_2: the species-1 soliton reaches all 3 = binom(3,1) weights", len(charges(3, 1)[0]) == comb(3, 1))
report("a_3: species 1 reaches 4 of 4, species 2 reaches 2 of binom(4,2) = 6 weights",
       len(charges(4, 1)[0]) == 4 == comb(4, 1) and len(charges(4, 2)[0]) == 2 and comb(4, 2) == 6)

# ---- figure ----
f = np.array([[1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]]) / 2.0    # orthonormal basis of the weight space of a_3
ay, ax_ = np.radians(28), np.radians(18)
Ry = np.array([[np.cos(ay), 0, np.sin(ay)], [0, 1, 0], [-np.sin(ay), 0, np.cos(ay)]])
Rx = np.array([[1, 0, 0], [0, np.cos(ax_), -np.sin(ax_)], [0, np.sin(ax_), np.cos(ax_)]])
view = lambda w: (Rx @ Ry @ (f @ np.array(w)))[:2]

fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.3))
for ax, a, title in ((axs[0], 1, r"species 1: $\Lambda^1\mathbb{C}^4$, 4 weights, all reached"),
                     (axs[1], 2, r"species 2: $\Lambda^2\mathbb{C}^4$, 6 weights, 2 reached")):
    ws = weights(4, a, sgn)
    ch, order = charges(4, a)
    P = {w: view(w) for w in ws}
    for w1, w2 in itertools.combinations(ws, 2):
        if np.isclose(np.sum((np.array(w1) - np.array(w2)) ** 2), 2):     # weights differ by one swap
            ax.plot(*zip(P[w1], P[w2]), color=style.GREY, lw=0.7, alpha=0.6)
    for w in ws:
        hit = w in ch
        ax.scatter(*P[w], s=70 if hit else 40, facecolor=style.RED if hit else "none",
                   edgecolor=style.RED if hit else style.GREY, lw=1.4, zorder=3)
    if len(order) == 2:                                                    # the two reached charges alternate with Im xi
        ax.annotate("", xy=P[order[1]], xytext=P[order[0]],
                    arrowprops=dict(arrowstyle="<->", color=style.RED, lw=1.4, shrinkA=6, shrinkB=6))
    ax.set_title(title, fontsize=8.5)
    ax.set_aspect("equal")
    ax.axis("off")
fig.tight_layout()
style.save(fig, "ch10_f14_topological_charges", style.OUT_PART3)
raise SystemExit(0 if style.all_ok() else 1)
