"""Figure F36 (chapter 28, sec-boundary-unitarity-crossing, eq-boundary-state): the boundary state.

Left: the Euclidean half-plane x <= 0 with Euclidean time y = it upwards; the boundary is a line at x = 0 and a
particle reflects off it.  Right: the same plane read the other way, with y as space and tau_E = -x >= 0 as
Euclidean time.  The boundary is now an initial condition at tau_E = 0, the boundary state |B>, and the reflected
world line becomes a pair of particles with opposite rapidities emitted from it (Ghoshal-Zamolodchikov, sec. 3).

Checks: the map (x, y) -> (y, -x) sends the half-plane x <= 0 to tau_E >= 0 and the boundary x = 0 to tau_E = 0;
it sends the incoming and outgoing halves of a reflected world line (slopes dy/dx = -c and +c) to two lines from
one boundary point with slopes d tau_E/dy = +1/c and -1/c, i.e. a pair with opposite rapidities.

Writes part7-half-line/figures/ch28_f36_boundary_state.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()


def rot(x, y):
    return y, -x


xs = np.linspace(-2.5, 0, 7)
c = 1.4                                            # |dy/dx| of the Euclidean world line

inc = (xs, c * xs)
out = (xs, -c * xs)                                # outgoing: y = -c x
ri, ro = rot(*inc), rot(*out)
report("the half-plane x <= 0 maps to tau_E >= 0 and x = 0 to tau_E = 0",
       all(rot(x, 0.3)[1] >= 0 for x in xs) and rot(0.0, 1.7)[1] == 0)
si = np.polyfit(ri[0], ri[1], 1)[0]
so = np.polyfit(ro[0], ro[1], 1)[0]
report("the reflected world line becomes a pair from one boundary point with opposite slopes (opposite rapidities)",
       abs(si + so) < 1e-12 and abs(abs(si) - 1 / c) < 1e-12 and rot(0.0, 0.0) == (0.0, -0.0))

fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.6, 2.9))
# left: x horizontal (x <= 0), y vertical
a1.axvspan(0, 0.8, color=style.LIGHTGREY, lw=0)
a1.axvline(0, color=style.GREY, lw=1.6)
a1.plot(*inc, color=style.BLUE, lw=1.6)
a1.plot(*out, color=style.RED, lw=1.6)
a1.annotate("", xy=(-3.3, 2.6), xytext=(-3.3, 0.6), arrowprops=dict(arrowstyle="->", color=style.GREY))
a1.text(-3.2, 2.7, r"time $y=it$", fontsize=8, color=style.GREY)
a1.text(0.08, -3.6, "boundary\ncondition", fontsize=8, color=style.GREY)
a1.text(-2.2, -3.6, r"$A(\theta)$ in", fontsize=8, color=style.BLUE)
a1.text(-2.2, 3.3, r"$A(-\theta)$ out", fontsize=8, color=style.RED)
a1.set_xlim(-3.6, 0.8)
a1.set_ylim(-4.2, 4.2)
a1.set_xlabel(r"space $x\leq0$")
a1.set_title("open channel: the boundary is a place", fontsize=9, color=style.GREY)
# right: y horizontal, tau_E = -x vertical
a2.axhspan(-0.8, 0, color=style.LIGHTGREY, lw=0)
a2.axhline(0, color=style.GREY, lw=1.6)
a2.plot(*ri, color=style.BLUE, lw=1.6)
a2.plot(*ro, color=style.RED, lw=1.6)
a2.text(3.0, -0.6, r"$|B\rangle$", fontsize=10, color=style.GREY)
a2.annotate("", xy=(-4.0, 3.3), xytext=(-4.0, 1.3), arrowprops=dict(arrowstyle="->", color=style.GREY))
a2.text(-3.9, 3.3, r"time $\tau_E=-x$", fontsize=8, color=style.GREY)
a2.text(-4.1, 1.0, r"$A(-\theta)$", fontsize=8, color=style.BLUE)
a2.text(2.7, 1.0, r"$A(\theta)$", fontsize=8, color=style.RED)
a2.set_xlim(-4.6, 4.6)
a2.set_ylim(-0.8, 3.6)
a2.set_xlabel(r"space $y$")
a2.set_title("closed channel: the boundary is a state", fontsize=9, color=style.GREY)
for ax in (a1, a2):
    ax.set_xticks([])
    ax.set_yticks([])
fig.tight_layout()
style.save(fig, "ch28_f36_boundary_state", style.OUT_PART7)
raise SystemExit(0 if style.all_ok() else 1)
