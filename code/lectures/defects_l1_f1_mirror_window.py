"""Lecture notes on defects, Lecture 1, figure 1 (fig-mirror-window): a mirror, a window, and a delta-function
defect that is both.

(a) A boundary: a particle of rapidity theta reflects with rapidity -theta.  (b) A purely transmitting defect:
the particle keeps its rapidity.  (c) The delta-function defect of a free boson, L = 1/2(d phi)^2 - 1/2 m^2 phi^2
- 1/2 g delta(x) phi^2: the incoming world line splits into a reflected and a transmitted one, drawn with widths
proportional to |R|^2 and |T|^2, R = g/(2ip-g), T = 2ip/(2ip-g), p = m sinh theta.

Checks: the world lines have slopes dx/dt = tanh(theta) and -tanh(theta); |R|^2 + |T|^2 = 1 at the plotted
rapidity and coupling; the plotted widths are proportional to |R|^2 and |T|^2.
Writes lecture-notes/defects/figures/l1_mirror_window.svg.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects" / "figures"

theta, m, g = 0.55, 1.0, 1.1
v = np.tanh(theta)
p = m * np.sinh(theta)
T = 2j * p / (2j * p - g)
R = g / (2j * p - g)
report("|R|^2 + |T|^2 = 1 at the plotted point", abs(abs(R) ** 2 + abs(T) ** 2 - 1) < 1e-14)
report("the defect both reflects and transmits here (both probabilities above 0.2)",
       abs(R) ** 2 > 0.2 and abs(T) ** 2 > 0.2)

fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.6))
t_hit = 0.0
ts_in = np.linspace(-1.6, t_hit, 50)
ts_out = np.linspace(t_hit, 1.6, 50)
lw0 = 2.4
slopes = []
for k, ax in enumerate(axes):
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.7, 1.7)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.annotate("", xy=(1.25, -1.6), xytext=(-1.25, -1.6),
                arrowprops=dict(arrowstyle="->", color=style.GREY, lw=0.7))
    ax.text(1.27, -1.6, "$x$", va="center", fontsize=8)
    ax.annotate("", xy=(-1.25, 1.65), xytext=(-1.25, -1.6),
                arrowprops=dict(arrowstyle="->", color=style.GREY, lw=0.7))
    ax.text(-1.22, 1.62, "$t$", va="top", fontsize=8)
    x_in = v * (ts_in - t_hit)
    ax.plot(x_in, ts_in, color=style.BLUE, lw=lw0)
    slopes.append(abs(np.polyfit(ts_in, x_in, 1)[0] - v) < 1e-12)
    if k == 0:
        ax.fill_betweenx([-1.55, 1.7], 0, 1.3, color=style.LIGHTGREY, lw=0)
        ax.plot([0, 0], [-1.55, 1.7], color=style.GREY, lw=2.0)
        x_out = -v * (ts_out - t_hit)
        ax.plot(x_out, ts_out, color=style.BLUE, lw=lw0)
        slopes.append(abs(np.polyfit(ts_out, x_out, 1)[0] + v) < 1e-12)
        ax.text(0.55, 0.9, "boundary", fontsize=8, ha="center")
        ax.text(-0.95, 1.15, r"$K(\theta)$", fontsize=8.5)
        ax.set_title("(a) a mirror", fontsize=8.5, color=style.GREY)
    elif k == 1:
        ax.plot([0, 0], [-1.55, 1.7], color=style.PURPLE, lw=2.0)
        x_out = v * (ts_out - t_hit)
        ax.plot(x_out, ts_out, color=style.BLUE, lw=lw0)
        slopes.append(abs(np.polyfit(ts_out, x_out, 1)[0] - v) < 1e-12)
        ax.text(0.08, -1.25, "defect", fontsize=8, color=style.PURPLE)
        ax.text(0.45, 1.15, r"$\mathbb{T}(\theta)$", fontsize=8.5)
        ax.set_title("(b) a window", fontsize=8.5, color=style.GREY)
    else:
        ax.plot([0, 0], [-1.55, 1.7], color=style.ORANGE, lw=2.0)
        wR, wT = lw0 * abs(R) ** 2, lw0 * abs(T) ** 2
        ax.plot(-v * (ts_out - t_hit), ts_out, color=style.BLUE, lw=wR)
        ax.plot(v * (ts_out - t_hit), ts_out, color=style.BLUE, lw=wT)
        report("plotted widths proportional to |R|^2 and |T|^2, summing to the incoming width",
               abs(wR + wT - lw0) < 1e-12)
        ax.text(-1.2, 1.15, rf"$|R|^2={abs(R)**2:.2f}$", fontsize=8)
        ax.text(0.35, 1.15, rf"$|\mathbb{{T}}|^2={abs(T)**2:.2f}$", fontsize=8)
        ax.text(0.08, -1.25, r"$g\,\delta(x)$", fontsize=8, color=style.ORANGE)
        ax.set_title("(c) both: a free field", fontsize=8.5, color=style.GREY)
report("world lines have slopes +-tanh(theta)", all(slopes))
fig.tight_layout(w_pad=1.0)
style.save(fig, "l1_mirror_window", OUT)
raise SystemExit(0 if style.all_ok() else 1)
