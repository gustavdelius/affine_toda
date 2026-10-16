"""Lecture notes on defects, Lecture 4, figure 3 (fig-l4-one-function): one function, read on four lines.

ln|z(u)| for the soliton delay z(u) = -coth(u/2) over the complex plane of u = theta - theta_d. Solitons of real
rapidity live on Im u = 0, where |z| > 1 and z has the pole of absorption at u = 0; breathers have constituents at
Im u = +-alpha; linear waves live on Im u = pi/2, where |z| = 1 and z is the transmission factor T; the odd sector
is Im u = pi, where z = -tanh(u/2), |z| < 1, and the zero sits at u = i pi.

Checks: |z| = 1 on Im u = +-pi/2; the pole at 0 and the zero at i pi; z(u + i pi) = 1/z(u).
Writes lecture-notes/defects/figures/l4_one_function.svg.
"""
import cmath
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects" / "figures"


def z(u):
    return -1 / np.tanh(u / 2)


us = np.linspace(-4, 4, 81)
report("|z| = 1 on Im u = pi/2 and Im u = -pi/2",
       max(abs(abs(z(u + 1j * s * math.pi / 2)) - 1) for u in us for s in (1, -1)) < 1e-12)
report("z has a pole at u = 0 and a zero at u = i pi",
       abs(1 / z(1e-9 + 0j)) < 1e-8 and abs(z(1j * math.pi + 1e-9)) < 1e-8)
report("z(u + i pi) = 1/z(u)",
       max(abs(z(u + 0.37j + 1j * math.pi) - 1 / z(u + 0.37j)) for u in us) < 1e-12)

fig, ax = plt.subplots(figsize=(5.6, 3.6))
re = np.linspace(-4, 4, 641)
im = np.linspace(-0.55 * math.pi, 1.55 * math.pi, 561)
RE, IM = np.meshgrid(re, im)
with np.errstate(divide="ignore", invalid="ignore"):
    L = np.log(np.abs(z(RE + 1j * IM)))
L = np.clip(L, -2.5, 2.5)
cs = ax.pcolormesh(re, im / math.pi, L, cmap="PuOr_r", vmin=-2.5, vmax=2.5, alpha=0.9, shading="auto",
                   rasterized=True)
ax.contour(RE, IM / math.pi, L, levels=[-2, -1, 1, 2], colors=style.GREY, linewidths=0.3)
cb = fig.colorbar(cs, ax=ax, shrink=0.85, pad=0.02, ticks=[-2, -1, 0, 1, 2])
cb.ax.tick_params(labelsize=7, colors=style.GREY)
cb.set_label(r"$\ln|z(u)|$", color=style.GREY)
lines = [(0.0, style.RED, "solitons: delayed, $|z|>1$", 0.07),
         (0.5, style.GREEN, "linear waves: $z=\\mathbb{T}$, $|z|=1$", 0.07),
         (1.0, style.BLUE, "odd sector: advanced, $|z|<1$", 0.07)]
for y, col, lab, dy in lines:
    ax.axhline(y, color=col, lw=1.4)
    ax.text(-3.85, y + dy, lab, color=col, fontsize=7.5, va="bottom")
ax.axhline(0.22, color=style.ORANGE, lw=1.0, ls="--")
ax.axhline(-0.22, color=style.ORANGE, lw=1.0, ls="--")
ax.text(3.85, 0.27, "breathers: $\\pm\\alpha$", color=style.ORANGE, fontsize=7.5, ha="right", va="bottom")
ax.plot(0, 0, marker="x", ms=7, mew=1.6, color="black", ls="none")
ax.plot(0, 1, marker="o", ms=6, mfc="white", mec="black", mew=1.2, ls="none")
ax.text(0.85, -0.13, "pole: absorption", fontsize=7.5, color=style.GREY)
ax.text(0.85, 0.87, "zero", fontsize=7.5, color=style.GREY)
ax.set_xlim(re[0], re[-1])
ax.set_ylim(im[0] / math.pi, im[-1] / math.pi)
ax.set_yticks([-0.5, 0, 0.5, 1, 1.5], [r"$-\frac{\pi}{2}$", "0", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$"])
ax.set_xlabel(r"$\operatorname{Re}u=\operatorname{Re}(\theta-\theta_d)$")
ax.set_ylabel(r"$\operatorname{Im}u$")
fig.tight_layout()
style.save(fig, "l4_one_function", OUT)
raise SystemExit(0 if style.all_ok() else 1)
