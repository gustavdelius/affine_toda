"""Lecture notes on defects, Lecture 2, figure 2 (fig-bianchi): Bianchi's square and the collision it describes.

(a) Bianchi's diagram: Phi_1 = B_{a1} Phi_0, Phi_2 = B_{a2} Phi_0, and Phi_12 = B_{a2} Phi_1 = B_{a1} Phi_2.
(b) The two-soliton solution Phi_12 obtained from the vacuum, a kink of rapidity theta1 = 0.9 and one of
theta2 = -0.4, drawn as the energy density in the (x, t) plane.  Dashed: the free world lines continued through the
collision.  Soliton 1 emerges displaced forward by -ln A/(m0 cosh theta1), A = tanh^2((theta1 - theta2)/2).

Checks: the plotted field solves sine-Gordon (finite differences); the ridge of soliton 1 at early and late times
sits on the two free world lines whose offset is -ln A/(m0 cosh theta1).
Writes lecture-notes/defects/figures/l2_bianchi.svg.
"""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects" / "figures"

th1, th2 = 0.9, -0.4
a1, a2 = -math.exp(-th1), -math.exp(-th2)
xi1, xi2 = 0.3, -0.2
c = (a1 + a2) / (a1 - a2)
A = ((a1 - a2) / (a1 + a2)) ** 2


def Yk(ak, xik, T, X):
    return np.exp(-0.5 * (ak + 1 / ak) * X + 0.5 * (1 / ak - ak) * T + xik)


def phi12(T, X):
    Y1, Y2 = Yk(a1, xi1, T, X), Yk(a2, xi2, T, X)
    return 4 * np.arctan2(c * (Y1 - Y2), 1 + Y1 * Y2)


# sine-Gordon residual by finite differences (m0 = 1), away from branch jumps
h = 1e-3
res = []
rng = np.random.default_rng(4)
for _ in range(20):
    T, X = rng.uniform(-6, 6), rng.uniform(-8, 8)
    f = lambda a_, b_: np.unwrap(np.array([phi12(a_, b_ - 2 * h), phi12(a_, b_ - h), phi12(a_, b_),
                                           phi12(a_, b_ + h), phi12(a_, b_ + 2 * h)]), period=2 * np.pi)
    fx = f(T, X)
    ft = np.unwrap(np.array([phi12(T - h, X), phi12(T, X), phi12(T + h, X)]), period=2 * np.pi)
    fxx = (fx[3] - 2 * fx[2] + fx[1]) / h ** 2
    ftt = (ft[2] - 2 * ft[1] + ft[0]) / h ** 2
    res.append(abs(ftt - fxx + np.sin(fx[2])))
report("plotted two-soliton field solves sine-Gordon (finite differences)", max(res) < 1e-4)

xs = np.linspace(-12, 12, 1201)
ts = np.linspace(-10, 10, 801)
TT, XX = np.meshgrid(ts, xs, indexing="ij")
P = np.unwrap(phi12(TT, XX), axis=1, period=2 * np.pi)
dx = xs[1] - xs[0]
dt = ts[1] - ts[0]
Px = np.gradient(P, dx, axis=1)
Pt = np.gradient(P, dt, axis=0)
dens = 0.5 * Px ** 2 + 0.5 * Pt ** 2 + (1 - np.cos(P))

v1, v2 = math.tanh(th1), math.tanh(th2)
# before the collision Y2 << 1 near soliton 1, and tan(Phi/4) ~ c Y1: its centre is where |c Y1| = 1;
# after it Y2 >> 1, tan(Phi/4) ~ -c/Y1: centre where |Y1| = |c|.  The shift of ln Y1 is 2 ln|c| = -ln A.
x1_in = lambda T: (-math.log(abs(c)) - xi1 - 0.5 * (1 / a1 - a1) * T) / (-0.5 * (a1 + 1 / a1))
x2_line = lambda T: (-xi2 - 0.5 * (1 / a2 - a2) * T) / (-0.5 * (a2 + 1 / a2))
shift = -math.log(A) / math.cosh(th1)
ok = True
for T, expect in ((-9.0, x1_in(-9.0)), (9.0, x1_in(9.0) + shift)):
    i = np.argmin(abs(ts - T))
    mask = np.abs(xs - x2_line(T)) > 4
    xpk = xs[mask][np.argmax(dens[i][mask])]
    ok &= abs(xpk - expect) < 0.05
report(f"ridge of soliton 1 follows the free world lines, displaced forward by -ln A/cosh(theta1) = {shift:.3f}", ok)

fig, (ax0, ax) = plt.subplots(1, 2, figsize=(7.0, 3.0), gridspec_kw=dict(width_ratios=[0.8, 1.3]))
# Bianchi square
pos = {"P0": (0, 0), "P1": (-1, 1), "P2": (1, 1), "P12": (0, 2)}
lab = {"P0": r"$\Phi_0$", "P1": r"$\Phi_1$", "P2": r"$\Phi_2$", "P12": r"$\Phi_{12}$"}
for k, (px, py) in pos.items():
    ax0.text(px, py, lab[k], ha="center", va="center", fontsize=10,
             bbox=dict(boxstyle="round,pad=0.25", fc="none", ec=style.GREY, lw=0.8))
for (p, q, txt, col) in (("P0", "P1", r"$a_1$", style.BLUE), ("P0", "P2", r"$a_2$", style.RED),
                         ("P1", "P12", r"$a_2$", style.RED), ("P2", "P12", r"$a_1$", style.BLUE)):
    (x0, y0), (x1, y1) = pos[p], pos[q]
    ax0.annotate("", xy=(x1 - 0.18 * (x1 - x0), y1 - 0.2 * (y1 - y0)), xytext=(x0 + 0.18 * (x1 - x0), y0 + 0.2 * (y1 - y0)),
                 arrowprops=dict(arrowstyle="->", color=col, lw=1.4))
    ax0.text((x0 + x1) / 2 + 0.18 * np.sign(x1 + x0 - 0.001 if x0 != x1 else 1), (y0 + y1) / 2, txt, color=col, fontsize=9,
             ha="center", va="center")
ax0.set_xlim(-1.6, 1.6)
ax0.set_ylim(-0.5, 2.5)
ax0.set_axis_off()
ax0.set_title("(a) Bianchi's square", fontsize=8.5, color=style.GREY)

im = ax.pcolormesh(xs, ts, dens, cmap="Blues", shading="auto", rasterized=True, vmin=0, vmax=np.percentile(dens, 99.5))
Tl = np.array([-10, 10])
ax.plot(x1_in(Tl), Tl, ls="--", color=style.ORANGE, lw=1.0)
ax.plot(x1_in(Tl) + shift, Tl, ls="--", color=style.ORANGE, lw=1.0)
ax.plot(x2_line(Tl), Tl, ls="--", color=style.GREY, lw=0.8)
ax.set_xlim(-12, 12)
ax.set_ylim(-10, 10)
ax.set_xlabel(r"$m_0x$")
ax.set_ylabel(r"$m_0t$")
ax.text(x1_in(4.5) + shift + 0.7, 4.5, r"soliton 1, $\theta_1$", color=style.ORANGE, fontsize=8, ha="left")
ax.text(x2_line(4.5) - 0.7, 4.5, r"soliton 2, $\theta_2$", color=style.GREY, fontsize=8, ha="right")
ax.annotate("", xy=(x1_in(-7.0) + shift, -7.0), xytext=(x1_in(-7.0), -7.0),
            arrowprops=dict(arrowstyle="->", color=style.ORANGE, lw=0.9))
ax.text(x1_in(-7.0) + shift + 0.6, -8.4, r"$\Delta x=-\ln A/(m_0\cosh\theta_1)$", color=style.ORANGE, fontsize=8)
ax.set_title(r"(b) $\Phi_{12}$ from the vacuum: energy density", fontsize=8.5, color=style.GREY)
fig.tight_layout(w_pad=1.0)
style.save(fig, "l2_bianchi", OUT)
raise SystemExit(0 if style.all_ok() else 1)
