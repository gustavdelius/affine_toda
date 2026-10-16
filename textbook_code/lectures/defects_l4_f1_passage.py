"""Lecture notes on defects, Lecture 4, figure 1 (fig-l4-passage): three things a window does to a soliton.

The exact solutions Phi = 4 arctan Y on x < 0, Psi = 4 arctan(z Y) on x > 0, Y = exp(m0 (x cosh theta - t sinh theta)),
z = -coth((theta - theta_d)/2), of the type-I_- defect with theta_d = 0.8 (m0 = 1), for
(a) theta = 0.4 < theta_d: the soliton passes, delayed;
(b) theta = 1.4 > theta_d: it comes out as an antisoliton, delayed less;
(c) theta = theta_d - 5e-4: it is held by the defect for a time ln|z|/(m0 sinh theta), then released.
Colour: the charge density d_x Phi, red for a kink and blue for an antikink. Dotted: the undelayed world line.

Checks: the plotted pairs satisfy the defect conditions at x = 0; the ridge of the outgoing field follows
x = t tanh theta - ln|z|/(m0 cosh theta); the sign of the outgoing charge density is the sign of z.
Writes lecture-notes/defects/figures/l4_passage.svg.
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

m0, thd = 1.0, 0.8
sig = math.exp(-thd)
cases = [(0.4, "(a) $\\theta=0.4<\\theta_d$: delayed"),
         (1.4, "(b) $\\theta=1.4>\\theta_d$: converted"),
         (thd - 5e-4, "(c) $\\theta\\approx\\theta_d$: held, then released")]


def zf(theta):
    return -1 / math.tanh((theta - thd) / 2)


def fields(theta, t, x):
    """Phi, Phi_t, Phi_x on the left (factor 1) or Psi, Psi_t, Psi_x on the right (factor z)."""
    out = {}
    for side, fac in (("L", 1.0), ("R", zf(theta))):
        Y = fac * np.exp(m0 * (x * math.cosh(theta) - t * math.sinh(theta)))
        g = 4 * Y / (1 + Y**2)
        out[side] = (4 * np.arctan(Y), -g * m0 * math.sinh(theta), g * m0 * math.cosh(theta))
    return out


# defect conditions at x = 0: Phi_x = Psi_t - B_Phi, Psi_x = Phi_t + B_Psi
res = 0.0
for theta, _ in cases:
    for t in np.linspace(-8, 18, 27):
        f = fields(theta, t, 0.0)
        ph, pt, px = f["L"]
        ps, qt, qx = f["R"]
        Bp = m0 * (sig * math.sin((ph + ps) / 2) + math.sin((ph - ps) / 2) / sig)
        Bq = m0 * (sig * math.sin((ph + ps) / 2) - math.sin((ph - ps) / 2) / sig)
        res = max(res, abs(px - (qt - Bp)), abs(qx - (pt + Bq)))
report(f"the plotted pairs satisfy the defect conditions at x = 0 (residual {res:.0e})", res < 1e-12)

ok_ridge, ok_sign = True, True
for theta, _ in cases:
    z = zf(theta)
    for t in (12.0, 16.0):
        xs = np.linspace(0, 14, 28001)
        dens = fields(theta, t, xs)["R"][2]
        xr = xs[np.argmax(np.abs(dens))]
        xp = t * math.tanh(theta) - math.log(abs(z)) / (m0 * math.cosh(theta))
        if xp > 0.5:
            ok_ridge &= abs(xr - xp) < 1e-3
            ok_sign &= np.sign(dens[np.argmax(np.abs(dens))]) == np.sign(z)
report("the ridge of the outgoing field follows x = t tanh theta - ln|z|/(m0 cosh theta)", ok_ridge)
report("the outgoing charge density has the sign of z: kink for theta < theta_d, antikink for theta > theta_d",
       ok_sign)

fig, axes = plt.subplots(1, 3, figsize=(7.0, 3.1), sharey=True)
xs = np.linspace(-7, 12, 571)
ts = np.linspace(-8, 18, 521)
TT, XX = np.meshgrid(ts, xs, indexing="ij")
for ax, (theta, title) in zip(axes, cases):
    f = fields(theta, TT, XX)
    dens = np.where(XX < 0, f["L"][2], f["R"][2]) / (2 * m0 * math.cosh(theta))
    ax.pcolormesh(xs, ts, dens, cmap="RdBu_r", vmin=-1, vmax=1, shading="auto", rasterized=True)
    ax.axvline(0, color=style.PURPLE, lw=1.4)
    tl = np.array([0.0, 18.0])
    ax.plot(tl * math.tanh(theta), tl, ls=":", color=style.GREY, lw=1.0)
    ax.set_xlim(xs[0], xs[-1])
    ax.set_ylim(ts[0], ts[-1])
    ax.set_xlabel(r"$m_0x$")
    ax.set_title(title, fontsize=8, color=style.GREY)
    z = zf(theta)
    ax.text(-6.6, 16.4, "in: soliton", color=style.RED, fontsize=7.5)
    ax.text(11.6, -7.2, "out: " + ("soliton" if z > 0 else "antisoliton"),
            color=style.RED if z > 0 else style.BLUE, fontsize=7.5, ha="right")
axes[0].set_ylabel(r"$m_0t$")
dt = math.log(abs(zf(cases[2][0]))) / (m0 * math.sinh(cases[2][0]))
axes[2].annotate("", xy=(-0.5, dt), xytext=(-0.5, 0.0),
                 arrowprops=dict(arrowstyle="<->", color=style.GREY, lw=0.8))
axes[2].text(-0.9, dt / 2, r"$\dfrac{\ln|z|}{m_0\sinh\theta}$", color=style.GREY, fontsize=9, va="center",
             ha="right")
fig.tight_layout(w_pad=0.6)
style.save(fig, "l4_passage", OUT)
raise SystemExit(0 if style.all_ok() else 1)
