"""Lecture notes on defects, Lecture 4, figure 2 (fig-l4-capture): a soliton near the defect rapidity is held.

The exact solution Phi = 4 arctan Y (x < 0), Psi = 4 arctan(z Y) (x > 0) for theta_d = 0.8 and
theta = theta_d - 5e-4 (m0 = beta = 1, M = 8), as in panel (c) of fig-l4-passage.
(a) The defect values Phi(0, t) and Psi(0, t): Phi drops from 2 pi to 0 when the soliton arrives, Psi only a time
ln|z|/(m0 sinh theta) later; in between the defect holds a jump of 2 pi.
(b) The energy budget in units of M: the field on x < 0, the defect B/beta^2, the field on x > 0, and the sum.
While the jump is held the defect carries almost exactly M cosh theta_d.

Checks: the half-way times differ by ln|z|/(m0 sinh theta); E + B/beta^2 is constant and equal to M cosh theta;
the maximum of B/beta^2 is within 2 |theta - theta_d| of M cosh theta_d.
Writes lecture-notes/defects/figures/l4_capture.svg.
"""
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects" / "figures"

m0, thd = 1.0, 0.8
eps = 5e-4
theta = thd - eps
sig = math.exp(-thd)
M = 8 * m0
z = -1 / math.tanh((theta - thd) / 2)
ch, sh = math.cosh(theta), math.sinh(theta)


def F(t, x, fac):
    Y = fac * np.exp(m0 * (x * ch - t * sh))
    g = 4 * Y / (1 + Y**2)
    return 4 * np.arctan(Y), -g * m0 * sh, g * m0 * ch


def energy(t, fac, a, b):
    def e(x):
        f, ft, fx = F(t, x, fac)
        return 0.5 * ft**2 + 0.5 * fx**2 + m0**2 * (1 - math.cos(f))
    return quad(e, a, b, limit=400, epsabs=1e-12)[0]


def Bd(ph, ps):
    return -2 * m0 * (sig * (math.cos((ph + ps) / 2) - 1) + (math.cos((ph - ps) / 2) - 1) / sig)


Phi0 = lambda t: F(t, 0.0, 1.0)[0]
Psi0 = lambda t: F(t, 0.0, z)[0]
t1 = brentq(lambda t: Phi0(t) - math.pi, -20, 20, xtol=1e-14)
t2 = brentq(lambda t: Psi0(t) - math.pi, -20, 40, xtol=1e-14)
dt = math.log(abs(z)) / (m0 * sh)
report(f"the half-way times differ by ln|z|/(m0 sinh theta) = {dt:.3f}/m0", abs((t2 - t1) - dt) < 1e-10)

ts = np.linspace(-6, 16, 177)
EL = np.array([energy(t, 1.0, -50, 0) for t in ts])
ER = np.array([energy(t, z, 0, 50) for t in ts])
BB = np.array([Bd(Phi0(t), Psi0(t)) for t in ts])
tot = EL + ER + BB
report("E + B/beta^2 is constant and equal to M cosh theta along the plotted solution",
       np.ptp(tot) < 1e-8 and abs(tot[0] - M * ch) < 1e-8)
report("the maximum of B/beta^2 is within 2 |theta - theta_d| of M cosh theta_d",
       abs(BB.max() / (M * math.cosh(thd)) - 1) < 2 * eps)

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.9), gridspec_kw=dict(width_ratios=[1.0, 1.15]))
tf = np.linspace(-6, 16, 1101)
ax.plot(tf, [Phi0(t) / math.pi for t in tf], color=style.BLUE, label=r"$\Phi(0,t)$")
ax.plot(tf, [Psi0(t) / math.pi for t in tf], color=style.RED, label=r"$\Psi(0,t)$")
ax.axvspan(t1, t2, color=style.LIGHTGREY, lw=0)
ax.annotate("", xy=(t2, 1.0), xytext=(t1, 1.0), arrowprops=dict(arrowstyle="<->", color=style.GREY, lw=0.8))
ax.text(0.5 * (t1 + t2), 1.12, "jump of $2\\pi$ held", color=style.GREY, fontsize=7.5, ha="center")
ax.set_yticks([0, 1, 2], ["0", r"$\pi$", r"$2\pi$"])
ax.set_xlim(tf[0], tf[-1])
ax.set_ylim(-0.1, 2.25)
ax.set_xlabel(r"$m_0t$")
ax.legend(loc="lower left", fontsize=7.5)
ax.set_title("(a) the defect values", fontsize=8.5, color=style.GREY)

ax2.plot(ts, EL / M, color=style.BLUE, label=r"field on $x<0$")
ax2.plot(ts, BB / M, color=style.PURPLE, label=r"defect, $\mathcal{B}/\beta^2$")
ax2.plot(ts, ER / M, color=style.RED, label=r"field on $x>0$")
ax2.plot(ts, tot / M, color=style.GREY, ls="--", lw=1.1, label="sum")
ax2.set_xlim(ts[0], ts[-1])
ax2.set_ylim(-0.05, 1.75)
ax2.set_xlabel(r"$m_0t$")
ax2.set_ylabel(r"energy$/M$")
ax2.legend(loc="upper right", fontsize=7, ncol=2)
ax2.set_title("(b) the energy budget", fontsize=8.5, color=style.GREY)
fig.tight_layout(w_pad=1.2)
style.save(fig, "l4_capture", OUT)
raise SystemExit(0 if style.all_ok() else 1)
