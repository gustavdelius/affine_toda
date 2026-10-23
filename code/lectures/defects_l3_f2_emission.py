"""Lecture notes on defects, Lecture 3, figure 2 (fig-l3-emission): a jump walks out as a soliton.

The exact solution Phi = 0 on x < 0, Psi = 4 arctan exp(m0 (x cosh theta_d - t sinh theta_d)) on x > 0 of the
frozen Baecklund defect with sigma = e^{-theta_d}, theta_d = 0.6 (m0 = beta = 1, M = 8).
(a) The energy density in the (x, t) plane; the dotted line continues the soliton's world line x = t tanh theta_d
into x < 0, where there is no field.
(b) The energy budget: the defect energy B(t), the bulk energy E(t) on x > 0, and their sum, in units of M.
The jump of 2 pi at early times stores M cosh theta_d, which leaves as the soliton.

Checks: the plotted pair satisfies the defect conditions at x = 0; E + B and P + P_d are constant along it;
B starts at M cosh theta_d and ends at 0.
Writes lecture-notes/defects/figures/l3_emission.svg.
"""
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects" / "figures"

m0, thd = 1.0, 0.6
sig = math.exp(-thd)
M = 8 * m0
ch, sh = math.cosh(thd), math.sinh(thd)


def Y(t, x):
    return np.exp(m0 * (x * ch - t * sh))


def Psi(t, x):
    return 4 * np.arctan(Y(t, x))


def Psi_t(t, x):
    y = Y(t, x)
    return 4 * y / (1 + y**2) * (-m0 * sh)


def Psi_x(t, x):
    y = Y(t, x)
    return 4 * y / (1 + y**2) * (m0 * ch)


def B_of(ph, ps):
    return -2 * m0 * (sig * (math.cos((ph + ps) / 2) - 1) + (math.cos((ph - ps) / 2) - 1) / sig)


def Pd_of(ph, ps):
    return 2 * m0 * (sig * (math.cos((ph + ps) / 2) - 1) - (math.cos((ph - ps) / 2) - 1) / sig)


def bulk(t):
    e = lambda x: 0.5 * Psi_t(t, x)**2 + 0.5 * Psi_x(t, x)**2 + m0**2 * (1 - math.cos(Psi(t, x)))
    p = lambda x: -Psi_t(t, x) * Psi_x(t, x)
    return quad(e, 0, 60, limit=400, epsabs=1e-12)[0], quad(p, 0, 60, limit=400, epsabs=1e-12)[0]


# defect conditions at x = 0 with Phi = 0: Phi_x = Psi_t - B_Phi, Psi_x = Phi_t + B_Psi
res = 0.0
for t in np.linspace(-10, 10, 21):
    ps = Psi(t, 0.0)
    Bphi = m0 * (sig * math.sin(ps / 2) + math.sin(-ps / 2) / sig)
    Bpsi = m0 * (sig * math.sin(ps / 2) - math.sin(-ps / 2) / sig)
    res = max(res, abs(0.0 - (Psi_t(t, 0.0) - Bphi)), abs(Psi_x(t, 0.0) - (0.0 + Bpsi)))
report(f"the plotted pair satisfies the defect conditions at x = 0 (residual {res:.0e})", res < 1e-12)

ts_b = np.linspace(-8, 14, 111)
Eb, Pb, Bd, Pdd = [], [], [], []
for t in ts_b:
    e_, p_ = bulk(t)
    Eb.append(e_)
    Pb.append(p_)
    Bd.append(B_of(0.0, Psi(t, 0.0)))
    Pdd.append(Pd_of(0.0, Psi(t, 0.0)))
Eb, Pb, Bd, Pdd = map(np.array, (Eb, Pb, Bd, Pdd))
report("E + B and P + P_d are constant along the plotted solution",
       np.ptp(Eb + Bd) < 1e-8 and np.ptp(Pb + Pdd) < 1e-8)
report("B starts at M cosh theta_d and ends at 0",
       abs(B_of(0.0, Psi(-40.0, 0.0)) - M * ch) < 1e-10 and abs(B_of(0.0, Psi(40.0, 0.0))) < 1e-10)

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.0), gridspec_kw=dict(width_ratios=[1.0, 1.1]))
xs = np.linspace(-6, 12, 721)
ts = np.linspace(-8, 14, 561)
TT, XX = np.meshgrid(ts, xs, indexing="ij")
dens = np.where(XX > 0, 0.5 * Psi_t(TT, XX)**2 + 0.5 * Psi_x(TT, XX)**2 + m0**2 * (1 - np.cos(Psi(TT, XX))), 0.0)
ax.pcolormesh(xs, ts, dens, cmap="Blues", shading="auto", rasterized=True, vmin=0, vmax=np.percentile(dens, 99.7))
ax.axvline(0, color=style.PURPLE, lw=1.4)
tl = np.array([-8, 14])
ax.plot(tl * math.tanh(thd), tl, ls=":", color=style.ORANGE, lw=1.0)
ax.text(0.25, -7.2, "defect", color=style.PURPLE, fontsize=8, ha="left")
ax.text(-5.6, 9.5, r"$\Phi=0$", color=style.GREY, fontsize=9)
ax.text(7.2, 9.0, "soliton,\n" + r"$\theta=\theta_d$", color=style.ORANGE, fontsize=8, ha="left")
ax.set_xlim(-6, 12)
ax.set_ylim(-8, 14)
ax.set_xlabel(r"$m_0x$")
ax.set_ylabel(r"$m_0t$")
ax.set_title("(a) energy density", fontsize=8.5, color=style.GREY)

ax2.plot(ts_b, Bd / M, color=style.PURPLE, label=r"defect, $\mathcal{B}/\beta^2$")
ax2.plot(ts_b, Eb / M, color=style.BLUE, label=r"field on $x>0$")
ax2.plot(ts_b, (Eb + Bd) / M, color=style.GREY, ls="--", lw=1.2, label="sum")
ax2.axhline(ch, color=style.LIGHTGREY, lw=0.6)
ax2.text(13.8, ch + 0.04, r"$M\cosh\theta_d$", color=style.GREY, fontsize=8, ha="right", va="bottom")
ax2.set_xlim(ts_b[0], ts_b[-1])
ax2.set_ylim(-0.05, 1.45)
ax2.set_xlabel(r"$m_0t$")
ax2.set_ylabel(r"energy$/M$")
ax2.legend(loc="center right", fontsize=7.5)
ax2.set_title("(b) the energy budget", fontsize=8.5, color=style.GREY)
fig.tight_layout(w_pad=1.2)
style.save(fig, "l3_emission", OUT)
raise SystemExit(0 if style.all_ok() else 1)
