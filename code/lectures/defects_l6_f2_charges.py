"""Lecture notes on defects, Lecture 6, figure 2 (fig-l6-charges): the defect holds every charge of a soliton.

The exact passage Phi = 4 arctan Y on x < 0, Psi = 4 arctan(z Y) on x > 0, Y = exp(m0 (x cosh theta - t sinh theta)),
z = -coth((theta - theta_d)/2), for theta_d = 0.8 and theta = theta_d - 0.01 (m0 = beta = 1): the soliton is held
by the defect for a time ln|z|/(m0 sinh theta) and then released. For the charges of spin -1, +1, -3, +3 of the
generating function (Lecture 6), the parts carried by the field on x < 0, by the defect and by the field on x > 0,
and their sum, in units of the total:
    spin -1:  int (1/2 (d_+Phi)^2 + V) dx,  defect B - P_d;
    spin +1:  int (1/2 (d_-Phi)^2 + V) dx,  defect B + P_d;
    spin -3:  (1/8 m0^3) int [(d_+^2 Phi)^2 - (d_+Phi)^4/4 + m0^2 cos Phi (d_+Phi)^2] dx,
              defect D_-3 = (sigma/2m0^2) cos(Sigma/2) d_+Phi d_+Psi + sigma^3 (cos(Sigma/2) - (2/3) cos^3(Sigma/2)) - sigma^3/3;
    spin +3:  the same with d_- and Delta, defect D_+3 = -(1/(2 m0^2 sigma)) cos(Delta/2) d_-Phi d_-Psi
              + sigma^-3 (cos(Delta/2) - (2/3) cos^3(Delta/2)) - sigma^-3/3,
with V = m0^2 (1 - cos Phi), B, P_d of Lecture 3, Sigma = Phi + Psi, Delta = Phi - Psi at x = 0. The defect contributions
vanish on even jumps.

Checks: the pair satisfies the defect conditions; each sum is constant to 1e-8; the totals are the charges of a
soliton of rapidity theta (8 m0 e^-+theta for spin -+1, computed by quadrature for spin -+3); half-way through the
hold each defect contribution is within 5 % of the total; the defect contributions at the odd jump (0, 2 pi) are the charges of a
soliton of rapidity theta_d: 8 m0 sigma, 8 m0/sigma, and the spin -+3 charges at theta_d, to 1e-8.
Writes lecture-notes/defects/figures/l6_charges.svg.
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
theta = thd - 0.01
sig = math.exp(-thd)
z = -1 / math.tanh((theta - thd) / 2)
ch, sh = math.cosh(theta), math.sinh(theta)


def jets(t, x, fac, th=None):
    """Phi and its light-cone derivatives d_+, d_+^2, d_-, d_-^2 for 4 arctan(fac Y)."""
    th = theta if th is None else th
    c_, s_ = math.cosh(th), math.sinh(th)
    Y = fac * math.exp(m0 * (x * c_ - t * s_))
    F = 4 * math.atan(Y)
    g1 = 4 * Y / (1 + Y * Y)                     # dF/d ln Y
    g2 = 4 * Y * (1 - Y * Y) / (1 + Y * Y)**2     # d^2F/d(ln Y)^2
    kp = m0 * (c_ - s_)                          # d_+ ln Y = m0 e^{-theta}
    km = -m0 * (c_ + s_)                         # d_- ln Y = -m0 e^{theta}
    return F, kp * g1, kp * kp * g2, km * g1, km * km * g2


def densities(t, x, fac, th=None):
    F, up, vp, um, vm = jets(t, x, fac, th)
    V = m0**2 * (1 - math.cos(F))
    return (0.5 * up**2 + V, 0.5 * um**2 + V,
            (vp**2 - up**4 / 4 + m0**2 * math.cos(F) * up**2) / (8 * m0**3),
            (vm**2 - um**4 / 4 + m0**2 * math.cos(F) * um**2) / (8 * m0**3))


def bulk(t, fac, side, k, th=None):
    a, b = (-40.0, 0.0) if side == 0 else (0.0, 40.0)
    return quad(lambda x: densities(t, x, fac, th)[k], a, b, limit=400, epsabs=1e-13)[0]


def defect_parts(P, Q, upL, upR, umL, umR):
    S, D = P + Q, P - Q
    B = -2 * m0 * (sig * (math.cos(S / 2) - 1) + (math.cos(D / 2) - 1) / sig)
    Pd = 2 * m0 * (sig * (math.cos(S / 2) - 1) - (math.cos(D / 2) - 1) / sig)
    cS, cD = math.cos(S / 2), math.cos(D / 2)
    Dm3 = sig / (2 * m0**2) * cS * upL * upR + sig**3 * (cS - 2 / 3 * cS**3) - sig**3 / 3
    Dp3 = -1 / (2 * m0**2 * sig) * cD * umL * umR + sig**-3 * (cD - 2 / 3 * cD**3) - sig**-3 / 3
    return B - Pd, B + Pd, Dm3, Dp3


def defect_parts_at(t):
    PL = jets(t, 0.0, 1.0)
    PR = jets(t, 0.0, z)
    return defect_parts(PL[0], PR[0], PL[1], PR[1], PL[3], PR[3])


# defect conditions
res = 0.0
for t in np.linspace(-6, 16, 23):
    F, up, _, um, _ = jets(t, 0.0, 1.0)
    G, wp, _, wm, _ = jets(t, 0.0, z)
    res = max(res, abs((up - wp) + 2 * sig * m0 * math.sin((F + G) / 2)),
              abs((um + wm) - 2 * m0 / sig * math.sin((F - G) / 2)))
report(f"the plotted pair satisfies the defect conditions (residual {res:.0e})", res < 1e-12)

ts = np.linspace(-6, 16, 177)
parts = np.zeros((4, 3, len(ts)))
for i, t in enumerate(ts):
    w = defect_parts_at(t)
    for k in range(4):
        parts[k, 0, i] = bulk(t, 1.0, 0, k)
        parts[k, 1, i] = w[k]
        parts[k, 2, i] = bulk(t, z, 1, k)
tot = parts.sum(axis=1)
report("each sum is constant along the passage (to 1e-8)", max(np.ptp(tot[k]) for k in range(4)) < 1e-8)
# the charges of a single soliton of rapidity th, by quadrature over the whole line
def soliton_charge(k, th):
    return quad(lambda x: densities(0.0, x, 1.0, th)[k], -40, 40, limit=400, epsabs=1e-13)[0]


Qsol = [soliton_charge(k, theta) for k in range(4)]
report("the totals are the charges of the incoming soliton (8 m0 e^-theta, 8 m0 e^theta, and the spin -+3 charges)",
       abs(tot[0, 0] - 8 * m0 * math.exp(-theta)) < 1e-8 and abs(tot[1, 0] - 8 * m0 * math.exp(theta)) < 1e-8
       and abs(tot[2, 0] - Qsol[2]) < 1e-8 and abs(tot[3, 0] - Qsol[3]) < 1e-8)
t_a = brentq(lambda t: jets(t, 0.0, 1.0)[0] - math.pi, -20, 20)
t_b = brentq(lambda t: jets(t, 0.0, z)[0] - math.pi, -20, 40)
wmid = defect_parts_at(0.5 * (t_a + t_b))
report("half-way through the hold each defect contribution is within 5 % of the total",
       all(abs(wmid[k] / tot[k, 0] - 1) < 0.05 for k in range(4)))
wodd = defect_parts(0.0, 2 * math.pi, 0.0, 0.0, 0.0, 0.0)
Qd = [soliton_charge(k, thd) for k in range(4)]
report("the odd jump (0, 2 pi) stores the charges of a soliton of rapidity theta_d, for spins -1, +1, -3, +3 "
       f"(spin -3: {Qd[2]:.6f} = -(2/3) sigma^3 = {-2 / 3 * sig**3:.6f})",
       abs(wodd[0] - 8 * m0 * sig) < 1e-12 and abs(wodd[1] - 8 * m0 / sig) < 1e-12
       and all(abs(wodd[k] - Qd[k]) < 1e-8 for k in range(4)))

titles = ["(a) spin $-1$: $E-P$", "(b) spin $+1$: $E+P$", "(c) spin $-3$", "(d) spin $+3$"]
fig, axes = plt.subplots(2, 2, figsize=(7.0, 4.5), sharex=True)
for k, ax in enumerate(axes.ravel()):
    ax.plot(ts, parts[k, 0] / tot[k], color=style.BLUE, label=r"field on $x<0$")
    ax.plot(ts, parts[k, 1] / tot[k], color=style.PURPLE, label="defect")
    ax.plot(ts, parts[k, 2] / tot[k], color=style.RED, label=r"field on $x>0$")
    ax.plot(ts, tot[k] / tot[k], color=style.GREY, ls="--", lw=1.1, label="sum")
    ax.axvspan(t_a, t_b, color=style.LIGHTGREY, lw=0)
    ax.set_title(titles[k], fontsize=8.5, color=style.GREY)
    ax.set_xlim(ts[0], ts[-1])
    if k < 2:
        ax.set_ylim(-0.08, 1.15)
    else:
        lo, hi = (parts[2:] / tot[2:, None, :]).min(), (parts[2:] / tot[2:, None, :]).max()
        ax.set_ylim(lo - 0.15, hi + 0.15)
for ax in axes[1]:
    ax.set_xlabel(r"$m_0t$")
for ax in axes[:, 0]:
    ax.set_ylabel("share of the charge")

axes[0, 0].legend(loc="center right", fontsize=7, ncol=1)
fig.tight_layout(h_pad=0.8, w_pad=1.0)
style.save(fig, "l6_charges", OUT)
raise SystemExit(0 if style.all_ok() else 1)
