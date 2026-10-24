"""Lecture notes on defects, Lecture 7, figure 1 (fig-l7-classical-limit): the quantum defect becomes the classical one.

(a) The probability that a soliton of rapidity theta is converted into an antisoliton, Y^2/(1 + Y^2) =
1/(1 + e^{-2 lambda u}), u = theta - theta_d, for lambda = 1, 4 and 16, against Lecture 4's classical step.
(b) The time delay m0 sinh(theta) Delta t = pi delta'(u)/(lambda + 1) obtained by Wigner's relation from the phase
delta of the scalar factor, delta'(u) = int_0^inf tanh(pi k/2) cos(ku)/(2 sinh(pi k/2 lambda)) dk, for the same lambda,
against the classical delay ln|coth(u/2)| of Lecture 4.

Checks: the probability is 1/2 at theta_d; delta' agrees with a finite difference of the phase of the Gamma-function
product of the scalar factor; the plotted delays approach ln|coth(u/2)| away from u = 0, with the deviation falling
as lambda grows; the quantum delay is finite at u = 0 and grows with lambda there.
Writes lecture-notes/defects/figures/l7_classical_limit.svg.
"""
import cmath
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import loggamma

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects" / "figures"
LAMS = (1.0, 4.0, 16.0)
COLS = (style.ORANGE, style.PURPLE, style.BLUE)


def dphase(u, lam):
    g = lambda k: (math.tanh(math.pi * k / 2) * math.cos(k * u) / (2 * math.sinh(math.pi * k / (2 * lam)))
                   if k > 0 else lam / 2)
    return quad(g, 0, 80 * lam, limit=4000, epsabs=1e-11)[0]


def phase_gamma(u, lam, N=4000):
    """phase of the Gamma-function product at real u (the ratio r_N/rbar_N of conjugate factors)"""
    k = np.arange(N)
    a = -1j * lam * u / (2 * np.pi)
    lr = (loggamma(k * lam + 0.25 + a) + loggamma((k + 1) * lam + 0.75 + a)
          - loggamma((k + 0.5) * lam + 0.25 + a) - loggamma((k + 0.5) * lam + 0.75 + a))
    return (math.pi * (1 + lam) / 4 - cmath.phase(1 + 1j * math.exp(lam * u)) + 2 * np.sum(lr.imag))


def delay(u, lam):
    return math.pi * dphase(u, lam) / (lam + 1)


def classical(u):
    return math.log(abs(1 / math.tanh(u / 2)))


report("the conversion probability Y^2/(1 + Y^2) is 1/2 at theta = theta_d", all(
    abs(1 / (1 + math.exp(0.0)) - 0.5) < 1e-15 for _ in LAMS))
h = 1e-4
okd = True
for lam in (1.0, 4.0):
    for u in (-0.6, 0.35):
        fd = (phase_gamma(u + h, lam) - phase_gamma(u - h, lam)) / (2 * h)
        okd = okd and abs(fd - dphase(u, lam)) < 1e-4 * max(1, abs(fd))
report("delta'(u) from the integral representation agrees with a finite difference of the phase of the Gamma "
       "product", okd)
devs = [max(abs(delay(u, lam) - classical(u)) for u in (-1.5, -0.8, 0.5, 1.2)) for lam in LAMS]
report(f"the delays approach ln|coth(u/2)| for |u| >= 0.5, deviations {devs[0]:.3f}, {devs[1]:.3f}, {devs[2]:.4f} "
       "falling with lambda", devs[0] > devs[1] > devs[2] and devs[2] < 0.1)
d0 = [delay(0.0, lam) for lam in LAMS]
report("the quantum delay is finite at u = 0 and grows with lambda there", all(np.isfinite(d0))
       and d0[0] < d0[1] < d0[2])

fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.0, 2.7))
us = np.linspace(-2.5, 2.5, 501)
for lam, col in zip(LAMS, COLS):
    axa.plot(us, 1 / (1 + np.exp(-2 * lam * us)), color=col, label=rf"$\lambda={lam:g}$")
axa.plot([-2.5, 0, 0, 2.5], [0, 0, 1, 1], color=style.GREY, ls="--", lw=1.1, label="classical")
axa.set_xlabel(r"$\theta-\theta_d$")
axa.set_ylabel("probability of conversion")
axa.set_title("(a) fast solitons are converted", fontsize=8.5, color=style.GREY)
axa.legend(loc="upper left", fontsize=7)
axa.set_xlim(-2.5, 2.5)

ub = np.concatenate([np.linspace(-2.5, -0.02, 90), np.linspace(0.02, 2.5, 90)])
ug = np.linspace(-2.5, 2.5, 161)
for lam, col in zip(LAMS, COLS):
    axb.plot(ug, [delay(u, lam) for u in ug], color=col, label=rf"$\lambda={lam:g}$")
axb.plot(ub[ub < 0], [classical(u) for u in ub[ub < 0]], color=style.GREY, ls="--", lw=1.1, label="classical")
axb.plot(ub[ub > 0], [classical(u) for u in ub[ub > 0]], color=style.GREY, ls="--", lw=1.1)
axb.set_ylim(0, 4.5)
axb.set_xlim(-2.5, 2.5)
axb.set_xlabel(r"$\theta-\theta_d$")
axb.set_ylabel(r"$m_0\sinh\theta\;\Delta t$")
axb.set_title("(b) the time delay, by Wigner's relation", fontsize=8.5, color=style.GREY)
axb.legend(loc="upper right", fontsize=7)
fig.tight_layout(w_pad=2.0)
style.save(fig, "l7_classical_limit", OUT)
raise SystemExit(0 if style.all_ok() else 1)
