"""Figure F20 (chapter 15, sec-tba, exr-tba-free, exr-tba-uv): the TBA plateau and c_eff.

Sinh-Gordon (S = f_{-B/2}) at B = 1, mass m = 1.  Kernel K = -i d ln S/d theta
= 2 b cosh(theta)/(sinh^2 theta + b^2), b = sin(pi B/2), with total integral 2 pi (N = 1).
(a) pseudo-energy eps(theta) from eq-tba at three values of L: a wide central region forms whose
    edges move out like ln(2/L). Since N = 1 there is no finite plateau value: eps decreases
    (x = e^-eps grows) towards the centre, as in exr-tba-uv.
(b) c_eff(L) = (3L/pi^2) int cosh(theta) ln(1+e^{-eps}) dtheta: free fermion (S = -1) goes to 1/2,
    sinh-Gordon to 1, but only logarithmically slowly.

Writes part4-real-coupling/figures/ch15_f20_tba.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
Bv = 1.0
b = np.sin(np.pi * Bv / 2)
N = 1201
T = 24.0
th = np.linspace(-T, T, N)
dth = th[1] - th[0]
K = lambda t: 2 * b * np.cosh(t) / (np.sinh(t) ** 2 + b**2)
Kmat = K(th[:, None] - th[None, :])


def solve(L, kernel=True):
    eps = L * np.cosh(th)
    if not kernel:
        return eps
    for _ in range(600):
        new = L * np.cosh(th) - Kmat @ np.log1p(np.exp(-eps)) * dth / (2 * np.pi)
        if np.max(np.abs(new - eps)) < 1e-11:
            eps = new
            break
        eps = 0.5 * eps + 0.5 * new
    return eps


def ceff(L, eps):
    return 3 * L / np.pi**2 * np.sum(np.cosh(th) * np.log1p(np.exp(-eps))) * dth


# ---- checks ---------------------------------------------------------------
report("kernel integral is 2 pi (N = 1)", abs(np.sum(K(th)) * dth - 2 * np.pi) < 1e-3)
Ls = np.logspace(-6, 1, 8)
cf = np.array([ceff(L, L * np.cosh(th)) for L in Ls])
report("free fermion: c_eff(L=1e-6) = 1/2 within 1e-3", abs(cf[0] - 0.5) < 1e-3)
report("free fermion: c_eff increases to 1/2 as L -> 0 (decreasing in L)", np.all(np.diff(cf) < 0))
cs = np.array([ceff(L, solve(L)) for L in Ls])
report("sinh-Gordon: c_eff(L) decreases with L and stays below 1", np.all(np.diff(cs) < 1e-9) and cs[0] < 1)
report("sinh-Gordon: c_eff(L=1e-6) > 0.75 and approaches 1 slowly (below the free value)", 0.75 < cs[0] < 1)
report("sinh-Gordon: c_eff(1e-6) > c_eff(1e-3) > c_eff(1) (monotone toward 1)",
       cs[0] > cs[np.argmin(abs(Ls - 1e-3))] > cs[np.argmin(abs(Ls - 1))])
report("sinh-Gordon: c_eff(L=10) < 0.01 (massive decoupling)", cs[-1] < 0.01)

# ---- figure ---------------------------------------------------------------
fig, (a, c) = plt.subplots(1, 2, figsize=(7.0, 2.9))
cols = [style.BLUE, style.GREEN, style.ORANGE]
for L, col in zip((1e-1, 1e-3, 1e-6), cols):
    e = solve(L)
    a.plot(th, e, color=col, label=rf"$L=10^{{{int(round(np.log10(L)))}}}$")
a.set_ylim(-6, 8)
a.set_xlim(-16, 16)
a.set_xlabel(r"$\theta$")
a.set_ylabel(r"$\varepsilon(\theta)$")
a.legend(fontsize=7.5, loc="upper center")
a.axhline(0, color=style.GREY, lw=0.5, ls=":")
a.set_title("(a) pseudo-energy", fontsize=9, color=style.GREY)
c.semilogx(Ls, cs, color=style.BLUE, label="sinh-Gordon, $B=1$")
c.semilogx(Ls, cf, color=style.ORANGE, label=r"free fermion, $S=-1$")
c.axhline(1, color=style.GREY, lw=0.6, ls=":")
c.axhline(0.5, color=style.GREY, lw=0.6, ls=":")
c.set_xlabel(r"$L$")
c.set_ylabel(r"$c_{\rm eff}(L)$")
c.set_ylim(0, 1.1)
c.legend(fontsize=7.5, loc="lower left")
c.set_title(r"(b) effective central charge", fontsize=9, color=style.GREY)
fig.tight_layout()
style.save(fig, "ch15_f20_tba", style.OUT_PART4)
raise SystemExit(0 if style.all_ok() else 1)
