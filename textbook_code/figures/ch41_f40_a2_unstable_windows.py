"""Figure F40 (chapter 41, sec-an-defect-bound-states, eq-dt-pole, eq-dt-window, prp-a2-unstable): the scalar-factor
poles of the a_2^(1) type-I_- defect and the windows in which they describe unstable states.

The species-a pole is at theta_a = theta_d + i pi/2 - i vartheta_a, vartheta_a = pi(h - a)/h + pi/(h omega)
(eq-dt-pole for a = 1, its conjugate for a = r), so Im(theta_a - theta_d)/pi = a/h - 1/2 - (1/omega)/h; for a_2 (h = 3)
species 1 is at -1/6 - 1/(3 omega), species 2 at 1/6 - 1/(3 omega).  With theta_d > 0 the energy
E_a = M_a cosh theta_a has Re E_a > 0 and Im E_a < 0 (an unstable state) exactly when -1/2 < Im(theta_a - theta_d)/pi <= 0,
the window eq-dt-window max(0, a - h/2) <= 1/omega < a: species 1 for omega > 1, species 2 for 1/2 < omega <= 2.
The coupling is omega = 4 pi/beta^2 - 1, so 1/omega = 0 is the classical limit and beta^2/pi = 4 x/(1 + x), x = 1/omega.

Checks: g (defect_common.an_scalar_g, r = 2) has a simple pole at theta_1 and gbar (an_scalar_gbar) at the species-2
position, for omega = 0.6, 1.37, 3; the sign test Re E > 0, Im E <= 0 on a grid of 1/omega reproduces eq-dt-window for
a = 1, 2 (and h = 4, 5, 6) and gives the a_2 windows omega > 1 and 1/2 < omega <= 2, i.e. beta^2 < 2 pi and
4 pi/3 <= beta^2 < 8 pi/3; as 1/omega -> 0 the species-1 pole tends to the classical absorption rapidity
theta_d - i(pi/2 - pi/h) (eq-ds-absorption).

Writes part8-defects/figures/ch41_f40_a2_unstable_windows.svg.
"""
import sys
from pathlib import Path

import numpy as np

import style
from style import plt, report

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from defect_common import an_scalar_g, an_scalar_gbar  # noqa: E402

style.setup()
H = 3


def y_pole(a, x, h=H):
    """Im(theta_a - theta_d)/pi as a function of x = 1/omega."""
    return a / h - 0.5 - x / h


def unstable(a, h, x, thd=0.4):
    E = np.cosh(thd + 1j * np.pi * y_pole(a, x, h))
    return E.real > 0 and E.imag <= 1e-14


# ---- checks ----------------------------------------------------------------
ok = True
out = []
for om in (0.6, 1.37, 3.0):
    D1 = 1j * np.pi * y_pole(1, 1 / om)                 # theta_1 - theta_d = i pi/2 - i(2 pi/3 + pi/(3 omega))
    D2 = 1j * np.pi * y_pole(2, 1 / om)
    e1, e2 = 1e-4, 1e-5
    p1 = abs(e1 * an_scalar_g(D1 + e1, 2, om) / (e2 * an_scalar_g(D1 + e2, 2, om)) - 1)
    p2 = abs(e1 * an_scalar_gbar(D2 + e1, 2, om) / (e2 * an_scalar_gbar(D2 + e2, 2, om)) - 1)
    ok &= p1 < 1e-3 and p2 < 1e-3
    out.append(f"{om}: {max(p1, p2):.0e}")
report("eq-dt-pole: g has a simple pole at theta_d - i pi/6 - i pi/(3 omega) (species 1) and gbar at "
       "theta_d + i pi/6 - i pi/(3 omega) (species 2), omega = " + ", ".join(out), ok)

xs = np.linspace(0, 4, 8001) + 1e-4 / np.e             # a grid that avoids the window ends themselves
ok = True
for h in (3, 4, 5, 6):
    for a in (1, h - 1):
        pred = np.array([max(0, a - h / 2) <= x < a for x in xs])
        ok &= np.array_equal(pred, np.array([unstable(a, h, x) for x in xs]))
s1 = xs[[unstable(1, 3, x) for x in xs]]
s2 = xs[[unstable(2, 3, x) for x in xs]]
report(f"Re E > 0, Im E <= 0 (theta_d = 0.4) reproduces eq-dt-window for a = 1, r and h = 3..6; a_2: species 1 for "
       f"0 <= 1/omega < {s1.max():.3f}, species 2 for {s2.min():.3f} <= 1/omega < {s2.max():.3f}",
       ok and s1.min() < 1e-3 and abs(s1.max() - 1) < 1e-3 and abs(s2.min() - 0.5) < 1e-3 and abs(s2.max() - 2) < 1e-3)
b2 = lambda x: 4 * x / (1 + x)                          # beta^2/pi from omega = 4 pi/beta^2 - 1
report("in beta^2: species 1 unstable for beta^2 < 2 pi, species 2 for 4 pi/3 <= beta^2 < 8 pi/3",
       abs(b2(1) - 2) < 1e-15 and abs(b2(0.5) - 4 / 3) < 1e-15 and abs(b2(2) - 8 / 3) < 1e-15)
report("classical limit 1/omega -> 0: species-1 pole -> theta_d - i(pi/2 - pi/h) (classical absorption, eq-ds-absorption), "
       "species-2 pole -> theta_d + i pi/6 (no classical counterpart in the window)",
       abs(y_pole(1, 0) + (0.5 - 1 / 3)) < 1e-15 and abs(y_pole(2, 0) - 1 / 6) < 1e-15)

# ---- figure ----------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.0, 2.9))
X1 = 3.0
x = np.linspace(0, X1, 601)
ax.axhspan(0, 0.5, color=style.GREY, alpha=0.12, lw=0)
ax.axhspan(-0.5, 0, color=style.GREEN, alpha=0.12, lw=0)
ax.axhline(0, color=style.GREY, lw=0.6)
ax.axhline(-0.5, color=style.GREY, lw=0.6)
ax.text(X1 - 0.05, 0.47, r"physical strip, ${\rm Im}\,E>0$", ha="right", va="top", fontsize=7.2)
ax.text(X1 - 0.05, -0.03, r"unstable: ${\rm Re}\,E>0$, ${\rm Im}\,E<0$", ha="right", va="top", fontsize=7.2,
        color=style.GREEN)
ax.text(X1 - 0.05, -0.53, r"${\rm Re}\,E<0$", ha="right", va="top", fontsize=7.2)
for a, c in ((1, style.BLUE), (2, style.RED)):
    y = y_pole(a, x)
    win = np.array([unstable(a, H, xx) for xx in x])
    ax.plot(x, y, color=c, lw=1.0, alpha=0.7)
    ax.plot(x, np.where(win, y, np.nan), color=c, lw=2.6, label=f"species {a}")
    ax.plot([0], [y_pole(a, 0)], "o", color=c, ms=4.5, zorder=5, clip_on=False)
for xv in (0.5, 1, 2):
    ax.axvline(xv, color=style.GREY, lw=0.7, ls=":")
ax.text(0.07, -0.035, "classical absorption", color=style.BLUE, va="top", fontsize=7.2)
ax.set_xlim(0, X1)
ax.set_ylim(-1.2, 0.5)
ax.set_xticks([0, 0.5, 1, 2, 3])
ax.set_xticklabels(["0", r"$\frac{1}{2}$", "1", "2", "3"])
ax.set_yticks([-1, -0.5, -1 / 6, 0, 1 / 6, 0.5])
ax.set_yticklabels(["$-1$", r"$-\frac{1}{2}$", r"$-\frac{1}{6}$", "0", r"$\frac{1}{6}$", r"$\frac{1}{2}$"])
ax.set_xlabel(r"$1/\omega$   (classical limit at 0)")
ax.set_ylabel(r"${\rm Im}(\theta_a-\theta_d)/\pi$")
top = ax.secondary_xaxis("top", functions=(lambda u: u, lambda u: u))
top.set_xticks([1 / 3, 0.5, 1, 2, 3])
top.set_xticklabels(["1", r"$\frac{4}{3}$", "2", r"$\frac{8}{3}$", "3"])
top.set_xlabel(r"$\beta^2/\pi$", fontsize=8)
ax.legend(loc="lower left", fontsize=7.5, handlelength=1.6)

fig.tight_layout()
style.save(fig, "ch41_f40_a2_unstable_windows", style.OUT_PART8)
raise SystemExit(0 if style.all_ok() else 1)
