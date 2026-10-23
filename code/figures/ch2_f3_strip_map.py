"""Figure F3 (chapter 2, sec-rapidities, exr-strip): the strip 0 <= Im theta <= pi and
the cut s-plane.

s(theta) = m_a^2 + m_b^2 + 2 m_a m_b cosh(theta) maps the open strip bijectively onto
the s-plane minus the two cuts (-inf, (m_a - m_b)^2] and [(m_a + m_b)^2, inf).
Lines of constant Re theta go to ellipses, lines of constant Im theta to hyperbolas,
all with foci at the two thresholds.

Writes part1-background/figures/ch2_f3_strip_map.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
ma, mb = 1.0, 0.6
A, B = ma**2 + mb**2, 2 * ma * mb
s_of = lambda th: A + B * np.cosh(th)
s_plus, s_minus = (ma + mb) ** 2, (ma - mb) ** 2
X = 2.3                                    # half-width of the strip shown

# ---- checks ---------------------------------------------------------------
rng = np.random.default_rng(1)
N = 5000
th = rng.uniform(-X, X, N) + 1j * rng.uniform(0.01, np.pi - 0.01, N)
t = np.arccosh((s_of(th) - A) / B + 0j)
cand = np.stack([t, -t])
inside = (cand.imag > 0) & (cand.imag < np.pi)
rec = np.where(inside[0], cand[0], cand[1])
report("exactly one preimage in the open strip for 5000 random points",
       bool((inside.sum(axis=0) == 1).all()))
report("the preimage is the original theta (injective on the strip)", np.allclose(rec, th))
x = np.linspace(-X, X, 400)
report("Im theta = 0 maps to [(m_a+m_b)^2, inf)", np.all(s_of(x) >= s_plus - 1e-12))
report("Im theta = pi maps to (-inf, (m_a-m_b)^2]", np.all(s_of(x + 1j * np.pi).real <= s_minus + 1e-12)
       and np.allclose(s_of(x + 1j * np.pi).imag, 0))
y = np.linspace(0, np.pi, 400)
sv = s_of(1j * y)
report("theta = i u (0<u<pi) maps onto the segment between the thresholds",
       np.allclose(sv.imag, 0) and sv.real.max() <= s_plus + 1e-12 and sv.real.min() >= s_minus - 1e-12)

# ---- figure ---------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.9), gridspec_kw={"width_ratios": [1, 1.25]})

# left: the strip
ax1.axhspan(0, np.pi, color=style.LIGHTGREY, lw=0)
for xv in (-1.5, -0.75, 0.75, 1.5):
    ax1.plot([xv, xv], [0, np.pi], color=style.GREY, lw=0.6, alpha=0.6)
for yv in np.pi * np.array([1, 2, 3, 4, 5]) / 6:
    ax1.plot([-X, X], [yv, yv], color=style.GREY, lw=0.6, alpha=0.6)
ax1.plot([-X, X], [0, 0], color=style.BLUE, lw=2.6)
ax1.plot([-X, X], [np.pi, np.pi], color=style.RED, lw=2.6)
ax1.plot([0, 0], [0, np.pi], color=style.GREEN, lw=2.6)
u0 = 1.9
ax1.plot([0], [u0], "*", color=style.GREEN, ms=10)
ax1.text(0.12, u0 + 0.08, r"$\theta=iu$", color=style.GREEN)
ax1.text(-X + 0.1, 0.13, r"$\mathrm{Im}\,\theta=0$: physical region", color=style.BLUE, fontsize=7.5)
ax1.text(-X + 0.1, np.pi - 0.3, r"$\mathrm{Im}\,\theta=\pi$: crossed channel", color=style.RED, fontsize=7.5)
ax1.text(0.1, 0.55, "bound\nstates", color=style.GREEN, fontsize=7.5)
ax1.set_xlim(-X - 0.2, X + 0.2)
ax1.set_ylim(-0.35, np.pi + 0.45)
ax1.set_xlabel(r"$\mathrm{Re}\,\theta$")
ax1.set_ylabel(r"$\mathrm{Im}\,\theta$")
ax1.set_yticks([0, np.pi / 2, np.pi])
ax1.set_yticklabels([r"$0$", r"$\pi/2$", r"$\pi$"])
ax1.set_title(r"(a) rapidity strip", fontsize=8.5, color=style.GREY)

# right: the s-plane
for xv in (-1.5, -0.75, 0.75, 1.5):
    w = s_of(xv + 1j * y)
    ax2.plot(w.real, w.imag, color=style.GREY, lw=0.6, alpha=0.6)
for yv in np.pi * np.array([1, 2, 3, 4, 5]) / 6:
    w = s_of(x + 1j * yv)
    ax2.plot(w.real, w.imag, color=style.GREY, lw=0.6, alpha=0.6)
smax = s_of(X)
ax2.plot([s_plus, smax], [0, 0], color=style.BLUE, lw=3.2)
ax2.plot([s_of(X + 1j * np.pi).real, s_minus], [0, 0], color=style.RED, lw=3.2)
ax2.plot([s_minus, s_plus], [0, 0], color=style.GREEN, lw=3.2)
ax2.plot([A + B * np.cos(u0)], [0], "*", color=style.GREEN, ms=10)
ax2.text(A + B * np.cos(u0), 0.5, r"$s=m_c^2$", color=style.GREEN, ha="center")
for val, lab in ((s_minus, r"$(m_a-m_b)^2$"), (s_plus, r"$(m_a+m_b)^2$")):
    ax2.plot([val], [0], "|", color=style.GREY, ms=14)
    ax2.text(val, -0.75, lab, ha="center", fontsize=8)
ax2.text(smax - 0.2, 0.28, "physical\ncut", color=style.BLUE, ha="right", fontsize=7.5)
ax2.text(s_of(X + 1j * np.pi).real + 0.2, 0.28, "crossed-channel\ncut", color=style.RED, fontsize=7.5)
ax2.set_xlim(s_of(X + 1j * np.pi).real - 0.3, smax + 0.3)
ax2.set_ylim(-2.7, 2.7)
ax2.set_aspect("equal")
ax2.set_xlabel(r"$\mathrm{Re}\,s$")
ax2.set_ylabel(r"$\mathrm{Im}\,s$")
ax2.set_title(r"(b) $s=m_a^2+m_b^2+2m_am_b\cosh\theta$: the cut plane", fontsize=8.5, color=style.GREY)
fig.tight_layout()
style.save(fig, "ch2_f3_strip_map", style.OUT_PART1)
raise SystemExit(0 if style.all_ok() else 1)
