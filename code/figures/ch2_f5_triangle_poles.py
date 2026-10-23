"""Figure F5 (chapter 2, sec-bootstrap-intro, sec-cdd; chapter 4, sec-lee-yang): bound
states and the sign of residues.

(a) The mass triangle: m_a, m_b, m_c with interior angle ubar = pi - u opposite m_c, so that
    m_c^2 = m_a^2 + m_b^2 + 2 m_a m_b cos u  (eq-bound-state-mass); the Lee-Yang case is
    equilateral, u = 2 pi/3.
(b),(c) Poles of the CDD factor f_x(theta) = (sinh t + i sin pi x)/(sinh t - i sin pi x) in the
    physical strip. The residue at i pi x is 2 i tan(pi x): i times a positive number for
    0<x<1/2 (an s-channel bound state in a unitary theory), i times a negative number for
    1/2<x<1 (scaling Lee-Yang, x = 2/3). Zeros sit at the negatives of the poles.

Writes part1-background/figures/ch2_f5_triangle_poles.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
f = lambda x, th: (np.sinh(th) + 1j * np.sin(np.pi * x)) / (np.sinh(th) - 1j * np.sin(np.pi * x))


def residue(x, th0, eps=1e-6):
    return eps * f(x, th0 + eps)


# ---- checks ---------------------------------------------------------------
for x in (0.25, 1 / 3, 2 / 3):
    r = residue(x, 1j * np.pi * x)
    report(f"x={x:.3f}: residue of f_x at i pi x equals 2i tan(pi x)",
           abs(r - 2j * np.tan(np.pi * x)) < 1e-4)
report("f_x depends on sin(pi x) only: f_{2/3} = f_{1/3}",
       abs(f(2 / 3, 0.37 + 0.4j) - f(1 / 3, 0.37 + 0.4j)) < 1e-12)
rng = np.random.default_rng(2)
th = rng.normal(size=20) + 1j * rng.uniform(0.1, 3.0, 20)
report("bootstrap f(t) = f(t - i pi/3) f(t + i pi/3) for x = 2/3",
       np.allclose(f(2 / 3, th), f(2 / 3, th - 1j * np.pi / 3) * f(2 / 3, th + 1j * np.pi / 3)))
ma, mb, u = 1.0, 0.75, 2.0
ubar = np.pi - u
mc = np.sqrt(ma**2 + mb**2 - 2 * ma * mb * np.cos(ubar))
report("law of cosines: m_c^2 = m_a^2 + m_b^2 + 2 m_a m_b cos u",
       abs(mc**2 - (ma**2 + mb**2 + 2 * ma * mb * np.cos(u))) < 1e-12)
report("Lee-Yang: u = 2 pi/3 and m_a = m_b gives m_c = m_a",
       abs(np.sqrt(2 * (1 + np.cos(2 * np.pi / 3))) - 1) < 1e-12)

# ---- figure ---------------------------------------------------------------
fig = plt.figure(figsize=(7.0, 3.1))
gs = fig.add_gridspec(1, 3, width_ratios=[1.6, 1, 1])
ax0 = fig.add_subplot(gs[0])
ax0.axis("off")
O = np.array([0.0, 0.0])
P = np.array([ma, 0.0])
Q = mb * np.array([np.cos(ubar), np.sin(ubar)])
ax0.fill([O[0], P[0], Q[0]], [O[1], P[1], Q[1]], color=style.LIGHTGREY)
ax0.plot([O[0], P[0], Q[0], O[0]], [O[1], P[1], Q[1], O[1]], color=style.GREY, lw=1.6)
ax0.text(0.5 * ma, -0.04, r"$m_a$", ha="center", va="top")
ax0.text(Q[0] * 0.5 - 0.06, Q[1] * 0.5 + 0.03, r"$m_b$", ha="right")
ax0.text(0.5 * (P[0] + Q[0]) + 0.05, 0.5 * (P[1] + Q[1]) + 0.03, r"$m_c$", color=style.GREEN)
arc = np.linspace(0, ubar, 40)
ax0.plot(0.17 * np.cos(arc), 0.17 * np.sin(arc), color=style.ORANGE)
ax0.text(0.2, 0.1, r"$\bar u$", color=style.ORANGE, fontsize=8.5)
ax0.text(0.45, 0.78, r"$m_c^2=m_a^2+m_b^2+2m_am_b\cos u$", ha="center", fontsize=7.5)
E = np.array([[1.45, 0.0], [2.05, 0.0], [1.75, 0.6 * np.sqrt(3) / 2]])
ax0.fill(E[:, 0], E[:, 1], color=style.LIGHTGREY)
ax0.plot(list(E[:, 0]) + [E[0, 0]], list(E[:, 1]) + [E[0, 1]], color=style.GREY, lw=1.4)
ax0.text(1.75, -0.04, "scaling Lee–Yang\n" r"$u=2\pi/3$", ha="center", va="top", fontsize=7.5)
ax0.set_xlim(-0.3, 2.2)
ax0.set_ylim(-0.3, 0.95)
ax0.set_aspect("equal")
ax0.set_title("(a) mass triangle", fontsize=9, color=style.GREY)


def strip_panel(ax, x, title, notes):
    ax.axhspan(0, np.pi, color=style.LIGHTGREY, lw=0)
    ax.axhline(0, color=style.BLUE, lw=2.2)
    ax.text(-0.9, 0.12, "physical", color=style.BLUE, fontsize=7)
    for yp in (np.pi * x, np.pi * (1 - x)):
        pos = residue(x, 1j * yp).imag > 0
        ax.plot([0], [yp], "X", color=style.GREEN if pos else style.RED, ms=9, zorder=5)
        ax.plot([0], [-yp], "o", mfc="none", mec=style.GREY, ms=6)
    for yp, txt, col in notes:
        ax.text(0.2, yp, txt, color=col, fontsize=7.5, va="center")
    ax.set_xlim(-1.0, 2.0)
    ax.set_ylim(-3.45, 3.45)
    ax.set_yticks([-np.pi, 0, np.pi])
    ax.set_yticklabels([r"$-\pi$", "0", r"$\pi$"])
    ax.set_xticks([])
    ax.set_ylabel(r"$\mathrm{Im}\,\theta$")
    ax.set_title(title, fontsize=9, color=style.GREY)


ax1 = fig.add_subplot(gs[1])
x1 = 0.25
strip_panel(ax1, x1, r"(b) unitary: $f_{1/4}$",
            [(np.pi * x1, r"$i\pi/4$" "\n" r"$B_1B_1\to B_2$", style.GREEN),
             (np.pi * (1 - x1), r"$3i\pi/4$" "\n" "crossed image", style.RED)])
ax2 = fig.add_subplot(gs[2])
x2 = 2 / 3
strip_panel(ax2, x2, r"(c) Lee–Yang: $f_{2/3}$",
            [(np.pi * x2, r"$2\pi i/3$" "\n" r"$1+1\to1$", style.RED),
             (np.pi * (1 - x2), r"$i\pi/3$" "\n" "crossed image", style.GREEN)])
fig.text(0.68, 0.0, r"green / red: sign of $\mathrm{Res}/i$;  circles: zeros", ha="center", fontsize=7.5,
         color=style.GREY)
fig.tight_layout(rect=(0, 0.04, 1, 1))
style.save(fig, "ch2_f5_triangle_poles", style.OUT_PART1)
raise SystemExit(0 if style.all_ok() else 1)
