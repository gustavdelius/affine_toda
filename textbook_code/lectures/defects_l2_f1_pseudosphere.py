"""Lecture notes on defects, Lecture 2, figure 1 (fig-pseudosphere): the static kink is the pseudosphere.

(a) The kink Phi = 4 arctan e^{m0 x} and its energy density (4 m0^2/beta^2) sech^2(m0 x).  (b) The pseudosphere
(sech s cos t, sech s sin t, s - tanh s), s = u + v = m0 x, t = u - v = m0 t_phys, with asymptotic lines u = const
and v = const drawn on it: these are the light rays x +- t = const of spacetime.  The rim s = 0 is the centre of
the kink, where the energy density is largest; far from the kink the surface shrinks to its axis.

Checks: the plotted surface has first fundamental form du^2 + 2 cos(Phi) du dv + dv^2 with Phi the kink (at sample
points, by finite differences); the energy density equals (4 m0^2/beta^2) radius^2; the drawn asymptotic lines are
curves u = const or v = const on the surface and have unit speed in the other coordinate.
Writes lecture-notes/defects/figures/l2_pseudosphere.svg.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects" / "figures"


def surf(U, V):
    s, t = U + V, U - V
    return np.array([np.cos(t) / np.cosh(s), np.sin(t) / np.cosh(s), s - np.tanh(s)])


def kink(s):
    return 4 * np.arctan(np.exp(s))


# checks
h = 1e-5
ok = True
rng = np.random.default_rng(3)
for _ in range(6):
    U, V = rng.normal(), rng.normal()
    ru = (surf(U + h, V) - surf(U - h, V)) / (2 * h)
    rv = (surf(U, V + h) - surf(U, V - h)) / (2 * h)
    ok &= abs(ru @ ru - 1) < 1e-8 and abs(rv @ rv - 1) < 1e-8 and abs(ru @ rv - np.cos(kink(U + V))) < 1e-8
report("plotted surface has I = du^2 + 2 cos(Phi) du dv + dv^2 with the kink Phi = 4 arctan e^{u+v}", ok)
xs = np.linspace(-4, 4, 801)
Px = np.gradient(kink(xs), xs)
dens = 0.5 * Px ** 2 + (1 - np.cos(kink(xs)))          # m0 = beta = 1
report("energy density = 4 sech^2(x) = 4 (radius of the pseudosphere)^2", np.max(np.abs(dens - 4 / np.cosh(xs) ** 2)) < 1e-3)

fig = plt.figure(figsize=(7.0, 3.3))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.15])
ax = fig.add_subplot(gs[0])
ax.plot(xs, kink(xs), color=style.BLUE, lw=1.8, label=r"$\Phi=4\arctan e^{m_0x}$")
ax.fill_between(xs, 0, dens, color=style.ORANGE, alpha=0.35, lw=0, label=r"energy density $\times\,\beta^2/m_0^2$")
ax.axvline(0, color=style.GREY, lw=0.6, ls=":")
ax.set_yticks([0, np.pi, 2 * np.pi])
ax.set_yticklabels(["0", r"$\pi$", r"$2\pi$"])
ax.set_xlabel(r"$m_0x$")
ax.set_xlim(-4, 4)
ax.set_ylim(0, 2 * np.pi + 0.3)
ax.legend(fontsize=7.5, loc="upper left")
ax.set_title("(a) the kink", fontsize=8.5, color=style.GREY)

ax3 = fig.add_subplot(gs[1], projection="3d")
S, T = np.meshgrid(np.linspace(-2.2, 2.2, 90), np.linspace(0, 2 * np.pi, 90))
X, Y, Z = np.cos(T) / np.cosh(S), np.sin(T) / np.cosh(S), S - np.tanh(S)
srf = ax3.plot_surface(X, Y, Z, color=style.BLUE, alpha=0.25, linewidth=0, shade=False, rasterized=True)
lines_ok = True
for c in np.linspace(-np.pi / 2, np.pi / 2, 4, endpoint=False) + np.pi / 8:
    ss = np.linspace(-2.2, 2.2, 300)
    # u = c: s = u + v, t = u - v  =>  t = 2c - s ; v = c: t = s - 2c
    for sign, col in ((1, style.RED), (-1, style.GREEN)):
        tt = (2 * c - ss) if sign == 1 else (ss - 2 * c)
        P = np.array([np.cos(tt) / np.cosh(ss), np.sin(tt) / np.cosh(ss), ss - np.tanh(ss)])
        ax3.plot(P[0], P[1], P[2], color=col, lw=0.7)
        # check: along the line the other coordinate runs at unit speed (|dr/dv| = 1 or |dr/du| = 1)
        Uc = (ss + tt) / 2
        Vc = (ss - tt) / 2
        k = 150
        if sign == 1:
            d = (surf(Uc[k], Vc[k] + h) - surf(Uc[k], Vc[k] - h)) / (2 * h)
            lines_ok &= abs(Uc[k] - c) < 1e-12
        else:
            d = (surf(Uc[k] + h, Vc[k]) - surf(Uc[k] - h, Vc[k])) / (2 * h)
            lines_ok &= abs(Vc[k] - c) < 1e-12
        lines_ok &= abs(np.linalg.norm(d) - 1) < 1e-8
report("drawn red and green curves are u = const and v = const, with unit speed along them", lines_ok)
tt = np.linspace(0, 2 * np.pi, 200)
ax3.plot(np.cos(tt), np.sin(tt), 0 * tt, color=style.ORANGE, lw=2.0)
ax3.set_box_aspect((1, 1, 1.5), zoom=1.25)
ax3.set_xlim(-1, 1)
ax3.set_ylim(-1, 1)
ax3.set_zlim(-1.4, 1.4)
ax3.set_axis_off()
ax3.view_init(elev=12, azim=35)
ax3.set_title("(b) the pseudosphere", fontsize=8.5, color=style.GREY, y=0.98)
ax3.text2D(0.80, 0.62, "rim:\ncentre of\nthe kink", transform=ax3.transAxes, fontsize=7.5, color=style.ORANGE)
ax3.text2D(0.02, 0.06, "red: $u$ = const, green: $v$ = const\n(light rays $x\\pm t$ = const)",
           transform=ax3.transAxes, fontsize=7.2)
fig.tight_layout(w_pad=0.5)
style.save(fig, "l2_pseudosphere", OUT)
raise SystemExit(0 if style.all_ok() else 1)
