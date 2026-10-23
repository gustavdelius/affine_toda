"""Lecture notes on defects, Lecture 1, figure 2 (fig-delta-poles): where the poles of the delta-function defect go.

T(theta) = 2ip/(2ip - g), p = m sinh theta, has its poles where sinh theta = -i g/(2m).  As g runs from 0 to
-2m a pole climbs the imaginary axis from 0 to i pi/2 (a bound state, in the physical strip 0 < Im theta < pi/2
of the defect, with energy m cos u); for g < -2m the pair leaves the axis along Im theta = pi/2 (imaginary
energy: the vacuum is unstable).  For 0 < g < 2m a pole descends from 0 to -i pi/2 (a virtual state), and for
g > 2m the pair moves along Im theta = -pi/2.  The partner poles at i pi - theta are drawn faintly.
Right panel: the energy m cos u of the bound state against -g/2m, compared with the bound state
sqrt(m^2 - g^2/4) of the equation -phi'' + g delta(x) phi = (E^2 - m^2) phi.

Checks: every plotted pole satisfies 2i m sinh theta = g; the bound-state energies agree with the
quantum-mechanical formula; the pole is in the physical strip exactly for -2m < g < 0.
Writes lecture-notes/defects/figures/l1_delta_poles.svg.
"""
import cmath
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects" / "figures"
m = 1.0


def poles(g):
    """The poles of T in -pi < Im theta <= pi: theta and i pi - theta with sinh theta = -i g/2m."""
    z = cmath.asinh(-1j * g / (2 * m))
    return [z, 1j * math.pi - z]


ok_pole, ok_strip = True, True
gs_neg = np.linspace(-3.2, 0, 400)
gs_pos = np.linspace(0, 3.2, 400)
traj = {"neg": [], "pos": []}
for key, gs in (("neg", gs_neg), ("pos", gs_pos)):
    for gg in gs:
        for z in poles(gg):
            ok_pole &= abs(2j * m * cmath.sinh(z) - gg) < 1e-9
        traj[key].append(poles(gg))
        z = poles(gg)[0]
        in_strip = 0 < z.imag < math.pi / 2 - 1e-12 and abs(z.real) < 1e-12
        ok_strip &= in_strip == (-2 * m < gg < 0)
report("plotted poles satisfy 2 i m sinh(theta) = g", ok_pole)
report("a pole lies in the physical strip 0 < Im theta < pi/2 exactly for -2m < g < 0", ok_strip)

ok_E = True
us, Es, Eqm, xs = [], [], [], []
for gg in np.linspace(-1.999, -0.001, 200):
    u = poles(gg)[0].imag
    E = m * math.cos(u)
    qm = math.sqrt(m ** 2 - gg ** 2 / 4)
    ok_E &= abs(E - qm) < 1e-12
    xs.append(-gg / (2 * m)); Es.append(E); Eqm.append(qm)
report("bound-state energy m cos u equals sqrt(m^2 - g^2/4)", ok_E)

fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.0, 2.9), gridspec_kw=dict(width_ratios=[1.45, 1]))
ax.axhspan(0, math.pi / 2, color=style.GREEN, alpha=0.12, lw=0)
ax.text(1.55, 0.78, "physical strip\nof the defect", fontsize=7.5, color=style.GREEN, ha="center")
for key, col in (("neg", style.BLUE), ("pos", style.RED)):
    arr = np.array(traj[key])
    # main pole, both branches off the axis are traced by +-Re
    for idx, alpha in ((0, 1.0), (1, 0.25)):
        zz = arr[:, idx]
        ax.plot(zz.real, zz.imag, color=col, lw=1.6, alpha=alpha)
        ax.plot(-zz.real, zz.imag, color=col, lw=1.6, alpha=alpha)
ax.axhline(0, color=style.GREY, lw=0.6)
ax.axhline(math.pi / 2, color=style.GREY, lw=0.4, ls=":")
ax.axhline(-math.pi / 2, color=style.GREY, lw=0.4, ls=":")
ax.set_xlim(-2.2, 2.2)
ax.set_ylim(-math.pi - 0.1, math.pi + 0.1)
ax.set_yticks([-math.pi, -math.pi / 2, 0, math.pi / 2, math.pi])
ax.set_yticklabels([r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
ax.set_xlabel(r"$\mathrm{Re}\,\theta$")
ax.set_ylabel(r"$\mathrm{Im}\,\theta$")
ax.text(0.07, 0.62, "bound\nstate", fontsize=7.5, color=style.BLUE)
ax.text(0.9, 1.72, r"$g<-2m$: unstable", fontsize=7.5, color=style.BLUE)
ax.text(0.07, -0.95, "virtual\nstate", fontsize=7.5, color=style.RED)
ax.text(0.9, -1.92, r"$g>2m$", fontsize=7.5, color=style.RED)
ax.set_title(r"(a) poles of $\mathbb{T}$ as $g$ varies", fontsize=8.5, color=style.GREY)

bx.plot(xs, Es, color=style.BLUE, lw=1.8, label=r"$m\cos u$ (pole)")
bx.plot(xs[::12], Eqm[::12], "o", ms=3, color=style.ORANGE, label=r"$\sqrt{m^2-g^2/4}$")
bx.set_xlabel(r"$-g/2m$")
bx.set_ylabel(r"bound-state energy$/m$")
bx.set_xlim(0, 1)
bx.set_ylim(0, 1.05)
bx.legend(fontsize=7.5, loc="lower left")
bx.set_title("(b) the bound state", fontsize=8.5, color=style.GREY)
fig.tight_layout(w_pad=1.5)
style.save(fig, "l1_delta_poles", OUT)
raise SystemExit(0 if style.all_ok() else 1)
