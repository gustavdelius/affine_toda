"""Lecture notes on defects, Lecture 3, figure 1 (fig-l3-landscape): the defect potential and what a jump stores.

(a) The frozen Baecklund potential B = -2 m0 [sigma (cos((Phi+Psi)/2) - 1) + sigma^{-1} (cos((Phi-Psi)/2) - 1)]
over the plane of defect values (Phi, Psi), for sigma = e^{-theta_d}, theta_d = 0.5. It is a sum of a function of
Phi + Psi and a function of Phi - Psi, a wave in field space. The pairs of bulk vacua (2 pi m, 2 pi n) are its
critical points: even jumps are minima, odd jumps maxima.
(b) The energy and momentum (B, P_d) stored by a static jump of 2 pi, in units of M = 8 m0/beta^2, for several
theta_d: they lie on the soliton mass shell E^2 - P^2 = M^2, at the point of rapidity theta_d.

Checks: the critical points and their type (Hessian); the plotted potential equals f(Phi+Psi) + g(Phi-Psi);
the stored (B, P_d) equal M (cosh theta_d, sinh theta_d).
Writes lecture-notes/defects/figures/l3_landscape.svg.
"""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects" / "figures"

m0 = 1.0
thd = 0.5
sig = math.exp(-thd)


def Bpot(ph, ps, sg=sig):
    return -2 * m0 * (sg * (np.cos((ph + ps) / 2) - 1) + (np.cos((ph - ps) / 2) - 1) / sg)


def Ppot(ph, ps, sg=sig):
    return 2 * m0 * (sg * (np.cos((ph + ps) / 2) - 1) - (np.cos((ph - ps) / 2) - 1) / sg)


# critical points and their type, by finite differences
h = 1e-4
ok = True
for m in range(-1, 2):
    for n in range(-1, 2):
        p0, q0 = 2 * math.pi * m, 2 * math.pi * n
        gx = (Bpot(p0 + h, q0) - Bpot(p0 - h, q0)) / (2 * h)
        gy = (Bpot(p0, q0 + h) - Bpot(p0, q0 - h)) / (2 * h)
        hxx = (Bpot(p0 + h, q0) - 2 * Bpot(p0, q0) + Bpot(p0 - h, q0)) / h**2
        hyy = (Bpot(p0, q0 + h) - 2 * Bpot(p0, q0) + Bpot(p0, q0 - h)) / h**2
        hxy = (Bpot(p0 + h, q0 + h) - Bpot(p0 + h, q0 - h) - Bpot(p0 - h, q0 + h) + Bpot(p0 - h, q0 - h)) / (4 * h**2)
        ev = np.linalg.eigvalsh(np.array([[hxx, hxy], [hxy, hyy]]))
        even = (n - m) % 2 == 0
        ok &= abs(gx) < 1e-8 and abs(gy) < 1e-8 and (all(ev > 0) if even else all(ev < 0))
report("pairs of bulk vacua are critical points of B: even jumps minima, odd jumps maxima", ok)

rng = np.random.default_rng(5)
pts = rng.uniform(-7, 7, size=(20, 2))
f = lambda sv: -2 * m0 * sig * (np.cos(sv / 2) - 1)
g = lambda dv: -2 * m0 / sig * (np.cos(dv / 2) - 1)
report("the plotted potential is f(Phi+Psi) + g(Phi-Psi)",
       max(abs(Bpot(p, q) - f(p + q) - g(p - q)) for p, q in pts) < 1e-12)

M = 8 * m0                                       # beta^2 = 1
thetas = [-1.0, 0.0, 1.0]
ok = True
for th in thetas:
    sg = math.exp(-th)
    Eb, Pb = Bpot(0.0, 2 * math.pi, sg), Ppot(0.0, 2 * math.pi, sg)
    ok &= abs(Eb - M * math.cosh(th)) < 1e-12 and abs(Pb - M * math.sinh(th)) < 1e-12
report("a static jump of 2 pi stores (B, P_d) = M (cosh theta_d, sinh theta_d) for the plotted theta_d", ok)

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.2), gridspec_kw=dict(width_ratios=[1.0, 1.0]))
L = 3.2 * math.pi
ph = np.linspace(-L, L, 401)
PH, PS = np.meshgrid(ph, ph, indexing="xy")
Z = Bpot(PH, PS) / m0
cs = ax.contourf(PH, PS, Z, levels=14, cmap="RdBu_r", alpha=0.85)
ax.contour(PH, PS, Z, levels=14, colors=style.GREY, linewidths=0.3)
for m in range(-1, 2):
    for n in range(-1, 2):
        p0, q0 = 2 * math.pi * m, 2 * math.pi * n
        even = (n - m) % 2 == 0
        ax.plot(p0, q0, marker="o", ms=5, mfc=style.BLUE if even else style.RED, mec="white", mew=0.6, ls="none")
# characteristic directions of the wave equation in field space
for c in (-2 * math.pi, 0.0, 2 * math.pi):
    ax.plot(ph, c - ph, ls=":", color=style.GREY, lw=0.8)
    ax.plot(ph, ph + c, ls=":", color=style.GREY, lw=0.8)
ax.set_xlim(-L, L)
ax.set_ylim(-L, L)
ticks = [-2 * math.pi, 0, 2 * math.pi]
ax.set_xticks(ticks, [r"$-2\pi$", "0", r"$2\pi$"])
ax.set_yticks(ticks, [r"$-2\pi$", "0", r"$2\pi$"])
ax.set_xlabel(r"$\Phi(0,t)$")
ax.set_ylabel(r"$\Psi(0,t)$")
ax.set_aspect("equal")
from matplotlib.lines import Line2D  # noqa: E402
handles = [Line2D([], [], marker="o", ms=5, mfc=style.BLUE, mec="white", ls="none", label="even jump: minimum"),
           Line2D([], [], marker="o", ms=5, mfc=style.RED, mec="white", ls="none", label="odd jump: maximum")]
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=7.5,
          handletextpad=0.2, columnspacing=1.0)
cb = fig.colorbar(cs, ax=ax, shrink=0.8, pad=0.03)
cb.ax.tick_params(labelsize=7, colors=style.GREY)
cb.set_label(r"$\mathcal{B}/m_0$", color=style.GREY)
ax.set_title(r"(a) the defect potential, $\theta_d=0.5$", fontsize=8.5, color=style.GREY)

thg = np.linspace(-1.75, 1.75, 200)
ax2.plot(np.sinh(thg), np.cosh(thg), color=style.GREY, lw=1.2)
pl = np.array([-2.9, 0.0, 2.9])
ax2.plot(pl, np.abs(pl), ls="--", color=style.GREY, lw=0.6)
for th, off in ((-1.0, (-46, 4)), (0.0, (-10, 9)), (1.0, (8, -10))):
    sg = math.exp(-th)
    Eb, Pb = Bpot(0.0, 2 * math.pi, sg) / M, Ppot(0.0, 2 * math.pi, sg) / M
    ax2.plot(Pb, Eb, marker="o", ms=5.5, color=style.RED, ls="none")
    ax2.annotate(rf"$\theta_d={th:g}$", (Pb, Eb), textcoords="offset points", xytext=off,
                 fontsize=7.5, color=style.RED)
ax2.plot(0, 0, marker="o", ms=5.5, color=style.BLUE, ls="none", clip_on=False)
ax2.annotate("no jump", (0, 0), textcoords="offset points", xytext=(6, 4), fontsize=7.5, color=style.BLUE)
ax2.text(0.0, 2.55, r"mass shell $E^2-P^2=M^2$", color=style.GREY, fontsize=8, ha="center")
ax2.set_xlim(-2.9, 2.9)
ax2.set_ylim(0, 2.9)
ax2.set_xlabel(r"stored momentum $/M$")
ax2.set_ylabel(r"stored energy $/M$")
ax2.set_title(r"(b) a jump of $2\pi$ is a frozen soliton", fontsize=8.5, color=style.GREY)
fig.tight_layout(w_pad=1.2)
style.save(fig, "l3_landscape", OUT)
raise SystemExit(0 if style.all_ok() else 1)
