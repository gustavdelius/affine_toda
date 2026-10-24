"""Figure F21 (chapter 16, sec-transmission-factors, sec-dhn-examples): zeros of X and the one-loop mass.

(a) The transmission factor X(z) = prod (z - z_j)/(z + z_j) of the sine-Gordon and phi^4 kinks: zeros at
    z_j = sin t_j (bound states, omega_j = m cos t_j), poles at -z_j.
(b) The contribution -(m/pi)(sin t - t cos t) of one bound state to Delta M (eq-ccg) as a function of t,
    with the bound states of the two kinks marked. The zero mode (t = pi/2) is the most expensive one.

Checks: eq-ccg equals eq-dhn-psi; Delta M = -m/pi (sine-Gordon) and m(1/(4 sqrt 3) - 3/(2 pi)) (phi^4);
the bound-state frequencies omega = m cos t_j are 0 and sqrt(3)/2 m.
Writes part5-semiclassics/figures/ch17_f21_bound_states.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
varpi = lambda t: (t * np.cos(t) - np.sin(t)) / np.pi
ccg = lambda ts: -sum(np.sin(t) - t * np.cos(t) for t in ts) / np.pi


def dM_psi(zeros, poles):                     # eq-dhn-psi with m = 1
    return 0.5 * (sum(varpi(np.arcsin(w)) for w in zeros) - sum(varpi(np.arcsin(w)) for w in poles))


kinks = {"sine-Gordon": [1.0], r"$\phi^4$": [1.0, 0.5]}
# ---- checks ---------------------------------------------------------------
for name, zs in kinks.items():
    ts = np.arcsin(zs)
    report(f"{name}: eq-ccg = eq-dhn-psi", abs(ccg(ts) - dM_psi(zs, [-z for z in zs])) < 1e-14)
report("sine-Gordon: Delta M = -m/pi", abs(ccg(np.arcsin([1.0])) + 1 / np.pi) < 1e-14)
report("phi^4: Delta M = m(1/(4 sqrt3) - 3/(2 pi))",
       abs(ccg(np.arcsin([1.0, 0.5])) - (1 / (4 * np.sqrt(3)) - 3 / (2 * np.pi))) < 1e-14)
report("phi^4: bound-state frequencies 0 and sqrt(3)/2 m",
       np.allclose([np.cos(np.arcsin(z)) for z in (1.0, 0.5)], [0, np.sqrt(3) / 2], atol=1e-12))

# ---- figure ---------------------------------------------------------------
fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 3.0), gridspec_kw={"width_ratios": [1.15, 1]})
a.axhline(0, color=style.GREY, lw=0.8)
for row, (name, zs) in enumerate(kinks.items()):
    y = 1.8 * (1 - row)
    for z in zs:
        a.plot([z], [y], "o", mfc="none", mec=style.BLUE, mew=1.6, ms=8)
        a.plot([-z], [y], "x", color=style.RED, mew=1.6, ms=7)
        a.text(z, y + 0.17, rf"$z={z:g}$" + ("\n(zero mode)" if z == 1 else r"" + "\n" + r"$\omega=\frac{\sqrt{3}}{2}m$"),
               ha="center", va="bottom", fontsize=7, color=style.GREY)
    a.text(-1.28, y, name, ha="right", va="center", fontsize=8, color=style.GREY)
a.hlines([0, 1.8], -1.2, 1.2, color=style.LIGHTGREY, lw=0.8)
a.plot([], [], "o", mfc="none", mec=style.BLUE, mew=1.6, label=r"zero $z_j=\sin t_j$")
a.plot([], [], "x", color=style.RED, mew=1.6, label=r"pole $-z_j$")
a.legend(loc="lower center", ncol=2, fontsize=7.5, bbox_to_anchor=(0.5, -0.48), handletextpad=0.2, columnspacing=1.0)
a.set_xlim(-1.2, 1.2)
a.set_ylim(-0.55, 2.95)
a.set_yticks([])
a.spines["left"].set_visible(False)
a.set_xlabel(r"$z=\kappa/m$", labelpad=1)
a.set_title(r"(a) zeros and poles of $X(z)$", fontsize=9, color=style.GREY)

t = np.linspace(0, np.pi / 2, 200)
b.plot(t, -(np.sin(t) - t * np.cos(t)) / np.pi, color=style.GREY, lw=1.4)
for ts, col, nm in ((np.pi / 2, style.RED, "zero mode"), (np.pi / 6, style.BLUE, r"$\omega=\frac{\sqrt{3}}{2}m$")):
    v = -(np.sin(ts) - ts * np.cos(ts)) / np.pi
    b.plot([ts], [v], "o", color=col, ms=5)
    b.annotate(nm, (ts, v), (ts - (0.45 if ts > 1 else -0.12), v + (0.05 if ts > 1 else -0.06)), fontsize=7.5, color=col)
b.set_xticks([0, np.pi / 6, np.pi / 2])
b.set_xticklabels(["0", r"$\pi/6$", r"$\pi/2$"])
b.set_xlabel(r"$t_j$ ($\omega_j=m\cos t_j$)")
b.set_ylabel(r"$\Delta M_j/m$")
b.set_title("(b) share of one bound state", fontsize=9, color=style.GREY)
fig.tight_layout()
style.save(fig, "ch17_f21_bound_states", style.OUT_PART5)
raise SystemExit(0 if style.all_ok() else 1)
