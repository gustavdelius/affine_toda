"""Figure F35 (chapter 28, sec-reflection-equation, eq-reflection-braid): the two orders of two reflections.

Two particles with rapidities theta_1 > theta_2 > 0 (velocities v1 > v2 > 0) move towards the boundary at x = 0,
particle 2 nearer to it. Each world line is straight until it meets the wall, where the particle reflects with
velocity -v_i. Depending on the initial positions, either
  (a) particle 2 reflects, the two scatter at rapidity difference theta_1 + theta_2, particle 1 reflects,
      and they scatter again at theta_1 - theta_2:  S(th1 - th2) K_2(th1) S(th1 + th2) K_2(th2), the left-hand
      side of eq-reflection-braid read from the right; or
  (b) they scatter first at theta_1 - theta_2, particle 1 reflects, they scatter at theta_1 + theta_2, and
      particle 2 reflects:  K_2(th2) S(th1 + th2) K_2(th1) S(th1 - th2), the right-hand side.
Factorization says the two amplitudes agree, which is the reflection equation (Ghoshal-Zamolodchikov Fig. 10).

Checks: the event order in each panel is the one stated, each crossing happens at x < 0, and the crossings are
between particles moving in opposite directions (theta_1 + theta_2) or in the same direction (theta_1 - theta_2)
as labelled.

Writes part7-half-line/figures/ch28_f35_reflection_equation.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
v = np.array([0.8, 0.4])                      # v1 > v2 > 0
colors = [style.BLUE, style.ORANGE]
TMAX = 13.0


def path(x0, i, t):
    """position of particle i at time t (reflection at x = 0)"""
    tw = -x0 / v[i]
    return np.where(t < tw, x0 + v[i] * t, -v[i] * (t - tw))


def events(x0):
    """time-ordered list of (t, x, kind, label)"""
    ev = []
    for i in range(2):
        ev.append((-x0[i] / v[i], 0.0, "K", i + 1))
    tt = np.linspace(0, TMAX, 260001)
    d = path(x0[0], 0, tt) - path(x0[1], 1, tt)
    sg = np.sign(d)
    sg[sg == 0] = 1                                   # a crossing exactly on a grid point counts once
    for k in np.nonzero(sg[1:] != sg[:-1])[0]:
        t_ = tt[k]
        # velocities at the crossing: same sign -> theta_1 - theta_2, opposite -> theta_1 + theta_2
        s1 = np.sign(path(x0[0], 0, t_ + 1e-6) - path(x0[0], 0, t_ - 1e-6))
        s2 = np.sign(path(x0[1], 1, t_ + 1e-6) - path(x0[1], 1, t_ - 1e-6))
        ev.append((t_, float(path(x0[0], 0, t_)), "S", "-" if s1 == s2 else "+"))
    return sorted(ev)


cfgA = np.array([-4.0, -1.0])
cfgB = np.array([-4.0, -3.0])
eA, eB = events(cfgA), events(cfgB)
kindsA = [(k, l) for _, _, k, l in eA]
kindsB = [(k, l) for _, _, k, l in eB]
report("(a): K(theta_2), S(theta_1+theta_2), K(theta_1), S(theta_1-theta_2) -- the left-hand side of eq-reflection-braid",
       kindsA == [("K", 2), ("S", "+"), ("K", 1), ("S", "-")])
report("(b): S(theta_1-theta_2), K(theta_1), S(theta_1+theta_2), K(theta_2) -- the right-hand side",
       kindsB == [("S", "-"), ("K", 1), ("S", "+"), ("K", 2)])
report("all scatterings happen in front of the wall (x < 0)",
       all(x < 0 for _, x, k, _ in eA + eB if k == "S"))

fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.2), sharey=True)
tt = np.linspace(0, TMAX, 600)
lab = {("K", 1): r"$K(\theta_1)$", ("K", 2): r"$K(\theta_2)$",
       ("S", "+"): r"$S(\theta_1{+}\theta_2)$", ("S", "-"): r"$S(\theta_1{-}\theta_2)$"}
for ax, x0, ev, title in ((axes[0], cfgA, eA, "(a) particle 2 reflects first"),
                          (axes[1], cfgB, eB, "(b) the particles scatter first")):
    ax.axvspan(0, 1.2, color=style.LIGHTGREY, lw=0)
    ax.axvline(0, color=style.GREY, lw=1.4)
    for i in range(2):
        ax.plot(path(x0[i], i, tt), tt, color=colors[i], lw=1.6)
        ax.text(x0[i], -0.75, rf"${i + 1}$", color=colors[i], ha="center", fontsize=10)
    for n, (t_, x_, k, l) in enumerate(ev, 1):
        ax.plot([x_], [t_], "o", color="white", mec=style.GREY, ms=6.5, zorder=5)
        if k == "K":
            ax.text(0.1, t_, lab[(k, l)], fontsize=8, va="center")
        else:
            ax.text(x_ - 0.25, t_, lab[(k, l)], fontsize=8, va="center", ha="right")
    ax.set_title(title, fontsize=9, color=style.GREY)
    ax.set_xlim(-6.2, 1.2)
    ax.set_ylim(-1.2, TMAX)
    ax.set_xlabel(r"$x$")
    ax.set_xticks([])
axes[0].set_ylabel(r"$t$")
axes[0].set_yticks([])
fig.text(0.5, -0.03, "factorization: the amplitudes of (a) and (b) are equal, which is the reflection equation",
         ha="center", fontsize=8.5, color=style.GREY)
fig.tight_layout()
style.save(fig, "ch28_f35_reflection_equation", style.OUT_PART7)
raise SystemExit(0 if style.all_ok() else 1)
