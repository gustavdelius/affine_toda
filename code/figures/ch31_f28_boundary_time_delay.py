"""Figure F28 (chapter 31, sec-images-boundary, eq-mirror-shift, eq-boundary-time-delay, eq-virtual-collision):
space-time picture of a soliton reflected as its antisoliton, for C = +1, Neumann, C = -1.

The reflected soliton is read from the whole-line mirror pair E_1 = exp(m_a(x cosh th - t sinh th) + xi_1) (incoming) and
E_2 = exp(m_a(x cosh th + t sinh th) + xi_2) (mirror), with xi_1 = 0 and xi_1 + xi_2 = S, S the shift of eq-mirror-shift.
Before the collision the mirror is seen with parameter xi_2 + ln A (A = 1 - s^2, eq-mirror-pair); after it, the outgoing
antisoliton with xi_2.  Hence the outgoing line is x = -v(t - Delta t) with Delta t = -S/(m_a sinh th), and the two incoming
lines meet at the virtual collision point (d, d/v) with d = -(S + ln A)/(2 m_a cosh th).

Checks: A = 1 - s^2 is the interaction coefficient of the pair; Delta t and d computed from S equal the closed forms of
eq-boundary-time-delay and eq-virtual-collision; Delta t_N = (Delta t_+ + Delta t_-)/2; signs; the incoming lines meet at (d, d/v).
The shifts S themselves are verified against the boundary condition in code/boundary_images.py.

Writes part7-half-line/figures/ch31_f28_boundary_time_delay.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
h, a, th = 3, 1, 0.6                                # a_2^(1), species 1; m = 1, m_a = 2 sin(pi a/h)
ma = 2 * np.sin(np.pi * a / h)
v = np.tanh(th)
s = np.sin(np.pi * a / h) / np.cosh(th)
A = 1 - s**2


def shift(C):
    return -np.log(A) if C == 0 else -2 * np.log(1 + C * s)


def delay(C):
    return -shift(C) / (ma * np.sinh(th))


def vcp(C):
    return -(shift(C) + np.log(A)) / (2 * ma * np.cosh(th))


# ---- checks ----------------------------------------------------------------
def interaction(a_, b_, th12, h_):
    A_, B_ = np.pi * a_ / h_, np.pi * b_ / h_
    return (np.cosh(th12) - np.cos(A_ - B_)) / (np.cosh(th12) - np.cos(A_ + B_))


ok = True
for h_, a_, t_ in ((3, 1, 0.6), (4, 1, 0.3), (4, 2, 1.1), (6, 2, 0.45), (6, 3, 0.9)):
    s_ = np.sin(np.pi * a_ / h_) / np.cosh(t_)
    ok &= abs(interaction(a_, h_ - a_, 2 * t_, h_) - (1 - s_**2)) < 1e-12
report("interaction coefficient A_{a,abar}(2 theta) = 1 - sin^2(pi a/h)/cosh^2(theta) (eq-mirror-pair), five cases", ok)
dt_closed = {C: (2 * np.log(1 + C * s) if C else np.log(1 - s**2)) / (ma * np.sinh(th)) for C in (1, 0, -1)}
d_closed = {C: C / (2 * ma * np.cosh(th)) * np.log((1 + s) / (1 - s)) for C in (1, 0, -1)}
report("Delta t from the shift S equals eq-boundary-time-delay for C = +1, 0, -1", all(abs(delay(C) - dt_closed[C]) < 1e-12 for C in (1, 0, -1)))
report("d from the shift S equals eq-virtual-collision for C = +1, 0, -1", all(abs(vcp(C) - d_closed[C]) < 1e-12 for C in (1, 0, -1)))
report("Delta t_N = (Delta t_{+1} + Delta t_{-1})/2", abs(dt_closed[0] - (dt_closed[1] + dt_closed[-1]) / 2) < 1e-12)
report("C = +1: delay and collision behind the wall (Delta t > 0, d > 0); Neumann: advance, d = 0; C = -1: advance, d < 0",
       dt_closed[1] > 0 and d_closed[1] > 0 and dt_closed[0] < 0 and abs(d_closed[0]) < 1e-14 and dt_closed[-1] < 0 and d_closed[-1] < 0)
ok = True
for C in (1, 0, -1):
    a1, a2 = 0.0, (shift(C) + np.log(A)) / (ma * np.cosh(th))      # incoming and mirror lines: x = v t - a1, x = -v t - a2
    tstar = (a1 - a2) / (2 * v)
    ok &= abs(v * tstar - a1 - vcp(C)) < 1e-12 and abs(tstar - vcp(C) / v) < 1e-12
report("the incoming soliton line and the mirror line meet at the virtual collision point (d, d/v)", ok)

# ---- figure ----------------------------------------------------------------
fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.9), sharey=True)
tt = np.linspace(-4.0, 4.0, 2)
for ax, C, title in zip(axs, (1, 0, -1), (r"$C=+1$ (attractive)", "Neumann", r"$C=-1$ (repulsive)")):
    dt, d = delay(C), vcp(C)
    xin = v * np.linspace(-4.5, 0, 2)                              # incoming soliton, x = v t for t < 0 (reaches x = 0 at t = 0)
    ax.axvspan(0, 1.9, color=style.LIGHTGREY, lw=0)
    ax.axvline(0, color=style.GREY, lw=1.4)
    t_in = np.linspace(-4.5, d / v + 0.0, 2) if True else None
    ax.plot(v * np.array([-4.5, 4.5]), [-4.5, 4.5], color=style.GREY, lw=.8, ls=":")             # undelayed incoming/outgoing reference
    ax.plot(-v * np.array([-4.5, 4.5]), [-4.5, 4.5], color=style.GREY, lw=.8, ls=":")
    # incoming soliton
    tmax = d / v if C else 0.0
    t1 = np.array([-4.5, tmax]); ax.plot(v * t1, t1, color=style.BLUE, lw=1.8)
    # mirror (whole-line picture) to the right of the wall, before the collision
    a2 = (shift(C) + np.log(A)) / (ma * np.cosh(th))
    tm = np.array([tmax, 4.5])
    ax.plot(-v * np.array([-4.5, tmax]) - a2, [-4.5, tmax], color=style.GREY, lw=1.2, ls="--")
    # outgoing antisoliton on x < 0
    t_out = np.linspace(dt if d > 1e-9 else d / v, 4.5, 2)         # for d < 0 the antisoliton emerges at the collision point
    ax.plot(-v * (t_out - dt), t_out, color=style.RED, lw=1.8)
    ax.plot([d], [d / v], "o", color=style.ORANGE, ms=5, zorder=5)
    ax.set_xlim(-3.0, 1.9)
    ax.set_ylim(-4.5, 4.5)
    ax.set_xlabel("$x$")
    ax.set_title(title + "\n" + rf"$\Delta t={dt:+.2f},\ d={d + 0.0:+.2f}$".replace("-0.00", "0").replace("+0.00", "0"), fontsize=8)
axs[0].set_ylabel("$t$")
axs[0].plot([], [], color=style.BLUE, label="soliton")
axs[0].plot([], [], color=style.RED, label="reflected antisoliton")
axs[0].plot([], [], color=style.GREY, ls="--", label="mirror (image)")
axs[0].plot([], [], color=style.GREY, ls=":", lw=.8, label="no delay")
axs[0].plot([], [], "o", color=style.ORANGE, ms=4, label="virtual collision")
fig.legend(*axs[0].get_legend_handles_labels(), loc="lower center", ncol=5, fontsize=7, columnspacing=1.0, handlelength=1.5, bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(rect=(0, 0.07, 1, 1))
style.save(fig, "ch31_f28_boundary_time_delay", style.OUT_PART7)
raise SystemExit(0 if style.all_ok() else 1)
