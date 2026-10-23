"""Lecture notes on defects, Lecture 1, figure 3 (fig-histories): two particles, two histories.

Particles a (rapidity theta1) and b (rapidity theta2 < theta1) arrive from the left, a behind b.  Their initial
positions are fixed by the times t_a, t_b at which they would reach the defect at x = 0.  Shifting b's initial
position (what the higher conserved charges do) changes the order of events.

Top row, out-state "b reflected, a transmitted":
  (i)  t_b < t_a: b reflects, meets a head-on at rapidity difference theta1 + theta2, a is transmitted;
  (ii) t_a < t_b: a overtakes b at rapidity difference theta1 - theta2, a is transmitted, b reflects.
Bottom row, out-state "both transmitted":
  (iii) t_a < t_b: a overtakes b on the left (theta1 - theta2), then a and b are transmitted;
  (iv)  t_b < t_a: b and a are transmitted, then a overtakes b on the right (theta1 - theta2).

Checks (computed from the world lines, not read off the labels): the event order in each panel; the rapidity
difference at each collision; on which side of the defect each collision happens.
Writes lecture-notes/defects/figures/l1_histories.svg.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects" / "figures"

th1, th2 = 1.0, 0.45
v1, v2 = np.tanh(th1), np.tanh(th2)
T0, T1 = -2.4, 2.6     # time window


def world(v, t_hit, fate):
    """Piecewise world line: x = v (t - t_hit) before, then +-v (t - t_hit) after (transmit / reflect)."""
    def x(t):
        t = np.asarray(t, dtype=float)
        after = v if fate == "T" else -v
        return np.where(t < t_hit, v * (t - t_hit), after * (t - t_hit))
    rap_after = (th1 if v == v1 else th2) * (1 if fate == "T" else -1)
    return x, rap_after


def events(t_a, t_b, fate_a, fate_b):
    xa, ra = world(v1, t_a, fate_a)
    xb, rb = world(v2, t_b, fate_b)
    ev = [(t_a, "a " + ("transmitted" if fate_a == "T" else "reflected"), 0.0, None),
          (t_b, "b " + ("transmitted" if fate_b == "T" else "reflected"), 0.0, None)]
    ts = np.linspace(T0, T1, 200001)
    d = xa(ts) - xb(ts)
    for i in np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]:
        tc = ts[i] - d[i] * (ts[i + 1] - ts[i]) / (d[i + 1] - d[i])
        # rapidities at the collision
        rA = th1 if tc < t_a else ra
        rB = th2 if tc < t_b else rb
        ev.append((tc, "a meets b", float(xa(tc)), abs(rA - rB)))
    return sorted(ev), xa, xb


panels = [
    ("(i) b meets the defect first", -0.3, -1.2, "T", "R"),
    ("(ii) the particles meet first", -0.6, 0.2, "T", "R"),
    ("(iii) they meet on the left", -0.6, 0.2, "T", "T"),
    ("(iv) they meet on the right", 0.3, -0.9, "T", "T"),
]
expect = [
    (["b reflected", "a meets b", "a transmitted"], th1 + th2, "left"),
    (["a meets b", "a transmitted", "b reflected"], th1 - th2, "left"),
    (["a meets b", "a transmitted", "b transmitted"], th1 - th2, "left"),
    (["b transmitted", "a transmitted", "a meets b"], th1 - th2, "right"),
]
fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.4))
for ax, (title, t_a, t_b, fa, fb), (order, drap, side) in zip(axes.flat, panels, expect):
    ev, xa, xb = events(t_a, t_b, fa, fb)
    names = [e[1] for e in ev]
    coll = [e for e in ev if e[1] == "a meets b"]
    ok = names == order and len(coll) == 1
    ok = ok and abs(coll[0][3] - drap) < 1e-9 and (coll[0][2] < 0) == (side == "left")
    report(f"{title}: event order {names}, collision at rapidity difference {coll[0][3]:.3f}, on the {side}", ok)
    ts = np.linspace(T0, T1, 400)
    ax.plot([0, 0], [T0, T1], color=style.PURPLE, lw=2.0)
    if "R" in (fa, fb):
        ax.plot([0, 0], [T0, T1], color=style.ORANGE, lw=2.0, ls=(0, (3, 3)))
    ax.plot(xa(ts), ts, color=style.BLUE, lw=1.8)
    ax.plot(xb(ts), ts, color=style.RED, lw=1.8)
    tc, _, xc, _ = coll[0]
    ax.plot([xc], [tc], "o", ms=5, color=style.GREY)
    lab = r"$S(\theta_1+\theta_2)$" if abs(coll[0][3] - (th1 + th2)) < 1e-9 else r"$S(\theta_1-\theta_2)$"
    # label in empty space: right of the defect and below the collision, or left of it and above
    txy = (0.15, tc - 0.45) if xc < 0 else (0.12, tc + 0.3)
    ax.annotate(lab, xy=(xc, tc), xytext=txy, fontsize=8, va="center",
                arrowprops=dict(arrowstyle="-", color=style.GREY, lw=0.6))
    tl = T0 + 0.55
    ax.text(float(xa(tl)) - 0.32, tl, "$a$", color=style.BLUE, fontsize=9, va="center")
    ax.text(float(xb(tl)) + 0.12, tl, "$b$", color=style.RED, fontsize=9, va="center")
    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(T0, T1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("$x$", labelpad=1)
    ax.set_ylabel("$t$", labelpad=1)
    ax.set_title(title, fontsize=8.5, color=style.GREY)
fig.text(0.5, 0.985, r"out-state: $b$ reflected, $a$ transmitted", ha="center", va="top", fontsize=8.5)
fig.text(0.5, 0.505, r"out-state: both transmitted", ha="center", va="top", fontsize=8.5)
fig.tight_layout(rect=(0, 0, 1, 0.97), h_pad=2.2, w_pad=1.5)
style.save(fig, "l1_histories", OUT)
raise SystemExit(0 if style.all_ok() else 1)
