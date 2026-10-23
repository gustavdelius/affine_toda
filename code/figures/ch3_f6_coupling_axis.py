"""Figure F6 (chapter 3, sec-sg-coleman and sec-sg-smatrix): the coupling axis of sine-Gordon.

x = beta_SG^2 / 8 pi.  With lambda = 8 pi / beta_SG^2 - 1 = 1/x - 1 and xi = pi / lambda =
pi x / (1 - x):
  * the breather B_n exists iff n < lambda, i.e. below x = 1/(n+1)  (beta_SG^2 = 8 pi/(n+1));
  * x = 1/2 (beta_SG^2 = 4 pi, lambda = 1, xi = pi): free-fermion point, last breather gone,
    first collapse threshold of the Coulomb gas;
  * x = 1 (beta_SG^2 = 8 pi): Coleman's bound, lambda -> 0, xi -> infinity;
  * RSOS points beta_SG^2/8pi = p/p': the scaling Lee-Yang model is (p,p') = (2,5), x = 2/5;
    the unitary minimal models M(p,p+1) sit at x = p/(p+1) in the repulsive regime.

Writes part1-background/figures/ch3_f6_coupling_axis.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
lam = lambda x: 1 / x - 1
xi = lambda x: np.pi * x / (1 - x)
n_breathers = lambda x: int(np.ceil(lam(x)) - 1) if lam(x) > 1 else 0   # integers n with n < lambda


def count_ticks(x):                      # thresholds 1/(n+1) lying above x
    return sum(1 for n in range(1, 60) if x < 1 / (n + 1) - 1e-12)


# ---- checks ---------------------------------------------------------------
report("number of breathers (n < lambda) = number of thresholds 8 pi/(n+1) above beta^2",
       all(n_breathers(x) == count_ticks(x) for x in np.linspace(0.02, 0.99, 400)))
report("free fermion x=1/2: lambda=1, xi=pi, no breather",
       abs(lam(0.5) - 1) < 1e-12 and abs(xi(0.5) - np.pi) < 1e-12 and n_breathers(0.5) == 0)
report("thresholds x=1/(n+1) have xi = pi/n", all(abs(xi(1 / (n + 1)) - np.pi / n) < 1e-12
                                                  for n in range(1, 8)))
report("scaling Lee-Yang (p,p')=(2,5): x=2/5, xi=2 pi/3, exactly one breather",
       abs(xi(0.4) - 2 * np.pi / 3) < 1e-12 and n_breathers(0.4) == 1)
report("unitary series M(p,p+1) lies in the repulsive regime x>1/2, no breathers",
       all(p / (p + 1) > 0.5 and n_breathers(p / (p + 1)) == 0 for p in range(2, 12)))

# ---- figure ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.2, 3.0))
ax.set_xlim(-0.14, 1.2)
ax.set_ylim(-1.5, 1.75)
ax.axis("off")
ax.plot([0, 1.0], [0, 0], color=style.GREY, lw=3)
ax.fill_between([0, 0.5], 0, 0.55, color=style.BLUE, alpha=0.16, lw=0)
ax.fill_between([0.5, 1.0], 0, 0.55, color=style.ORANGE, alpha=0.16, lw=0)
ax.fill_between([1.0, 1.12], 0, 0.55, facecolor="none", edgecolor=style.GREY, hatch="///", lw=0)
ax.text(0.25, 0.43, "attractive", ha="center", color=style.BLUE, fontsize=8)
ax.text(0.75, 0.43, "repulsive: no breathers", ha="center", color=style.ORANGE, fontsize=8)
ax.text(1.06, 0.43, "massless", ha="center", color=style.GREY, fontsize=7)
for n in range(1, 6):                       # number of breathers on each interval
    lo, hi = 1 / (n + 2), 1 / (n + 1)
    ax.text(0.5 * (lo + hi), 0.17, str(n), ha="center", color=style.BLUE, fontsize=8)
ax.text(-0.02, 0.17, "breathers:", ha="right", color=style.BLUE, fontsize=7.5)
for n in range(1, 7):                       # thresholds
    xt = 1 / (n + 1)
    ax.plot([xt, xt], [-0.1, 0.1], color=style.GREY, lw=1.2)
    if n > 4:                               # too close together to label
        continue
    ax.text(xt, -0.25, [r"$4\pi$", r"$8\pi/3$", r"$2\pi$", r"$8\pi/5$"][n - 1], ha="center", fontsize=7.5)
    ax.text(xt, -0.72, r"$\pi$" if n == 1 else rf"$\pi/{n}$", ha="center", fontsize=7.5)
    ax.text(xt, -1.18, str(n), ha="center", fontsize=7.5)
for xt, lab, lab2, lab3 in ((0.0, r"$0$", r"$0$", r"$\infty$"), (1.0, r"$8\pi$", r"$\infty$", "0")):
    ax.plot([xt, xt], [-0.1, 0.1], color=style.GREY, lw=1.2)
    ax.text(xt, -0.25, lab, ha="center", fontsize=7.5)
    ax.text(xt, -0.72, lab2, ha="center", fontsize=7.5)
    ax.text(xt, -1.18, lab3, ha="center", fontsize=7.5)
ax.text(1.14, -0.25, r"$\beta_{\rm SG}^2$", ha="left", fontsize=8.5)
ax.text(1.14, -0.72, r"$\xi$", ha="left", fontsize=8.5)
ax.text(1.14, -1.18, r"$\lambda$", ha="left", fontsize=8.5)
ax.plot([0.5, 0.5], [-0.1, 1.15], "--", color=style.GREY, lw=1.0)
ax.text(0.5, 1.22, "free fermion", ha="center", fontsize=7.5)
ax.plot([1.0, 1.0], [-0.1, 1.15], "--", color=style.GREY, lw=1.0)
ax.text(1.0, 1.22, "Coleman's bound", ha="center", fontsize=7.5)
ax.plot([0.4], [0.0], "*", color=style.RED, ms=11, zorder=5)
ax.text(0.4, 0.78, "Lee–Yang\n" r"$\mathcal{M}(2,5)$", ha="center", color=style.RED, fontsize=7.5)
for p in range(2, 10):
    ax.plot([p / (p + 1)], [0.0], "D", color=style.GREEN, ms=4, zorder=5)
ax.text(0.79, 0.78, r"unitary $\mathcal{M}(p,p+1)$", ha="center", color=style.GREEN, fontsize=7.5)
ax.text(-0.14, 1.62, r"$x=\beta_{\rm SG}^2/8\pi$;  $B_n$ exists for $n<\lambda=8\pi/\beta_{\rm SG}^2-1$",
        ha="left", fontsize=8)
fig.tight_layout()
style.save(fig, "ch3_f6_coupling_axis", style.OUT_PART1)
raise SystemExit(0 if style.all_ok() else 1)
