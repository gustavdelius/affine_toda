"""Figure F23 (chapter 23, sec-cn-parameters, sec-bound-states-exact-masses): pole ladders of c_n^(1).

(a) Poles of S_{nn} in the strip variable t (physical strip 0 < t < T) for c_3^(1) at omega = 2.37:
    bound-state poles a_j = j omega + (j-1)/2, spaced by omega + 1/2, with a_n = T, and the breather
    tower T - k, spaced by 1. The crossed images t -> T - t are drawn in the second row.
(b) Fusion masses M_{n-j}/M_n = 2 cos(pi a_j / 2T) = 2 sin((n-j) pi / H) against omega, with H = 2n + B.

Checks: a_n = T; T - a_j = (n-j)(omega + 1/2); 2T = H (omega + 1/2) with H = 2n + B, B = -1/(omega+1/2)
= -2 beta^2/(8 pi - beta^2); the fusion mass 2 cos(pi a_j/2T) equals 2 sin((n-j) pi/H) (n = 2..6, random omega);
the c_3 poles at omega and 2 omega + 1/2 (table of sec-c3-spectrum).
Writes part6-soliton-smatrices/figures/ch24_f23_cn_ladders.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
aj = lambda j, w: j * w + (j - 1) / 2
Tn = lambda n, w: n * w + (n - 1) / 2

# ---- checks ---------------------------------------------------------------
rng = np.random.default_rng(1)
ok = [True] * 4
for _ in range(30):
    w = rng.uniform(0.3, 6)
    b2 = 4 * np.pi / (w + 1)
    for n in range(2, 7):
        T = Tn(n, w)
        H = 2 * n - 1 / (w + 0.5)
        ok[0] &= abs(aj(n, w) - T) < 1e-12
        ok[1] &= all(abs(T - aj(j, w) - (n - j) * (w + 0.5)) < 1e-12 for j in range(1, n + 1))
        ok[2] &= abs(2 * T / (w + 0.5) - H) < 1e-12 and abs(-1 / (w + 0.5) + 2 * b2 / (8 * np.pi - b2)) < 1e-12
        ok[3] &= all(abs(2 * np.cos(np.pi * aj(j, w) / (2 * T)) - 2 * np.sin((n - j) * np.pi / H)) < 1e-12 for j in range(1, n))
report("a_n = T", ok[0])
report("T - a_j = (n-j)(omega + 1/2)", ok[1])
report("2T = H(omega + 1/2) with H = 2n + B, B = -1/(omega+1/2) = -2 beta^2/(8 pi - beta^2)", ok[2])
report("fusion mass 2cos(pi a_j/2T) = 2 sin((n-j) pi/H), n = 2..6", ok[3])
w0 = 2.37
report("c_3: poles 3+3->2 at omega and 3+3->1 at 2 omega + 1/2",
       abs(aj(1, w0) - w0) < 1e-12 and abs(aj(2, w0) - (2 * w0 + 0.5)) < 1e-12)

# ---- figure ---------------------------------------------------------------
n, w = 3, w0
T = Tn(n, w)
fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.9), gridspec_kw={"width_ratios": [1.25, 1]})
ks = np.arange(1, int(np.floor(T)) + 1)
for y, f, tag in ((1.0, lambda t: t, "direct channel"), (0.0, lambda t: T - t, r"crossed, $t\to T-t$")):
    a.hlines(y, 0, T, color=style.GREY, lw=0.8)
    a.plot(f(T - ks), np.full(ks.shape, y), "|", color=style.ORANGE, ms=9, mew=1.4)
    for j in range(1, n):
        a.plot([f(aj(j, w))], [y], "o", color=style.BLUE, ms=6)
    a.text(-0.15, y, tag, ha="right", va="center", fontsize=7.5, color=style.GREY)
for j, lab in ((1, r"$3{+}3\to2$"), (2, r"$3{+}3\to1$")):
    a.text(aj(j, w), 1.18, rf"$a_{j}$" + "\n" + lab, ha="center", va="bottom", fontsize=7, color=style.BLUE)
a.annotate("", (aj(2, w), 0.78), (aj(1, w), 0.78), arrowprops=dict(arrowstyle="<->", color=style.BLUE, lw=0.8))
a.text((aj(1, w) + aj(2, w)) / 2, 0.66, r"$\omega+\frac{1}{2}$", ha="center", va="top", fontsize=7.5, color=style.BLUE)
a.plot([], [], "o", color=style.BLUE, ms=6, label=r"bound states $a_j$")
a.plot([], [], "|", color=style.ORANGE, ms=9, mew=1.4, label=r"breather tower $T-k$ (spacing 1)")
a.legend(loc="lower center", ncol=2, fontsize=7.5, bbox_to_anchor=(0.5, -0.55), columnspacing=1.0, handletextpad=0.2)
a.axvline(T, color=style.GREY, lw=0.6, ls=":")
a.set_xlim(-0.1, T + 0.3)
a.set_ylim(-0.35, 1.75)
a.set_yticks([])
a.spines["left"].set_visible(False)
a.set_xticks([0, T])
a.set_xticklabels(["0", r"$T$"])
a.set_xlabel(r"$t$ (physical strip $0<t<T$), $\omega=2.37$")
a.set_title(r"(a) poles of $S_{33}$, $c_3^{(1)}$", fontsize=9, color=style.GREY)

ww = np.linspace(0.3, 6, 300)
for aa, col in ((2, style.BLUE), (1, style.GREEN)):
    H = 2 * n - 1 / (ww + 0.5)
    b.plot(ww, 2 * np.sin(aa * np.pi / H), color=col)
    b.axhline(2 * np.sin(aa * np.pi / (2 * n)), color=col, lw=0.6, ls=":")
    for w1 in (2.37, 1.61):
        b.plot([w1], [2 * np.sin(aa * np.pi / (2 * n - 1 / (w1 + 0.5)))], "o", color=col, ms=4)
    b.text(6.0, 2 * np.sin(aa * np.pi / (2 * n)) - 0.03, rf"$M_{aa}/M_3$", color=col, ha="right", va="top", fontsize=7.5)
b.set_ylim(0.85, 2.0)
b.set_xlabel(r"$\omega=4\pi/\beta^2-1$")
b.set_ylabel(r"$2\sin(a\pi/H)$")
b.set_title(r"(b) soliton mass ratios, $H=6+B$", fontsize=9, color=style.GREY)
fig.tight_layout()
style.save(fig, "ch24_f23_cn_ladders", style.OUT_PART6)
raise SystemExit(0 if style.all_ok() else 1)
