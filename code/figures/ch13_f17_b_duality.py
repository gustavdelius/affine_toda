"""Figure F17 (chapter 13, sec-blocks and sec-strong-weak): the parameter B and the block {x}.

(a) B(beta) = (1/2pi) beta^2/(1+beta^2/4pi) rises from 0 to 2 and satisfies B(4pi/beta) = 2 - B(beta),
    so beta^2 = 4 pi is the self-dual point B = 1.
(b) Poles and zeros of {x} = (x-1)(x+1)/((x-1+B)(x+1-B)) on the line Im theta = pi*y/h (y is a position
    in units of pi/h) for x = 5: the poles at x -+ 1 are fixed, the zeros at x-1+B and x+1-B move with B,
    meet at x at B = 1, and exchange places under B -> 2 - B.

Writes part4-real-coupling/figures/ch13_f17_b_duality.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
B = lambda b: b**2 / (2 * np.pi) / (1 + b**2 / (4 * np.pi))
h = 12


def blk(x, th):                      # (x) with Coxeter number h
    return np.sinh(th / 2 + 1j * np.pi * x / (2 * h)) / np.sinh(th / 2 - 1j * np.pi * x / (2 * h))


def cdd(x, Bv, th):                  # {x}
    return blk(x - 1, th) * blk(x + 1, th) / (blk(x - 1 + Bv, th) * blk(x + 1 - Bv, th))


# ---- checks ---------------------------------------------------------------
bs = np.linspace(0.05, 30, 400)
report("B(4 pi/beta) = 2 - B(beta)", np.allclose(B(4 * np.pi / bs), 2 - B(bs)))
report("B monotonic from 0 to 2", np.all(np.diff(B(bs)) > 0) and B(1e-6) < 1e-10 and abs(B(1e6) - 2) < 1e-6)
report("B = 1 at beta^2 = 4 pi", abs(B(np.sqrt(4 * np.pi)) - 1) < 1e-12)
th = 0.3 + 0.7j
report("{x} invariant under B -> 2 - B", abs(cdd(5, 0.3, th) - cdd(5, 1.7, th)) < 1e-12)
report("{x} = 1 at B = 0", abs(cdd(5, 0.0, th) - 1) < 1e-12)
report("{x} is a phase for real theta", abs(abs(cdd(5, 0.4, 0.8)) - 1) < 1e-12)
eps = 1e-7
for Bv in (0.4, 1.3):                 # poles at x -+ 1 fixed, zeros at x-1+B, x+1-B
    x = 5
    pole_ok = all(abs(cdd(x, Bv, 1j * np.pi * (x + s) / h + eps)) > 1e5 for s in (-1, 1))
    zero_ok = all(abs(cdd(x, Bv, 1j * np.pi * y / h + eps)) < 1e-4 for y in (x - 1 + Bv, x + 1 - Bv))
    report(f"B={Bv}: simple poles at x+-1, zeros at x-1+B and x+1-B (h={h}, x={x})", pole_ok and zero_ok)

# ---- figure ---------------------------------------------------------------
fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.9), gridspec_kw={"width_ratios": [1.1, 1]})
bb = np.linspace(0, 12, 600)
a.plot(bb, B(bb), color=style.BLUE)
a.axhline(2, color=style.GREY, lw=0.6, ls=":")
a.plot([np.sqrt(4 * np.pi)], [1], "o", color=style.RED)
a.annotate(r"self-dual: $\beta^2=4\pi$, $B=1$", (np.sqrt(4 * np.pi), 1), (4.6, 0.45), fontsize=7.5,
           color=style.RED, arrowprops=dict(arrowstyle="-", color=style.RED, lw=0.7))
a.plot([1.2, 4 * np.pi / 1.2], [B(1.2), B(4 * np.pi / 1.2)], "s", color=style.ORANGE, ms=4)
a.text(0.2, 1.72, r"$\beta\leftrightarrow4\pi/\beta$: $B\leftrightarrow2-B$", fontsize=7.5, color=style.ORANGE)
a.set_xlabel(r"$\beta$")
a.set_ylabel(r"$B(\beta)$")
a.set_ylim(0, 2.15)
a.set_title("(a) coupling parameter", fontsize=9, color=style.GREY)

x = 5
bv = np.linspace(0, 2, 200)
b.axvline(0, color=style.GREY, lw=0)
for s in (-1, 1):
    b.plot(bv, np.full_like(bv, x + s), color=style.GREEN, lw=2.2)
b.plot(bv, x - 1 + bv, color=style.RED, lw=1.5, ls="--")
b.plot(bv, x + 1 - bv, color=style.RED, lw=1.5, ls="--")
b.axvline(1, color=style.GREY, lw=0.6, ls=":")
b.text(1.97, x + 1.07, r"poles $x\pm1$ (fixed)", color=style.GREEN, fontsize=7.5, ha="right", va="bottom")
b.text(1.05, x + 0.1, r"zeros meet at $B=1$", color=style.RED, fontsize=7.5)
b.text(0.12, x + 0.55, r"zero $x+1-B$", color=style.RED, fontsize=7.5, rotation=-34, rotation_mode="anchor")
b.text(0.12, x - 0.7, r"zero $x-1+B$", color=style.RED, fontsize=7.5, rotation=34, rotation_mode="anchor")
b.set_xlabel(r"$B$")
b.set_ylabel(r"position $y$ on $\mathrm{Im}\,\theta=\pi y/h$")
b.set_yticks([x - 1, x, x + 1])
b.set_yticklabels([r"$x-1$", r"$x$", r"$x+1$"])
b.set_xlim(0, 2)
b.set_ylim(x - 1.4, x + 1.5)
b.set_title(r"(b) poles and zeros of $\{x\}$", fontsize=9, color=style.GREY)
fig.tight_layout()
style.save(fig, "ch13_f17_b_duality", style.OUT_PART4)
raise SystemExit(0 if style.all_ok() else 1)
