"""Lecture notes on defects, Lecture 9, figure 1 (fig-l9-fold-fuse): folding the line and fusing a defect.

(a) A defect at x = 0 (purple) with Phi on x < 0 and Psi on x > 0. A particle of rapidity theta = 0.6 arrives from the
left and is transmitted, amplitude T(theta); a particle of rapidity -theta' (theta' = 0.4) arrives from the right and
is transmitted, amplitude Tbar(theta').
(b) Folded: x -> -x on the right, Psi(x) -> Psi(-x). Both fields live on x <= 0, the defect is a boundary, and the
transmitted pieces of the world lines (dashed, the copy V') come back from it: the doubled reflection matrix is
off-diagonal, K = [[0, Tbar], [T, 0]].
(c) Fusing: the defect at x = -a in front of a boundary at x = 0 with reflection factor K. A particle crosses the
defect, is reflected, and crosses back, with amplitude Tbar(theta) K(theta) T(theta), whatever a is; as a -> 0 the
composite is a boundary with K_f = Tbar K T.

Checks: the mirror map sends the transmitted segment, of velocity tanh theta on x > 0, to one of velocity -tanh theta
on x < 0, which starts where the incoming segment ends; in (c) the reflected segment has velocity -tanh theta and the
segments join at the defect and at the boundary; the amplitude along the path of (c), with the Lee-Yang factors at
b = 1.37 and theta = 0.6, is K_b(theta) for a = 0.5, 1, 2.
Writes lecture-notes/defects/figures/l9_fold_fuse.svg.
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
I = 1j
PI = math.pi
TH, THP = 0.6, 0.4
v, vp = math.tanh(TH), math.tanh(THP)


def blk(x, th):
    return np.sinh(th / 2 + I * PI * x / 12) / np.sinh(th / 2 - I * PI * x / 12)


def f23(th):
    return (np.sinh(th) + I * math.sin(2 * PI / 3)) / (np.sinh(th) - I * math.sin(2 * PI / 3))


def K1(th):
    return blk(1, th) * blk(3, th) / blk(4, th)


def Kb(th, b):
    return K1(th) * blk(1 + b, th) * blk(5 - b, th) / (blk(1 - b, th) * blk(5 + b, th))


# world lines of (a): particle 1 from (-3, -3/v... ) crossing x = 0 at t = -0.6; particle 2 crossing at t = 0.8
t1, t2 = -0.6, 0.8
seg1_in = np.array([[-2.6, t1 - 2.6 / v], [0, t1]])           # x < 0, velocity v
seg1_out = np.array([[0, t1], [2.6, t1 + 2.6 / v]])           # x > 0, velocity v
seg2_in = np.array([[2.6, t2 - 2.6 / vp], [0, t2]])           # x > 0, velocity -vp
seg2_out = np.array([[0, t2], [-2.6, t2 + 2.6 / vp]])         # x < 0, velocity -vp


def mirror(seg):
    return seg * np.array([-1, 1])


def vel(seg):
    return (seg[1, 0] - seg[0, 0]) / (seg[1, 1] - seg[0, 1])


ok = abs(vel(mirror(seg1_out)) + v) < 1e-12 and np.allclose(mirror(seg1_out)[0], seg1_in[1])
ok = ok and abs(vel(mirror(seg2_in)) - vp) < 1e-12 and mirror(seg2_in)[1, 0] == 0
report("the mirror map sends the transmitted segment (velocity tanh theta, x > 0) to velocity -tanh theta on x < 0, "
       "joined to the incoming segment at the boundary", ok)

# (c): defect at x = -a, boundary at x = 0
ok = True
for a in (0.5, 1.0, 2.0):
    tA = -1.0  # crossing the defect
    tB = tA + a / v  # reaching the boundary
    tC = tB + a / v  # crossing back
    seg_ref = np.array([[0, tB], [-a, tC]])
    ok = ok and abs(vel(seg_ref) + v) < 1e-12
    amp = f23(TH + I * PI * (3 - 1.37) / 6) * Kb(TH, 0) * f23(TH - I * PI * (3 - 1.37) / 6)  # Tbar K_(1) T
    ok = ok and abs(amp - Kb(TH, 1.37)) < 1e-12
report("in (c) the reflected segment has velocity -tanh theta; the amplitude along the path, Tbar K_(1) T for "
       "Lee-Yang at b = 1.37, is K_b(theta) for a = 0.5, 1, 2", ok)

# figure
fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.6), gridspec_kw=dict(width_ratios=[2.0, 1.25, 1.25]))
for ax in axs:
    ax.set_xticks([])
    ax.set_yticks([])
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
    ax.set_ylim(-4.6, 4.8)
ax = axs[0]
ax.set_xlim(-2.8, 2.8)
ax.axvline(0, color=style.PURPLE, lw=2.4)
ax.plot(*seg1_in.T, color=style.GREEN, lw=1.6)
ax.plot(*seg1_out.T, color=style.GREEN, lw=1.6)
ax.plot(*seg2_in.T, color=style.ORANGE, lw=1.6)
ax.plot(*seg2_out.T, color=style.ORANGE, lw=1.6)
ax.text(0.15, t1 - 0.55, r'$\mathbb{T}(\theta)$', fontsize=8.5, color=style.GREEN)
ax.text(-1.1, t2 + 0.15, r'$\bar{\mathbb{T}}(\theta^\prime)$', fontsize=8.5, color=style.ORANGE)
ax.text(-1.6, 4.2, r'$\Phi$', fontsize=10, ha='center')
ax.text(1.6, 4.2, r'$\Psi$', fontsize=10, ha='center')
ax.annotate('', xy=(2.7, -4.3), xytext=(-2.7, -4.3), arrowprops=dict(arrowstyle='->', color=style.GREY, lw=0.8))
ax.text(2.55, -4.15, r'$x$', fontsize=8.5)
ax.set_title('(a) a defect', fontsize=8.5, color=style.GREY)
ax = axs[1]
ax.set_xlim(-2.8, 0.6)
ax.axvspan(0, 0.6, color=style.LIGHTGREY, lw=0)
ax.axvline(0, color=style.PURPLE, lw=2.4)
ax.plot(*seg1_in.T, color=style.GREEN, lw=1.6)
ax.plot(*mirror(seg1_out).T, color=style.GREEN, lw=1.6, ls='--')
ax.plot(*mirror(seg2_in).T, color=style.ORANGE, lw=1.6, ls='--')
ax.plot(*seg2_out.T, color=style.ORANGE, lw=1.6)
ax.text(-2.75, 0.2, r'$\Phi(x),\ \Psi(-x)$', fontsize=8.5)
ax.set_title(r'(b) folded: $K=\left(\genfrac{}{}{0}{}{0\ \ \bar{\mathbb{T}}}{\mathbb{T}\ \ 0}\right)$', fontsize=8.5,
             color=style.GREY)
ax = axs[2]
ax.set_xlim(-2.8, 0.6)
ax.axvspan(0, 0.6, color=style.LIGHTGREY, lw=0)
ax.axvline(0, color=style.GREY, lw=2.4)
a = 1.0
ax.axvline(-a, color=style.PURPLE, lw=2.4)
tA = -2.2
tB, tC = tA + a / v, tA + 2 * a / v
ax.plot([-a - 1.7, -a], [tA - 1.7 / v, tA], color=style.GREEN, lw=1.6)
ax.plot([-a, 0], [tA, tB], color=style.GREEN, lw=1.6)
ax.plot([0, -a], [tB, tC], color=style.GREEN, lw=1.6)
ax.plot([-a, -a - 1.7], [tC, tC + 1.7 / v], color=style.GREEN, lw=1.6)
ax.text(-a + 0.06, tA - 0.45, r'$\mathbb{T}$', fontsize=8.5, color=style.GREEN)
ax.text(0.06, tB - 0.2, r'$K$', fontsize=8.5, color=style.GREY)
ax.text(-a + 0.06, tC + 0.15, r'$\bar{\mathbb{T}}$', fontsize=8.5, color=style.GREEN)
ax.annotate('', xy=(-0.05, -4.0), xytext=(-a + 0.05, -4.0), arrowprops=dict(arrowstyle='<->', color=style.GREY, lw=0.8))
ax.text(-a / 2, -4.5, r'$a\to0$', fontsize=8.5, ha='center')
ax.set_title(r'(c) fused: $K_{\rm f}=\bar{\mathbb{T}}\,K\,\mathbb{T}$', fontsize=8.5, color=style.GREY)
fig.tight_layout(w_pad=0.8)
style.save(fig, "l9_fold_fuse", OUT)
raise SystemExit(0 if style.all_ok() else 1)
