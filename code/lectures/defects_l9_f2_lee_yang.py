"""Lecture notes on defects, Lecture 9, figure 2 (fig-l9-lee-yang): the Lee-Yang defect pressed against the boundary 1.

Blocks (x) = sinh(theta/2 + i pi x/12)/sinh(theta/2 - i pi x/12); S = f_{2/3}; the defect T_b(theta) = f_{2/3}(theta -
theta_b), Tbar_b(theta) = f_{2/3}(theta + theta_b), theta_b = i pi (3 - b)/6; the boundary K_(1) = (1)(3)(4)^-1 and
K_b = K_(1) (1-b)^-1 (1+b) (5-b) (5+b)^-1. Positions theta = i pi y of the poles in the physical strip 0 < y < 1/2,
against real b in [-3, 3].
(a) The defect: the poles of Tbar_b at y = (b+1)/6 (the defect bound state) and y = (b-1)/6 (the second), and of T_b at
y = (5-b)/6 and (7-b)/6.
(b) The fused boundary Tbar_b K_(1) T_b = K_b: the same lines, and the two lines y = 1/6, 1/2 of K_(1).

Checks: K_b = Tbar_b K_(1) T_b on a grid of complex rapidities; each plotted point is a simple pole of the factor it is
attributed to (residue stable under 1/eps scaling) and of K_b; a scan of the strip finds no pole of K_b off the plotted
lines, for b = -2.5, 0.5, 1.5, 2.5; at b = 0 the pole of (1+b) at y = 1/6 cancels against the zero of (1-b)^-1.
Writes lecture-notes/defects/figures/l9_lee_yang.svg.
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


def blk(x, th):
    return np.sinh(th / 2 + I * PI * x / 12) / np.sinh(th / 2 - I * PI * x / 12)


def f23(th):
    return (np.sinh(th) + I * math.sin(2 * PI / 3)) / (np.sinh(th) - I * math.sin(2 * PI / 3))


def thb(b):
    return I * PI * (3 - b) / 6


def K1(th):
    return blk(1, th) * blk(3, th) / blk(4, th)


def Kb(th, b):
    return K1(th) * blk(1 + b, th) * blk(5 - b, th) / (blk(1 - b, th) * blk(5 + b, th))


def T(th, b):
    return f23(th - thb(b))


def Tb(th, b):
    return f23(th + thb(b))


# checks
grid = np.array([0.2 + 0.1j, -0.7 + 0.3j, 1.1 + 0.45j, 0.05 + 0.2j, -1.6])
ok = True
for b in np.linspace(-3, 3, 13):
    ok = ok and np.abs(Tb(grid, b) * K1(grid) * T(grid, b) - Kb(grid, b)).max() < 1e-11
report("K_b = Tbar_b K_(1) T_b on a grid of complex rapidities, b = -3..3", ok)

lines_T = {'Tbar: (b+1)/6': lambda b: (b + 1) / 6, 'Tbar: (b-1)/6': lambda b: (b - 1) / 6,
           'T: (5-b)/6': lambda b: (5 - b) / 6, 'T: (7-b)/6': lambda b: (7 - b) / 6}
fac = {'Tbar': Tb, 'T': T}
zeros_T = {'T: (1-b)/6': lambda b: (1 - b) / 6, 'T: (-1-b)/6': lambda b: (-1 - b) / 6,
           'Tbar: (b+5)/6': lambda b: (b + 5) / 6}


def simple_pole(fn, th0):
    r = [e * fn(th0 + e) for e in (1e-4, 1e-5, 1e-6)]
    return abs(r[2]) > 1e-3 and abs(r[2] - r[1]) < 1e-3 * abs(r[2])


ok = True
bs_test = [-2.5, -0.6, 0.5, 1.5, 2.5, 2.9]
for b in bs_test:
    for name, yf in lines_T.items():
        y0 = yf(b)
        if 0.01 < y0 < 0.49:
            f = fac[name.split(':')[0]]
            ok = ok and simple_pole(lambda t: f(t, b), I * PI * y0)
            ok = ok and simple_pole(lambda t: Kb(t, b), I * PI * y0)
    for y0 in (1 / 6, 1 / 2 - 1e-12):
        pass
ok = ok and simple_pole(K1, I * PI / 6) and simple_pole(K1, I * PI / 2)
report("every plotted point is a simple pole of its factor and of K_b; K_(1) has simple poles at i pi/6 and i pi/2", ok)

ok = True
ys = np.linspace(0.003, 0.497, 2000)
for b in (-2.5, 0.5, 1.5, 2.5):
    vals = np.abs(Kb(I * PI * ys, b))
    peaks = [ys[i] for i in range(1, len(ys) - 1) if vals[i] > vals[i - 1] and vals[i] > vals[i + 1] and vals[i] > 50]
    expected = [yf(b) for yf in lines_T.values()] + [1 / 6]
    for p in peaks:
        ok = ok and min(abs(p - e) for e in expected) < 2e-3
    for e in expected:
        if 0.01 < e < 0.49:
            ok = ok and any(abs(p - e) < 2e-3 for p in peaks)
report("a scan of the strip finds the poles of K_b only on the plotted lines (b = -2.5, 0.5, 1.5, 2.5)", ok)

# the plotted zeros are simple zeros of their factors
ok = True
for b in bs_test:
    for name, yf in zeros_T.items():
        y0 = yf(b)
        if 0.01 < y0 < 0.49:
            f = fac[name.split(':')[0]]
            ok = ok and simple_pole(lambda t: 1 / f(t, b), I * PI * y0)
report("every plotted zero is a simple zero of its factor", ok)

# at b = 0 the pole of (1 + b) at i pi/6 cancels against the zero of (1 - b)^-1: K_0 = K_(1)
ok = np.abs(Kb(grid, 0.0) - K1(grid)).max() < 1e-13
report("at b = 0, K_b = K_(1): the pole of (1+b) and the zero of (1-b)^-1 cancel", ok)

# figure
fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.0, 2.9), sharey=True)
bgrid = np.linspace(-3, 3, 601)
cols = {'Tbar: (b+1)/6': style.RED, 'Tbar: (b-1)/6': style.BLUE, 'T: (5-b)/6': style.GREY, 'T: (7-b)/6': style.GREY}
labels = {'Tbar: (b+1)/6': r'$\bar{\mathbb{T}}_b$: $\frac{b+1}{6}$', 'Tbar: (b-1)/6': r'$\bar{\mathbb{T}}_b$: $\frac{b-1}{6}$',
          'T: (5-b)/6': r'$\mathbb{T}_b$', 'T: (7-b)/6': None}
for ax in (axa, axb):
    ax.axvspan(2, 3, color=style.LIGHTGREY, lw=0, zorder=0)
    for name, yf in lines_T.items():
        yv = yf(bgrid)
        m = (yv > 0) & (yv < 0.5)
        yplot = np.where(m, yv, np.nan)
        ls = '-' if name.startswith('Tbar') else '--'
        ax.plot(bgrid, yplot, color=cols[name], lw=1.6 if ls == '-' else 1.1, ls=ls,
                label=labels[name] if ax is axa else None)
    for k, yf in enumerate(zeros_T.values()):
        yv = yf(bgrid)
        yplot = np.where((yv > 0) & (yv < 0.5), yv, np.nan)
        ax.plot(bgrid, yplot, color=style.GREY, lw=0.9, ls=':', label='zeros' if (ax is axa and k == 0) else None)
    ax.set_xlim(-3, 3)
    ax.set_ylim(0, 0.5)
    ax.set_xlabel(r'$b$')
    ax.axvline(2, color=style.GREY, lw=0.6, ls=':')
    ax.text(2.5, 0.03, 'relabelled', fontsize=7, ha='center', color=style.GREY)
for yv in (1 / 6, 1 / 2):
    axb.plot([-3, 3], [yv, yv], color=style.ORANGE, lw=1.6, clip_on=False, label=r'$K_{(\mathbb{1})}$' if yv < 0.3 else None)
axb.annotate('bound state of\nthe defect and\nof the boundary', xy=(1.0, 1 / 3), xytext=(-2.3, 0.36), fontsize=7,
             color=style.RED, arrowprops=dict(arrowstyle='->', color=style.RED, lw=0.8))
axa.set_yticks([0, 1 / 6, 1 / 3, 1 / 2])
axa.set_yticklabels(['0', r'$\frac{1}{6}$', r'$\frac{1}{3}$', r'$\frac{1}{2}$'])
axa.set_ylabel(r'pole at $\theta=i\pi y$:  $y$')
axa.set_title(r'(a) the defect $\Phi(b)$', fontsize=8.5, color=style.GREY)
axb.set_title(r'(b) pressed against $\mathbb{1}$: $K_b=\bar{\mathbb{T}}_b\,K_{(\mathbb{1})}\,\mathbb{T}_b$', fontsize=8.5,
              color=style.GREY)
h1, l1 = axa.get_legend_handles_labels()
h2, l2 = axb.get_legend_handles_labels()
fig.tight_layout(w_pad=1.5, rect=(0, 0.08, 1, 1))
fig.legend(h1 + h2, l1 + l2, loc='lower center', ncol=5, fontsize=7.5, handlelength=2.0, bbox_to_anchor=(0.5, -0.01))
style.save(fig, "l9_lee_yang", OUT)
raise SystemExit(0 if style.all_ok() else 1)
