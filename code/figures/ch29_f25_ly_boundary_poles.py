"""Figure F25 (chapter 28, sec-boundary-lee-yang, eq-ly-reflection): poles and zeros of the boundary
Lee-Yang reflection factor K_b in the physical strip 0 < Im theta < pi/2, as functions of b.

K_b = (1)(3)(4)^{-1}(1-b)^{-1}(1+b)(5-b)(5+b)^{-1}, with the blocks (x) = sinh(th/2 + i pi x/12)/sinh(th/2 - i pi x/12).
A block (x) has its pole at th = i pi x/6 and its zero at th = -i pi x/6 (mod 2 pi i).  The pole of (1+b) enters the strip
at b = -1, leaves it at b = 2 where the first bound-state energy M cos((b+1) pi/6) vanishes; (b-1) is the second one.
The factor (1-b)^{-1} has its zero at th = i(1-b) pi/6, in the strip for -2 < b < 1; at b = 0 it cancels the pole of
(1+b), so that K_0 = K_(1) has a simple pole at i pi/6.  At b = 2 the tracks 1+b and 5-b both reach i pi/2 (a triple
pole with the fixed pole of (3)), and the track b-1 reaches i pi/6 (a double pole with the fixed pole of (1)).

Writes part7-half-line/figures/ch29_f25_ly_boundary_poles.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
FACT = lambda b: [(1, 1), (3, 1), (4, -1), (1 - b, -1), (1 + b, 1), (5 - b, 1), (5 + b, -1)]   # (x, exponent)


def blk(x, th):
    return np.sinh(th / 2 + 1j * np.pi * x / 12) / np.sinh(th / 2 - 1j * np.pi * x / 12)


def K(b, th):
    out = 1.0 + 0j
    for x, e in FACT(b):
        out = out * blk(x, th) ** e
    return out


def net_orders(b):
    """net order of the pole (+) / zero (-) at theta = i pi y/6, y mod 12, from the factor list"""
    d = {}
    for x, e in FACT(b):
        for y, o in ((x % 12, e), ((-x) % 12, -e)):      # block: pole at +x, zero at -x
            d[round(y, 9)] = d.get(round(y, 9), 0) + o
    return {y: o for y, o in d.items() if o != 0}


def strip_poles(b):
    return sorted((y, o) for y, o in net_orders(b).items() if 0 < y <= 3 + 1e-9)


# ---- checks ----------------------------------------------------------------
eps = 1e-7
ok_all = True
for b in (0.0, 0.7, 1.5, 2.6, 3.4, 4.5):
    for y, o in strip_poles(b):
        th0 = 1j * np.pi * y / 6
        v = abs(K(b, th0 + eps)) * eps if o > 0 else abs(K(b, th0 + eps)) / eps
        ok_all &= (o in (1, -1)) and (1e-3 < v < 1e3)
report("every pole/zero predicted by the factor list in the strip is a simple pole/zero of K_b (b = 0, .7, 1.5, 2.6, 3.4, 4.5)", ok_all)
report("b = 0: simple poles at i pi/6 and i pi/2 only (one-point pole and double reflection)",
       [(round(y, 6), o) for y, o in strip_poles(0.0) if o > 0] == [(1.0, 1), (3.0, 1)])
th = 0.3 + 0.9j
report("K_{b=0} equals K_(1) = (1)(3)(4)^{-1}",
       abs(K(0.0, th) - blk(1, th) * blk(3, th) / blk(4, th)) < 1e-12)
bb = np.linspace(-2.5, 5.5, 3201)
y1 = 1 + bb
inside = (y1 > 0) & (y1 < 3)
report("the pole of (1+b) is in the strip exactly for -1 < b < 2", abs(bb[inside].min() + 1) < 4e-3 and abs(bb[inside].max() - 2) < 4e-3)
report("first gap M cos((b+1) pi/6): M cos(pi/6) at b = 0 and 0 at b = 2",
       abs(np.cos(np.pi / 6) - np.sqrt(3) / 2) < 1e-14 and abs(np.cos(3 * np.pi / 6)) < 1e-15)
inside2 = (bb - 1 > 0) & (bb - 1 < 3)
report("the pole of (1-b)^{-1} at i(b-1) pi/6 is in the strip exactly for 1 < b < 4", abs(bb[inside2].min() - 1) < 4e-3 and abs(bb[inside2].max() - 4) < 4e-3)
scan_ok = True
ys = np.linspace(0.02, 3.3, 2400)               # theta = i pi y/6 on the imaginary axis (a little past the strip edge y = 3)
for b in (-0.5, 0.5, 1.5, 2.5, 4.5):
    pred = [y for y, o in strip_poles(b) if o > 0]
    v = np.abs(K(b, 1j * np.pi * ys / 6 + 1e-3))
    w = 40
    peaks = [ys[i] for i in range(w, len(ys) - w) if v[i] >= v[i - 1] and v[i] >= v[i + 1]
             and v[i] > 3 * min(v[i - w], v[i + w]) and ys[i] < 3.05]
    scan_ok &= all(any(abs(f - p) < 0.05 for p in pred) for f in peaks) and all(any(abs(f - p) < 0.05 for f in peaks) for p in pred)
report("b = 0: the zero of (1-b)^{-1} cancels the pole of (1+b) at i pi/6; b = 2: triple pole at i pi/2, double pole at i pi/6",
       strip_poles(0.0)[0] == (1, 1) and dict(strip_poles(2.0)) == {1: 2, 3: 3}
       and 1e-3 < abs(K(2.0, 1j * np.pi / 2 + 1e-4)) * 1e-12 < 1e3)
report("a scan of K_b along the imaginary axis finds large values exactly at the predicted poles (b = -.5, .5, 1.5, 2.5, 4.5)", scan_ok)

# ---- figure ----------------------------------------------------------------
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.1), gridspec_kw={"width_ratios": [1.35, 1]})
b = np.linspace(-1.5, 5.0, 2001)
tracks = {  # label, y(b) in units i pi/6, colour, kind
    "1+b": (lambda b: 1 + b, style.RED), "b-1": (lambda b: b - 1, style.BLUE),
    "5-b": (lambda b: 5 - b, style.GREEN), "7-b": (lambda b: 7 - b, style.PURPLE),
}
ax.axhspan(0, 0.5, color=style.LIGHTGREY, lw=0)
ax.axhline(1 / 6, color=style.GREY, lw=1.6, alpha=.8)
ax.axhline(1 / 2, color=style.GREY, lw=1.6, alpha=.8)
ax.text(5.05, 1 / 6, r"$i\pi/6$", va="center", fontsize=7.5)
ax.text(5.05, 1 / 2, r"$i\pi/2$", va="center", fontsize=7.5)
names = {"1+b": r"$i(1+b)\pi/6$ (first bound state)", "b-1": r"$i(b-1)\pi/6$ (second)",
         "5-b": r"$i(5-b)\pi/6$", "7-b": r"$i(7-b)\pi/6$"}
for lab, (f, c) in tracks.items():
    y = f(b) / 6
    y = np.where((y > 0) & (y <= .5), y, np.nan)
    ax.plot(b, y, color=c, lw=1.8, label=names[lab])
yz = (1 - b) / 6
ax.plot(b, np.where((yz > 0) & (yz <= .5), yz, np.nan), color=style.RED, lw=1.2, ls="--", label=r"zero at $i(1-b)\pi/6$")
for bx, ls in ((0, ":"), (2, "--")):
    ax.axvline(bx, color=style.GREY, lw=.8, ls=ls)
ax.text(0.03, 0.51, r"$b=0$", fontsize=7.5)
ax.text(2.03, 0.51, r"$h_{\rm crit}$", fontsize=7.5)
ax.set_xlim(-1.5, 5.0)
ax.set_ylim(0, .54)
ax.set_xlabel(r"$b$")
ax.set_ylabel(r"pole position  ${\rm Im}\,\theta/\pi$")
ax.set_title("poles of $K_b$ in the physical strip", fontsize=8.5)
ax.legend(fontsize=6.8, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, columnspacing=1.0, handlelength=1.4)

e1 = np.where((bb > -1) & (bb < 2), np.cos((bb + 1) * np.pi / 6), np.nan)
e2 = np.where((bb > 1) & (bb < 4), np.cos((bb - 1) * np.pi / 6), np.nan)
ax2.plot(bb, e1, color=style.RED, label=r"first: $\cos((b+1)\pi/6)$")
ax2.plot(bb, e2, color=style.BLUE, label=r"second: $\cos((b-1)\pi/6)$")
ax2.axvline(2, color=style.GREY, lw=.8, ls="--")
ax2.set_xlim(-1.5, 5.0)
ax2.set_ylim(0, 1.02)
ax2.set_xlabel(r"$b$")
ax2.set_ylabel(r"$e_\beta-e_0\;(\,M\,)$")
ax2.set_title("boundary bound-state gaps", fontsize=8.5)
ax2.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=1)
fig.tight_layout()
style.save(fig, "ch29_f25_ly_boundary_poles", style.OUT_PART7)
raise SystemExit(0 if style.all_ok() else 1)
