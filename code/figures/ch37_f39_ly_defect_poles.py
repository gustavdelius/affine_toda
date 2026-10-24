"""Figure F39 (chapter 36, sec-defect-lee-yang, eq-ly-transmission, prp-ly-transmission): poles and zeros of the
Lee-Yang defect transmission factors in the physical strip 0 < Im theta < pi/2, against the defect parameter b.

TT_b(theta) = f_{2/3}(theta - i pi(3-b)/6), TTbar_b(theta) = f_{2/3}(theta + i pi(3-b)/6), f_x = (sinh theta + i sin pi x)/
(sinh theta - i sin pi x).  On the imaginary axis theta = i pi y:
  TT_b:    poles y = (5-b)/6, (7-b)/6;  zeros y = (1-b)/6, -(1+b)/6     (mod 2)
  TTbar_b: poles y = (b+1)/6, (b-1)/6;  zeros y = (b+5)/6, (b+7)/6     (mod 2)
The pole of TTbar_b at (b+1)/6 is the defect bound state |1^+> for -1 < b < 2, the one at (b-1)/6 the second bound
state for 1 < b < 2 (Bajnok-Simon, sec. 6).  The fundamental ground-state range is -3 <= b <= 2; the trivial defect 1
has TT = TTbar = 1 and no poles or zeros.

Checks: the closed forms agree with Bajnok-Simon's blocks [b+1][b-1] and [5-b][-5-b]; the predicted points are simple
poles and simple zeros; a scan along the imaginary axis of the strip finds large and small values exactly at the
predicted points (b = -2.5, ..., 2.5); the strip ranges of each line; the excited defect has TT_{b+2}, TTbar_{b+2}
(eq-defect-bootstrap) and energy m sin((b+2) pi/6) = m sin(b pi/6) + m cos u; the gap m cos u vanishes at b = 2.

Writes part8-defects/figures/ch37_f39_ly_defect_poles.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
S3 = np.sin(2 * np.pi / 3)
f23 = lambda z: (np.sinh(z) + 1j * S3) / (np.sinh(z) - 1j * S3)
TT = lambda b, z: f23(z - 1j * np.pi * (3 - b) / 6)
TTb = lambda b, z: f23(z + 1j * np.pi * (3 - b) / 6)
LINES = {  # (factor, kind): list of y(b) = Im theta / pi
    ("T", "pole"): [lambda b: (5 - b) / 6, lambda b: (7 - b) / 6],
    ("T", "zero"): [lambda b: (1 - b) / 6, lambda b: -(1 + b) / 6],
    ("Tb", "pole"): [lambda b: (b + 1) / 6, lambda b: (b - 1) / 6],
    ("Tb", "zero"): [lambda b: (b + 5) / 6, lambda b: (b + 7) / 6],
}
FUN = {"T": TT, "Tb": TTb}


def in_strip(y):
    return 1e-9 < y < 0.5 - 1e-9


# ---- checks ----------------------------------------------------------------
brk = lambda x, z: 1j * np.sinh(z / 2 + 1j * np.pi * x / 12) / np.sinh(z / 2 + 1j * np.pi * x / 12 - 1j * np.pi / 2)
pts = np.array([0.37 + 0.11j, -0.8 + 0.3j, 1.3 - 0.2j])
e = max(abs(TT(b, z) - brk(b + 1, z) * brk(b - 1, z)) + abs(TTb(b, z) - brk(5 - b, z) * brk(-5 - b, z))
        for b in (0.4, -1.3, 1.7, -2.6) for z in pts)
report(f"eq-ly-transmission = Bajnok-Simon [b+1][b-1] and [5-b][-5-b] [{e:.0e}]", e < 1e-12)

ok = True
for b in (-2.7, -1.3, 0.4, 1.6, 2.6):
    for (fac, kind), ys in LINES.items():
        for yf in ys:
            z0, eps = 1j * np.pi * yf(b), 1e-6
            val = FUN[fac](b, z0 + eps)
            if kind == "pole":         # simple pole: eps * TT(z0 + eps) tends to a finite non-zero residue
                r1, r2 = eps * val, 0.5 * eps * FUN[fac](b, z0 + 0.5 * eps)
                ok &= abs(r1) > 1e-3 and abs(r1 / r2 - 1) < 1e-4
            else:                      # simple zero
                ok &= abs(val) < 1e-4 and abs(FUN[fac](b, z0 + 2 * eps) / val - 2) < 1e-4
report("the predicted points are simple poles and simple zeros of TT_b and TTbar_b (b = -2.7, -1.3, 0.4, 1.6, 2.6)", ok)

ok = True
ys = np.linspace(0.002, 0.498, 2481)
for b in (-2.5, -1.5, -0.5, 0.5, 1.5, 2.5):
    for fac in ("T", "Tb"):
        v = np.log(np.abs(FUN[fac](b, 1j * np.pi * ys + 1e-4)))
        w = 30
        found = {"pole": [ys[i] for i in range(w, len(ys) - w) if v[i] == v[i - w:i + w].max() and v[i] > 4],
                 "zero": [ys[i] for i in range(w, len(ys) - w) if v[i] == v[i - w:i + w].min() and v[i] < -4]}
        for kind in ("pole", "zero"):
            pred = [yf(b) % 2 for yf in LINES[(fac, kind)] if in_strip(yf(b) % 2)]
            ok &= len(found[kind]) == len(pred) and all(any(abs(f - p) < 2e-3 for f in found[kind]) for p in pred)
report("a scan of TT_b and TTbar_b along i(0, pi/2) finds poles and zeros exactly at the predicted points "
       "(b = -2.5, -1.5, ..., 2.5)", ok)

bb = np.linspace(-3, 3, 6001)
rng = lambda yf: (bb[[in_strip(yf(b)) for b in bb]].min(), bb[[in_strip(yf(b)) for b in bb]].max())
r1, r2 = rng(LINES[("Tb", "pole")][0]), rng(LINES[("Tb", "pole")][1])
report(f"in -3 <= b <= 3 the TTbar pole (b+1)/6 is in the strip for {r1[0]:.3f} < b < {r1[1]:.3f}, (b-1)/6 for "
       f"{r2[0]:.3f} < b <= 3, and TT has a pole in the strip only for b > 2",
       abs(r1[0] + 1) < 2e-3 and abs(r1[1] - 2) < 2e-3 and abs(r2[0] - 1) < 2e-3
       and all(not in_strip(yf(b)) for yf in LINES[("T", "pole")] for b in bb[bb <= 2]))

ok = True
S = lambda z: f23(z)
for b in (-0.5, 0.6, 1.4):
    u = (b + 1) * np.pi / 6
    ok &= max(abs(S(z + 1j * u) * TT(b, z) - TT(b + 2, z)) + abs(S(z - 1j * u) * TTb(b, z) - TTb(b + 2, z))
              for z in pts) < 1e-12
    ok &= abs(np.sin(b * np.pi / 6) + np.cos(u) - np.sin((b + 2) * np.pi / 6)) < 1e-14
report("excited defect (eq-defect-bootstrap, u -> -u): TT_b S(theta + iu) = TT_{b+2}, TTbar_b S(theta - iu) = "
       "TTbar_{b+2}; m sin(b pi/6) + m cos u = m sin((b+2) pi/6); the gap m cos u vanishes at b = 2",
       ok and abs(np.cos(3 * np.pi / 6)) < 1e-15)

# ---- figure ----------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.0), sharey=True)
B0, B1 = -3.0, 3.0
b = np.linspace(B0, B1, 3001)
col = {("Tb", 0): style.RED, ("Tb", 1): style.BLUE}
for ax, fac, title in zip(axs, ("T", "Tb"), (r"$\mathbb{T}_b(\theta)$, from the left",
                                               r"$\bar{\mathbb{T}}_b(\theta)$, from the right")):
    ax.axvspan(2, B1, color=style.LIGHTGREY, lw=0)
    ax.axvline(2, color=style.GREY, lw=0.8, ls="--")
    for kind in ("pole", "zero"):
        for k, yf in enumerate(LINES[(fac, kind)]):
            y = np.array([yf(x) % 2 for x in b])
            y = np.where([in_strip(v) for v in y], y, np.nan)
            if kind == "zero":
                ax.plot(b, y, color=style.GREY, lw=1.2, ls=(0, (3, 2)))
                continue
            c = col.get((fac, k), style.GREY)
            if (fac, k) in col:     # bound-state pole: coloured inside the fundamental range, grey beyond b = 2
                yin = np.where(b <= 2, y, np.nan)
                ax.plot(b, np.where(b >= 2, y, np.nan), color=style.GREY, lw=1.8)
                ax.plot(b, yin, color=c, lw=2.0)
            else:
                ax.plot(b, y, color=c, lw=1.8)
    ax.set_xlim(B0, B1)
    ax.set_ylim(0, 0.5)
    ax.set_xlabel(r"$b$")
    ax.set_title(title, fontsize=8.5)
    ax.text(2.5, 0.485 if fac == "Tb" else 0.015, "beyond the\nground-state\nrange", ha="center",
            va="top" if fac == "Tb" else "bottom", fontsize=6.8)
axs[0].set_ylabel(r"${\rm Im}\,\theta/\pi$ (physical strip)")
axs[0].set_yticks([0, 1 / 6, 1 / 3, 1 / 2])
axs[0].set_yticklabels(["0", r"$\frac{1}{6}$", r"$\frac{1}{3}$", r"$\frac{1}{2}$"])
a0 = axs[0]
a0.text(-0.9, 0.36, r"zero $\frac{1-b}{6}$", fontsize=7.5, rotation=0)
a0.text(-2.95, 0.05, r"zero $-\frac{1+b}{6}$", fontsize=7.5)
a0.text(2.12, 0.25, r"pole" + "\n" + r"$\frac{5-b}{6}$", fontsize=7.5)
a1 = axs[1]
a1.text(-0.7, 0.24, r"$|1^+\rangle$", color=style.RED, fontsize=8)
a1.text(1.15, 0.035, "2nd", color=style.BLUE, fontsize=7.5)
a1.text(-2.95, 0.2, r"zero $\frac{b+5}{6}$", fontsize=7.5)
handles = [plt.Line2D([], [], color=style.RED, lw=2.0), plt.Line2D([], [], color=style.BLUE, lw=2.0),
           plt.Line2D([], [], color=style.GREY, lw=1.8), plt.Line2D([], [], color=style.GREY, lw=1.2, ls=(0, (3, 2)))]
labels = [r"pole: defect bound state $|1^+\rangle$ at $\frac{b+1}{6}$",
          r"pole: second bound state at $\frac{b-1}{6}$", "other pole", "zero"]
fig.legend(handles, labels, loc="lower center", ncol=2, fontsize=7.2, bbox_to_anchor=(0.5, -0.01),
           columnspacing=1.5, handlelength=2.0)
fig.tight_layout(rect=(0, 0.12, 1, 1))
style.save(fig, "ch37_f39_ly_defect_poles", style.OUT_PART8)
raise SystemExit(0 if style.all_ok() else 1)
