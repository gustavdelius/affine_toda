"""Figure F18 (chapter 13, sec-dorey-formula, exr-e8-poles): poles of the e8 S-matrix S_11.

S_11 = {1}{11}{19}{29}, h = 30.  Poles in the physical strip 0 < Im theta < pi sit at
theta = i pi y/30 for y = 2, 10, 12, 18, 20, 28.  The s-channel poles y = 20, 12, 2 are the
fusings 1+1 -> 1, 2, 3 (u = 2pi/3, 2pi/5, pi/15) with m_2 = 2 m_1 cos(pi/5), m_3 = 2 m_1 cos(pi/30);
the others y = 10, 18, 28 = 30 - (20, 12, 2) are their crossed images.  The figure shows the poles
and zeros (B = 0.5) and the residue sign.

Writes part4-real-coupling/figures/ch13_f18_e8_poles.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
h = 30
blk = lambda x, th: np.sinh(th / 2 + 1j * np.pi * x / (2 * h)) / np.sinh(th / 2 - 1j * np.pi * x / (2 * h))
cdd = lambda x, Bv, th: blk(x - 1, th) * blk(x + 1, th) / (blk(x - 1 + Bv, th) * blk(x + 1 - Bv, th))
S11 = lambda Bv, th: cdd(1, Bv, th) * cdd(11, Bv, th) * cdd(19, Bv, th) * cdd(29, Bv, th)

# ---- checks ---------------------------------------------------------------
Bv = 0.5
eps = 1e-7
poles = [2, 10, 12, 18, 20, 28]
res = {}
for y in poles:
    th0 = 1j * np.pi * y / h
    res[y] = eps * S11(Bv, th0 + eps)
report("S_11 has simple poles at y = 2,10,12,18,20,28 (nonzero finite residue)",
       all(1e-3 < abs(r) < 1e3 for r in res.values()))
others = [y for y in range(1, 30) if y not in poles]
report("no pole at the other integer positions y = 1..29",
       all(abs(eps * S11(Bv, 1j * np.pi * y / h + eps)) < 1e-4 for y in others))
report("residue/i is positive at the s-channel poles y = 20, 12, 2",
       all((res[y] / 1j).real > 0 and abs(res[y].real) < 1e-6 for y in (20, 12, 2)))
report("residue/i is negative at the crossed poles y = 10, 18, 28",
       all((res[y] / 1j).real < 0 for y in (10, 18, 28)))
m1 = 1.0
for y, ratio, name in ((20, 1.0, "m_1"), (12, 2 * np.cos(np.pi / 5), "m_2"), (2, 2 * np.cos(np.pi / 30), "m_3")):
    u = np.pi * y / h
    mc = np.sqrt(2 * m1**2 * (1 + np.cos(u)))
    report(f"fusing angle {y}pi/30 gives {name}/m_1 = {ratio:.4f}", abs(mc - ratio) < 1e-12)
t = 0.3 + 0.9j
report("S_11 unitarity S(t) S(-t) = 1", abs(S11(Bv, t) * S11(Bv, -t) - 1) < 1e-10)
report("S_11 crossing S(i pi - t) = S(t)", abs(S11(Bv, 1j * np.pi - t) - S11(Bv, t)) < 1e-10)

# ---- figure ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.0, 2.6))
ax.axhspan(0, np.pi, color=style.LIGHTGREY, lw=0)
ax.axhline(0, color=style.BLUE, lw=2.2)
names = {20: r"$1+1\to1$", 12: r"$1+1\to2$", 2: r"$1+1\to3$"}
for y in poles:
    s_chan = y in names
    ax.plot([y], [0], "X", color=style.GREEN if s_chan else style.RED, ms=9, zorder=5)
    if s_chan:
        ax.text(y, 0.42, names[y], ha="center", fontsize=8, color=style.GREEN)
    else:
        ax.text(y, -0.32, rf"crossed, $y={y}$", ha="center", fontsize=7.5, color=style.RED)
zeros = []
for x in (1, 11, 19, 29):
    zeros += [x - 1 + Bv, x + 1 - Bv]
zeros = [z for z in zeros if 0 < z < 30]
ax.plot(zeros, [0] * len(zeros), "o", mfc="none", mec=style.GREY, ms=6)
ax.set_xlim(-1, 31)
ax.set_ylim(-0.5, 1.2)
ax.set_xticks(range(0, 31, 5))
ax.set_xlabel(r"position $y$ on the imaginary axis, $\theta=i\pi y/30$")
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.text(30.7, 1.05, r"green: bound state, red: crossed;  circles: zeros at $B=0.5$", ha="right", fontsize=7.5)
ax.set_title(r"$e_8^{(1)}$: poles of $S_{11}=\{1\}\{11\}\{19\}\{29\}$ in the physical strip", fontsize=9, color=style.GREY)
fig.tight_layout()
style.save(fig, "ch13_f18_e8_poles", style.OUT_PART4)
raise SystemExit(0 if style.all_ok() else 1)
