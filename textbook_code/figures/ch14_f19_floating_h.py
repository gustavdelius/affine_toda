"""Figure F19 (chapter 14, sec-floating-coxeter, tbl-floating-h): the floating Coxeter number.

(a) H(B) = h + (h_dual - h) B/2 for the four non-simply-laced families (eq-floating-h), with the
    values of tbl-floating-h: c_n (n=3): 2n+B; b_n (n=3): 2n-B/2; g_2: 6+3B; f_4: 12+3B.
(b) c_n^(1), n = 2: the mass ratio m_1/m_2 = sin(pi/H)/sin(2pi/H) = 1/(2 cos(pi/H)) as B runs from
    0 (classical, H = 2n) to 2 (H = 2n+2, the masses of d_{n+1}^(2)); dots mark the endpoints.

Writes part4-real-coupling/figures/ch14_f19_floating_h.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
fam = {  # name: (h, h_dual, H(B) from tbl-floating-h)
    r"$c_3^{(1)}$": (6, 8, lambda B: 6 + B),
    r"$b_3^{(1)}$": (6, 5, lambda B: 6 - B / 2),
    r"$g_2^{(1)}$": (6, 12, lambda B: 6 + 3 * B),
    r"$f_4^{(1)}$": (12, 18, lambda B: 12 + 3 * B),
}
Bs = np.linspace(0, 2, 201)

# ---- checks ---------------------------------------------------------------
for name, (h, hd, H) in fam.items():
    report(f"{name}: H(0)=h={h}, H(2)=h_dual={hd}", abs(H(0) - h) < 1e-12 and abs(H(2) - hd) < 1e-12)
    report(f"{name}: H = h + (h_dual - h) B/2 for all B", np.allclose(H(Bs), h + (hd - h) * Bs / 2))
n = 2
H = lambda B: 2 * n + B
ratio = lambda B: np.sin(np.pi / H(B)) / np.sin(2 * np.pi / H(B))
report("c_2: classical m_1/m_2 = sin(pi/4)/sin(pi/2) at B=0", abs(ratio(0) - np.sin(np.pi / 4)) < 1e-12)
report("c_2: m_1/m_2 at B=2 equals sin(pi/6)/sin(pi/3) (d_3^(2) masses)", abs(ratio(2) - np.sin(np.pi / 6) / np.sin(np.pi / 3)) < 1e-12)
report("c_2: m_1/m_2 = 1/(2cos(pi/H))", np.allclose(ratio(Bs), 1 / (2 * np.cos(np.pi / H(Bs)))))
report("c_2: ratio decreases monotonically with B", np.all(np.diff(ratio(Bs)) < 0))
# one-loop: d/dB ratio at 0 = (pi/(2 n H^2... ) check numerically against the exact expansion
d = (ratio(1e-6) - ratio(0)) / 1e-6
report("c_2: slope at B=0 is -pi sin(pi/4) /(2 (2n)^2 cos^2(pi/4)) (one-loop drift)",
       abs(d - (-np.pi * np.sin(np.pi / 4) / (2 * (2 * n) ** 2 * np.cos(np.pi / 4) ** 2))) < 1e-4)

# ---- figure ---------------------------------------------------------------
fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.9))
cols = [style.BLUE, style.ORANGE, style.GREEN, style.PURPLE]
for (name, (h, hd, Hf)), c in zip(fam.items(), cols):
    a.plot(Bs, Hf(Bs), color=c)
    a.plot([0, 2], [h, hd], "o", color=c, ms=3.5)
    a.text(2.04, hd, name, color=c, fontsize=8, va="center")
a.set_xlim(0, 2.45)
a.set_xlabel(r"$B$")
a.set_ylabel(r"$H$")
a.set_xticks([0, 1, 2])
a.set_title("(a) floating Coxeter number", fontsize=9, color=style.GREY)
b.plot(Bs, ratio(Bs), color=style.BLUE)
b.plot([0, 2], [ratio(0), ratio(2)], "o", color=style.BLUE, ms=4)
b.text(0.12, ratio(0) + 0.002, r"$c_2^{(1)}$ classical, $H=4$", fontsize=7.5, va="center")
b.text(1.9, ratio(2) - 0.004, r"$d_3^{(2)}$ classical, $H=6$", fontsize=7.5, ha="right", va="top")
b.set_xlabel(r"$B$")
b.set_ylabel(r"$m_1/m_2$")
b.set_xticks([0, 1, 2])
b.set_title(r"(b) mass drift in $c_2^{(1)}$", fontsize=9, color=style.GREY)
fig.tight_layout()
style.save(fig, "ch14_f19_floating_h", style.OUT_PART4)
raise SystemExit(0 if style.all_ok() else 1)
