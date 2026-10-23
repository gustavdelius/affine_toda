"""Figure F2 (introduction): real versus imaginary coupling.

Panel (a): the a_1^(1) potentials in units m_0^2/beta^2 (= 1), as functions of
beta*phi: V_real = cosh(beta phi) - 1 (one vacuum) and Re V_imag = 1 - cos(beta phi)
(vacua at 2 pi k / beta; the kink runs between neighbours).
Panel (b): a_2^(1), Re V_imag = sum_j (1 - cos(beta alpha_j . phi)) in the plane of
the two fields: the degenerate vacua form a lattice, the solitons run between them.

Writes figures/intro_f2_potentials.svg.   Run: python3 intro_f2_potentials.py
"""
import numpy as np

import style
from style import plt, report

style.setup()

# ---- checks ---------------------------------------------------------------
a1, a2 = np.array([np.sqrt(2), 0.0]), np.array([-1 / np.sqrt(2), np.sqrt(1.5)])
alpha = [a1, a2, -(a1 + a2)]                       # a_2^(1): alpha_0 = -(alpha_1 + alpha_2)
report("a_2^(1): sum of the affine simple roots is zero", np.allclose(sum(alpha), 0))
report("a_2^(1): all simple roots have squared length 2",
       all(abs(a @ a - 2) < 1e-12 for a in alpha))

A = np.array([a1, a2])
def vacuum(k):                                     # alpha_j . phi in 2 pi Z for all j
    return 2 * np.pi * np.linalg.solve(A, np.array(k, float))
def re_v(p):
    return sum(1 - np.cos(a @ p) for a in alpha)

ks = [(i, j) for i in range(-2, 3) for j in range(-2, 3)]
report("Re V_imag vanishes at every vacuum 2 pi A^{-1} k",
       all(abs(re_v(vacuum(k))) < 1e-12 for k in ks))
rng = np.random.default_rng(0)
pts = rng.uniform(-6, 6, (500, 2))
vac = np.array([vacuum(k) for k in ks])
away = [p for p in pts if np.min(np.linalg.norm(vac - p, axis=1)) > 0.3]
report("Re V_imag > 0 away from the vacua", all(re_v(p) > 1e-6 for p in away))
report("V_real >= 0 with its only zero at phi = 0",
       np.all(np.cosh(np.linspace(-5, 5, 2001)) - 1 >= 0))

# ---- figure ---------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.2),
                               gridspec_kw={"width_ratios": [1.25, 1]})
x = np.linspace(-7.8, 7.8, 900)
ax1.plot(x, np.cosh(x) - 1, color=style.BLUE, label=r"real coupling: $\cosh\beta_{\rm SG}\phi-1$")
ax1.plot(x, 1 - np.cos(x), color=style.ORANGE, label=r"imaginary coupling: $1-\cos\beta_{\rm SG}\phi$")
ax1.plot([0], [0], "o", color=style.BLUE, ms=6)
for k in (-1, 0, 1):
    ax1.plot([2 * np.pi * k], [0], "o", color=style.ORANGE, ms=6)
ax1.annotate("", xy=(2 * np.pi - 0.15, 0.32), xytext=(0.15, 0.32),
             arrowprops=dict(arrowstyle="->", color=style.GREEN, lw=1.6))
ax1.text(np.pi, 0.5, "kink", ha="center", color=style.GREEN)
ax1.set_ylim(-0.12, 2.9)
ax1.set_xlim(-7.8, 7.8)
ax1.set_xticks([-2 * np.pi, 0, 2 * np.pi])
ax1.set_xticklabels([r"$-2\pi$", r"$0$", r"$2\pi$"])
ax1.set_xlabel(r"$\beta_{\rm SG}\phi$")
ax1.set_ylabel(r"$\beta_{\rm SG}^2V/m_0^2$")
ax1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=7.5)
ax1.set_title(r"(a)  $a_1^{(1)}$: one vacuum versus infinitely many", fontsize=9, color=style.GREY)

g = np.linspace(-9.5, 9.5, 400)
X, Y = np.meshgrid(g, g)
Z = sum(1 - np.cos(a[0] * X + a[1] * Y) for a in alpha)
cs = ax2.contourf(X, Y, Z, levels=np.linspace(0, 4.5, 10), cmap="viridis_r", alpha=0.9)
for c in cs.collections:
    c.set_rasterized(True)   # keep the SVG small; matplotlib 3.6 ignores rasterized= in contourf
ax2.plot(vac[:, 0], vac[:, 1], "o", color="white", mec=style.GREY, ms=4.5)
ax2.set_xlim(-9.5, 9.5)
ax2.set_ylim(-9.5, 9.5)
ax2.set_aspect("equal")
ax2.set_xlabel(r"$\beta\phi_1$")
ax2.set_ylabel(r"$\beta\phi_2$")
ax2.set_title(r"(b)  $a_2^{(1)}$: $\mathrm{Re}\,V_{\rm imag}$, vacua as dots", fontsize=10,
              color=style.GREY)
cb = fig.colorbar(cs, ax=ax2, fraction=0.046, pad=0.03)
cb.outline.set_edgecolor(style.GREY)
cb.ax.tick_params(colors=style.GREY, labelsize=8)
fig.tight_layout()
style.save(fig, "intro_f2_potentials", style.OUT_INTRO)
raise SystemExit(0 if style.all_ok() else 1)
