"""Lectures on integrable defects (defects-gemini), Lecture 1, figure 2 (fig-lec1-ising):
The Free-Fermion Escape Hatch: the Ising defect (sec-lec1-ising, eq-lec1-ising-rt).

In the off-critical Ising model (S = -1), the mixed-channel obstructions S(theta_1 + theta_2) - S(theta_1 - theta_2) = 0
and 1 - S(theta_1 - theta_2) S(theta_1 + theta_2) = 0 vanish identically, allowing reflection R(theta) and transmission
T(theta) to coexist across a line of modified energy bonds of strength g.

Left panel (a): Transmission probability |T(theta)|^2 (solid) and reflection probability |R(theta)|^2 (dashed) against
rapidity theta for defect couplings g = 0.5, g = 2(sqrt(2)-1) (where sin^2 chi = cos^2 chi = 1/2), and g = 1.4.
Right panel (b): Asymptotic transmission |T(infinity)|^2 = cos^2 chi = ((1 - g^2/4)/(1 + g^2/4))^2 and the imaginary
part of the pole rapidity theta_pole = i chi = -i arcsin(g/(1 + g^2/4)) as g varies from -2 (free boundary) through 0
(transparent line) to +2 (fixed boundary).

Checks:
  1. Single-particle unitarity R(theta) R(-theta) + T(theta) T(-theta) = 1 and R(theta) T(-theta) + T(theta) R(-theta) = 0
     (eq-lec1-ising-unitarity) and |R(theta)|^2 + |T(theta)|^2 = 1 for real theta.
  2. Separate crossing symmetries T(i*pi - theta) = T(theta), R(i*pi - theta) = -R(theta).
  3. Special limits: g = 0 gives T = 1, R = 0; g = +-2 gives T = 0 and R = -+ i cosh(theta)/(sinh(theta) +- i).
  4. Mixed Channels A and B (eq-lec1-mixed-a, eq-lec1-mixed-b) vanish identically for S = -1 with R*T != 0.

Writes lecture-notes/defects-gemini/figures/l1_ising.svg.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects-gemini" / "figures"


def chi_from_g(g):
    g = np.asarray(g, dtype=float)
    s = -g / (1.0 + 0.25 * g**2)
    c = (1.0 - 0.25 * g**2) / (1.0 + 0.25 * g**2)
    return s, c


def ising_rt(theta, g):
    s, c = chi_from_g(g)
    sh = np.sinh(theta)
    ch = np.cosh(theta)
    denom = sh - 1j * s
    T = c * sh / denom
    R = 1j * s * ch / denom
    return T, R


# ---- Checks ----
test_thetas = np.array([0.23, 0.85, 1.42, 0.61 + 0.33j, -0.44 + 0.19j])
test_gs = [-1.7, -0.6, 0.4, 2.0 * (np.sqrt(2.0) - 1.0), 1.5]

ok_unit = True
ok_cross = True
for g in test_gs:
    T_p, R_p = ising_rt(test_thetas, g)
    T_m, R_m = ising_rt(-test_thetas, g)
    u1 = R_p * R_m + T_p * T_m
    u2 = R_p * T_m + T_p * R_m
    ok_unit &= np.allclose(u1, 1.0) and np.allclose(u2, 0.0)

    # Real theta probability conservation
    th_real = np.linspace(0.05, 3.0, 50)
    Tr, Rr = ising_rt(th_real, g)
    ok_unit &= np.allclose(np.abs(Tr) ** 2 + np.abs(Rr) ** 2, 1.0)

    # Crossing relations
    T_shift, R_shift = ising_rt(0.5j * np.pi - test_thetas, g)
    T_pi, R_pi = ising_rt(1j * np.pi - test_thetas, g)
    # Note: sinh(i*pi/2 - th) = i cosh(th), cosh(i*pi/2 - th) = -i sinh(th)
    ok_cross &= np.allclose(T_pi, T_p) and np.allclose(R_pi, -R_p)

report(
    "Ising defect amplitudes satisfy unitarity R(th)R(-th)+T(th)T(-th)=1, R(th)T(-th)+T(th)R(-th)=0 and |R|^2+|T|^2=1",
    ok_unit,
)
report(
    "Ising defect amplitudes satisfy crossing T(i*pi - th) = T(th) and R(i*pi - th) = -R(th)",
    ok_cross,
)

# Check limits g = 0 and g = +-2
T0, R0 = ising_rt(test_thetas, 0.0)
Tp2, Rp2 = ising_rt(test_thetas, 2.0)
Tm2, Rm2 = ising_rt(test_thetas, -2.0)
report(
    "special limits: g=0 is transparent (T=1, R=0); g=+-2 decouples into two half-line boundaries (T=0, |R|=1)",
    np.allclose(T0, 1.0)
    and np.allclose(R0, 0.0)
    and np.allclose(Tp2, 0.0)
    and np.allclose(Tm2, 0.0)
    and np.allclose(Rp2, -1j * np.cosh(test_thetas) / (np.sinh(test_thetas) + 1j))
    and np.allclose(Rm2, 1j * np.cosh(test_thetas) / (np.sinh(test_thetas) - 1j)),
)

# Check Mixed Channels A and B for S = -1
S_ising = lambda th: -1.0
th1_v, th2_v = 1.15, 0.42
T1_v, R1_v = ising_rt(th1_v, 0.9)
T2_v, R2_v = ising_rt(th2_v, 0.9)
mixed_a = R2_v * T1_v * (S_ising(th1_v + th2_v) - S_ising(th1_v - th2_v))
mixed_b = R1_v * T2_v * (1.0 - S_ising(th1_v - th2_v) * S_ising(th1_v + th2_v))
report(
    "free-fermion loophole: at S = -1, Mixed Channels A and B vanish identically with R*T != 0",
    abs(R1_v * T2_v) > 0.1 and abs(mixed_a) == 0.0 and abs(mixed_b) == 0.0,
)

# ---- Figure ----
fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.0, 2.85), gridspec_kw=dict(wspace=0.32))

# Panel (a): |T(theta)|^2 and |R(theta)|^2 vs theta
thetas = np.linspace(0.0, 2.6, 300)
g_list = [
    (0.45, style.BLUE, r"$g=0.45$"),
    (2.0 * (np.sqrt(2.0) - 1.0), style.PURPLE, r"$g=2(\sqrt{2}-1)$"),
    (1.35, style.RED, r"$g=1.35$"),
]

for g_val, col, lab in g_list:
    T_curve, R_curve = ising_rt(thetas, g_val)
    axa.plot(thetas, np.abs(T_curve) ** 2, color=col, lw=1.8, label=lab)
    axa.plot(thetas, np.abs(R_curve) ** 2, color=col, lw=1.4, ls="--")

axa.axhline(0.5, color=style.LIGHTGREY, lw=0.8, ls=":")
axa.text(2.52, 0.94, r"solid: $|\mathbb{T}(\theta)|^2$, dashed: $|R(\theta)|^2$", fontsize=7.3, color=style.GREY, ha="right", va="top")
axa.set_xlim(0.0, 2.6)
axa.set_ylim(-0.02, 1.04)
axa.set_xlabel(r"rapidity $\theta$")
axa.set_ylabel(r"probability $|\mathbb{T}(\theta)|^2,\; |R(\theta)|^2$")
axa.set_title(r"(a) Partial transmission & reflection ($S=-1$)", fontsize=8.5, color=style.GREY)
axa.legend(loc="center right", bbox_to_anchor=(1.0, 0.43), fontsize=7.4, handlelength=1.8)

# Panel (b): interpolation from g = -2 (free boundary) to g = +2 (fixed boundary)
gs = np.linspace(-3.2, 3.2, 400)
s_chi, c_chi = chi_from_g(gs)
T_inf_sq = c_chi**2
chi_over_pi = np.arcsin(s_chi) / np.pi  # in [-1/2, 1/2]

axb.axvspan(-2.0, 0.0, color=style.GREEN, alpha=0.07, lw=0)
axb.axhline(0.0, color=style.LIGHTGREY, lw=0.8)
axb.axvline(-2.0, color=style.LIGHTGREY, lw=0.8, ls=":")
axb.axvline(2.0, color=style.LIGHTGREY, lw=0.8, ls=":")

axb.plot(gs, T_inf_sq, color=style.BLUE, lw=1.9, label=r"$|\mathbb{T}(\infty)|^2 = \cos^2\chi$")
axb.plot(gs, chi_over_pi, color=style.ORANGE, lw=1.8, label=r"pole $\operatorname{Im}\theta_{\rm pole}/\pi = \chi/\pi$")
axb.plot([-2.0, 2.0], [0.0, 0.0], "o", color=style.RED, ms=4.8, zorder=5)
axb.plot([0.0], [1.0], "o", color=style.BLUE, ms=4.8, zorder=5)

axb.annotate(
    r"free ($g=-2$)" + "\n" + r"$\mathbb{T}=0$",
    xy=(-2.0, 0.0),
    xytext=(-3.1, -0.30),
    fontsize=7.3,
    color=style.RED,
    arrowprops=dict(arrowstyle="->", color=style.RED, lw=0.7, shrinkB=4),
)
axb.annotate(
    r"fixed ($g=+2$)" + "\n" + r"$\mathbb{T}=0$",
    xy=(2.0, 0.0),
    xytext=(1.05, 0.55),
    fontsize=7.3,
    color=style.RED,
    arrowprops=dict(arrowstyle="->", color=style.RED, lw=0.7, shrinkB=4),
)
axb.text(-1.05, 0.86, r"bound state" + "\n" + r"($\chi>0$)", fontsize=7.3, color=style.GREEN, ha="center")

axb.set_xlim(-3.2, 3.2)
axb.set_ylim(-0.58, 1.08)
axb.set_xlabel(r"defect bond coupling $g$")
axb.set_ylabel(r"$|\mathbb{T}(\infty)|^2$ and pole $\operatorname{Im}\theta_{\rm pole}/\pi$")
axb.set_title(r"(b) Transparent ($g=0$) to boundaries ($g=\pm2$)", fontsize=8.5, color=style.GREY)
axb.legend(loc="lower left", fontsize=7.2, handlelength=1.8)

fig.subplots_adjust(left=0.08, right=0.98, bottom=0.17, top=0.89, wspace=0.32)
style.save(fig, "l1_ising", OUT)
raise SystemExit(0 if style.all_ok() else 1)
