"""Lecture notes on defects, Lecture 5, figure 2 (fig-l5-rotation): from the identity to the rotation.

For a_3^(1) (h = 4), the phases that the type-I_- window gives to what crosses it, against theta - theta_d:
(a) linear waves of particle k = 1, 2, 3, T_k = (e^{-(theta - theta_d)} + e^{i gamma_k})/(e^{-(theta - theta_d)} +
    e^{-i gamma_k}), gamma_k = pi k/h, whose phase rises continuously from 0 to 2 gamma_k, passing gamma_k at
    theta = theta_d;
(b) solitons of species a = 1, 2, 3, z_a(u) = (e^{-u} - i e^{i gamma_a})/(e^{-u} - i e^{-i gamma_a}),
    u = theta - theta_d, whose phase arg(e^{-2u} + e^{2 i gamma_a}) goes from 0 to 2 gamma_a (a < h/2), jumps by pi
    (a = h/2, the sine-Gordon conversion), or goes to 2 gamma_a - 2 pi (a > h/2).
In both panels the end values are the eigenvalues omega^k = e^{2 pi i k/h} of the rotation pi: fast particles and
fast solitons see the window as the symmetry line psi = pi phi.

Checks: |T_k| = 1 on the plotted range; the end values are omega^k; arg T_k(theta_d) = gamma_k;
z_a(u) = T_a(u + theta_d - i pi/2); the jump of arg z_2 by pi at u = 0.
Writes lecture-notes/defects/figures/l5_rotation.svg.
"""
import cmath
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects" / "figures"

h = 4
thd = 0.0
u = np.linspace(-6, 6, 1201)


def T(uu, k):
    g = math.pi * k / h
    x = np.exp(-uu)
    return (x + np.exp(1j * g)) / (x + np.exp(-1j * g))


def z(uu, a):
    g = math.pi * a / h
    x = np.exp(-uu)
    return (x - 1j * np.exp(1j * g)) / (x - 1j * np.exp(-1j * g))


def phase_T(uu, k):
    """Continuous phase of T_k: 2 arg(e^{-u} + e^{i gamma_k}), in (0, 2 gamma_k)."""
    return 2 * np.angle(np.exp(-uu) + np.exp(1j * math.pi * k / h))


def phase_z(uu, a):
    return np.angle(np.exp(-2 * uu) + np.exp(2j * math.pi * a / h))


ok = all(np.max(np.abs(np.abs(T(u, k)) - 1)) < 1e-13 for k in (1, 2, 3))
ok &= all(np.max(np.abs(np.exp(1j * phase_T(u, k)) - T(u, k))) < 1e-12 for k in (1, 2, 3))
report("|T_k| = 1 on the plotted range, and the plotted phase is the phase of T_k", ok)
ends = all(abs(T(np.array([30.0]), k)[0] - cmath.exp(2j * math.pi * k / h)) < 1e-10 and
           abs(z(np.array([30.0]), k)[0] - cmath.exp(2j * math.pi * k / h)) < 1e-10 and
           abs(T(np.array([-30.0]), k)[0] - 1) < 1e-10 for k in (1, 2, 3))
ends &= all(abs(phase_T(np.array([0.0]), k)[0] - math.pi * k / h) < 1e-12 for k in (1, 2, 3))
report("the end values are 1 and omega^k for waves and solitons; arg T_k(theta_d) = gamma_k", ends)
cont = all(np.max(np.abs(z(uu := np.array([-1.3, 0.2, 0.9]), a) -
                         T(uu - 0.5j * math.pi, a))) < 1e-12 for a in (1, 2, 3))
jump = abs(phase_z(np.array([1e-6]), 2)[0] - phase_z(np.array([-1e-6]), 2)[0]) % (2 * math.pi)
cont &= abs(jump - math.pi) < 1e-6
report("z_a(u) = T_a(u + theta_d - i pi/2), and the phase of the self-conjugate soliton jumps by pi at u = 0", cont)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.8), sharey=True)
cols = {1: style.BLUE, 2: style.PURPLE, 3: style.RED}
for k in (1, 2, 3):
    ax1.plot(u, phase_T(u, k) / math.pi, color=cols[k], label=f"$k={k}$")
    ax1.axhline(2 * k / h, color=cols[k], ls=":", lw=0.8)
    ax1.plot(0, k / h, marker="o", ms=3.5, color=cols[k])
u_neg, u_pos = u[u < 0], u[u > 0]
ax2.plot(u, phase_z(u, 1) / math.pi, color=cols[1], label="$a=1$")
ax2.plot(u_neg, phase_z(u_neg, 2) / math.pi, color=cols[2], label="$a=2$")
ax2.plot(u_pos, phase_z(u_pos, 2) / math.pi, color=cols[2])
ax2.plot(u, phase_z(u, 3) / math.pi, color=cols[3], label="$a=3$")
ax2.text(5.8, 0.92, r"$\omega^2$", color=cols[2], fontsize=8, ha="right", va="top")
ax2.text(5.8, 0.42, r"$\omega$", color=cols[1], fontsize=8, ha="right", va="top")
ax2.text(5.8, -0.36, r"$\omega^3$", color=cols[3], fontsize=8, ha="right", va="bottom")
for y, k in ((0.5, 1), (1.0, 2), (-0.5, 3)):
    ax2.axhline(y, color=cols[k], ls=":", lw=0.8)
for ax in (ax1, ax2):
    ax.axhline(0, color=style.GREY, lw=0.6)
    ax.axvline(0, color=style.LIGHTGREY, lw=0.8)
    ax.set_xlim(u[0], u[-1])
    ax.set_xlabel(r"$\theta-\theta_d$")
    ax.legend(loc="upper left", fontsize=7)
ax1.set_ylim(-0.7, 1.7)
ax1.set_yticks([-0.5, 0, 0.5, 1, 1.5], [r"$-\frac{\pi}{2}$", "0", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$"])
ax1.set_ylabel("phase")
ax1.set_title(r"(a) linear waves: $\arg\mathbb{T}_k$", fontsize=8.5, color=style.GREY)
ax2.set_title(r"(b) solitons: $\arg z_a$", fontsize=8.5, color=style.GREY)
ax1.text(5.8, 1.42, r"$\omega^3$", color=cols[3], fontsize=8, ha="right", va="top")
ax1.text(5.8, 0.92, r"$\omega^2$", color=cols[2], fontsize=8, ha="right", va="top")
ax1.text(5.8, 0.42, r"$\omega$", color=cols[1], fontsize=8, ha="right", va="top")
fig.tight_layout(w_pad=0.8)
style.save(fig, "l5_rotation", OUT)
raise SystemExit(0 if style.all_ok() else 1)
