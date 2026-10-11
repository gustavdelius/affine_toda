"""Figure F9 (chapter 4, sec-krein-spaces, prp-krein-spectrum item 4): why the Hilbert norm
can grow while the evolution stays unitary for the Krein form.

(a) H Hermitian (H = sigma_y): U(t) = e^{-iHt} is a rotation, orbits are circles and the
    Hilbert norm is conserved.
(b) H = i sigma_x is J-self-adjoint for J = sigma_z but not Hermitian: U(t) = e^{-iHt} =
    e^{sigma_x t} is a hyperbolic rotation. Orbits are hyperbolas x^2 - y^2 = const (the Krein
    form [v,v] = <v, J v>), the lines y = +-x are neutral, and the Hilbert norm of (1,0) grows as
    cosh 2t.  (c) the two norms along the orbit.

Writes part1-background/figures/ch4_f9_krein_geometry.svg.
"""
import numpy as np
from scipy.linalg import expm

import style
from style import plt, report

style.setup()
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], complex)
J = sz
Ha, Hb = sy, 1j * sx

# ---- checks ---------------------------------------------------------------
report("(a) H = sigma_y is Hermitian", np.allclose(Ha, Ha.conj().T))
report("(b) H = i sigma_x is not Hermitian", not np.allclose(Hb, Hb.conj().T))
report("(b) H is J-self-adjoint: H^dagger = J H J", np.allclose(Hb.conj().T, J @ Hb @ J))
ok_u = ok_k = ok_h = ok_ha = True
v0 = np.array([1.0, 0.0], complex)
for t in np.linspace(0, 3, 13):
    Ub = expm(-1j * Hb * t)
    ok_u &= np.allclose(Ub, np.cosh(t) * np.eye(2) + np.sinh(t) * sx)
    v = Ub @ v0
    ok_k &= abs(v.conj() @ J @ v - 1) < 1e-12                            # Krein norm conserved
    ok_h &= abs((v.conj() @ v).real - np.cosh(2 * t)) < 1e-12            # Hilbert norm = cosh 2t
    va = expm(-1j * Ha * t) @ v0
    ok_ha &= abs((va.conj() @ va).real - 1) < 1e-12
report("(b) e^{-iHt} = cosh t + sigma_x sinh t", ok_u)
report("(b) Krein norm [v,v] = 1 is conserved", ok_k)
report("(b) Hilbert norm of (1,0) equals cosh 2t", ok_h)
report("(a) Hilbert norm is conserved for Hermitian H", ok_ha)

# ---- figure ---------------------------------------------------------------
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(7.0, 2.8), gridspec_kw={"width_ratios": [1, 1, 1.1]})
L = 2.3
ph = np.linspace(0, 2 * np.pi, 400)
for ax in (ax1, ax2):
    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
for r0 in (0.6, 1.1, 1.6, 2.0):
    ax1.plot(r0 * np.cos(ph), r0 * np.sin(ph), color=style.BLUE, lw=1.4)
ax1.annotate("", xy=(0.0, 1.05), xytext=(1.05, 0.2), arrowprops=dict(arrowstyle="->", color=style.GREY, lw=1.0,
                                                                       connectionstyle="arc3,rad=0.4"))
ax1.set_title(r"(a) Hermitian: circles" "\n" r"$x^2+y^2$ conserved", fontsize=8.5, color=style.GREY)
tt = np.linspace(-3.2, 3.2, 1200)
for c in (0.25, 1.0, 2.2):                       # positive cone: x^2 - y^2 = c
    for sg in (+1, -1):
        ax2.plot(sg * np.sqrt(c) * np.cosh(tt), np.sqrt(c) * np.sinh(tt), color=style.BLUE, lw=1.4)
for c in (0.25, 1.0, 2.2):                       # negative cone: x^2 - y^2 = -c
    for sg in (+1, -1):
        ax2.plot(np.sqrt(c) * np.sinh(tt), sg * np.sqrt(c) * np.cosh(tt), color=style.RED, lw=1.4)
ax2.plot([-L, L], [-L, L], "--", color=style.GREY, lw=1.0)
ax2.plot([-L, L], [L, -L], "--", color=style.GREY, lw=1.0)
ax2.text(1.5, 0.2, r"$[v,v]>0$", color=style.BLUE, fontsize=8)
ax2.text(-0.55, 1.95, r"$[v,v]<0$", color=style.RED, fontsize=8)
ax2.text(1.35, 1.9, "neutral", color=style.GREY, fontsize=7.5, rotation=40)
ax2.set_title(r"(b) Krein-unitary: hyperbolas" "\n" r"$x^2-y^2$ conserved", fontsize=8.5, color=style.GREY)
t = np.linspace(0, 2.0, 300)
ax3.plot(t, np.cosh(2 * t), color=style.PURPLE, label=r"Hilbert norm $\langle v,v\rangle=\cosh2t$")
ax3.plot(t, np.ones_like(t), color=style.BLUE, label=r"Krein norm $[v,v]=1$")
ax3.set_xlabel(r"$t$")
ax3.legend(loc="upper left", fontsize=7.5)
ax3.set_title(r"(c) along the orbit of $(1,0)$", fontsize=8.5, color=style.GREY)
fig.tight_layout()
style.save(fig, "ch4_f9_krein_geometry", style.OUT_PART1)
raise SystemExit(0 if style.all_ok() else 1)
