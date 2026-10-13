"""Figure F41 (chapter 42, sec-jump-coideal, prp-diagonal-borel, tbl-dbd-symmetries): the symmetries of the doubled
boundary of sine-Gordon lines, and how the type-I_- jump defect interpolates between two invertible lines.

Top: as theta_d -> +infinity the transmission matrix tends, on each defect charge n, to a twist line, whose doubled
boundary keeps the twisted diagonal (1 (x) t_n) Delta(U); the transparent line is the twist line with t = 1.  As
theta_d -> -infinity it tends to a charge conjugation dressed by a shift of the defect charge; the charge-conjugation
line C itself keeps (1 (x) tau) Delta(U).  At finite theta_d the doubled boundary keeps one Borel half of each:
Delta(U_q(b^-)) on the Weyl module rho and (1 (x) tau) Delta(U_q(b^+)) on the second module rho~
(eq-dbd-second-module); no twisted diagonal survives, and the other halves do not return in the limits on all of W
(prp-diagonal-borel (5): the sector twist t_n depends on n, and the charge shift cannot be absorbed into W).
The arrows mark matrix limits, not restored symmetries.  Bottom: the probability that a
soliton is converted, |Phi C_n|^2 = Y^2/(1 + Y^2), Y = e^{lambda(theta - theta_d)}, against lambda(theta_d - theta).

Checks: the transmission matrix of defect_common.sg_tmatrix has |C_n/A_n| = Y on every label, and with the scalar
factor (eq-sgd-scalar) |Phi A_n|^2 + |Phi C_n|^2 = 1 and |Phi C_n|^2 = Y^2/(1 + Y^2) (lambda = 1.37, three theta_d);
the limits of prp-defect-twist-limit: for theta - theta_d -> -infinity the V -> V block on a fixed label is diagonal,
diag(A_n, D_n), and conjugates pi(e_1), pi(e_0), pi(k_1) into t pi(e_1), t^-1 pi(e_0), pi(k_1) with t = A_n/D_n
(a torus twist t_1 = 1/t_0, depending on n); for theta - theta_d -> +infinity the block is antidiagonal and
conjugates pi(e_1) into a multiple of pi(e_0) = pi(tau(e_1)) (charge conjugation, as a matrix on each total-charge
sector, not a line with W = C); the twist line diag(e^{i alpha},
e^{-i alpha}) gives t_1 = 1/t_0 = e^{2 i alpha}, and C = sigma_x gives tau exactly.

Writes part8-defects/figures/ch42_f41_doubled_symmetries.svg.
"""
import sys
from pathlib import Path

import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

import style
from style import plt, report

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from boundary_common import q_phys, sl2_spin  # noqa: E402
from defect_common import sg_tmatrix, tmatrix_elements  # noqa: E402

style.setup()
LAM = 1.37
q = q_phys(LAM)

# ---- checks ----------------------------------------------------------------
ok_ratio, ok_prob = True, True
for thd in (-0.8, 0.3, 1.1):
    th = 0.2
    Y = np.exp(LAM * (th - thd))
    T0, V, W = sg_tmatrix(th, LAM, q, thd, 4, scalar=False)
    A, B, C, D = tmatrix_elements(T0, V, W)
    ok_ratio &= all(abs(abs(C[n] / A[n]) - Y) < 1e-12 for n in C)
    T1, _, _ = sg_tmatrix(th, LAM, q, thd, 4, scalar=True)
    A1, _, C1, _ = tmatrix_elements(T1, V, W)
    ok_prob &= all(abs(abs(A1[n])**2 + abs(C1[n])**2 - 1) < 1e-8 and abs(abs(C1[n])**2 - Y**2 / (1 + Y**2)) < 1e-8
                   for n in C1)
report("sg_tmatrix (eq-sg-transmission-matrix): |C_n/A_n| = Y = e^{lambda(theta - theta_d)} on every label", ok_ratio)
report("with the scalar factor (eq-sgd-scalar): |Phi A_n|^2 + |Phi C_n|^2 = 1 and the conversion probability "
       "|Phi C_n|^2 = Y^2/(1 + Y^2) (lambda = 1.37, theta_d = -0.8, 0.3, 1.1)", ok_prob)

Vy = sl2_spin(1, q, np.exp(LAM * 0.2))
e1, e0, k1 = Vy[('e', 1)], Vy[('e', 0)], Vy[('k', 1)]
ok = True
for n in (-2, 0, 3):
    for s in (0.0, 0.3):
        An, Dn = q**-n, q**(n + s)
        M = np.diag([An, Dn])
        Mi = np.linalg.inv(M)
        t = An / Dn
        ok &= np.allclose(M @ e1 @ Mi, t * e1) and np.allclose(M @ e0 @ Mi, e0 / t) and np.allclose(M @ k1 @ Mi, k1)
        ok &= abs(abs(t) - 1) < 1e-12
report("theta - theta_d -> -inf: diag(A_n, D_n) conjugates pi(e_1), pi(e_0), pi(k_1) to t pi(e_1), t^-1 pi(e_0), pi(k_1), "
       "t = q^{-2n-s} a phase: a twist line with t_1 = 1/t_0 depending on the defect charge n", ok)
ok = True
for n in (-2, 0, 3):
    for s in (0.0, 0.3):
        Cn, Bn1 = -1j * q**-n, -1j * q**(n + 1 + s)          # v_+ -> v_- (n -> n+1) and v_- -> v_+ (n+1 -> n)
        M = np.array([[0, Bn1], [Cn, 0]])
        Mi = np.linalg.inv(M)
        X = M @ e1 @ Mi
        c = X[1, 0] / e0[1, 0] if abs(e0[1, 0]) > 0 else X[0, 1] / e0[0, 1]
        ok &= np.allclose(X, c * e0) and np.allclose(M @ k1 @ Mi, np.linalg.inv(k1))
report("theta - theta_d -> +inf: the antidiagonal block conjugates pi(e_1) to a multiple of pi(e_0) = pi(tau(e_1)) and "
       "pi(k_1) to pi(k_0): charge conjugation dressed by a torus twist", ok)
al = 0.4
Mt, Cx = np.diag([np.exp(1j * al), np.exp(-1j * al)]), np.array([[0, 1], [1, 0]])
report("twist line diag(e^{i alpha}, e^{-i alpha}): t_1 = 1/t_0 = e^{2 i alpha}; C = sigma_x maps pi(e_1), pi(k_1) "
       "exactly to pi(e_0), pi(k_0) (the diagram automorphism tau)",
       np.allclose(Mt @ e1 @ np.linalg.inv(Mt), np.exp(2j * al) * e1) and np.allclose(Mt @ e0 @ np.linalg.inv(Mt), np.exp(-2j * al) * e0)
       and np.allclose(Cx @ e1 @ Cx, e0) and np.allclose(Cx @ k1 @ Cx, Vy[('k', 0)]))

# ---- figure ----------------------------------------------------------------
fig = plt.figure(figsize=(7.0, 3.2))
ax = fig.add_axes([0.0, 0.42, 1.0, 0.58])
ax.set_xlim(-4.4, 4.4)
ax.set_ylim(0.85, 3.0)
ax.axis("off")


def box(xc, yc, w, h, ec, fc="none", lw=1.0, ls="-"):
    ax.add_patch(FancyBboxPatch((xc - w / 2, yc - h / 2), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                                ec=ec, fc=fc, lw=lw, ls=ls))


def half(xc, yc, w, h, text, col, kept):
    box(xc, yc, w, h, col if kept else style.GREY, fc=(col + "26") if kept else "none", lw=1.4 if kept else 0.8,
        ls="-" if kept else (0, (3, 2)))
    ax.text(xc, yc, text, ha="center", va="center", fontsize=7.5, color=col if kept else style.GREY)


# left end: charge conjugation
xl, xr = -3.1, 3.1
ax.text(xl, 2.75, r"$\theta_d\to-\infty$: $C$ with a charge shift", ha="center", va="center", fontsize=8)
ax.text(xl, 2.38, r"line $C$: $(1\otimes\tau)\Delta(U)$", ha="center", va="center", fontsize=8)
half(xl - 0.6, 1.8, 1.15, 0.42, r"$\mathfrak{b}^+$ half", style.ORANGE, True)
half(xl + 0.6, 1.8, 1.15, 0.42, r"$\mathfrak{b}^-$ half", style.ORANGE, False)
# right end: twist lines
ax.text(xr, 2.75, r"$\theta_d\to+\infty$: twist lines $e^{i\alpha_n\cdot\mu}$", ha="center", va="center", fontsize=8)
ax.text(xr + 0.3, 2.38, r"$(1\otimes t_n)\Delta(U)$, each $n$", ha="center", va="center", fontsize=8)
ax.text(xr, 1.42, r"$t=1$: the transparent line, $\Delta(U)$", ha="center", va="center", fontsize=7.2)
half(xr - 0.6, 1.8, 1.15, 0.42, r"$\mathfrak{b}^-$ half", style.BLUE, True)
half(xr + 0.6, 1.8, 1.15, 0.42, r"$\mathfrak{b}^+$ half", style.BLUE, False)
# centre: the jump defect
ax.text(0, 2.75, r"type-I$_-$ jump, finite $\theta_d$", ha="center", va="center", fontsize=8)
half(0, 2.12, 2.5, 0.42, r"$\Delta(U_q(\mathfrak{b}^-))$ on $\rho$", style.BLUE, True)
half(0, 1.55, 2.5, 0.42, r"$(1\otimes\tau)\Delta(U_q(\mathfrak{b}^+))$ on $\tilde\rho$", style.ORANGE, True)
ax.text(0, 1.08, "no twisted diagonal, no coideal with\nboth halves; dashed halves stay lost", ha="center", va="center",
        fontsize=7.2)
ax.add_patch(FancyArrowPatch((1.3, 2.12), (xr - 1.2, 1.86), arrowstyle="->", mutation_scale=9, color=style.BLUE, lw=1.0,
                             connectionstyle="arc3,rad=-0.15"))
ax.add_patch(FancyArrowPatch((-1.3, 1.55), (xl + 0.03, 1.57), arrowstyle="->", mutation_scale=9, color=style.ORANGE,
                             lw=1.0, connectionstyle="arc3,rad=-0.25"))
ax.text(1.5, 2.36, r"matrix limit $\theta_d\to+\infty$", color=style.BLUE, ha="center", va="bottom", fontsize=7)
ax.text(-3.0, 1.22, r"matrix limit $\theta_d\to-\infty$", color=style.ORANGE, ha="center", va="top", fontsize=7)

# bottom: conversion probability
ax2 = fig.add_axes([0.1, 0.12, 0.82, 0.25])
u = np.linspace(-4, 4, 400)
ax2.plot(u, 1 / (1 + np.exp(2 * u)), color=style.PURPLE)
ax2.set_xlim(-4, 4)
ax2.set_ylim(-0.05, 1.05)
ax2.set_yticks([0, 0.5, 1])
ax2.set_xlabel(r"$\lambda(\theta_d-\theta)$")
ax2.set_ylabel("converted", fontsize=7.5)
ax2.text(-3.9, 0.45, r"$|\Phi C_n|^2=\frac{Y^2}{1+Y^2}$, $Y=e^{\lambda(\theta-\theta_d)}$", fontsize=7.5, va="center")
ax2.text(3.9, 0.85, "diagonal: twist", ha="right", fontsize=7.2, color=style.BLUE)
ax2.text(-3.9, 0.85, "off-diagonal: conjugation", ha="left", fontsize=7.2, color=style.ORANGE)
style.save(fig, "ch42_f41_doubled_symmetries", style.OUT_PART8)
raise SystemExit(0 if style.all_ok() else 1)
