"""Figure F38 (chapter 36, sec-doubling-trick, eq-doubled-k, prp-doubling): the doubling trick.

Left: the defect theory on the line, phi on x < 0 and psi on x > 0, with two world lines: a phi-particle of rapidity
theta transmitted to the right (TT(theta)), and a psi-particle of rapidity -theta' transmitted to the left
(TTbar(theta')).  Right: the doubled theory on x <= 0 with Phi(x) = (phi(x), psi(-x)); the psi-pieces of the world
lines are mirrored, x -> -x, rapidity theta -> -theta, and the defect becomes a boundary whose reflection matrix is
off-diagonal, K = (0, TTbar P; P TT, 0) on (V + V') (x) W.  The mirrored copy V' has opposite chirality: right-movers
become left-movers, and its bulk S-matrix is the parity image P Scheck P.

Checks: (1) x -> -x maps a plane wave of rapidity theta to one of rapidity -theta, and a right-mover f(t - x) to a
left-mover, i.e. d_+ <-> d_- (sympy); the mirrored world lines of the figure are the images of the psi-pieces;
(2) for the sine-Gordon Rcheck (boundary_common.rcheck, Part VI coproduct), P Rcheck P intertwines the opposite
coproduct and not the coproduct; (3) for the Lee-Yang defect (eq-ly-transmission) and the classical sine-Gordon
factor (eq-shg-linear-transmission), the doubled K satisfies boundary unitarity K(theta) K(-theta) = 1 and the scalar
form of boundary crossing-unitarity, K_{V'V}(i pi/2 - theta) = K_{VV'}(i pi/2 + theta), exactly because TT satisfies
defect unitarity and crossing; a modification that keeps unitarity and breaks crossing fails the second.
The full matrix statement (prp-doubling item 1) is checked in code/defect_axioms.py (d).

Writes part8-defects/figures/ch37_f38_doubling_trick.svg.
"""
import sys
from pathlib import Path

import numpy as np
import sympy as sp

import style
from style import plt, report

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from boundary_common import cop, q_phys, rcheck, sl2_spin  # noqa: E402

style.setup()
TH, THP = 0.8, 0.55                    # rapidity of the phi-particle, and theta' of the psi-particle
T1, T2 = -0.9, 0.9                     # times at which they reach the defect
TMIN, TMAX = -2.2, 2.2

# ---- checks ----------------------------------------------------------------
x, t, th, m = sp.symbols('x t theta m', real=True)
f = sp.Function('f')
wave = sp.exp(-sp.I * m * (sp.cosh(th) * t - sp.sinh(th) * x))
mirrored = wave.subs(x, -x)
report("x -> -x maps the plane wave of rapidity theta to the plane wave of rapidity -theta",
       sp.simplify(mirrored - wave.subs(th, -th)) == 0)
g = f(t - x)
gm = g.subs(x, -x)                    # f(t + x)
dp = lambda h: sp.diff(h, t) + sp.diff(h, x)
dm = lambda h: sp.diff(h, t) - sp.diff(h, x)
report("x -> -x turns a right-mover f(t - x) (d_+ f = 0) into a left-mover (d_- = 0), i.e. d_+ <-> d_-",
       sp.simplify(dp(g)) == 0 and sp.simplify(dm(gm)) == 0 and sp.simplify(dp(gm)) != 0)


def wl(th_, t0, t_from, t_to, x0=0.0):
    """World line x = x0 + tanh(th_) (t - t0) between two times."""
    tt = np.array([t_from, t_to])
    return x0 + np.tanh(th_) * (tt - t0), tt


x_psi1, t_psi1 = wl(TH, T1, T1, TMAX)                  # transmitted psi-particle, rapidity theta, x > 0
x_psi2, t_psi2 = wl(-THP, T2, TMIN, T2)                # incoming psi-particle, rapidity -theta', x > 0
x_m1, t_m1 = wl(-TH, T1, T1, TMAX)                     # mirrored pieces in the doubled picture
x_m2, t_m2 = wl(THP, T2, TMIN, T2)
report("the dashed world lines of the doubled picture are the mirror images x -> -x of the psi-pieces, with rapidities "
       "theta -> -theta (right-mover <-> left-mover)",
       np.allclose(x_m1, -x_psi1) and np.allclose(x_m2, -x_psi2) and np.all(x_psi1 >= 0) and np.all(x_m1 <= 0)
       and np.sign(np.diff(x_psi1)[0]) == -np.sign(np.diff(x_m1)[0]))

q = q_phys(1.37)
A, B = sl2_spin(1, q, np.exp(1.37 * 0.4)), sl2_spin(1, q, np.exp(1.37 * -0.3))
R, nul = rcheck(A, B)                                  # A (x) B -> B (x) A, intertwines Delta
P = np.zeros((4, 4))
for i in range(2):
    for j in range(2):
        P[j * 2 + i, i * 2 + j] = 1
PRP = P @ R @ P                                        # B (x) A -> A (x) B


def res(M, X_in, X_out):
    return np.abs(M @ X_in - X_out @ M).max()


r_op = max(res(PRP, P @ cop(A, B, k, j) @ P, P @ cop(B, A, k, j) @ P) for j in (0, 1) for k in 'efk')
r_co = max(res(PRP, cop(B, A, k, j), cop(A, B, k, j)) for j in (0, 1) for k in 'efk')
report(f"sine-Gordon: P Rcheck P intertwines the opposite coproduct [res {r_op:.0e}] and not the coproduct "
       f"[res {r_co:.2f}]: the mirrored copy V' has the opposite chirality", nul == 1 and r_op < 1e-10 and r_co > 1e-3)

S_ly = lambda z: (np.sinh(z) + 1j * np.sin(2 * np.pi / 3)) / (np.sinh(z) - 1j * np.sin(2 * np.pi / 3))
Tcl = lambda thd: (lambda z: -1j * np.sinh((z - thd) / 2 - 1j * np.pi / 4) / np.sinh((z - thd) / 2 + 1j * np.pi / 4))
cases = {}
for b in (0.4, -1.3, 1.7 + 0.5j):
    cases[f"Lee-Yang b = {b}"] = (lambda z, b=b: S_ly(z - 1j * np.pi * (3 - b) / 6),
                                  lambda z, b=b: S_ly(z + 1j * np.pi * (3 - b) / 6))
cases["classical sine-Gordon, theta_d = 0.3"] = (Tcl(0.3), lambda z: Tcl(0.3)(1j * np.pi - z))
pts = np.array([0.37 + 0.11j, -0.8 + 0.3j, 1.3 - 0.2j, 0.05j])


def Kd(T, Tb, z):
    return np.array([[0, Tb(z)], [T(z), 0]])


def doubled(T, Tb):
    eu = max(np.abs(Kd(T, Tb, z) @ Kd(T, Tb, -z) - np.eye(2)).max() for z in pts)
    ec = max(abs(Kd(T, Tb, 1j * np.pi / 2 - z)[1, 0] - Kd(T, Tb, 1j * np.pi / 2 + z)[0, 1]) for z in pts)
    du = max(abs(T(z) * Tb(-z) - 1) for z in pts)
    dc = max(abs(T(z) - Tb(1j * np.pi - z)) for z in pts)
    return eu, ec, du, dc


ok = True
for name, (T, Tb) in cases.items():
    eu, ec, du, dc = doubled(T, Tb)
    ok &= max(eu, ec, du, dc) < 1e-12
report("doubled K = (0, TTbar; TT, 0): K(theta) K(-theta) = 1 and K_{V'V}(i pi/2 - theta) = K_{VV'}(i pi/2 + theta) "
       "(boundary unitarity and crossing-unitarity) for " + ", ".join(cases), ok)
blk = lambda xx: (lambda z: np.sinh(z / 2 + 1j * np.pi * xx / 12) / np.sinh(z / 2 - 1j * np.pi * xx / 12))
T0, Tb0 = cases["Lee-Yang b = 0.4"]
eu, ec, du, dc = doubled(lambda z: T0(z) * blk(1.3)(z), lambda z: Tb0(z) * blk(1.3)(z))
report(f"control: TT, TTbar times a block (1.3) keep unitarity [{max(eu, du):.0e}] and break doubled crossing "
       f"[{ec:.2f}] together with defect crossing [{dc:.2f}]", max(eu, du) < 1e-12 and ec > 1e-3 and dc > 1e-3)

# ---- figure ----------------------------------------------------------------
v1, v2 = np.tanh(TH), np.tanh(THP)
tc = (v2 * T2 + v1 * T1) / (v1 + v2)                   # the two psi-particles meet at (xc, tc) (bulk scattering)
xc = v1 * (tc - T1)
fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.0, 3.1), gridspec_kw=dict(width_ratios=[1.5, 1.0], wspace=0.5))
XL = 2.2
for ax in (axa, axb):
    ax.set_ylim(TMIN, TMAX)
    ax.set_xticks([0])
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(r"$x$", labelpad=1)
axa.set_ylabel(r"$t$")

# (a) the defect theory
axa.set_xlim(-XL, XL)
axa.axvspan(-XL, 0, color=style.BLUE, alpha=0.07, lw=0)
axa.axvspan(0, XL, color=style.RED, alpha=0.07, lw=0)
axa.axvline(0, color=style.PURPLE, lw=2.0)
axa.plot(*wl(TH, T1, TMIN, T1), color=style.BLUE)
axa.plot(x_psi1, t_psi1, color=style.RED)
axa.plot(x_psi2, t_psi2, color=style.RED)
axa.plot(*wl(-THP, T2, T2, TMAX), color=style.BLUE)
for tt in (T1, T2):
    axa.plot([0], [tt], "o", color=style.PURPLE, ms=4, zorder=5)
axa.plot([xc], [tc], "o", color=style.ORANGE, ms=4, zorder=5)
axa.text(xc + 0.12, tc, r"Š", color=style.ORANGE, ha="left", va="center", fontsize=8.5,
         fontstyle="italic")  # S-caron = check S (mathtext misplaces the accent)
axa.text(0.1, T1 - 0.05, r"$\mathbb{T}(\theta)$", ha="left", va="top", fontsize=8)
axa.text(-0.1, T2 + 0.05, r"$\bar{\mathbb{T}}(\theta')$", ha="right", va="bottom", fontsize=8)
axa.text(-1.5, TMIN + 0.12, r"$\vec\phi$ on $x<0$", color=style.BLUE, ha="center", va="bottom", fontsize=8)
axa.text(0.72, TMIN + 0.12, r"$\vec\psi$ on $x>0$", color=style.RED, ha="center", va="bottom", fontsize=8)
axa.text(-0.85, -1.75, r"$\theta$", color=style.BLUE, ha="right", fontsize=8)
axa.text(1.25, 0.75, r"$\theta$", color=style.RED, ha="right", fontsize=8)
axa.text(1.15, -0.75, r"$-\theta'$", color=style.RED, ha="left", fontsize=8)
axa.text(-0.55, 1.6, r"$-\theta'$", color=style.BLUE, ha="right", fontsize=8)
axa.set_title(r"(a) defect on the line", fontsize=8)

# (b) the doubled theory on the half-line
axb.set_xlim(-XL, 0.75)
axb.axvspan(-XL, 0, color=style.GREY, alpha=0.05, lw=0)
axb.axvline(0, color=style.PURPLE, lw=2.0)
for yy in np.linspace(TMIN, TMAX - 0.16, 19):
    axb.plot([0, 0.14], [yy, yy + 0.16], color=style.PURPLE, lw=0.8)
axb.plot(*wl(TH, T1, TMIN, T1), color=style.BLUE)
axb.plot(x_m1, t_m1, color=style.RED, ls="--")
axb.plot(x_m2, t_m2, color=style.RED, ls="--")
axb.plot(*wl(-THP, T2, T2, TMAX), color=style.BLUE)
for tt in (T1, T2):
    axb.plot([0], [tt], "o", color=style.PURPLE, ms=4, zorder=5)
axb.plot([-xc], [tc], "o", color=style.ORANGE, ms=4, zorder=5)
axb.text(-xc - 0.12, tc, r"$P$Š$P$", color=style.ORANGE, ha="right", va="center", fontsize=8.5,
         fontstyle="italic")
axb.text(0.2, T1, r"$P\,\mathbb{T}$", ha="left", va="center", fontsize=8)
axb.text(0.2, T2, r"$\bar{\mathbb{T}}P$", ha="left", va="center", fontsize=8)
axb.text(-1.0, 0.45, r"$-\theta$", color=style.RED, ha="right", fontsize=8)
axb.text(-1.55, -1.95, r"$\theta'$", color=style.RED, ha="right", fontsize=8)
axb.set_title(r"(b) doubled theory on $x\leq0$", fontsize=8)

# the mirror arrow between the panels; legend and reflection matrix below panel (b)
fig.subplots_adjust(left=0.04, right=0.99, bottom=0.2, top=0.92)
pa, pb = axa.get_position(), axb.get_position()
xm, ym = (pa.x1 + pb.x0) / 2, (pa.y0 + pa.y1) / 2
fig.text(xm, ym + 0.05, r"$x\to-x$" + "\n" + r"on $\vec\psi$", ha="center", va="bottom", fontsize=8, color=style.GREY)
fig.patches.append(plt.matplotlib.patches.FancyArrowPatch((xm - 0.04, ym), (xm + 0.04, ym), transform=fig.transFigure,
                                                          arrowstyle="->", mutation_scale=10, color=style.GREY))
axa.plot([], [], color=style.BLUE, label=r"copy $V$: $\vec\phi(x)$")
axa.plot([], [], color=style.RED, ls="--", label=r"copy $V'$: $\vec\psi(-x)$, $\theta\to-\theta$")
fig.legend(*axa.get_legend_handles_labels(), loc="lower left", bbox_to_anchor=(pa.x0, 0.0), ncol=2, fontsize=7.5,
           handlelength=2.0, columnspacing=1.5)
fig.text(pb.x0 + 0.5 * (pb.x1 - pb.x0), 0.005,
         r"$K=\left(\genfrac{}{}{0}{0}{0}{P\,\mathbb{T}}\ \ \genfrac{}{}{0}{0}{\bar{\mathbb{T}}P}{0}\right)$",
         ha="center", va="bottom", fontsize=9.5)
style.save(fig, "ch37_f38_doubling_trick", style.OUT_PART8)
raise SystemExit(0 if style.all_ok() else 1)
