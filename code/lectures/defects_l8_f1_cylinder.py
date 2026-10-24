"""Lecture notes on defects, Lecture 8, figure 1 (fig-l8-cylinder): two ways to put a defect on a cylinder.

Space is a circle, time runs upwards, and three solitons (green) wind around the cylinder with rapidities
theta_j = -0.5, 0.1, 0.6, drawn with slopes tanh(theta_j).
(a) The defect along time (purple), at a point of the circle: one more site of the soliton chain. The transfer
matrix t(theta) is an auxiliary soliton (orange, dashed) that runs once around the circle and crosses the three
solitons with R and the defect with the L-operator of Lecture 7.
(b) The defect along space: at one instant it runs around the circle and crosses each soliton once, with L_j^-1; its
state |n> propagates along it and closes on itself. The train argument moves the auxiliary soliton of t(theta) through
it.

Checks: the defect conditions turned through 90 degrees are Baecklund's B_a with a = -i sigma (sympy); for the three
drawn solitons the operator <0| L_3^-1 L_2^-1 L_1^-1 |0>, in the order in which the defect meets them, commutes with
t_omega(theta) in every sector (omega = q^{-m/2}); the reversed order and the product without inverses do not.
Writes lecture-notes/defects/figures/l8_cylinder.svg.
"""
import cmath
import math
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects" / "figures"

# ---------------------------------------------------------------- checks
m0, sig = sp.symbols('m0 sigma', positive=True)
Ph, Ps, a_ = sp.symbols('Phi Psi a')
Pt, Px, St, Sx = sp.symbols('Phi_t Phi_x Psi_t Psi_x')
B = -2 * m0 * (sig * (sp.cos((Ph + Ps) / 2) - 1) + (sp.cos((Ph - Ps) / 2) - 1) / sig)
# turned conditions (Phi before the line): Phi_t = Psi_x - i B_Phi, Psi_t = Phi_x + i B_Psi
turned = {Pt: Sx - sp.I * sp.diff(B, Ph), St: Px + sp.I * sp.diff(B, Ps)}
bt = sp.solve([(Pt + Px) - (St + Sx) - 2 * a_ * m0 * sp.sin((Ph + Ps) / 2),
               (Pt - Px) + (St - Sx) + 2 * m0 / a_ * sp.sin((Ph - Ps) / 2)], [Pt, St], dict=True)[0]
diffs = [sp.simplify(((bt[k] - turned[k]).subs(a_, -sp.I * sig)).rewrite(sp.exp)) for k in (Pt, St)]
report("the defect turned through 90 degrees imposes Baecklund's B_a with a = -i sigma", diffs == [0, 0])

lam = math.sqrt(5) + 0.3
q = -cmath.exp(-1j * math.pi * lam)
THS = np.array([-0.5, 0.1, 0.6])
N = len(THS)
YS = np.exp(lam * THS)


def site(Y, K, inverse=True, s=0.37):
    ns = np.arange(-K, K + 1)
    d = 2 * K + 1
    A = np.zeros((2, 2, d, d), dtype=complex)
    if inverse:
        den = 1 + Y ** 2
        A[0, 0] = np.diag(q ** ns.astype(complex)) / den
        A[1, 1] = np.diag(q ** (-ns - s).astype(complex)) / den
    else:
        A[0, 0] = np.diag(q ** (-ns).astype(complex))
        A[1, 1] = np.diag(q ** (ns + s).astype(complex))
    for i, n in enumerate(ns):
        if i + 1 < d:
            A[1, 0][i + 1, i] = (1j * Y * q ** (-n - 1 - s) / (1 + Y ** 2)) if inverse else -1j * Y * q ** (-n)
        if i >= 1:
            A[0, 1][i - 1, i] = (1j * Y * q ** (n - 1) / (1 + Y ** 2)) if inverse else -1j * Y * q ** (n + s)
    return A


def apply_sites(order, inverse=True, yd=cmath.exp(lam * (0.21 + 0.05j))):
    """<0| A_{order[-1]} ... A_{order[0]} |0> on V^3, slot j = soliton j."""
    K = N + 1
    d = 2 * K + 1
    out = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for a in range(2 ** N):
        ain = [(a >> (N - 1 - k)) & 1 for k in range(N)]
        vec = {tuple(ain): np.eye(d, dtype=complex)[:, K]}
        for j in order:
            A = site(YS[j] / yd, K, inverse)
            new = {}
            for conf, w in vec.items():
                for b2 in (0, 1):
                    blk = A[b2, conf[j]]
                    if blk.any():
                        c2 = list(conf)
                        c2[j] = b2
                        new[tuple(c2)] = new.get(tuple(c2), 0) + blk @ w
            vec = new
        for conf, w in vec.items():
            out[int(''.join(map(str, conf)), 2), a] = w[K]
    return out


def t_op(y, om):
    T = np.eye(2 ** (N + 1), dtype=complex)
    for j in range(N):
        z = y / YS[j]
        al = q * (z ** 2 - 1) / (q ** 2 * z ** 2 - 1)
        be = z * (q ** 2 - 1) / (q ** 2 * z ** 2 - 1)
        Op = np.zeros_like(T)
        for v in range(2 ** (N + 1)):
            bits = [(v >> (N - k)) & 1 for k in range(N + 1)]  # bits[0] = auxiliary
            aa, bb = bits[0], bits[j + 1]
            if aa == bb:
                Op[v, v] += 1
            else:
                Op[v, v] += al
                b2 = list(bits)
                b2[0], b2[j + 1] = bb, aa
                Op[int(''.join(map(str, b2)), 2), v] += be
        T = Op @ T
    h = 2 ** N
    return om * T[:h, :h] + T[h:, h:] / om


chg = np.array([sum(1 if ((v >> (N - 1 - k)) & 1) == 0 else -1 for k in range(N)) for v in range(2 ** N)])
D_drawn = apply_sites([0, 1, 2])
D_rev = apply_sites([2, 1, 0])
D_noinv = apply_sites([0, 1, 2], inverse=False)
yth = cmath.exp(lam * (0.4 + 0.3j))
good, bad_rev, bad_noinv = 0.0, np.inf, np.inf
for m in range(-N, N + 1, 2):
    idx = np.where(chg == m)[0]
    t = t_op(yth, q ** (-m / 2))[np.ix_(idx, idx)]

    def comm(Dx):
        Ds = Dx[np.ix_(idx, idx)]
        return np.abs(Ds @ t - t @ Ds).max() / (np.abs(Ds).max() * np.abs(t).max())
    good = max(good, comm(D_drawn))
    if len(idx) > 1:
        bad_rev, bad_noinv = min(bad_rev, comm(D_rev)), min(bad_noinv, comm(D_noinv))
report(f"the drawn order with inverses commutes with t_omega in every sector ({good:.0e})", good < 1e-13)
report(f"the reversed order ({bad_rev:.1e}) and the product without inverses ({bad_noinv:.1e}) do not",
       bad_rev > 1e-3 and bad_noinv > 1e-3)

# ---------------------------------------------------------------- figure
R, b, H = 1.0, 0.26, 3.0
PHI0 = np.array([3.55, 4.55, 5.6])          # angles of the solitons at mid-height
SLOPE = np.tanh(THS) / R * 0.9              # d(angle)/dt


def pt(phi, y):
    return R * np.cos(phi), y + b * np.sin(phi)


def ellipse(ax, y, color, lw, ls_front='-', ls_back=(0, (2, 2)), alpha_back=0.6, z=2):
    ph = np.linspace(np.pi, 2 * np.pi, 100)
    ax.plot(*pt(ph, y), color=color, lw=lw, ls=ls_front, zorder=z)
    ph = np.linspace(0, np.pi, 100)
    ax.plot(*pt(ph, y), color=color, lw=lw * 0.8, ls=ls_back, alpha=alpha_back, zorder=z - 1)


def cylinder(ax):
    ax.plot([-R, -R], [0, H], color=style.GREY, lw=0.9)
    ax.plot([R, R], [0, H], color=style.GREY, lw=0.9)
    ellipse(ax, 0, style.GREY, 0.9)
    ph = np.linspace(0, 2 * np.pi, 200)
    ax.plot(*pt(ph, H), color=style.GREY, lw=0.9)


def worldlines(ax):
    ys = np.linspace(0, H, 400)
    for p0, k in zip(PHI0, SLOPE):
        ph = p0 + k * (ys - H / 2)
        x, y = pt(ph, ys)
        front = np.sin(ph) < 0
        for mask, ls, al in ((front, '-', 1.0), (~front, (0, (1.5, 1.5)), 0.45)):
            xm, ym = np.where(mask, x, np.nan), np.where(mask, y, np.nan)
            ax.plot(xm, ym, color=style.GREEN, lw=1.3, ls=ls, alpha=al, zorder=3)


fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.0, 3.3))
for ax in (axa, axb):
    cylinder(ax)
    worldlines(ax)
    ax.set_xlim(-1.9, 1.9)
    ax.set_ylim(-0.45, H + 0.45)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.annotate('', xy=(-1.55, 2.4), xytext=(-1.55, 0.6),
                arrowprops=dict(arrowstyle='->', color=style.GREY, lw=0.8))
    ax.text(-1.68, 1.5, 'time', rotation=90, ha='center', va='center', fontsize=7.5)

# (a) the defect along time
phd = 4.95
xd = R * np.cos(phd)
axa.plot([xd, xd], [b * np.sin(phd), H + b * np.sin(phd)], color=style.PURPLE, lw=2.4, zorder=4)
axa.text(xd + 0.06, H + 0.12, 'defect', color=style.PURPLE, fontsize=8, ha='left')
ellipse(axa, 1.55, style.ORANGE, 1.4, ls_front=(0, (4, 2)), z=5)
axa.text(R + 0.06, 1.55, r'$t(\theta)$', color=style.ORANGE, fontsize=8.5, va='center')
axa.text(0, -0.42, '(a) along time: one more site', ha='center', fontsize=8.5)

# (b) the defect along space
yD = 1.25
ellipse(axb, yD, style.PURPLE, 2.4, z=5)
for j, (p0, k) in enumerate(zip(PHI0, SLOPE)):
    # crossing of world line j with the defect: solve y + b sin(phi(y)) = yD + b sin(phi) at the same angle -> y = yD
    ph = p0 + k * (yD - H / 2)
    x, y = pt(ph, yD)
    if np.sin(ph) < 0:
        axb.plot([x], [y], 'o', color=style.PURPLE, ms=3.5, zorder=6)
        axb.text(x, y - 0.2, rf'$L_{j + 1}^{{-1}}$', color=style.PURPLE, fontsize=7.5, ha='center', va='top')
axb.text(R + 0.06, yD, r'$D(\theta_d)$', color=style.PURPLE, fontsize=8.5, va='center')
ellipse(axb, 2.45, style.ORANGE, 1.4, ls_front=(0, (4, 2)), z=5)
axb.text(R + 0.06, 2.45, r'$t(\theta)$', color=style.ORANGE, fontsize=8.5, va='center')
axb.annotate('', xy=(1.32, 1.5), xytext=(1.32, 2.25),
             arrowprops=dict(arrowstyle='->', color=style.ORANGE, lw=0.8))
axb.text(0, -0.42, '(b) along space: an operator', ha='center', fontsize=8.5)
fig.tight_layout(w_pad=0.5)
style.save(fig, "l8_cylinder", OUT)
raise SystemExit(0 if style.all_ok() else 1)
