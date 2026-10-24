"""Lecture notes on defects, Lecture 6, figure 1 (fig-l6-holonomy): parallel transport around a defect.

The exact passage Phi = 4 arctan Y on x < 0, Psi = 4 arctan(z Y) on x > 0, Y = exp(m0 (x cosh theta - t sinh theta)),
z = -coth((theta - theta_d)/2), for theta = 0.6, theta_d = 1 (m0 = 1): the soliton passes, delayed.
The Lax pair of the lectures, A_+ = (i/4) d_+Phi h + (m0 lambda/2)(e^{i Phi/2} e + e^{-i Phi/2} f),
A_- = -(i/4) d_-Phi h - (m0/2 lambda)(e^{i Phi/2} f + e^{-i Phi/2} e), a_t = (A_+ + A_-)/2, a_x = (A_+ - A_-)/2,
with the defect matrix L_d = e^{i Delta h/4} + (sigma/lambda) e^{-i Sigma h/4} s at x = 0, lambda = 0.8 + 0.5 i.
(a) The charge density d_x Phi, the defect, the undelayed world line (dotted), and a rectangle from t_1 to t_2
between x = -6 and x = 6.
(b) The deviation from 1 of the holonomy H around the rectangle, crossing the defect with L_d^{-1} at t_1 and L_d at
t_2, as a function of t_2: with L_d, with L_d replaced by 1, and with sigma replaced by 1.3 sigma.

Checks: the plotted pair satisfies the defect conditions at x = 0; with L_d the holonomy is 1 to 1e-9 for every
t_2; without L_d, or with the wrong sigma, it deviates by more than 0.1 once the soliton has reached the defect.
Writes lecture-notes/defects/figures/l6_holonomy.svg.
"""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects" / "figures"

m0, thd, theta = 1.0, 1.0, 0.6
sig = math.exp(-thd)
z = -1 / math.tanh((theta - thd) / 2)
ch, sh = math.cosh(theta), math.sinh(theta)
lam = 0.8 + 0.5j
xm, xp, t1 = -6.0, 6.0, -5.0
h = np.diag([1.0, -1.0]).astype(complex)
e = np.array([[0, 1], [0, 0]], complex)
f = e.T.copy()
s = e + f


def fields(t, x, side):
    fac = 1.0 if side == 0 else z
    Y = fac * np.exp(m0 * (x * ch - t * sh))
    g = 4 * Y / (1 + Y**2)
    return 4 * np.arctan(Y), -g * m0 * sh, g * m0 * ch


def X(F):
    return np.exp(1j * F / 2) * e + np.exp(-1j * F / 2) * f


def Ym(F):
    return np.exp(1j * F / 2) * f + np.exp(-1j * F / 2) * e


def a_x(t, x, side):
    F, Ft, Fx = fields(t, x, side)
    return 1j / 4 * Ft * h + m0 / 4 * (lam * X(F) + Ym(F) / lam)


def a_t(t, x, side):
    F, Ft, Fx = fields(t, x, side)
    return 1j / 4 * Fx * h + m0 / 4 * (lam * X(F) - Ym(F) / lam)


def rk4(Afun, a, b, n):
    """Solve dM/ds = -A(s) M from a to b, M(a) = 1."""
    M = np.eye(2, dtype=complex)
    hs = (b - a) / n
    for k in range(n):
        s0 = a + k * hs
        k1 = -Afun(s0) @ M
        k2 = -Afun(s0 + hs / 2) @ (M + hs / 2 * k1)
        k3 = -Afun(s0 + hs / 2) @ (M + hs / 2 * k2)
        k4 = -Afun(s0 + hs) @ (M + hs * k3)
        M = M + hs / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return M


def Ld(t, fac=1.0):
    P = fields(t, 0.0, 0)[0]
    Q = fields(t, 0.0, 1)[0]
    S, D = P + Q, P - Q
    return np.diag([np.exp(1j * D / 4), np.exp(-1j * D / 4)]) + fac * sig / lam * np.diag(
        [np.exp(-1j * S / 4), np.exp(1j * S / 4)]) @ s


# defect conditions at x = 0
res = 0.0
for t in np.linspace(-8, 12, 21):
    ph, pt, px = fields(t, 0.0, 0)
    ps, qt, qx = fields(t, 0.0, 1)
    Bp = m0 * (sig * math.sin((ph + ps) / 2) + math.sin((ph - ps) / 2) / sig)
    Bq = m0 * (sig * math.sin((ph + ps) / 2) - math.sin((ph - ps) / 2) / sig)
    res = max(res, abs(px - (qt - Bp)), abs(qx - (pt + Bq)))
report(f"the plotted pair satisfies the defect conditions at x = 0 (residual {res:.0e})", res < 1e-12)

nx, nt_step = 2400, 80
Tm1 = rk4(lambda xv: a_x(t1, xv, 0), xm, 0.0, nx)
Tp1 = rk4(lambda xv: a_x(t1, xv, 1), 0.0, xp, nx)
t2s = np.linspace(t1 + 0.25, 10.0, 61)
Um = np.eye(2, dtype=complex)
Up = np.eye(2, dtype=complex)
tprev = t1
dev = {"L": [], "one": [], "wrong": []}
for t2 in t2s:
    Um = rk4(lambda tv: a_t(tv, xm, 0), tprev, t2, nt_step) @ Um
    Up = rk4(lambda tv: a_t(tv, xp, 1), tprev, t2, nt_step) @ Up
    tprev = t2
    Tm2 = rk4(lambda xv: a_x(t2, xv, 0), xm, 0.0, nx)
    Tp2 = rk4(lambda xv: a_x(t2, xv, 1), 0.0, xp, nx)
    for mode, fac in (("L", 1.0), ("one", None), ("wrong", 1.3)):
        L1 = np.eye(2) if fac is None else Ld(t1, fac)
        L2 = np.eye(2) if fac is None else Ld(t2, fac)
        H = (np.linalg.inv(Um) @ np.linalg.inv(Tm2) @ L2 @ np.linalg.inv(Tp2) @ Up @ Tp1 @ np.linalg.inv(L1) @ Tm1)
        dev[mode].append(np.linalg.norm(H - np.eye(2), 2))
dev = {k: np.array(v) for k, v in dev.items()}
t_arrive = -math.log(1.0) / (m0 * sh)       # the soliton's centre (|Y| = 1) reaches x = 0 at t = 0
report(f"with L_d the holonomy is 1 for every t_2 (max deviation {dev['L'].max():.0e})", dev["L"].max() < 1e-9)
after = t2s > t_arrive
report("without L_d, or with 1.3 sigma, the holonomy deviates by more than 0.1 after the soliton reaches the defect",
       dev["one"][after].min() > 0.1 and dev["wrong"][after].min() > 0.1)
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.0), gridspec_kw=dict(width_ratios=[1.0, 1.2]))
xs = np.linspace(-8, 8, 481)
ts = np.linspace(-8, 11, 571)
TT, XX = np.meshgrid(ts, xs, indexing="ij")
dens = np.where(XX < 0, fields(TT, XX, 0)[2], fields(TT, XX, 1)[2]) / (2 * m0 * ch)
ax.pcolormesh(xs, ts, dens, cmap="RdBu_r", vmin=-1, vmax=1, shading="auto", rasterized=True)
ax.axvline(0, color=style.PURPLE, lw=1.4)
tl = np.array([0.0, 11.0])
ax.plot(tl * math.tanh(theta), tl, ls=":", color=style.GREY, lw=1.0)
t2d = 6.0
kw = dict(color=style.GREY, lw=1.1)
ax.plot([xm, xp], [t1, t1], **kw)
ax.plot([xm, xp], [t2d, t2d], **kw)
ax.plot([xm, xm], [t1, t2d], **kw)
ax.plot([xp, xp], [t1, t2d], **kw)
arr = dict(arrowstyle="-|>", color=style.GREY, lw=1.1, mutation_scale=9)
ax.annotate("", xy=(-2.6, t1), xytext=(-3.4, t1), arrowprops=arr)
ax.annotate("", xy=(3.4, t1), xytext=(2.6, t1), arrowprops=arr)
ax.annotate("", xy=(xp, 0.9), xytext=(xp, 0.1), arrowprops=arr)
ax.annotate("", xy=(2.6, t2d), xytext=(3.4, t2d), arrowprops=arr)
ax.annotate("", xy=(-3.4, t2d), xytext=(-2.6, t2d), arrowprops=arr)
ax.annotate("", xy=(xm, 0.1), xytext=(xm, 0.9), arrowprops=arr)
ax.plot([0, 0], [t1, t2d], "o", color=style.PURPLE, ms=4)
ax.text(0.35, t1 - 1.0, r"$L_d^{-1}$", color=style.PURPLE, fontsize=8.5)
ax.text(-1.3, t2d + 0.45, r"$L_d$", color=style.PURPLE, fontsize=8.5)
ax.text(-3.0, t1 + 0.4, r"$T_-$", color=style.GREY, fontsize=8)
ax.text(3.0, t1 + 0.4, r"$T_+$", color=style.GREY, fontsize=8)
ax.set_xlim(xs[0], xs[-1])
ax.set_ylim(ts[0], ts[-1])
ax.set_xlabel(r"$m_0x$")
ax.set_ylabel(r"$m_0t$")
ax.set_title("(a) a loop that crosses the defect", fontsize=8.5, color=style.GREY)

ax2.semilogy(t2s, np.maximum(dev["L"], 1e-16), color=style.PURPLE, label=r"with $L_d$")
ax2.semilogy(t2s, dev["one"], color=style.ORANGE, label=r"$L_d\to1$")
ax2.semilogy(t2s, dev["wrong"], color=style.TEAL, ls="--", label=r"$\sigma\to1.3\,\sigma$")
ax2.axvline(t_arrive, color=style.GREY, lw=0.8, ls=":")
ax2.text(t_arrive + 0.2, 3e-15, "soliton reaches\nthe defect", color=style.GREY, fontsize=7)
ax2.set_xlim(t2s[0], t2s[-1])
ax2.set_ylim(1e-16, 300)
ax2.set_xlabel(r"$m_0t_2$")
ax2.set_ylabel(r"$\Vert H-1\Vert$")
ax2.legend(loc="center right", fontsize=7.5)
ax2.set_title(r"(b) the holonomy, $\lambda=0.8+0.5i$", fontsize=8.5, color=style.GREY)
fig.tight_layout(w_pad=1.2)
style.save(fig, "l6_holonomy", OUT)
raise SystemExit(0 if style.all_ok() else 1)
