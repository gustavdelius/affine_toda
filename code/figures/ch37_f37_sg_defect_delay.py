"""Figure F37 (chapter 37, sec-sg-defect-solitons, eq-sg-delay, prp-sg-delay): a sine-Gordon soliton through a
type-I_- jump defect.

Left: the shift ln|z| of the soliton's position constant against theta - theta_d, with
z_{I-}(theta) = (e^{-theta} + sigma)/(e^{-theta} - sigma) = -coth((theta - theta_d)/2), sigma = e^{-theta_d}.
The soliton keeps its character for theta < theta_d (z > 0) and emerges as an antisoliton for theta > theta_d (z < 0);
ln|z| diverges at theta = theta_d, where the soliton is absorbed.  Right: world lines x = v t (incoming, x < 0) and
x = v t + Delta x (outgoing, x > 0), v = tanh theta, Delta x = -ln|z|/(m0 cosh theta), Delta t = ln|z|/(m0 sinh theta),
for theta_d = 1 and theta = 0.5, 1.6 (m0 = 1).

Checks: the two forms of z agree; |z| > 1 for real theta != theta_d and z < 0 iff theta > theta_d; ln|z| ~ ln(2/|u|) at
u = theta - theta_d -> 0; the soliton on x < 0 and the soliton with E -> zE on x > 0 satisfy both defect conditions
(defect_common.defect_condition_residual, as in defect_classical_sg.py check 4) for the plotted rapidities and a
complex position constant; for the real kink the centre (field = pi/beta_SG) of the outgoing soliton is at
x = v t + Delta x, and for theta > theta_d the outgoing field is an antisoliton; at theta = theta_d the configuration
(soliton, 2 pi/beta_SG) satisfies the defect conditions (absorption).

Writes part8-defects/figures/ch37_f37_sg_defect_delay.svg.
"""
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.optimize import brentq

import style
from style import plt, report

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from defect_common import defect_condition_residual  # noqa: E402

style.setup()
M0, BSG = 1.0, 1.0                     # m0 = 2m, beta_SG = sqrt2 beta (book normalization of defect_common)
M_BOOK, BETA_BOOK = M0 / 2, BSG / np.sqrt(2)
THD = 1.0
THS = (0.5, 1.6)


def z_closed(th, thd):
    return -1 / np.tanh((th - thd) / 2)


def z_sigma(th, thd):
    s = np.exp(-thd)
    return (np.exp(-th) + s) / (np.exp(-th) - s)


def dx(th, thd):
    return -np.log(abs(z_closed(th, thd))) / (M0 * np.cosh(th))


# ---- checks ----------------------------------------------------------------
u = np.concatenate([np.linspace(-6, -1e-3, 400), np.linspace(1e-3, 6, 400)])
for thd in (THD, -0.7):
    th = u + thd
    report(f"z_I- = (e^-theta + sigma)/(e^-theta - sigma) = -coth((theta - theta_d)/2) (eq-sg-delay), theta_d = {thd}",
           np.allclose(z_sigma(th, thd), z_closed(th, thd), rtol=1e-10))
zz = z_closed(u + THD, THD)
report("|z| > 1 for real theta != theta_d, and z < 0 exactly for theta > theta_d (conversion, prp-sg-delay)",
       np.all(np.abs(zz) > 1) and np.all((zz < 0) == (u > 0)))
us = np.array([1e-3, 1e-4, -1e-5])
report("ln|z| = ln(2/|theta - theta_d|) + O(u^2) as theta -> theta_d (absorption)",
       np.allclose(np.log(np.abs(z_closed(us + THD, THD))), np.log(2 / np.abs(us)), atol=1e-6))

mp.mp.dps = 30


def sol(x, t, th, xi, z=1):
    E = z * mp.e ** (M0 * (x * mp.cosh(th) - t * mp.sinh(th)) + xi)
    return 2j / BSG * mp.log((1 - E) / (1 + E))


def cond_res(fphi, fpsi, t0, thd):
    p, q = fphi(0, t0), fpsi(0, t0)
    px, qx = mp.diff(lambda x: fphi(x, t0), 0), mp.diff(lambda x: fpsi(x, t0), 0)
    pt, qt = mp.diff(lambda t: fphi(0, t), t0), mp.diff(lambda t: fpsi(0, t), t0)
    c = lambda v: np.array([complex(v)])
    return defect_condition_residual(c(p), c(q), c(px), c(qx), c(pt), c(qt), float(np.exp(-thd)), 1,
                                     beta=BETA_BOOK, m=M_BOOK, coupling='imag', typ=1)


worst = 0.0
for th in THS:
    z = mp.mpf(z_closed(th, THD))
    for xi in (mp.mpc('0.2', '0.9'), mp.mpc(0, mp.pi / 2)):
        for t0 in (-0.6, 0.0, 0.4, 1.3):
            worst = max(worst, cond_res(lambda x, t: sol(x, t, th, xi), lambda x, t: sol(x, t, th, xi, z), t0, THD))
report(f"soliton on x < 0 and E -> z E on x > 0 satisfy both type-I_- defect conditions, theta = {THS}, theta_d = {THD}, "
       f"two position constants, four times [res {worst:.0e}]", worst < 1e-10)

ok = True
XI_R = mp.mpc(0, mp.pi / 2)                           # real kink: E = i e^X, phi = (4/beta_SG) arctan(z e^X)


def kink(x, t, th, z=1.0):
    return 4 / BSG * np.arctan(z * np.exp(M0 * (x * np.cosh(th) - t * np.sinh(th))))


for th in THS:
    z = z_closed(th, THD)
    v = np.tanh(th)
    for t0 in (4.0, 6.0):
        x0 = v * t0 + dx(th, THD)                     # predicted centre of the outgoing soliton
        for x in (x0 - 1.0, x0, x0 + 0.7):            # arctan form = Hirota form with xi = i pi/2 (away from the cut)
            ok &= abs(complex(sol(x, t0, th, XI_R, mp.mpf(z))) - kink(x, t0, th, z)) < 1e-12
        jump_in = kink(x0 + 60, t0, th) - kink(x0 - 60, t0, th)
        jump_out = kink(x0 + 60, t0, th, z) - kink(x0 - 60, t0, th, z)
        xs = brentq(lambda x: kink(x, t0, th, z) - jump_out / 2, x0 - 3, x0 + 3, xtol=1e-13)
        ok &= abs(xs - x0) < 1e-8 and abs(jump_in - 2 * np.pi / BSG) < 1e-12
        ok &= abs(jump_out - (1 if th < THD else -1) * jump_in) < 1e-12
report("real kink: the outgoing centre is at x = v t + Delta x with Delta x = -ln|z|/(m0 cosh theta); the outgoing field "
       "jumps by 2 pi/beta_SG with the incoming sign for theta < theta_d and the opposite sign for theta > theta_d", ok)

worst = 0.0
for t0 in (-0.5, 0.3, 1.0):
    xi = mp.mpc('0.1', '0.7')
    worst = max(worst, cond_res(lambda x, t: sol(x, t, THD, xi), lambda x, t: 2 * mp.pi / BSG, t0, THD))
report(f"theta = theta_d: (soliton, 2 pi/beta_SG) satisfies the defect conditions (absorption) [res {worst:.0e}]",
       worst < 1e-12)

# ---- figure ----------------------------------------------------------------
fig = plt.figure(figsize=(7.0, 2.75))
gs = fig.add_gridspec(1, 3, width_ratios=[1.55, 1, 1], wspace=0.32)
ax = fig.add_subplot(gs[0])
ul, ur = np.linspace(-4, -0.012, 400), np.linspace(0.012, 4, 400)
ax.plot(ul, np.log(np.abs(z_closed(ul, 0))), color=style.BLUE)
ax.plot(ur, np.log(np.abs(z_closed(ur, 0))), color=style.RED)
ax.axvline(0, color=style.GREY, lw=0.8, ls=":")
ax.set_xlim(-4, 4)
ax.set_ylim(0, 4.2)
ax.set_xlabel(r"$\theta-\theta_d$")
ax.set_ylabel(r"$\ln|z|=m_0\cosh\theta\,|\Delta x|$")
ax.text(-2.5, 2.0, "soliton\nstays soliton\n" + r"$z>0$", color=style.BLUE, ha="center", fontsize=7.5)
ax.text(2.5, 2.0, "soliton becomes\nantisoliton\n" + r"$z<0$", color=style.RED, ha="center", fontsize=7.5)
ax.text(0.12, 3.95, r"absorption at $\theta=\theta_d$", color=style.GREY, ha="left", va="top", fontsize=7.5)
for th, mk, lab in zip(THS, ("o", "s"), ("(b)", "(c)")):
    ax.plot([th - THD], [np.log(abs(z_closed(th, THD)))], mk, color=style.ORANGE, ms=4, zorder=5)
    ax.text(th - THD + (-0.15 if th < THD else 0.15), np.log(abs(z_closed(th, THD))) + 0.12, lab, color=style.ORANGE,
            ha="right" if th < THD else "left", va="bottom", fontsize=7.5)
ax.set_title(r"(a) shift $\ln|z|$, $z=-\coth\frac{\theta-\theta_d}{2}$", fontsize=8)

T = np.array([-5.0, 7.0])
for k, (th, mk) in enumerate(zip(THS, ("o", "s"))):
    a = fig.add_subplot(gs[k + 1])
    v, d = np.tanh(th), dx(th, THD)
    dt = -d / v
    a.axvspan(-3, 0, color=style.BLUE, alpha=0.06, lw=0)
    a.axvspan(0, 3, color=style.RED, alpha=0.06, lw=0)
    a.axvline(0, color=style.PURPLE, lw=1.6)
    tin = np.array([T[0], 0.0])
    a.plot(v * tin, tin, color=style.BLUE)
    tout = np.array([dt, T[1]])
    out_col = style.BLUE if th < THD else style.RED
    a.plot(v * tout + d, tout, color=out_col)
    a.plot(v * np.array([0, T[1]]), [0, T[1]], color=style.GREY, lw=0.8, ls=":")      # undelayed continuation
    ta = min(3.2, 1.8 / v)                                                              # height of the Delta x arrow
    a.annotate("", xy=(v * ta + d, ta), xytext=(v * ta, ta),
               arrowprops=dict(arrowstyle="<->", color=style.ORANGE, lw=1.0, shrinkA=0, shrinkB=0))
    a.text(v * ta + 0.12, ta - 0.15, r"$\Delta x$", color=style.ORANGE, ha="left", va="top", fontsize=7.5)
    a.plot([0], [0], mk, color=style.ORANGE, ms=4, zorder=5)
    a.set_xlim(-3, 3)
    a.set_ylim(*T)
    a.set_xlabel(r"$x$")
    a.text(-0.5, T[0] + 0.3, r"$\phi$", ha="center", va="bottom", color=style.BLUE)
    a.text(0.5, T[0] + 0.3, r"$\psi$", ha="center", va="bottom", color=style.RED)
    kind = "soliton" if th < THD else "antisoliton"
    a.set_title(rf"({'bc'[k]}) $\theta={th}$: {kind} out" + "\n" + rf"$\Delta x={d:.2f},\ \Delta t={dt:.2f}$", fontsize=8)
    if k == 0:
        a.set_ylabel(r"$t$")
fig.subplots_adjust(left=0.08, right=0.99, bottom=0.16, top=0.86)
style.save(fig, "ch37_f37_sg_defect_delay", style.OUT_PART8)
raise SystemExit(0 if style.all_ok() else 1)
