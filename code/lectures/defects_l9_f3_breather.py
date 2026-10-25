"""Lecture notes on defects, Lecture 9, figure 3 (fig-l9-breather): the classical breather and the quantum breathers.

(a) A classical breather crosses the type-I defect (m0 = 1, theta_d = 0.15). Its constituents have rapidities
theta pm i alpha, theta = 0.55, alpha = 0.7, so Y_{1,2} = exp(x cosh(theta pm i alpha) - t sinh(theta pm i alpha) + c),
c = 0.2 pm 0.4 i, and Phi = 4 arctan[(Y_1 + Y_2)/(1 + tan(alpha)^2 Y_1 Y_2)] on x < 0; on x > 0, Psi is the same with
Y_{1,2} -> z(theta pm i alpha - theta_d) Y_{1,2}, z(u) = -coth(u/2). Shown: the field, with the line of the envelope's
centre (|Y_1| = 1 on the left, |z Y_1| = 1 on the right) and its continuation without the defect (dashed).
(b) The time delay of the breather, m0 cos(alpha) sinh(theta) Delta t, against u = theta - theta_d: the classical value
ln|z(u + i alpha)| = (1/2) ln[(cosh u + cos alpha)/(cosh u - cos alpha)] (black) and the delays of the quantum breathers
B_n given by Wigner's relation, (pi/(2(lambda + 1))) d delta_n/d theta with e^{i delta_n} = T_{B_n} = prod_k
Tcl(theta - theta_d + i pi (n + 1 - 2k)/(2 lambda)), for n = 3, 6, 24 at lambda = n pi/(pi - 2 alpha), so that the
constituents of B_n sit at theta pm i alpha.
(c) The shift of the internal clock: arg z(u + i alpha) (black) and delta_{n+1} - delta_n for the same n.

Checks: Psi and Phi satisfy the defect conditions at x = 0 (mpmath, 30 digits) and the conjugate factors do not; the
drawn centre lines are where |Y_1| = 1 and |z Y_1| = 1; the quantum delays are computed both from the product formula
by differentiation and from the sum of the derivatives of its factors; their deviation from the classical delay and the
deviation of the clock shift fall as n grows.
Writes lecture-notes/defects/figures/l9_breather.svg.
"""
import cmath
import math
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects" / "figures"
I = 1j
PI = math.pi
TH, THD, AL = 0.55, 0.15, 0.7
C = 0.2 + 0.4j


def zcl(u):
    return -np.cosh(u / 2) / np.sinh(u / 2)


def Tcl(u):
    return zcl(u + I * PI / 2)


def breather_np(x, t, z1=1.0, z2=1.0):
    r1, r2 = TH + I * AL, TH - I * AL
    Y1 = z1 * np.exp(x * np.cosh(r1) - t * np.sinh(r1) + C)
    Y2 = z2 * np.exp(x * np.cosh(r2) - t * np.sinh(r2) + np.conj(C))
    return np.real(4 * np.arctan((Y1 + Y2) / (1 + math.tan(AL) ** 2 * Y1 * Y2)))


z1 = complex(zcl(TH + I * AL - THD))
z2 = complex(zcl(TH - I * AL - THD))

# checks: the defect conditions (mpmath)
mp.mp.dps = 30
m_th, m_thd, m_al = mp.mpf(TH), mp.mpf(THD), mp.mpf(AL)
m_c = mp.mpc(C.real, C.imag)
sig = mp.e ** (-m_thd)


def breather_mp(x, t, zz1=1, zz2=1):
    r1, r2 = m_th + 1j * m_al, m_th - 1j * m_al
    Y1 = zz1 * mp.e ** (x * mp.cosh(r1) - t * mp.sinh(r1) + m_c)
    Y2 = zz2 * mp.e ** (x * mp.cosh(r2) - t * mp.sinh(r2) + mp.conj(m_c))
    return 4 * mp.atan((Y1 + Y2) / (1 + mp.tan(m_al) ** 2 * Y1 * Y2))


mz1 = -mp.coth((m_th + 1j * m_al - m_thd) / 2)
mz2 = -mp.coth((m_th - 1j * m_al - m_thd) / 2)


def resid(a1, a2):
    worst = mp.mpf(0)
    for t in (mp.mpf(-4), mp.mpf('-0.5'), mp.mpf('2.5'), mp.mpf(6)):
        Ph = lambda x, tt: breather_mp(x, tt)
        Ps = lambda x, tt: breather_mp(x, tt, a1, a2)
        Pt, Px = mp.diff(lambda tt: Ph(0, tt), t), mp.diff(lambda x: Ph(x, t), 0)
        St, Sx = mp.diff(lambda tt: Ps(0, tt), t), mp.diff(lambda x: Ps(x, t), 0)
        S_, D_ = Ph(0, t) + Ps(0, t), Ph(0, t) - Ps(0, t)
        worst = max(worst, abs((Pt + Px) - (St + Sx) + 2 * sig * mp.sin(S_ / 2)),
                    abs((Pt - Px) + (St - Sx) - 2 / sig * mp.sin(D_ / 2)))
    return worst


rg, rb = resid(mz1, mz2), resid(mp.conj(mz1), mp.conj(mz2))
report(f"the plotted breather satisfies the defect conditions at x = 0 (residual {float(rg):.0e}); the conjugate "
       f"factors fail ({float(rb):.0e})", rg < 1e-20 and rb > 1e-3)

# the centre lines: x = v t - Re(c)/(cos alpha cosh theta) on the left, with Re c -> Re c + ln|z1| on the right
v = math.tanh(TH)
den = math.cos(AL) * math.cosh(TH)
xc_left = lambda t: v * t - C.real / den
xc_right = lambda t: v * t - (C.real + math.log(abs(z1))) / den
ok = True
for t in (-6.0, 3.0):
    r1 = TH + I * AL
    Y1l = abs(cmath.exp(xc_left(t) * cmath.cosh(r1) - t * cmath.sinh(r1) + C))
    Y1r = abs(z1 * cmath.exp(xc_right(t) * cmath.cosh(r1) - t * cmath.sinh(r1) + C))
    ok = ok and abs(Y1l - 1) < 1e-12 and abs(Y1r - 1) < 1e-12
report("the centre lines are |Y_1| = 1 on the left and |z Y_1| = 1 on the right", ok)

# quantum delays


def TBn(th, n, lam):
    return np.prod([Tcl(th - THD + I * PI * (n + 1 - 2 * k) / (2 * lam)) for k in range(1, n + 1)])


def dTcl_log(u):
    """d/du ln Tcl(u) = d/du ln z(u + i pi/2), with d/du ln z(u) = -1/sinh(u)"""
    return -1 / np.sinh(u + I * PI / 2)


def delay_q(th, n, lam):
    ks = np.arange(1, n + 1)
    dlog = np.sum(dTcl_log(th - THD + I * PI * (n + 1 - 2 * ks) / (2 * lam)))
    return (PI / (2 * (lam + 1))) * dlog.imag  # d delta/d theta = Im d ln T/d theta


def delay_cl(u):
    return 0.5 * np.log((np.cosh(u) + math.cos(AL)) / (np.cosh(u) - math.cos(AL)))


def clock_q(th, n, lam):
    return np.angle(TBn(th, n + 1, lam) / TBn(th, n, lam))


def clock_cl(u):
    return np.angle(zcl(u + I * AL))


NS = (3, 6, 24)
ok = True
for n in NS:
    lam = n * PI / (PI - 2 * AL)
    h = 1e-6
    for th in (0.3, 1.2):
        num = (np.angle(TBn(th + h, n, lam) / TBn(th - h, n, lam)) / (2 * h)) * PI / (2 * (lam + 1))
        ok = ok and abs(num - delay_q(th, n, lam)) < 1e-6
report("the quantum delays from the product by differentiation agree with the sum over its factors", ok)
us = np.linspace(-3, 3, 241)
devs, cdevs = [], []
for n in NS:
    lam = n * PI / (PI - 2 * AL)
    dq = np.array([delay_q(u + THD, n, lam) for u in us])
    devs.append(np.abs(dq - delay_cl(us)).max())
    cq = np.array([clock_q(u + THD, n, lam) for u in us])
    cdevs.append(np.abs(np.angle(np.exp(I * (cq - clock_cl(us))))).max())
ok = all(devs[i + 1] < devs[i] for i in range(len(NS) - 1)) and all(cdevs[i + 1] < cdevs[i] for i in range(len(NS) - 1))
report(f"the deviations from the classical delay ({', '.join(f'{d:.3f}' for d in devs)}) and clock shift "
       f"({', '.join(f'{d:.3f}' for d in cdevs)}) fall as n grows", ok and devs[-1] < 0.1 and cdevs[-1] < 0.1)

# figure
fig = plt.figure(figsize=(7.0, 2.7))
gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1, 1])
axa, axb, axc = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
xs = np.linspace(-9, 9, 361)
ts = np.linspace(-10, 12, 441)
Xg, Tg = np.meshgrid(xs, ts)
F = np.where(Xg < 0, breather_np(Xg, Tg), breather_np(Xg, Tg, z1, z2))
F = np.mod(F + PI, 2 * PI) - PI
axa.pcolormesh(Xg, Tg, F, cmap='RdBu_r', vmin=-PI, vmax=PI, shading='auto', rasterized=True)
axa.axvline(0, color=style.PURPLE, lw=2.0)
tl = np.linspace(-10, 12, 50)
xl, xr = xc_left(tl), xc_right(tl)
axa.plot(np.where(xl < 0, xl, np.nan), tl, color=style.GREY, lw=0.9)
axa.plot(np.where(xl >= 0, xl, np.nan), tl, color=style.GREY, lw=0.9, ls='--')
axa.plot(np.where(xr > 0, xr, np.nan), tl, color='k', lw=0.9)
axa.set_xlim(-9, 9)
axa.text(-4.5, 10.3, r'$\Phi$', fontsize=9, ha='center', color=style.GREY)
axa.text(6.8, -8.8, r'$\Psi$', fontsize=9, ha='center', color=style.GREY)
axa.set_ylim(-10, 12)
axa.set_xlabel(r'$m_0x$')
axa.set_ylabel(r'$m_0t$')
axa.set_title('(a) a classical breather crosses', fontsize=8.5, color=style.GREY)
cols = (style.ORANGE, style.GREEN, style.BLUE)
for n, cc in zip(NS, cols):
    lam = n * PI / (PI - 2 * AL)
    axb.plot(us, [delay_q(u + THD, n, lam) for u in us], color=cc, lw=1.2, label=f'$n={n}$')
    axc.plot(us, [clock_q(u + THD, n, lam) for u in us], color=cc, lw=1.2, label=f'$n={n}$')
axb.plot(us, delay_cl(us), color='k', lw=1.0, ls='--', label='classical')
axc.plot(us, clock_cl(us), color='k', lw=1.0, ls='--', label='classical')
axb.set_xlabel(r'$\theta-\theta_d$')
axc.set_xlabel(r'$\theta-\theta_d$')
axb.set_title(r'(b) delay $m_0\cos\alpha\,\sinh\theta\,\Delta t$', fontsize=8.5, color=style.GREY)
axc.set_title(r'(c) clock shift $\delta_{n+1}-\delta_n$', fontsize=8.5, color=style.GREY)
axb.legend(fontsize=7, loc='upper right', bbox_to_anchor=(1.08, 1.02), handlelength=1.4)
axc.set_yticks([0, PI / 2, PI])
axc.set_yticklabels(['0', r'$\frac{\pi}{2}$', r'$\pi$'])
fig.tight_layout(w_pad=1.2)
style.save(fig, "l9_breather", OUT)
raise SystemExit(0 if style.all_ok() else 1)
