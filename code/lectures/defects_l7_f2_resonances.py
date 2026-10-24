"""Lecture notes on defects, Lecture 7, figure 2 (fig-l7-resonances): absorption becomes a ladder of resonances.

(a) The plane of u = theta - theta_d for lambda = 4.6. The scalar factor Phi is analytic and free of zeros in the strip
-pi/(2 lambda) < Im u < pi + pi/(2 lambda) (shaded). Below it, crossing gives simple poles of Phi, and so of the
transmission matrix, at u_n = -i(2n + 1) pi/(2 lambda), n = 0, 1, ... with n < lambda; the one nearest the real axis is
the resonance. Inside the strip the matrix part degenerates (det = 0) at u = +i(2n + 1) pi/(2 lambda). The classical
defect matrix degenerates at u = 0 and u = i pi.
(b) The widths Gamma_n = 2M sinh(theta_d) sin((2n+1) pi/(2 lambda)) of the poles, in units of m0 sinh(theta_d), against
b^2 = beta^2/(8 pi) = 1/(lambda + 1). As b^2 -> 0 they tend to 1, 3, 5, 7: the widths (2n + 1) hbar gamma of the
resonances of an inverted oscillator with the growth rate gamma = m0 sinh(theta_d) of the odd jump.

Checks: the plotted poles are poles of the Gamma-function product (1/eps scaling); at the plotted zeros of the
determinant the scalar factor is finite and 1 + Y^2 = 0; no further pole of Phi between them on the imaginary axis;
the plotted widths are -2 Im M cosh(theta_d + u_n); they tend to 2n + 1.
Writes lecture-notes/defects/figures/l7_resonances.svg.
"""
import cmath
import math
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.special import loggamma

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
mp.mp.dps = 30
OUT = style.ROOT / "lecture-notes" / "defects" / "figures"
LAM = 4.6


def lnPhi(u, lam, N=300, nmax=12):
    """the Gamma-function product of the scalar factor, with the factors k >= N summed by Stirling's series"""
    a = -1j * lam * u / (2 * np.pi)
    k = np.arange(N)

    def lr(aa):
        return (loggamma(k * lam + 0.25 + aa) + loggamma((k + 1) * lam + 0.75 + aa)
                - loggamma((k + 0.5) * lam + 0.25 + aa) - loggamma((k + 0.5) * lam + 0.75 + aa))

    s = np.sum(lr(a) - lr(-a))
    am = mp.mpc(a)
    tail = mp.mpc(0)
    for n in range(2, nmax + 1):
        C = mp.mpc(0)
        for sj, cj in zip((1, 1, -1, -1), (0.25, lam + 0.75, lam / 2 + 0.25, lam / 2 + 0.75)):
            C += sj * (mp.bernpoly(n + 1, cj + am) - mp.bernpoly(n + 1, cj - am))
        tail += C * (-1) ** (n + 1) / mp.mpf(n * (n + 1)) * mp.zeta(n, N) / mp.mpf(lam) ** n
    return 1j * np.pi * (1 + lam) / 4 - np.log(1 + 1j * np.exp(lam * u)) + s + complex(tail)


def Phi(u, lam):
    return np.exp(lnPhi(u, lam))


poles = [-1j * (2 * n + 1) * math.pi / (2 * LAM) for n in range(10) if n < LAM]
zeros = [1j * (2 * n + 1) * math.pi / (2 * LAM) for n in range(10) if (2 * n + 1) / (2 * LAM) < 1 + 1 / (2 * LAM)]
okp = True
for un in poles:
    r1, r2 = [eps * Phi(un + eps, LAM) for eps in (1e-4, 1e-5)]
    okp = okp and abs(r1 - r2) / abs(r2) < 2e-3 and abs(r2) > 1e-6
report(f"the {len(poles)} plotted poles u_n = -i(2n+1) pi/(2 lambda) are simple poles of the Gamma product", okp)
okz = all(abs(1 + cmath.exp(2 * LAM * z)) < 1e-12 and abs(Phi(z + 1e-7, LAM)) < 1e3 for z in zeros)
report(f"at the {len(zeros)} plotted zeros 1 + Y^2 = 0 and Phi is finite", okz)
# between the poles on the imaginary axis |Phi| stays bounded
mids = [-1j * n * math.pi / LAM for n in range(1, len(poles))]
report("no further pole of Phi between neighbouring poles on the imaginary axis",
       all(abs(Phi(m, LAM)) < 1e3 for m in mids))
thd = 0.8
okw = True
for lam in (LAM, 9.3):
    M = (lam + 1) / math.pi
    for n in range(4):
        G = 2 * (lam + 1) * math.sin((2 * n + 1) * math.pi / (2 * lam)) / math.pi
        E = M * cmath.cosh(thd - 1j * (2 * n + 1) * math.pi / (2 * lam))
        okw = okw and abs(-2 * E.imag / math.sinh(thd) - G) < 1e-12
report("the plotted widths are -2 Im M cosh(theta_d + u_n), in units of m0 sinh theta_d", okw)
lim = [2 * (1 / 1e-5) * math.sin((2 * n + 1) * math.pi / (2 * (1 / 1e-5 - 1))) / math.pi for n in range(4)]
report("as b^2 -> 0 the widths tend to 1, 3, 5, 7", all(abs(lim[n] - (2 * n + 1)) < 1e-3 for n in range(4)))

fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.0, 3.3), gridspec_kw=dict(width_ratios=[1.0, 1.35]))
lo, hi = -math.pi / (2 * LAM), math.pi + math.pi / (2 * LAM)
axa.axhspan(lo / math.pi, hi / math.pi, color=style.LIGHTGREY, lw=0)
axa.axhline(0, color=style.GREY, lw=0.6)
xq, xc = -0.3, 0.5
axa.plot([xq] * len(poles), [p.imag / math.pi for p in poles], "x", color=style.RED, ms=6, mew=1.6,
         label=r"poles of $\mathbb{T}$")
axa.plot([xq] * len(zeros), [z.imag / math.pi for z in zeros], "o", color=style.BLUE, ms=4.5, mfc="none", mew=1.3,
         label=r"$\det=0$ (rank one)")
axa.plot([xc, xc], [0, 1], "s", color=style.GREY, ms=5, mfc="none", mew=1.2, label=r"$\det L_d=0$")
axa.annotate("resonance", xy=(xq + 0.05, poles[0].imag / math.pi), xytext=(xq + 0.2, -0.42), fontsize=7,
             color=style.RED, arrowprops=dict(arrowstyle="-", color=style.RED, lw=0.6))
axa.text(0.97, 0.5, r"$\Phi$: no poles," + "\nno zeros", fontsize=7, color=style.GREY, va="center", ha="right")
axa.set_xlim(-1, 1)
axa.set_ylim(-1.25, 1.3)
axa.set_xticks([xq, xc])
axa.set_xticklabels([rf"quantum, $\lambda={LAM}$", "classical"])
axa.set_ylabel(r"$\mathrm{Im}\,(\theta-\theta_d)/\pi$")
axa.set_title(r"(a) on the line $\mathrm{Re}\,\theta=\theta_d$", fontsize=8.5, color=style.GREY)
axa.legend(loc="lower left", fontsize=6.5, bbox_to_anchor=(-0.02, -0.02))

b2 = np.linspace(1e-3, 0.5, 400)
lam_b = 1 / b2 - 1
cols = (style.RED, style.ORANGE, style.PURPLE, style.BLUE)
for n in range(4):
    ok = 2 * n + 1 < lam_b + 1e-12
    G = 2 * (lam_b + 1) * np.sin((2 * n + 1) * np.pi / (2 * lam_b)) / np.pi
    axb.plot(b2[ok], G[ok], color=cols[n], label=rf"$n={n}$")
    axb.axhline(2 * n + 1, color=cols[n], ls=":", lw=0.9)
axb.set_xlim(0, 0.5)
axb.set_ylim(0, 8.2)
axb.set_xlabel(r"$b^2=\beta^2/8\pi=1/(\lambda+1)$")
axb.set_ylabel(r"$\Gamma_n/(m_0\sinh\theta_d)$")
axb.set_title(r"(b) widths of the poles; dotted: $2n+1$", fontsize=8.5, color=style.GREY)
axb.legend(loc="upper right", fontsize=7, ncol=2)
fig.tight_layout(w_pad=2.0)
style.save(fig, "l7_resonances", OUT)
raise SystemExit(0 if style.all_ok() else 1)
