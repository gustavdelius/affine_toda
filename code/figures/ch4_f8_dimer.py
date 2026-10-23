"""Figure F8 (chapter 4, sec-dimer): the PT-symmetric dimer.

H = [[a + i g, s], [s, a - i g]], J = sigma_x, Theta = sigma_x K.  Eigenvalues a +- sqrt(s^2 - g^2)
with eigenvectors x_pm = (s, lambda_pm - a - i g).  For |g| < s both are real and the Krein
norms [x_pm, x_pm] = +-2 s sqrt(s^2 - g^2) have opposite signs (blue: positive, red: negative);
at |g| = s the eigenvalues collide in a Jordan block (exceptional point); for |g| > s they form a
complex-conjugate pair of neutral eigenvectors (purple).

Writes part1-background/figures/ch4_f8_dimer.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
a, s = 0.0, 1.0
J = np.array([[0, 1], [1, 0]], complex)


def H(g):
    return np.array([[a + 1j * g, s], [s, a - 1j * g]])


# ---- checks ---------------------------------------------------------------
ok_alg = ok_krein = ok_neutral = ok_sym = True
for g in np.linspace(0, 2.0, 41):
    if abs(g - s) < 1e-9:
        continue
    h = H(g)
    ok_sym &= np.allclose(J @ h.conj().T @ J, h)                        # J H^dagger J = H
    lam = np.linalg.eigvals(h)
    r = np.sqrt(complex(s**2 - g**2))
    closed = np.array([a + r, a - r])
    # compare as sets (sorting complex numbers is unstable when the real parts agree to rounding)
    ok_alg &= all(np.min(np.abs(lam - c)) < 1e-9 for c in closed) and \
        all(np.min(np.abs(closed - e)) < 1e-9 for e in lam)
    for sign in (+1, -1):
        lm = a + sign * r
        x = np.array([s, lm - a - 1j * g])
        kn = (x.conj() @ J @ x)                                            # [x, x] = <x, J x>
        if g < s:
            ok_krein &= abs(kn - sign * 2 * s * np.sqrt(s**2 - g**2)) < 1e-9
        else:
            ok_neutral &= abs(kn) < 1e-9                                  # neutral when g > s
report("J H^dagger J = H for every g", ok_sym)
report("eigenvalues equal a +- sqrt(s^2 - g^2)", ok_alg)
report("Krein norms are +-2 s sqrt(s^2-g^2) of opposite sign for g<s", ok_krein)
report("both eigenvectors are neutral for g>s", ok_neutral)
h = H(s) - a * np.eye(2)
report("g=s: H-a is nilpotent and nonzero (Jordan block of size 2)",
       np.allclose(h @ h, 0) and np.linalg.norm(h) > 0.5)

# ---- figure ---------------------------------------------------------------
g = np.linspace(0, 2.0, 801)
r = np.sqrt((s**2 - g**2).astype(complex))
lp, lm = a + r, a - r
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.8, 2.9))
un = g < s
bk = g > s
for ax, fn, ylab in ((ax1, np.real, r"$\mathrm{Re}\,\lambda$"), (ax2, np.imag, r"$\mathrm{Im}\,\lambda$")):
    ax.plot(g[un], fn(lp[un]), color=style.BLUE, lw=2.4)
    ax.plot(g[un], fn(lm[un]), color=style.RED, lw=2.4)
    ax.plot(g[bk], fn(lp[bk]), color=style.PURPLE, lw=2.4)
    ax.plot(g[bk], fn(lm[bk]), color=style.PURPLE, lw=2.4)
    ax.axvline(s, color=style.GREY, ls="--", lw=1.0)
    ax.axvspan(s, 2.0, color=style.PURPLE, alpha=0.08, lw=0)
    ax.set_xlabel(r"$g/s$")
    ax.set_ylabel(ylab)
    ax.set_xlim(0, 2.0)
ax1.plot([s], [a], "o", color=style.GREY, ms=7, zorder=5)
ax1.text(0.06, 0.42, r"$\lambda_+$:  $[x_+,x_+]>0$", color=style.BLUE)
ax1.text(0.06, -0.50, r"$\lambda_-$:  $[x_-,x_-]<0$", color=style.RED)
ax1.text(1.5, 0.5, "complex pair,\nneutral\neigenvectors", color=style.PURPLE, ha="center", fontsize=8)
ax1.text(1.03, -0.12, "exceptional\npoint", color=style.GREY, fontsize=7.5, va="top")
ax1.set_title(r"(a) real parts", fontsize=9, color=style.GREY)
ax2.text(0.04, 0.07, r"unbroken: real", color=style.GREY, fontsize=8)
ax2.set_ylim(-1.15, 1.15)
ax2.set_title(r"(b) imaginary parts", fontsize=9, color=style.GREY)
fig.tight_layout()
style.save(fig, "ch4_f8_dimer", style.OUT_PART1)
raise SystemExit(0 if style.all_ok() else 1)
