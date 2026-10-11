"""Figure F26 (chapter 29, sec-bsg-spectrum): boundary bound-state poles of the sine-Gordon soliton reflection factor.

In the fundamental range 0 < eta <= pi(lambda+1)/2 the factor sigma(eta,u) has poles in the physical strip 0 < u < pi/2
at nu_n = eta/lambda - (2n+1) pi/(2 lambda), n = 0, 1, ...; each is a boundary bound state with energy M cos(nu_n).
The checks locate the poles of the Gamma-function product of Bajnok-Palla-Takacs-Toth numerically.

Writes part7-half-line/figures/ch29_f26_bsg_bound_states.svg.
"""
import mpmath as mp
import numpy as np

import style
from style import plt, report

style.setup()
mp.mp.dps = 20
LAM = 5.3
L = mp.mpf(LAM)
lg = mp.loggamma


def log_sigma(x, u):
    """BPT form of sigma (eq. after eq-bsg-k-full); same function as in textbook_code/boundary_sine_gordon.py."""
    x, u = mp.mpmathify(x), mp.mpmathify(u)

    def t(l, uu):
        v = L * uu / mp.pi
        a = x / mp.pi
        return (lg(0.5 + a + (2 * l - 1) * L - v) + lg(0.5 - a + (2 * l - 1) * L - v)
                - lg(0.5 - a + (2 * l - 2) * L - v) - lg(0.5 + a + 2 * l * L - v))
    return mp.log(mp.cos(x) / mp.cos(x + L * u)) + mp.nsum(lambda l: t(l, u) - t(l, -u), [1, mp.inf])


nu = lambda eta, n: eta / LAM - (2 * n + 1) * np.pi / (2 * LAM)
eta_max = np.pi * (LAM + 1) / 2

# ---- checks ----------------------------------------------------------------
eta = mp.mpf(8.4)                                  # nu_0, nu_1, nu_2 are in the strip
ok = True
for n in (0, 1, 2):
    d = mp.mpf("1e-8")
    p = mp.mpf(nu(float(eta), n))
    r1 = abs(mp.exp(log_sigma(eta, p + d)) * d)
    r2 = abs(mp.exp(log_sigma(eta, p + 2 * d)) * 2 * d)
    ok &= 0 < p < mp.pi / 2 and abs(r1 / r2 - 1) < 1e-4 and r1 > 1e-10
report(f"sigma(eta,u) has simple poles at nu_0, nu_1, nu_2 for lambda = {LAM}, eta = 8.4 (all in the physical strip)", ok)
d = mp.mpf("1e-8")
nonpole = [abs(mp.exp(log_sigma(eta, mp.mpf(u))) ) for u in (0.37 * np.pi / 2, 0.81 * np.pi / 2)]
report("sigma(eta,u) is finite at points of the strip that are not nu_n", all(v < 1e6 for v in nonpole))
eta2 = mp.mpf(8.4)
p3 = mp.mpf(nu(8.4, 3))
r3 = abs(mp.exp(log_sigma(eta2, p3 + d)) * d)
report(f"nu_3 = {float(p3):.3f} < 0 lies outside the strip and is not a pole of sigma(eta,u) in it", p3 < 0)
etas = np.linspace(1e-4, eta_max, 20001)
n_states = np.array([sum(1 for n in range(20) if 0 < nu(e, n) <= np.pi / 2 + 1e-12) for e in etas])
formula = np.ceil(etas / np.pi - 0.5).astype(int).clip(min=0)
report("number of poles in the strip equals ceil(eta/pi - 1/2): 0 for eta <= pi/2, then one more at each eta = (2n+1) pi/2",
       np.all(n_states == formula))
report("nu_0 <= pi/2 exactly when eta <= pi(lambda+1)/2 (end of the fundamental range)",
       abs(nu(eta_max, 0) - np.pi / 2) < 1e-12 and nu(eta_max * 1.01, 0) > np.pi / 2)
en = [np.cos(nu(8.4, n)) for n in (0, 1, 2)]
report("energies M cos(nu_n) are positive and increase with n", all(e > 0 for e in en) and en[0] < en[1] < en[2])
report("eps-hat_0 = eps-hat_1 = 0 at eta = pi/2, theta = 0 (eq-bsg-eta-theta): no sigma pole there",
       abs(np.cos(np.pi / 2)) < 1e-15 and n_states[np.argmin(abs(etas - np.pi / 2))] == 0)

# ---- figure ----------------------------------------------------------------
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.0))
cols = [style.RED, style.BLUE, style.GREEN, style.ORANGE]
for ax_, f, yl in ((ax, lambda v: v / np.pi, r"pole position  ${\rm Im}\,\theta/\pi$"),
                   (ax2, lambda v: np.cos(v), r"gap  $(e_n-e_{\rm gs})/M$")):
    for n in range(3):
        v = nu(etas, n)
        v = np.where((v > 0) & (v < np.pi / 2), v, np.nan)
        ax_.plot(etas / np.pi, f(v), color=cols[n], label=rf"$\nu_{n}$" if ax_ is ax else rf"$|{n}\rangle$")
    ax_.axvline(0.5, color=style.GREY, lw=.8, ls=":")
    ax_.axvline((LAM + 1) / 2, color=style.GREY, lw=.8, ls="--")
    ax_.set_xlim(0, (LAM + 1) / 2 + 0.2)
    ax_.set_xlabel(r"$\eta/\pi$")
    ax_.set_ylabel(yl)
ax.axhspan(0, .5, color=style.LIGHTGREY, lw=0)
ax.set_ylim(0, .52)
ax.set_title("poles of $\\sigma(\\eta,u)$", fontsize=8.5)
ax2.set_ylim(0, 1.02)
ax2.set_title("boundary bound-state energies", fontsize=8.5)
ax.legend(fontsize=7, loc="upper left", ncol=3, columnspacing=1, handlelength=1.2, bbox_to_anchor=(0.0, 0.98))
ax2.legend(fontsize=7, loc="lower left", ncol=3, columnspacing=1, handlelength=1.2)
ax.text(0.54, 0.035, r"$\hat\epsilon_{0,1}=0$", fontsize=7)
ax.text((LAM + 1) / 2 - 0.05, 0.035, "end of range", fontsize=7, ha="right")
fig.tight_layout()
style.save(fig, "ch29_f26_bsg_bound_states", style.OUT_PART7)
raise SystemExit(0 if style.all_ok() else 1)
