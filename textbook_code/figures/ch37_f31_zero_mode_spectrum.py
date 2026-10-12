"""Figure F31 (chapter 37, sec-solvable-check): spectrum of the a_2 zero-mode reduction
H = |p|^2/2 - g sum_j e^{i alpha_j.x} on the field-space torus, in the complex energy plane, for the physical
sign of the potential and for the opposite sign.

Plane-wave basis p = c_1 alpha_1 + c_2 alpha_2, |c_i| <= P, |alpha|^2 = 4 (the normalisation of the chapter), g = 1.
The potential shifts p by alpha_1, alpha_2, alpha_0 = -alpha_1 - alpha_2, so H is real and not symmetric.

Checks: spectrum closed under complex conjugation; lowest level real with smallest real part (physical sign);
Z(T) = Tr exp(-T H) > 0 for T = 0.5, 2, 8 (physical sign); the opposite sign has Z(2) < 0 and a complex lowest
pair; the spectrum does not change with P (8 versus 12); complex levels exist among the lowest 40.
Writes part8-foundations/figures/ch37_f31_zero_mode_spectrum.svg."""
import itertools

import numpy as np
import style
from style import plt, report

style.setup()
G = 2 * np.array([[2, -1], [-1, 2]])          # Gram matrix of alpha_1, alpha_2 for |alpha|^2 = 4
SHIFTS = [(1, 0), (0, 1), (-1, -1)]


def spectrum(g, P, sign=1):
    keys = list(itertools.product(range(-P, P + 1), repeat=2))
    idx = {c: i for i, c in enumerate(keys)}
    H = np.zeros((len(keys), len(keys)))
    for c, i in idx.items():
        v = np.array(c)
        H[i, i] = 0.5 * v @ G @ v
        for s in SHIFTS:
            c2 = (c[0] + s[0], c[1] + s[1])
            if c2 in idx:
                H[idx[c2], i] -= sign * g
    ev = np.linalg.eigvals(H)
    return ev[np.argsort(ev.real, kind="stable")]


def conj_closed(ev, n, tol=1e-6):
    low = ev[:n]
    return all(np.abs(low - np.conj(z)).min() < tol for z in low[:n - 2])  # skip the cut-off edge


ev = spectrum(1.0, 12)
ev8 = spectrum(1.0, 8)
evm = spectrum(1.0, 12, sign=-1)
low = ev[:40]
report("physical sign: lowest 40 levels closed under complex conjugation", conj_closed(ev, 40))
report("physical sign: lowest level real, and no level has smaller real part",
       abs(ev[0].imag) < 1e-9 and np.all(ev.real >= ev[0].real - 1e-12))
report("physical sign: Z(T) = Tr exp(-T H) > 0 for T = 0.5, 2, 8",
       all(np.exp(-T * ev).sum().real > 0 for T in (0.5, 2, 8)))
report("physical sign: complex levels among the lowest 40", (np.abs(low.imag) > 1e-6).sum() > 0)
report("physical sign: lowest 10 levels unchanged between P = 8 and P = 12", np.abs(ev[:10] - ev8[:10]).max() < 1e-6)
Zm = np.exp(-2 * evm).sum().real
report("opposite sign: Z(2) < 0", Zm < 0)
report("opposite sign: lowest level is a complex pair", abs(evm[0].imag) > 1e-3 and abs(evm[0] - np.conj(evm[1])) < 1e-6)

fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.1), sharex=True, sharey=True)
for ax, e, title, col in ((axs[0], ev, "physical sign of the potential", style.BLUE),
                          (axs[1], evm, "opposite sign of the potential", style.RED)):
    e = e[:40]
    ax.axhline(0, color=style.LIGHTGREY, lw=0.8)
    ax.scatter(e.real, e.imag, s=16, color=col, zorder=3, label="levels")
    ax.scatter(e[:1].real, e[:1].imag, s=48, facecolor="none", edgecolor=style.GREEN, lw=1.6, zorder=4)
    ax.set_title(title, fontsize=8.5, color=style.GREY, pad=10)
    ax.set_ylim(-1.7, 1.9)
    ax.set_xlabel(r"$\mathrm{Re}\,E$")
axs[0].set_ylabel(r"$\mathrm{Im}\,E$")
axs[0].annotate("real vacuum,\nsmallest $\\mathrm{Re}\\,E$", xy=(ev[0].real, 0), xytext=(7.0, 0.95), fontsize=7.5,
                color=style.GREEN, arrowprops=dict(arrowstyle="->", color=style.GREEN, lw=0.9))
axs[1].annotate("lowest level is\na complex pair", xy=(evm[0].real, evm[0].imag), xytext=(7.0, 1.0), fontsize=7.5,
                color=style.RED, arrowprops=dict(arrowstyle="->", color=style.RED, lw=0.9))
fig.tight_layout()
style.save(fig, "ch37_f31_zero_mode_spectrum", style.OUT_PART8)
raise SystemExit(0 if style.all_ok() else 1)
