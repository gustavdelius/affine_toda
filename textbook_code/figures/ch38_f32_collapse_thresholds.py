"""Figure F32 (chapter 38, sec-collapse-thresholds, sec-threshold-sharp): where the Coulomb-gas continuum limit holds.
Left: the first collapse threshold beta_UV^2 of several affine algebras on the axis of beta^2/pi, with the regions
of Theorem below-2pi (nearest-neighbour bound), Theorem window (cluster bound) and Proposition uv-divergence
(counterterms needed).  Right: the factor of the continued bulk free energy that carries the poles, for h = 2, 3, 6;
its poles sit at the neutral thresholds beta^2 = 4 pi (1 - 1/(k h)) and none at 2 pi for h >= 3.

Checks: beta_UV^2 from the definition (minimum over clusters of up to 7 simple affine roots) equals 2, 8/3, 3, 10/3,
4, 4 (units of pi) for a_1, a_2, a_3, d_4, c_2, g_2, and 4(1-1/h) for a_4, a_5; the null vector sum n_j alpha_j
vanishes; the cluster formula gives 4 pi / |alpha_i.alpha_j| for a pair; beta_UV^2 > 2 pi except for a_1; the
numerical poles of the Gamma-function factor are at kappa = 1 - 1/(k h), the h = 2 poles at odd xi.
Writes part8-foundations/figures/ch38_f32_collapse_thresholds.svg."""
import itertools

import numpy as np
from scipy.special import gamma as Gam

import style
from style import plt, report

style.setup()


def cyclic(n):
    """Gram matrix of the affine simple roots of a_n^(1) (n >= 2), |alpha|^2 = 2."""
    h = n + 1
    G = 2 * np.eye(h)
    for j in range(h):
        G[j, (j + 1) % h] -= 1
        G[(j + 1) % h, j] -= 1
    return G, np.ones(h)


def chain(lengths, links):
    G = np.diag(lengths).astype(float)
    for i, j, v in links:
        G[i, j] = G[j, i] = v
    return G


ALG = {
    "a_1": (np.array([[2, -2], [-2, 2.0]]), np.array([1, 1.0])),
    "a_2": cyclic(2), "a_3": cyclic(3), "a_4": cyclic(4), "a_5": cyclic(5),
    "d_4": (chain([2] * 5, [(0, 2, -1), (1, 2, -1), (2, 3, -1), (2, 4, -1)]), np.array([1, 1, 2, 1, 1.0])),
    "c_2": (chain([2, 1, 2], [(0, 1, -1), (1, 2, -1)]), np.array([1, 2, 1.0])),
    "g_2": (chain([2, 2, 2 / 3], [(0, 1, -1), (1, 2, -1)]), np.array([1, 2, 3.0])),
}


def beta_uv_over_pi(G, size_max=7):
    best = np.inf
    for size in range(2, size_max + 1):
        for S in itertools.combinations_with_replacement(range(len(G)), size):
            idx = list(S)
            den = sum(G[i, i] for i in idx) - G[np.ix_(idx, idx)].sum()      # sum |a|^2 - |sum a|^2
            if den > 1e-12:
                best = min(best, 8 * (size - 1) / den)
    return best


uv = {name: beta_uv_over_pi(G) for name, (G, n) in ALG.items()}
for name, (G, n) in ALG.items():
    report(f"{name}: sum n_j alpha_j = 0", abs(n @ G @ n) < 1e-12)
want = {"a_1": 2, "a_2": 8 / 3, "a_3": 3, "d_4": 10 / 3, "c_2": 4, "g_2": 4}
report("beta_UV^2/pi from the cluster definition = 2, 8/3, 3, 10/3, 4, 4 for a_1, a_2, a_3, d_4, c_2, g_2",
       all(abs(uv[k] - v) < 1e-9 for k, v in want.items()))
report("beta_UV^2 = 4 pi (1 - 1/h) for a_4 (h = 5) and a_5 (h = 6)",
       abs(uv["a_4"] - 4 * (1 - 1 / 5)) < 1e-9 and abs(uv["a_5"] - 4 * (1 - 1 / 6)) < 1e-9)
Gc = ALG["c_2"][0]
report("pair rule: 8 pi/(-2 alpha_i.alpha_j) = 4 pi/|alpha_i.alpha_j| = 4 pi for the adjacent roots of c_2",
       abs(8 / (-2 * Gc[0, 1]) - 4) < 1e-9 and abs(4 / abs(Gc[0, 1]) - 4) < 1e-9)
report("beta_UV^2 > 2 pi for every algebra but a_1", all(uv[k] > 2 + 1e-9 for k in uv if k != "a_1"))

# bulk free-energy factor
def Gfun(k, h):
    a = 1 / (h * (1 - k)); c = k / (h * (1 - k)); B = -k / (1 - k)
    return Gam(1 - a) * Gam(1 + c) / (B * Gam(a) * Gam(1 - c))

kap = np.linspace(0.02, 0.975, 400001)
ok_poles = True
for h in (2, 3, 4, 6):
    v = np.abs(Gfun(kap, h))
    peaks = [kap[i] for i in range(1, len(v) - 1) if v[i] > v[i - 1] and v[i] > v[i + 1] and v[i] > 1e3]
    allpred = [1 - 1 / (k * h) for k in range(1, 60)]
    ok_poles &= all(min(abs(p - a) for a in allpred) < 2e-4 for p in peaks)       # every blow-up is a predicted pole
    ok_poles &= all(min(abs(p - a) for p in peaks) < 2e-4 for a in allpred[:3])   # the first three are present
report("poles of the continued bulk energy at kappa = 1 - 1/(k h), h = 2, 3, 4, 6, and no others", ok_poles)
report("h = 2: poles at xi = kappa/(1-kappa) = 1, 3, 5 (odd)", all(abs((1 - 1 / (2 * k)) / (1 / (2 * k)) - (2 * k - 1)) < 1e-9 for k in (1, 2, 3)))

# ---- figure ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.2), gridspec_kw={"width_ratios": [1.15, 1]})
rows = [("a_1", 2.0), ("a_2", uv["a_2"]), ("a_3", uv["a_3"]), ("d_4", uv["d_4"]),
        ("e_6", 4 * 11 / 12), ("e_8", 4 * 29 / 30), ("c_2, g_2", 4.0)]
for i, (lab, u) in enumerate(rows[::-1]):
    y = i
    ax1.barh(y, min(u, 2), left=0, height=0.55, color=style.GREEN, alpha=0.85)
    if u > 2:
        ax1.barh(y, u - 2, left=2, height=0.55, color=style.BLUE, alpha=0.85)
    ax1.barh(y, 4.6 - u, left=u, height=0.55, color=style.RED, alpha=0.75)
    ax1.text(-0.12, y, "$" + lab.replace(", ", ",\\,") + "$", ha="right", va="center", fontsize=8.5)
    ax1.text(u + 0.06, y, f"{u:.2f}".rstrip("0").rstrip("."), ha="left", va="center", fontsize=6.5, color="white")
ax1.axvline(2, color=style.GREY, lw=0.8, ls=":")
ax1.axvline(4, color=style.GREY, lw=0.8, ls=":")
ax1.set_xlim(0, 4.6); ax1.set_ylim(-0.6, len(rows) - 0.4)
ax1.set_yticks([])
ax1.spines["left"].set_visible(False)
ax1.set_xticks([0, 1, 2, 3, 4]); ax1.set_xlabel(r"$\beta^2/\pi$")
ax1.text(1.0, len(rows) - 0.25, "nearest-\nneighbour", ha="center", va="bottom", fontsize=7, color=style.GREEN)
ax1.text(3.0, len(rows) - 0.25, "cluster\nbound", ha="center", va="bottom", fontsize=7, color=style.BLUE)
ax1.text(4.3, len(rows) - 0.25, "diverges", ha="center", va="bottom", fontsize=7, color=style.RED)
ax1.set_title(r"first collapse threshold $\beta^2_{\rm UV}$", fontsize=8.5, pad=22)

cols = {2: style.PURPLE, 3: style.BLUE, 6: style.ORANGE}
for h, c in cols.items():
    k = np.linspace(0.1, 0.9, 40000)
    ax2.semilogy(4 * k, np.abs(Gfun(k, h)), color=c, lw=1.3, label=f"$h={h}$")
    for m in (1, 2, 3):
        ax2.axvline(4 * (1 - 1 / (m * h)), color=c, lw=0.6, ls=":")
ax2.axvline(2, color=style.GREY, lw=0.8, ls="--")
ax2.text(2.04, 0.13, r"$2\pi$", fontsize=7.5, color=style.GREY)
ax2.set_ylim(0.1, 3e3); ax2.set_xlim(0.4, 3.6)
ax2.set_xlabel(r"$\beta^2/\pi$"); ax2.set_ylabel(r"$\Gamma$-function factor of $f$")
ax2.legend(loc="upper left", fontsize=7.5, ncol=1, handlelength=1.2, columnspacing=1.0)
ax2.set_title(r"poles at $\beta^2=4\pi(1-\frac{1}{kh})$", fontsize=8.5)
fig.tight_layout()
style.save(fig, "ch38_f32_collapse_thresholds", style.OUT_PART8)
raise SystemExit(0 if style.all_ok() else 1)
