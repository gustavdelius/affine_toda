"""Figure F16 (chapter 11, sec-theta-classical; chapter 10, sec-breathers): breathers and Theta.

(a) The rapidity plane.  Theta maps the data (a, theta, xi) of a soliton to (abar, conj(theta), conj(xi)), i.e. it
    reflects the rapidity in the real axis and exchanges a with abar = h - a.  A neutral breather is a pair
    (a, theta0 + i u), (abar, theta0 - i u): the two constituents are exchanged, so the data are Theta-invariant.
(b) The regular neutral breathers of a_5^(1) (h = 6).  The breather of species a has mass 2 M_a cos u, but it is
    regular only where A_{a abar}(2iu) = 1 - sin^2(pi a/h)/cos^2 u < 0, i.e. |pi/2 - pi a/h| < u < pi/2; so its mass
    lies in (0, 2 M_a sin(pi a/h)).  For 0 < A < 1 some tau_j vanishes, whatever the phases.  Only a = h/2 = 3
    fills (0, 2 M_a).

Checks: M_a cosh(theta0 + iu) + M_abar cosh(theta0 - iu) = 2 M_a cosh(theta0) cos u and the same with sinh give real
E and P; E^2 - P^2 = (2 M_a cos u)^2; (a, theta) -> (h - a, conj theta) maps the pair to itself; M_a = M_{h-a};
A_{a,h-a}(2iu) from eq-interaction-coefficient equals 1 - sin^2(pi a/h)/cos^2 u and is negative exactly in the band;
in the band, xi_2 = conj(xi_1) + i pi gives Re tau_j >= 1 on a grid (regular); below it (0 < A < 1) a zero of tau_0
is found by root-finding; the band ends at 2 M_a sin(pi a/h), which is 2 M_a only for a = h/2.
Writes part3-classical/figures/ch11_f16_breather_theta.svg.
"""
import itertools
import warnings

import numpy as np
from scipy.optimize import fsolve

import style
from style import plt, report

style.setup()
warnings.simplefilter("ignore", RuntimeWarning)
h, beta2 = 6, 1.0
m = lambda a: 2 * np.sin(np.pi * a / h)
M = lambda a: 2 * h * m(a) / beta2
A = lambda a, b, th: (np.cosh(th) - np.cos(np.pi * (a - b) / h)) / (np.cosh(th) - np.cos(np.pi * (a + b) / h))
band = lambda a: abs(np.pi / 2 - np.pi * a / h)          # regular for band(a) < u < pi/2


def taus(sols, x, t):
    """tau_j of eq-n-soliton for a_{h-1}^(1), sols = [(a, theta, xi), ...]."""
    om = np.exp(2j * np.pi / h)
    out = np.zeros((h, np.size(x)), complex)
    for r in range(len(sols) + 1):
        for S in itertools.combinations(range(len(sols)), r):
            e = np.ones(np.size(x), complex)
            for i in S:
                a, th, xi = sols[i]
                e = e * np.exp(m(a) * (x * np.cosh(th) - t * np.sinh(th)) + xi)
            for i1, i2 in itertools.combinations(S, 2):
                e = e * A(sols[i1][0], sols[i2][0], sols[i1][1] - sols[i2][1])
            out += np.outer([om ** (j * sum(sols[i][0] for i in S)) for j in range(h)], e)
    return out


rng = np.random.default_rng(5)
ok_re = ok_inv = ok_mass = True
for _ in range(50):
    a = int(rng.integers(1, h))
    th0, u = rng.normal(), rng.uniform(band(a) + 0.02, np.pi / 2 - 0.02)
    E = M(a) * np.cosh(th0 + 1j * u) + M(h - a) * np.cosh(th0 - 1j * u)
    P = M(a) * np.sinh(th0 + 1j * u) + M(h - a) * np.sinh(th0 - 1j * u)
    ok_re &= abs(E.imag) < 1e-9 and abs(P.imag) < 1e-9 and np.isclose(E.real, 2 * M(a) * np.cosh(th0) * np.cos(u))
    ok_mass &= np.isclose(E.real**2 - P.real**2, (2 * M(a) * np.cos(u)) ** 2)
    pair = {(a, th0 + 1j * u), (h - a, th0 - 1j * u)}
    mapped = {(h - aa, np.conj(t)) for aa, t in pair}
    ok_inv &= all(any(x[0] == y[0] and abs(x[1] - y[1]) < 1e-12 for y in pair) for x in mapped)
report("breather E and P are real and E = 2 M_a cosh(theta0) cos u", ok_re)
report("E^2 - P^2 = (2 M_a cos u)^2", ok_mass)
report("Theta: (a, theta) -> (h - a, conj theta) maps the breather pair to itself", ok_inv)
report("M_a = M_{h-a} for a_5^(1)", all(np.isclose(M(a), M(h - a)) for a in range(1, h)))

ok_A = True
for a in (1, 2, 3):
    for u in np.linspace(0.03, np.pi / 2 - 0.03, 40):
        f = 1 - np.sin(np.pi * a / h) ** 2 / np.cos(u) ** 2
        ok_A &= abs(A(a, h - a, 2j * u).real / f - 1) < 1e-12 and ((f < 0) == (u > band(a)))
report("A_{a,h-a}(2iu) = 1 - sin^2(pi a/h)/cos^2 u, negative exactly for |pi/2 - pi a/h| < u < pi/2", ok_A)

ok_reg = True
for a in (1, 2, 3):
    u = (band(a) + np.pi / 2) / 2
    xg, tg = np.meshgrid(np.linspace(-12, 12, 241), np.linspace(0, 2 * np.pi / (m(a) * np.sin(u)), 121))
    T = taus([(a, 1j * u, 0.4 + 0.3j), (h - a, -1j * u, 0.4 - 0.3j + 1j * np.pi)], xg.ravel(), tg.ravel())
    ok_reg &= T.real.min() > 0.999
report("in the band, xi_2 = conj(xi_1) + i pi gives Re tau_j >= 1: a regular breather, for a = 1, 2, 3", ok_reg)

a, u = 1, 0.5 * band(1)
sols = [(a, 1j * u, 0.2 + 0.1j), (h - a, -1j * u, 0.3 + 0.9j)]
g = lambda v: [taus(sols, np.array([v[0]]), v[1])[0, 0].real, taus(sols, np.array([v[0]]), v[1])[0, 0].imag]
res = min((np.abs(g(fsolve(g, [x0, t0], xtol=1e-13))).max() for x0 in np.linspace(-4, 4, 9) for t0 in (0.5, 2, 4)))
report(f"below the band (a = 1, u = {u:.3f}, A = {A(1, 5, 2j * u).real:.3f}): tau_0 has a zero (residual {res:.1e})",
       0 < A(1, 5, 2j * u).real < 1 and res < 1e-10)
tops = [2 * np.cos(band(a)) for a in (1, 2, 3)]
report("band top 2 M_a cos(|pi/2 - pi a/h|) = 2 M_a sin(pi a/h); equal to 2 M_a only for a = h/2",
       np.allclose(tops, [2 * np.sin(np.pi * a / h) for a in (1, 2, 3)]) and np.isclose(tops[2], 2)
       and tops[0] < 2 and tops[1] < 2)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.1), gridspec_kw=dict(width_ratios=[1, 1.25]))
u0 = 0.9
ax1.axhline(0, color=style.GREY, lw=0.8)
ax1.set_xlim(-1.1, 1.6)
ax1.set_ylim(-1.5, 1.5)
ax1.scatter([0.4, 0.4], [u0, -u0], s=60, color=[style.RED, style.BLUE], zorder=3)
ax1.annotate(r"$(a,\ \theta_0+iu)$", (0.4, u0), (0.55, u0 + 0.05), fontsize=8, color=style.RED)
ax1.annotate(r"$(\bar a,\ \theta_0-iu)$", (0.4, -u0), (0.55, -u0 - 0.2), fontsize=8, color=style.BLUE)
ax1.annotate("", xy=(0.4, -u0 + 0.12), xytext=(0.4, u0 - 0.12),
             arrowprops=dict(arrowstyle="<->", color=style.GREY, lw=1.2))
ax1.text(-0.15, 0.1, r"$\Theta$", fontsize=11, color=style.GREY, ha="right")
ax1.set_xlabel(r"$\mathrm{Re}\,\theta$")
ax1.set_ylabel(r"$\mathrm{Im}\,\theta$")
ax1.set_yticks([-u0, 0, u0])
ax1.set_yticklabels([r"$-u$", "0", r"$u$"])
ax1.set_xticks([0.4])
ax1.set_xticklabels([r"$\theta_0$"])
ax1.set_title(r"$\Theta$: $(a,\theta)\mapsto(\bar a,\bar\theta)$ exchanges the pair", fontsize=8.5)

labels = {1: r"$a=1,\ \bar a=5$", 2: r"$a=2,\ \bar a=4$", 3: r"$a=\bar a=3$"}
ulab = {1: r"$u>\pi/3$", 2: r"$u>\pi/6$", 3: r"all $u$"}
for k, a in enumerate((1, 2, 3)):
    y, top = 2 - k, tops[k]
    ax2.barh(y, top, height=0.5, color=style.RED, alpha=0.75, lw=0)
    if top < 2 - 1e-9:
        ax2.barh(y, 2 - top, left=top, height=0.5, color="none", edgecolor=style.GREY, hatch="////", lw=0.6)
    ax2.text(top / 2, y, ulab[a], fontsize=7.5, color="white", ha="center", va="center")
ax2.text(1.5, 2.42, r"$\tau_j$ vanishes", fontsize=7.5, color=style.GREY, ha="center")
ax2.set_yticks([2, 1, 0])
ax2.set_yticklabels([labels[1], labels[2], labels[3]], fontsize=8)
ax2.set_ylim(-0.5, 2.75)
ax2.set_xlim(0, 2)
ax2.set_xticks([0, 1, np.sqrt(3), 2])
ax2.set_xticklabels(["0", "1", r"$\sqrt{3}$", "2"])
ax2.set_xlabel(r"breather mass $2M_a\cos u$, in units of $M_a$")
ax2.set_title(r"regular neutral breathers of $a_5^{(1)}$", fontsize=8.5)
fig.tight_layout()
style.save(fig, "ch11_f16_breather_theta", style.OUT_PART3)
raise SystemExit(0 if style.all_ok() else 1)
