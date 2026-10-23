"""Figure F16 (chapter 11, sec-theta-classical, and chapter 10, sec-breathers): breathers and Theta.

(a) The rapidity plane.  Theta maps the data (a, theta, xi) of a soliton to (abar, conj(theta), conj(xi)), i.e. it
    reflects the rapidity in the real axis and exchanges a with abar = h - a.  A breather is a pair (a, theta0 + i u),
    (abar, theta0 - i u): the two constituents are exchanged, so the data are Theta-invariant.
(b) The classical breather mass 2 M_a cos u as u runs over (0, pi/2): a continuous band; quantization (Parts V, VI)
    selects discrete values.

Checks: M_a cosh(theta0 + i u) + M_abar cosh(theta0 - i u) = 2 M_a cosh(theta0) cos u and the same with sinh give
real E and P; E^2 - P^2 = (2 M_a cos u)^2; the map (a, theta) -> (h - a, conj theta) sends the breather pair to itself;
M_a = M_{h-a} for a_n^(1), M_a = 2 h m_a/beta^2 with m_a = 2 m sin(pi a/h).
Writes part3-classical/figures/ch11_f16_breather_theta.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
h, beta2 = 5, 1.0
m = lambda a: 2 * np.sin(np.pi * a / h)
M = lambda a: 2 * h * m(a) / beta2
rng = np.random.default_rng(5)
ok_re = ok_inv = ok_mass = True
for _ in range(50):
    a = int(rng.integers(1, h))
    th0, u = rng.normal(), rng.uniform(0.05, np.pi / 2 - 0.05)
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
report("M_a = M_{h-a} for a_4^(1)", all(np.isclose(M(a), M(h - a)) for a in range(1, h)))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.1), gridspec_kw=dict(width_ratios=[1, 1.1]))
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
u = np.linspace(0, np.pi / 2, 200)
ax2.plot(u, np.cos(u) * 2, color=style.RED)
ax2.fill_between(u, 0, 2 * np.cos(u), color=style.RED, alpha=0.12, lw=0)
ax2.set_xlabel(r"$u$")
ax2.set_ylabel(r"$M_{\rm breather}/M_a$")
ax2.set_xticks([0, np.pi / 4, np.pi / 2])
ax2.set_xticklabels(["0", r"$\pi/4$", r"$\pi/2$"])
ax2.set_title(r"classical breather mass $2M_a\cos u$: a continuous band", fontsize=8.5)
fig.tight_layout()
style.save(fig, "ch11_f16_breather_theta", style.OUT_PART3)
raise SystemExit(0 if style.all_ok() else 1)
