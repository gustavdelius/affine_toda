"""Figure F7 (chapter 3, sec-sg-classical and sec-sg-semiclassics): kink and breather.

Units beta_SG = m_0 = 1.  (a) kink and antikink phi_pm = +-4 arctan e^x and (b) the energy
density (4 m_0^2/beta^2) sech^2 x, whose integral is M_cl = 8 m_0/beta^2.  (c) the breather
at its turning point, phi_b = 4 arctan(tan u / cosh(x sin u)), for three values of u.
(d) classical breather energy E_b = 2 M sin u against the DHN levels m_n = 2 M sin(n xi/2),
which lie on the same curve at u = n xi/2 (example xi = pi/4: three breathers).

Writes part1-background/figures/ch3_f7_kink_breather.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
x = np.linspace(-45, 45, 180001)
dx = x[1] - x[0]
kink = lambda x, s=+1: s * 4 * np.arctan(np.exp(x))
breather = lambda x, t, u: 4 * np.arctan(np.tan(u) * np.sin(t * np.cos(u)) / np.cosh(x * np.sin(u)))

# ---- checks (same equation and energies as textbook_code/sine_gordon.py) ---
h = 1e-3
xs = np.linspace(-6, 6, 13)
res_k = (kink(xs + h) - 2 * kink(xs) + kink(xs - h)) / h**2 - np.sin(kink(xs))        # phi'' = sin phi
report("kink solves the static equation phi'' = sin phi (finite-difference residual)",
       np.max(np.abs(res_k)) < 1e-4)
dens = 4 / np.cosh(x) ** 2
report("kink energy: integral of (4/beta^2) sech^2 x = 8 = M_cl", abs(np.sum(dens) * dx - 8) < 1e-4)
kx = np.gradient(kink(x), dx)
report("energy density (1/2)phi'^2 + (1 - cos phi) equals 4 sech^2 x",
       np.max(np.abs(0.5 * kx**2 + 1 - np.cos(kink(x)) - dens)) < 1e-4)
for u in (0.4, 0.9, 1.3):
    tt = np.array([0.0, 0.7, 2.0])
    X, T = np.meshgrid(np.linspace(-5, 5, 9), tt)
    phi = lambda xx, t: breather(xx, t, u)
    res = ((phi(X, T + h) - 2 * phi(X, T) + phi(X, T - h)) / h**2
           - (phi(X + h, T) - 2 * phi(X, T) + phi(X - h, T)) / h**2 + np.sin(phi(X, T)))
    report(f"breather u={u}: solves phi_tt - phi_xx + sin phi = 0", np.max(np.abs(res)) < 5e-4)
    phit = (breather(x, h, u) - breather(x, -h, u)) / (2 * h)               # phi = 0 at t = 0
    E = np.sum(0.5 * phit**2) * dx
    report(f"breather u={u}: energy at t=0 equals 2 M_cl sin u = 16 sin u", abs(E - 16 * np.sin(u)) < 2e-3)

# ---- figure ---------------------------------------------------------------
fig, axs = plt.subplots(2, 2, figsize=(6.8, 5.2))
(a, b), (c, d) = axs
xx = np.linspace(-9, 9, 800)
a.plot(xx, kink(xx), color=style.BLUE, label="kink")
a.plot(xx, 2 * np.pi - kink(xx), color=style.ORANGE, label="antikink")
a.axhline(0, color=style.GREY, lw=0.6)
a.axhline(2 * np.pi, color=style.GREY, lw=0.6)
a.set_yticks([0, np.pi, 2 * np.pi])
a.set_yticklabels([r"$0$", r"$\pi$", r"$2\pi$"])
a.set_xlabel(r"$m_0x$")
a.set_ylabel(r"$\beta_{\rm SG}\,\phi$")
a.legend(loc="center right", fontsize=8)
a.set_title("(a) kink and antikink", fontsize=9, color=style.GREY)
b.fill_between(xx, 4 / np.cosh(xx) ** 2, color=style.BLUE, alpha=0.35)
b.plot(xx, 4 / np.cosh(xx) ** 2, color=style.BLUE)
b.set_xlabel(r"$m_0x$")
b.set_ylabel(r"$\beta_{\rm SG}^2\,\mathcal{E}/m_0^2$")
b.text(1.8, 2.6, r"$\int\mathcal{E}\,dx=\frac{8m_0}{\beta_{\rm SG}^2}$", color=style.BLUE)
b.set_title("(b) energy density", fontsize=9, color=style.GREY)
for u, col in ((0.4, style.TEAL), (0.9, style.PURPLE), (1.3, style.RED)):
    tmax = (np.pi / 2) / np.cos(u)                       # sin(m_0 t cos u) = 1
    c.plot(xx, breather(xx, tmax, u), color=col, label=rf"$u={u}$")
c.set_xlabel(r"$m_0x$")
c.set_ylabel(r"$\beta_{\rm SG}\,\phi_b$")
c.legend(loc="upper right", fontsize=8)
c.set_title("(c) breathers at their turning point", fontsize=9, color=style.GREY)
uu = np.linspace(0, np.pi / 2, 300)
d.plot(uu, np.sin(uu), color=style.GREY)
xi0 = np.pi / 4
for n in (1, 2, 3):
    un = n * xi0 / 2
    d.plot([un], [np.sin(un)], "o", color=style.GREEN, ms=7)
    d.plot([0, un], [np.sin(un)] * 2, ":", color=style.GREEN, lw=1.0)
    d.text(0.02, np.sin(un) + 0.025, rf"$m_{n}$", color=style.GREEN, fontsize=9)
d.set_xticks([0, np.pi / 4, np.pi / 2])
d.set_xticklabels([r"$0$", r"$\pi/4$", r"$\pi/2$"])
d.set_xlabel(r"$u$")
d.set_ylabel(r"$E_b/2M$")
d.text(0.55, 0.20, r"classical: $2M\sin u$" "\n" r"DHN: $2M\sin\frac{n\xi}{2}$," "\n" r"$\xi=\pi/4$", fontsize=8)
d.set_title("(d) continuum versus tower", fontsize=9, color=style.GREY)
fig.tight_layout()
style.save(fig, "ch3_f7_kink_breather", style.OUT_PART1)
raise SystemExit(0 if style.all_ok() else 1)
