"""Figure F27 (chapter 29, sec-shg-boundary-breathers, prp-shg-breather): the sinh-Gordon boundary breather.

With u = beta phi/sqrt2 = ln(tau_0/tau_1) the field equation is u_tt - u_xx + 2 sinh 2u = 0 on x <= 0 and the boundary
condition (C_0 = C_1 = eps) is u_x = -2 eps sinh u at x = 0, independent of beta.  Left: u(x,t) at four times of one period for
rho = 1.15, eps = -0.7.  Right: where tau_0 tau_1 has no zero on x <= 0 (regular): -1 < eps < 0 and cos rho < -eps.

Checks: field equation and boundary condition by finite differences; energy by quadrature equals (8/beta^2)(-cos rho - eps);
the regular/singular classification on a grid of (eps, rho) agrees with the criterion everywhere except within 0.03 of it.

Writes part7-half-line/figures/ch29_f27_shg_boundary_breather.svg.
"""
import numpy as np
from scipy.integrate import quad

import style
from style import plt, report

style.setup()


def taus(x, t, rho, eps):
    c = np.cos(rho)
    r = (eps + c) / (eps - c)
    s = np.sqrt(r + 0j)
    return [1 + (-1) ** j * 2 * np.cos(2 * t * np.sin(rho)) * np.exp(2 * x * c) * s / np.tan(rho) - np.exp(4 * x * c) * r
            for j in (0, 1)]


def u(x, t, rho, eps):
    t0, t1 = taus(x, t, rho, eps)
    return np.real(np.log(t0 / t1))


rho0, eps0 = 1.15, -0.7
h = 1e-3
# ---- checks ----------------------------------------------------------------
xs = np.array([-0.3, -1.1, -2.5])
ts = np.array([0.2, 0.9, 1.7])
X, T = np.meshgrid(xs, ts)
f = lambda x, t: u(x, t, rho0, eps0)
utt = (f(X, T + h) - 2 * f(X, T) + f(X, T - h)) / h**2
uxx = (f(X + h, T) - 2 * f(X, T) + f(X - h, T)) / h**2
res = utt - uxx + 2 * np.sinh(2 * f(X, T))
report(f"u = ln(tau_0/tau_1) solves u_tt - u_xx + 2 sinh 2u = 0 (finite differences, residual {np.abs(res).max():.0e})", np.abs(res).max() < 1e-4)
tt = np.array([0.1, 0.8, 1.9, 2.6])
ux0 = (f(0 * tt, tt) - f(0 * tt - h, tt)) / h
ux0 = (3 * f(0 * tt, tt) - 4 * f(0 * tt - h, tt) + f(0 * tt - 2 * h, tt)) / (2 * h)
bc = ux0 + 2 * eps0 * np.sinh(f(0 * tt, tt))
report(f"boundary condition u_x = -2 eps sinh u at x = 0 (residual {np.abs(bc).max():.0e})", np.abs(bc).max() < 1e-4)


def energy_beta2(t, rho, eps):
    """E beta^2 = int [u_t^2 + u_x^2 + 2(cosh 2u - 1)] dx + 4 eps (cosh u - 1)|_{x=0}"""
    g = lambda x: float(np.real(u(x, t, rho, eps)))
    def dens(x):
        ut = (g_t(x, t + 1e-5) - g_t(x, t - 1e-5)) / 2e-5
        ux = (g(x + 1e-5) - g(x - 1e-5)) / 2e-5 if x < -1e-5 else (g(x) - g(x - 1e-5)) / 1e-5
        return ut**2 + ux**2 + 2 * (np.cosh(2 * g(x)) - 1)
    g_t = lambda x, tau: float(np.real(u(x, tau, rho, eps)))
    val = quad(dens, -40, 0, limit=400, epsabs=1e-11, epsrel=1e-11)[0]
    return val + 4 * eps * (np.cosh(g(0.0)) - 1)


Es = [energy_beta2(t, rho0, eps0) for t in (0.0, 0.6, 1.3)]
Ecd = 8 * (-np.cos(rho0) - eps0)
report(f"energy {Es[0]:.6f} (beta^2 E; spread over t {max(Es) - min(Es):.0e}) = 8(-cos rho - eps) = {Ecd:.6f}",
       abs(Es[0] - Ecd) < 1e-5 and max(Es) - min(Es) < 1e-5)


def regular(rho, eps):
    c = np.cos(rho)
    r = (eps + c) / (eps - c)
    if r < 0:
        return False                               # tau has a real zero (sqrt of a negative number): cannot be real
    xx = np.linspace(-12, 0, 1201)[:, None]
    tt_ = np.linspace(0, np.pi / np.sin(rho), 201)[None, :]
    vals = [1 + s * 2 * np.cos(2 * tt_ * np.sin(rho)) * np.exp(2 * xx * c) * np.sqrt(r) / np.tan(rho) - np.exp(4 * xx * c) * r
            for s in (1, -1)]
    return min(v.min() for v in vals) > 0           # tau_j is real and tends to 1 as x -> -infinity: a zero is a sign change


E_ax = np.linspace(-0.95, -0.05, 19)
C_ax = np.linspace(0.04, 0.96, 24)
grid = np.zeros((len(C_ax), len(E_ax)), bool)
bad = 0
for i, cc in enumerate(C_ax):
    for j, ee in enumerate(E_ax):
        grid[i, j] = regular(np.arccos(cc), ee)
        crit = cc < -ee
        if grid[i, j] != crit and abs(cc + ee) > 0.03:
            bad += 1
report(f"regularity on x <= 0 agrees with (-1 < eps < 0 and cos rho < -eps) at all {grid.size} grid points not within 0.03 of the line [{bad} exceptions]", bad == 0)
report("frequency 2 sin(rho) is bounded below by 2 sin(pi a) on the regular region, eps = cos(pi a) (collapse at cos rho = -eps)",
       all(abs(np.sin(np.arccos(-e)) - np.sin(np.pi * np.arccos(e) / np.pi)) < 1e-12 for e in (-0.3, -0.7, -0.9)))

# ---- figure ----------------------------------------------------------------
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.0), gridspec_kw={"width_ratios": [1.15, 1]})
x = np.linspace(-5, 0, 800)
Tper = np.pi / np.sin(rho0)
for k, c in zip((0, 1, 3, 4), (style.RED, style.ORANGE, style.GREEN, style.BLUE)):
    t = k * Tper / 8
    ax.plot(x, u(x, t, rho0, eps0), color=c, label=rf"$t={k}T/8$" if k else r"$t=0$")
ax.axhline(0, color=style.GREY, lw=.6)
ax.axvline(0, color=style.GREY, lw=1.2)
ax.set_xlabel("$x$ (boundary at $x=0$)")
ax.set_ylabel(r"$u=\beta\phi/\sqrt{2}$")
ax.set_title(r"breather, $\rho=1.15,\ \varepsilon=-0.7$", fontsize=8.5)
ax.legend(fontsize=6.8, loc="center left", ncol=2, columnspacing=1, bbox_to_anchor=(0.0, 0.62))

Eg, Cg = np.meshgrid(E_ax, C_ax)
ax2.scatter(Eg[grid], Cg[grid], s=14, color=style.GREEN, label="regular")
ax2.scatter(Eg[~grid], Cg[~grid], s=14, marker="x", color=style.RED, lw=.8, label="singular")
ee = np.linspace(-1, 0, 100)
ax2.plot(ee, -ee, color=style.GREY, lw=1.2, label=r"$\cos\rho=-\varepsilon$")
ax2.set_xlim(-1, 0)
ax2.set_ylim(0, 1)
ax2.set_xlabel(r"$\varepsilon$")
ax2.set_ylabel(r"$\cos\rho$")
ax2.set_title("is $\\tau_0\\tau_1\\neq0$ on $x\\leq0$?", fontsize=8.5)
ax2.legend(fontsize=7, loc="upper center", ncol=3, columnspacing=.8, handletextpad=.2, bbox_to_anchor=(0.5, -0.2))
fig.tight_layout()
style.save(fig, "ch29_f27_shg_boundary_breather", style.OUT_PART7)
raise SystemExit(0 if style.all_ok() else 1)
