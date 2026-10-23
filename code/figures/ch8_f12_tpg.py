"""Figure F12 (chapter 8, sec-tpg): two tensor-product graphs.

(a) U_q(sl2 hat), spin 1 (x) spin 1: components L(2), L(1), L(0) of parity +,-,+; the edge rule
    rho_mu/rho_nu = <(C(nu)-C(mu))/4> gives rho_2 = 1, rho_1 = <2>, rho_0 = <2><1> (exr-spin1).
(b) U_q(d_{n+1}^(2)), spinor (x) spinor, n = 3: the chain Lambda^3 - Lambda^2 - Lambda^1 - Lambda^0 of
    U_q(b_3), with rho_{Lambda^{n-j}} = prod_{i<=j} <i/2>_{(-1)^i}; <a>_- = <a>, <a>_+ = (1+zq^{2a})/(z+q^{2a}).

Checks: C(j) = 2j(j+1) from lambda.(lambda+2 rho) with lambda = j alpha; the rule applied along the edges of
(a) reproduces <2> and <2><1>; every edge factor and every rho_mu of (a) and (b) is unitary,
rho(z) rho(1/z) = 1, since <a>(1/z) = 1/<a>(z).  The same eigenvalues of (a) are obtained from Jimbo's equations
in code/spin1_tpg.py.
Writes part2-algebra/figures/ch8_f12_tpg.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
rng = np.random.default_rng(1)


def br(a, z, q):
    return (1 - z * q ** (2 * a)) / (z - q ** (2 * a))


def brp(a, z, q):
    return (1 + z * q ** (2 * a)) / (z + q ** (2 * a))


# ---- checks ---------------------------------------------------------------------------------------------
alpha, rho_w = np.array([np.sqrt(2.0)]), np.array([np.sqrt(2.0) / 2])    # sl2: alpha^2 = 2, rho = alpha/2
C = lambda j: (j * alpha) @ (j * alpha + 2 * rho_w)
report("sl2: C(j) = 2 j (j+1) from lambda.(lambda + 2 rho)", all(np.isclose(C(j), 2 * j * (j + 1)) for j in (0, .5, 1, 2)))
ok_rule = ok_unit = True
for _ in range(20):
    z = complex(rng.normal(), rng.normal())
    q = np.exp(1j * rng.uniform(0.1, 3.0))
    # edges (mu -> nu): L(2) - L(1) - L(0); rho_mu/rho_nu = <(C(nu)-C(mu))/4>
    r2 = 1.0
    r1 = r2 / br((C(1) - C(2)) / 4, z, q)
    r0 = r1 / br((C(0) - C(1)) / 4, z, q)
    ok_rule &= np.isclose(r1, br(2, z, q)) and np.isclose(r0, br(2, z, q) * br(1, z, q))
    for rho in (lambda zz: br(2, zz, q), lambda zz: br(2, zz, q) * br(1, zz, q)):
        ok_unit &= np.isclose(rho(z) * rho(1 / z), 1)
report("(a) the edge rule gives rho_1 = <2>, rho_0 = <2><1>", ok_rule)
report("(a) rho_1 and rho_0 are unitary: rho(z) rho(1/z) = 1", ok_unit)
ok_b = True
for _ in range(20):
    z = complex(rng.normal(), rng.normal())
    q = np.exp(1j * rng.uniform(0.1, 3.0))
    facs = [br(0.5, z, q), brp(1.0, z, q), br(1.5, z, q)]              # i = 1 (-), 2 (+), 3 (-)
    facs_inv = [br(0.5, 1 / z, q), brp(1.0, 1 / z, q), br(1.5, 1 / z, q)]
    ok_b &= all(np.isclose(f * g, 1) for f, g in zip(facs, facs_inv))
report("(b) the factors <1/2>, <1>_+, <3/2> are unitary, hence all rho_{Lambda^k}", ok_b)

from math import comb
dims = [comb(7, k) for k in (3, 2, 1, 0)]                      # Lambda^k of C^7
par = "+--+"
sym = sum(d for d, p in zip(dims, par) if p == "+")
anti = sum(d for d, p in zip(dims, par) if p == "-")
report("(b) 8 x 8 = 36 + 28: Lambda^3 + Lambda^0 symmetric (35+1), Lambda^2 + Lambda^1 antisymmetric (21+7)",
       sum(dims) == 64 and sym == 36 and anti == 28)
edge_plus = [par[i] == par[i + 1] for i in range(3)]
report("(b) the <.>_+ factor sits exactly on the edge joining equal parities (twisted rule, only edge 2)",
       edge_plus == [False, True, False])

# ---- figure ---------------------------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.0, 3.0), gridspec_kw={"height_ratios": [1, 1]})
PCOL = {"+": style.BLUE, "-": style.ORANGE}


def chain(ax, names, parity, edges, rhos, title):
    n = len(names)
    xs = np.linspace(1.4, 8.6, n)
    ax.set_xlim(0, 10)
    ax.set_ylim(-1.0, 1.0)
    ax.axis("off")
    for i in range(n - 1):
        ax.plot(xs[i:i + 2], [0, 0], color=style.GREY, lw=1.4, zorder=1)
        ax.text((xs[i] + xs[i + 1]) / 2, 0.2, edges[i], ha="center", va="bottom", fontsize=9.5, color=style.GREY)
    for x, name, p, r in zip(xs, names, parity, rhos):
        ax.add_patch(plt.Circle((x, 0), 0.34, fc=PCOL[p], ec="none", zorder=3, transform=ax.transData))
        ax.text(x, 0, "$+$" if p == "+" else "$-$", ha="center", va="center", fontsize=11, color="white", fontweight="bold", zorder=4)
        ax.text(x, 0.45, name, ha="center", va="bottom", fontsize=9.5, color=style.GREY)
        ax.text(x, -0.45, r, ha="center", va="top", fontsize=9.5, color=style.GREY)
    ax.set_title(title, fontsize=9, color=style.GREY, loc="left")
    ax.set_aspect("equal")


chain(ax1, [r"$L(2)$", r"$L(1)$", r"$L(0)$"], "+-+",
      [r"$\rho_2/\rho_1=\langle-2\rangle$", r"$\rho_1/\rho_0=\langle-1\rangle$"],
      [r"$\rho_2=1$", r"$\rho_1=\langle2\rangle$", r"$\rho_0=\langle2\rangle\langle1\rangle$"],
      r"(a) $U_q(\widehat{sl}_2)$, spin $1\otimes1$")
chain(ax2, [r"$\Lambda^3$", r"$\Lambda^2$", r"$\Lambda^1$", r"$\Lambda^0$"], "+--+",
      [r"$\langle\frac{1}{2}\rangle$", r"$\langle1\rangle_+$", r"$\langle\frac{3}{2}\rangle$"],
      [r"$\rho=1$", r"$\langle\frac{1}{2}\rangle$", r"$\langle\frac{1}{2}\rangle\langle1\rangle_+$",
       r"$\langle\frac{1}{2}\rangle\langle1\rangle_+\langle\frac{3}{2}\rangle$"],
      r"(b) $U_q(d_4^{(2)})$, spinor $\otimes$ spinor, under $U_q(b_3)$")
fig.tight_layout()
style.save(fig, "ch8_f12_tpg", style.OUT_PART2)
raise SystemExit(0 if style.all_ok() else 1)
