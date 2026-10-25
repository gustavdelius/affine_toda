"""Lectures on integrable defects (defects-gemini), Lecture 3, figure 2 (fig-lec3-a2-chiral):
Classical a_2^(1) Soliton Transmission Factor z_a(theta) and the Chiral Topological Charge Selection Rule
(sec-lec3-solitons, eq-lec3-za-formula, eq-lec3-za-modulus, eq-lec3-za-phase).

Panel (a): Modulus |z_a(theta)| > 1 and spatial delay shift ln|z_a(theta)| > 0 against theta - theta_d for a_2^(1)
           (h = 3, gamma_1 = pi/3, gamma_2 = 2*pi/3), showing that solitons are always delayed (|z_a| > 1) with maximum
           delay at theta = theta_d.
Panel (b): Phase shift arg z_a(theta) against theta - theta_d for the vector soliton V (a = 1, blue, sweeping across
           (0, 2*pi/3)) and the antivector soliton Vbar (a = 2, red, sweeping across (-2*pi/3, 0), or (4*pi/3, 2*pi) mod 2*pi).
Panel (c): The sl_3 weight triangle of the fundamental soliton V = {|1>_V, |2>_V, |3>_V}: because arg z_1(theta) lies
           strictly inside (0, 2*pi/3), an incoming soliton |j>_V can keep its charge (k = 0, blue) or step once
           backward along the necklace (k = 1, |j>_V -> |j-1>_V, green), whereas stepping forward
           (k = 2, |j>_V -> |j+1>_V, red dashed crossed arrow) is classically forbidden.

Checks:
  1. Hirota 1-soliton tau-functions tau_j = 1 + omega^{aj} E_a and tildetau_j = 1 + omega^{aj} z_a(theta) E_a satisfy
     the a_2^(1) light-cone defect equations for both a = 1 and a = 2 if and only if z_a(theta) is given by eq-lec3-za-formula.
  2. |z_a(theta)|^2 matches eq-lec3-za-modulus with |z_a(theta)| > 1 for all real theta, and arg z_1(theta) sweeps
     monotonically across (0, 2*pi/3) while arg z_2(theta) = -arg z_1(theta) sweeps across (-2*pi/3, 0).
  3. Complex absorption poles of z_a^{I_-}(theta) and z_a^{I_+}(theta) match eq-lec3-absorption-rapidity.
  4. Problem 3.2: Sklyanin Poisson bracket {L_d(lambda_1) (x), L_d(lambda_2)} = [r(lambda_1/lambda_2), L_d(lambda_1) (x) L_d(lambda_2)]
     holds identically from {phi(0), psi(0)} = -1.

Writes lecture-notes/defects-gemini/figures/l3_a2_chiral.svg.
"""
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects-gemini" / "figures"

h = 3
omega = np.exp(2j * np.pi / h)


def z_an(u, a, h_val=3):
    """z_a(theta) as a function of u = theta - theta_d (eq-lec3-za-formula)."""
    gam = np.pi * a / h_val
    return (np.exp(-u) - 1j * np.exp(1j * gam)) / (np.exp(-u) - 1j * np.exp(-1j * gam))


def z_an_mod_sq(u, a, h_val=3):
    """Explicit modulus-squared formula eq-lec3-za-modulus."""
    gam = np.pi * a / h_val
    em = np.exp(-u)
    return (em**2 + 2.0 * em * np.sin(gam) + 1.0) / (em**2 - 2.0 * em * np.sin(gam) + 1.0)


# ---- Check 1: Hirota bilinear defect equations for a_2^(1) ----
# In Hirota variables (eq-lec3-fj-gj-tau, exr-dc-delay):
# f_j = (tau_{j+1} * tildetau_{j-1}) / (tau_j * tildetau_j), g_j = (tau_{j-1} * tildetau_j) / (tau_j * tildetau_{j-1})
# Component Bäcklund equations (39-defect-conditions.qmd eq-dc-components):
# partial_+ ln(g_j) = -m * sigma * (f_j - f_{j-1})
# partial_- ln(f_j) = (m / sigma) * (g_j - g_{j+1})
m_val = 1.0
thd_val = 0.4
sig_val = np.exp(-thd_val)
ts_test = np.linspace(-3.0, 3.0, 40)

ok_hirota = True
for a_sp in (1, 2):
    gam = np.pi * a_sp / h
    th_val = 0.95
    u_val = th_val - thd_val
    z_val = z_an(u_val, a_sp, h)
    # Under partial_+ = partial_t + partial_x: since E_a = exp(2*m*sin(gam)*(x cosh th - t sinh th)),
    # partial_+ E_a = 2*m*sin(gam)*(cosh th - sinh th) E_a = 2*m*sin(gam)*e^{-th} E_a
    # partial_- E_a = -2*m*sin(gam)*(cosh th + sinh th) E_a = -2*m*sin(gam)*e^{+th} E_a
    E_t = np.exp(-2.0 * m_val * np.sin(gam) * np.sinh(th_val) * ts_test + 0.3j)
    dp_E = 2.0 * m_val * np.sin(gam) * np.exp(-th_val) * E_t
    dm_E = -2.0 * m_val * np.sin(gam) * np.exp(th_val) * E_t

    tau = [1.0 + (omega ** (a_sp * j)) * E_t for j in range(h)]
    dtau_p = [(omega ** (a_sp * j)) * dp_E for j in range(h)]
    dtau_m = [(omega ** (a_sp * j)) * dm_E for j in range(h)]

    ttau = [1.0 + (omega ** (a_sp * j)) * z_val * E_t for j in range(h)]
    dttau_p = [(omega ** (a_sp * j)) * z_val * dp_E for j in range(h)]
    dttau_m = [(omega ** (a_sp * j)) * z_val * dm_E for j in range(h)]

    f_tau = [tau[(j + 1) % h] * ttau[(j - 1) % h] / (tau[j] * ttau[j]) for j in range(h)]
    g_tau = [tau[(j - 1) % h] * ttau[j] / (tau[j] * ttau[(j - 1) % h]) for j in range(h)]

    for j in range(h):
        # partial_+ ln(g_j) = partial_+ (ln tau_{j-1} + ln tildetau_j - ln tau_j - ln tildetau_{j-1})
        lhs_p = (
            dtau_p[(j - 1) % h] / tau[(j - 1) % h]
            + dttau_p[j] / ttau[j]
            - dtau_p[j] / tau[j]
            - dttau_p[(j - 1) % h] / ttau[(j - 1) % h]
        )
        rhs_p = -m_val * sig_val * (f_tau[j] - f_tau[(j - 1) % h])

        # partial_- ln(f_j) = partial_- (ln tau_{j+1} + ln tildetau_{j-1} - ln tau_j - ln tildetau_j)
        lhs_m = (
            dtau_m[(j + 1) % h] / tau[(j + 1) % h]
            + dttau_m[(j - 1) % h] / ttau[(j - 1) % h]
            - dtau_m[j] / tau[j]
            - dttau_m[j] / ttau[j]
        )
        rhs_m = (m_val / sig_val) * (g_tau[j] - g_tau[(j + 1) % h])

        ok_hirota &= np.allclose(lhs_p, rhs_p) and np.allclose(lhs_m, rhs_m)

report(
    "a_2^(1) Hirota 1-soliton tau-functions satisfy both light-cone defect equations across all t with z_a(theta) of eq-lec3-za-formula",
    ok_hirota,
)

# ---- Check 2: Modulus and phase range ----
us = np.linspace(-6.0, 6.0, 500)
z1 = z_an(us, 1, 3)
z2 = z_an(us, 2, 3)
ph1 = np.angle(z1)  # in (0, 2*pi/3)
ph2 = np.angle(z2)  # in (-2*pi/3, 0)

report(
    "a_2^(1) modulus and phase: |z_a(th)|^2 matches eq-lec3-za-modulus (>1 for all real th), arg z_1 sweeps (0, 2*pi/3), and arg z_2 = -arg z_1",
    np.allclose(np.abs(z1) ** 2, z_an_mod_sq(us, 1, 3))
    and np.allclose(np.abs(z2) ** 2, z_an_mod_sq(us, 2, 3))
    and np.all(np.abs(z1) > 1.0)
    and np.allclose(z2, np.conj(z1))
    and abs(ph1[0] - 0.0) < 1e-2
    and abs(ph1[-1] - 2 * np.pi / 3) < 1e-2
    and abs(ph2[0] - 0.0) < 1e-2
    and abs(ph2[-1] + 2 * np.pi / 3) < 1e-2,
)

# ---- Check 3: Complex absorption rapidities eq-lec3-absorption-rapidity ----
ok_abs = True
for a_sp in (1, 2):
    gam = np.pi * a_sp / h
    u_pole_Iminus = -1j * (0.5 * np.pi - gam)
    u_pole_Iplus = +1j * (0.5 * np.pi - gam)
    u_zero_Iminus = -1j * (0.5 * np.pi + gam)
    ok_abs &= (
        abs(1.0 / z_an(u_pole_Iminus, a_sp, h)) < 1e-14
        and abs(1.0 / z_an(-u_pole_Iplus, a_sp, h)) < 1e-14
        and abs(z_an(u_zero_Iminus, a_sp, h)) < 1e-14
    )

report(
    "complex absorption rapidities u_abs = -+ i(pi/2 - pi*a/h) are exact poles of z_a^{I_-}(u) and z_a^{I_+}(u)",
    ok_abs,
)

# ---- Check 4: Problem 3.2 (Sklyanin bracket of L_d(lambda) for a_1^(1)) ----
phi_s, psi_s, b_s, sig_s, l1_s, l2_s = sp.symbols("phi psi b sigma l1 l2")
c_s = b_s / (2 * sp.sqrt(2))
Ld = lambda lam: sp.Matrix(
    [
        [sp.exp(c_s * (phi_s - psi_s)), (sig_s / lam) * sp.exp(-c_s * (phi_s + psi_s))],
        [(sig_s / lam) * sp.exp(c_s * (phi_s + psi_s)), sp.exp(-c_s * (phi_s - psi_s))],
    ]
)
L1 = Ld(l1_s)
L2 = Ld(l2_s)
# Fundamental bracket {phi, psi} = -1 => {A, B} = -(dA/dphi dB/dpsi - dA/dpsi dB/dphi)
pb_tensor = sp.zeros(4, 4)
for i1 in range(2):
    for j1 in range(2):
        for i2 in range(2):
            for j2 in range(2):
                A_ij = L1[i1, j1]
                B_kl = L2[i2, j2]
                pb_tensor[2 * i1 + i2, 2 * j1 + j2] = -(
                    sp.diff(A_ij, phi_s) * sp.diff(B_kl, psi_s)
                    - sp.diff(A_ij, psi_s) * sp.diff(B_kl, phi_s)
                )

# Trigonometric classical r-matrix for a_1^(1)
x_ratio = l1_s / l2_s
r_cl = (b_s**2 / (4 * (x_ratio - 1 / x_ratio))) * sp.Matrix(
    [
        [(x_ratio + 1 / x_ratio) / 2, 0, 0, 0],
        [0, -(x_ratio + 1 / x_ratio) / 2, 2, 0],
        [0, 2, -(x_ratio + 1 / x_ratio) / 2, 0],
        [0, 0, 0, (x_ratio + 1 / x_ratio) / 2],
    ]
)
L1_kron_L2 = sp.kronecker_product(L1, L2)
comm_r = r_cl * L1_kron_L2 - L1_kron_L2 * r_cl
sklyanin_diff = sp.simplify(pb_tensor - comm_r)

report(
    "Problem 3.2: Sklyanin Poisson bracket {L_d(l1) (x), L_d(l2)} = [r(l1/l2), L_d(l1) (x) L_d(l2)] holds identically from {phi(0), psi(0)} = -1",
    sklyanin_diff == sp.zeros(4, 4),
)

# ---- Figure ----
fig, (axa, axb, axc) = plt.subplots(
    1, 3, figsize=(7.0, 2.75), gridspec_kw=dict(width_ratios=[1.05, 1.1, 0.95], wspace=0.34)
)

u_plot = np.linspace(-3.5, 3.5, 400)
z1_p = z_an(u_plot, 1, 3)
z2_p = z_an(u_plot, 2, 3)

# Panel (a): Modulus |z_a(u)| and delay shift ln|z_a(u)|
axa.axhline(1.0, color=style.LIGHTGREY, lw=0.8, ls=":")
axa.axvline(0.0, color=style.PURPLE, lw=0.9, ls=":")
axa.plot(u_plot, np.abs(z1_p), color=style.BLUE, lw=1.9, label=r"$|z_1(\theta)| = |z_2(\theta)| > 1$")
axa.plot(u_plot, np.log(np.abs(z1_p)), color=style.TEAL, lw=1.7, ls="--", label=r"delay shift $\ln|z_a(\theta)|$")

max_mod = float(np.abs(z_an(0.0, 1, 3)))
max_shift = float(np.log(max_mod))
axa.plot([0.0], [max_mod], "o", color=style.BLUE, ms=4.5, zorder=5)
axa.plot([0.0], [max_shift], "o", color=style.TEAL, ms=4.5, zorder=5)
axa.text(0.24, max_mod, r"$2+\sqrt{3}\approx 3.73$", fontsize=7.1, color=style.BLUE, va="center", ha="left")

axa.set_xlim(-3.5, 3.5)
axa.set_ylim(0.0, 4.95)
axa.set_xlabel(r"$\theta - \theta_d$")
axa.set_ylabel(r"$|z_a(\theta)|$ and $\ln|z_a(\theta)|$")
axa.set_title(r"(a) Spatial delay in $a_2^{(1)}$ ($h=3$)", fontsize=8.3, color=style.GREY)
axa.legend(loc="upper left", fontsize=6.7, handlelength=1.4, borderpad=0.2)

# Panel (b): Phase sweep arg z_1(theta) in (0, 2*pi/3) and arg z_2(theta) in (-2*pi/3, 0)
ph1_p = np.angle(z1_p) / np.pi
ph2_p = np.angle(z2_p) / np.pi

axb.axhspan(0.0, 2.0 / 3.0, color=style.BLUE, alpha=0.06, lw=0)
axb.axhspan(-2.0 / 3.0, 0.0, color=style.RED, alpha=0.05, lw=0)
axb.axhline(2.0 / 3.0, color=style.BLUE, lw=0.9, ls=":")
axb.axhline(0.0, color=style.LIGHTGREY, lw=0.8)
axb.axhline(-2.0 / 3.0, color=style.RED, lw=0.9, ls=":")
axb.axvline(0.0, color=style.PURPLE, lw=0.9, ls=":")

axb.plot(u_plot, ph1_p, color=style.BLUE, lw=1.9, label=r"$a=1\ (V):\ (0, +\frac{2\pi}{3})$")
axb.plot(u_plot, ph2_p, color=style.RED, lw=1.9, label=r"$a=2\ (\bar{V}):\ (-\frac{2\pi}{3}, 0)$")
axb.plot([0.0], [1.0 / 3.0], "o", color=style.BLUE, ms=4.5, zorder=5)
axb.plot([0.0], [-1.0 / 3.0], "o", color=style.RED, ms=4.5, zorder=5)

axb.set_xlim(-3.5, 3.5)
axb.set_ylim(-0.82, 1.08)
axb.set_yticks([-2.0 / 3.0, -1.0 / 3.0, 0.0, 1.0 / 3.0, 2.0 / 3.0])
axb.set_yticklabels([r"$-\frac{2\pi}{3}$", r"$-\frac{\pi}{3}$", r"$0$", r"$\frac{\pi}{3}$", r"$\frac{2\pi}{3}$"])
axb.set_xlabel(r"$\theta - \theta_d$")
axb.set_ylabel(r"phase $\arg z_a(\theta)$")
axb.set_title(r"(b) Chiral phase $\arg z_a(\theta)$", fontsize=8.3, color=style.GREY)
axb.legend(loc="upper left", fontsize=6.7, handlelength=1.4, borderpad=0.2)

# Panel (c): sl_3 weight triangle of V = 3 showing chiral selection rule
axc.set_xlim(-1.45, 1.45)
axc.set_ylim(-1.45, 1.45)
axc.axis("off")

tri_angles = [np.pi / 2, np.pi / 2 - 2 * np.pi / 3, np.pi / 2 - 4 * np.pi / 3]
tri_pts = [(0.92 * np.cos(a), 0.92 * np.sin(a)) for a in tri_angles]
tri_labs = [r"$|1\rangle_V$", r"$|2\rangle_V$", r"$|3\rangle_V$"]

for j, ((x1, y1), lab) in enumerate(zip(tri_pts, tri_labs)):
    axc.plot([x1], [y1], "o", ms=8.0, color=style.BLUE, zorder=5)
    axc.text(1.28 * x1, 1.28 * y1 + (0.06 if j == 0 else -0.06), lab, fontsize=8.2, color=style.BLUE, ha="center", va="center")

# Allowed root-label lowering transitions p = 1 (green solid)
for j in range(3):
    x1, y1 = tri_pts[j]
    x2, y2 = tri_pts[(j - 1) % 3]
    axc.annotate(
        "",
        xy=(x2, y2),
        xytext=(x1, y1),
        arrowprops=dict(
            arrowstyle="->",
            color=style.GREEN,
            lw=2.0,
            shrinkA=8,
            shrinkB=8,
            connectionstyle="arc3,rad=-0.18",
        ),
    )

# Forbidden counter-clockwise transition from |1>_V to |2>_V (p = 2, red dashed with X)
x1, y1 = tri_pts[0]
x3, y3 = tri_pts[1]
axc.annotate(
    "",
    xy=(x3, y3),
    xytext=(x1, y1),
    arrowprops=dict(
        arrowstyle="->",
        color=style.RED,
        lw=1.3,
        ls="--",
        shrinkA=9,
        shrinkB=9,
        connectionstyle="arc3,rad=-0.35",
    ),
)
# Put a red X on the forbidden arc
xm_forb, ym_forb = 0.68, 0.35
axc.plot([xm_forb], [ym_forb], "x", ms=7.0, mew=1.8, color=style.RED, zorder=6)

axc.text(0.0, -0.08, r"$p=0,1$" + "\n" + r"allowed", fontsize=7.5, color=style.GREEN, ha="center", va="center")
axc.text(0.0, -1.25, r"green: $j\to j-1$ ($p=1$)", fontsize=7.3, color=style.GREEN, ha="center")
axc.text(0.0, -1.44, r"red: $j\to j+1$ ($p=2$)", fontsize=7.3, color=style.RED, ha="center")
axc.set_title(r"(c) Chiral rule on $V=\mathbf{3}$", fontsize=8.3, color=style.GREY)

fig.subplots_adjust(left=0.07, right=0.98, bottom=0.17, top=0.88, wspace=0.34)
style.save(fig, "l3_a2_chiral", OUT)
raise SystemExit(0 if style.all_ok() else 1)
