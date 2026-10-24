"""Lectures on integrable defects (defects-gemini), Lecture 2, figure 1 (fig-lec2-soliton-defect):
The Four Physical Wonders of a Sine-Gordon Soliton Crossing a Defect (sec-lec2-soliton-defect, eq-lec2-z-delay).

Panel (a): The Bäcklund delay factor z_{I_-}(theta) = -coth((theta - theta_d)/2), the spatial shift
           ln|z_{I_-}(theta)|, and the "half-soliton" relation z_{odd}(theta)^2 = tanh^2((theta - theta_d)/2) = A_{11}(theta - theta_d).
Panel (b): Slow soliton delay (theta = 0.45 < theta_d = 0.9): emerges as a soliton (positive topological charge density,
           blue ridge) delayed by Delta x = -ln|z|/(m_0 cosh theta) < 0.
Panel (c): Fast soliton-to-antisoliton conversion / topological NOT gate (theta = 1.45 > theta_d = 0.9): since z < 0,
           the incoming soliton (blue ridge) emerges on x > 0 as an antisoliton (red trough), leaving a zero-energy
           even jump psi(0) - phi(0) = 4*pi/beta_SG on the defect.
Panel (d): Near-critical soliton capture (theta = theta_d - 5e-4): the defect holds the odd jump 2*pi/beta_SG for a long
           duration Delta t = ln|z|/(m_0 sinh theta), storing the exact soliton energy B_d = M cosh theta_d while total
           energy E_L + B_d + E_R = M cosh theta stays constant.

Checks:
  1. Momentum conservation functional identity (1/2) f'(p) g'(q) = V(phi) - V(psi) holds for sine-Gordon and fails for phi^4.
  2. Exact Hirota pairs (phi, psi) satisfy the light-cone defect conditions (eq-lec2-sg-lc) at x = 0 across all times t.
  3. Outgoing topological charge Q_R = (beta_SG / 2*pi) int_0^inf dx partial_x psi at late time is +1 for theta < theta_d
     and -1 for theta > theta_d, and the ridge displacement matches Delta x = -ln|z|/(m_0 cosh theta).
  4. Total energy E_L + B_d + E_R and total momentum P_L + P_d + P_R are conserved across the capture passage, the plateau
     of B_d matches M cosh theta_d, and z_{odd}(theta)^2 = A_{11}(theta - theta_d).

Writes lecture-notes/defects-gemini/figures/l2_soliton_defect.svg.
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

m0 = 1.0
beta = 1.0  # units where beta_SG = 1, so classical soliton mass M = 8*m0/beta^2 = 8
M_sol = 8.0 * m0 / beta**2
thd = 0.85
sigma = np.exp(-thd)


def z_delay(th, th_d=thd):
    return (np.exp(-th) + np.exp(-th_d)) / (np.exp(-th) - np.exp(-th_d))


def sg_soliton_pair(x, t, th, th_d=thd):
    """Exact Hirota 1-soliton pair (phi on x<=0, psi on x>=0) and their x, t derivatives."""
    z = z_delay(th, th_d)
    ch = np.cosh(th)
    sh = np.sinh(th)
    # E = +i * exp(X) so phi = +(4/beta) arctan(exp(X)), increasing from 0 to 2*pi/beta (Q = +1)
    X = m0 * (x * ch - t * sh)
    E = 1j * np.exp(X)
    Ep = z * E

    phi = np.real((2j / beta) * np.log((1.0 - E) / (1.0 + E)))
    psi = np.real((2j / beta) * np.log((1.0 - Ep) / (1.0 + Ep)))

    # Analytical derivatives from d/dX [(2i/beta) ln((1-E)/(1+E))] = -(4i/beta) E / (1 - E^2)
    dphi_dX = np.real((-4j / beta) * E / (1.0 - E**2))
    dpsi_dX = np.real((-4j / beta) * Ep / (1.0 - Ep**2))

    phi_x = m0 * ch * dphi_dX
    phi_t = -m0 * sh * dphi_dX
    psi_x = m0 * ch * dpsi_dX
    psi_t = -m0 * sh * dpsi_dX
    return phi, psi, phi_x, phi_t, psi_x, psi_t


def B_defect(phi0, psi0, th_d=thd):
    sig = np.exp(-th_d)
    return -(2.0 * m0 / beta**2) * (
        sig * (np.cos(0.5 * beta * (phi0 + psi0)) - 1.0)
        + (1.0 / sig) * (np.cos(0.5 * beta * (phi0 - psi0)) - 1.0)
    )


def P_defect(phi0, psi0, th_d=thd):
    sig = np.exp(-th_d)
    return (2.0 * m0 / beta**2) * (
        sig * (np.cos(0.5 * beta * (phi0 + psi0)) - 1.0)
        - (1.0 / sig) * (np.cos(0.5 * beta * (phi0 - psi0)) - 1.0)
    )


# ---- Check 1: Momentum factorization identity ----
p_sym, q_sym = sp.symbols("p q", real=True)
phi_s, psi_s = p_sym + q_sym, p_sym - q_sym
V_sg = 1 - sp.cos(phi_s) - (1 - sp.cos(psi_s))
fact_sg = 2 * sp.sin(p_sym) * sp.sin(q_sym)
V_phi4 = phi_s**4 - psi_s**4
# A function W(p,q) factors as F(p)*G(q) iff d^2/dp dq ln(W(p,q)) = 0
cross_phi4 = sp.simplify(sp.diff(sp.log(V_phi4), p_sym, q_sym))
cross_sg = sp.simplify(sp.diff(sp.log(fact_sg), p_sym, q_sym))
report(
    "momentum factorization: V(p+q)-V(p-q) factors as F(p)G(q) for sine-Gordon and fails for phi^4",
    sp.simplify(V_sg - fact_sg) == 0 and cross_sg == 0 and cross_phi4 != 0,
)

# ---- Check 2: Defect conditions at x = 0 for slow, fast, and near-critical solitons ----
ts_check = np.linspace(-6.0, 6.0, 200)
max_res = 0.0
for th_test in (0.45, 1.45, thd - 5e-4):
    phi0, psi0, phi_x0, phi_t0, psi_x0, psi_t0 = sg_soliton_pair(0.0, ts_check, th_test)
    # Light-cone conditions eq-lec2-sg-lc
    res_plus = (phi_t0 + phi_x0 - (psi_t0 + psi_x0)) + (2.0 * sigma * m0 / beta) * np.sin(0.5 * beta * (phi0 + psi0))
    res_minus = (phi_t0 - phi_x0 + (psi_t0 - psi_x0)) - (2.0 * m0 / (sigma * beta)) * np.sin(0.5 * beta * (phi0 - psi0))
    max_res = max(max_res, np.max(np.abs(res_plus)), np.max(np.abs(res_minus)))

report(
    f"exact Hirota 1-soliton pairs satisfy the sine-Gordon defect conditions at x=0 across all regimes [max res {max_res:.1e}]",
    max_res < 1e-12,
)

# ---- Check 3: Delay shift and NOT-gate charge conversion ----
xs_right = np.linspace(0.0, 35.0, 5000)
dx_grid = xs_right[1] - xs_right[0]
t_late = 18.0

_, _, _, _, psi_x_slow, _ = sg_soliton_pair(xs_right, t_late, 0.45)
_, _, _, _, psi_x_fast, _ = sg_soliton_pair(xs_right, t_late, 1.45)

Q_slow = np.trapz(psi_x_slow, xs_right) / (2.0 * np.pi / beta)
Q_fast = np.trapz(psi_x_fast, xs_right) / (2.0 * np.pi / beta)

x_peak_slow = xs_right[np.argmax(np.abs(psi_x_slow))]
x_exp_slow = np.tanh(0.45) * t_late - np.log(abs(z_delay(0.45))) / (m0 * np.cosh(0.45))
x_peak_fast = xs_right[np.argmax(np.abs(psi_x_fast))]
x_exp_fast = np.tanh(1.45) * t_late - np.log(abs(z_delay(1.45))) / (m0 * np.cosh(1.45))

report(
    f"slow soliton (th<thd) transmits with Q={Q_slow:+.3f} and shift dx={x_peak_slow - np.tanh(0.45)*t_late:.3f}; "
    f"fast soliton (th>thd) converts to antisoliton Q={Q_fast:+.3f}",
    abs(Q_slow - 1.0) < 1e-3
    and abs(Q_fast - (-1.0)) < 1e-3
    and abs(x_peak_slow - x_exp_slow) < 2 * dx_grid
    and abs(x_peak_fast - x_exp_fast) < 2 * dx_grid,
)

# ---- Check 4: Capture energy/momentum conservation and half-soliton theorem ----
th_cap = thd - 5e-4
ts_cap = np.linspace(-5.0, 24.5, 160)
xs_left = np.linspace(-25.0, 0.0, 1500)
xs_r_cap = np.linspace(0.0, 35.0, 2000)

E_L_list, E_R_list, B_d_list, P_tot_list = [], [], [], []
for tc in ts_cap:
    phi_l, _, phi_xl, phi_tl, _, _ = sg_soliton_pair(xs_left, tc, th_cap)
    _, psi_r, _, _, psi_xr, psi_tr = sg_soliton_pair(xs_r_cap, tc, th_cap)
    phi0, psi0, _, _, _, _ = sg_soliton_pair(0.0, tc, th_cap)

    e_l = np.trapz(0.5 * phi_tl**2 + 0.5 * phi_xl**2 + (m0**2 / beta**2) * (1.0 - np.cos(beta * phi_l)), xs_left)
    e_r = np.trapz(0.5 * psi_tr**2 + 0.5 * psi_xr**2 + (m0**2 / beta**2) * (1.0 - np.cos(beta * psi_r)), xs_r_cap)
    p_l = np.trapz(-phi_tl * phi_xl, xs_left)
    p_r = np.trapz(-psi_tr * psi_xr, xs_r_cap)

    E_L_list.append(e_l)
    E_R_list.append(e_r)
    B_d_list.append(B_defect(phi0, psi0))
    P_tot_list.append(p_l + p_r + P_defect(phi0, psi0))

E_L_arr = np.array(E_L_list)
E_R_arr = np.array(E_R_list)
B_d_arr = np.array(B_d_list)
E_tot_arr = E_L_arr + E_R_arr + B_d_arr
P_tot_arr = np.array(P_tot_list)

th_grid = np.linspace(-2.0, 2.5, 100)
z_odd_sq = (1.0 / z_delay(th_grid)) ** 2
A11_grid = np.tanh(0.5 * (th_grid - thd)) ** 2

report(
    "capture & half-soliton: E_tot=M cosh(th), P_tot=M sinh(th) conserved to 1e-4, odd jump stores (M cosh thd, M sinh thd), and z_odd^2 = A_11",
    np.allclose(E_tot_arr, M_sol * np.cosh(th_cap), rtol=1e-4)
    and np.allclose(P_tot_arr, M_sol * np.sinh(th_cap), rtol=1e-4)
    and abs(B_defect(0.0, 2.0 * np.pi / beta) - M_sol * np.cosh(thd)) < 1e-12
    and abs(P_defect(0.0, 2.0 * np.pi / beta) - M_sol * np.sinh(thd)) < 1e-12
    and abs(B_defect(0.0, 4.0 * np.pi / beta)) < 1e-12
    and np.allclose(z_odd_sq, A11_grid),
)

# ---- Figure ----
fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.0))
(axa, axb), (axc, axd) = axes

# Panel (a): Delay factor and half-soliton relation
u_left = np.linspace(-2.5, -0.06, 200)
u_right = np.linspace(0.06, 2.5, 200)
u_all = np.linspace(-2.5, 2.5, 400)

axa.axhline(0.0, color=style.LIGHTGREY, lw=0.8)
axa.axvline(0.0, color=style.PURPLE, lw=1.1, ls=":")

axa.plot(u_left, np.log(np.abs(1.0 / np.tanh(0.5 * u_left))), color=style.BLUE, lw=1.9, label=r"$\ln|z_{\rm I_-}|$ ($z>0$, soliton)")
axa.plot(u_right, np.log(np.abs(1.0 / np.tanh(0.5 * u_right))), color=style.RED, lw=1.9, label=r"$\ln|z_{\rm I_-}|$ ($z<0$, NOT gate)")
axa.plot(u_all, np.tanh(0.5 * u_all) ** 2, color=style.GREEN, lw=1.6, ls="--", label=r"$z_{\rm odd}^2 = A_{11}(\theta-\theta_d)$")

# Mark the rapidities of panels (b) and (c)
for th_m, col_m, lab_m in [(0.45, style.BLUE, "(b)"), (1.45, style.RED, "(c)")]:
    u_m = th_m - thd
    y_m = np.log(abs(z_delay(th_m)))
    axa.plot([u_m], [y_m], "o", color=col_m, ms=4.8, zorder=5)
    axa.text(u_m + (0.12 if u_m > 0 else -0.12), y_m + 0.15, lab_m, color=col_m, fontsize=8, ha="left" if u_m > 0 else "right")

axa.set_xlim(-2.5, 2.5)
axa.set_ylim(-0.15, 3.3)
axa.set_xlabel(r"relative rapidity $\theta - \theta_d$")
axa.set_ylabel(r"shift $\ln|z_{\rm I_-}|$ and $z_{\rm odd}^2$")
axa.set_title(r"(a) Bäcklund delay & half-soliton identity", fontsize=8.5, color=style.GREY)
axa.legend(loc="upper left", fontsize=7.0, handlelength=1.6)

# Helper for spacetime panels (b) and (c)
x_sp = np.linspace(-4.5, 4.5, 220)
t_sp = np.linspace(-5.0, 5.5, 220)
XX, TT = np.meshgrid(x_sp, t_sp)

for ax, th_p, title_p in [
    (axb, 0.45, r"(b) Slow soliton ($\theta<\theta_d$): delayed kink"),
    (axc, 1.45, r"(c) Fast soliton ($\theta>\theta_d$): NOT gate ($s\to\bar{s}$)"),
]:
    _, _, phi_x_g, _, psi_x_g, _ = sg_soliton_pair(XX, TT, th_p)
    rho = np.where(XX <= 0.0, phi_x_g, psi_x_g) / (2.0 * m0 / beta)  # normalized in [-1, 1]

    ax.imshow(
        rho,
        extent=[x_sp[0], x_sp[-1], t_sp[0], t_sp[-1]],
        origin="lower",
        aspect="auto",
        cmap="coolwarm",
        vmin=-1.0,
        vmax=1.0,
         alpha=0.85,
        rasterized=True,
    )
    ax.axvline(0.0, color=style.PURPLE, lw=2.0)
    # Undelayed continuation on x > 0
    t_pos = np.linspace(0.0, 5.5, 100)
    v_p = np.tanh(th_p)
    ax.plot(v_p * t_pos, t_pos, color=style.GREY, lw=1.2, ls=":")
    # Delayed worldline on x > 0
    dx_p = -np.log(abs(z_delay(th_p))) / (m0 * np.cosh(th_p))
    dt_p = np.log(abs(z_delay(th_p))) / (m0 * np.sinh(th_p))
    t_del = np.linspace(dt_p, 5.5, 100)
    ax.plot(v_p * t_del + dx_p, t_del, color=style.RED if th_p < thd else style.BLUE, lw=1.2, ls="--")

    ax.annotate(
        r"$\Delta t = \frac{\ln|z|}{m_0\sinh\theta}$",
        xy=(0.0, dt_p),
        xytext=(-4.1, 2.8),
        fontsize=7.6,
        color=style.GREY,
        arrowprops=dict(arrowstyle="->", color=style.GREY, lw=0.7, shrinkB=3),
    )
    if th_p > thd:
        ax.text(1.8, -3.5, r"outgoing $\bar{s}$" + "\n" + r"($z<0$)", fontsize=7.8, color=style.BLUE, ha="center")
    else:
        ax.text(1.8, -3.5, r"outgoing $s$" + "\n" + r"($z>0$)", fontsize=7.8, color=style.RED, ha="center")

    ax.set_xlim(-4.5, 4.5)
    ax.set_ylim(-5.0, 5.5)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$t$")
    ax.set_title(title_p, fontsize=8.5, color=style.GREY)

# Panel (d): Capture at theta = theta_d - 5e-4
dt_hold = np.log(abs(z_delay(th_cap))) / (m0 * np.sinh(th_cap))
axd.axvspan(0.0, dt_hold, color=style.PURPLE, alpha=0.08, lw=0)
axd.plot(ts_cap, E_L_arr / M_sol, color=style.BLUE, lw=1.7, label=r"$E_L / M$ ($x<0$)")
axd.plot(ts_cap, B_d_arr / M_sol, color=style.PURPLE, lw=1.9, label=r"$\mathcal{B}_d / M$ (defect)")
axd.plot(ts_cap, E_R_arr / M_sol, color=style.TEAL, lw=1.7, label=r"$E_R / M$ ($x>0$)")
axd.plot(ts_cap, E_tot_arr / M_sol, color=style.GREY, lw=1.2, ls="--", label=r"$(E_{\rm bulk}+\mathcal{B}_d)/M = \cosh\theta_d$")

axd.text(0.5 * dt_hold, 0.85, r"odd jump $\frac{2\pi}{\beta_{\rm SG}}$" + "\n" + r"$\mathcal{B}_d = M\cosh\theta_d$",
         fontsize=7.3, color=style.PURPLE, ha="center")
axd.set_xlim(-5.0, 24.5)
axd.set_ylim(-0.05, 1.68)
axd.set_xlabel(r"time $t$")
axd.set_ylabel(r"energy (in units of $M$)")
axd.set_title(r"(d) Soliton capture & energy conservation ($\theta\approx\theta_d$)", fontsize=8.5, color=style.GREY)
axd.legend(loc="center right", bbox_to_anchor=(1.0, 0.45), fontsize=6.8, handlelength=1.4)

fig.tight_layout(h_pad=1.8, w_pad=1.4)
style.save(fig, "l2_soliton_defect", OUT)
raise SystemExit(0 if style.all_ok() else 1)
