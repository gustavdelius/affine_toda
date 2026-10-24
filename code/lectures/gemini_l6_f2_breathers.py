"""Lectures on integrable defects (defects-gemini), Lecture 6, figure 2 (fig-lec6-breathers):
The Miraculous Breather Cancellation, Higher Sine-Gordon Breathers, and Universal a_r^(1)
Toda Particle Transmission Phases (sec-lec6-breathers, eq-lec6-b1-transmission,
eq-lec6-ar-particle-t, eq-lec6-crossing-check, eq-lec6-fusing-check, exr-lec6-higher-breathers).

Panel (a): Demonstration of the B_1 breather cancellation across three couplings
           lambda = sqrt(2), 2*sqrt(2), 4*sqrt(2): the constituent soliton scalar-factor
           phases arg[Phi(theta_1)*Phi(theta_2)] (dashed curves, strongly lambda-dependent)
           combine with the singlet matrix eigenvalue factor to collapse onto the single
           coupling-independent classical curve arg T_{B_1}(theta) = arg T^{cl}(theta - theta_d)
           (bold red curve), while the moduli |Phi(theta_1)*Phi(theta_2)| * |lambda_singlet| = 1.
Panel (b): Universal real-coupling a_3^(1) (h = 4) Toda particle transmission phases
           arg T_a(theta)/pi for species a = 1 (gamma_1 = pi/4), a = 2 (gamma_2 = pi/2, equal
           to B_1), and a = 3 (gamma_3 = 3*pi/4), alongside the coupling-dependent higher
           sine-Gordon breather phases arg T_{B_2}(theta)/pi and arg T_{B_3}(theta)/pi at
           lambda = 3.5 (Problem 6.1).

Checks:
  1. Sine-Gordon B_1 bootstrap: On the singlet |B_1> in V(theta - i*u/2) (x) V(theta + i*u/2)
     (u = pi - pi/lambda), the two-soliton product (T_1 (x) 1)(1 (x) T_2) acts as 1_W times
     T^{cl}(theta - theta_d), independent of n, s, and lambda.
  2. Problem 6.1: Higher breather transmission factors T_{B_m}(theta) (m = 2, 3) satisfy
     unitarity |T_{B_m}| = 1, crossing T_{B_m}(theta)*T_{B_m}(theta + i*pi) = 1, and vary
     non-trivially with lambda.
  3. Higher-rank a_2^(1) (h = 3) and a_3^(1) (h = 4) lightest breathers: Fusing V and Vbar on
     the U_q(a_r^(1)) singlet at theta -+ i*(pi/2 - pi/(h*omega)) gives 1_W times T_1(theta)
     and T_r(theta).
  4. Universal a_r^(1) particle factors T_a(theta) (h = 3, 4, 5) satisfy unitarity, crossing
     T_a(-theta)^{-1} = T_{h-a}(i*pi - theta), and bulk fusing T_a(theta + i*gamma_b)*T_b(theta - i*gamma_a)
     = T_{a+b}(theta).

Writes lecture-notes/defects-gemini/figures/l6_breathers.svg.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from boundary_common import an_vector, cop, prop_residual, q_phys  # noqa: E402
from defect_common import (  # noqa: E402
    an_scalar_constant,
    an_scalar_g,
    an_scalar_gbar,
    an_T_vector,
    an_weyl_module,
    invariant_vector,
    sg_tmatrix,
    sg_transmission_scalar,
)

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects-gemini" / "figures"

thd = 0.25


def T_particle(D, a, h):
    """Universal type-I_- particle transmission factor (eq-lec6-ar-particle-t)."""
    gam_a = np.pi * a / h
    return (np.exp(-D) + np.exp(1j * gam_a)) / (np.exp(-D) + np.exp(-1j * gam_a))


# ---- Check 1: Sine-Gordon B_1 singlet bootstrap cancels n, s, and lambda ----
dev_b1 = 0.0
fit_b1 = 0.0
for lam_v in (np.sqrt(2.0), 2.0 * np.sqrt(2.0), 3.0 * np.sqrt(2.0)):
    qq = q_phys(lam_v)
    u_fus = np.pi * (1.0 - 1.0 / lam_v)
    for th_v in (-0.45, 0.35):
        t1, t2 = th_v - 0.5j * u_fus, th_v + 0.5j * u_fus
        rel_sv, w_sing = invariant_vector(np.exp(lam_v * t1), np.exp(lam_v * t2), qq)
        for s_v in (0.0, 0.35):
            # Include q^{-s/2} normalization so det(M_n) = 1 + Y^2 for general s
            T1, _, W = sg_tmatrix(t1, lam_v, qq, thd, 5, s=s_v, scalar=True)
            T2 = sg_tmatrix(t2, lam_v, qq, thd, 5, s=s_v, scalar=True)[0]
            T1 = (qq ** (-0.5 * s_v)) * T1
            T2 = (qq ** (-0.5 * s_v)) * T2
            L12 = np.kron(T1, np.eye(2)) @ np.kron(np.eye(2), T2)
            for n0 in (-1, 0, 1):
                e_n = np.zeros(W["dim"])
                e_n[list(W["n"]).index(n0)] = 1.0
                res_p, c_val = prop_residual(L12 @ np.kron(w_sing, e_n), np.kron(e_n, w_sing))
                fit_b1 = max(fit_b1, rel_sv, res_p)
                dev_b1 = max(dev_b1, abs(c_val - T_particle(th_v - thd, 1, 2)))

report(
    f"sine-Gordon B_1 singlet bootstrap equals T^cl(theta - theta_d) across 3 couplings, 2 s-values, and 3 defect labels "
    f"[max dev {dev_b1:.1e}, prop {fit_b1:.1e}]",
    dev_b1 < 1e-8 and fit_b1 < 1e-10,
)


# ---- Check 2: Problem 6.1 (Higher breathers T_{B_m} for m = 2, 3) ----
def T_Bm(D, m, lam_v):
    out = 1.0 + 0.0j
    for k in range(1, m + 1):
        shift = (1j * np.pi / lam_v) * (0.5 * (m + 1) - k)
        out = out * T_particle(D + shift, 1, 2)
    return out


u_grid = np.linspace(-2.0, 2.0, 30)
bm_unit_err = 0.0
bm_cross_err = 0.0
for m_val in (2, 3):
    vals_r = np.array([T_Bm(u, m_val, 3.5) for u in u_grid])
    vals_sh = np.array([T_Bm(u + 1j * np.pi, m_val, 3.5) for u in u_grid])
    bm_unit_err = max(bm_unit_err, np.max(np.abs(np.abs(vals_r) - 1.0)))
    bm_cross_err = max(bm_cross_err, np.max(np.abs(vals_r * vals_sh - 1.0)))

# Verify non-trivial lambda-dependence for m = 2 (unlike m = 1)
diff_m2_lam = abs(T_Bm(0.5, 2, 2.5) - T_Bm(0.5, 2, 4.5))
diff_m1_lam = abs(T_Bm(0.5, 1, 2.5) - T_Bm(0.5, 1, 4.5))

report(
    f"Problem 6.1: higher breather factors T_{{B_m}} (m=2,3) satisfy unitarity [{bm_unit_err:.1e}] and crossing "
    f"[{bm_cross_err:.1e}], and depend on lambda for m>=2 [|Delta T_B2|={diff_m2_lam:.3f} vs |Delta T_B1|={diff_m1_lam:.1e}]",
    bm_unit_err < 1e-14 and bm_cross_err < 1e-14 and diff_m2_lam > 0.1 and diff_m1_lam < 1e-15,
)


# ---- Check 3: Higher-rank a_2^(1) and a_3^(1) lightest breather singlet bootstrap ----
def ar_singlet_bootstrap(r_val, om_val, th_val, cL, cR, R2=14):
    h_val = r_val + 1
    qq = q_phys(om_val)
    mu_val = -np.exp(-0.5j * np.pi * om_val) / (qq - 1.0 / qq)
    W = an_weyl_module(r_val, -1, qq, np.exp(om_val * thd), R2, s=np.linspace(-0.2, 0.2, h_val), mu=mu_val)
    dW = W["dim"]
    ub = 0.5 * np.pi - np.pi / (h_val * om_val)
    thL, thR = th_val - 1j * ub, th_val + 1j * ub

    A_rep = an_vector(r_val, qq, np.exp(om_val * thL), conj=cL)
    B_rep = an_vector(r_val, qq, np.exp(om_val * thR), conj=cR)
    M_sing = np.vstack(
        [cop(A_rep, B_rep, kd, j) - (np.eye(h_val * h_val) if kd == "k" else 0) for j in range(h_val) for kd in "efk"]
    )
    _, sv, Vh = np.linalg.svd(M_sing, full_matrices=False)
    inv_vec = Vh[-1].conj()

    c_norm = an_scalar_constant(r_val, om_val, qq)

    def T_full(t_arg, conj_flag):
        M_mat = an_T_vector(W, np.exp(om_val * (t_arg - thd - 0.5j * np.pi)), conj_flag)
        sc = (
            an_scalar_gbar(t_arg - thd, r_val, om_val, 400, c_norm)
            if conj_flag
            else an_scalar_g(t_arg - thd, r_val, om_val, 400, c_norm)
        )
        return sc * M_mat

    TL, TR = T_full(thL, cL), T_full(thR, cR)
    k_deep = np.where(W["deep"])[0][0]
    v = np.kron(inv_vec, np.eye(dW)[k_deep]).reshape(h_val, h_val * dW)
    v = (TR @ v.T).T.reshape(h_val * dW, h_val)
    v = TL @ v
    res_p, ratio_val = prop_residual(v.reshape(-1), np.kron(np.eye(dW)[k_deep], inv_vec))
    return sv[-1] / sv[0], res_p, ratio_val


ar_dev = 0.0
ar_res = 0.0
for r_v in (2, 3):
    h_v = r_v + 1
    for (cL, cR), k_exp in (((True, False), 1), ((False, True), r_v)):
        sv_r, res_r, val_r = ar_singlet_bootstrap(r_v, np.sqrt(2.0), 0.42, cL, cR, R2=14)
        ar_res = max(ar_res, sv_r, res_r)
        ar_dev = max(ar_dev, abs(val_r - T_particle(0.42 - thd, k_exp, h_v)))

report(
    f"a_2^(1) and a_3^(1) V (x) Vbar singlet bootstrap reproduces T_1(theta) and T_r(theta) times 1_W "
    f"[max dev {ar_dev:.1e}, singlet/prop res {ar_res:.1e}]",
    ar_dev < 1e-7 and ar_res < 1e-12,
)

# ---- Check 4: Unitarity, crossing, and bulk fusing bootstrap for a_r^(1) particle factors ----
toda_err = 0.0
for h_v in (3, 4, 5):
    for a in range(1, h_v):
        # Unitarity and crossing
        u_vals = np.linspace(-2.0, 2.0, 20)
        Ta_u = T_particle(u_vals, a, h_v)
        Tbar_a = 1.0 / T_particle(-u_vals - 2.0 * thd, a, h_v)
        Ta_cross = T_particle(1j * np.pi - u_vals - 2.0 * thd, h_v - a, h_v)
        toda_err = max(toda_err, np.max(np.abs(np.abs(Ta_u) - 1.0)), np.max(np.abs(Tbar_a - Ta_cross)))
        # Bulk fusing bootstrap: T_a(D + i*gamma_b) * T_b(D - i*gamma_a) = T_{a+b}(D)
        for b in range(1, h_v - a):
            gam_a = np.pi * a / h_v
            gam_b = np.pi * b / h_v
            fused_ab = T_particle(u_vals + 1j * gam_b, a, h_v) * T_particle(u_vals - 1j * gam_a, b, h_v)
            target_ab = T_particle(u_vals, a + b, h_v)
            toda_err = max(toda_err, np.max(np.abs(fused_ab - target_ab)))

report(
    f"a_r^(1) particle factors T_a(theta) (h=3,4,5) satisfy unitarity, crossing, and bulk fusing bootstrap [err {toda_err:.1e}]",
    toda_err < 1e-14,
)

# ---- Figure ----
fig, (axa, axb) = plt.subplots(
    1, 2, figsize=(7.0, 2.95), gridspec_kw=dict(width_ratios=[1.08, 1.12], wspace=0.28)
)

# Panel (a): Miraculous Breather Cancellation across couplings lambda = sqrt(2), 2*sqrt(2), 4*sqrt(2)
u_plot = np.linspace(-2.2, 2.2, 45)
axa.axhline(0.0, color=style.LIGHTGREY, lw=0.8)
axa.axhline(1.0, color=style.LIGHTGREY, lw=0.8)
axa.axvline(0.0, color=style.PURPLE, lw=0.9, ls=":")

for lam_c, col_c, ls_c, lab_c in [
    (np.sqrt(2.0), style.BLUE, "--", r"constituents $\arg[\Phi(\theta_1)\Phi(\theta_2)]$ ($\lambda=\sqrt{2}$)"),
    (2.0 * np.sqrt(2.0), style.TEAL, "-.", r"constituents $\arg[\Phi(\theta_1)\Phi(\theta_2)]$ ($\lambda=2\sqrt{2}$)"),
    (4.0 * np.sqrt(2.0), style.ORANGE, ":", r"constituents $\arg[\Phi(\theta_1)\Phi(\theta_2)]$ ($\lambda=4\sqrt{2}$)"),
]:
    u_f = np.pi * (1.0 - 1.0 / lam_c)
    phi_prod = np.array(
        [sg_transmission_scalar(u - 0.5j * u_f, lam_c) * sg_transmission_scalar(u + 0.5j * u_f, lam_c) for u in u_plot]
    )
    ph_const = np.unwrap(np.angle(phi_prod)) / np.pi
    axa.plot(u_plot, np.clip(ph_const, -0.2, 1.75), color=col_c, lw=1.45, ls=ls_c, label=lab_c)

# Universal collapsed B_1 phase curve arg T_{B_1}(u) / pi
u_dense = np.linspace(-2.2, 2.2, 240)
ph_b1 = np.unwrap(np.angle(T_particle(u_dense, 1, 2))) / np.pi
axa.plot(
    u_dense,
    ph_b1,
    color=style.RED,
    lw=2.3,
    label=r"fused singlet $\arg\mathbb{T}_{B_1}(\theta) = \arg\mathbb{T}^{\rm cl}$ (all $\lambda$!)",
)

axa.annotate(
    r"all $\lambda, n, s$ cancel" + "\n" + r"on singlet $|B_1\rangle$!",
    xy=(0.0, 0.5),
    xytext=(-1.15, 0.68),
    fontsize=6.6,
    color=style.RED,
    ha="center",
    arrowprops=dict(arrowstyle="->", color=style.RED, lw=0.8, shrinkB=4),
)

axa.set_xlim(-2.2, 2.2)
axa.set_ylim(-0.08, 1.82)
axa.set_xlabel(r"relative rapidity $\theta - \theta_d$")
axa.set_ylabel(r"transmission phase / $\pi$")
axa.set_title(r"(a) Singlet cancellation onto $\mathbb{T}_{B_1} = \mathbb{T}^{\rm cl}$", fontsize=8.3, color=style.GREY)
axa.legend(loc="upper left", fontsize=5.8, handlelength=1.4, borderpad=0.16, labelspacing=0.16)

# Panel (b): Universal a_3^(1) (h = 4) Toda Particle Phases & Higher Breathers B_2, B_3
axb.axhline(0.0, color=style.LIGHTGREY, lw=0.8)
axb.axhline(1.0, color=style.LIGHTGREY, lw=0.8)
axb.axhline(1.5, color=style.LIGHTGREY, lw=0.8, ls=":")
axb.axvline(0.0, color=style.PURPLE, lw=0.9, ls=":")

for a_sp, col_sp, lw_sp, lab_sp in [
    (1, style.BLUE, 1.8, r"$a_3^{(1)}$ species $a=1$ ($\gamma_1 = \pi/4$)"),
    (2, style.RED, 2.0, r"$a_3^{(1)}$ species $a=2$ ($\gamma_2 = \pi/2$, $= B_1$)"),
    (3, style.TEAL, 1.8, r"$a_3^{(1)}$ species $a=3$ ($\gamma_3 = 3\pi/4$)"),
]:
    ph_a = np.unwrap(np.angle(T_particle(u_dense, a_sp, 4))) / np.pi
    axb.plot(u_dense, ph_a, color=col_sp, lw=lw_sp, label=lab_sp)
    axb.plot([0.0], [a_sp / 4.0], "o", color=col_sp, ms=4.3, zorder=6)

# Also plot higher sine-Gordon breathers B_2 and B_3 at lambda = 3.5 (rescaled by 1/m so asymptotic jump is 1*pi)
ph_b2 = np.unwrap(np.angle([T_Bm(u, 2, 3.5) for u in u_dense])) / np.pi
ph_b3 = np.unwrap(np.angle([T_Bm(u, 3, 3.5) for u in u_dense])) / np.pi
axb.plot(u_dense, ph_b2 / 2.0, color=style.ORANGE, lw=1.35, ls="--", label=r"$\frac{1}{2}\arg\mathbb{T}_{B_2}(\theta)/\pi$ ($\lambda=3.5$)")
axb.plot(u_dense, ph_b3 / 3.0, color=style.PURPLE, lw=1.35, ls=":", label=r"$\frac{1}{3}\arg\mathbb{T}_{B_3}(\theta)/\pi$ ($\lambda=3.5$)")

axb.set_xlim(-2.2, 2.2)
axb.set_ylim(-0.08, 1.65)
axb.set_xlabel(r"relative rapidity $\theta - \theta_d$")
axa_ylabel = r"transmission phase $\arg\mathbb{T}_a(\theta) / \pi$"
axb.set_ylabel(axa_ylabel)
axb.set_title(r"(b) $a_3^{(1)}$ Toda particles ($h=4$) & higher $B_m$", fontsize=8.3, color=style.GREY)
axb.legend(loc="upper left", fontsize=5.8, handlelength=1.4, borderpad=0.18, labelspacing=0.18)

fig.subplots_adjust(left=0.08, right=0.98, bottom=0.17, top=0.88, wspace=0.28)
style.save(fig, "l6_breathers", OUT)
raise SystemExit(0 if style.all_ok() else 1)
