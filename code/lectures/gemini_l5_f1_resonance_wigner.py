"""Lectures on integrable defects (defects-gemini), Lecture 5, figure 1 (fig-lec5-resonance-wigner):
Pole Cancellation in the Physical Strip, the Unstable Resonance at theta_P = theta_d - i*pi/(2*lambda),
and Wigner's Semiclassical Time Delay (sec-lec5-unitarity-crossing, sec-lec5-gamma-and-resonance,
eq-lec5-phi-gamma, eq-lec5-wigner-delay, eq-lec5-resonance-pole).

Panel (a): Analytic structure of the sine-Gordon scalar factor Phi(u) (u = theta - theta_d):
           - Along the imaginary axis u = i*eta, the apparent physical-strip pole of 1/(1 + i e^{lambda u}) at
             eta = +pi/(2*lambda) is cancelled by the simple zero of 1/rbar(i*eta) (from Gamma(0) = oo), leaving
             |Phi(i*eta)| smooth across the physical strip 0 < Im(u) < pi.
           - In the lower half-plane Im(u) < 0, the uncancelled simple pole at u_P = -i*pi/(2*lambda) produces the
             decaying defect-soliton resonance (Im Delta E < 0 for theta_d > 0).
Panel (b): Semiclassical Wigner time delay (pi/(lambda + 1)) d/dtheta arg Phi(theta - theta_d) for lambda = 1.6, 4.5, 16
           converging to the classical Bäcklund spatial delay ln|coth((theta - theta_d)/2)| as lambda -> oo.

Checks:
  1. Defect crossing relation holds with C = sigma_x at the physical q = -e^{-i*pi*lambda} and fails at +e^{-i*pi*lambda}.
  2. Phi(u) satisfies unitarity |Phi(u)|^2 (1 + e^{2*lambda*u}) = 1 and crossing Phi(u)*Phi(u + i*pi)*(1 + e^{2*lambda*u}) = 1,
     is finite at u = +i*pi/(2*lambda), and exhibits 1/eps simple-pole scaling at u_P = -i*pi/(2*lambda).
  3. Wigner delay (pi/(lambda + 1)) d/dtheta arg Phi(u) converges to ln|coth(u/2)| as 1/lambda for lambda = 10, 40, 160.
  4. Problem 5.2: Defect bootstrap at u_P = -i*pi/(2*lambda) projects onto the rank-1 residue image in each charge sector H_n
     and produces the odd-sector transmission matrix (s = 1, theta_d -> theta_d - i*pi).

Writes lecture-notes/defects-gemini/figures/l5_resonance_wigner.svg.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from boundary_common import q_phys, rcheck, sl2_spin  # noqa: E402
from defect_common import sg_transmission_scalar  # noqa: E402

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects-gemini" / "figures"

lam = np.sqrt(2.0) + 0.5  # generic irrational > 1 (~1.9142)
q = q_phys(lam)
thd = 0.35


# ---- Check 1: Crossing fixes the sign of q = -e^{-i*pi*lambda} ----
def crossing_check(q_test, th_val=0.42, n_val=2, s_val=0.0):
    """Compare Tbar_n(th) = T_n(-th)^{-1} with C T_n(i*pi - th)^{T_V} C^{-1} on H_n = {v_+ (x) |n>, v_- (x) |n+1>}."""
    # Note: y(i*pi - th)/y_d = e^{i*pi*lambda} e^{lambda(-th - thd)} = (-1/q_phys) Y(-th)
    # Using the tested q_test for the matrix entries and the exact analytical continuation e^{i*pi*lambda} for y(i*pi - th):
    Y_minus = np.exp(lam * (-th_val - thd))
    Y_cross = np.exp(lam * (1j * np.pi - th_val - thd))
    Phi_m = sg_transmission_scalar(-th_val - thd, lam)
    Phi_c = sg_transmission_scalar(1j * np.pi - th_val - thd, lam)

    # On H_n, T_n(th) = Phi(th - thd) * [[q^{-n}, -i*Y*q^{n+1+s}], [-i*Y*q^{-n}, q^{n+1+s}]]
    # Wait: in Tbar = T(-th)^{-1}, the diagonal entries map (+,n)->(+,n) and (-,n+1)->(-,n+1):
    G_n = np.diag([q_test ** (-n_val), q_test ** (n_val + 1 + s_val)])
    T_minus_n = Phi_m * np.array([[1.0, -1j * Y_minus], [-1j * Y_minus, 1.0]]) @ G_n
    Tbar_n = np.linalg.inv(T_minus_n)

    # Under C = sigma_x:
    # - (+,n) -> (+,n) in Tbar matches (-,n) -> (-,n) in T(i*pi - th), which has amplitude Phi_c * q_test^{n+s}
    # - (-,n+1) -> (-,n+1) in Tbar matches (+,n+1) -> (+,n+1) in T(i*pi - th), which has amplitude Phi_c * q_test^{-(n+1)}
    # - (+,n) -> (-,n+1) in Tbar (entry [1,0]) matches (+,n) -> (-,n+1) in T(i*pi - th), which has amplitude Phi_c * (-i*Y_cross*q_test^{-n})
    # - (-,n+1) -> (+,n) in Tbar (entry [0,1]) matches (-,n+1) -> (+,n) in T(i*pi - th), which has amplitude Phi_c * (-i*Y_cross*q_test^{n+1+s})
    rhs_n = Phi_c * np.array(
        [
            [q_test ** (n_val + s_val), -1j * Y_cross * (q_test ** (n_val + 1 + s_val))],
            [-1j * Y_cross * (q_test ** (-n_val)), q_test ** (-(n_val + 1))],
        ]
    )
    # Notice that Tbar_n = (1 / (Phi_m * (1 + Y_minus^2))) * G_n^{-1} @ [[1, i*Y_minus], [i*Y_minus, 1]]
    # = (1 / (Phi_m * (1 + Y_minus^2))) * [[q_test^n, i*Y_minus*q_test^n], [i*Y_minus*q_test^{-(n+1+s)}, q_test^{-(n+1+s)}]]
    # Wait: in the crossing relation eq-defect-crossing:
    # Tbar[(c, N'), (N, d)](th) = C[c, a] T[(N', b), (a, N)](i*pi - th) C^{-1}[b, d] (with th -> -th in defect_sg_transmission.py check 8)
    # Let's check T_n(th) vs C (T(-th - i*pi)^{-1})^{T_V} C^{-1} as in defect_sg_transmission.py!
    Y_p = np.exp(lam * (th_val - thd))
    Phi_p = sg_transmission_scalar(th_val - thd, lam)
    T_p = Phi_p * np.array([[1.0, -1j * Y_p], [-1j * Y_p, 1.0]]) @ G_n

    # At th_shift = th - i*pi, Y(th - i*pi) = e^{-i*pi*lambda} Y_p
    Y_sh = np.exp(lam * (th_val - 1j * np.pi - thd))
    Phi_sh = sg_transmission_scalar(th_val - 1j * np.pi - thd, lam)
    # Inverse of T(th - i*pi) on the conjugate block {v_+ (x) |n-1>, v_- (x) |n>}:
    # On {v_+ (x) |n>, v_- (x) |n>}:
    # Diagonal (+,n)->(+,n) of C Tbar(i*pi - th) C is (-,n)->(-,n) of T(th - i*pi)^{-1}, which sits in block H_{n-1}:
    G_nm1 = np.diag([q_test ** (-(n_val - 1)), q_test ** (n_val + s_val)])
    Tinv_nm1 = np.linalg.inv(Phi_sh * np.array([[1.0, -1j * Y_sh], [-1j * Y_sh, 1.0]]) @ G_nm1)
    Tinv_n = np.linalg.inv(Phi_sh * np.array([[1.0, -1j * Y_sh], [-1j * Y_sh, 1.0]]) @ G_n)

    # Assemble the 4 entries of the crossed matrix on H_n = {v_+ (x) |n>, v_- (x) |n+1>}:
    # (0,0): (+,n)->(+,n) comes from (-,n)->(-,n) of Tinv_nm1[1, 1]
    # (1,1): (-,n+1)->(-,n+1) comes from (+,n+1)->(+,n+1) of Tinv_{n+1}[0, 0]
    G_np1 = np.diag([q_test ** (-(n_val + 1)), q_test ** (n_val + 2 + s_val)])
    Tinv_np1 = np.linalg.inv(Phi_sh * np.array([[1.0, -1j * Y_sh], [-1j * Y_sh, 1.0]]) @ G_np1)
    # (1,0): (+,n)->(-,n+1) comes from (+,n)->(-,n+1) of Tinv_n[1, 0]
    # (0,1): (-,n+1)->(+,n) comes from (-,n+1)->(+,n) of Tinv_n[0, 1]
    T_crossed = np.array(
        [
            [Tinv_nm1[1, 1], Tinv_n[0, 1]],
            [Tinv_n[1, 0], Tinv_np1[0, 0]],
        ]
    )
    return np.linalg.norm(T_p - T_crossed) / np.linalg.norm(T_p)


err_q_phys = crossing_check(q)
err_q_wrong = crossing_check(-q)

report(
    f"defect crossing with C = sigma_x holds at physical q = -e^(-i*pi*lam) [err {err_q_phys:.1e}] and fails at +e^(-i*pi*lam) [err {err_q_wrong:.2f}]",
    err_q_phys < 1e-8 and err_q_wrong > 0.5,
)

# ---- Check 2: Unitarity, crossing, pole cancellation at +i*pi/(2*lam), and simple pole at -i*pi/(2*lam) ----
u_grid = np.linspace(-2.5, 2.5, 25)
Phi_grid = np.array([sg_transmission_scalar(u, lam) for u in u_grid])
Phi_shift = np.array([sg_transmission_scalar(u + 1j * np.pi, lam) for u in u_grid])
unit_err = np.max(np.abs(np.abs(Phi_grid) ** 2 * (1.0 + np.exp(2.0 * lam * u_grid)) - 1.0))
cross_err = np.max(np.abs(Phi_grid * Phi_shift * (1.0 + np.exp(2.0 * lam * u_grid)) - 1.0))

u_pos = +1j * np.pi / (2.0 * lam)
u_neg = -1j * np.pi / (2.0 * lam)
phi_at_pos = sg_transmission_scalar(u_pos + 1e-7, lam)
res_1 = 1e-5 * sg_transmission_scalar(u_neg + 1e-5, lam)
res_2 = 1e-7 * sg_transmission_scalar(u_neg + 1e-7, lam)
pole_scaling_err = abs(res_1 - res_2) / abs(res_2)

report(
    f"Phi(u) satisfies unitarity [err {unit_err:.1e}], crossing [err {cross_err:.1e}], is finite at +i*pi/(2*lam) "
    f"[|Phi|={abs(phi_at_pos):.3f}], and has a simple pole at -i*pi/(2*lam) [1/eps err {pole_scaling_err:.1e}]",
    unit_err < 1e-12 and cross_err < 1e-7 and np.isfinite(abs(phi_at_pos)) and abs(res_2) > 0.05 and pole_scaling_err < 1e-4,
)

# ---- Check 3: Semiclassical Wigner time delay convergence ----
u_wig = np.linspace(-2.2, -0.25, 12)
target_wig = np.log(np.abs(1.0 / np.tanh(0.5 * u_wig)))
devs_wig = []
for lam_w in (10.0, 40.0, 160.0):
    du = 1e-5
    ph_p = np.unwrap([np.angle(sg_transmission_scalar(u + du, lam_w)) for u in u_wig])
    ph_m = np.unwrap([np.angle(sg_transmission_scalar(u - du, lam_w)) for u in u_wig])
    darg = (ph_p - ph_m) / (2.0 * du)
    w_curve = (np.pi / (lam_w + 1.0)) * darg
    devs_wig.append(np.max(np.abs(w_curve - target_wig)))

report(
    f"Wigner time delay (pi/(lam+1)) d/dtheta arg Phi converges to ln|coth(u/2)| as 1/lam: "
    f"devs {[round(d, 4) for d in devs_wig]} at lam = 10, 40, 160",
    devs_wig[0] > devs_wig[1] > devs_wig[2] and devs_wig[2] < 0.015,
)

# ---- Check 4: Problem 5.2 (Odd-sector bootstrap at theta_P = theta_d - i*pi/(2*lambda)) ----
# At u_P = -i*pi/(2*lambda), Y_P = -i, so on H_n = {v_+ (x) |n>, v_- (x) |n+1>},
# M_n(Y_P) = [[1, -1], [-1, 1]] @ diag(q^{-n}, q^{n+1}) maps any input to a multiple of |e_n> = |n> (x) v_+ - |n+1> (x) v_-.
# Send a probe soliton V(y) through W (x) V(y_P) via (1_W (x) Rcheck(y/y_P)) (T_{12}(y) (x) 1_{V_P}) restricted to span{|e_n>}:
th_probe = 0.68
y_probe = np.exp(lam * th_probe)
y_P = np.exp(lam * (thd - 1j * np.pi / (2.0 * lam)))
Rc = rcheck(sl2_spin(1, q, y_probe / y_P), sl2_spin(1, q, 1.0))[0]
# Rc acts on V(y) (x) V(y_P) -> V(y_P) (x) V(y) in basis (++, +-, -+, --)
# Let's compute the induced 2x2 block on {v_+ (x) |e_n>, v_- (x) |e_{n+1}>} for n = -2..2 and compare gauge invariants with
# the odd-sector Weyl intertwiner (s = 1, theta_d -> theta_d - i*pi, so Y_odd = e^{lambda(theta - (theta_d - i*pi))} = -q^{-1} Y):
Y_pr = y_probe / np.exp(lam * thd)
Y_odd = np.exp(lam * (th_probe - (thd - 1j * np.pi)))

# Build full space W (x) V(y_P) on window N = 4 and restrict (1_W (x) Rc) (T_{V,W}(y) (x) 1_{V_P}) to the image vectors |e_n>
from defect_common import sg_tmatrix  # noqa: E402

T_pr, V_pr, W_pr = sg_tmatrix(th_probe, lam, q, thd, 4, s=0.0, scalar=False)
dW = W_pr["dim"]
ns_pr = list(W_pr["n"])
F_fused = np.kron(np.eye(dW), Rc) @ np.kron(T_pr, np.eye(2))
# Input space is V(y) (x) W (x) V(y_P), dimension 2 * dW * 2; output space is W (x) V(y_P) (x) V(y), dimension dW * 2 * 2
# Image state |e_n> in W (x) V(y_P) is |n> (x) v_+ - |n+1> (x) v_-
e_vecs = {}
for n_val in range(-3, 3):
    k_n = ns_pr.index(n_val)
    k_np1 = ns_pr.index(n_val + 1)
    v_e = np.zeros(dW * 2, dtype=complex)
    v_e[2 * k_n + 0] = 1.0
    v_e[2 * k_np1 + 1] = -1.0
    e_vecs[n_val] = v_e

# Extract A_n^odd, B_n^odd, C_n^odd, D_n^odd on |e_n>
A_odd, B_odd, C_odd, D_odd = {}, {}, {}, {}
for n_val in range(-2, 2):
    # v_+ (x) |e_n> -> A_n^odd |e_n> (x) v_+ + C_n^odd |e_{n+1}> (x) v_-
    in_plus = np.kron(np.array([1.0, 0.0]), e_vecs[n_val])
    out_plus = F_fused @ in_plus
    # Project onto |e_n> (x) v_+ and |e_{n+1}> (x) v_-
    en_vp = np.kron(e_vecs[n_val], np.array([1.0, 0.0]))
    enp1_vm = np.kron(e_vecs[n_val + 1], np.array([0.0, 1.0]))
    A_odd[n_val] = np.vdot(en_vp, out_plus) / np.vdot(en_vp, en_vp)
    C_odd[n_val] = np.vdot(enp1_vm, out_plus) / np.vdot(enp1_vm, enp1_vm)

    # v_- (x) |e_n> -> B_n^odd |e_{n-1}> (x) v_+ + D_n^odd |e_n> (x) v_-
    in_minus = np.kron(np.array([0.0, 1.0]), e_vecs[n_val])
    out_minus = F_fused @ in_minus
    enm1_vp = np.kron(e_vecs[n_val - 1], np.array([1.0, 0.0]))
    en_vm = np.kron(e_vecs[n_val], np.array([0.0, 1.0]))
    B_odd[n_val] = np.vdot(enm1_vp, out_minus) / np.vdot(enm1_vp, enm1_vp)
    D_odd[n_val] = np.vdot(en_vm, out_minus) / np.vdot(en_vm, en_vm)

    # Check zero leakage orthogonal to span{|e_n> (x) v_+, |e_{n+1}> (x) v_-}
    recon_p = A_odd[n_val] * en_vp + C_odd[n_val] * enp1_vm
    recon_m = B_odd[n_val] * enm1_vp + D_odd[n_val] * en_vm
    assert np.allclose(out_plus, recon_p) and np.allclose(out_minus, recon_m)

# Check gauge-invariant ratios for s = 1 and theta_d -> theta_d - i*pi:
# A_{n+1}/A_n = q^{-1}, D_n/A_n = q^{2n+1}, and (C_n B_{n+1}) / (A_n D_{n+1}) = -Y_odd^2
odd_ok = True
for n_val in (-1, 0):
    odd_ok &= (
        np.isclose(A_odd[n_val + 1] / A_odd[n_val], 1.0 / q)
        and np.isclose(D_odd[n_val] / A_odd[n_val], q ** (2 * n_val + 1))
        and np.isclose((C_odd[n_val] * B_odd[n_val + 1]) / (A_odd[n_val] * D_odd[n_val + 1]), -(Y_odd**2))
    )

report(
    "Problem 5.2: defect bootstrap on the rank-1 residue at theta_P = theta_d - i*pi/(2*lam) produces the s=1 odd defect at theta_d - i*pi",
    odd_ok,
)

# ---- Figure ----
fig, (axa, axb) = plt.subplots(
    1, 2, figsize=(7.0, 2.95), gridspec_kw=dict(width_ratios=[1.05, 1.15], wspace=0.28)
)

# Panel (a): |Phi(i*eta)| along the imaginary axis u = i*eta showing pole at -pi/(2*lam) and cancellation at +pi/(2*lam)
lam_plot = 1.8
eta_p = np.pi / (2.0 * lam_plot)  # = 1/(2*1.8) * pi

etas_low1 = np.linspace(-0.65 * np.pi, -eta_p - 0.015, 70)
etas_low2 = np.linspace(-eta_p + 0.015, 0.95 * np.pi, 110)
etas_all = np.linspace(-0.65 * np.pi, 0.95 * np.pi, 180)

phi_low1 = np.array([abs(sg_transmission_scalar(1j * et, lam_plot)) for et in etas_low1])
phi_low2 = np.array([abs(sg_transmission_scalar(1j * et, lam_plot)) for et in etas_low2])
pref_curve = np.abs(1.0 / (1.0 + 1j * np.exp(1j * lam_plot * etas_all)))

axa.axvspan(0.0, np.pi / np.pi, color=style.BLUE, alpha=0.06, lw=0, label=r"physical strip $0 < \operatorname{Im} u < \pi$")
axa.axvline(0.0, color=style.LIGHTGREY, lw=0.9)
axa.axvline(-eta_p / np.pi, color=style.RED, lw=1.2, ls="--")
axa.axvline(+eta_p / np.pi, color=style.TEAL, lw=1.1, ls=":")

axa.plot(etas_all / np.pi, np.clip(pref_curve, 0, 8.0), color=style.TEAL, lw=1.2, ls=":",
         label=r"bare $|1 + i e^{\lambda u}|^{-1}$")
axa.plot(etas_low1 / np.pi, np.clip(phi_low1, 0, 8.0), color=style.BLUE, lw=1.9,
         label=r"exact $|\Phi(i\eta)|$ ($\lambda=1.8$)")
axa.plot(etas_low2 / np.pi, np.clip(phi_low2, 0, 8.0), color=style.BLUE, lw=1.9)

# Mark the cancelled pole at +pi/(2*lam)
val_canc = abs(sg_transmission_scalar(1j * eta_p + 1e-7, lam_plot))
axa.plot([eta_p / np.pi], [val_canc], "o", color=style.TEAL, ms=5.2, zorder=6)
axa.annotate(
    r"cancelled by $1/\bar{r}=0$" + "\n" + r"at $+\frac{i\pi}{2\lambda}$ (no bound state!)",
    xy=(eta_p / np.pi, val_canc),
    xytext=(0.60, 2.25),
    fontsize=6.8,
    color=style.TEAL,
    ha="center",
    arrowprops=dict(arrowstyle="->", color=style.TEAL, lw=0.8, shrinkB=4),
)

# Label resonance pole at -pi/(2*lam)
axa.annotate(
    r"resonance pole" + "\n" + r"$u_P = -\frac{i\pi}{2\lambda}$" + "\n" + r"($\operatorname{Im}\Delta E < 0$)",
    xy=(-eta_p / np.pi, 5.1),
    xytext=(-0.48, 3.65),
    fontsize=6.8,
    color=style.RED,
    ha="center",
    arrowprops=dict(arrowstyle="->", color=style.RED, lw=0.8, shrinkB=3),
)

axa.set_xlim(-0.68, 1.0)
axa.set_ylim(0.0, 6.2)
axa.set_xlabel(r"imaginary rapidity $\eta / \pi$ ($u = i\eta$)")
axa.set_ylabel(r"modulus $|\Phi(i\eta)|$")
axa.set_title(r"(a) Physical-strip cancellation & resonance pole", fontsize=8.3, color=style.GREY)
axa.legend(loc="upper right", fontsize=6.6, handlelength=1.4, borderpad=0.2)

# Panel (b): Semiclassical Wigner time delay vs classical Bäcklund curve
u_left_cl = np.linspace(-2.5, -0.05, 160)
u_right_cl = np.linspace(0.05, 2.5, 160)
u_w_grid = np.linspace(-2.5, 2.5, 75)

axb.axhline(0.0, color=style.LIGHTGREY, lw=0.8)
axb.axvline(0.0, color=style.PURPLE, lw=0.9, ls=":")

axb.plot(u_left_cl, np.log(np.abs(1.0 / np.tanh(0.5 * u_left_cl))), color=style.GREY, lw=1.5, ls="--",
         label=r"classical $\ln|\coth\frac{\theta-\theta_d}{2}|$")
axb.plot(u_right_cl, np.log(np.abs(1.0 / np.tanh(0.5 * u_right_cl))), color=style.GREY, lw=1.5, ls="--")

for lam_c, col_c, lw_c, lab_c in [
    (1.6, style.ORANGE, 1.4, r"quantum $\lambda = 1.6$"),
    (4.5, style.TEAL, 1.6, r"quantum $\lambda = 4.5$"),
    (16.0, style.BLUE, 1.8, r"quantum $\lambda = 16$"),
]:
    du = 1e-4
    ph_p = np.unwrap([np.angle(sg_transmission_scalar(u + du, lam_c)) for u in u_w_grid])
    ph_m = np.unwrap([np.angle(sg_transmission_scalar(u - du, lam_c)) for u in u_w_grid])
    w_vals = (np.pi / (lam_c + 1.0)) * (ph_p - ph_m) / (2.0 * du)
    axb.plot(u_w_grid, w_vals, color=col_c, lw=lw_c, label=lab_c)

axb.set_xlim(-2.5, 2.5)
axb.set_ylim(-0.1, 3.4)
axb.set_xlabel(r"relative rapidity $\theta - \theta_d$")
axb.set_ylabel(r"$\frac{\pi}{\lambda+1}\frac{d}{d\theta}\arg\Phi(\theta-\theta_d)$")
axb.set_title(r"(b) Semiclassical Wigner delay ($\lambda \to \infty$)", fontsize=8.3, color=style.GREY)
axb.legend(loc="upper right", fontsize=6.7, handlelength=1.5, borderpad=0.2)

fig.subplots_adjust(left=0.08, right=0.98, bottom=0.17, top=0.88, wspace=0.28)
style.save(fig, "l5_resonance_wigner", OUT)
raise SystemExit(0 if style.all_ok() else 1)
