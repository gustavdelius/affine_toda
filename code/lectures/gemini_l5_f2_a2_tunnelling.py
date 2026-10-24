"""Lectures on integrable defects (defects-gemini), Lecture 5, figure 2 (fig-lec5-a2-tunnelling):
Chiral Multi-Branch Quantum Tunnelling in a_2^(1) and Species-Dependent Resonance Windows
(sec-lec5-a2-tunnelling, prp-lec5-a2-matrices, eq-lec5-vbar-from-fusion, exr-lec5-omega-trick).

Panel (a): Normalized branch probabilities p_j(theta) = |X|^{2j} / (1 + |X|^2 + |X|^4) (j = 0, 1, 2)
           for anti-soliton Vbar transmission in a_2^(1) (h = 3), where |X| = exp(omega*(theta - theta_d)),
           at quantum coupling omega = 1.2 (solid) vs semiclassical omega = 4.5 (dashed).
           The j = 1 quantum tunnelling branch (absent classically) peaks at p_1(theta_d) = 1/3 and is
           exponentially suppressed as exp(-omega |theta - theta_d|) away from theta_d as omega -> oo.
Panel (b): Imaginary parts of the a_2^(1) defect resonance poles Im(theta_a - theta_d)/pi
           = 1/2 - ((3-a)/3 + 1/(3*omega)) as functions of 1/omega for species a = 1 (V) and a = 2 (Vbar),
           highlighting the lower half-plane decaying-resonance strip Im(theta_a - theta_d) in (-pi/2, 0)
           where Re(Delta E_a) > 0 and Im(Delta E_a) < 0 (unstable on 1/omega in (0, 1) for a = 1 and
           (1/2, 2) for a = 2, i.e. max(0, a - 3/2) < 1/omega < a).

Checks:
  1. Closed-form intertwiners M_V (2 branches) and M_{Vbar} (3 branches) in a_2^(1) commute with U_q(b^-)
     on all interior labels of the Weyl module W (R2 = 18).
  2. Problem 5.1: On W, the adjacent-node clock relation f_j f_{j+1} = q^-1 f_{j+1} f_j holds on interior
     labels, so Omega_j = Delta(f_j) Delta(f_{j+1}) - q^-1 Delta(f_{j+1}) Delta(f_j) has zero pure-defect
     1 (x) W^2 leg, while pi_V(f_j) pi_V(f_{j+1}) = 0 on V (and pi_{Vbar}(f_{j+1}) pi_{Vbar}(f_j) = 0 on Vbar).
  3. Bootstrap fusion of M_V(theta - i*pi/3) and M_V(theta + i*pi/3) on the unique submodule
     Vbar(y) subset V(theta - i*pi/3) (x) V(theta + i*pi/3) reproduces the 3-branch matrix M_{Vbar}(theta)
     with ratio 1 on all deep lattice labels.
  4. The j = 1 quantum tunnelling probability p_1(theta) = |X|^2 / (1 + |X|^2 + |X|^4) peaks at 1/3 at
     theta = theta_d and decays exponentially (< 1e-4 at |theta - theta_d| = 1.5 for omega = 7);
     scalar factors g(D) and gbar(D) exhibit simple poles at D_a = i*pi/2 - i*pi*((3-a)/3 + 1/(3*omega)).

Writes lecture-notes/defects-gemini/figures/l5_a2_tunnelling.svg.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from style import plt, report  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from boundary_common import an_vector, cop, intertwiners, prop_residual, q_phys  # noqa: E402
from defect_common import (  # noqa: E402
    an_scalar_constant,
    an_scalar_g,
    an_scalar_gbar,
    an_T_vector,
    an_tensor_rep,
    an_weyl_module,
)

style.setup()
OUT = style.ROOT / "lecture-notes" / "defects-gemini" / "figures"

r = 2
h = r + 1  # h = 3 for a_2^(1)
omega = np.sqrt(2.0)  # generic irrational
q = q_phys(omega)
thd = 0.25
yd = np.exp(omega * thd)
mu = -np.exp(-0.5j * np.pi * omega) / (q - 1.0 / q)

# ---- Check 1: M_V (2 branches) and M_{Vbar} (3 branches) intertwine U_q(b^-) on interior labels ----
W = an_weyl_module(r, -1, q, yd, R2=18, mu=mu, orient=1)
dW = W["dim"]
inter_mask = W["interior"]
rows_ok = np.repeat(inter_mask, h)
cols_ok = np.tile(inter_mask, h)
blk = np.ix_(rows_ok, cols_ok)

th_test = 0.62
y_test = np.exp(omega * th_test)
X_test = np.exp(omega * (th_test - thd - 0.5j * np.pi))

V_rep = an_vector(r, q, y_test, conj=False)
Vbar_rep = an_vector(r, q, y_test, conj=True)
T_V = an_T_vector(W, X_test, conj=False)
T_Vbar = an_T_vector(W, X_test, conj=True)

err_inter_V = 0.0
err_inter_Vbar = 0.0
for j in range(h):
    for kd in ("k", "f"):
        err_inter_V = max(
            err_inter_V,
            np.max(np.abs((T_V @ cop(V_rep, W, kd, j) - cop(W, V_rep, kd, j) @ T_V)[blk])),
        )
        err_inter_Vbar = max(
            err_inter_Vbar,
            np.max(np.abs((T_Vbar @ cop(Vbar_rep, W, kd, j) - cop(W, Vbar_rep, kd, j) @ T_Vbar)[blk])),
        )

report(
    f"a_2^(1) intertwiners M_V (2 branches) and M_Vbar (3 branches) commute with U_q(b^-) on interior labels "
    f"[err_V {err_inter_V:.1e}, err_Vbar {err_inter_Vbar:.1e}]",
    err_inter_V < 1e-12 and err_inter_Vbar < 1e-12,
)

# ---- Check 2: Problem 5.1 (Adjacent-node clock relation & quadratic Omega_j trick on V (x) W and Vbar (x) W) ----
blk_W = np.ix_(inter_mask, inter_mask)
clock_err = 0.0
omega_trick_err = 0.0
v_adj_norm = 0.0
vbar_adj_norm = 0.0

for j in range(h):
    jp1 = (j + 1) % h
    fj_W = W[("f", j)]
    fjp1_W = W[("f", jp1)]
    clock_err = max(clock_err, np.max(np.abs((fj_W @ fjp1_W - (1.0 / q) * (fjp1_W @ fj_W))[blk_W])))

    fj_V = V_rep[("f", j)]
    fjp1_V = V_rep[("f", jp1)]
    kj_V = V_rep[("k", j)]
    kjp1_V = V_rep[("k", jp1)]
    v_adj_norm = max(v_adj_norm, np.linalg.norm(fj_V @ fjp1_V))
    vbar_adj_norm = max(vbar_adj_norm, np.linalg.norm(Vbar_rep[("f", jp1)] @ Vbar_rep[("f", j)]))

    # In Part VI coproduct Delta(f_j) = f_j (x) k_j^-1 + 1 (x) f_j:
    # Delta(f_j) Delta(f_{j+1}) - q^-1 Delta(f_{j+1}) Delta(f_j) kills the pure 1 (x) (f_j f_{j+1} - q^-1 f_{j+1} f_j) leg!
    Df_j = cop(V_rep, W, "f", j)
    Df_jp1 = cop(V_rep, W, "f", jp1)
    Omega_j = Df_j @ Df_jp1 - (1.0 / q) * (Df_jp1 @ Df_j)

    # Exact expansion with pure-W term removed (and fj_V @ fjp1_V = 0 on V):
    kinv_j_W = np.linalg.inv(W[("k", j)])
    kinv_jp1_W = np.linalg.inv(W[("k", jp1)])
    rhs_Omega = (
        -(1.0 / q) * np.kron(fjp1_V @ fj_V, kinv_jp1_W @ kinv_j_W)
        + np.kron(fj_V, kinv_j_W @ fjp1_W)
        - (1.0 / q) * np.kron(fj_V, fjp1_W @ kinv_j_W)
        + np.kron(fjp1_V, fj_W @ kinv_jp1_W)
        - (1.0 / q) * np.kron(fjp1_V, kinv_jp1_W @ fj_W)
    )
    omega_trick_err = max(omega_trick_err, np.max(np.abs((Omega_j - rhs_Omega)[blk])))

report(
    f"Problem 5.1: W satisfies f_j f_{{j+1}} = q^-1 f_{{j+1}} f_j [err {clock_err:.1e}], pi_V(f_j)pi_V(f_{{j+1}}) = 0 "
    f"[{v_adj_norm:.1e}], and Omega_j kills the pure defect leg [err {omega_trick_err:.1e}]",
    clock_err < 1e-12 and v_adj_norm < 1e-14 and vbar_adj_norm < 1e-14 and omega_trick_err < 1e-12,
)

# ---- Check 3: Bootstrap fusion Vbar(y) subset V(theta - i*pi/3) (x) V(theta + i*pi/3) reproduces M_{Vbar}(theta) ----
shifts = [(-1.0 / 3.0), (+1.0 / 3.0)]
V1 = an_vector(r, q, y_test * np.exp(1j * np.pi * omega * shifts[0]))
V2 = an_vector(r, q, y_test * np.exp(1j * np.pi * omega * shifts[1]))
P12 = an_tensor_rep(V1, V2)
null_phi, _ = intertwiners([(Vbar_rep[(kd, j)], P12[(kd, j)]) for j in range(h) for kd in "efk"])
Phi_emb = null_phi[0]

Ts_shift = [an_T_vector(W, np.exp(omega * (th_test + 1j * np.pi * sh - thd - 0.5j * np.pi))) for sh in shifts]
deep_idx = np.where(W["deep"])[0]
fusion_res = 0.0
fusion_ratios = []
for i in range(h):
    for k in deep_idx:
        v = np.kron(Phi_emb[:, i], np.eye(dW)[k])
        for p in range(r - 1, -1, -1):
            v = np.einsum("ij,ajb->aib", Ts_shift[p], v.reshape(h**p, h * dW, h ** (r - 1 - p))).reshape(-1)
        col = T_Vbar[:, i * dW + k]
        target = sum(
            col[k2 * h + j] * np.kron(np.eye(dW)[k2], Phi_emb[:, j])
            for k2 in range(dW)
            for j in range(h)
            if col[k2 * h + j] != 0
        )
        res_ik, c_ik = prop_residual(v, target)
        fusion_res = max(fusion_res, res_ik)
        fusion_ratios.append(c_ik)

ratio_err = max(abs(c - 1.0) for c in fusion_ratios)
report(
    f"bootstrap fusion on Vbar subset V(theta-i*pi/3) (x) V(theta+i*pi/3) (dim={len(null_phi)}) reproduces 3-branch "
    f"M_Vbar(theta) on all {len(deep_idx)} deep labels [res {fusion_res:.1e}, ratio-1 err {ratio_err:.1e}]",
    len(null_phi) == 1 and fusion_res < 1e-12 and ratio_err < 1e-12,
)


# ---- Check 4: Exponential suppression of j=1 tunnelling branch & resonance poles of g, gbar ----
def branch_probs(u, om):
    X2 = np.exp(2.0 * om * u)
    denom = 1.0 + X2 + X2**2
    return 1.0 / denom, X2 / denom, (X2**2) / denom


_, p1_peak, _ = branch_probs(0.0, omega)
_, p1_tail, _ = branch_probs(1.5, 7.0)

# Pole check for g(D) and gbar(D) at D_a = i*pi/2 - i*pi*((3-a)/3 + 1/(3*omega)) (eq-dt-pole, Delta = theta_d + i*pi/2)
c_norm = an_scalar_constant(r, omega, q)
D1_pole = 0.5j * np.pi - 1j * np.pi * (2.0 / 3.0 + 1.0 / (3.0 * omega))
D2_pole = 0.5j * np.pi - 1j * np.pi * (1.0 / 3.0 + 1.0 / (3.0 * omega))
e1, e2 = 1e-4, 1e-5
pole1_err = abs(e1 * an_scalar_g(D1_pole + e1, r, omega, c=c_norm) / (e2 * an_scalar_g(D1_pole + e2, r, omega, c=c_norm)) - 1.0)
pole2_err = abs(e1 * an_scalar_gbar(D2_pole + e1, r, omega, c=c_norm) / (e2 * an_scalar_gbar(D2_pole + e2, r, omega, c=c_norm)) - 1.0)

report(
    f"j=1 tunnelling branch peaks at p_1(0) = {p1_peak:.4f} (= 1/3), decays to {p1_tail:.1e} at |u|=1.5 (omega=7), "
    f"and g, gbar have simple poles at D_1, D_2 [1/eps errs {pole1_err:.1e}, {pole2_err:.1e}]",
    abs(p1_peak - 1.0 / 3.0) < 1e-14 and p1_tail < 1e-4 and pole1_err < 1e-3 and pole2_err < 1e-3,
)

# ---- Figure ----
fig, (axa, axb) = plt.subplots(
    1, 2, figsize=(7.0, 2.95), gridspec_kw=dict(width_ratios=[1.12, 1.08], wspace=0.28)
)

# Panel (a): Relative branch probabilities p_0, p_1, p_2 of Vbar in a_2^(1)
u_vals = np.linspace(-2.2, 2.2, 300)
p0_q, p1_q, p2_q = branch_probs(u_vals, 1.2)
p0_sc, p1_sc, p2_sc = branch_probs(u_vals, 4.5)

axa.axhline(0.0, color=style.LIGHTGREY, lw=0.8)
axa.axhline(1.0, color=style.LIGHTGREY, lw=0.8)
axa.axvline(0.0, color=style.PURPLE, lw=0.9, ls=":")

# Shade the quantum tunnelling branch j=1 for omega = 1.2
axa.fill_between(u_vals, 0.0, p1_q, color=style.RED, alpha=0.14)

axa.plot(u_vals, p0_q, color=style.BLUE, lw=1.7, label=r"$j=0$ (diagonal $\bar{v}_i \to \bar{v}_i$)")
axa.plot(u_vals, p2_q, color=style.TEAL, lw=1.7, label=r"$j=2$ (classical convert $\bar{v}_i \to \bar{v}_{i+2}$)")
axa.plot(u_vals, p1_q, color=style.RED, lw=2.0, label=r"$j=1$ (quantum tunnelling $\bar{v}_i \to \bar{v}_{i+1}$)")

# Semiclassical omega = 4.5 curves in dashed
axa.plot(u_vals, p0_sc, color=style.BLUE, lw=1.2, ls="--")
axa.plot(u_vals, p2_sc, color=style.TEAL, lw=1.2, ls="--")
axa.plot(u_vals, p1_sc, color=style.RED, lw=1.35, ls="--", label=r"$j=0,1,2$ at $\omega = 4.5$ (dashed)")

axa.plot([0.0], [1.0 / 3.0], "o", color=style.RED, ms=4.8, zorder=6)
axa.annotate(
    r"quantum tunnelling" + "\n" + r"$p_1(0) = \frac{1}{3}$" + "\n" + r"($\propto e^{-\omega|\theta-\theta_d|}$)",
    xy=(0.0, 1.0 / 3.0),
    xytext=(-1.25, 0.44),
    fontsize=6.8,
    color=style.RED,
    ha="center",
    arrowprops=dict(arrowstyle="->", color=style.RED, lw=0.8, shrinkB=4),
)

axa.set_xlim(-2.2, 2.2)
axa.set_ylim(-0.04, 1.40)
axa.set_xlabel(r"relative rapidity $\theta - \theta_d$")
axa.set_ylabel(r"branch weight $p_j(\theta) = |T_{\bar{V},j}|^2 / \sum_k |T_{\bar{V},k}|^2$")
axa.set_title(r"(a) $\bar{V}$ three-branch weights ($\omega=1.2$ vs $4.5$)", fontsize=8.3, color=style.GREY)
axa.legend(loc="upper left", fontsize=6.1, handlelength=1.4, borderpad=0.18, labelspacing=0.18)

# Panel (b): Resonance pole trajectories Im(theta_a - theta_d)/pi vs 1/omega in a_2^(1) (eq-dt-pole, eq-dt-window)
inv_om = np.linspace(0.0, 2.4, 240)
# Exact trajectories from eq-dt-pole: Im(theta_a - theta_d)/pi = 1/2 - (3-a)/3 - 1/(3*omega)
# a = 1 (V):    Im(theta_1 - theta_d)/pi = -1/6 - 1/(3*omega) -> in (-1/2, 0) for 1/omega in (0, 1)
# a = 2 (Vbar): Im(theta_2 - theta_d)/pi = +1/6 - 1/(3*omega) -> in (-1/2, 0) for 1/omega in (1/2, 2)
y_a1 = -1.0 / 6.0 - inv_om / 3.0
y_a2 = +1.0 / 6.0 - inv_om / 3.0

axb.axhspan(-0.5, 0.0, color=style.RED, alpha=0.10, lw=0, label=r"decaying resonance ($\operatorname{Re}\Delta E > 0, \operatorname{Im}\Delta E < 0$)")
axb.axhspan(0.0, 0.35, color=style.BLUE, alpha=0.06, lw=0, label=r"physical strip ($\operatorname{Im} u > 0$)")

axb.axhline(0.0, color=style.GREY, lw=0.9)
axb.axhline(-0.5, color=style.RED, lw=0.9, ls="--")

# Plot a=1 (V) trajectory and highlight unstable window 1/omega in (0, 1)
axb.plot(inv_om, y_a1, color=style.BLUE, lw=1.5, ls="--", label=r"$a=1$ ($V$): $-\frac{1}{6} - \frac{1}{3\omega}$")
mask_a1 = inv_om <= 1.0
axb.plot(inv_om[mask_a1], y_a1[mask_a1], color=style.BLUE, lw=3.0)
axb.plot([0.0, 1.0], [-1.0 / 6.0, -0.5], "o", color=style.BLUE, ms=4.5, zorder=6)

# Plot a=2 (Vbar) trajectory and highlight unstable window 1/omega in (0.5, 2.0)
axb.plot(inv_om, y_a2, color=style.RED, lw=1.5, ls="--", label=r"$a=2$ ($\bar{V}$): $+\frac{1}{6} - \frac{1}{3\omega}$")
mask_a2 = (inv_om >= 0.5) & (inv_om <= 2.0)
axb.plot(inv_om[mask_a2], y_a2[mask_a2], color=style.RED, lw=3.0)
axb.plot([0.5, 2.0], [0.0, -0.5], "o", color=style.RED, ms=4.5, zorder=6)

axb.annotate(
    r"$V$ unstable $\frac{1}{\omega} \in (0, 1)$",
    xy=(0.30, -0.27),
    xytext=(0.56, -0.23),
    fontsize=6.5,
    color=style.BLUE,
    ha="left",
    arrowprops=dict(arrowstyle="->", color=style.BLUE, lw=0.8, shrinkB=3),
)
axb.annotate(
    r"$\bar{V}$ unstable $\frac{1}{\omega} \in (\frac{1}{2}, 2)$",
    xy=(1.25, -0.25),
    xytext=(1.55, -0.10),
    fontsize=6.5,
    color=style.RED,
    ha="center",
    arrowprops=dict(arrowstyle="->", color=style.RED, lw=0.8, shrinkB=3),
)

axb.set_xlim(0.0, 2.4)
axb.set_ylim(-1.12, 0.35)
axb.set_xlabel(r"inverse coupling $1 / \omega$")
axb.set_ylabel(r"pole location $\operatorname{Im}(\theta_a - \theta_d) / \pi$")
axb.set_title(r"(b) $a_2^{(1)}$ resonance windows vs $1/\omega$", fontsize=8.3, color=style.GREY)
axb.legend(loc="lower left", fontsize=6.0, handlelength=1.4, borderpad=0.2, labelspacing=0.2)

fig.subplots_adjust(left=0.08, right=0.98, bottom=0.17, top=0.88, wspace=0.28)
style.save(fig, "l5_a2_tunnelling", OUT)
raise SystemExit(0 if style.all_ok() else 1)
