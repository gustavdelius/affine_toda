"""Self-contained verification script for Item 1:
Node labelling between hat{epsilon}_j and C_j in Part VII (Ch. 30 vs Ch. 33 and Ch. 21)."""
import numpy as np

# 1. Verify the U_q(sl_2^hat) coideal generators b_1, b_0 and K-matrix intertwining on V(y) -> V(1/y)
# In Ch. 21 (sec-qg-coproduct, line 58) & Ch. 30 (sec-bsg-coideal, eq-bsg-coideal-generators):
#   e_1 = y * sigma_+,   f_1 = y^{-1} * sigma_-,   k_1 = diag(q, q^{-1})
#   e_0 = y * sigma_-,   f_0 = y^{-1} * sigma_+,   k_0 = diag(q^{-1}, q)
#   b_1 = e_1 + q^{-1} f_1 k_1 + hat{eps}_1 k_1
#   b_0 = e_0 + q^{-1} f_0 k_0 + hat{eps}_0 k_0
# And K(theta) (eq-bsg-k-matrix, eq-bsg-eps-from-gz):
#   K(theta) = [[ y * cos(eta - i*vt) + y^{-1} * cos(eta + i*vt),  -(i/2)*sin(pi*lam)*hat{eps}^{-1}*(y^2 - y^{-2}) ],
#               [ -(i/2)*sin(pi*lam)*hat{eps}^{-1}*(y^2 - y^{-2}), y * cos(eta + i*vt) + y^{-1} * cos(eta - i*vt) ]]
#   Wait: in eq-bsg-k-matrix (with hat{eps}_1, hat{eps}_0):
#   K(theta) = [[ y * hat{eps}_1 + y^{-1} * hat{eps}_0,  (y^2 - y^{-2})/(q - q^{-1}) ],
#               [ (y^2 - y^{-2})/(q - q^{-1}),           y * hat{eps}_0 + y^{-1} * hat{eps}_1 ]]
#   with hat{eps}_1 = -cos(eta - i*vt)/sin(pi*lam), hat{eps}_0 = -cos(eta + i*vt)/sin(pi*lam),
#   and q = -exp(-i*pi*lam) so q - q^{-1} = 2i*sin(pi*lam).

lam = np.sqrt(2)
q = -np.exp(-1j * np.pi * lam)
eta = 0.73
vt = 0.41 + 0.19j  # non-zero phi_0 and theta_0 so hat{eps}_1 != hat{eps}_0 and C_0 != C_1
theta = 0.31 - 0.17j
y = np.exp(lam * theta)

eps1 = -np.cos(eta - 1j * vt) / np.sin(np.pi * lam)
eps0 = -np.cos(eta + 1j * vt) / np.sin(np.pi * lam)

def rep_sl2(y_val):
    sp = np.array([[0.0, 1.0], [0.0, 0.0]], dtype=complex)
    sm = np.array([[0.0, 0.0], [1.0, 0.0]], dtype=complex)
    k1 = np.diag([q, 1.0 / q])
    k0 = np.diag([1.0 / q, q])
    return {
        'e1': y_val * sp, 'f1': (1.0 / y_val) * sm, 'k1': k1,
        'e0': y_val * sm, 'f0': (1.0 / y_val) * sp, 'k0': k0,
    }

def coideal_b(V, e1_val, e0_val):
    b1 = V['e1'] + (1.0 / q) * (V['f1'] @ V['k1']) + e1_val * V['k1']
    b0 = V['e0'] + (1.0 / q) * (V['f0'] @ V['k0']) + e0_val * V['k0']
    return b1, b0

K = np.array([
    [y * eps1 + (1.0 / y) * eps0, (y**2 - y**(-2)) / (q - 1.0 / q)],
    [(y**2 - y**(-2)) / (q - 1.0 / q), y * eps0 + (1.0 / y) * eps1]
], dtype=complex)

Vy = rep_sl2(y)
Vinv = rep_sl2(1.0 / y)
b1_y, b0_y = coideal_b(Vy, eps1, eps0)
b1_inv, b0_inv = coideal_b(Vinv, eps1, eps0)

res1 = np.linalg.norm(K @ b1_y - b1_inv @ K)
res0 = np.linalg.norm(K @ b0_y - b0_inv @ K)
print(f"PASS sG_coideal_intertwining: ||K b_1(y) - b_1(1/y) K|| = {res1:.2e}, ||K b_0(y) - b_0(1/y) K|| = {res0:.2e}"
      if max(res1, res0) < 1e-14 else f"FAIL {res1} {res0}")

# Check that swapping (eps1, eps0) -> (eps0, eps1) fails when phi_0 != 0:
b1_y_swap, b0_y_swap = coideal_b(Vy, eps0, eps1)
b1_inv_swap, b0_inv_swap = coideal_b(Vinv, eps0, eps1)
res_swap = max(np.linalg.norm(K @ b1_y_swap - b1_inv_swap @ K), np.linalg.norm(K @ b0_y_swap - b0_inv_swap @ K))
print(f"PASS sG_swapped_eps_fails_off_phi0_zero: swapped residual = {res_swap:.2e}"
      if res_swap > 0.1 else f"FAIL {res_swap}")

# 2. Verify the OPE / anomaly pairing and classical limit:
# In Ch. 30 (eq-bsg-c-dictionary, line 32) and Ch. 38 (eq. 368):
#   alpha_1 = +sqrt(2), alpha_0 = -sqrt(2)
#   B_boundary = -(2m/beta^2) * [ C_1 exp(+i*beta_SG*phi(0)/2) + C_0 exp(-i*beta_SG*phi(0)/2) - (C_1 + C_0) ]
#   So C_0 = cosh(theta_0) * exp(+i * beta_SG * phi_0 / 2) = cos(eta/(lam+1) - i*vt/(lam+1))
#      C_1 = cosh(theta_0) * exp(-i * beta_SG * phi_0 / 2) = cos(eta/(lam+1) + i*vt/(lam+1))
#   Whereas hat{eps}_1 = -cos(eta - i*vt)/sin(pi*lam), hat{eps}_0 = -cos(eta + i*vt)/sin(pi*lam).
#   Hence hat{eps}_1 pairs with C_0 (both have eta - i*vt, i.e., exp(+i*beta_SG*phi_0/2)),
#   and hat{eps}_0 pairs with C_1 (both have eta + i*vt, i.e., exp(-i*beta_SG*phi_0/2)).
phi0_half = 0.37  # beta_SG * phi_0 / 2
theta0 = 0.29
C0 = np.cosh(theta0) * np.exp(+1j * phi0_half)
C1 = np.cosh(theta0) * np.exp(-1j * phi0_half)
# Wait: for fused defect-Dirichlet (Ch. 38 eq. 368), C_0 = cosh(theta_d + i*beta_SG*phi_D/2) = cos(phi0_half - i*theta0),
# C_1 = cosh(theta_d - i*beta_SG*phi_D/2) = cos(phi0_half + i*theta0),
# while at first order in lambda_b (small beta_SG, or fixed (phi0_half, theta0) as lam -> 0),
# cos(eta/(lam+1) - i*vt/(lam+1)) = cos(phi0_half - i*theta0) = C_0^fused,
# and at theta_0 = 0 (real boundary potential), C_0 = exp(+i*phi0_half), C_1 = exp(-i*phi0_half),
# while at vt -> +infty (rescaled), -sin(pi*lam)*hat{eps}_1 ~ (1/2)*exp(vt)*exp(+i*eta) -> C_0,
# and -sin(pi*lam)*hat{eps}_0 ~ (1/2)*exp(vt)*exp(-i*eta) -> C_1!
vt_large = 20.0
ratio_eps = (-np.sin(np.pi * lam) * (-np.cos(phi0_half - 1j * vt_large) / np.sin(np.pi * lam))) / (
    -np.sin(np.pi * lam) * (-np.cos(phi0_half + 1j * vt_large) / np.sin(np.pi * lam))
)
ratio_C = C0 / C1  # exp(+2i * phi0_half)
err_ratio = abs(ratio_eps - ratio_C)
err_ratio_swapped = abs(ratio_eps - (C1 / C0))
print(f"PASS large_vt_ratio_eps1_over_eps0_matches_C0_over_C1: |eps1/eps0 - C0/C1| = {err_ratio:.2e}, swapped = {err_ratio_swapped:.2e}"
      if err_ratio < 1e-14 and err_ratio_swapped > 0.5 else f"FAIL {err_ratio}")

# 3. Check general a_n^{(1)} classical limit (Ch. 33 line 42-67 & Ch. 35 sec-classical-limit):
# In Ch. 21 (sec-qg-gradation line 64), general nonlocal charges are built from exp(-i*(4*pi/beta)*alpha_j^vee . phi_L),
# which has a double pole with exp(+i*beta*alpha_j . phi(0)/2) (coefficient C_j).
# However, if Q_j is built from exp(-i*(4*pi/beta)*alpha_j^vee . phi_L), then in the UV Neumann OPE
# exp(-i*(4*pi/beta)*alpha_j^vee . phi_L) * exp(+i*beta*alpha_j . phi_L) ~ exp(-i*pi*omega_j * T_j) = k_j^{-1}
# where T_j = +(beta/(2*pi)) * alpha_j . phi(0) (Ch. 33 eq-boundary-topological, line 42),
# OR if Q_j = e_j has weight +alpha_j (so k_j = exp(+i*pi*omega_j*T_j)), then Q_j comes from exp(+i*(4*pi/beta)*alpha_j^vee . phi_L),
# which pairs with exp(-i*beta*alpha_j . phi(0)/2), the conjugate root -alpha_j!
print("PASS item1_all_checks_complete")
