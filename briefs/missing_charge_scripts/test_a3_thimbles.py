"""Analyze a_3^(1) exact partition functions, critical points, and continuum profiles in reachable (1,3) vs missing (1,2) sectors."""
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../textbook_code/thimbles')))

import numpy as np
import mpmath as mp
from lie import An
from model import Chain
from newton_batch import find_cps

g = An(3)  # a_3^(1), h=4, rank=3

# 1. Exact analytical verification of Hirota bilinear equation for
#    tau_j(x, t) = 1 + i^j C1 e^{x - i t} + (-1)^j C4 e^{2x}
C1 = 1.5 + 0.8j
C4 = 0.7 - 1.1j
xs = np.linspace(-3, 3, 11)[:, None]
ts = np.linspace(-2, 2, 11)[None, :]
j = np.arange(4)[:, None, None]
A = C1 * np.exp(xs - 1j * ts)
B = C4 * np.exp(2.0 * xs)

tau = 1.0 + (1j ** j) * A + ((-1.0) ** j) * B
tau_x = (1j ** j) * A + 2.0 * ((-1.0) ** j) * B
tau_xx = (1j ** j) * A + 4.0 * ((-1.0) ** j) * B
tau_t = -1j * (1j ** j) * A
tau_tt = -(1j ** j) * A

# Bilinear equation: tau * (tau_xx - tau_tt) - (tau_x^2 - tau_t^2) = tau^2 - tau_{j-1} tau_{j+1}
lhs_bilin = tau * (tau_xx - tau_tt) - (tau_x ** 2 - tau_t ** 2)
rhs_bilin = tau ** 2 - np.roll(tau, 1, axis=0) * np.roll(tau, -1, axis=0)
print("1a. Exact Hirota bilinear residual:", np.max(np.abs(lhs_bilin - rhs_bilin)))

# Exact PDE residual for u = i sum_j alpha_j log tau_j:
d2_log_tau = (tau * (tau_xx - tau_tt) - (tau_x ** 2 - tau_t ** 2)) / (tau ** 2)
uxx_minus_utt = np.einsum('ja,jxy->axy', 1j * g.alpha, d2_log_tau)
log_tau = np.log(tau)
u = np.einsum('ja,jxy->axy', 1j * g.alpha, log_tau)
exp_ia_u = np.exp(np.einsum('ja,axy->jxy', 1j * g.alpha, u))
gradW = np.einsum('ja,jxy->axy', -1j * g.alpha, exp_ia_u)
print("1b. Exact analytical PDE residual |u_xx - u_tt - grad W(u)|:", np.max(np.abs(uxx_minus_utt - gradW)))

# 2. Exact Coulomb-gas partition function Z_1(lambda) for all 6 weights of Lambda^2(C^4) in a_3^(1)
def Z_L1_coulomb_general(g_alg, kappa, lam, hbar, M=14):
    r = g_alg.r
    h = g_alg.h
    Gram = g_alg.alpha @ g_alg.alpha.T
    al_lam = np.round(g_alg.alpha @ lam).astype(int)
    k = kappa / hbar
    fac = np.array([float(mp.factorial(i)) for i in range(M + 1)])
    pw = np.array([k ** i / fac[i] for i in range(M + 1)])
    grid = np.meshgrid(*([np.arange(M + 1)] * h), indexing='ij')
    m = np.stack([x.ravel() for x in grid], axis=-1)
    q2 = np.einsum('ni,ij,nj->n', m, Gram, m)
    parity = (m @ al_lam) % 2
    signs = np.where(parity == 0, 1.0, -1.0)
    weights = np.prod(pw[m], axis=-1)
    s = np.sum(weights * signs * np.exp(-hbar * q2 / 4.0))
    lam2 = float(lam @ lam)
    pref = (2 * np.pi * hbar) ** (-r) * (np.pi * hbar) ** (r / 2.0) * np.exp(-np.pi ** 2 * lam2 / hbar - h * k)
    return pref * s

w_reps = g.fundamental_weights_rep(2)
kappa = 1.0
for hbar in (1.0, 0.5, 0.25):
    print(f"\n2. Exact Z_1(lambda) for all 6 weights of Lambda^2 in a_3^(1) (kappa={kappa}, hbar={hbar}):")
    for key, lam in sorted(w_reps.items()):
        Z_ex = Z_L1_coulomb_general(g, kappa, lam, hbar, M=14)
        print(f"   sector {key}: Z_1 = {Z_ex:.12e}")

# 3. Analyze critical points and saddle expansion for L=1 in sector (1,3) and missing sector (1,2)
for sec_name, key in [("reachable (1,3)", (1, 3)), ("missing (1,2)", (1, 2))]:
    lam = w_reps[key]
    ch = Chain(g, 1, kappa, lam)
    cps, cnt = find_cps(ch, nstart=15000, re_spread=3.5, im_spread=1.5)
    S = np.array([ch.action(u) for u in cps])
    o = np.argsort(S.real)
    print(f"\n3. Sector {sec_name} (L=1, kappa={kappa}): found {len(cps)} critical points")
    for idx in o[:8]:
        U = cps[idx]
        mu, V = ch.tangent(U)
        Z1_05 = ch.one_loop(U, 0.5, V, mu)
        c1 = ch.two_loop_c1(U)
        Z2_05 = Z1_05 * (1.0 + 0.5 * c1)
        ue = g.eps_coords(U[0]) / (2 * np.pi)
        print(f"   S = {S[idx].real:9.5f} {S[idx].imag:+9.5f}j | u1/(2pi) = {np.round(ue, 4)} | Z2(hbar=0.5) = {Z2_05.real:+.6e}{Z2_05.imag:+.2e}j")
