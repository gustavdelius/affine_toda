"""Check L=1 saddle decomposition in missing sector (1,2) vs reachable (1,3), and solve the continuum BVP in V_(1,3) and V_(1,2)."""
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../textbook_code/thimbles')))

import numpy as np
import mpmath as mp
from scipy.optimize import root
from lie import An
from model import Chain
from newton_batch import find_cps
from flows import intersection_number
from flows2 import realmin_ReS

g = An(3)
w_reps = g.fundamental_weights_rep(2)

def Z_L1_coulomb_general(g_alg, kappa, lam, hbar, M=16):
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

# 1. Complete L=1 Picard-Lefschetz check for missing (1,2)
kappa = 1.0
lam = w_reps[(1, 2)]
ch0 = Chain(g, 1, kappa, lam)
rmin = realmin_ReS(ch0, nstart=60)
cps, cnt = find_cps(ch0, nstart=15000, re_spread=3.5, im_spread=1.5)
S = np.array([ch0.action(u) for u in cps])
o = np.argsort(S.real)
cps = [cps[i] for i in o]
S = S[o]
print(f"=== 1. Sector missing (1,2): L=1, kappa={kappa}, min_R Re(S) = {rmin:.6f}, Imax = {ch0.Imax:.4f} ===")
contrib = []
for idx, (U, s) in enumerate(zip(cps[:8], S[:8])):
    if s.real < rmin - 1e-7 or s.real > rmin + 6.0:
        n, sols = 0, []
        mu, V = ch0.tangent(U)
    else:
        n, sols, mu, V = intersection_number(ch0, U, ndir=600, seed=idx + 1, delta=0.01, tmax=8.0)
    ue = g.eps_coords(U[0]) / (2 * np.pi)
    print(f"  cp {idx:2d}: S={s.real:9.5f}{s.imag:+9.5f}j  n={n:+d}  #hits={len(sols)}  mu={np.round(mu, 3)}  u1/2pi={np.round(ue, 4)}")
    if n != 0:
        contrib.append((U, s, n, mu, V))

print("  --- Saddle expansion vs exact Coulomb-gas Z_1 for missing (1,2) ---")
for hbar in (1.0, 0.5, 0.25, 0.15):
    Zex = Z_L1_coulomb_general(g, kappa, lam, hbar, M=16 if hbar >= 0.25 else 20)
    tot1 = 0j
    tot2 = 0j
    for (U, s, n, mu, V) in contrib:
        Z1 = ch0.one_loop(U, hbar, V, mu)
        c1 = ch0.two_loop_c1(U)
        Z2 = Z1 * (1.0 + hbar * c1)
        tot1 += n * Z1
        tot2 += n * Z2
    rel1 = min(abs(tot1.real - Zex), abs(-tot1.real - Zex)) / abs(Zex)
    rel2 = min(abs(tot2.real - Zex), abs(-tot2.real - Zex)) / abs(Zex)
    print(f"  hbar={hbar:4.2f}: Z_exact={Zex:.8e} | sum n Z_1={abs(tot1.real):+.8e} (rel.err {rel1:.2e}) | sum n Z_2={abs(tot2.real):+.8e} (rel.err {rel2:.2e})")

# 2. Continuum scaling with vacuum boundary condition gamma(L+1) = pi in V_(1,3)
print("\n=== 2. Continuum scaling in V_(1,3) with gamma(L+1) = pi ===")
ell = 16.0
for L in (4, 8, 16, 32, 64, 128, 256):
    a = ell / (L + 1)
    kap = a ** 2
    xgrid = np.linspace(0, ell, L + 2)[1:-1] - ell / 2.0
    gam0 = 2.0 * np.arctan(np.exp(2.0 * xgrid))
    def grad_13(gam):
        g_ext = np.concatenate([[0.0], gam, [np.pi]])
        return (4.0 * (2 * gam - g_ext[:-2] - g_ext[2:]) + kap * 8.0 * np.sin(2.0 * gam)) / kap
    res13 = root(grad_13, gam0, tol=1e-12)
    gam = res13.x
    g_ext = np.concatenate([[0.0], gam, [np.pi]])
    S13 = np.sum(2.0 * np.diff(g_ext) ** 2) + kap * np.sum(4.0 * (1.0 - np.cos(2.0 * gam)))
    print(f"  L={L:3d} (a={a:.4f}): S_cont(1,3)={S13 / a:10.6f} (exact M_2 = 16.0)")
