"""a_2^(1), L=1: critical points, n_sigma, and decomposition check Z_exact = sum n_sigma Z_sigma."""
import numpy as np, time, sys, json
from lie import An
from model import Chain
from newton_batch import find_cps
from flows2 import intersection_number_n2_refined, realmin_ReS, floor_of
from thimble_exact import thimble_integral
from exactZ import Z_grid
g = An(2)
kappa = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
sector = sys.argv[2] if len(sys.argv) > 2 else 'w1'
lams = {'w1': g.fundamental_weights_rep(1)[(1,)], '0': np.zeros(2), 'root': g.alpha[1],
        'sym2': 2*g.fundamental_weights_rep(1)[(1,)]}
lam = lams[sector]
ch0 = Chain(g, 1, kappa, lam)
rmin = realmin_ReS(ch0, nstart=40)
cps, cnt = find_cps(ch0, nstart=20000, re_spread=4.5, im_spread=2.0)
S = np.array([ch0.action(u) for u in cps]); o = np.argsort(S.real); cps = [cps[i] for i in o]; S = S[o]
print("kappa=%g sector=%s lam=%s  #crit found=%d  min_R ReS=%.5f  Imax=%.2f" % (kappa, sector, np.round(lam, 4), len(cps), rmin, ch0.Imax))
out = {}
for eps in (0.0, 0.02, -0.02):
    ch = Chain(g, 1, kappa, lam, eps=eps)
    fl = floor_of(ch, rmin)
    rows = []
    for U, s in zip(cps, S):
        if abs(s.imag) > ch0.Imax + 0.5 or s.real > rmin + 40: continue
        if s.real < rmin - 1e-9 and eps == 0.0:
            n, sols = 0, []
            mu, V = ch.tangent(U)
        else:
            n, sols, mu, V, h0 = intersection_number_n2_refined(ch, U, fl)
        rows.append((U, s, n, mu, V, sols))
    out[eps] = rows
    print("eps=%+.2f" % eps)
    for (U, s, n, mu, V, sols) in rows:
        ue = g.eps_coords(U[0]) / (2*np.pi)
        print("  S=%9.5f%+9.5fi  n=%+d  #hits=%d  u1/2pi(eps-basis)=%s" % (s.real, s.imag, n, len(sols), np.round(ue, 4)))
    sys.stdout.flush()
# decomposition check at eps=0 and eps=0.02
for eps in (0.0, 0.02):
    ch = Chain(g, 1, kappa, lam, eps=eps)
    rows = [r for r in out[eps] if r[2] != 0]
    print("decomposition check eps=%+.2f, contributing:" % eps, [(round(r[1].real, 4), round(r[1].imag, 4), r[2]) for r in rows])
    for hbar in (2.0, 1.0, 0.5, 0.25):
        Zex = Z_grid(ch, hbar, h=0.02); Zex2 = Z_grid(ch, hbar, h=0.03)
        tot_ex = 0; tot1 = 0; tot2 = 0; parts = []
        for (U, s, n, mu, V, sols) in rows:
            Z1 = ch.one_loop(U, hbar, V, mu); c1 = ch.two_loop_c1(U)
            Z2 = Z1*(1 + hbar*c1)
            ratio = max(mu)/min(mu)
            if ratio < 1.6:
                Zs, dIm = thimble_integral(ch, U, hbar, V, mu, nrho=32, nphi=64, dt=0.005)
                exact = True
            else:
                Zs, exact = Z2, False
            tot_ex += n*Zs; tot1 += n*Z1; tot2 += n*Z2
            parts.append((s, n, Zs, exact))
        err = abs(tot_ex - Zex)/abs(Zex)
        print("  hbar=%.2f Z_exact=%.10e%+.2ei (grid h=.02 vs .03 rel diff %.1e) | sum n Z_sigma (exact thimbles where available, else 2-loop) rel.err=%.2e | all 1-loop %.2e | all 2-loop %.2e" % (
            hbar, Zex.real, Zex.imag, abs(Zex - Zex2)/abs(Zex), err, abs(tot1 - Zex)/abs(Zex), abs(tot2 - Zex)/abs(Zex)))
        # residual test for non-exact (complex pair): Z_exact - exact parts vs 2-loop of others
        ex_part = sum(n*Zs for (s, n, Zs, e) in parts if e)
        approx = [(s, n, Zs) for (s, n, Zs, e) in parts if not e]
        if approx:
            resid = Zex - ex_part
            pred = sum(n*Zs for (s, n, Zs) in approx)
            pred1 = sum(n*ch.one_loop(U, hbar, V, mu) for (U, s, n, mu, V, sols) in rows if max(mu)/min(mu) >= 1.6)
            print("      residual Z_exact - exact thimbles = %.6e%+.2ei ; 2-loop prediction from n_sigma of the others = %.6e%+.2ei (rel.err %.2e; 1-loop rel.err %.2e); residual/Z_exact=%.2e" % (
                resid.real, resid.imag, pred.real, pred.imag, abs(resid - pred)/abs(resid), abs(resid - pred1)/abs(resid), abs(resid)/abs(Zex)))
        sys.stdout.flush()
