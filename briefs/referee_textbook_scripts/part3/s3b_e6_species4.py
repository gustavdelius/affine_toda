"""Rerun e6 species 4 (heaviest self-conjugate) at higher precision: the dps-30 run did not truncate."""
import numpy as np, mpmath as mp
exec(open('s3_soliton_mass.py').read().split('for typ, r in')[0].replace("delta_mp(Al, a, dps=40)", "delta_mp(Al, a, dps=80, iters=10)")
     .replace("mp.mp.dps = 30\n    mu", "mp.mp.dps = 60\n    mu").replace("MpCtx(30), Kmax=14, tol=mp.mpf(10)**-22", "MpCtx(60), Kmax=14, tol=mp.mpf(10)**-35"))
Al = Algebra('e', 6)
E, pred, formula, degs, n, mn, xi = energy(Al, 4)
print(f'   e6 species 4: deg tau_j={degs} (n_j={n}), E={mp.nstr(E, 16)}, 2h m_a/beta^2={mp.nstr(pred, 16)}')
report('e6 species 4: E = 2 h m_a/beta^2 and deg tau_j = n_j', abs(E - pred) < 1e-15*abs(pred) and degs == n)
