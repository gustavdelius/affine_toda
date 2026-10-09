"""a_3^(1): charges reached by static multi-solitons of distinct species, and by 1+1 'charged breathers'
(two species-1 solitons with rapidities theta = +-iu) at t = 0. Energies: static sum of M_a; breather 2 M_1 cos u."""
import numpy as np, collections, itertools, sys
from lie import An
from charges import charge_label, rep_name
g = An(3)
rng = np.random.default_rng(2)
def scan(species, thetas, nsamp, re_range=4.0):
    labs = collections.Counter()
    for k in range(nsamp):
        xis = [rng.uniform(-re_range, re_range) + 1j*rng.uniform(0, 2*np.pi) for _ in species]
        lab = charge_label(g, species, thetas, xis, X=20.0, N=8001)
        labs[rep_name(g, lab)] += 1
    return labs
M = {a: 16*np.sin(np.pi*a/4) for a in (1,2,3)}  # soliton energies in units of sqrt(kappa) (u-variables)
print("soliton actions (units sqrt(kappa)/hbar):", M)
for sp in ([2], [1,3], [1,2], [2,3], [1,2,3]):
    labs = scan(sp, [0.0]*len(sp), 600)
    print("static species", sp, "E =", round(sum(M[a] for a in sp),3), dict(labs)); sys.stdout.flush()
for u in (np.pi/4 - 0.15, np.pi/4 - 0.05, np.pi/4 + 0.05, np.pi/4+0.15):
    labs = scan([1,1], [1j*u, -1j*u], 600)
    print("1+1 breather u=%.3f (fusing angle pi/4) E=2M1 cos u=%.3f"%(u, 2*M[1]*np.cos(u)), dict(labs)); sys.stdout.flush()
