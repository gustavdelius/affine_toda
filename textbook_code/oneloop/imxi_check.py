"""The one-loop mass does not depend on Im xi (sec-which-saddles-contribute, Chapter 17).

Run from this directory:  python3 imxi_check.py     (about 15 min)

For the g_2 triple composite and the b_3 spinor composite (the folded solitons of lattice.py), evaluate the
one-loop energy (8.1) directly from the eigenvalues of the Fourier-grid fluctuation operator at several values
of Im xi: the most regular one, and values 1.0 and 0.3 away from the nearest Im xi at which some tau_j vanishes
on the real line. Each must agree with the closed-form Delta M of prp-one-loop-mass.

Near a singular Im xi the soliton is sharply peaked, so the grid must be fine: N = 256 points on L = 48 is off by
4% at distance 1.0 and wrong in sign at 0.3, while N = 1024 agrees to 1e-3 at both.
"""
import numpy as np
from lattice import build, spectrum
from analysis import dspin, D4_Z3
from fold import folded

L, N, TOL = 48.0, 1024, 2e-3
cases = [("g_2 triple composite", 'd', 4, D4_Z3, 0), ("b_3 spinor composite", 'd', 4, dspin(4), 1)]

for name, alg, r, perm, a0 in cases:
    B = build(alg, r, perm, a0)
    F = folded(alg, r, perm, name, verbose=False)
    closed = [s for s in F['sols'] if a0 in s['cons']][0]['dM']
    sing = []                                      # Im xi at which some tau_j has a zero on the real line
    for j in range(len(B['n'])):
        cj = B['c'][j][:B['deg'][j]+1]
        if len(cj) > 1: sing += list(np.angle(np.roots(cj[::-1])))
    a = min(sing)
    print(f"{name}: closed-form Delta M = {closed:.6f}; singular Im xi at {np.round(sorted(set(np.round(sing, 8))), 4)}")
    for d in [None, 1.0, 0.3]:
        xi = None if d is None else 1j*(a + d)
        lam, lam0, dM, xiu, mn = spectrum(B, L=L, Npts=N, xi=xi)
        ok = abs(dM - closed) < TOL
        where = "most regular" if d is None else f"{d} from singular"
        print(f"{'PASS' if ok else 'FAIL'} {name}, Im xi = {xiu.imag:+.4f} ({where}): "
              f"Delta M = {dM.real:.6f} {dM.imag:+.1e}i, |diff| = {abs(dM - closed):.1e}", flush=True)
