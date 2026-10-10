"""Lie duality for twisted algebras (10-solitons.qmd:24): fold by a centre element g of the simply-laced parent.
Folded particles = M^2 eigenvectors in the g-invariant subspace; folded single solitons = parent single solitons of
character chi_g(a)=1 (non-resonance => g maps soliton (a,xi) to (a, xi + i arg chi), ch.18 step 4), each with
M = 2 h_parent m_a/beta^2.  Compare soliton-mass ratios with the folded particle ratios and with the dual's particles."""
import numpy as np
from lie import Algebra
def report(n, ok): print(('PASS ' if ok else 'FAIL ') + n)
def fold_masses(Al, g):
    Ms = np.array([Al.alpha[j] - Al.alpha[g[j]] for j in range(Al.r + 1)])
    u, s, vt = np.linalg.svd(Ms); V = vt[int((s > 1e-9).sum()):].T
    return np.sort(np.sqrt(np.abs(np.linalg.eigvalsh(V.T @ Al.M0 @ V))))
dual_particles = {   # particle masses of the dual untwisted theories (algebra_data.py closed forms)
    'g2': np.array([np.sqrt(2), np.sqrt(6)]),
    'c3': 2*np.sin(np.pi*np.arange(1, 4)/6), 'b3': np.sort(np.r_[2*np.sqrt(2)*np.sin(np.pi*np.arange(1, 3)/6), np.sqrt(2)]),
    'c2': 2*np.sin(np.pi*np.arange(1, 3)/4),
}
cases = [('d4^(3) <- e6, Z3', ('e', 6), 'g2'), ('d3^(2) <- d4, (01)(34)', ('d', 4), 'c2'),
         ('d4^(2) <- d5, (01)(45)', ('d', 5), 'c3'), ('a5^(2) <- d6, flip', ('d', 6), 'b3'), ('e6^(2) <- e7, Z2', ('e', 7), None)]
ok_all = True
for name, (typ, r), dual in cases:
    Al = Algebra(typ, r)
    for g in Al.centre:
        if g[0] == 0: continue
        if typ == 'd' and r in (4, 5) and not (g[0] == 1): continue          # the (0 1)(..) element
        if typ == 'd' and r == 6 and g[0] == 1: continue                        # the flip moving 0 to a far fork node
        inv = sorted(a for a in Al.mass if abs(Al.chi(g, a) - 1) < 1e-9)
        sol = np.sort([Al.mass[a] for a in inv]); part = fold_masses(Al, g)
        same = len(sol) == len(part) and np.allclose(sol, part)
        ok_all &= same
        line = f'   {name} (g={g}): invariant species {inv}; soliton ratios {np.round(sol/sol[0], 5)}; folded particle ratios {np.round(part/part[0], 5)}'
        if dual is not None:
            dp = np.sort(dual_particles[dual]); line += f'; dual {dual}^(1) particle ratios {np.round(dp/dp[0], 5)}'
        print(line)
        break
report('twisted foldings: single-soliton mass ratios = particle mass ratios of the same (twisted) theory', ok_all)
