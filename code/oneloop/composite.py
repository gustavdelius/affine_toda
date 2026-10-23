"""Hirota form (H) for folded composite solitons: build the sigma-invariant coincident multisoliton by the
tau-recursion, then the Hirota fluctuation in each folded channel, and compare its transmission factor with
the product of parent X_ab (the product rule used in fold.py)."""
import numpy as np, sys
from toda import soliton, fluct, degree, perm_matrix
from fold import load, sigma_action, orbits
pi = np.pi

def composite_check(alg, r, perm, name, zs=(0.31+0.22j, -0.47+0.61j, 0.83-0.35j)):
    D = load(alg, r); C, n, h = D['C'], D['n'], D['h']
    masses, vecs, table = D['masses'], D['vecs'], D['table']
    P = perm_matrix(perm)
    act = sigma_action(P, vecs); orbs = [o for o, chi in orbits(act) if abs(chi - 1) < 1e-8]
    print(f"\n== {name} ({alg}_{r}^(1)): invariant orbits {orbs}")
    worst = 0
    for so in orbs:
        if len(so) == 1: continue
        mu = masses[so[0]]
        G = [np.eye(len(n))]; g = P.copy()
        while not np.allclose(g, np.eye(len(n))): G.append(g); g = g @ P
        k0 = np.argmax(np.abs(vecs[so[0]]))
        delta = sum(gm @ vecs[so[0]] for gm in G); delta = delta/delta[k0] if abs(delta[k0]) > 1e-9 else delta
        Kmax = int(len(so)*max(n)) + 3
        c = soliton(C, n, mu, delta, Kmax=2*Kmax)
        deg = degree(c)
        expect = [int(len(so)*nj) for nj in n]
        print(f"  composite {so}: tau degrees {deg} (expected |orbit| n_j = {expect})")
        for co in orbs:
            u = sum(gm @ vecs[co[0]] for gm in G)
            mb = masses[co[0]]
            # product rule prediction
            from collections import Counter
            cnt = Counter()
            for cc in so:
                for x, o in table[(cc, co[0])]: cnt[x] += o
            for z in zs:
                p = fluct(C, n, c, mu, mb, u, z*mb)
                Dp = degree(p, tol=1e-8)
                ratio = np.array([p[j, deg[j]]/c[j, deg[j]] for j in range(len(n))])
                X = np.vdot(u, ratio)/np.vdot(u, u)
                res = np.linalg.norm(ratio - X*u)/np.linalg.norm(ratio)
                Xp = np.prod([(z - np.cos(pi*x/h))**o for x, o in cnt.items()])
                worst = max(worst, abs(X/Xp - 1))
                term = all(Dp[j] <= deg[j] for j in range(len(n)))
            print(f"     channel {co}: fluctuation degree <= tau degree: {term};  |X_Hirota/X_product - 1| <= {abs(X/Xp-1):.1e}; diag resid {res:.1e}")
    print(f"  worst deviation from product rule: {worst:.1e}")
    return worst

if __name__ == "__main__":
    from analysis import dspin, E6_Z2, D4_Z3, arefl
    composite_check('d', 4, dspin(4), "b_3")
    composite_check('d', 5, dspin(5), "b_4")
    composite_check('d', 6, dspin(6), "b_5")
    composite_check('d', 4, D4_Z3, "g_2")
    composite_check('e', 6, E6_Z2, "f_4")
    composite_check('a', 7, arefl(7), "c_4")
