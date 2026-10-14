"""Folded solitons: transmission factors as products of parent X_ab, one-loop masses by the counting rule.

Parent data from xtable.run (Hirota-computed X_ab roots on the grid z = cos(pi x/h)).
"""
import numpy as np, pickle, os, sys
from collections import Counter, defaultdict
from massformula import psi
from toda import perm_matrix
pi = np.pi
HERE = os.path.dirname(os.path.abspath(__file__))

def load(alg, r):
    fn = os.path.join(HERE, f"X_{alg}{r}.pkl")
    if os.path.exists(fn):
        return pickle.load(open(fn, "rb"))
    from xtable import run
    out = run(alg, r, verbose=False)
    pickle.dump(out, open(fn, "wb"))
    return out

def sigma_action(P, vecs, tol=1e-8):
    """for each species a: (image species a', phase chi) with P v_a = chi v_a'."""
    act = []
    for a, v in enumerate(vecs):
        w = P @ v
        found = None
        for b, u in enumerate(vecs):
            chi = np.vdot(u, w)/np.vdot(u, u)
            if np.linalg.norm(w - chi*u) < tol*np.linalg.norm(w):
                found = (b, chi); break
        assert found is not None, "species not mapped to a species"
        act.append(found)
    return act

def orbits(act):
    seen, orbs = set(), []
    for a in range(len(act)):
        if a in seen: continue
        orb = [a]; b = act[a][0]
        while b != a:
            orb.append(b); b = act[b][0]
        # total phase around the orbit (for the invariant combination to exist it must be 1)
        chi = 1.0+0j
        for c in orb: chi *= act[c][1]
        seen.update(orb); orbs.append((orb, chi))
    return orbs

def mult(roots_list):
    cnt = Counter()
    for roots in roots_list:
        for x, o in roots: cnt[x] += o
    return sorted((x, o) for x, o in cnt.items() if o != 0)

def folded(alg, r, perm, name, beta2_factor=1.0, verbose=True):
    """beta2_factor: beta_parent^2 = beta2_factor * beta_folded^2 (1 for direct folds with a fixed long node)."""
    D = load(alg, r)
    h, masses, vecs, table = D['h'], D['masses'], D['vecs'], D['table']
    P = perm_matrix(perm)
    act = sigma_action(P, vecs)
    orbs = orbits(act)
    inv = [(o, chi) for o, chi in orbs if abs(chi - 1) < 1e-8]
    sol_orbs = [o for o, chi in inv]          # folded solitons = invariant composites of the orbit
    ch_orbs = [o for o, chi in inv]           # folded channels = invariant combinations
    excluded = [(o, chi) for o, chi in orbs if abs(chi - 1) >= 1e-8]
    res = []
    for so in sol_orbs:
        chans = []
        for co in ch_orbs:
            # transmission of the folded channel: product over constituents, same for each orbit member
            prods = [mult([table[(c, b)] for c in so]) for b in co]
            assert all(p == prods[0] for p in prods), "orbit members disagree"
            chans.append((co, masses[co[0]], prods[0]))
        Mcl = 2*h*sum(masses[c] for c in so)*beta2_factor   # coefficient of 1/beta_folded^2 ... see below
        # classical mass = (2h/beta_p^2) sum m_c = (2h/(f beta^2)) sum m_c ; we store coefficient of 1/beta^2
        Mcl = 2*h*sum(masses[c] for c in so)/beta2_factor
        dM = 0.5*sum(mb*sum(o*psi(np.cos(pi*x/h)) for x, o in roots) for co, mb, roots in chans).real
        # MacKay-Watts style count: each bound-state zero counted once (order capped at 1), poles ignored
        res.append(dict(cons=so, Mcl=Mcl, dM=dM, chans=chans))
    if verbose:
        print(f"\n=== {name}: folded from {alg}_{r}^(1), h_parent = {h}, beta_p^2 = {beta2_factor} beta^2")
        print("  species masses:", np.round(masses, 5))
        print("  sigma action (image, phase):", [(b, np.round(chi, 3)) for b, chi in act])
        print("  invariant orbits:", [o for o in sol_orbs], " excluded:", [(o, np.round(chi, 3)) for o, chi in excluded])
    return dict(name=name, h=h, masses=masses, sols=res, ch_orbs=ch_orbs, act=act, D=D, beta2_factor=beta2_factor)

def describe(F, prec=6):
    h = F['h']
    for s in F['sols']:
        print(f"  soliton {s['cons']}: M_cl = {s['Mcl']:.6f}/beta^2,  Delta M = {s['dM']:.8f},  Delta M/M_cl = {s['dM']/s['Mcl']:.8f} beta^2")
        for co, mb, roots in s['chans']:
            bs = [(x, o) for x, o in roots if x < h/2 - 1e-9]       # Re z > 0: bound-state half plane
            thr = [(x, o) for x, o in roots if abs(x - h/2) < 1e-9]
            txt = " ".join(f"{x:g}{'+' if o > 0 else '-'}{abs(o) if abs(o) > 1 else ''}" for x, o in roots)
            print(f"     channel {co} (m={mb:.5f}): roots x[ord] = {txt}   bound half-plane: {bs}  at threshold z=0: {thr}")
