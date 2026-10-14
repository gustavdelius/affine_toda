"""All foldings: one-loop soliton masses by the counting rule, MacKay-Watts-style single count,
comparison with one-loop particle masses (Lie-dual real coupling, and same theory continued)."""
import numpy as np, json, sys
from fold import folded
from particles import folded_particles
from massformula import psi
pi = np.pi

def drev(r):
    p = [0]*(r+1); p[0], p[1], p[r-1], p[r] = r-1, r, 0, 1
    for j in range(2, r-1): p[j] = r - j
    return p
def drho(r):
    p = [0]*(r+1); p[0], p[r], p[1], p[r-1] = r, 1, r-1, 0
    for j in range(2, r-1): p[j] = r - j
    return p
def ds1(r):
    p = list(range(r+1)); p[0], p[1] = 1, 0; p[r-1], p[r] = r, r-1
    return p
def dspin(r):
    p = list(range(r+1)); p[r-1], p[r] = r, r-1; return p
def arefl(r, shift=0):          # a_r^(1): j -> shift - j mod r+1
    N = r + 1; return [(shift - j) % N for j in range(N)]

E6_Z2 = [0, 6, 2, 5, 4, 3, 1]
E6_Z3 = [1, 6, 3, 5, 4, 2, 0]
E7_Z2 = [7, 6, 2, 5, 4, 3, 1, 0]
D4_Z3 = [0, 3, 2, 4, 1]

def soliton_data(F):
    h = F['h']; out = []
    for s in F['sols']:
        extra = 0.0; mults = []
        for co, mb, roots in s['chans']:
            for x, o in roots:
                if o > 1 and x < h/2 - 1e-9:
                    w = np.cos(pi*x/h); om = mb*np.sqrt(1 - w*w)
                    extra += (o - 1)*0.5*om
                    mults.append(dict(channel=co, m=mb, x=x, order=o, omega=om, lam=om**2))
        out.append(dict(cons=s['cons'], Mcl=s['Mcl'], r=s['dM']/s['Mcl'], rMW=(s['dM'] - extra)/s['Mcl'],
                        dM=s['dM'], mults=mults))
    return out

def thresholds(F):
    return sorted(set(round(mb**2, 10) for co, mb, roots in F['sols'][0]['chans']))

def report(name, F, part_dual, part_same, f=1.0):
    """part_dual / part_same: (masses, dlog m/dbeta_p^2) of Lie-dual and same theory (beta_p = parent coupling)."""
    S = soliton_data(F)
    Ms = np.array([s['Mcl'] for s in S])
    order = np.argsort(Ms); S = [S[i] for i in order]; Ms = Ms[order]
    res = dict(name=name, solitons=[], h_parent=F['h'])
    print(f"\n##### {name}  (beta_parent^2 = {f} beta^2)")
    thr = thresholds(F)
    print("  folded channel masses^2 (thresholds):", np.round(thr, 6))
    for s in S:
        print(f"  soliton {s['cons']}: M_cl = {s['Mcl']:.6f}/beta^2, ratio to lightest {s['Mcl']/Ms[0]:.6f}, dM = {s['dM']:.8f},"
              f" dM/M_cl = {s['r']:.8f} beta^2 (single count: {s['rMW']:.8f})")
        for mu in s['mults']:
            at_thr = [t for t in thr if abs(t - mu['lam']) < 1e-8]
            below = mu['lam'] < thr[0] - 1e-9
            print(f"      zero of order {mu['order']} in channel {mu['channel']} (m={mu['m']:.5f}) at x={mu['x']:g}: lambda={mu['lam']:.6f}, "
                  + ("AT THRESHOLD" if at_thr else ("isolated (below all thresholds)" if below else "embedded")))
        res['solitons'].append(dict(cons=s['cons'], Mcl=s['Mcl'], ratio=s['Mcl']/Ms[0], r=s['r'], rMW=s['rMW'],
                                    mults=[dict(order=m['order'], lam=m['lam'], x=m['x'], m=m['m']) for m in s['mults']]))
    # match classical ratios with particle masses
    for label, P, sign in (("Lie-dual real coupling, beta_r^2 = +beta^2", part_dual, +1), ("same theory continued, beta_r^2 = -beta^2", part_same, -1)):
        if P is None: continue
        m, dl = P
        dl = dl*f        # per folded beta^2
        mi = np.argsort(m); m = m[mi]; dl = dl[mi]
        # assign: soliton i <-> particle with matching ratio
        rat_s = Ms/Ms[0]; rat_p = m/m[0]
        ok = len(rat_s) == len(rat_p) and np.allclose(rat_s, rat_p, atol=1e-6)
        if not ok:
            print(f"  [{label}] classical ratios do not match: solitons {np.round(rat_s,5)}, particles {np.round(rat_p,5)}")
            continue
        # within groups of degenerate classical mass, the soliton-particle pairing is fixed by labels;
        # we take the pairing that matches (and say so)
        import itertools
        groups = []; i = 0
        while i < len(m):
            g = [i]
            while g[-1] + 1 < len(m) and abs(m[g[-1]+1] - m[i]) < 1e-9: g.append(g[-1] + 1)
            groups.append(g); i = g[-1] + 1
        if any(len(g) > 1 for g in groups):
            best = None
            for perms in itertools.product(*[list(itertools.permutations(g)) for g in groups]):
                idx = [k for p in perms for k in p]
                dd = dl[idx]
                dev = max(abs((S[i]['r'] - S[0]['r']) - sign*(dd[i] - dd[0])) for i in range(len(S)))
                if best is None or dev < best[0]: best = (dev, idx)
            dl = dl[best[1]]
            print(f"  [{label}] degenerate classical masses: pairing within degenerate groups chosen by match: {best[1]}")
        dev_rule = [(S[i]['r'] - S[0]['r']) - sign*(dl[i] - dl[0]) for i in range(len(S))]
        dev_MW = [(S[i]['rMW'] - S[0]['rMW']) - sign*(dl[i] - dl[0]) for i in range(len(S))]
        print(f"  [{label}] particle dlog(m_a/m_1) per beta^2: {np.round(sign*(dl - dl[0]), 8)}")
        print(f"      soliton  dlog(M_a/M_1) counting rule : {np.round([S[i]['r'] - S[0]['r'] for i in range(len(S))], 8)}  -> max dev {max(abs(np.array(dev_rule))):.1e}")
        print(f"      soliton  dlog(M_a/M_1) single count  : {np.round([S[i]['rMW'] - S[0]['rMW'] for i in range(len(S))], 8)}  -> max dev {max(abs(np.array(dev_MW))):.1e}")
        res[label] = dict(part=list(sign*(dl - dl[0])), rule=[S[i]['r'] - S[0]['r'] for i in range(len(S))],
                          single=[S[i]['rMW'] - S[0]['rMW'] for i in range(len(S))],
                          dev_rule=float(max(abs(np.array(dev_rule)))), dev_single=float(max(abs(np.array(dev_MW)))),
                          universal_rule=[S[i]['r'] - sign*dl[i] for i in range(len(S))])
        print(f"      M_a^qu/m_a^qu - classical, per beta^2 (should be a-independent): {np.round([S[i]['r'] - sign*dl[i] for i in range(len(S))], 8)}")
    return res

if __name__ == "__main__":
    which = sys.argv[1:] or ['all']
    allres = {}
    def want(k): return 'all' in which or k in which
    if want('cn'):
        for n in [2, 3, 4, 5]:
            F = folded('a', 2*n-1, arefl(2*n-1), f"c_{n}", verbose=False)
            allres[f"c_{n}^(1)"] = report(f"c_{n}^(1) from a_{2*n-1}^(1)", F,
                   folded_particles('d', n+1, dspin(n+1)) if n >= 3 else None, folded_particles('a', 2*n-1, arefl(2*n-1)))
    if want('bn'):
        for n in [3, 4, 5, 6, 7]:
            F = folded('d', n+1, dspin(n+1), f"b_{n}", verbose=False)
            allres[f"b_{n}^(1)"] = report(f"b_{n}^(1) from d_{n+1}^(1)", F,
                   folded_particles('a', 2*n-1, arefl(2*n-1)), folded_particles('d', n+1, dspin(n+1)))
    if want('g2'):
        F = folded('d', 4, D4_Z3, "g2", verbose=False)
        P = folded_particles('d', 4, D4_Z3)
        allres["g_2^(1)"] = report("g_2^(1) from d_4^(1)", F, P, P)
    if want('f4'):
        F = folded('e', 6, E6_Z2, "f4", verbose=False)
        P = folded_particles('e', 6, E6_Z2)
        allres["f_4^(1)"] = report("f_4^(1) from e_6^(1)", F, P, P)
    if want('tw'):
        for n in [2, 3, 4]:
            F = folded('d', 2*n, drev(2*n), "a2n-1", verbose=False)
            P = folded_particles('d', 2*n, drev(2*n))
            allres[f"a_{2*n-1}^(2)"] = report(f"a_{2*n-1}^(2) from d_{2*n}^(1)", F, None, P)
        for n in [3, 4, 5, 6]:
            F = folded('d', n+2, ds1(n+2), "dn+1", verbose=False)
            P = folded_particles('d', n+2, ds1(n+2))
            allres[f"d_{n+1}^(2)"] = report(f"d_{n+1}^(2) from d_{n+2}^(1)", F, None, P)
        for n in [1, 2, 3]:
            F = folded('d', 2*n+2, drho(2*n+2), "a2n", verbose=False)
            P = folded_particles('d', 2*n+2, drho(2*n+2))
            allres[f"a_{2*n}^(2)"] = report(f"a_{2*n}^(2) from d_{2*n+2}^(1) (Z4)", F, None, P)
        F = folded('e', 6, E6_Z3, "d43", verbose=False); P = folded_particles('e', 6, E6_Z3)
        allres["d_4^(3)"] = report("d_4^(3) from e_6^(1) (Z3)", F, None, P)
        F = folded('e', 7, E7_Z2, "e62", verbose=False); P = folded_particles('e', 7, E7_Z2)
        allres["e_6^(2)"] = report("e_6^(2) from e_7^(1)", F, None, P)
    if want('indirect'):
        for n in [1, 2, 3, 4]:
            F = folded('a', 2*n, arefl(2*n), "a2n ind", verbose=False)
            P = folded_particles('a', 2*n, arefl(2*n))
            allres[f"a_{2*n}^(2) ind"] = report(f"a_{2*n}^(2) from a_{2*n}^(1) (indirect)", F, None, P)
        for n in [2, 3, 4]:
            F = folded('a', 2*n+1, arefl(2*n+1, shift=-1), "dn+1 ind", beta2_factor=2.0, verbose=False)
            P = folded_particles('a', 2*n+1, arefl(2*n+1, shift=-1))
            allres[f"d_{n+1}^(2) ind"] = report(f"d_{n+1}^(2) from a_{2*n+1}^(1) (indirect)", F, None, P, f=2.0)
    json.dump(allres, open("results_" + "_".join(which) + ".json", "w"), indent=1, default=float)
