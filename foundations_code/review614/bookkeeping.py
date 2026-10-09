"""Test of Step 3 of the proof of Theorem 6.14 on random multi-scale clusters.

Rescaled units (ell = 1), Dirichlet wall at y = 0, lattice cutoff b.
For each random 1-connected cluster Q we compute
  F  = actual product of m-th-neighbour factors and intra-unit pair factors of (6.5)^p,
       for candidate B (units = maximal proper tight subsets of Q) and, if |Q|<=m,
       candidate A (Q inside one unit, pair scales min(1, sqrt(d_k d_l)) up to factor 6);
  M  = prod_j hat s_j^{-e_j} with the exponents placed by the rules of Step 3;
and check: facts (a),(b),(c); E_j > 0 at every node and E_j >= eta at active nodes;
log(F/M) bounded by a constant (independent of scale separation).
"""
import itertools, sys
import numpy as np
from collections import defaultdict

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)

# ---------- root data ----------
def a4():
    e = np.eye(5)
    al = [e[i]-e[i+1] for i in range(4)]
    al0 = -sum(al)
    return [al0]+al, [1, 1, 1, 1, 1]

def c2():
    e = np.eye(2)
    a1 = (e[0]-e[1])/np.sqrt(2)  # short, |.|^2 = 1
    a2 = np.sqrt(2)*e[1]         # long,  |.|^2 = 2
    a0 = -(2*a1+a2)
    return [a0, a1, a2], [1, 2, 1]

def g2():
    # G2 with long roots |.|^2=2: short simple a1, long simple a2
    # use standard realization in plane
    a1 = np.array([1.0, 0.0])*np.sqrt(2/3)
    a2 = np.array([-1.5, np.sqrt(3)/2])*np.sqrt(2/3)
    theta = 3*a1+2*a2
    return [-theta, a1, a2], [1, 3, 2]

ALG = {'a4': (a4, 0.79), 'c2': (c2, 0.76), 'g2': (g2, 0.76)}

def eta_of(roots, kappa, p, m, ext=None):
    best = np.inf
    idx = list(range(len(roots)))
    for size in range(2, m+1):
        for S in itertools.combinations_with_replacement(idx, size):
            v = [roots[i] for i in S]
            D = sum(x@x for x in v) - (sum(v)@sum(v))
            best = min(best, 2*(size-1) - p*kappa*max(D, 0))
        if ext is not None and size-1 >= 1:
            for S in itertools.combinations_with_replacement(idx, size-1):
                v = [roots[i] for i in S]+[ext]
                D = sum(x@x for x in v) - (sum(v)@sum(v))
                best = min(best, 2*(size-1) - p*kappa*max(D, 0))
    return best

# ---------- configurations ----------
def gen(n, center, scale):
    if n == 1:
        return [center]
    if rng.random() < 0.3:
        pts = []
        for _ in range(n):
            r = scale*np.sqrt(rng.random()); th = 2*np.pi*rng.random()
            pts.append(center+r*np.array([np.cos(th), np.sin(th)]))
        return pts
    k = 2 if (n == 2 or rng.random() < 0.6) else 3
    k = min(k, n)
    cuts = np.sort(rng.choice(np.arange(1, n), size=k-1, replace=False))
    sizes = np.diff(np.concatenate([[0], cuts, [n]]))
    pts = []
    for s in sizes:
        r = scale*np.sqrt(rng.random()); th = 2*np.pi*rng.random()
        c = center+r*np.array([np.cos(th), np.sin(th)])
        sub = scale*10**(-rng.uniform(0.0, 3.5))
        pts += gen(int(s), c, sub)
    return pts

def dendrogram(P):
    n = len(P)
    pairs = []
    for i in range(n):
        for j in range(i+1, n):
            pairs.append((np.linalg.norm(P[i]-P[j]), i, j))
    pairs.sort()
    comp = {i: i for i in range(n)}       # union-find parent
    top = {i: ('leaf', i) for i in range(n)}
    def find(x):
        while comp[x] != x:
            comp[x] = comp[comp[x]]; x = comp[x]
        return x
    nodes = []   # dicts: s, children, S, parent
    members = {i: {i} for i in range(n)}
    lca = {}
    for d, i, j in pairs:
        ri, rj = find(i), find(j)
        if ri == rj:
            continue
        jid = len(nodes)
        nodes.append({'s': d, 'children': [top[ri], top[rj]], 'S': members[ri] | members[rj], 'parent': None, 'edge': (i, j)})
        for c in (top[ri], top[rj]):
            if c[0] == 'node':
                nodes[c[1]]['parent'] = jid
        for k in members[ri]:
            for l in members[rj]:
                lca[(min(k, l), max(k, l))] = jid
        comp[ri] = rj
        members[rj] = members[ri] | members[rj]
        top[rj] = ('node', jid)
    return nodes, lca

def run_one(alg, p=1.0, Ktight=21.0, use_ext=False, caseA_scale=1.0):
    mk, kappa = ALG[alg]
    roots, marks = mk()
    gbar = 2*kappa
    m = int(np.floor(kappa/(1-kappa)))+1
    ext = None
    if use_ext:
        # random external weight satisfying (6.2) for clusters up to m-1 roots
        for _ in range(100):
            w = rng.normal(size=roots[0].shape)*rng.uniform(0.1, 1.5)
            ok = True
            idx = list(range(len(roots)))
            for size in range(1, m):
                for S in itertools.combinations_with_replacement(idx, size):
                    v = [roots[i] for i in S]
                    D = sum(x@x for x in v)+w@w-(w+sum(v))@(w+sum(v))
                    if not (2*size > kappa*D):
                        ok = False; break
                if not ok: break
            gam_ext = kappa*(w@w)
            if ok and (2-gbar)*m > gam_ext:
                ext = w; break
        if ext is None:
            return None
    eta = eta_of(roots, kappa, p, m, ext)
    n = int(rng.integers(2, 2*m+2))
    P = gen(n, np.zeros(2), 1.0)
    P = [np.array(x) for x in P]
    # rescale so that max MST edge <= 1
    nodes, lca = dendrogram(P)
    smax = max(nd['s'] for nd in nodes)
    fac = smax*rng.uniform(1.0, 3.0)
    P = [x/fac for x in P]
    # place near the wall y=0
    ymin = min(x[1] for x in P)
    dmin = 10**rng.uniform(-5, 0.7)
    P = [x+np.array([0, dmin-ymin]) for x in P]
    nodes, lca = dendrogram(P)
    b = 10**rng.uniform(-8, np.log10(1/6.0)-0.01)
    hat = lambda s: max(s, b)
    # species
    if use_ext:
        spec = [None]+list(rng.integers(0, len(roots), size=n-1))
        vec = [ext]+[roots[s] for s in spec[1:]]
    else:
        spec = list(rng.integers(0, len(roots), size=n))
        vec = [roots[s] for s in spec]
    gam = [kappa*(v@v) for v in vec]
    Dkl = lambda k, l: -2*(vec[k]@vec[l])
    dist = lambda k, l: np.linalg.norm(P[k]-P[l])
    dwall = [x[1] for x in P]
    # ---- fact (a)
    for (k, l), j in lca.items():
        r = dist(k, l); sj = nodes[j]['s']; Sj = len(nodes[j]['S'])
        assert sj <= r*(1+1e-12) and r <= (Sj-1)*sj*(1+1e-12), 'fact (a) fails'
    # ---- tight sets by brute force (proper subsets of Q), fact (b)
    allidx = list(range(n))
    def diam(S):
        return max([dist(k, l) for k, l in itertools.combinations(S, 2)], default=0.0)
    def dout(S):
        rest = [x for x in allidx if x not in S]
        return min(dist(k, l) for k in S for l in rest)
    node_of = {frozenset(nd['S']): j for j, nd in enumerate(nodes)}
    tight = []
    for size in range(1, n):
        for S in itertools.combinations(allidx, size):
            if size == 1 or Ktight*max(diam(S), 2*b) <= dout(S):
                tight.append(frozenset(S))
                if size >= 2:
                    assert frozenset(S) in node_of, 'fact (b): tight set not a node'
                    j = node_of[frozenset(S)]
                    par = nodes[j]['parent']
                    assert par is not None and abs(nodes[par]['s']-dout(S)) < 1e-12, 'fact (b): R_B not parent edge'
    # units (candidate B)
    small = [S for S in tight if len(S) <= m]
    unitsB = [S for S in small if not any(S < T for T in small)]
    assert sorted(sum([sorted(S) for S in unitsB], [])) == allidx
    # ---- m-th neighbour distances and i(k), fact (c)
    delta = []
    for k in allidx:
        ds = sorted(dist(k, l) for l in allidx if l != k)
        delta.append(ds[m-1] if len(ds) >= m else np.inf)
    ik = []
    for k in allidx:
        # ancestors of leaf k
        j = None
        for jj, nd in enumerate(nodes):
            if k in nd['S'] and len(nd['S']) >= m+1:
                if j is None or len(nodes[jj]['S']) < len(nodes[j]['S']):
                    j = jj
        ik.append(j)
        if j is not None:
            assert delta[k] >= nodes[j]['s']*(1-1e-12), 'fact (c) fails'
    def h(j): return len(nodes[j]['S'])-1
    results = []
    candidates = [('B', unitsB)]
    if n <= m:
        candidates.append(('A', [frozenset(allidx)]))
    for cand, units in candidates:
        # actual F
        logF = 0.0
        for k in allidx:
            logF += p*gam[k]*max(0.0, -np.log(hat(delta[k])))
        for B in units:
            if cand == 'B':
                RB = dout(B) if len(B) < n else np.inf
                rhoB = max(0.0, RB/3-b)
            else:
                rhoB = caseA_scale
            for k, l in itertools.combinations(sorted(B), 2):
                s = min(rhoB, 1.0, np.sqrt(dwall[k]*dwall[l]))
                logF += p*kappa*Dkl(k, l)*max(0.0, np.log(s/hat(dist(k, l))))
        # exponents
        e = defaultdict(float)
        for k in allidx:
            if ik[k] is not None:
                e[ik[k]] += p*gam[k]
        active_nodes = set()
        for B in units:
            if len(B) < 2: continue
            if cand == 'B':
                RB = dout(B); rhoB = max(0.0, RB/3-b)
            else:
                rhoB = 1.0
            inB = [j for j, nd in enumerate(nodes) if nd['S'] <= B]
            dj = {j: min(dwall[k] for k in nodes[j]['S']) for j in inB}
            sig = {j: min(rhoB, 1.0, dj[j]) for j in inB}
            act = {j for j in inB if hat(nodes[j]['s']) < sig[j]}
            # down-closedness check
            for j in act:
                for jj in inB:
                    if nodes[jj]['S'] <= nodes[j]['S']:
                        assert jj in act, 'active set not down-closed'
            active_nodes |= act
            for j in act:
                Dj = sum(Dkl(k, l) for (k, l), jj in lca.items() if jj == j)
                e[j] += p*kappa*Dj
            tops = [j for j in act if nodes[j]['parent'] not in act]
            for js in tops:
                S = nodes[js]['S']
                DS = sum(Dkl(k, l) for k, l in itertools.combinations(sorted(S), 2))
                if DS < 0:
                    e[js] -= p*kappa*DS
                else:
                    q = nodes[js]['parent']
                    if q is None:
                        pass  # cluster root: dropped
                    else:
                        e[q] -= p*kappa*DS
        logM = sum(-e[j]*np.log(hat(nodes[j]['s'])) for j in range(len(nodes)))
        # cumulative exponents
        Cj = {}
        for j, nd in enumerate(nodes):
            Cj[j] = sum(e[jj] for jj, nd2 in enumerate(nodes) if nd2['S'] <= nd['S'])
        Ej = {j: 2*h(j)-Cj[j] for j in range(len(nodes))}
        minE = min(Ej.values())
        minEact = min([Ej[j] for j in active_nodes], default=np.inf)
        spread = np.log(max(nd['s'] for nd in nodes)/max(min(nd['s'] for nd in nodes), b))
        results.append(dict(cand=cand, n=n, logFM=logF-logM, minE=minE, minEact=minEact,
                            eta=eta, spread=spread, nact=len(active_nodes)))
    return results

if __name__ == '__main__':
    out = []
    bad = 0
    for alg in ['a4', 'c2', 'g2']:
        for use_ext in [False, True]:
            for caseA_scale in [1.0, 1/6]:
                for trial in range(1500):
                    r = run_one(alg, use_ext=use_ext, caseA_scale=caseA_scale)
                    if r is None: continue
                    for x in r:
                        x['alg'] = alg; x['ext'] = use_ext
                        out.append(x)
                        if x['minE'] <= 0 or x['minEact'] < x['eta']-1e-9:
                            bad += 1
                            print('BAD', x)
    logFM = np.array([x['logFM'] for x in out]); n = np.array([x['n'] for x in out])
    spread = np.array([x['spread'] for x in out]); nact = np.array([x['nact'] for x in out])
    print('samples', len(out), 'bad E', bad)
    print('max log(F/M)', logFM.max(), ' max log(F/M)/n', (logFM/n).max())
    print('samples with active nodes', (nact > 0).sum())
    for lo, hi in [(0, 3), (3, 6), (6, 10), (10, 20), (20, 60)]:
        sel = (spread >= lo) & (spread < hi)
        if sel.any():
            print(f'scale spread in [{lo},{hi}): count {sel.sum()}, max log(F/M) {logFM[sel].max():.3f}')
