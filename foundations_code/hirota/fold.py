"""Folded (non-simply-laced and twisted) theories as sigma-invariant sectors of simply-laced parents.

Folded soliton of species a: the coincident static sigma-orbit of parent soliton a (first-order data
sum_{a' in orbit} delta^{a'}); for central sigma the orbit is {a} and only species with character 1 occur.
Folded channel b: sigma-invariant parent first-order data sum_k delta^b o sigma^k.

For each folded soliton and channel: exact fluctuation solution (Hirota), check of (H) in the invariant
sector, X_b in closed form (fit to cos(pi q/h) factors), comparison with prod_{a' in orbit} X_{a'b}.
"""
import sys, json
import numpy as np, mpmath as mp
from lie import Algebra, delta_mp
from hirota import MpCtx, soliton_deg, fluct_deg
from compare_dorey import X_dorey

DPS = 30

FOLDS = {
    # name: (parent, sigma as dict on affine nodes (others fixed), kind)
    'c2(1)<a3': ('a3', 'refl_a'), 'c3(1)<a5': ('a5', 'refl_a'), 'c4(1)<a7': ('a7', 'refl_a'),
    'a2(2)<a2': ('a2', 'refl_a'), 'a4(2)<a4': ('a4', 'refl_a'), 'a6(2)<a6': ('a6', 'refl_a'),
    'b3(1)<d4': ('d4', {3: 4, 4: 3}), 'b4(1)<d5': ('d5', {4: 5, 5: 4}), 'b5(1)<d6': ('d6', {5: 6, 6: 5}),
    'b6(1)<d7': ('d7', {6: 7, 7: 6}),
    'g2(1)<d4': ('d4', {1: 3, 3: 4, 4: 1}),
    'f4(1)<e6': ('e6', {1: 6, 6: 1, 3: 5, 5: 3}),
    # twisted: central sigma
    'd3(2)<d4': ('d4', {0: 1, 1: 0, 3: 4, 4: 3}), 'd4(2)<d5': ('d5', {0: 1, 1: 0, 4: 5, 5: 4}),
    'd5(2)<d6': ('d6', {0: 1, 1: 0, 5: 6, 6: 5}),
    'a3(2)<d4': ('d4', 'refl_d'), 'a5(2)<d6': ('d6', 'refl_d'),
    'd4(3)<e6': ('e6', 'z3_e6'),
    'e6(2)<e7': ('e7', 'z2_e7'),
}


def make_sigma(A, spec):
    r = A.r
    if spec == 'refl_a':
        h = r + 1
        return tuple((-j) % h for j in range(r + 1))
    if spec == 'refl_d':
        # d_r reflection with fixed middle node, centre element for r even: 0<->r, 1<->r-1, j<->r-j
        perm = list(range(r + 1))
        perm[0], perm[r], perm[1], perm[r - 1] = r, 0, r - 1, 1
        for j in range(2, r - 1):
            perm[j] = r - j
        g = tuple(perm)
        assert g in A.autos and g in A.centre, g
        return g
    if spec in ('z3_e6', 'z2_e7'):
        g = [c for c in A.centre if c[0] != 0][0]
        return tuple(g)
    perm = list(range(r + 1))
    for i, j in spec.items():
        perm[i] = j
    g = tuple(perm)
    assert g in A.autos, (spec, A.autos)
    return g


def compose_pow(g, k):
    p = list(range(len(g)))
    for _ in range(k):
        p = [g[i] for i in p]
    return p


def order(g):
    k = 1
    while compose_pow(g, k) != list(range(len(g))):
        k += 1
    return k


class Fold:
    def __init__(self, name):
        self.name = name
        pname, spec = FOLDS[name]
        self.A = A = Algebra(pname[0], int(pname[1:]))
        self.sigma = make_sigma(A, spec)
        self.central = self.sigma in A.centre
        self.ord = order(self.sigma)
        self.ctx = MpCtx(DPS)
        self.d, self.lam = {}, {}
        for b in A.mass:
            self.d[b], self.lam[b] = delta_mp(A, b, dps=DPS)[:2]
        # invariant data for each species
        self.inv = {}
        for b in A.mass:
            v = [mp.mpc(0)] * (A.r + 1)
            for k in range(self.ord):
                p = compose_pow(self.sigma, k)
                for j in range(A.r + 1):
                    v[j] += self.d[b][p[j]]
            if max(abs(x) for x in v) > 1e-10:
                self.inv[b] = [x / v[0] for x in v] if abs(v[0]) > 1e-10 else v
        # orbits of species: species b, b' in same orbit iff inv vectors proportional
        self.orbits = []
        for b in sorted(self.inv):
            for O in self.orbits:
                c = O[0]
                if abs(self.lam[c] - self.lam[b]) < 1e-12:
                    # proportional?
                    ratio = None
                    ok = True
                    for j in range(A.r + 1):
                        x, y = self.inv[b][j], self.inv[c][j]
                        if abs(y) > 1e-12:
                            q = x / y
                            if ratio is None:
                                ratio = q
                            elif abs(q - ratio) > 1e-12:
                                ok = False
                        elif abs(x) > 1e-12:
                            ok = False
                    if ok:
                        O.append(b)
                        break
            else:
                self.orbits.append([b])
        # parent species belonging to each orbit = those whose inv vector is nonzero and proportional;
        # the full orbit members: species b' with delta^{b'} = delta^b o sigma^k
        self.members = {}
        for O in self.orbits:
            b = O[0]
            mem = set()
            for k in range(self.ord):
                p = compose_pow(self.sigma, k)
                v = [self.d[b][p[j]] for j in range(A.r + 1)]
                for c in A.mass:
                    if abs(self.lam[c] - self.lam[b]) < 1e-12 and max(abs(v[j] - self.d[c][j] * v[0]) for j in range(A.r + 1)) < 1e-12:
                        mem.add(c)
            self.members[b] = sorted(mem)

    def soliton(self, b):
        A = self.A
        mem = self.members[b]
        degs = [len(mem) * int(x) for x in A.n]
        first = [sum(self.d[c][j] for c in mem) for j in range(A.r + 1)]
        t, dr, res = soliton_deg(A, first, mp.sqrt(self.lam[b]), degs=degs, ctx=self.ctx)
        inv_err = max(abs(t[self.sigma[j]][p] - t[j][p]) for j in range(A.r + 1) for p in range(len(t[j])))
        return t, float(dr), float(res), float(inv_err)


def analyse(name, verbose=True):
    F = Fold(name)
    A = F.A
    h = A.h
    rng = np.random.default_rng(2)
    zs = [mp.mpc(*x) for x in (rng.normal(size=(3 * (h + 2), 2)) * 1.5)]
    cq = [mp.cos(mp.pi * q / h) for q in range(h + 1)]
    Xd = X_dorey(A)
    sols = [O[0] for O in F.orbits]
    out = {'name': name, 'parent': A.name, 'h': h, 'central': F.central,
           'species': {str(O[0]): {'members': F.members[O[0]], 'mass': A.mass[O[0]]} for O in F.orbits}, 'X': {}}
    if verbose:
        print(f"== {name}: parent {A.name}, sigma={F.sigma}, central={F.central}, folded species (orbit members):",
              {O[0]: F.members[O[0]] for O in F.orbits})
    worst = 0
    for a in sols:
        t, dr, res, inv_err = F.soliton(a)
        D = [len(tj) - 1 for tj in t]
        if verbose:
            print(f"  soliton {a} (orbit {F.members[a]}): deg={D} dropped={dr:.1e} resid={res:.1e} sigma-inv err={inv_err:.1e}")
        # degenerate invariant channels: group folded channels by mass
        for b in sols:
            degen = [c for c in sols if abs(A.mass[c] - A.mass[b]) < 1e-9]
            Bm = mp.matrix([[F.inv[c][j] for c in degen] for j in range(A.r + 1)])
            Xs, mix = [], 0
            for z in zs:
                s, dr2, res2 = fluct_deg(A, t, mp.sqrt(F.lam[a]), F.inv[b], F.lam[b], z, ctx=F.ctx)
                w = [s[j][D[j]] / t[j][D[j]] for j in range(A.r + 1)]
                coef = mp.lu_solve(Bm.T * Bm, Bm.T * mp.matrix(w))
                resv = max(abs(sum(Bm[j, i] * coef[i] for i in range(len(degen))) - w[j]) for j in range(A.r + 1))
                Xb = coef[degen.index(b)]
                m_ = max([abs(coef[i]) for i, c in enumerate(degen) if c != b] + [resv]) / abs(Xb)
                mix = max(mix, float(m_), float(dr2), float(res2))
                Xs.append(Xb)
            M = np.array([[float(mp.log(abs(z - c))) for c in cq] + [1.0] for z in zs])
            y = np.array([float(mp.log(abs(X))) for X in Xs])
            sol, *_ = np.linalg.lstsq(M, y, rcond=None)
            N = {q: int(round(sol[q])) for q in range(h + 1) if int(round(sol[q]))}
            err = 0
            for z, X in zip(zs, Xs):
                pred = mp.mpf(1)
                for q, m in N.items():
                    pred *= (z - cq[q]) ** m
                err = max(err, float(abs(pred / X - 1)))
            # product formula over the orbit
            P = {}
            for a2 in F.members[a]:
                for q, m in Xd[(a2, b)].items():
                    P[q] = P.get(q, 0) + m
            P = {q: m for q, m in P.items() if m}
            worst = max(worst, mix, err)
            out['X'][f"{a},{b}"] = {'N': {str(q): m for q, m in N.items()}, 'product_ok': P == N, 'mix': mix, 'fit': err}
            if verbose:
                print(f"    channel {b}: X exps {dict(sorted(N.items()))} fit={err:.1e} (H)-mix/resid={mix:.1e} product-formula={'OK' if P == N else 'FAIL ' + str(P)}")
    out['worst'] = worst
    return out, F


if __name__ == '__main__':
    names = sys.argv[1:] or list(FOLDS)
    allout = {}
    for nm in names:
        o, F = analyse(nm)
        allout[nm] = o
        print(nm, 'WORST', o['worst'], flush=True)
    json.dump(allout, open('data/folds_' + ('all' if not sys.argv[1:] else '_'.join(x.split('<')[0] for x in names)) + '.json', 'w'), indent=1)
