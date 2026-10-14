"""Exact fluctuation solutions about every static single soliton: checks of (H) and closed-form X_ab.
Usage: python3 fluct_scan.py d4 d5 e6 ...   Writes results to data/X_<alg>.json"""
import sys, json, os, time
import numpy as np, mpmath as mp
from lie import Algebra, delta_mp
from hirota import MpCtx, soliton_deg, fluct_deg

DPS = int(os.environ.get("DPS", "30"))
os.makedirs('data', exist_ok=True)


def degenerate_basis(A, b, deltas):
    return [c for c in A.mass if abs(A.mass[c] - A.mass[b]) < 1e-8]


def analyse(name):
    typ, r = name[0], int(name[1:])
    A = Algebra(typ, r)
    ctx = MpCtx(DPS)
    h = A.h
    deltas, lams = {}, {}
    for a in A.mass:
        d, lam, _ = delta_mp(A, a, dps=DPS)
        deltas[a], lams[a] = d, lam
    cq = [mp.cos(mp.pi * q / h) for q in range(h + 1)]
    rng = np.random.default_rng(1)
    zs = [mp.mpc(*x) for x in (rng.normal(size=(3 * (h + 2), 2)) * 1.5)]
    out = {'alg': name, 'h': h, 'mass': {str(a): A.mass[a] for a in A.mass}, 'conj': {str(a): A.conj[a] for a in A.mass},
           'pairs': {}}
    worst = {'dropped': 0, 'resid': 0, 'eig': 0, 'mix': 0, 'fit': 0}
    for a in sorted(A.mass):
        t, dr, res = soliton_deg(A, deltas[a], mp.sqrt(lams[a]), ctx=ctx)
        assert dr < 1e-20 and res < 1e-20
        D = [len(tj) - 1 for tj in t]
        for b in sorted(A.mass):
            degen = degenerate_basis(A, b, deltas)
            Bmat = mp.matrix([[deltas[c][j] for c in degen] for j in range(r + 1)])
            Xs, mixes = [], []
            for z in zs:
                s, dr2, res2 = fluct_deg(A, t, mp.sqrt(lams[a]), deltas[b], lams[b], z, ctx=ctx)
                worst['dropped'] = max(worst['dropped'], float(dr2))
                worst['resid'] = max(worst['resid'], float(res2))
                w = [s[j][D[j]] / t[j][D[j]] for j in range(r + 1)]
                # eigen check: NC w = m_b^2 w
                NCw = [sum(int(A.NC[i, j]) * w[j] for j in range(r + 1)) for i in range(r + 1)]
                wn = max(abs(x) for x in w)
                ev = max(abs(NCw[i] - lams[b] * w[i]) for i in range(r + 1)) / (wn * lams[b]) if wn > 0 else 0
                worst['eig'] = max(worst['eig'], float(ev))
                # decompose w in degenerate basis (least squares)
                coef = mp.lu_solve(Bmat.T * Bmat, Bmat.T * mp.matrix(w))
                resv = max(abs(sum(Bmat[j, i] * coef[i] for i in range(len(degen))) - w[j]) for j in range(r + 1))
                Xb = coef[degen.index(b)]
                mix = max([abs(coef[i]) for i, c in enumerate(degen) if c != b] + [0]) / max(abs(Xb), mp.mpf(10) ** (-DPS))
                mixes.append(float(max(mix, resv / max(wn, 1e-300))))
                Xs.append(Xb)
            worst['mix'] = max(worst['mix'], max(mixes))
            # fit log|X| = sum_q N_q log|z - c_q| + const
            M = np.array([[float(mp.log(abs(z - c))) for c in cq] + [1.0] for z in zs])
            y = np.array([float(mp.log(abs(X))) for X in Xs])
            sol, *_ = np.linalg.lstsq(M, y, rcond=None)
            N = np.round(sol[:-1]).astype(int)
            # verify exactly with integer exponents, const fixed by X(z->inf)=1 -> const 0
            err = 0
            for z, X in zip(zs, Xs):
                pred = mp.mpf(1)
                for q, c in enumerate(cq):
                    if N[q]:
                        pred *= (z - c) ** int(N[q])
                err = max(err, abs(pred / X - 1))
            worst['fit'] = max(worst['fit'], float(err))
            out['pairs'][f"{a},{b}"] = {'N': {str(q): int(N[q]) for q in range(h + 1) if N[q]}, 'fit_err': float(err),
                                       'mix': max(mixes), 'D': D}
            print(f"{name} a={a} b={b}: zeros/poles q:N = {{{', '.join(f'{q}:{N[q]}' for q in range(h+1) if N[q])}}} fit={float(err):.1e} mix={max(mixes):.1e}", flush=True)
    out['worst'] = worst
    print(name, 'WORST', worst)
    json.dump(out, open(f'data/X_{name}.json', 'w'), indent=1)


if __name__ == '__main__':
    for nm in sys.argv[1:]:
        t0 = time.time()
        analyse(nm)
        print('time', time.time() - t0)
