"""Factories for fluctuation sectors: parent single solitons and folded (invariant) sectors."""
import numpy as np
import mpmath as mp
from lie import Algebra, delta_mp
from hirota import MpCtx, soliton_deg
from sector import Sector
from compare_dorey import X_dorey

_cache = {}


def algebra(nm):
    if nm not in _cache:
        A = Algebra(nm[0], int(nm[1:]))
        A.Xd = X_dorey(A)
        _cache[nm] = A
    return _cache[nm]


def parent_sector(nm, a, dps=30, xi=None):
    A = algebra(nm)
    ctx = MpCtx(dps)
    chans = []
    for b in sorted(A.mass):
        d, lam, _ = delta_mp(A, b, dps=dps)
        chans.append((b, d, lam))
    da = [c for c in chans if c[0] == a][0]
    t, dr, res = soliton_deg(A, da[1], mp.sqrt(da[2]), ctx=ctx)
    S = Sector(A, t, mp.sqrt(da[2]), chans, xi=xi, dps=dps, label=f"{nm} soliton {a}")
    S.exps = {b: A.Xd[(a, b)] for b in sorted(A.mass)}
    return S


def fold_sector(foldname, a, dps=30, xi=None):
    import fold as F_
    F_.DPS = dps
    F = F_.Fold(foldname)
    A = F.A
    A.Xd = X_dorey(A)
    t, dr, res, inv_err = F.soliton(a)
    assert dr < 1e-15 and res < 1e-15 and inv_err < 1e-15
    sols = [O[0] for O in F.orbits]
    chans = [(b, F.inv[b], F.lam[b]) for b in sols]
    # invariant subspace projector
    g = F.sigma
    R = A.autoR[A.autos.index(g)]
    # invariant subspace of R (orthogonal map): eigenvectors with eigenvalue 1
    w, V = np.linalg.eig(R)
    inv = V[:, np.abs(w - 1) < 1e-8]
    Q, _ = np.linalg.qr(np.real(inv) if np.allclose(inv.imag, 0) else inv)
    # real orthonormal basis
    Pr = np.linalg.svd(np.hstack([inv.real, inv.imag]))[0][:, :inv.shape[1]]
    S = Sector(A, t, mp.sqrt(F.lam[a]), chans, P=Pr, xi=xi, dps=dps, label=f"{foldname} soliton {a}")
    exps = {}
    for b in sols:
        P = {}
        for a2 in F.members[a]:
            for q, m in A.Xd[(a2, b)].items():
                P[q] = P.get(q, 0) + m
        exps[b] = {q: m for q, m in P.items() if m}
    S.exps = exps
    S.fold = F
    return S


def Xclosed(S, nm):
    h = S.A.h
    N = S.exps[nm]
    def X(z):
        out = 1.0 + 0j
        for q, m in N.items():
            out *= (z - np.cos(np.pi * q / h)) ** m
        return out
    return X


def net_isolated(S, tol=1e-9):
    """net count (zeros - poles) of prod_b X_b at energies below the lowest threshold (z in (0,1])"""
    h = S.A.h
    mmin2 = min(S.mb2.values())
    ent = []
    for nm, d, lb in S.channels:
        for q, m in S.exps[nm].items():
            if 2 * q < h:
                lam = S.mb2[nm] * np.sin(np.pi * q / h) ** 2
                ent.append((lam, m, nm, q))
    ent.sort()
    groups = []
    for e in ent:
        if groups and abs(groups[-1][0] - e[0]) < tol:
            groups[-1][1] += e[1]
            groups[-1][2].append(e)
        else:
            groups.append([e[0], e[1], [e]])
    iso = [(g[0], g[1], g[2]) for g in groups if g[0] < mmin2 - 1e-9]
    emb = [(g[0], g[1], g[2]) for g in groups if g[0] >= mmin2 - 1e-9]
    return iso, emb, mmin2
