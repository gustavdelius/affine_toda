"""Closed-form transmission factors X_ab(z), z = kappa_b/m_b, from the (Hirota-verified) block formula, and
the zero/pole bookkeeping of Proposition 8.5: net order of prod_b X_b at each energy lambda.

For folded theories, X for a folded channel is the product over the soliton orbit (see fold.py)."""
import numpy as np
from fractions import Fraction
from lie import Algebra
from compare_dorey import X_dorey


def Xfun(N, h):
    """callable X(z) from exponent dict q -> N_q"""
    cq = {q: np.cos(np.pi * q / h) for q in N}
    def X(z):
        z = np.asarray(z, complex)
        out = np.ones_like(z)
        for q, m in N.items():
            out = out * (z - cq[q]) ** m
        return out
    return X


def zeros_poles(N, h, mb2):
    """list of (lambda, q, order, kind) for z = cos(pi q/h) >= 0, i.e. bound-state half-plane and threshold"""
    out = []
    for q, m in N.items():
        c = np.cos(np.pi * q / h)
        if q * 2 < h:
            lam = mb2 * (1 - c * c)
            out.append((lam, q, m, 'bs'))
        elif q * 2 == h:
            out.append((mb2, q, m, 'thr'))
    return out


def net_count(entries, tol=1e-9):
    """entries: list of (lambda, b, q, order) -> dict of distinct lambda -> (net, details)"""
    groups = []
    for e in sorted(entries):
        for g in groups:
            if abs(g[0][0] - e[0]) < tol:
                g.append(e)
                break
        else:
            groups.append([e])
    return [(g[0][0], sum(e[3] for e in g), g) for g in groups]


def soliton_table(A, a, Xd=None, channels=None, product_over=None):
    """bound-state-region bookkeeping for soliton a (or orbit product_over) of algebra A"""
    Xd = Xd or X_dorey(A)
    chans = channels or sorted(A.mass)
    orbit = product_over or [a]
    entries, thr = [], []
    Ns = {}
    for b in chans:
        N = {}
        for a2 in orbit:
            for q, m in Xd[(a2, b)].items():
                N[q] = N.get(q, 0) + m
        N = {q: m for q, m in N.items() if m}
        Ns[b] = N
        for lam, q, m, kind in zeros_poles(N, A.h, A.mass[b] ** 2):
            if kind == 'bs':
                entries.append((lam, b, q, m))
            else:
                thr.append((b, q, m))
    return Ns, net_count(entries), thr
