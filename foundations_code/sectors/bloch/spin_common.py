"""Spinor-soliton R-matrices (c_n^(1): alg='c'; a_{2n-1}^(2): alg='a') from ../spinor_scan.py, as a module."""
import sys, os, math, itertools, numpy as np
SECTORS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, SECTORS)
from lib import *
def make(N, alg):
    D = 2**N; Wn = spin_weights(N); P = flip(D, D)
    roots = {i: np.eye(N)[i-1] - np.eye(N)[i] for i in range(1, N)}; roots[N] = np.eye(N)[N-1]
    def rhos(x, q):
        if alg == 'c':
            out = [1.0]
            for j in range(1, N+1): out.append(out[-1]*br(x, 2*j, (-1)**(j+1), q))
        else:
            out = []
            for k in range(N+1):
                v = 1.0
                for i in range(1, (k+1)//2 + 1): v = v*br(x, 4*k - 8*i + 6, 1, q)
                out.append(v)
        return out
    def setup(q):
        e, f, k, qi = spinor_rep(N, q, 1.0, alg)
        comps = components(e, f, k, Wn, range(1, N+1), roots)
        byd = {d: Pm for d, Pm, _ in comps}
        Pk = [byd[math.comb(2*N+1, N-kk)] for kk in range(N+1)]
        return lambda x: P @ sum(r*Pm for r, Pm in zip(rhos(x, q), Pk))   # particle-labelled R = P Rcheck
    return D, Wn, setup
def weights_multi(Wn, Nn):
    """total weight vector and weight of particle 1 for each basis state of V^{x Nn}"""
    D = Wn.shape[0]; tot = []; first = []
    for c in itertools.product(range(D), repeat=Nn):
        tot.append(sum(Wn[i] for i in c)); first.append(Wn[c[0]])
    return np.array(tot), np.array(first)
