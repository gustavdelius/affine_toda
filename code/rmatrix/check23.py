import numpy as np
from spinor import spin_rep, Rgen, Ws, q, n
from dtw import rep, weights
Wv = weights(n); B = np.load('B29.npy'); a, b = q, 1/q; I8 = np.eye(8)
cache = {}
def R13(x):
    k = complex(np.round(x, 14))
    if k not in cache: cache[k] = Rgen(rep(n, q, x), spin_rep(1.0), Wv, Ws)[0]
    return cache[k]
def rank_at(fun, z, D):
    dims = []
    for d in (1e-5, 1e-7):
        R = fun(z*(1+d)); R = R/np.linalg.norm(R); s = np.linalg.svd(R, compute_uv=False); dims.append(int(np.sum(s > 1e-3*s[0])))
    return dims[1] if dims[0] == dims[1] else None
# (i) R_13 on the real axis, both signs
hits = []
for l in np.arange(-7, 7.5, 0.5):
    for sg, lab in ((1, '+'), (-1, '-')):
        d = rank_at(R13, sg*q**l, 64)
        if d is not None and d < 64: hits.append(f"{lab}q^{l:+.1f}->{d}")
print("R_13 real axis:", "; ".join(hits) if hits else "none")
# (ii) R_23 by fusion: W(a,b) (x) S(c) -> S(c) (x) W ;  z = 1/c
def R23(z):
    c = 1/z
    M = np.kron(R13(a/c), I8) @ np.kron(I8, R13(b/c))          # V(a)V(b)S(c) -> S(c)V(a)V(b)
    img = M @ np.kron(B, I8); coef, *_ = np.linalg.lstsq(np.kron(I8, B), img, rcond=None)
    return coef
c0 = R23(1.31+0.27j)
res = np.abs(np.kron(I8, B) @ c0 - (np.kron(R13(a*(1.31+0.27j)), I8) @ np.kron(I8, R13(b*(1.31+0.27j))) @ np.kron(B, I8))).max()
print(f"R_23 (29 x 8 = 232) by fusion; image in S (x) W to {res:.1e}")
hits = []
for l in np.arange(-6, 7, 1.0):
    for ph, lab in ((1, '+'), (-1, '-'), (1j, '+i'), (-1j, '-i')):
        d = rank_at(R23, ph*q**l, 232)
        if d is not None and d < 232: hits.append(f"{lab}q^{l:+.0f}->{d}")
print("R_23 special points:", "; ".join(hits) if hits else "none")
