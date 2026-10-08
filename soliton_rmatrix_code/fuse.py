import numpy as np
from dtw import rep, weights
from rsolve import intertwiner
q, n = 0.7, 3
W8 = weights(n); N = 2*n+2; nodes = range(n+1); I8 = np.eye(N)
_cache = {}
def R11(x):
    key = complex(np.round(x, 15))
    if key not in _cache:
        Rc, null, _, _ = intertwiner(rep(n, q, x), rep(n, q, 1.0), nodes, W8, W8)
        _cache[key] = Rc/Rc[0, 0]
    return _cache[key]
def subspace_basis(Rx, wts):
    """weight-adapted orthonormal basis of the image of Rx"""
    U, s, _ = np.linalg.svd(Rx); r = np.sum(s > 1e-8*s[0]); P = U[:, :r] @ U[:, :r].conj().T
    cols, cw = [], []
    for w in np.unique(np.round(wts, 6), axis=0):
        sel = np.where(np.all(np.isclose(wts, w), axis=1))[0]
        if not len(sel): continue
        sub = P[:, sel]; u, ss, _ = np.linalg.svd(sub); k = np.sum(ss > 1e-8)
        for t in range(k): cols.append(u[:, t]); cw.append(w)
    return np.array(cols).T, np.array(cw)
# fused multiplet: image of R11(x*) with x* = q^-2 ; lives in V(a) (x) V(b), a/b = q^2
a, b = q, 1/q
w64 = (W8[:, None, :] + W8[None, :, :]).reshape(N*N, -1)
B, wB = subspace_basis(R11(b/a), w64)          # image of R11(x1/x2) : V(x1)(x)V(x2) -> V(x2)(x)V(x1), x1=b, x2=a
print("fused multiplet dimension:", B.shape[1])
def R12(z):
    M = np.kron(I8, R11(z/b)) @ np.kron(R11(z/a), I8)        # V(z)(x)V(a)(x)V(b) -> V(a)(x)V(b)(x)V(z)
    inp = np.kron(I8, B); outb = np.kron(B, I8)
    img = M @ inp; coef, *_ = np.linalg.lstsq(outb, img, rcond=None)
    resid = np.abs(outb @ coef - img).max()/np.abs(img).max()
    return coef, resid
R, resid = R12(1.37+0.21j)
print(f"R_12 is {R.shape[0]}x{R.shape[1]}; image lies in W (x) V to relative accuracy {resid:.1e}")
def image_dim(zs):
    out = []
    for d in (1e-5, 1e-7):
        R, _ = R12(zs*(1+d)); R = R/np.linalg.norm(R); s = np.linalg.svd(R, compute_uv=False); out.append(int(np.sum(s > 1e-3*s[0])))
    return out
D = 8*B.shape[1]
print(f"special points of R_12 (generic rank {D}):")
for l in np.arange(-9, 9.5, 0.5):
    for sgn in (+1, -1):
        dims = image_dim(sgn*q**l)
        if dims[1] < D and dims[0] == dims[1]:
            print(f"   z* = {'+' if sgn>0 else '-'}q^{l:+.1f}: image dimension {dims[1]}")
np.save('B29.npy', B); np.save('wB29.npy', wB)
