import numpy as np, time
from dtw import rep, weights
from rsolve import intertwiner
q, n = 0.7, 3; W8 = weights(n); B = np.load('B29.npy'); a, b = q, 1/q
cache = {}
def R11(x):
    k = complex(np.round(x, 14))
    if k not in cache:
        Rc = intertwiner(rep(n, q, x), rep(n, q, 1.0), range(n+1), W8, W8)[0]; cache[k] = Rc/np.abs(Rc).max()
    return cache[k]
Win = np.einsum('ia,jb->ijab', B, B).reshape(8**4, 29*29)      # basis of W (x) W inside V^(x)4, index order (1,2,3,4)
def apply_local(vecs, M, pos):
    # vecs: (8,8,8,8,K); apply 64x64 M on tensor factors (pos, pos+1)
    v = np.moveaxis(vecs, [pos, pos+1], [0, 1]); sh = v.shape
    v = (M @ v.reshape(64, -1)).reshape(sh); return np.moveaxis(v, [0, 1], [pos, pos+1])
def R22(z):
    v = Win.reshape(8, 8, 8, 8, -1)
    v = apply_local(v, R11(z*b/a), 1)      # [za, zb, a, b] -> [za, a, zb, b]
    v = apply_local(v, R11(z*b/b), 2)      #                -> [za, a, b, zb]
    v = apply_local(v, R11(z*a/a), 0)      #                -> [a, za, b, zb]
    v = apply_local(v, R11(z*a/b), 1)      #                -> [a, b, za, zb]
    img = v.reshape(8**4, -1); coef, *_ = np.linalg.lstsq(Win, img, rcond=None)
    return coef, np.abs(Win @ coef - img).max()/np.abs(img).max()
if __name__ == "__main__":
    C, res = R22(1.21+0.33j)
    print(f"R_22 by fusion: image in W (x) W to {res:.1e}")
