import numpy as np
from toda import *
pi = np.pi
def coxeter(alg, r):
    return {'a': r+1, 'd': 2*r-2, 'e': {6: 12, 7: 18, 8: 30}.get(r)}[alg]

def fitX(C, n, c, ma, mb, vb, h, npts=60, seed=0, half=False):
    rng = np.random.default_rng(seed)
    zs = rng.uniform(-1.5, 1.5, npts) + 1j*rng.uniform(0.05, 1.5, npts)*rng.choice([-1, 1], npts)
    Xs = []; res = 0
    for z in zs:
        X, r_, D, Dp, p = Xfactor(C, n, c, ma, mb, vb, z)
        Xs.append(X); res = max(res, r_)
    Xs = np.array(Xs)
    grid = np.arange(0, 2*h+1)/2.0 if half else np.arange(0, h+1)*1.0   # x grid
    cs = np.cos(pi*grid/h)
    A = np.log(np.abs(zs[:, None] - cs[None, :]))
    e, *_ = np.linalg.lstsq(A, np.log(np.abs(Xs)), rcond=None)
    ei = np.round(e).astype(int)
    pred = np.prod((zs[:, None] - cs[None, :])**ei[None, :], axis=1)
    err = np.max(np.abs(pred/Xs - 1))
    roots = [(grid[k], ei[k]) for k in range(len(grid)) if ei[k] != 0]
    return roots, err, res, np.max(np.abs(e-ei))

def run(alg, r, verbose=True):
    h = coxeter(alg, r)
    sp, C, n, autos = species(alg, r)
    sols = []
    for ma, va in sp:
        k0 = np.argmax(np.abs(va)); va = va/va[0] if abs(va[0]) > 1e-8 else va/va[k0]
        c = soliton(C, n, ma, va)
        sols.append((ma, va, c, degree(c)))
    table = {}
    for ia, (ma, va, c, D) in enumerate(sols):
        for ib, (mb, vb, _, _) in enumerate(sols):
            roots, err, res, dev = fitX(C, n, c, ma, mb, vb, h)
            table[(ia, ib)] = roots
            if verbose:
                print(f"  sol {ia} (m={ma:.5f}, deg {D}) ch {ib} (m={mb:.5f}): roots x[ord] =",
                      " ".join(f"{x:g}[{o:+d}]" for x, o in roots), f" fit err {err:.1e}, diag resid {res:.1e}, int dev {dev:.1e}")
    return sols, table, h, C, n, autos

if __name__ == "__main__":
    import sys
    alg, r = sys.argv[1], int(sys.argv[2])
    run(alg, r)
