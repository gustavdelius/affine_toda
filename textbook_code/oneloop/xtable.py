import numpy as np, sys, pickle
from toda import *
from xfit import fitX, coxeter
def run(alg, r, verbose=True):
    h = coxeter(alg, r)
    sp, C, n, autos = species(alg, r)
    sols = []
    for ma, va in sp:
        k0 = np.argmax(np.abs(va)); va = va/va[k0]
        c = soliton(C, n, ma, va)
        sols.append((ma, va, c, degree(c)))
    table = {}; worst = 0
    for ia, (ma, va, c, D) in enumerate(sols):
        for ib, (mb, vb, _, _) in enumerate(sols):
            roots, err, res, dev = fitX(C, n, c, ma, mb, vb, h)
            table[(ia, ib)] = roots; worst = max(worst, err, res)
    if verbose:
        print(f"{alg}_{r}^(1), h = {h}; masses:", np.round([s[0] for s in sols], 5))
        print("  tau degrees:", [s[3] for s in sols], " (Kac labels", list(n.astype(int)), ")")
        for ia in range(len(sols)):
            print(f"  sol {ia}: " + " | ".join(" ".join(f"{x:g}{'+' if o>0 else '-'}{abs(o) if abs(o)>1 else ''}" for x, o in table[(ia, ib)]) for ib in range(len(sols))))
        print("  worst fit/diagonality error:", f"{worst:.1e}")
    return dict(alg=alg, r=r, h=h, masses=[s[0] for s in sols], vecs=[s[1] for s in sols],
                taus=[s[2] for s in sols], degs=[s[3] for s in sols], table=table, C=C, n=n)
if __name__ == "__main__":
    alg, r = sys.argv[1], int(sys.argv[2])
    out = run(alg, r)
    pickle.dump(out, open(f"X_{alg}{r}.pkl", "wb"))
