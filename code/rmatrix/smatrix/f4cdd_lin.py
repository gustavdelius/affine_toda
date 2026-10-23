import numpy as np
from collections import defaultdict
exec(open('f4particle2.py').read().split("Y = lambda U: U*T/H")[0])
ours = Amp.blocks(BBs, BBsign); c11 = cds_blocks(1, 1, H)
need = c11*Amp(list(ours.den.elements()), list(ours.num.elements()), ours.sign)   # target: CDS / ours
w = omega
def to_ab(y):
    for a in range(-24, 25):
        b = y - a*w
        if abs(b*4 - round(b*4)) < 1e-6 and abs(b) < 20: return canon(a, round(b*4)/4)
    raise ValueError(y)
def canon(a, b):
    while a > 5: a -= 12; b -= 3
    while a < -6: a += 12; b += 3
    return (a, b)
def sh(cell, da, db): return canon(cell[0] + da, cell[1] + db)
target = defaultdict(int)
for c in need.num.elements(): target[to_ab(c)] += 1
for c in need.den.elements(): target[to_ab(c)] -= 1
print("target (CDS/ours) exponents at lattice positions (a, b) <-> a*w + b:", dict(target))
sig = (6, 0.5)                                    # T - 1 = 6w + 1/2
tri = (2, 0.5)                                    # T/3 = 2w + 1/2
def solve(bmax=8.0, with_boot=True):
    cells = [(a, b) for a in range(-6, 6) for b in np.arange(-bmax, bmax + 0.01, 0.25)]
    idx = {(a, round(b, 2)): i for i, (a, b) in enumerate(cells)}
    rows_cells = [(a, b) for a in range(-6, 6) for b in np.arange(-bmax - 5, bmax + 5.01, 0.25)]
    A1 = np.zeros((len(rows_cells), len(cells))); A2 = np.zeros_like(A1); rhs = np.zeros(len(rows_cells))
    for r, (a, b) in enumerate(rows_cells):
        key = (a, round(b, 2)); rhs[r] = target.get(key, 0)
        for (da, db), cf in (((0, 0), 2), (sig, 1), ((-sig[0], -sig[1]), 1)):
            src = sh(key, -da, -db); j = idx.get((src[0], round(src[1], 2)))
            if j is not None: A1[r, j] += cf
        for (da, db), cf in (((0, 0), 1), (tri, -1), ((-tri[0], -tri[1]), -1)):
            src = sh(key, -da, -db); j = idx.get((src[0], round(src[1], 2)))
            if j is not None: A2[r, j] += cf
    missing = [k for k in target if k not in {(a, round(b, 2)) for a, b in rows_cells}]
    A = np.vstack([A1, A2]) if with_boot else A1
    y = np.concatenate([rhs, np.zeros(len(rows_cells))]) if with_boot else rhs
    sol, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = np.abs(A @ sol - y).max()
    return sol, res, cells, missing
for wb in (False, True):
    for bmax in (6.0, 8.0):
        sol, res, cells, missing = solve(bmax, wb)
        print(f"breather condition{' + self-fusion bootstrap' if wb else ''}, window |b| <= {bmax}: best residual {res:.3f}" + (f" (target cells outside rows: {missing})" if missing else ""))
print("--- larger windows")
for bmax in (12.0, 16.0, 24.0):
    for wb in (False, True):
        sol, res, cells, missing = solve(bmax, wb)
        print(f"   |b| <= {bmax:4.0f}, bootstrap {wb}: residual {res:.4f}")
sol, res, cells, missing = solve(24.0, True)
big = sorted([(abs(v), c, round(v, 3)) for v, c in zip(sol, cells) if abs(v) > 0.05], reverse=True)[:24]
print("largest exponents of the least-squares Z (|b|<=24, with bootstrap):")
for _, c, v in big: print(f"   n({c[0]:+d} w {c[1]:+.2f}) = {v:+.3f}")
