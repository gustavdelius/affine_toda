import numpy as np, time
exec(open('f4gamma.py').read().split("Ta = (6, F(9, 2)); s =")[0])
thr = np.linspace(0.06, 4.0, 60); h = 1e-5
def dlog(fun):
    return np.array([((fun(t + h)/fun(t - h))).__abs__() and np.angle(fun(t + h)/fun(t - h))/(2*h) for t in thr])
def solve(label, Tf, ours_ex, cds_ex):
    t0 = time.time()
    T = val(Tf, w0); bet = np.pi*(T - 1)/T
    ours_y = [val(y, w0) for y in ours_ex[0]]; cds_y = [val(y, w0) for y in cds_ex[0]]
    need = lambda th: blocks_val(cds_y, cds_ex[1], T, th)/blocks_val(ours_y, ours_ex[1], T, th)
    rhs = dlog(need)
    grid = sorted({round(al*w0 + b2/2, 9) for al in range(0, 7) for b2 in range(-6, 13) if 0 < al*w0 + b2/2 < T/2})
    cols = []
    for e in grid:
        Fp = Fmaker([e, T - e], T, J=1200)
        W = lambda th, Fp=Fp: Fp(th + 1j*bet)*Fp(th)**2*Fp(th - 1j*bet)
        cols.append(dlog(W))
    M = np.array(cols).T
    sol, *_ = np.linalg.lstsq(M, rhs, rcond=None)
    res = np.abs(M @ sol - rhs).max()/np.abs(rhs).max()
    nz = sorted([(abs(v), e, v) for e, v in zip(grid, sol) if abs(v) > 0.05], reverse=True)[:8]
    print(f"{label}: {len(grid)} candidate pairs; best least-squares relative residual {res:.3e}  [{time.time()-t0:.0f}s]")
    print("   largest multiplicities:", [(round(e, 3), round(v, 3)) for _, e, v in nz])
    # sanity: the needed factor itself is far from trivial
    print(f"   size of the needed correction: max |d log W_needed / d theta| = {np.abs(rhs).max():.3f}")
Ta = (6, F(9, 2)); s = (F(3), (F(9, 2) - 1)/2)
oursA = normalize([(y[0] + s[0], y[1] + s[1]) for y in SBa] + [(y[0] - s[0], y[1] - s[1]) for y in SBa], 1, Ta)
solve("(a)", Ta, oursA, cds_exact(Ta, 1))
Tb = (6, F(-3, 2))
solve("(b)", Tb, normalize(BBb, int(sgn_b), Tb), cds_exact(Tb, 0))
