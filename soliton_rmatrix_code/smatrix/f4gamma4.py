import numpy as np, time
from collections import defaultdict
exec(open('f4gamma.py').read().split("Ta = (6, F(9, 2)); s =")[0])
def analyse(label, Tf, ours_ex, cds_ex):
    t0 = time.time(); T = val(Tf, w0); bet = np.pi*(T - 1)/T
    # exact needed exponents on the lattice (a, b): CDS / ours
    def canon(a, b):
        a, b = F(a), F(b)
        while val((a, b), w0) > T + 1e-9: a, b = a - 2*Tf[0], b - 2*F(Tf[1])
        while val((a, b), w0) <= -T + 1e-9: a, b = a + 2*Tf[0], b + 2*F(Tf[1])
        return (a, b)
    need = defaultdict(int)
    for y in cds_ex[0]: need[canon(*y)] += 1; need[canon(-y[0], -y[1])] -= 1
    for y in ours_ex[0]: need[canon(*y)] -= 1; need[canon(-y[0], -y[1])] += 1
    need = {k: v for k, v in need.items() if v}
    grid = [(al, F(b2, 2)) for al in range(0, 7) for b2 in range(-6, 13) if 0 < val((al, F(b2, 2)), w0) < T/2]
    cols = []; checked = 0; devs = []
    for (al, bb) in grid:
        e = val((al, bb), w0); Fp = Fmaker([e, T - e], T, J=3000)
        Wf = lambda th: Fp(th + 1j*bet)*Fp(th)**2*Fp(th - 1j*bet)
        if checked < 1:
            th = 0.3 + 0.2j; t = T*th/(1j*np.pi)
            alt = Fp(th)**2/(Fp(1j*np.pi*(t + 1)/T)*Fp(1j*np.pi*(t - 1)/T))
            print(f"   identity check W = F(t)^2/(F(t+1)F(t-1)): ratio {Wf(th)/alt:.6f}"); checked += 1
        cands = set()
        for s_ in (1, -1):
            for base in ((al, bb), (Tf[0] - al, F(Tf[1]) - bb)):
                for k in range(-3, 4): cands.add(canon(s_*base[0], s_*base[1] + k))
        for k in range(-3, 4): cands.add(canon(0, k))
        col = defaultdict(int)
        for c in cands:
            tc = val(c, w0)
            if abs(tc) < 1e-9 or abs(abs(tc) - T) < 1e-9: continue
            n_ = [abs(Wf(1j*np.pi*(tc + d)/T)) for d in (1e-4, 1e-6)]
            od = np.log10(n_[1]/n_[0])/2
            if abs(od) > 0.5: col[c] += int(round(od))
        # store as block exponents: a pole of order o at t = y means block (y)^o ; zeros at y mean block (-y)
        blocks = defaultdict(int)
        for c, o in col.items():
            if o > 0: blocks[c] += o
            else: blocks[canon(-c[0], -c[1])] += -o
        cols.append(blocks)
        mdl = lambda th: np.prod([(np.sinh(th/2 + 1j*np.pi*val(k, w0)/(2*T))/np.sinh(th/2 - 1j*np.pi*val(k, w0)/(2*T)))**v for k, v in blocks.items()]) if blocks else 1.0
        dev = max(abs(Wf(th)/mdl(th) - 1) for th in (0.31+0.17j, -0.4+0.5j))
        devs.append(dev)
    keys = sorted(set(need) | set(k for c in cols for k in c), key=lambda y: val(y, w0))
    M = np.array([[c.get(k, 0) for c in cols] for k in keys], float)
    rhs = np.array([need.get(k, 0) for k in keys], float)
    sol, *_ = np.linalg.lstsq(M, rhs, rcond=None)
    res = np.abs(M @ sol - rhs).max()
    print(f"{label}: {len(grid)} pairs, {len(keys)} lattice positions; exact block problem residual {res:.3f}, rank {np.linalg.matrix_rank(M)}  [{time.time()-t0:.0f}s]")
    print(f"   every W_e reproduced by its identified block product: max deviation {max(devs):.1e} (truncation level)")
    ex = cols[grid.index((1, F(1, 2)))] if (1, F(1, 2)) in grid else cols[0]
    print(f"   example: W_e for e = w + 1/2 has blocks {sorted([(str(k[0]) + 'w' + ('+' if k[1] >= 0 else '') + str(k[1]), v) for k, v in ex.items()], key=lambda z: z[0])}")
    return res
Ta = (6, F(9, 2)); s = (F(3), (F(9, 2) - 1)/2)
oursA = normalize([(y[0] + s[0], y[1] + s[1]) for y in SBa] + [(y[0] - s[0], y[1] - s[1]) for y in SBa], 1, Ta)
analyse("(a)", Ta, oursA, cds_exact(Ta, 1))
Tb = (6, F(-3, 2))
analyse("(b)", Tb, normalize(BBb, int(sgn_b), Tb), cds_exact(Tb, 0))
