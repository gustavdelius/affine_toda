"""Consolidated checks for the A_2^(1) RSOS (IRF) soliton S-matrix in the symmetric gauge."""
import qsign_patch
from irf_a2 import *
import itertools, sys
from math import gcd
sys.path.insert(0, '/home/gustav/Git/affine_toda/foundations_code')
rng = np.random.default_rng(7)
which = sys.argv[1]

def sigma(m, A, B, u):
    br = m.br
    return br(1+u)*br(1-u)/br(1)**2 if A == B else br(1.5+u)*br(1.5-u)/br(1)**2

if which == 'ybe':
    for (p, pp) in [(6, 7), (7, 9), (7, 12)]:
        m = A2RSOS(p, pp); W = make_weights(m)
        ey = max(ybe_error(m, W, ts, 0.31+0.2j, -0.47+0.1j)[0]/ybe_error(m, W, ts, 0.31+0.2j, -0.47+0.1j)[1]
                 for ts in itertools.product([0, 1], repeat=3))
        eu = 0
        for u in [0.3+0.2j, 0.77j]:
            for a in m.H:
                for c in m.H:
                    for A, B in itertools.product([0, 1], repeat=2):
                        ins, outs, M = block(m, W, A, B, a, c, u)
                        if not ins: continue
                        _, _, M2 = block(m, W, B, A, a, c, -u)
                        eu = max(eu, np.abs(M2 @ M - sigma(m, A, B, u)*np.eye(len(ins))).max())
        print(f"W3({p},{pp}): YBE (8 type triples) rel err {ey:.1e}; braiding unitarity err {eu:.1e}", flush=True)

if which == 'ha':
    for (p, pp) in [(4, 5), (5, 6), (6, 7), (7, 8), (10, 11), (4, 7), (5, 9), (7, 13), (5, 7), (7, 9), (8, 11), (10, 13), (7, 12)]:
        m = A2RSOS(p, pp); W = make_weights(m)
        ha = {}; un = {}; bad = 0
        for th in [0.4, 1.3, -2.1]:
            u = uof(th)
            for a in m.H:
                for c in m.H:
                    for A, B in itertools.product([0, 1], repeat=2):
                        ins, outs, M = block(m, W, A, B, a, c, u)
                        if not ins: continue
                        _, _, M2 = block(m, W, B, A, a, c, -u)   # reversed process at -theta
                        key = 'ss' if A == B else 'sa'
                        ha[key] = max(ha.get(key, 0), np.abs(M - M2.conj().T).max())
                        un[key] = max(un.get(key, 0), np.abs(M @ M.conj().T - sigma(m, A, B, u)*np.eye(len(outs))).max())
        r = pp - p
        neg = sum(m.qdim(h) < 0 for h in m.H)
        print(f"W3({p},{pp}) r={r} r mod p={r % p} lam={r/p:.3f} breathers={'yes' if 3*r > p else 'no'} "
              f"neg qdims {neg}/{len(m.H)}: HA err 33:{ha['ss']:.1e} 33b:{ha['sa']:.1e}; unitarity err 33:{un['ss']:.1e} 33b:{un['sa']:.1e}", flush=True)

if which == 'tm':
    TS = [(0, 1), (0, 0, 0), (0, 1, 0, 1), (0, 0, 1, 1), (0, 0, 0, 1, 1, 1)]
    for (p, pp) in [(4, 5), (5, 6), (6, 7), (7, 8), (4, 7), (5, 9), (5, 7), (7, 9), (8, 11), (7, 12)]:
        f = Fast(A2RSOS(p, pp)); out = []
        for types in TS:
            P, _ = f.per(types)
            if len(P) > 500: continue
            worst = 0
            for t in range(8):
                ev = np.linalg.eigvals(f.transfer(types, rng.normal(size=len(types))*1.5))
                worst = max(worst, np.abs(np.abs(ev)-1).max())
            out.append(f"{''.join('3' if t == 0 else 'b' for t in types)}[{len(P)}]:{worst:.0e}")
        print(f"W3({p},{pp}):", " ".join(out), flush=True)

if which == 'vertex':
    from krein_bethe import R, embed
    for (p, pp) in [(4, 5), (6, 7), (7, 9), (5, 7)]:
        lam = (pp-p)/p; q = qsign_patch.QS*np.exp(1j*np.pi*lam); worst = 0
        for t in range(20):
            th = rng.normal(size=3)*1.5; x = np.exp(1.5*lam*th)
            Tm = np.eye(27, dtype=complex)
            for j in (1, 2): Tm = embed(R(3, x[0]/x[j], q), 0, j, 3, 3) @ Tm
            ids = [np.ravel_multi_index(c, [3]*3) for c in itertools.permutations(range(3))]
            ev = np.linalg.eigvals(Tm[np.ix_(ids, ids)]); worst = max(worst, np.abs(np.abs(ev)-1).max())
        print(f"vertex U_q(sl3^), q={'-' if qsign_patch.QS < 0 else ''}e^(i pi lam), lam={lam:.3f} [W3({p},{pp}) coupling]: 3 solitons, all-colour sector max||s|-1| = {worst:.3f}", flush=True)

if which == 'signs':
    bad = []
    for p in range(4, 61):
        for r in range(1, 2*p):
            if gcd(r, p) != 1: continue
            sg = [np.sign(np.sin(np.pi*r*x/p)) for x in range(1, p)]
            C = all(sg[x-1] == sg[x+1] for x in range(1, p-2))
            if C != (r % p in (1, p-1)): bad.append((p, r))
    print("condition sigma(x+2)=sigma(x) on [1,p-1]  <=>  r = +-1 mod p, for 4<=p<=60, 1<=r<2p: exceptions", bad)
