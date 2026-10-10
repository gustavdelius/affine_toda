import sys, os, itertools, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'r3_L_bicharacter.py')).read().split('for n, alg in')[0])
for n, alg in ((2, 'c'), (2, 'a'), (3, 'c')):
    D, Wn, bits, L = Lmatrix(n, alg, 1.61)
    print(alg, n, 'bits', bits.tolist()); print(L)
    # blocks by total weight: within each block, values of L
    cnt = 0
    for Mf in itertools.product((0, 1), repeat=n*n):
        M = np.array(Mf).reshape(n, n)
        for u in itertools.product((0, 1), repeat=n):
            for v in itertools.product((0, 1), repeat=n):
                u_, v_ = np.array(u), np.array(v); ok = True; tot = {}
                for a in range(D):
                    for b in range(D):
                        val = L[a, b]*(-1)**(bits[a]@M@bits[b] + u_@bits[a] + v_@bits[b])
                        if tot.setdefault(tuple(bits[a]+bits[b]), val) != val: ok = False; break
                    if not ok: break
                if ok:
                    cnt += 1
                    if cnt <= 3: print('  fit M', M.tolist(), 'u', u, 'v', v)
    print('  general fits:', cnt)
