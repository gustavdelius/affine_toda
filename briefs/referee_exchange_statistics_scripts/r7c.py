import sys, os, itertools, numpy as np
here = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(here, 'r3_L_bicharacter.py')).read().split('for n, alg in')[0])
for n, alg in ((3, 'c'), (3, 'a')):
    for om in (0.3, 0.7, 1.61):
        try:
            D, Wn, bits, L = Lmatrix(n, alg, om)
        except Exception as e:
            print(alg, n, om, 'no diagonal +-1 similarity:', e); continue
        Pi = np.array([(-1)**b.sum() for b in bits])
        # test L ~ Pi x 1 modulo function of total weight
        ok = True; tot = {}
        for a in range(D):
            for b in range(D):
                val = L[a, b]*Pi[a]
                if tot.setdefault(tuple(bits[a]+bits[b]), val) != val: ok = False
        # general fit search
        cnt = 0; first = None
        for u in itertools.product((0, 1), repeat=n):
            for v in itertools.product((0, 1), repeat=n):
                for Mf in itertools.product((0, 1), repeat=n*n):
                    M = np.array(Mf).reshape(n, n); u_, v_ = np.array(u), np.array(v); good = True; t2 = {}
                    for a in range(D):
                        for b in range(D):
                            val = L[a, b]*(-1)**(bits[a]@M@bits[b] + u_@bits[a] + v_@bits[b])
                            if t2.setdefault(tuple(bits[a]+bits[b]), val) != val: good = False; break
                        if not good: break
                    if good:
                        cnt += 1
                        if first is None: first = (M.tolist(), u, v)
        print(f"{alg}_{n} om={om}: L ~ Pi x 1: {ok}; general fits {cnt}; first {first}", flush=True)
