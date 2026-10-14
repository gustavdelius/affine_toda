import numpy as np
import b2_dn_classical as b2
sm = b2.load_smatrix()
n = 5; h, bl, pfm, _ = sm.setup('d', n)
al, mk = b2.dn_roots(n); masses = np.sqrt(np.linalg.eigvalsh((al.T*mk) @ al)); mmin = masses.min()
for a in range(n):
    print(a, round(pfm[a]*mmin, 4), bl[(a, 2)] if a == 2 else '', flush=True)
a = 2
print('blocks of node 2 with all b:', {b: bl[(a, b)] for b in range(n)})
for e in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
    v = np.prod([sm.S(bl[(a, b)], h, 1.0, 1j*np.pi*3/h + e) for b in range(n)])
    print(f"eps {e:.0e}: |prod_b S_2b(i 3pi/8 + eps, B=1)| = {abs(v):.6e}", flush=True)
