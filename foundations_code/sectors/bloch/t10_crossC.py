"""c_n spinor: charge-conjugation intertwiner C: V*(x) = V(x t) (antipode as in soliton_rmatrix_code/smatrix/crossing.py); symmetry of C."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')); from lib import spinor_rep
om = float(sys.argv[1]) if len(sys.argv) > 1 else 2.37; q = np.exp(-1j*np.pi*om); x = 0.83+0.29j
for n in (1, 2, 3, 4):
    D = 2**n
    def dual(x):
        e, f, k, _ = spinor_rep(n, q, x, 'c'); out = []
        for i in range(n+1):
            ki = np.linalg.inv(k[i]); out += [(-ki@e[i]).T, (-f[i]@k[i]).T, ki.T]
        return out
    def plain(x):
        e, f, k, _ = spinor_rep(n, q, x, 'c'); out = []
        for i in range(n+1): out += [e[i], f[i], k[i]]
        return out
    for ph in (1, -1):
        hit = None
        for l in range(-2*n-4, 2*n+5):
            t = ph*q**l; A, B = dual(x), plain(x*t)
            M = np.vstack([np.kron(np.eye(D), Ag.T) - np.kron(Bg, np.eye(D)) for Ag, Bg in zip(A, B)])
            s = np.linalg.svd(M, compute_uv=False)
            if s[-1] < 1e-9*s[0]: hit = (l, t, M); break
        if hit: break
    l, t, M = hit
    C = np.linalg.svd(M)[2][-1].conj().reshape(D, D)       # C A = B C  convention of crossing.py
    Mx = np.linalg.inv(C) @ C.T
    off = np.abs(Mx - np.diag(np.diag(Mx))).max()
    d = np.diag(Mx)
    print(f"n={n}: t={'+' if ph>0 else '-'}q^{l}; C^-1 C^T diagonal (offdiag {off:.0e}); diag/q-powers:",
          np.round(d.real, 3))
    # q -> 1 limit: repeat at q close to 1
