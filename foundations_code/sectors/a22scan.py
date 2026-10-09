import numpy as np
from a22 import rep, W
from rsolve import intertwiner
P = np.zeros((9, 9))
for i in range(3):
    for j in range(3): P[j*3+i, i*3+j] = 1
# closed form via projectors: compute R at 3 generic x and get projectors
def projectors(q):
    x0 = 1.37+0.41j
    R = intertwiner(rep(x0, q), rep(1.0, q), (0, 1), W, W)[0]; R = R/R[0, 0]
    ev = np.linalg.eigvals(R); grp = []
    for v in ev:
        if not any(abs(v-g) < 1e-7 for g in grp): grp.append(v)
    Pk = []
    for g in grp:
        M = np.eye(9, dtype=complex)
        for h in grp:
            if h != g: M = M @ (R - h*np.eye(9))/(g - h)
        Pk.append(M)
    return R, grp, Pk
lxs = np.linspace(0, 30, 601)
for om in np.arange(0.02, 2.0, 0.04):
    q = np.exp(-1j*np.pi*om)
    out = []
    for sgn in (1, -1):
        d = 0; arg = None
        for lx in lxs:
            x = sgn*np.exp(lx)
            R = intertwiner(rep(x, q), rep(1.0, q), (0, 1), W, W)[0]; R = R/R[0, 0]
            e = np.linalg.eigvals(P @ R); dd = np.abs(np.abs(e) - 1).max()
            if dd > d: d, arg = dd, lx
        out.append((d, arg))
    print(f"omega {om:.2f}: x>0 maxdev {out[0][0]:.2e} at logx {out[0][1]}, x<0 maxdev {out[1][0]:.2e} at logx {out[1][1]}", flush=True)
