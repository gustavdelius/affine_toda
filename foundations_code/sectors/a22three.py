import numpy as np
from lib import *
from a22 import rep, W
from rsolve import intertwiner
P9 = flip(3, 3)
def setup(q):
    e, f, k, qi = rep(1.0, q)
    comps = components({1: e[1]}, {1: f[1]}, {1: k[1]}, W, [1], {1: np.array([1.0])})
    byd = {d: Pm for d, Pm, _ in comps}
    def Rpl(x):
        Rc = byd[5] + br(x, 4, 1, q)*byd[3] + br(x, 6, -1, q)*byd[1]
        return P9 @ Rc
    return Rpl
q = np.exp(-1j*0.77); Rpl = setup(q); x0 = 1.9+0.3j
Rc = intertwiner(rep(x0, q), rep(1.0, q), (0, 1), W, W)[0]; Rc /= Rc[0, 0]
print("closed form check", np.abs(P9 @ Rc - Rpl(x0)).max())
lx = np.concatenate([np.linspace(0.002, 8, 400), np.linspace(8, 40, 65)])
for xi_pi in (0.3, 0.35, 0.4, 0.43, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.3, 1.34, 1.4, 1.5, 2.0, 3.0, 6.0):
    xi = (xi_pi + 0.0013)*np.pi; qTW = 1j*np.exp(1j*np.pi**2/(3*xi)); q = np.conj(qTW); Rpl = setup(q)
    two = maxdev_scan(Rpl, lx, (1,)); twom = maxdev_scan(Rpl, lx, (-1,))
    three, kr = scan3(Rpl, 3, L=8, n=400)
    c = np.cos(2*np.pi**2/(3*xi))
    print(f"xi={xi_pi:5.2f} pi  cos(2pi^2/3xi)={c:+.3f}:  K0K0 2-body {two[0]:.1e} [logx {two[2]:.2f}];  R(-x) 2-body {twom[0]:.1e};  K0K0K0 3-body {three[0]:.1e} (Krein err {kr:.0e})", flush=True)
