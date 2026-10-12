"""Part (b): boundary unitarity K(theta)K(-theta)=1 and boundary crossing-unitarity k(theta) = S(2 theta) k(-theta) for the full soliton
reflection matrix on excited states S = |n1,...,nk>:  K^S = [R_0 sigma(x) sigma(i vt)/(cos x cosh vt)] M(x,vt) prod_odd a_ni(eta) prod_even a_ni(etabar),
x = eta (k even) or etabar (k odd, with s <-> sbar).  Book S_0 = a_S0 (validated against the integral representation in boundary_sine_gordon.py)."""
import sys, numpy as np
sys.path.insert(0, '/home/gustav/Git/affine_toda/textbook_code')
from boundary_common import q_phys, sl2_spin, rcheck
from p12b_common import *
def report(name, ok): print(('PASS ' if ok else 'FAIL ') + name, flush=True)
Xm = np.array([[0, 1], [1, 0]])
def Kstate(S, theta, eta, vt, lam):
    k = len(S); L = mp.mpf(lam); x = eta if k % 2 == 0 else etabar(eta, lam)
    u = -1j*mp.mpmathify(theta)
    P = lambda s: mp.cos(L*u)*mp.cos(x)*mp.cosh(vt) - s*mp.sin(L*u)*mp.sin(x)*mp.sinh(vt)
    Q0 = -mp.sin(L*u)*mp.cos(L*u)
    pref = mp.exp(log_R0(u, lam) + log_sigma(x, u, lam) + log_sigma(1j*vt, u, lam))/(mp.cos(x)*mp.cosh(vt))
    for i, n in enumerate(S, start=1):
        pref *= a_n(n, eta if i % 2 == 1 else etabar(eta, lam), u, lam)
    M = np.array([[complex(P(1)), complex(Q0)], [complex(Q0), complex(P(-1))]])
    M = complex(pref)*M
    return M if k % 2 == 0 else Xm @ M @ Xm
def kvec(K): return (K @ Xm).T.reshape(4)
res = []
for lam, eta, vt, th in ((3.3, 5.5, 0.3, 0.29), (6.3, 5.0, 0.8, 0.41)):
    q = q_phys(lam); Rsg = lambda z: rcheck(sl2_spin(1, q, z), sl2_spin(1, q, 1.0))[0]
    for S in [(), (0,), (3,), (0, 5), (2, 1), (1, 3, 0), (2, 0, 1, 4)]:
        K1, K2 = Kstate(S, th, eta, vt, lam), Kstate(S, -th, eta, vt, lam)
        devU = np.abs(K1 @ K2 - np.eye(2)).max()
        km = kvec(Kstate(S, 1j*mp.pi/2 - th, eta, vt, lam)); kp = kvec(Kstate(S, 1j*mp.pi/2 + th, eta, vt, lam))
        S2 = complex(a_S0(-1j*(2*th), lam))*Rsg(np.exp(2*lam*th))
        devC = np.abs(km - S2 @ kp).max()/np.abs(km).max()
        res.append((lam, eta, S, devU, devC))
        print(lam, eta, S, f'unitarity {devU:.1e}  crossing-unitarity {devC:.1e}', flush=True)
report('unitarity and crossing-unitarity of K^S for 14 states (depth <= 4) at 2 parameter points [max %.0e]' % max(max(r[3], r[4]) for r in res),
       max(max(r[3], r[4]) for r in res) < 1e-8)
