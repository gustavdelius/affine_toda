"""Part (b): breathers B_1, B_2 on the excited states.  Full matrix bootstrap  v = (1 x K^S(theta_a)) S(theta_a+theta_b) (1 x K^S(theta_b)) w,
w the U_q-invariant (breather) vector in V(theta_a) x V(theta_b), theta_{a,b} = theta -+ i u_m/2, u_m = pi - m pi/lambda, against
   K_{B_m}^S(u) = K_{B_m}(x_S; u) prod_{i odd} b_{m n_i}(eta;u) prod_{i even} b_{m n_i}(etabar;u),
   b_{mk}(x;u) = a_k(x; u+pi/2-m pi/(2 lam)) a_k(x; u-pi/2+m pi/(2 lam))
(BPTT (3.9) and the last formula of sec. 3.7); x_S = eta (even depth), etabar (odd depth); K_{B_m}(x;u) = Ghoshal's B_m factor.
Also checks that the vector v is proportional to w (the breather decouples from the soliton multiplet)."""
import sys, numpy as np
sys.path.insert(0, '/home/gustav/Git/affine_toda/textbook_code')
from boundary_common import q_phys, sl2_spin, rcheck, cop
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
    M = complex(pref)*np.array([[complex(P(1)), complex(Q0)], [complex(Q0), complex(P(-1))]])
    return M if k % 2 == 0 else Xm @ M @ Xm

def ghoshal_breather(n, u, x, vt, lam):
    L = mp.mpf(lam)
    pre = ((-1)**(n+1)*mp.cos(u/2 + n*PI/(4*L))*mp.cos(u/2 - PI/4 - n*PI/(4*L))*mp.sin(u/2 + PI/4)
           / (mp.cos(u/2 - n*PI/(4*L))*mp.cos(u/2 + PI/4 + n*PI/(4*L))*mp.sin(u/2 - PI/4)))
    for l in range(1, n):
        pre *= (mp.sin(u + l*PI/(2*L))*mp.cos(u/2 - PI/4 - l*PI/(4*L))**2 / (mp.sin(u - l*PI/(2*L))*mp.cos(u/2 + PI/4 + l*PI/(4*L))**2))
    def Sx(y):
        if n == 1: return (mp.cos(y/L) - mp.sin(u))/(mp.cos(y/L) + mp.sin(u))
        c1, c2 = mp.cos(y/L - PI/(2*L)), mp.cos(y/L + PI/(2*L))
        return (mp.sin(u) - c1)*(mp.sin(u) - c2)/((mp.sin(u) + c1)*(mp.sin(u) + c2))
    return pre*Sx(x)*Sx(1j*vt)

def inv_vec(q, ya, yb):
    A, B = sl2_spin(1, q, ya), sl2_spin(1, q, yb)
    M = np.vstack([cop(A, B, kd, j) - (np.eye(4) if kd == 'k' else 0) for j in (0, 1) for kd in 'efk'])
    _, s, Vh = np.linalg.svd(M); w = Vh[-1].conj()
    return s[-1]/s[0], w/w[1]

devs = []; ncase = 0
for lam, eta, vt in ((3.3, 5.5, 0.3), (3.7, 4.1, 0.0)):
    q = q_phys(lam); Rsg = lambda z: rcheck(sl2_spin(1, q, z), sl2_spin(1, q, 1.0))[0]
    L = mp.mpf(lam); etab = etabar(eta, lam)
    for m in (1, 2):
        um = np.pi - m*np.pi/lam; th = 0.3
        ta, tb = th - 1j*um/2, th + 1j*um/2
        r_in, w = inv_vec(q, np.exp(lam*ta), np.exp(lam*tb)); assert r_in < 1e-10
        S2 = complex(a_S0(-1j*(ta + tb), lam))*Rsg(np.exp(lam*(ta + tb)))
        uth = -1j*mp.mpf(th)
        for S in [(), (0,), (2,), (4,), (5,), (1, 2), (0, 3), (3, 1), (1, 2, 0), (2, 0, 3)]:
            Ka, Kb = Kstate(S, ta, eta, vt, lam), Kstate(S, tb, eta, vt, lam)
            v = np.kron(np.eye(2), Ka) @ S2 @ np.kron(np.eye(2), Kb) @ w
            RB = v[1]
            xS = eta if len(S) % 2 == 0 else etab
            pred = ghoshal_breather(m, uth, xS, vt, lam)
            for i, n in enumerate(S, start=1):
                xx = eta if i % 2 == 1 else etab
                pred *= a_n(n, xx, uth + PI/2 - m*PI/(2*L), lam)*a_n(n, xx, uth - PI/2 + m*PI/(2*L), lam)
            d = max(np.linalg.norm(v - RB*w)/np.linalg.norm(v), abs(RB/complex(pred) - 1))
            devs.append(d); ncase += 1
            print(lam, m, S, f'{d:.1e}', flush=True)
report(f'B_1, B_2 on excited states up to depth 3, labels up to 5 ({ncase} cases): full matrix bootstrap gives K_(B_m)(x_S) prod b_(m n_i) and decouples the soliton multiplet [{max(devs):.0e}]',
       max(devs) < 1e-9)
