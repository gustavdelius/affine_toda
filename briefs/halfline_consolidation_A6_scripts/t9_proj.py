"""d_n: (i) unitarity factor of the normalized K-hat = K_poly/((w^2-1)(y^h - sig)(w y)^(h/2)), sig = s(-1)^n;
(ii) singlet projection of the lowest breather: (1 x K(theta_a)) Rcheck (1 x K(theta_b)) on the B_1 singlet of
V(y_a) x V(y_b), theta_{a,b} = theta -+ i u/2, u = pi - pi/T; for eps = 0 (K = J) and s = +-1 (K-hat)."""
import sys
import numpy as np
from t7_cross import Kmat, n, N
from qk import dn_vector, q_phys, rcheck_w, cop, nodes_of, kron

h = 2*n - 2


def Khat(name, om, th, sig):
    w = np.exp(-0.5j*np.pi*om); y = np.exp(om*th)
    K = Kmat(name, om, th)
    if name == '0':
        return K
    return K/((w**2 - 1)*(y**h - sig)*(w*y)**(h//2))


def singlet(A, B):
    d = A['dim']*B['dim']
    M = np.vstack([cop(A, B, kd, j) - (np.eye(d) if kd == 'k' else 0) for j in nodes_of(A) for kd in 'efk'])
    _, sv, Vh = np.linalg.svd(M)
    wv = Vh[-1].conj()
    return wv/wv[int(np.argmax(np.abs(wv) > 1e-8))], sv[-1]/sv[0], sv[-2]/sv[0]


if __name__ == '__main__':
    S_ = lambda z: np.sin(np.pi*z)
    C_ = lambda z: np.cos(np.pi*z)
    for om, mu in ((0.53, 0.11 + 0.19j), (0.61, -0.17 + 0.23j), (1.37, 0.05 - 0.12j)):
        T = (n - 1)*om
        q = q_phys(om)
        for name, s in (('s+', 1), ('s-', -1)):
            sig = s*(-1)**n
            th = 0.37 - 0.21j
            mth = -1j*T*th/np.pi
            U = (Khat(name, om, -th, sig) @ Khat(name, om, th, sig))
            Upred = sig*(np.cos(2*np.pi*mth) - sig*np.cos(np.pi*T))*(np.cos(2*np.pi*mth) - sig*np.cos(np.pi*om)) \
                / (2*np.sin(np.pi*om/2)**2*(np.cos(2*np.pi*mth) - sig))
            print(f'om={om} {name} sig={sig}: K(-th)K(th) scalar [{np.abs(U - U[0, 0]*np.eye(N)).max()/abs(U[0, 0]):.1e}], '
                  f'U/Upred - 1 = {abs(U[0, 0]/Upred - 1):.1e}')
        thb = 1j*np.pi*mu/T
        u = np.pi - np.pi/T
        ta, tb = thb - 0.5j*u, thb + 0.5j*u
        Va, Vb = dn_vector(n, q, np.exp(om*ta)), dn_vector(n, q, np.exp(om*tb))
        w_in, sv1, sv2 = singlet(Va, Vb)
        Wa = dn_vector(n, q, np.exp(-om*tb))
        w_out, sv1o, _ = singlet(Wa, dn_vector(n, q, np.exp(-om*ta)))
        R = rcheck_w(Va, Wa)
        print(f'  singlet in V x V at u = pi - pi/T: smallest sv {sv1:.1e} (next {sv2:.1e}); out {sv1o:.1e}')
        for name, s in (('0', 0), ('s+', 1), ('s-', -1)):
            sig = s*(-1)**n
            v = kron(np.eye(N), Khat(name, om, ta, sig)) @ (R @ (kron(np.eye(N), Khat(name, om, tb, sig)) @ w_in))
            c = np.vdot(w_out, v)/np.vdot(w_out, w_out)
            res = np.linalg.norm(v - c*w_out)/np.linalg.norm(v)
            # candidates
            cands = {}
            if name != '0':
                cands['c_+an'] = -C_(mu)*C_(mu + om/2)*C_(mu - T)/(S_(om/2)**2*C_(mu - om/2))
                cands['c_-an'] = S_(mu)*S_(mu + om/2)*S_(mu - T)/(S_(om/2)**2*S_(mu - om/2))
            print(f'  {name}: eigen-residual {res:.1e}, c = {c:.10f}; ' + ', '.join(f'{k}={v_:.10f} ratio {c/v_:.10f}'
                                                                              for k, v_ in cands.items()))
