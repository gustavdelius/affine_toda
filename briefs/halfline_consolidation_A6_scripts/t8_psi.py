"""Fit the crossing function of the d_n K with s = +-1 to closed forms."""
import sys
import numpy as np
from t7_cross import Kmat, Piso, n, N
from qk import dn_vector, q_phys, chevalley, rcheck_w

h = 2*n - 2


def psi_num(name, om, th, Pi):
    q = q_phys(om)
    km = (Pi @ Kmat(name, om, 1j*np.pi/2 - th)).T.reshape(N*N)
    kp = (Pi @ Kmat(name, om, 1j*np.pi/2 + th)).T.reshape(N*N)
    R = rcheck_w(chevalley(dn_vector(n, q, 1.0), np.exp(2*om*th)), chevalley(dn_vector(n, q, 1.0), 1.0))
    v = R @ kp
    return np.vdot(v, km)/np.vdot(v, v)


S = lambda z: np.sin(np.pi*z)
Cc = lambda z: np.cos(np.pi*z)
for om in (0.37, 0.61, 1.29):
    T = (n - 1)*om
    P, _ = Piso(om)
    Pi = np.linalg.inv(P)
    w = np.exp(-0.5j*np.pi*om)
    for name, s in (('s+', 1), ('s-', -1)):
        sig = s*(-1)**n
        for th in (0.31, 0.83 + 0.1j, -0.4 + 0.25j):
            ps = psi_num(name, om, th, Pi)
            mu = -1j*T*th/np.pi
            # normalization n(theta) = (y^h - sig)(w y)^(h/2), y = e^{om theta}
            nf = lambda t: (np.exp(h*om*t) - sig)*(w*np.exp(om*t))**(h/2)
            ps_hat = ps*nf(1j*np.pi/2 + th)/nf(1j*np.pi/2 - th)
            cand = {
                'a_n Psi_+': S(om/2 - mu)/S(om/2 + mu),
                'a_n Psi_-': Cc(mu - om/2)/Cc(mu + om/2),
            }
            print(f'om={om} {name} sig={sig} th={th}: Psi_hat = {ps_hat:.10f}; ' +
                  '; '.join(f'{k}: {v:.10f}' for k, v in cand.items()))
