"""Fit the singlet coefficient c_sigma(mu) of d_n to a product of sines (least squares on log|c|)."""
import sys
import itertools
import numpy as np
from t9_proj import Khat, singlet, n, N
from qk import dn_vector, q_phys, rcheck_w, kron

h = 2*n - 2


def cval(name, om, mu, sig):
    T = (n - 1)*om
    q = q_phys(om)
    thb = 1j*np.pi*mu/T
    u = np.pi - np.pi/T
    ta, tb = thb - 0.5j*u, thb + 0.5j*u
    Va, Vb = dn_vector(n, q, np.exp(om*ta)), dn_vector(n, q, np.exp(om*tb))
    w_in = singlet(Va, Vb)[0]
    Wa = dn_vector(n, q, np.exp(-om*tb))
    w_out = singlet(Wa, dn_vector(n, q, np.exp(-om*ta)))[0]
    R = rcheck_w(Va, Wa)
    v = kron(np.eye(N), Khat(name, om, ta, sig)) @ (R @ (kron(np.eye(N), Khat(name, om, tb, sig)) @ w_in))
    return np.vdot(w_out, v)/np.vdot(w_out, w_out)


if __name__ == '__main__':
    rng = np.random.default_rng(3)
    for om in (0.53, 0.71):
        T = (n - 1)*om
        shifts = {'0': 0, 'h': 0.5, 'o/2': om/2, '-o/2': -om/2, 'h+o/2': 0.5 + om/2, 'h-o/2': 0.5 - om/2,
                  'T': T, 'T/2': T/2, 'h+T/2': 0.5 + T/2, 'T+o/2': T + om/2, 'T-o/2': T - om/2,
                  'h+T': 0.5 + T, 'T/2+o/2': T/2 + om/2, 'T/2-o/2': T/2 - om/2, 'h+T/2+o/2': 0.5 + T/2 + om/2,
                  'h+T/2-o/2': 0.5 + T/2 - om/2, 'o': om, 'h+o': 0.5 + om}
        keys = list(shifts)
        mus = [complex(a, b) for a, b in rng.uniform(-0.45, 0.45, size=(40, 2))]
        for name, s in (('s+', 1), ('s-', -1)):
            sig = s*(-1)**n
            cs = np.array([cval(name, om, m, sig) for m in mus])
            Aux = np.array([[1.0] + [np.log(abs(np.sin(np.pi*(m + shifts[k])))) for k in keys] for m in mus])
            sol, res, rk, sv = np.linalg.lstsq(Aux, np.log(np.abs(cs)), rcond=None)
            print(f'om={om} {name}: rank {rk}/{Aux.shape[1]}, residual {res}')
            print('   ', {k: round(v, 3) for k, v in zip(['const'] + keys, sol) if abs(v) > 1e-6})
