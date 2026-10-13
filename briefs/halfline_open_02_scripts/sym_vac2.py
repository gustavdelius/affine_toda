"""Uniform (C = +1) d_n vacua in the automorphism-invariant sector, closed form / quadrature profiles, and the
linearized reflection about them compared with the B -> 0 limit of the eps = 0 amplitude K_1.
d_4: phi = f (e2 - e3), u = e^f, u' = (u-1) sqrt(2u(u+2)), u(0) = 9/8 (Tzitzeica).
d_5: phi = f (e2 - e4), tanh(f/4) = e^{2(x - x0)}, e^{f(0)} = 4 - 2 sqrt 2 (sinh-Gordon)."""
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from dn_classical import droots, masses, linear_K


class Prof:
    def __init__(self, n, X=16.0):
        E = np.eye(n); self.n = n
        if n == 5:
            self.v = E[1] - E[3]
            f0 = np.log(4 - 2*np.sqrt(2)); x0 = -0.5*np.log(np.tanh(f0/4))
            self.f = lambda x: 4*np.arctanh(np.exp(2*(x - x0)))
        elif n == 4:
            self.v = E[1] - E[2]
            u0 = 9/8
            ds = np.linspace(np.log(u0 - 1), np.log(u0 - 1) - 40, 400)       # log(u - 1) grid
            xs = []
            for d in ds:
                u = 1 + np.exp(d)
                # x(u) = - int_u^{u0} du/((u-1) sqrt(2u(u+2))) ; in variable s = log(u-1): du/(u-1) = ds
                xs.append(-quad(lambda s: 1/np.sqrt(2*(1 + np.exp(s))*(3 + np.exp(s))), d, np.log(u0 - 1),
                                epsabs=1e-14, epsrel=1e-13)[0])
            xs = np.array(xs)[::-1]; ds = ds[::-1]
            sp = CubicSpline(xs, ds)
            self.f = lambda x: np.log(1 + np.exp(sp(x)))
        self.t = np.array([-X, 0.0])

    def sol(self, x):
        return np.concatenate([self.f(x)*self.v, np.zeros(self.n)])


def bl(x, th, h):
    return np.sinh(th/2 + 1j*np.pi*x/(2*h))/np.sinh(th/2 - 1j*np.pi*x/(2*h))


def K1cl(th, n):
    h = 2*n - 2
    return -1/(bl(1, th, h)*bl(h - 1, th, h))*bl(h/2 - 1, th, h)*bl(h/2 + 1, th, h)/bl(h/2, th, h)**2


if __name__ == '__main__':
    import sys
    np.set_printoptions(linewidth=180, precision=8, suppress=True)
    for n in [int(a) for a in sys.argv[1:]] or (5, 4):
        m, U = masses(n)
        P = Prof(n)
        al, nm = droots(n)
        phi0 = P.sol(0.0)[:n]
        bc = 0  # BC check with numerical derivative
        dphi = (P.sol(0.0)[:n] - P.sol(-1e-6)[:n])/1e-6
        print(n, 'BC residual', np.abs(dphi + (al.T*np.sqrt(nm)) @ np.exp(al @ phi0/2)).max(), flush=True)
        for th in (0.6, 0.3):
            Om = m[0]*np.cosh(th)
            K, op = linear_K(n, np.ones(n + 1), P, Om)
            print(n, 'theta', th, 'open', op, '\n K =', K, '\n eig', np.linalg.eigvals(K), ' K1cl =', K1cl(th, n),
                  ' a_n-type -1/((1)(h-1)) =', -1/(bl(1, th, 2*n - 2)*bl(2*n - 3, th, 2*n - 2)), flush=True)
