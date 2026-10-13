"""Vacuum checks (mp precision) for d_4, d_5 C = +1 invariant kinks: field equation, BC, regularity, energy."""
import mpmath as mp
import numpy as np
from cl_lin import Kink, droots

mp.mp.dps = 30
for n in (4, 5):
    K = Kink(n)
    A, nj = droots(n)
    Am = mp.matrix(A.tolist()); v = mp.matrix(K.v.tolist())
    c = [sum(Am[j, k]*v[k] for k in range(n)) for j in range(n + 1)]
    f = K.f_mp
    # field equation f'' v = sum n_j a_j e^{c_j f}: project on v (others vanish identically: check)
    worst = 0
    for x in (mp.mpf('-0.3'), mp.mpf('-2.1'), mp.mpf('-5'), mp.mpf(0)):
        d2 = mp.diff(f, x, 2)
        rhs = [sum(nj[j]*Am[j, k]*mp.exp(c[j]*f(x)) for j in range(n + 1)) for k in range(n)]
        worst = max(worst, max(abs(d2*v[k] - rhs[k]) for k in range(n)))
    # BC
    x0 = mp.mpf(0)
    d1 = mp.diff(f, x0)
    bc = [d1*v[k] + sum(mp.sqrt(nj[j])*Am[j, k]*mp.exp(c[j]*f(x0)/2) for j in range(n + 1)) for k in range(n)]
    # energy
    def dens(x):
        fx = f(x)
        return mp.diff(f, x)**2*(v.T*v)[0]/2 + sum(nj[j]*(mp.exp(c[j]*fx) - 1) for j in range(n + 1))
    bulk = mp.quad(dens, [-mp.inf, -5, -1, 0])
    bnd = 2*sum(mp.sqrt(nj[j])*(mp.exp(c[j]*f(x0)/2) - 1) for j in range(n + 1))
    # regularity: f finite and real on a grid
    xs = np.linspace(-20, 0, 2001)
    fx = K.f(xs)
    print(f'd_{n}: c_j = {[int(round(float(t))) for t in c]}, e^f(0) = {mp.nstr(mp.exp(f(x0)), 20)}, '
          f'field eq residual {mp.nstr(worst, 3)}, BC residual {mp.nstr(max(abs(t) for t in bc), 3)}, '
          f'f regular on [-20,0]: {np.all(np.isfinite(fx))}, min f {fx.min():.2e} max f {fx.max():.6f}')
    print(f'   energy (m/beta^2): bulk {mp.nstr(bulk, 15)}, boundary {mp.nstr(bnd, 15)}, total {mp.nstr(bulk + bnd, 15)}')
