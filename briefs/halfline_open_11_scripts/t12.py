import numpy as np
from scal import *
for (n, l) in [(1, 1), (2, 1), (2, 0), (3, 1), (3, 2)]:
    om, eta = 0.37, 0.3+0.4j
    T, eps, epsb, lam = params(n, om, l, eta)
    la = lambda m: logRt(m, n, om) + logg_a(m, n, om, l, eta)
    lb = lambda m: logRt(m, n, om) + logg_b(m, n, om, l, eta)
    Ub = lambda m: lam**2*U(m, eta + 1j*np.pi*(epsb-eps))
    print(n, l)
    for mu in (0.11+0.05j, -0.2+0.13j, 0.3-0.1j):
        pa = np.exp(la(mu)+la(-mu))*U(mu, eta)
        pb = np.exp(lb(mu)+lb(-mu))*Ub(mu)
        q1 = np.exp(la(T/2-mu) - lb(T/2+mu))/Y(mu, n, om)
        q2 = np.exp(lb(T/2-mu) - la(T/2+mu))/Y(mu, n, om)
        print('   mu', mu, 'Pa', np.round(pa, 6), 'Pb', np.round(pb, 6), 'Q1', np.round(q1, 6), 'Q2', np.round(q2, 6))
