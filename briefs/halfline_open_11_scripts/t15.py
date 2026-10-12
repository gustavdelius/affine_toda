import numpy as np
from scal import *
from gz import logR0, logsigma
for (n, l) in [(1, 1), (2, 1), (2, 0), (3, 1), (3, 2), (3, 0)]:
    om, eta = 0.37, 0.3+0.4j
    T, eps, epsb, lam = params(n, om, l, eta)
    Ub = lambda m: lam**2*U(m, eta + 1j*np.pi*(epsb-eps))
    out = []
    for mu in (0.11+0.05j, -0.2+0.13j, 0.3-0.1j):
        a, b = amps(mu, n, om, l, eta); am, bm = amps(-mu, n, om, l, eta)
        u1 = a*am*U(mu, eta); u2 = b*bm*Ub(mu)
        a1, b1 = amps(T/2-mu, n, om, l, eta); a2, b2 = amps(T/2+mu, n, om, l, eta)
        c1 = a1/(Y(mu, n, om)*b2); c2 = b1/(Y(mu, n, om)*a2)
        line = [u1, u2, c1, c2]
        if n == 1:
            xi = np.pi/2 - np.pi*eps; u = 1j*np.pi*mu/T  # u = -i theta, theta = i pi mu/T
            gzv = np.exp(logR0(u, om) + logsigma(xi, u, om))/np.cos(xi)
            line.append(a/gzv)
        out.append(np.round(line, 6))
    print(n, l, out)
