import numpy as np
from scal import *
for (n, l) in [(1, 1), (2, 1), (2, 0), (3, 1), (3, 2)]:
    om, eta = 0.37, 0.3+0.4j
    T, eps, epsb, lam = params(n, om, l, eta)
    ca = 1j; cb = 1j/lam; ratio = (cb/ca)**2
    print(n, l)
    for mu in (0.11+0.05j, -0.2+0.13j, 0.3-0.1j, 0.05-0.4j):
        qR = np.exp(logR(T/2-mu, n, om) - logR(T/2+mu, n, om))/Y(mu, n, om)
        # gamma crossing residual: gamma(T/2-mu)/gammabar(T/2+mu) vs (cb/ca) sin pi(eps+T/2-mu)/sin pi(epsb+T/2+mu)
        Z = (cb/ca)*np.sin(np.pi*(eps+T/2-mu))/np.sin(np.pi*(epsb+T/2+mu))
        qg = np.exp(logga(T/2-mu, n, om, l, eta, ratio) - loggb(T/2+mu, n, om, l, eta, ratio))/Z
        print('   mu', mu, 'QR', np.round(qR, 6), ' Qg', np.round(qg, 6))
