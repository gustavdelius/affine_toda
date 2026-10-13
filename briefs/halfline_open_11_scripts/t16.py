import numpy as np
from scal import *
from gz import logR0, logsigma
from t14 import ratfit
n, l, om = 1, 1, 0.37
T = om
for eta in (0.3+0.4j, -0.2+0.9j):
    T, eps, epsb, lam = params(n, om, l, eta)
    xi = np.pi/2 - np.pi*eps
    mus = [0.05*k - 0.3 + 0.07j*((k*k) % 5 - 2) for k in range(1, 14)]
    rat = []
    for mu in mus:
        a, b = amps(mu, n, om, l, eta)
        u = np.pi*mu/T
        gzv = np.exp(logR0(u, om) + logsigma(xi, u, om))/np.cos(xi)
        rat.append(a/gzv)
    print(eta, np.round(rat[:4], 6))
    ws = [np.exp(1j*np.pi*m/T) for m in mus]
    for d in range(0, 5):
        num, den = ratfit(rat[:2*d+2], ws[:2*d+2], d)
        res = max(abs(np.polyval(num[::-1], w)/np.polyval(den[::-1], w) - q) for w, q in zip(ws, rat))
        if res < 1e-6:
            f = lambda w: np.round(T*np.log(w+0j)/(1j*np.pi), 5)
            print('  deg', d, 'zeros', [f(z) for z in np.roots(num[::-1])], 'poles', [f(p) for p in np.roots(den[::-1])], 'lead', num[-1]/den[-1] if d else num[0])
            break
