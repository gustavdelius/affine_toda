import numpy as np, sys
from breather import KB
from t18 import lsfit
om = 0.37; l = 1; eta = 0.3+0.4j
n = int(sys.argv[1])
T = (n+1)*om/2; x = (T-1)/2
eps = eta/(1j*np.pi); A = eps + x; Ab = A + om*l - T
mus = [0.03*k - 0.33 + 0.045j*((k*k) % 7 - 3) for k in range(1, 23)]
vals = []
for m in mus:
    th = 1j*np.pi*m/T
    E = (np.sinh(th) - 1j*np.sin(np.pi*A/T))/(np.sinh(th) + 1j*np.sin(np.pi*Ab/T))
    vals.append(KB(m, n, om, l, eta)/E)
ws = [np.exp(1j*np.pi*m/T) for m in mus]
g = lambda w: np.round(T*np.log(w+0j)/(1j*np.pi), 4)
for d in range(0, 10):
    num, den, res = lsfit(vals, ws, d)
    if res < 1e-5:
        print(n, 'T', T, 'deg', d, 'res %.1e' % res, 'zeros', sorted([g(z) for z in np.roots(num[::-1])], key=lambda z: z.real),
              'poles', sorted([g(p) for p in np.roots(den[::-1])], key=lambda z: z.real), 'lead', np.round(num[-1]/den[-1], 5))
        break
else:
    print(n, 'no fit; last res', res)
