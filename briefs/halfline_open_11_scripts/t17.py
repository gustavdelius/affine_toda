import numpy as np
from breather import KB
from scal import blk, params
om = 0.37; n = 1; l = 1; T = om
def ghoshal(u, lam, eta):
    pre = (np.cos(u/2+np.pi/(4*lam))*np.cos(u/2-np.pi/4-np.pi/(4*lam))*np.sin(u/2+np.pi/4)/
           (np.cos(u/2-np.pi/(4*lam))*np.cos(u/2+np.pi/4+np.pi/(4*lam))*np.sin(u/2-np.pi/4)))
    return pre*(np.cos(eta/lam)-np.sin(u))/(np.cos(eta/lam)+np.sin(u))
cdd = lambda m: blk(-T/2, m, T)**2
for eta in (0.3+0.4j, -0.1+0.7j):
    eps = params(n, om, l, eta)[1]; xi = np.pi/2 - np.pi*eps
    for mu in (0.11+0.05j, -0.2+0.13j, 0.3-0.1j):
        u = np.pi*mu/T
        k1 = KB(mu, n, om, l, eta, cdd=cdd); k2 = KB(mu, n, om, l, eta, cdd=cdd, conj=True)
        print(eta, mu, 'KB/Ghoshal', np.round(k1/ghoshal(u, om, xi), 6), np.round(k2/ghoshal(u, om, xi), 6), ' without cdd', np.round(KB(mu, n, om, l, eta)/ghoshal(u, om, xi), 6))
