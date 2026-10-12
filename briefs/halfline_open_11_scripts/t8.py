import numpy as np
from sp_core import *
from boundary_common import q_phys
om = 0.37; q = q_phys(om)
# Rbar transmission entries
for n in (2, 3):
    N = n+1
    Rb = R(Vy(n, q, 0.3+0.1j, om, True), Vy(n, q, -0.2, om))
    vals = [Rb[j*N+i, i*N+j] for i in range(N) for j in range(N) if i != j]
    print(n, 'Rbar V transmission entries', np.round(vals, 6))
# fit Phi1 as rational in x=exp(2 T th)
def ratfit(f, xs, dn, dd):
    # f = P(x)/Q(x), Q monic leading? use Q(0)=1 normalisation: P - f*Q = 0
    M = []; b = []
    for x, y in zip(xs, f):
        row = [x**k for k in range(dn+1)] + [-y*x**k for k in range(1, dd+1)]
        M.append(row); b.append(y)
    sol, *_ = np.linalg.lstsq(np.array(M), np.array(b), rcond=None)
    return sol[:dn+1], np.concatenate([[1], sol[dn+1:]])
for n in (1, 2, 3):
    N = n+1; T = N*om/2
    for l in range(0, N):
        eta = 0.3+0.4j
        ths = [0.05*k + 0.01 + 0.03j*((k*k) % 7) for k in range(1, 16)]
        P1 = []; P2 = []
        for th in ths:
            c1, c2, r = Phis(n, om, l, eta, th); P1.append(c1); P2.append(c2)
        xs = [np.exp(2*T*th) for th in ths]
        for name, P in (('Phi1', P1), ('Phi2', P2)):
            for d in (1, 2, 3):
                num, den = ratfit(P[:2*d+3], xs[:2*d+3], d, d)
                res = max(abs(np.polyval(num[::-1], x)/np.polyval(den[::-1], x) - p) for x, p in zip(xs, P))
                if res < 1e-8:
                    zn = np.roots(num[::-1]); zd = np.roots(den[::-1])
                    # express roots x0 = exp(2 eta') ... print log(x0)/2 relative to eta and i pi om
                    f = lambda z: np.round(np.log(z+0j)/2, 5)
                    print(n, l, name, 'deg', d, 'zeros log/2:', [f(z) for z in zn], 'poles log/2:', [f(z) for z in zd], 'lead', np.round(num[-1]/den[-1], 5))
                    break
            else:
                print(n, l, name, 'no fit')
