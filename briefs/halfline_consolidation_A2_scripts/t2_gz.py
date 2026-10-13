"""Task 2, n = 1: SP minimal solution against Ghoshal-Zamolodchikov Dirichlet, with the book's R_0, sigma (ch. 29,
eq-bsg-k-full).  lambda = omega, u = -i theta = pi mu/omega.  Run: python3 -u t2_gz.py > t2_gz.log"""
import mpmath as mp
from sp_mp import SP, logprod, blk
mp.mp.dps = 30


def R0(u, lam):
    P = 4*lam
    t = lambda s: [(1, -2*lam*s*u/mp.pi), (1, -4*lam + 1 - 2*lam*s*u/mp.pi),
                   (-1, -3*lam - 2*lam*s*u/mp.pi), (-1, -lam + 1 - 2*lam*s*u/mp.pi)]
    return mp.exp(logprod(P, t(1) + [(-c, a) for c, a in t(-1)]))


def sigma(x, u, lam):
    P = 2*lam
    h = mp.mpf(1)/2
    t = lambda s: [(1, h + x/mp.pi - lam - lam*s*u/mp.pi), (1, h - x/mp.pi - lam - lam*s*u/mp.pi),
                   (-1, h - x/mp.pi - 2*lam - lam*s*u/mp.pi), (-1, h + x/mp.pi - lam*s*u/mp.pi)]
    return mp.cos(x)/mp.cos(x + lam*u)*mp.exp(logprod(P, t(1) + [(-c, a) for c, a in t(-1)]))


lam = mp.mpf('0.37')
# (i) Dirichlet limit of eq-bsg-k-full: sigma(i vartheta, u) * P0^pm / cosh(vartheta) -> cos(eta_B +- lam u) ?
print('(i) sigma(i*vartheta, u) as vartheta grows (u = 0.3+0.2i):')
u = mp.mpc('0.3', '0.2')
for vt in (5, 10, 20, 40):
    print('   vartheta', vt, mp.nstr(sigma(1j*vt, u, lam), 15))

# (ii) unitarity and crossing of the GZ Dirichlet scalar R0 sigma(xi,u)/cos xi as printed in the book
xi = mp.mpc('0.7', '0.3')
for u in (mp.mpc('0.3', '0.2'), mp.mpc('-0.5', '0.1')):
    g = lambda uu: R0(uu, lam)*sigma(xi, uu, lam)/mp.cos(xi)
    print('(ii) R0 unitarity', mp.nstr(R0(u, lam)*R0(-u, lam) - 1, 3),
          ' sigma(x,u)sigma(x,-u) identity', mp.nstr(sigma(xi, u, lam)*sigma(xi, -u, lam)
                                                    * mp.cos(xi + lam*u)*mp.cos(xi - lam*u)/mp.cos(xi)**2 - 1, 3),
          ' sigma even in x', mp.nstr(sigma(xi, u, lam)/sigma(-xi, u, lam) - 1, 3))

# (iii) a_min K_V against GZ Dirichlet, n = 1, l = 1, xi = pi eps - pi/2 (from eps_1/eps_0 = -e^{2 eta})
print('(iii) n = 1, l = 1: P^+/(a K_V00) and P^-/(a K_V11) divided by (-T/2)_I^2')
for om in ('0.37', '0.83'):
    lam = mp.mpf(om)
    for eta in ('0.41+0.23j', '-0.3+1.1j', '0.05+1.4j'):
        S = SP(1, lam, 1, mp.mpmathify(eta))
        xi = mp.pi*S.eps - mp.pi/2
        worst = []
        for mu in (mp.mpc('0.11', '0.19'), mp.mpc('-0.23', '0.07'), mp.mpc('0.2', '-0.15')):
            u = mp.pi*mu/lam
            th = 1j*u
            gz = R0(u, lam)*sigma(xi, u, lam)/mp.cos(xi)
            KV00 = mp.sinh(S.eta + S.T*th)
            KV11 = mp.sinh(S.eta - S.T*th)
            a = S.a(mu)
            sA = blk(-S.T/2, mu, S.T)**2
            worst.append((gz*mp.cos(xi + lam*u)/(a*KV00)/sA, gz*mp.cos(xi - lam*u)/(a*KV11)/sA))
        print(f'   om={om} eta={eta}: ' + ', '.join(f'({mp.nstr(p, 12)}, {mp.nstr(m_, 12)})' for p, m_ in worst))
