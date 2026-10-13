"""Gandenberger-type closed forms of the s = +-1 scalar factors for all n (check against A_d rho sigma, -A_d rho_-).
A^(+) = a0p a1 (T/2 - om/2)_I (T/2 + om/2)_I,   A^(-) = -i (-T/2)_I a0m a1,  a1 = DG (3.9) without prefactor,
a0p = sin(pi om/2)/sin pi(mu - T/2) prod_{j>=1} G_j(mu, om/2-3T/2) G_j(mu, 1-om/2-T/2) / [G_j(mu, om/2-T/2) G_j(mu, 1-om/2-3T/2)]
a0m = sin(pi om/2)/cos pi(mu - T/2) prod_{j>=1} G_j(mu, om/2-3T/2+1/2) G_j(mu, 1/2-om/2-T/2) / [G_j(mu, om/2-T/2+1/2) G_j(mu, 1/2-om/2-3T/2)]
G_j(mu, a) = Gamma(mu + j h om + a)/Gamma(-mu + j h om + a), T = h om/2."""
import numpy as np
from scipy.special import loggamma as lg
from rhonp import rho_plus, blk
from rhom import rho_minus
from gz import rich, A1
from boundary_common import report
def a0(mu, n, om, s, K=4000):
    T = (n + 1)*om/2; sh = 0 if s == 1 else 0.5; one = 1 if s == 1 else 0.5
    G = lambda j, a: lg(mu + 2*j*T + a) - lg(-mu + 2*j*T + a)
    part = lambda KK: np.sum(G(np.arange(1, KK + 1), om/2 - 1.5*T + sh) + G(np.arange(1, KK + 1), one - om/2 - T/2)
                             - G(np.arange(1, KK + 1), om/2 - T/2 + sh) - G(np.arange(1, KK + 1), one - om/2 - 1.5*T))
    pre = np.sin(np.pi*om/2)/(np.sin(np.pi*(mu - T/2)) if s == 1 else np.cos(np.pi*(mu - T/2)))
    return pre*np.exp(rich(part, K))
w = 0
for om in (0.37, 0.9):
    for n in (1, 2, 3, 4, 5):
        T = (n + 1)*om/2
        for mu in (0.13 + 0.21j, -0.3 + 0.1j):
            a1 = A1(mu, n, om)/blk(-T/2, mu, T)
            Ap = A1(mu, n, om)*rho_plus(mu, n, om)*blk(T/2, mu, T)**2*blk(T/2 - om/2, mu, T)*blk(T/2 + om/2, mu, T)
            Am = -A1(mu, n, om)*rho_minus(mu, n, om)
            w = max(w, abs(Ap/(a0(mu, n, om, 1)*a1*blk(T/2 - om/2, mu, T)*blk(T/2 + om/2, mu, T)) - 1),
                    abs(Am/(-1j*blk(-T/2, mu, T)*a0(mu, n, om, -1)*a1) - 1))
report(f'closed forms A^(+) = a0p a1 (T/2-om/2)(T/2+om/2), A^(-) = -i(-T/2) a0m a1, n = 1..5 [{w:.1e}]', w < 1e-6)
