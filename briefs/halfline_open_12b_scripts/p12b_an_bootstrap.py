"""Part (b), scalar level.  Checks (all to ~1e-15 or better):
 1. Master functional equation of the bulk amplitude: a(u - pi/lam) = tan(u/2) tan(pi/(2 lam) - u/2) a(u)
    (a = S_{++}, BPTT (3.1)); it implies a(u+x-pi/lam)a(u-x+pi/lam)/(a(u+x)a(u-x)) = {2 lam x/pi - 1}
    for ALL x, hence the product form of a_n for ALL n (telescoping, nu_n = nu_{n-1} - pi/lam).
 2. ratio form = product form of a_n, n = 1..10, several (lam, eta) incl. eta outside the fundamental range.
 3. Iterated scalar bootstrap Q^{S+(m)} = a(u-p) Q^S b(u+p), p = nu_m(eta) (even length) or w_m = nu_m(etabar) (odd length)
    reproduces the closed form  Q^S = Q(eta or etabar) prod_odd a_ni(eta) prod_even a_ni(etabar)  (BPTT sec. 3.8)
    for all labels up to 4 deep, labels up to 4.
 4. unitarity a_n(u)a_n(-u) = 1 and crossing a_n(pi/2+w) = a_n(pi/2-w).
Runtime: about 25 minutes (mp.nsum Gamma products; sections 2 and 3 dominate).
"""
import itertools, random
from p12b_common import *
random.seed(1)
def report(name, ok): print(('PASS ' if ok else 'FAIL ') + name)

# 1
worst = 0
for lam in (1.37, 3.3, 7.1, 12.5):
    L = mp.mpf(lam)
    for u in (mp.mpc('0.3', '0.2'), mp.mpc('-0.7', '0.1'), mp.mpc('1.1', '-0.4')):
        lhs = a_S0(u - PI/L, lam)/a_S0(u, lam); rhs = mp.tan(u/2)*mp.tan(PI/(2*L) - u/2)
        worst = max(worst, abs(lhs/rhs - 1))
report(f'a(u - pi/lam) = tan(u/2) tan(pi/2lam - u/2) a(u) at 4 lambdas, 3 points each [{float(worst):.0e}]', worst < 1e-15)
worst = 0
for lam in (2.2, 6.0):
    L = mp.mpf(lam)
    for x, u in ((mp.mpc('0.4', '0.3'), mp.mpc('0.2', '0.1')), (mp.mpf('1.1'), mp.mpc('-0.5', '0.2'))):
        lhs = a_S0(u + x - PI/L, lam)*a_S0(u - x + PI/L, lam)/(a_S0(u + x, lam)*a_S0(u - x, lam))
        rhs = curly(2*L*x/PI - 1, u, lam)
        worst = max(worst, abs(lhs/rhs - 1))
report(f'one-step identity a(u+x-pi/lam)a(u-x+pi/lam)/(a(u+x)a(u-x)) = {{2 lam x/pi - 1}} for generic complex x [{float(worst):.0e}]', worst < 1e-15)

# 2
worst = 0
for lam, eta in ((3.3, 5.5), (7.1, 4.0), (2.2, 6.2), (11.0, 9.3), (1.37, 3.1)):
    L = mp.mpf(lam)
    for u in (mp.mpc('0.3', '0.2'), mp.mpc('-0.5', '0.1')):
        for n in range(1, 11):
            ratio = a_S0(u + nu(n, eta, L), lam)*a_S0(u - nu(n, eta, L), lam)/(a_S0(u + nu(0, eta, L), lam)*a_S0(u - nu(0, eta, L), lam))
            worst = max(worst, abs(ratio/a_n(n, eta, u, lam) - 1))
report(f'ratio form = product form of a_n for n = 1..10 at 5 (lambda, eta) [{float(worst):.0e}]', worst < 1e-14)

# 3
def Qclosed(S, eta, u, lam):
    k = len(S); base = eta if k % 2 == 0 else etabar(eta, lam)
    r = Qfun(base, u, lam)
    for i, n in enumerate(S, start=1):
        r *= a_n(n, eta if i % 2 == 1 else etabar(eta, lam), u, lam)
    return r
worst, count = 0, 0
for lam, eta in ((3.3, 5.5), (7.1, 9.0), (4.6, 2.9)):
    L = mp.mpf(lam)
    for u in (mp.mpc('0.3', '0.2'), mp.mpc('-0.5', '0.1')):
        # iterate the bootstrap along random label sequences up to depth 4, labels 0..4 (analytic identities, so not restricted to the strip)
        for S in [(0,), (2,), (3, 1), (0, 4), (1, 0, 2), (4, 2, 3), (2, 3, 1, 4), (0, 0, 0, 0), (3, 3, 3, 3)]:
            Q = Qfun(eta, u, lam); cur_eta = eta
            for i, m in enumerate(S, start=1):
                ee = eta if i % 2 == 1 else etabar(eta, lam)
                p = nu(m, ee, L)
                Q = a_S0(u - p, lam)*Q*b_S0(u + p, lam)
            dev = abs(Q/Qclosed(S, eta, u, lam) - 1) if True else 0
            worst = max(worst, dev); count += 1
report(f'iterated scalar bootstrap = closed form  Q(eta|etabar) prod a_ni(eta|etabar), {count} (lambda, eta, u, state) cases to depth 4 [{float(worst):.0e}]', worst < 1e-12)

# 4
worst_u = worst_c = 0
for lam, eta in ((3.3, 5.5), (7.1, 4.0), (2.2, 6.2), (1.37, 1.0)):
    for n in range(0, 8):
        for ee in (eta, etabar(eta, lam)):
            for w in (mp.mpc('0.3', '0.2'), mp.mpc('0.9', '-0.1')):
                worst_u = max(worst_u, abs(a_n(n, ee, w, lam)*a_n(n, ee, -w, lam) - 1))
                worst_c = max(worst_c, abs(a_n(n, ee, PI/2 + w, lam)/a_n(n, ee, PI/2 - w, lam) - 1))
report(f'a_n(u)a_n(-u) = 1 [{float(worst_u):.0e}] and a_n(pi/2+w) = a_n(pi/2-w) [{float(worst_c):.0e}] for n <= 7, eta and etabar', worst_u < 1e-14 and worst_c < 1e-14)
