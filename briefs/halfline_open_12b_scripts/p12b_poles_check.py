"""Part (b): checks that the exact pole-order counter of p12b_poles.py (Gamma-function arguments of R_0, sigma and sines of the a_n blocks)
agrees with numerically measured orders of  Q^S(u) = Q(x_S, u) prod a_ni  near every nu_j, w_j in the physical strip
(order = log2 |f(u0+d)/f(u0+2d)|, d = 1e-6; 20-digit arithmetic)."""
import math
from p12b_common import *
from p12b_poles import total_events
def report(name, ok): print(('PASS ' if ok else 'FAIL ') + name, flush=True)
bad = 0; n = 0
for lam, eta, states in ((6.3, 5.0, [(), (0,), (1,), (0, 5), (0, 5, 1)]),
                         (7.0, 0.9*math.pi*8/2, [(0,), (1,), (0, 1), (0, 1, 1), (2, 3)])):
    Lm = mp.mpf(lam); em = mp.mpf(eta); eb = etabar(em, Lm)
    for S in states:
        def f(u):
            k = len(S); r = Qfun(em if k % 2 == 0 else eb, u, Lm)
            for i, nn in enumerate(S, 1): r *= a_n(nn, em if i % 2 == 1 else eb, u, Lm)
            return r
        ev = total_events(S, eta, lam)
        pos = sorted({round(float(nu(j, em, Lm)), 9) for j in range(-6, 8)} | {round(float(nu(j, eb, Lm)), 9) for j in range(-2, 12)})
        for p in pos:
            if not (0 < p < math.pi/2): continue
            d = mp.mpf('1e-6')
            num = float(mp.log(abs(f(mp.mpf(p) + d)/f(mp.mpf(p) + 2*d)))/mp.log(2))
            ex = ev.get(round(p, 9), 0.0)
            # candidate positions may be off by 1e-9 rounding; recompute exactly
            n += 1
            if abs(num - ex) > 0.02: bad += 1; print('MISMATCH', lam, S, p, num, ex, flush=True)
report(f'exact pole-order counter agrees with numerically measured orders at {n} (state, position) pairs [{bad} mismatches]', bad == 0)
