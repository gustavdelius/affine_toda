"""Exact (symbolic-in-position) pole/zero bookkeeping of the soliton reflection factors on the excited boundary states.
For a real u in the open physical strip 0 < u < pi/2 we count, with multiplicity, the poles (+) and zeros (-) of
   Q^S(u) = R_0(u) sigma(x_S,u) Q_0(u) prod_odd a_ni(eta,u) prod_even a_ni(etabar,u),     x_S = eta (|S| even), etabar (|S| odd)
from the Gamma-function arguments of R_0 and sigma (BPT/BPTT (3.2)) and from the zeros of the sines in the blocks.
(sigma(i vt,u) has no pole or zero in the strip [BPTT]; P^pm differ from Q only by the factor P0/Q0, which is entire, so P can
only have extra ZEROS, never extra poles: see check at the bottom.)
"""
import math
from collections import defaultdict

TOL = 1e-9

def gamma_events(factors, lam, lmax=60):
    """factors: list of (alpha, s, e): Gamma(alpha + s*lam*u/pi)^e.  Return dict position->order contribution in strip."""
    ev = defaultdict(float)
    for alpha, s, e in factors:
        if abs(s) < 1e-14: continue
        for m in range(0, 4000):
            # alpha + s*lam*u/pi = -m  ->  u = pi(-m - alpha)/(s lam)
            u = math.pi*(-m - alpha)/(s*lam)
            if 0 < u < math.pi/2: ev[round(u, 9)] += e
            # as m grows u moves monotonically; stop once outside on the far side
            if (s > 0 and u < 0) or (s < 0 and u > math.pi/2): 
                if (s > 0 and u < 0) or (s < 0 and u > math.pi/2): break
    return ev

def sigma_factors(x, lam, lmax=60):
    a = x/math.pi; f = []
    for l in range(1, lmax):
        # E(u): numerators Gamma(1/2 +- a + (2l-1) lam - v), denominators Gamma(1/2 -+ a ...), v = lam u/pi ; /E(-u) flips exponent and s
        for s, sign in ((-1, 1), (+1, -1)):
            f += [(0.5 + a + (2*l-1)*lam, s, +sign), (0.5 - a + (2*l-1)*lam, s, +sign),
                  (0.5 - a + (2*l-2)*lam, s, -sign), (0.5 + a + 2*l*lam, s, -sign)]
    return f

def R0_factors(lam, lmax=60):
    f = []
    for l in range(1, lmax):
        for s, sign in ((-2, 1), (+2, -1)):
            # Gamma(4 l lam - 2v) etc, 2v = 2 lam u/pi -> s*lam*u/pi with s = -+2
            f += [(4*l*lam, s, +sign), (4*lam*(l-1) + 1, s, +sign), ((4*l-3)*lam, s, -sign), ((4*l-1)*lam + 1, s, -sign)]
    return f

def trig_events(lam, x, n_odd_args):
    """1/cos(x + lam u) in sigma:  poles at x + lam u = pi/2 + k pi."""
    ev = defaultdict(float)
    for k in range(-400, 400):
        u = (math.pi/2 + k*math.pi - x)/lam
        if 0 < u < math.pi/2: ev[round(u, 9)] += 1
    return ev

def block_events(xx, ee):
    """(x)^ee, (x) = sin(u/2 + pi x/2)/sin(u/2 - pi x/2): pole at u = pi x + 2 pi m, zero at u = -pi x + 2 pi m."""
    ev = defaultdict(float)
    for m in range(-6, 7):
        up = math.pi*xx + 2*math.pi*m; uz = -math.pi*xx + 2*math.pi*m
        if 0 < up < math.pi/2: ev[round(up, 9)] += ee
        if 0 < uz < math.pi/2: ev[round(uz, 9)] -= ee
    return ev

def an_events(n, eta, lam):
    ev = defaultdict(float)
    for l in range(1, n+1):
        y = 2*(eta/math.pi - l)
        A, B = (y+1)/(2*lam), (y-1)/(2*lam)
        for xx, ee in ((A, 1), (B, 1), (A-1, -1), (B+1, -1)):
            for k, v in block_events(xx, ee).items(): ev[k] += v
    return ev

def add(d, e):
    for k, v in e.items(): d[k] += v

def total_events(S, eta, lam, lmax=60):
    k = len(S); etab = math.pi*(lam+1) - eta
    x = eta if k % 2 == 0 else etab
    ev = defaultdict(float)
    add(ev, gamma_events(R0_factors(lam, lmax), lam))
    add(ev, gamma_events(sigma_factors(x, lam, lmax), lam))
    add(ev, trig_events(lam, x, 0))
    # Q_0 = -sin(lam u) cos(lam u): zeros
    for kk in range(1, 400):
        for u in (kk*math.pi/(2*lam),):
            if 0 < u < math.pi/2: ev[round(u, 9)] -= 1
    for i, n in enumerate(S, start=1):
        add(ev, an_events(n, eta if i % 2 == 1 else etab, lam))
    # merge near-equal keys
    out = {}
    for u in sorted(ev):
        if out and abs(u - max(out)) < 1e-7: out[max(out)] += ev[u]
        else: out[u] = ev[u]
    return {u: o for u, o in out.items() if abs(o) > 1e-9}

def label_position(u, eta, lam):
    nu_ = lambda j: eta/lam - (2*j+1)*math.pi/(2*lam)
    etab = math.pi*(lam+1) - eta
    w_ = lambda j: etab/lam - (2*j+1)*math.pi/(2*lam)    # w_j = nu_j(etabar)
    out = []
    for j in range(-8, 40):
        if abs(nu_(j) - u) < 1e-7: out.append(f'nu_{j}')
        if abs(w_(j) - u) < 1e-7: out.append(f'w_{j}')
    for n in range(1, 40):
        if abs(n*math.pi/(2*lam) - u) < 1e-7: out.append(f'u_{n}')
        if abs(math.pi/2 - n*math.pi/(2*lam) - u) < 1e-7: out.append(f'pi/2-u_{n}')
    return '='.join(out) if out else '??'

if __name__ == '__main__':
    import sys
    lam, eta = 6.3, 5.0
    for S in [(), (0,), (1,), (3,), (2, 0)]:
        print(S, {label_position(u, eta, lam): o for u, o in total_events(S, eta, lam).items()})
