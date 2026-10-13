"""Boundary bootstrap for the vector-like d_n particles: K_c(th) = K_a(th + i ub_ac^b) K_b(th - i ub_bc^a)
S_ab(2th + i(ub_ac^b - ub_bc^a)), fusings detected from the book's S-matrix poles (setup('d', n)); K_a built from K_1 along a
chain and checked at every other fusing and in crossing-unitarity. Applied to the book's uniform K_1 (validation) and to
the s = +-1 amplitudes of this work (with and without the CDD factor)."""
import itertools
import mpmath as mp
from t17_cont import K1, bb, setup

mp.mp.dps = 30


def Smp(bl, h, B, th):
    r = mp.mpf(1)
    for x, k in bl.items():
        r *= (bb(x - 1, th, h)*bb(x + 1, th, h)/(bb(x - 1 + B, th, h)*bb(x + 1 - B, th, h)))**k
    return r


def fusings(n):
    h, bl, mass, _ = setup('d', n)
    vec = list(range(n - 2))              # nodes 0..n-3 = particles 1..n-2
    out = []
    for a, b, c in itertools.product(vec, vec, vec):
        ma, mb, mc = mass[a], mass[b], mass[c]
        if not (abs(ma - mb) < mc < ma + mb):
            continue
        cu = (mc**2 - ma**2 - mb**2)/(2*ma*mb)
        u = mp.acos(cu)
        x = u*h/mp.pi
        if abs(x - mp.nint(x)) > 1e-9:
            continue
        u = mp.pi*mp.nint(x)/h
        eps = mp.mpf('1e-12')
        res = eps*Smp(bl[(a, b)], h, mp.mpf('0.37'), 1j*u + eps)
        if abs(res) < 1e-6 or abs(res) > 1e6 or not (mp.re(-1j*res) > 0 and abs(mp.im(-1j*res)) < 1e-6*abs(res)):
            continue              # keep only simple forward-channel poles (-i Res > 0)
        al = mp.atan2(mb*mp.sin(u), ma + mb*mp.cos(u))   # theta_a = th + i al, theta_b = th - i (u - al)
        out.append((a, b, c, al, u - al))
    return h, bl, out


def check(n, Kfun, B):
    h, bl, fus = fusings(n)
    K = {0: Kfun}
    # build K_c along the chain using fusings with a = 0 (particle 1)
    changed = True
    while changed:
        changed = False
        for a, b, c, al, be in fus:
            if a in K and b in K and c not in K:
                K[c] = (lambda t, a=a, b=b, al=al, be=be: K[a](t + 1j*al)*K[b](t - 1j*be)
                        * Smp(bl[(a, b)], h, B, 2*t + 1j*(al - be)))
                changed = True
    dev = 0
    ths = [mp.mpc('0.31', '0.17'), mp.mpc('-0.7', '0.05')]
    for a, b, c, al, be in fus:
        for t in ths:
            lhs = K[c](t)
            rhs = K[a](t + 1j*al)*K[b](t - 1j*be)*Smp(bl[(a, b)], h, B, 2*t + 1j*(al - be))
            dev = max(dev, abs(lhs/rhs - 1))
    for a in K:
        for t in ths:
            dev = max(dev, abs(K[a](t)*K[a](t + 1j*mp.pi)/Smp(bl[(a, a)], h, B, 2*t) - 1), abs(K[a](t)*K[a](-t) - 1))
    return len(fus), sorted(K), dev


if __name__ == '__main__':
    B = mp.mpf('0.43')
    def K1u(t, n):
        h = 2*n - 2
        b = lambda x: bb(x, t, h)
        return b(h - B/2)*b(h/2 - 1)*b(h/2 + 1 - B/2)/(b(h - 1)*b(1 - B/2)*b(h/2)*b(h/2 - B/2))
    for n in (4, 5, 6, 7):
        r = check(n, lambda t: K1u(t, n), B)
        print(f'd_{n} uniform (book eq-dn-uniform-amplitude): {r[0]} fusings, particles {r[1]}, max dev {mp.nstr(r[2], 3)}')
        for sig in (1, -1):
            for cdd in ((), (0, 1, n - 2)):
                r = check(n, lambda t: K1(sig, t, B, n, cdd), B)
                print(f'd_{n} sig={sig:+d} CDD k={cdd}: {r[0]} fusings, particles {r[1]}, max dev {mp.nstr(r[2], 3)}',
                      flush=True)
