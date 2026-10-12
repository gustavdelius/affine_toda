"""Particle-level (real coupling) checks for d_n^(1) boundary amplitudes, book blocks (eq-blocks):
(x) = sinh(th/2 + i pi x/2h)/sinh(th/2 - i pi x/2h), {x} = (x-1)(x+1)/((x-1+B)(x+1-B)), h = 2n-2.
d_n S-matrix (BCDS): vector-like particles a, b <= n-2: S_ab = prod_{p=|a-b|+1, step 2}^{a+b-1} {p}{h-p};
spinors s, s' (n-1, n): S_as = prod_{k=1}^{a} {h/2 - a - 1 + 2k} (a <= n-2);
S_ss = S_s's' = prod_{k} {4k-1}... built from the standard formulas below (checked by unitarity/crossing/bootstrap)."""
from fractions import Fraction as Fr
from collections import Counter
import mpmath as mp

mp.mp.dps = 30
I = mp.mpc(0, 1)
PI = mp.pi


def blk(x, th, h):
    return mp.sinh(th/2 + I*PI*x/(2*h))/mp.sinh(th/2 - I*PI*x/(2*h))


class Prod:
    """sign * prod (c + d B)^e"""
    def __init__(self, h, terms=None, sign=1):
        self.h, self.sign, self.t = h, sign, Counter()
        for (c, d), e in (terms.items() if isinstance(terms, dict) else (terms or [])):
            self.t[(Fr(c), Fr(d))] += e

    def mul(self, x, e=1):
        self.t[(Fr(x[0]), Fr(x[1]))] += e
        return self

    def __mul__(self, o):
        r = Prod(self.h, dict(self.t), self.sign*o.sign)
        for k, e in o.t.items():
            r.t[k] += e
        return r

    def inv(self):
        return Prod(self.h, {k: -e for k, e in self.t.items()}, self.sign)

    def dual(self):
        r = Prod(self.h, sign=self.sign)
        for (c, d), e in self.t.items():
            r.t[(c + 2*d, -d)] += e
        return r

    def __call__(self, th, B):
        v = mp.mpf(self.sign)
        for (c, d), e in self.t.items():
            if e:
                v *= blk(c + d*mp.mpf(B), th, self.h)**e
        return v

    def poles(self, B):
        h, acc = self.h, Counter()
        for (c, d), e in self.t.items():
            x = (c + d*B) % (2*h)
            if x > h:
                x -= 2*h
            if x != 0 and x != h:
                acc[x] += e
                acc[-x] -= e
        return {p: o for p, o in sorted(acc.items()) if 0 < p < Fr(h, 2) and o != 0}

    def __repr__(self):
        s = '-' if self.sign < 0 else ''
        return s + ''.join(f'({c}{"+" if d >= 0 else ""}{d}B)' + (f'^{e}' if e != 1 else '')
                           for (c, d), e in sorted(self.t.items()) if e)


def curly(p, h):
    return Prod(h, {(p - 1, 0): 1, (p + 1, 0): 1, (-p + 1, -1): 1, (-p - 1, 1): 1})


def S_dn(a, b, n):
    """d_n S-matrix; a, b in 1..n (n-1, n spinors)."""
    h = 2*n - 2
    r = Prod(h)
    sp = (n - 1, n)
    if a not in sp and b not in sp:
        for p in range(abs(a - b) + 1, a + b, 2):
            r = r*curly(p, h)*curly(h - p, h)
        return r
    if a in sp and b not in sp:
        a, b = b, a
    if a not in sp:            # vector-like a with spinor b
        for k in range(1, a + 1):
            r = r*curly(h//2 - a - 1 + 2*k, h)
        return r
    # spinor-spinor: S_ss = prod_{k=1}^{[ (n)/2 ]} {4k - 1 - ...}; use standard BCDS:
    # n even: S_ss = prod_{k=1}^{n/2} {4k-3}, S_ss' = prod_{k=1}^{n/2-1} {4k-1}; n odd: S_ss' = prod_{k=1}^{(n-1)/2} {4k-3},
    # S_ss = prod_{k=1}^{(n-1)/2} {4k-1}
    same = (a == b)
    if n % 2 == 0:
        ks = range(1, n//2 + 1) if same else range(1, n//2)
        off = 3 if same else 1
    else:
        ks = range(1, (n - 1)//2 + 1)
        off = 1 if same else 3
    for k in ks:
        r = r*curly(4*k - off, h)
    return r


def close(f, g, tol=mp.mpf('1e-20')):
    return abs(f - g) < tol*(1 + abs(g))


THS = [mp.mpc('0.37', '0.11'), mp.mpc('-0.83', '0.29'), mp.mpc('1.21', '-0.47')]
BS = [mp.mpf('0.317'), mp.mpf('1.43')]


def bulk_checks(n):
    h = 2*n - 2
    ok = True
    for a in range(1, n + 1):
        for b in range(1, n + 1):
            S = S_dn(a, b, n)
            for th in THS:
                for B in BS:
                    ok &= close(S(th, B)*S(-th, B), 1)
                    # crossing: S_ab(i pi - th) = S_{a bbar}(th); for d_n: bar = id except n odd swaps spinors
                    bb = b if (n % 2 == 0 or b < n - 1) else (2*n - 1 - b)
                    ok &= close(S(I*PI - th, B), S_dn(a, bb, n)(th, B))
    return ok
