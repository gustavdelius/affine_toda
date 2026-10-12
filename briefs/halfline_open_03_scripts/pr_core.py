"""Particle reflection amplitudes of real-coupling a_n^(1) theory on the half-line (Part VII, ch. 35:
"""
from fractions import Fraction as Fr
from collections import Counter
import mpmath as mp

mp.mp.dps = 30
I = mp.mpc(0, 1)
PI = mp.pi


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


def blk(x, th, h):
    """Book block (x), eq-blocks."""
    return mp.sinh(th/2 + I*PI*x/(2*h))/mp.sinh(th/2 - I*PI*x/(2*h))


class Prod:
    """sign * prod (c + d B)^e, stored as {(c, d): e} with Fractions c, d."""
    def __init__(self, h, terms=None, sign=1):
        self.h, self.sign, self.t = h, sign, Counter()
        for (c, d), e in (terms.items() if isinstance(terms, dict) else (terms or [])):
            self.t[(Fr(c), Fr(d))] += e

    def mul(self, x, e=1):          # multiply by (x)^e, x = (c, d) meaning c + d B
        self.t[(Fr(x[0]), Fr(x[1]))] += e
        return self

    def __mul__(self, o):
        r = Prod(self.h, dict(self.t), self.sign*o.sign)
        for k, e in o.t.items():
            r.t[k] += e
        return r

    def dual(self):                 # B -> 2 - B
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
        """Net order at theta = i pi p/h for 0 < p < h/2 (physical strip), B a Fraction; >0 pole, <0 zero."""
        h, acc = self.h, Counter()
        for (c, d), e in self.t.items():
            x = (c + d*B) % (2*h)
            if x > h:
                x -= 2*h
            if x != 0 and x != h:
                acc[x] += e
                acc[-x] -= e
        return {p: o for p, o in sorted(acc.items()) if 0 < p < Fr(h, 2) and o != 0}


def curly(p, h):
    """{p} = (p-1)(p+1)(-p+1-B)(-p-1+B)."""
    return Prod(h, {(p - 1, 0): 1, (p + 1, 0): 1, (-p + 1, -1): 1, (-p - 1, 1): 1})


def S(a, b, n):
    h = n + 1
    if a + b > h:
        a, b = h - a, h - b
    r = Prod(h)
    for p in range(abs(a - b) + 1, a + b, 2):
        r = r*curly(p, h)
    return r


def K_uniform(a, n):           # DG (3.16)
    r = Prod(n + 1)
    for c in range(1, a + 1):
        r.mul((c - 1, 0)).mul((c - n - 1, 0)).mul((-c, Fr(1, 2))).mul((-c - n, Fr(-1, 2)))
    return r


def K_neumann(a, n):           # DG (3.20)
    r = Prod(n + 1)
    for c in range(1, a + 1):
        r.mul((c - 1, 0)).mul((c - n - 1, 0)).mul((-c + 1, Fr(-1, 2))).mul((-c + n + 1, Fr(1, 2)))
    return r


def M_sol(a, b, n):            # DG (5.4)
    h = n + 1
    r = Prod(h)
    for k in range(1, a + 1):
        r.mul((a + b + Fr(h, 2) - 2*k, Fr(1, 2))).mul((b - a - Fr(h, 2) + 2*k, Fr(-1, 2)))
        r.mul((a + b - Fr(h, 2) - 2*k, Fr(1, 2)), -1).mul((b - a + Fr(h, 2) + 2*k, Fr(-1, 2)), -1)
    return r


THS = [mp.mpc('0.37', '0.11'), mp.mpc('-0.83', '0.29'), mp.mpc('1.21', '-0.47')]
BS = [mp.mpf('0.317'), mp.mpf('1.43')]
TOL = mp.mpf('1e-20')


def close(f, g):
    return max(abs(f(th, B) - g(th, B)) for th in THS for B in BS) < TOL


def axioms(K, n, label, conj=True):
    """K(a) -> amplitude; checks unitarity, crossing-unitarity (book and DG forms), K_a = K_{h-a}, bootstrap."""
    h = n + 1
    ok_u = all(close(lambda t, B: K(a)(t, B)*K(a)(-t, B), lambda t, B: 1) for a in range(1, h))
    ok_c = all(close(lambda t, B: K(a)(I*PI/2 - t, B), lambda t, B: S(a, h - a, n)(2*t, B)*K(h - a)(I*PI/2 + t, B))
               for a in range(1, h))
    ok_dg = all(close(lambda t, B: K(a)(t, B)*K(h - a)(t + I*PI, B), lambda t, B: S(a, a, n)(2*t, B))
                for a in range(1, h))
    ok_cc = all(close(K(a), K(h - a)) for a in range(1, h)) if conj else True
    ok_b = True
    for a in range(1, h):
        for b in range(1, h):
            if a + b < h:
                c, ua, ub = a + b, PI*b/h, PI*a/h
            elif a + b > h:
                c, ua, ub = a + b - h, PI*(h - b)/h, PI*(h - a)/h
            else:
                continue
            ok_b &= close(K(c), lambda t, B: K(a)(t + I*ua, B)*S(a, b, n)(2*t + I*ua - I*ub, B)*K(b)(t - I*ub, B))
    return ok_u, ok_c, ok_dg, ok_cc, ok_b


