"""Common tools for problem 8 (residues of boundary Coleman-Thun diagrams).

Blocks (book eq-blocks): (x)(theta) = sinh(theta/2 + i pi x/2h)/sinh(theta/2 - i pi x/2h), pole at theta = i pi x/h.
A Prod is sign * prod (c + d B)^e over blocks, kept as a list so that coinciding blocks are handled exactly.
Laurent data at theta0 = i pi p/h are computed exactly from the block structure (order of the pole, leading
coefficient), then evaluated with mpmath.

Conventions (book ch. 28, GZ): near a boundary bound-state pole  K ~ (i/2) g^2/(theta - i v);
bulk: S ~ i Gamma^2/(theta - i u).
"""
from fractions import Fraction as Fr
from collections import Counter
import mpmath as mp

mp.mp.dps = 40
I = mp.mpc(0, 1)
PI = mp.pi


def blk(x, th, h):
    return mp.sinh(th/2 + I*PI*x/(2*h))/mp.sinh(th/2 - I*PI*x/(2*h))


class Prod:
    """sign * prod_k (c_k + d_k B)^{e_k}; terms: list of ((c, d), e)."""

    def __init__(self, h, terms=None, sign=1):
        self.h, self.sign = h, sign
        self.t = Counter()
        for (c, d), e in (terms.items() if isinstance(terms, dict) else (terms or [])):
            self.t[(Fr(c), Fr(d))] += e

    def copy(self):
        r = Prod(self.h, sign=self.sign)
        r.t = Counter(self.t)
        return r

    def mul(self, x, e=1):
        self.t[(Fr(x[0]), Fr(x[1]))] += e
        return self

    def __mul__(self, o):
        r = self.copy()
        r.sign *= o.sign
        for k, e in o.t.items():
            r.t[k] += e
        return r

    def inv(self):
        r = Prod(self.h, sign=self.sign)
        for k, e in self.t.items():
            r.t[k] = -e
        return r

    def __call__(self, th, B):
        v = mp.mpf(self.sign)
        for (c, d), e in self.t.items():
            if e:
                v *= blk(c + d*mp.mpf(B), th, self.h)**e
        return v

    # -------- exact local analysis at theta0 = i pi p / h (p a Fraction or mpf) with B a Fraction
    def order_at(self, p, B):
        """Net order (pole > 0, zero < 0) at theta = i pi p/h, exact for Fraction B, p."""
        h = self.h
        o = 0
        for (c, d), e in self.t.items():
            if not e:
                continue
            x = c + d*B
            if (x - p) % (2*h) == 0 and (x % h) != 0:
                o += e        # pole of (x)
            if (x + p) % (2*h) == 0 and (x % h) != 0:
                o -= e        # zero of (x)
        return o

    def laurent(self, p, B, Bnum=None):
        """Return (order, leading coefficient) of the expansion in (theta - theta0) at theta0 = i pi p/h.
        K ~ coef * (theta - theta0)^(-order).  Exact block-by-block: each block is factored into its local
        singular part times a regular value."""
        h = self.h
        Bn = mp.mpf(B.numerator)/B.denominator if Bnum is None else Bnum
        pn = mp.mpf(p.numerator)/p.denominator
        th0 = I*PI*pn/h
        coef = mp.mpf(self.sign)
        order = 0
        for (c, d), e in self.t.items():
            if not e:
                continue
            x = c + d*B
            xn = mp.mpf(c.numerator)/c.denominator + mp.mpf(d.numerator)/d.denominator*Bn
            num = mp.sinh(th0/2 + I*PI*xn/(2*h))
            den = mp.sinh(th0/2 - I*PI*xn/(2*h))
            ispole = (x - p) % (2*h) == 0 and (x % h) != 0
            iszero = (x + p) % (2*h) == 0 and (x % h) != 0
            if ispole and iszero:
                raise ValueError('degenerate block')
            if ispole:
                # den ~ cosh(th0/2 - i pi x/2h) * (theta - th0)/2
                dden = mp.cosh(th0/2 - I*PI*xn/(2*h))/2
                coef *= (num/dden)**e
                order += e
            elif iszero:
                dnum = mp.cosh(th0/2 + I*PI*xn/(2*h))/2
                coef *= (dnum/den)**e
                order -= e
            else:
                coef *= (num/den)**e
        return order, coef

    def strip_poles(self, B, half=True):
        """Physical-strip poles 0 < p < h/2 (boundary) of the product, as {p: order}, exact."""
        h = self.h
        cand = set()
        for (c, d), e in self.t.items():
            if e:
                x = (c + d*B) % (2*h)
                for y in (x, (-x) % (2*h)):
                    cand.add(y)
        top = Fr(h, 2) if half else Fr(h)
        res = {}
        for y in sorted(cand):
            if 0 < y < top:
                o = self.order_at(y, B)
                if o != 0:
                    res[y] = o
        return res


def curly(p, h):
    """{p} = (p-1)(p+1)/((p-1+B)(p+1-B))."""
    return Prod(h, [((p - 1, 0), 1), ((p + 1, 0), 1), ((p - 1, 1), -1), ((p + 1, -1), -1)])


def S(a, b, n):
    h = n + 1
    if a + b > h:
        a, b = h - a, h - b
    r = Prod(h)
    for p in range(abs(a - b) + 1, a + b, 2):
        r = r*curly(p, h)
    return r


def K_uniform(a, n):
    r = Prod(n + 1)
    for c in range(1, a + 1):
        r.mul((c - 1, 0)).mul((c - n - 1, 0)).mul((-c, Fr(1, 2))).mul((-c - n, Fr(-1, 2)))
    return r


def M_sol(a, b, n):
    h = n + 1
    r = Prod(h)
    for k in range(1, a + 1):
        r.mul((a + b + Fr(h, 2) - 2*k, Fr(1, 2))).mul((b - a - Fr(h, 2) + 2*k, Fr(-1, 2)))
        r.mul((a + b - Fr(h, 2) - 2*k, Fr(1, 2)), -1).mul((b - a + Fr(h, 2) + 2*k, Fr(-1, 2)), -1)
    return r


def K1nm(nn, mm):
    """a_2, DG (6.7): K_1 on b_{n,m}."""
    return K_uniform(1, 2)*Prod(3, [((Fr(1, 2), -Fr(2*nn + 1, 2)), 1), ((Fr(5, 2), Fr(2*mm + 1, 2)), 1),
                                    ((Fr(1, 2), -Fr(2*nn - 1, 2)), 1), ((Fr(5, 2), Fr(2*mm - 1, 2)), 1),
                                    ((Fr(3, 2), Fr(2*nn - 1, 2)), 1), ((Fr(3, 2), -Fr(2*mm - 1, 2)), 1),
                                    ((-Fr(5, 2), Fr(2*nn + 1, 2)), 1), ((-Fr(1, 2), -Fr(2*mm + 1, 2)), 1)])


def Ka2(a, nn, mm):
    return K1nm(nn, mm) if a == 1 else K1nm(mm, nn)


def F(k, B):
    d = PI*B/6
    return mp.sin(k*d)*mp.cos(PI/6 - k*d)/mp.sin(d)


def E(nn, mm, B):
    return F(nn, B) + F(mm, B)


def fnum(x):
    return mp.mpf(x.numerator)/x.denominator if isinstance(x, Fr) else mp.mpf(x)
