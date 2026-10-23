"""Grassmann algebra, (1,1) superspace and supermatrices for the checks of the
chapter on affine Toda theories based on Lie superalgebras (sec-super-toda).

This is a local helper. The other scripts in code/super_toda/ import it; no
other script in the book does. Run it on its own to check the algebra:

    python3 code/super_toda/superspace.py

Conventions (those of Evans-Hollowood 1992 and Evans-Madsen 1997):

* Light-cone coordinates x^+- = (t +- x)/2, so that d_+- = d_t +- d_x is the
  derivative d/dx^+-.  These are the book's d_+- of sec-lagrangians.
* Real Grassmann coordinates theta^+, theta^-; D_+- = d/dtheta^+- - i theta^+- d_+-,
  so D_+-^2 = -i d_+- and {D_+, D_-} = 0.  Q_+- = d/dtheta^+- + i theta^+- d_+-
  anticommute with D_+- and obey Q_+-^2 = i d_+-.
* A scalar superfield is Phi = phi + i theta^+ psi_+ + i theta^- psi_- + i theta^+ theta^- sigma.

Representation.  A Grassmann element is a dict {monomial: coefficient}.  A
monomial is a tuple of generator labels in increasing order of their creation
index; a product is reordered to that order with the sign of the permutation.
d/dtheta acts from the left.  Coefficients are sympy expressions.

Fields live in a jet space.  The bosonic jets phi[a]_{k,l} = d_+^k d_-^l phi^a
and sigma[a]_{k,l} are sympy symbols; the fermionic jets psi_+[a]_{k,l},
psi_-[a]_{k,l} are Grassmann generators.  Total derivatives d_+- act on
coefficients by the chain rule and on generators by raising the jet index.

Supermatrices.  An element of the Grassmann envelope of gl(m|n) is a dict
{(r, s): Grassmann element}, meaning sum_rs c_rs (x) E_rs with c_rs to the
left of the matrix unit.  Matrix units of parity p(r)+p(s) anticommute with odd
coefficients, so (c (x) E)(d (x) F) = (-1)^{|E||d|} cd (x) EF.
"""
import itertools
import sys
import sympy as sp

I = sp.I
_GEN_ORDER = {}      # label -> creation index


def _register(label):
    if label not in _GEN_ORDER:
        _GEN_ORDER[label] = len(_GEN_ORDER)
    return label


def _sort_sign(seq):
    """Sign of the permutation sorting seq by creation index; 0 if repeated."""
    seq = list(seq)
    if len(set(seq)) < len(seq):
        return 0, None
    keys = [_GEN_ORDER[s] for s in seq]
    sign = 1
    for i in range(len(keys)):            # bubble sort, counting swaps
        for j in range(len(keys) - 1 - i):
            if keys[j] > keys[j + 1]:
                keys[j], keys[j + 1] = keys[j + 1], keys[j]
                seq[j], seq[j + 1] = seq[j + 1], seq[j]
                sign = -sign
    return sign, tuple(seq)


class G:
    """Element of a Grassmann algebra with sympy coefficients."""

    def __init__(self, terms=None):
        self.t = {}
        if terms:
            for m, c in terms.items():
                if c != 0:
                    self.t[m] = self.t.get(m, 0) + c

    @staticmethod
    def scalar(c):
        return G({(): sp.sympify(c)})

    @staticmethod
    def gen(label):
        _register(label)
        return G({(label,): sp.Integer(1)})

    def copy(self):
        return G(dict(self.t))

    def __add__(self, o):
        o = _g(o)
        r = dict(self.t)
        for m, c in o.t.items():
            r[m] = r.get(m, 0) + c
        return G(r)

    __radd__ = __add__

    def __neg__(self):
        return G({m: -c for m, c in self.t.items()})

    def __sub__(self, o):
        return self + (-_g(o))

    def __rsub__(self, o):
        return _g(o) + (-self)

    def __mul__(self, o):
        o = _g(o)
        r = {}
        for m1, c1 in self.t.items():
            for m2, c2 in o.t.items():
                s, m = _sort_sign(m1 + m2)
                if s:
                    r[m] = r.get(m, 0) + s * c1 * c2
        return G(r)

    def __rmul__(self, o):
        return _g(o) * self

    def map(self, f):
        return G({m: f(c) for m, c in self.t.items()})

    def simplify(self):
        r = {}
        for m, c in self.t.items():
            c = sp.expand(sp.simplify(c))
            if c != 0:
                r[m] = c
        return G(r)

    def is_zero(self, simplify=True):
        return all(expr_is_zero(c) if simplify else c == 0 for c in self.t.values())

    def even(self):
        return G({m: c for m, c in self.t.items() if len(m) % 2 == 0})

    def odd(self):
        return G({m: c for m, c in self.t.items() if len(m) % 2 == 1})

    def body(self):
        return self.t.get((), sp.Integer(0))

    def coeff(self, mono):
        """Coefficient of a monomial given as a tuple of labels (any order)."""
        s, m = _sort_sign(mono)
        if not s:
            return sp.Integer(0)
        return s * self.t.get(m, sp.Integer(0))

    def __repr__(self):
        return ' + '.join(f'({c})*{"*".join(m) or "1"}' for m, c in self.t.items()) or '0'


def _g(x):
    return x if isinstance(x, G) else G.scalar(x)


def expr_is_zero(c, trials=3):
    """Zero test for a sympy expression in exponentials, hyperbolic functions and
    polynomials: expand after rewriting in exponentials; if sympy cannot decide,
    evaluate at `trials` random rational points and require |value| < 1e-12."""
    c = sp.expand(sp.sympify(c).rewrite(sp.exp))
    if c == 0:
        return True
    c2 = sp.simplify(c)
    if c2 == 0:
        return True
    syms = sorted(c.free_symbols, key=lambda q: q.name)
    import random
    rnd = random.Random(1)
    for _ in range(trials):
        val = c.subs({q: sp.Rational(rnd.randint(-9, 9), rnd.randint(5, 13)) for q in syms})
        if abs(complex(sp.N(val, 30))) > 1e-12:
            return False
    return True


def d_theta(x, label):
    """Left derivative d/d(label) of a Grassmann element."""
    r = {}
    for m, c in x.t.items():
        if label in m:
            i = m.index(label)
            mm = m[:i] + m[i + 1:]
            r[mm] = r.get(mm, 0) + (-1) ** i * c
    return G(r)


# ------------------------------------------------------------------ jet space
TP, TM = G.gen('thp'), G.gen('thm')        # theta^+, theta^- (created first)


class Jets:
    """Jet variables of r scalar superfields and total derivatives."""

    def __init__(self, r, name='f'):
        self.r = r
        self.name = name

    def phi(self, a, k=0, l=0):
        return sp.Symbol(f'{self.name}phi{a}_{k}_{l}')

    def sig(self, a, k=0, l=0):
        return sp.Symbol(f'{self.name}sig{a}_{k}_{l}')

    def psi(self, s, a, k=0, l=0):
        return G.gen(f'{self.name}psi{s}{a}_{k}_{l}')

    def _shift_symbol(self, s, sign):
        n = s.name
        if not n.startswith(self.name):
            return sp.Integer(0)
        body = n[len(self.name):]
        if body.startswith('phi') or body.startswith('sig'):
            kind, rest = body[:3], body[3:]
            a, k, l = rest.split('_')
            k, l = int(k), int(l)
            if sign > 0:
                k += 1
            else:
                l += 1
            return sp.Symbol(f'{self.name}{kind}{a}_{k}_{l}')
        return sp.Integer(0)

    def _shift_gen(self, label, sign):
        if not label.startswith(self.name + 'psi'):
            return None
        rest = label[len(self.name) + 3:]
        s, rest = rest[0], rest[1:]
        a, k, l = rest.split('_')
        k, l = int(k), int(l)
        if sign > 0:
            k += 1
        else:
            l += 1
        return G.gen(f'{self.name}psi{s}{a}_{k}_{l}')

    def dcoef(self, c, sign):
        """Total derivative of a sympy coefficient."""
        c = sp.sympify(c)
        res = sp.Integer(0)
        for s in c.free_symbols:
            sh = self._shift_symbol(s, sign)
            if sh != 0:
                res += sp.diff(c, s) * sh
        return res

    def d(self, x, sign):
        """Total derivative d_+ (sign=+1) or d_- (sign=-1) of a Grassmann element."""
        x = _g(x)
        out = G()
        for m, c in x.t.items():
            dc = self.dcoef(c, sign)
            if dc != 0:
                out = out + G({m: dc})
            for i, lab in enumerate(m):
                sh = self._shift_gen(lab, sign)
                if sh is not None:
                    left = G({m[:i]: sp.Integer(1)})
                    right = G({m[i + 1:]: sp.Integer(1)})
                    # reorder: left * sh * right keeps the position of the generator
                    out = out + c * (left * sh * right)
        return out

    def D(self, x, sign):
        """Superspace derivative D_+- = d/dtheta^+- - i theta^+- d_+-."""
        th, lab = (TP, 'thp') if sign > 0 else (TM, 'thm')
        return d_theta(x, lab) - I * (th * self.d(x, sign))

    def Q(self, x, sign):
        th, lab = (TP, 'thp') if sign > 0 else (TM, 'thm')
        return d_theta(x, lab) + I * (th * self.d(x, sign))

    def superfield(self, a):
        """Phi^a = phi + i th+ psi_+ + i th- psi_- + i th+ th- sigma."""
        return (G.scalar(self.phi(a)) + I * (TP * self.psi('p', a))
                + I * (TM * self.psi('m', a)) + I * (TP * TM) * self.sig(a))

    def component(self, x, which):
        """Theta components: x = x0 + th+ x+ + th- x- + th+ th- x2."""
        if which == '0':
            return G({m: c for m, c in x.t.items() if 'thp' not in m and 'thm' not in m})
        if which == '+':
            y = d_theta(x, 'thp')
            return G({m: c for m, c in y.t.items() if 'thm' not in m})
        if which == '-':
            y = d_theta(x, 'thm')
            return G({m: c for m, c in y.t.items() if 'thp' not in m})
        if which == '2':
            y = d_theta(d_theta(x, 'thp'), 'thm')
            return G({m: c for m, c in y.t.items() if 'thp' not in m and 'thm' not in m})
        raise ValueError(which)


def sexp(x):
    """exp of a Grassmann element: exp(body) * sum_k n^k/k! with n nilpotent."""
    x = _g(x)
    b = x.body()
    n = x - G.scalar(b)
    out, term, k = G.scalar(1), G.scalar(1), 0
    while True:
        k += 1
        term = term * n
        if not term.t:
            break
        out = out + term.map(lambda c, k=k: c / sp.factorial(k))
    return out.map(lambda c: sp.exp(b) * c)


def berezin(x):
    """Component Lagrangian int d^2theta x := i (d/dth- d/dth+ x)|_{th=0}, so that
    int d^2theta th+ th- = i.  The factor i makes the action of real superfields
    real, since (th+ th-)^* = th- th+ = -th+ th-.  With this choice the
    Evans-Madsen action gives the kinetic term +1/2 d_+phi.d_-phi
    (checked in components.py)."""
    y = I * d_theta(d_theta(x, 'thp'), 'thm')
    return G({m: c for m, c in y.t.items() if 'thp' not in m and 'thm' not in m})


# ------------------------------------------------------------- supermatrices
class SMat:
    """Element sum_rs c_rs (x) E_rs of the Grassmann envelope of gl(m|n).
    parity: list of row parities (0 even, 1 odd)."""

    def __init__(self, parity, entries=None):
        self.p = list(parity)
        self.e = {}
        if entries:
            for k, v in entries.items():
                v = _g(v)
                if v.t:
                    self.e[k] = v

    def unit_parity(self, r, s):
        return (self.p[r] + self.p[s]) % 2

    def __add__(self, o):
        r = dict(self.e)
        for k, v in o.e.items():
            r[k] = r.get(k, G()) + v
        return SMat(self.p, r)

    def __neg__(self):
        return SMat(self.p, {k: -v for k, v in self.e.items()})

    def __sub__(self, o):
        return self + (-o)

    def scale(self, c):
        """Left multiplication by a Grassmann element c: c (sum d E) = sum cd E."""
        c = _g(c)
        return SMat(self.p, {k: c * v for k, v in self.e.items()})

    def __mul__(self, o):
        r = {}
        for (a, b), c in self.e.items():
            for (b2, d), dd in o.e.items():
                if b != b2:
                    continue
                pe = self.unit_parity(a, b)
                term = c * (dd.even() + (-1) ** pe * dd.odd())
                r[(a, d)] = r.get((a, d), G()) + term
        return SMat(self.p, r)

    def map(self, f):
        return SMat(self.p, {k: f(v) for k, v in self.e.items()})

    def is_zero(self):
        return all(v.is_zero() for v in self.e.values())


def unit(parity, r, s, c=1):
    return SMat(parity, {(r, s): _g(c)})


def gbracket(X, Y, px, py):
    """Graded bracket [X, Y} = XY - (-1)^{px py} YX of homogeneous elements
    of parities px, py (total parity: coefficient plus matrix unit)."""
    return X * Y - (Y * X).scale(sp.Integer((-1) ** (px * py)))


# ------------------------------------------------- component Lagrangians
def dot(Gram, u, v):
    return (sp.Matrix(u).T * Gram * sp.Matrix(v))[0, 0]


def superspace_lagrangian(J, Gram, roots, odd, mu):
    """The Evans-Madsen superspace Lagrangian for root vectors (coordinate
    columns), a Gram matrix, parities odd[j] and coefficients mu[j]."""
    r = Gram.shape[0]
    Phi = [J.superfield(a) for a in range(r)]
    Dp = [J.D(f, +1) for f in Phi]
    Dm = [J.D(f, -1) for f in Phi]
    L = G()
    for a in range(r):
        for b in range(r):
            if Gram[a, b] != 0:
                L = L + (Dp[a] * Dm[b]).map(lambda c, g=Gram[a, b]: I / 2 * g * c)
    for al, o, m in zip(roots, odd, mu):
        x = G()
        w = list(Gram * sp.Matrix(al))
        for a in range(r):
            if w[a] != 0:
                x = x + Phi[a].map(lambda c, k=w[a]: k * c)
        e = sexp(x)
        L = L + (e.map(lambda c: m * c) if o else (TP * TM) * e.map(lambda c: I * m * c))
    return L


def eliminate_sigma(J, Lc, r):
    """Solve the (linear) auxiliary-field equations and substitute."""
    sig = [J.sig(a) for a in range(r)]
    body = Lc.body()
    sol = sp.solve([sp.diff(body, s) for s in sig], sig, dict=True)[0]
    return Lc.map(lambda c: sp.expand(c.subs(sol))), sol


def psi_lin(J, s, Gram, v):
    """v.psi_s = sum_ab v^a Gram_ab psi_s^b as a Grassmann element."""
    w = list(Gram * sp.Matrix(v))
    out = G()
    for a in range(len(w)):
        if w[a] != 0:
            out = out + J.psi(s, a).map(lambda c, k=w[a]: k * c)
    return out


def expected_components(J, Gram, roots, odd, mu):
    """eq-super-components (Evans-Madsen (3.4)) built directly."""
    r = Gram.shape[0]
    phi = sp.Matrix([J.phi(a) for a in range(r)])
    dpp = sp.Matrix([J.phi(a, 1, 0) for a in range(r)])
    dmp = sp.Matrix([J.phi(a, 0, 1) for a in range(r)])
    Lb = (dpp.T * Gram * dmp)[0, 0] / 2
    for al, o, m in zip(roots, odd, mu):
        if not o:
            Lb -= m * sp.exp(dot(Gram, al, phi))
    for (al, o, m), (bl, p, n) in itertools.product(zip(roots, odd, mu), repeat=2):
        if o and p:
            Lb -= sp.Rational(1, 2) * m * n * dot(Gram, al, bl) * sp.exp(dot(Gram, al, phi) + dot(Gram, bl, phi))
    out = G.scalar(Lb)
    for a in range(r):
        for b in range(r):
            if Gram[a, b] != 0:
                out = out + (J.psi('p', a) * J.psi('p', b, 0, 1)).map(lambda c, g=Gram[a, b]: I / 2 * g * c)
                out = out + (J.psi('m', a) * J.psi('m', b, 1, 0)).map(lambda c, g=Gram[a, b]: I / 2 * g * c)
    for al, o, m in zip(roots, odd, mu):
        if o:
            out = out + (psi_lin(J, 'p', Gram, al) * psi_lin(J, 'm', Gram, al)).map(
                lambda c, m=m, al=al: I * m * sp.exp(dot(Gram, al, phi)) * c)
    return out


def components(Gram, roots, odd, mu, name='f'):
    J = Jets(Gram.shape[0], name)
    L = superspace_lagrangian(J, Gram, roots, odd, mu)
    Lc = berezin(L)
    Lc, sol = eliminate_sigma(J, Lc, Gram.shape[0])
    return J, Lc, sol


def potential(J, Lc, r):
    """U(phi) = -(Lagrangian at zero fermions and zero derivatives)."""
    body = Lc.body()
    subs = {J.phi(a, 1, 0): 0 for a in range(r)}
    subs.update({J.phi(a, 0, 1): 0 for a in range(r)})
    return -sp.expand(body.subs(subs))


def mass_matrices(Gram, roots, odd, mu):
    """Boson mass^2 matrix Gram^{-1} Hess U and fermion mass matrix
    Gram^{-1} K, K = sum_odd mu_j (Gram alpha_j)(Gram alpha_j)^T, at phi = 0."""
    r = Gram.shape[0]
    J, Lc, _ = components(Gram, roots, odd, mu)
    U = potential(J, Lc, r)
    phis = [J.phi(a) for a in range(r)]
    z = {p: 0 for p in phis}
    grad = sp.Matrix([sp.diff(U, p) for p in phis]).subs(z)
    H = sp.Matrix(r, r, lambda a, b: sp.diff(U, phis[a], phis[b])).subs(z)
    K = sp.zeros(r, r)
    for al, o, m in zip(roots, odd, mu):
        if o:
            w = Gram * sp.Matrix(al)
            K += m * w * w.T
    return sp.simplify(grad), sp.simplify(Gram.inv() * H), sp.simplify(Gram.inv() * K), U, J, Lc


# ------------------------------------------------------------------ self-test
def _report(name, ok, fails):
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok:
        fails.append(name)


def self_test():
    fails = []
    a, b, c = G.gen('ta'), G.gen('tb'), G.gen('tc')
    _report('Grassmann generators anticommute and square to zero',
            (a * b + b * a).is_zero() and (a * a).is_zero()
            and ((a * b) * c - a * (b * c)).is_zero(), fails)
    y = G.scalar(sp.Symbol('u')) + b * c
    _report('d/dtheta is an odd derivation acting from the left',
            (d_theta(a * b * c, 'tb') + a * c).is_zero()
            and (d_theta(a * y, 'ta') - y).is_zero()
            and (d_theta(y * a, 'ta') - y).is_zero(), fails)
    J = Jets(1)
    Phi = J.superfield(0)
    DDp = J.D(J.D(Phi, +1), +1) + I * J.d(Phi, +1)
    DDm = J.D(J.D(Phi, -1), -1) + I * J.d(Phi, -1)
    anti = J.D(J.D(Phi, -1), +1) + J.D(J.D(Phi, +1), -1)
    _report('D_+^2 = -i d_+, D_-^2 = -i d_-, {D_+, D_-} = 0 on a superfield',
            DDp.is_zero() and DDm.is_zero() and anti.is_zero(), fails)
    QD = J.Q(J.D(Phi, +1), +1) + J.D(J.Q(Phi, +1), +1)
    QD2 = J.Q(J.D(Phi, -1), +1) + J.D(J.Q(Phi, +1), -1)
    QQ = J.Q(J.Q(Phi, +1), +1) - I * J.d(Phi, +1)
    _report('{Q_+, D_+} = {Q_+, D_-} = 0 and Q_+^2 = i d_+',
            QD.is_zero() and QD2.is_zero() and QQ.is_zero(), fails)
    th = TP * TM
    _report('berezin(th+ th-) = i and berezin(1) = berezin(th+) = 0',
            sp.simplify(berezin(th).body() - I) == 0 and berezin(G.scalar(1)).is_zero()
            and berezin(TP).is_zero(), fails)
    f = sp.Symbol('fphi0_0_0')
    e = sexp(G.scalar(f) + a * b + c)
    _report('sexp agrees with the Taylor series on a nilpotent element',
            (e - (G.scalar(1) + a * b + c + a * b * c).map(lambda q: q * sp.exp(f))
             ).is_zero(), fails)
    # supermatrices: gl(1|1) graded Jacobi and supertrace of brackets
    par = [0, 1]
    E11, E22 = unit(par, 0, 0), unit(par, 1, 1)
    E12, E21 = unit(par, 0, 1), unit(par, 1, 0)
    anticomm = E12 * E21 + E21 * E12
    _report('gl(1|1): {E12, E21} = E11 + E22 for odd matrix units',
            (anticomm - E11 - E22).is_zero(), fails)
    # envelope: odd coefficients times odd units give even elements whose
    # commutator is [a (x) X, b (x) Y] = (-1)^{|X||b|} ab (x) [X, Y}
    X, Y = E12.scale(a), E21.scale(b)
    lhs = X * Y - Y * X
    rhs = (E11 + E22).scale(-(a * b))
    _report('envelope: [a E12, b E21] = -a b {E12, E21} with odd a, b',
            (lhs - rhs).is_zero(), fails)
    # graded Jacobi on random homogeneous elements of gl(2|1)
    par3 = [0, 0, 1]
    elems = [(unit(par3, 0, 1), 0), (unit(par3, 0, 2), 1), (unit(par3, 2, 1), 1),
             (unit(par3, 1, 2), 1), (unit(par3, 2, 2), 0)]
    ok = True
    for (A, pa), (B, pb), (C, pc) in itertools.product(elems, repeat=3):
        t1 = gbracket(A, gbracket(B, C, pb, pc), pa, (pb + pc) % 2)
        t2 = gbracket(gbracket(A, B, pa, pb), C, (pa + pb) % 2, pc)
        t3 = gbracket(B, gbracket(A, C, pa, pc), pb, (pa + pc) % 2).scale(
            sp.Integer((-1) ** (pa * pb)))
        if not (t1 - t2 - t3).is_zero():
            ok = False
    _report('gl(2|1): graded Jacobi identity on matrix units', ok, fails)
    return fails


if __name__ == '__main__':
    fails = self_test()
    sys.exit(1 if fails else 0)
