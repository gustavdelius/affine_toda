"""Particle reflection amplitudes of real-coupling a_n^(1) theory on the half-line (Part VII, ch. 35:
sec-particle-reflection).  Blocks (eq-blocks): (x) = sinh(theta/2 + i pi x/2h)/sinh(theta/2 - i pi x/2h),
{x} = (x-1)(x+1)/((x-1+B)(x+1-B)), h = n+1; S_ab = prod_{p=|a-b|+1, step 2}^{a+b-1} {p} (eq-an-smatrix).
DG = Delius-Gandenberger hep-th/9904002, G98 = Gandenberger hep-th/9806003.

Checks:
 1. Block identities (including the doubling identity (x)(2 theta) = -(x/2)(x/2+h)) and the bulk S-matrix
    (unitarity, crossing) used below.
 2. Uniform (+...+) boundary, DG (3.16) (sec-uniform-boundary): iterating the bootstrap (3.15) from K_1 (3.12)
    gives (3.16); unitarity (eq-boundary-unitarity); crossing-unitarity in the book form
    eq-boundary-crossing-diagonal and in DG's form (3.17); K_a = K_{h-a}; the bulk boundary bootstrap
    eq-boundary-bootstrap-bulk at every fusing a+b -> c of a_n (both a+b < h and a+b > h); B -> 0 limit
    -1/((a)(h-a)) = eq-linear-reflection at C = +1 (and its closed form); equality with the conjecture of
    Corrigan-Dorey-Rietdijk-Sasaki hep-th/9404108 (5.16) with <x> replaced by <x~>; n = 1..7.
 3. Breather amplitude at imaginary coupling, DG (3.11): the closed form, continued by omega -> -2/B
    (the bulk continuation beta^2 -> -beta^2), is DG (3.12).  (The singlet projection that produces (3.11) from
    the soliton K involves the matrix part of S_VV and is not checked here; the product of the three scalar
    factors A_1 F_11 A_1 alone is not equal to (3.11).)  For n = 1, (3.11) is Ghoshal's K_B1 (eq-bsg-b1-reflection)
    at eta = pi/2, vartheta = 0, the point eps-hat = 0.
 4. n = 1: DG (3.16) is the sinh-Gordon factor eq-shg-reflection at E = B/2, F = 0 (the point eps-hat = 0,
    eta = pi/2, theta-param = 0 of Ghoshal's B_1 factor), DG (3.20) is E = 1 - B/2 (Neumann); E = 0 differs only
    at O(B^2).
 5. Neumann amplitudes DG (3.20) (sec-neumann-duality): unitarity, crossing-unitarity, bootstrap, K_a = K_{h-a},
    B -> 0 limit 1, (3.20) = (3.16) at B -> 2-B, neither is self-dual; O(beta^2) agreement with Kim
    hep-th/9506031 (33), (39) for a_3 and with G98 (6.11) for a_2; G98's K^(-) (from the s = -1 soliton K)
    satisfies the axioms, has classical limit 1, and is dual to the CDRS minimal conjecture (5.7) for C = -1.
 6. Solitonic boundaries DG (5.4)-(5.6) (sec-solitonic-boundaries): K_a^(b) = K_a M_a^(b) satisfies unitarity,
    crossing-unitarity and the bootstrap, M_a^(b) = M_a^(h-b), K_a^(b) = K_{h-a}^(b); the B -> 0 limit of
    M_a^(b) is the four-block product (5.7); DG (6.3), (6.18), (6.19), (6.25) are special cases, and so is
    (6.26) once the block (-9/2+B/2) missing from its second line is restored (as printed it violates
    crossing-unitarity);
    the CDD factors k_c (3.25) are unitary, satisfy the boundary CDD crossing condition and tend to 1.
 7. a_2 excited boundary states, DG (6.7) (sec-boundary-spectrum-a2): unitarity, crossing-unitarity with
    K_2^{b_nm} = K_1^{b_mn}, the bulk bootstrap on every state, the boundary bootstrap eq-boundary-bootstrap
    at the poles that create b_{n+1,m}, b_{n,m+1} and b_{-n-1,n+1}; (6.6), (6.9), (6.10) are special cases;
    energies E_{nm} from the bootstrap; pole table and spectrum at B = 9/20 (tbl-a2-boundary-states).
    a_3 (b = 2): DG (6.11) needs (-3+(2n1+1)B/2) in line (d) (printed with a minus sign) to reduce to (5.5)
    and to satisfy crossing-unitarity.
Runtime: about a minute.
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


# ---------------------------------------------------------------- 1. blocks and bulk S-matrix
ok = True
for h in (2, 3, 5):
    for x in (mp.mpf('0.3'), mp.mpf('1.7')):
        for th in THS:
            ok &= abs(blk(x, -th, h)*blk(x, th, h) - 1) < TOL and abs(blk(x + 2*h, th, h) - blk(x, th, h)) < TOL
            ok &= abs(blk(h, th, h) + 1) < TOL and abs(blk(x, I*PI - th, h) + blk(h - x, th, h)) < TOL
            ok &= abs(blk(x, 2*th, h) + blk(x/2, th, h)*blk(x/2 + h, th, h)) < TOL
report('blocks: (x)(-theta) = 1/(x), (x+2h) = (x), (h) = -1, (x)(i pi - theta) = -(h-x)(theta), '
       '(x)(2 theta) = -(x/2)(x/2+h)', ok)
ok = True
for n in (2, 3, 5):
    h = n + 1
    for a in range(1, h):
        for b in range(1, h):
            ok &= close(lambda t, B: S(a, b, n)(t, B)*S(a, b, n)(-t, B), lambda t, B: 1)
            ok &= close(lambda t, B: S(a, b, n)(I*PI - t, B), S(a, h - b, n))
            ok &= close(S(a, b, n), S(h - a, h - b, n))
report('a_n S-matrix (eq-an-smatrix): unitarity, crossing S_ab(i pi - theta) = S_{a,h-b}, C-invariance, n = 2, 3, 5', ok)

# ---------------------------------------------------------------- 2. uniform (+...+) boundary, DG (3.12)-(3.19)
ok_iter = True
for n in range(1, 8):
    h = n + 1
    K1 = Prod(h, {(-n, 0): 1, (n + 1, Fr(-1, 2)): 1, (-1, Fr(1, 2)): 1})          # DG (3.12)
    ok_iter &= close(K1, K_uniform(1, n))
    for a in range(2, h):                                                          # DG (3.15)
        f = lambda t, B, a=a: (K1(t + I*PI*(a - 1)/h, B) *
                               mp.fprod([K1(t - I*PI*(a + 1 - 2*c)/h, B)*S(a - c, 1, n)(2*t + I*PI*(3*c - a - 1)/h, B)
                                         for c in range(1, a)]))
        ok_iter &= close(f, K_uniform(a, n))
report('DG (3.15) iterated from K_1 = (3.12) gives the closed form (3.16), n = 1..7', ok_iter)
res = [axioms(lambda a, n=n: K_uniform(a, n), n, 'uniform') for n in range(1, 8)]
for i, nm in enumerate(['boundary unitarity', 'crossing-unitarity eq-boundary-crossing-diagonal',
                        'crossing-unitarity in DG form (3.17)', 'K_a = K_{h-a} (3.18)',
                        'bulk boundary bootstrap eq-boundary-bootstrap-bulk at every fusing']):
    report(f'DG (3.16): {nm}, n = 1..7', all(r[i] for r in res))
ok_cl = ok_lin = True
for n in range(1, 8):
    h = n + 1
    for a in range(1, h):
        for th in THS:
            Kcl = -1/(blk(a, th, h)*blk(h - a, th, h))
            ok_cl &= abs(K_uniform(a, n)(th, 0) - Kcl) < TOL
            s = mp.sin(PI*a/h)
            for C in (1, -1):
                ok_lin &= abs((mp.sinh(th) - I*C*s)/(mp.sinh(th) + I*C*s) + (blk(a, th, h)*blk(h - a, th, h))**(-C)) < TOL
report('DG (3.19): B -> 0 limit of (3.16) is -1/((a)(h-a)) = eq-linear-reflection at C = +1, n = 1..7', ok_cl)
report('eq-linear-reflection: (sinh th - iC sin(pi a/h))/(sinh th + iC sin(pi a/h)) = -[(a)(h-a)]^(-C), C = +-1',
       ok_lin)


def cdrs(a, n, tilde):          # CDRS (5.16); <x> = (x+1/2)/(x-1/2+B/2), <x~> = (x-1/2)/(x+1/2-B/2)
    h = n + 1

    def br(x, t, B):
        if tilde:
            return blk(x - mp.mpf(1)/2, t, h)/blk(x + mp.mpf(1)/2 - B/2, t, h)
        return blk(x + mp.mpf(1)/2, t, h)/blk(x - mp.mpf(1)/2 + B/2, t, h)
    return lambda t, B: mp.fprod([br(a - k - mp.mpf(1)/2, t, B)/br(h - a + k + mp.mpf(1)/2, t, B) for k in range(a)])


report('CDRS hep-th/9404108 (5.16) with <x> -> <x~> (their conjecture for A = +2, i.e. C = +1) equals DG (3.16), '
       'n = 1..7', all(close(cdrs(a, n, True), K_uniform(a, n)) for n in range(1, 8) for a in range(1, n + 1)))
report('... and CDRS (5.16) itself (their A = -2, C = -1 conjecture) is not (3.16) (control)',
       not close(cdrs(1, 2, False), K_uniform(1, 2)))
print('     physical-strip poles of DG (3.16) at B = 3/10 (theta = i pi p/h: order):',
      {n: {a: {str(p): o for p, o in K_uniform(a, n).poles(Fr(3, 10)).items() if o > 0} for a in range(1, (n + 3)//2)}
       for n in (2, 3, 4, 5)})
ok = all(all(p == Fr(k) and o == 1 for p, o in [(p, o) for p, o in K_uniform(a, n).poles(Fr(3, 10)).items() if o > 0]
             for k in [p]) and sorted(p for p, o in K_uniform(a, n).poles(Fr(3, 10)).items() if o > 0)
         == [Fr(k) for k in range(1, a) if Fr(k) < Fr(n + 1, 2)]
         for n in range(2, 8) for a in range(1, (n + 3)//2))
report('DG (6.1): the physical-strip poles of K_a (a <= h/2) are simple, at theta = i pi p/h, p = 1..a-1, '
       'and coupling-independent, n = 2..7', ok)

# ---------------------------------------------------------------- 3. breather bootstrap at imaginary coupling
def blkI(a, mu, h, w):
    """DG (3.4): (a)_I = sin(pi(mu+a)/(h w))/sin(pi(mu-a)/(h w)), w = DG's lambda = book omega."""
    return mp.sin(PI*(mu + a)/(h*w))/mp.sin(PI*(mu - a)/(h*w))


ok = True
for n in range(1, 8):
    h = n + 1
    for B in BS:
        wc = -2/B                                    # continuation beta^2 -> -beta^2: omega -> -2/B
        for th in THS:
            mu = -I*h*wc*th/(2*PI)
            v = blkI(-n*wc/2, mu, h, wc)*blkI(mp.mpf(1)/2 - h*wc/2, mu, h, wc)*blkI(-wc/2 - mp.mpf(1)/2, mu, h, wc)
            ok &= abs(v - K_uniform(1, n)(th, B)) < TOL
report('(a)_I -> (x) with x = 2a/omega and omega -> -2/B: the closed form (3.11) continues to DG (3.12), n = 1..7', ok)

lam_sg = mp.mpf('2.3')                     # n = 1, imaginary coupling: omega = sine-Gordon lambda, u = -i theta
err = 0
for th in THS:
    u, mu = -I*th, -I*2*lam_sg*th/(2*PI)
    dg = blkI(-lam_sg/2, mu, 2, lam_sg)*blkI(mp.mpf(1)/2 - lam_sg, mu, 2, lam_sg)*blkI(-lam_sg/2 - mp.mpf(1)/2, mu, 2, lam_sg)
    eta = PI/2                             # eps-hat = 0: (eta, vartheta) = (pi/2, 0)
    gh = (mp.cos(u/2 + PI/(4*lam_sg))*mp.cos(u/2 - PI/4 - PI/(4*lam_sg))*mp.sin(u/2 + PI/4)
          / (mp.cos(u/2 - PI/(4*lam_sg))*mp.cos(u/2 + PI/4 + PI/(4*lam_sg))*mp.sin(u/2 - PI/4))
          * (mp.cos(eta/lam_sg) - mp.sin(u))/(mp.cos(eta/lam_sg) + mp.sin(u)) * (1 - mp.sin(u))/(1 + mp.sin(u)))
    err = max(err, abs(dg - gh))
report(f'n = 1: DG (3.11) is Ghoshal\'s K_B1 (eq-bsg-b1-reflection) at eta = pi/2, vartheta = 0 [{mp.nstr(err, 2)}]',
       err < TOL)

# ---------------------------------------------------------------- 4. n = 1: Ghoshal / Corrigan-Delius
def K_shg(E, F):               # eq-shg-reflection, h = 2
    return lambda t, B: (blk(1, t, 2)*blk(2 - B/2, t, 2)*blk(1 + B/2, t, 2) /
                         (blk(1 - E(B), t, 2)*blk(1 + E(B), t, 2)*blk(1 - F(B), t, 2)*blk(1 + F(B), t, 2)))


zero = lambda B: 0
report('n = 1: DG (3.16) = eq-shg-reflection at E = B/2, F = 0 (eps-hat = 0: eta = pi/2, vartheta = 0)',
       close(K_uniform(1, 1), K_shg(lambda B: B/2, zero)))
report('n = 1: DG (3.20) = eq-shg-reflection at E = 1 - B/2, F = 0 (Neumann)',
       close(K_neumann(1, 1), K_shg(lambda B: 1 - B/2, zero)))
th = THS[0]
rat = [abs(K_shg(zero, zero)(th, B) - K_uniform(1, 1)(th, B))/B**2 for B in (mp.mpf('1e-3'), mp.mpf('1e-4'))]
report(f'n = 1: E = 0 (Chenaghlou-Corrigan at C_0 = C_1 = 1) and E = B/2 differ at O(B^2) only '
       f'[|diff|/B^2 = {mp.nstr(rat[0], 4)}, {mp.nstr(rat[1], 4)}]', abs(rat[0] - rat[1]) < 1e-2*rat[1])

# ---------------------------------------------------------------- 5. Neumann amplitudes DG (3.20), duality, O(beta^2)
res = [axioms(lambda a, n=n: K_neumann(a, n), n, 'Neumann') for n in range(1, 8)]
report('DG (3.20): unitarity, crossing-unitarity (both forms), K_a = K_{h-a}, bulk bootstrap at every fusing, n = 1..7',
       all(all(r) for r in res))
report('DG (3.20) = DG (3.16) with B -> 2-B, block by block, n = 1..7',
       all(close(K_neumann(a, n), K_uniform(a, n).dual()) for n in range(1, 8) for a in range(1, n + 1)))
report('DG (3.21): B -> 0 limit of (3.20) is 1, n = 1..7',
       all(abs(K_neumann(a, n)(th, 0) - 1) < TOL for n in range(1, 8) for a in range(1, n + 1) for th in THS))
report('neither (3.16) nor (3.20) is self-dual under B -> 2-B (n = 2, 3)',
       all(not close(K_uniform(a, n), K_uniform(a, n).dual()) for n in (2, 3) for a in range(1, n + 1)))
beta = mp.mpf('1e-6')
Bb = beta**2/(2*PI)/(1 + beta**2/(4*PI))
r2 = mp.sqrt(2)
kim = {1: lambda t: I/16*(mp.sinh(t)/(mp.cosh(t) + 1/r2) - mp.sinh(t)/(mp.cosh(t) - 1)),
       2: lambda t: I/16*(mp.sinh(t)/mp.cosh(t) - mp.sinh(t)/(mp.cosh(t) - 1) - mp.sinh(t)/(mp.cosh(t) - 1/r2)
                          + mp.sinh(t)/(mp.cosh(t) + 1/r2))}
err = max(abs((K_neumann(a, 3)(t, Bb) - 1)/beta**2 - kim[a](t)) for a in (1, 2) for t in THS)
report(f'a_3: O(beta^2) term of DG (3.20) = Kim hep-th/9506031 (33), (39) (Neumann, one loop) [{mp.nstr(err, 2)}]',
       err < 1e-8)
err_u = min(abs((K_uniform(1, 3)(t, Bb) - K_uniform(1, 3)(t, 0))/beta**2 - kim[1](t)) for t in THS)
report('... while the uniform amplitude (3.16) does not reproduce it (control)', err_u > 1e-3)
g98p = lambda t: -mp.mpf(1)/4*mp.cos(t/(2*I))/mp.sin(3*t/(2*I))
err = max(abs((K_neumann(1, 2)(t, Bb) - 1)/beta**2 - g98p(t)) for t in THS)
report(f'a_2: O(beta^2) term of DG (3.20) = G98 (6.11) for K^(+) [{mp.nstr(err, 2)}]', err < 1e-8)
Kminus = Prod(3, {(1, 0): 1, (-1, Fr(1, 2)): 1, (3, Fr(-1, 2)): 1}, sign=-1)       # G98 (6.4), s = -1
g98m = lambda t: -mp.mpf(1)/4*mp.sin(t/(2*I))/mp.cos(3*t/(2*I))
ok = (all(axioms(lambda a: Kminus, 2, 'K-')) and abs(Kminus(th, 0) - 1) < TOL
      and max(abs((Kminus(t, Bb) - 1)/beta**2 - g98m(t)) for t in THS) < 1e-8)
report('G98 K^(-) = -(1)(-1+B/2)(3-B/2) (s = -1): axioms with K_2 = K_1, B -> 0 limit 1, O(beta^2) term of G98 (6.11)',
       ok)
cdrs57 = Prod(3, {(1, 0): 1, (2, Fr(1, 2)): 1, (0, Fr(1, 2)): -1}, sign=-1)          # CDRS (5.7), C = -1
report('the dual of K^(-) is the CDRS minimal conjecture (5.7) -(1)(2+B/2)/(B/2) for C = -1', close(Kminus.dual(), cdrs57))

# ---------------------------------------------------------------- 6. solitonic boundaries DG (5.4)-(5.7)
ok_ax, ok_sym, ok_cl, ok_cl_swap = True, True, True, True
for n in range(2, 7):
    h = n + 1
    for b in range(1, h):
        r = axioms(lambda a, b=b, n=n: K_uniform(a, n)*M_sol(a, b, n), n, 'sol')
        ok_ax &= all(r)
        ok_sym &= all(close(M_sol(a, b, n), M_sol(a, h - b, n)) for a in range(1, h))
        for a in range(1, h):
            four = lambda t, a=a, b=b: (blk(a + b - mp.mpf(h)/2, t, h)*blk(-a - b - mp.mpf(h)/2, t, h)
                                        * blk(a - b + mp.mpf(h)/2, t, h)*blk(b - a + mp.mpf(h)/2, t, h))
            ok_cl &= all(abs(M_sol(a, b, n)(t, 0) - four(t)) < TOL for t in THS)
report('DG (5.5): K_a^(b) = K_a M_a^(b) satisfies unitarity, crossing-unitarity (both forms), K_a^(b) = K_{h-a}^(b) '
       'and the bulk bootstrap at every fusing, n = 2..6, all b', ok_ax)
report('DG (5.6): M_a^(b) = M_a^(h-b), n = 2..6', ok_sym)
report('DG (5.7): B -> 0 limit of M_a^(b) = (a+b-h/2)(-a-b-h/2)(a-b+h/2)(b-a+h/2), n = 2..6', ok_cl)
ok = close(K_uniform(1, 2)*M_sol(1, 1, 2),
           K_uniform(1, 2)*Prod(3, {(Fr(1, 2), Fr(-1, 2)): 1, (Fr(3, 2), Fr(-1, 2)): 1, (Fr(3, 2), Fr(1, 2)): 1,
                                    (Fr(5, 2), Fr(1, 2)): 1}))
p4 = lambda d: Prod(5, d)
ok &= close(M_sol(1, 1, 4), p4([((Fr(5, 2), Fr(-1, 2)), 1), ((Fr(5, 2), Fr(1, 2)), 1), ((Fr(-1, 2), Fr(-1, 2)), 1),
                                ((Fr(-9, 2), Fr(1, 2)), 1)]))                               # DG (6.18)
ok &= close(M_sol(2, 1, 4), p4({(Fr(1, 2), Fr(-1, 2)): 1, (Fr(9, 2), Fr(1, 2)): 1, (Fr(3, 2), Fr(-1, 2)): 1,
                                (Fr(7, 2), Fr(1, 2)): 1}))                                  # DG (6.19)
ok &= close(M_sol(1, 2, 4), M_sol(2, 1, 4))                                                 # DG (6.25)
ok &= close(M_sol(2, 2, 4), M_sol(2, 1, 4)*p4([((Fr(5, 2), Fr(-1, 2)), 1), ((Fr(5, 2), Fr(1, 2)), 1),
                                                 ((Fr(-1, 2), Fr(-1, 2)), 1), ((Fr(-9, 2), Fr(1, 2)), 1)]))  # (6.26)
report('DG (6.3) (a_2), (6.18), (6.19), (6.25) and (6.26) with the block (-9/2+B/2) restored (a_4) are DG (5.5)', ok)
K26 = lambda a: (K_uniform(a, 4)*M_sol(2, 1, 4)*p4([((Fr(5, 2), Fr(-1, 2)), 1), ((Fr(5, 2), Fr(1, 2)), 1),
                                                   ((Fr(-1, 2), Fr(-1, 2)), 1)]) if a in (2, 3) else K_uniform(a, 4)*M_sol(a, 2, 4))
report('DG (6.26) as printed (three blocks in its second line) violates crossing-unitarity (correction)',
       not axioms(K26, 4, '')[1])
print('     physical-strip poles of M_a^(b) at B = 3/10:',
      {f'a_{n},b={b}': {a: {str(p): o for p, o in M_sol(a, b, n).poles(Fr(3, 10)).items()} for a in range(1, (n + 3)//2)}
       for n, b in ((2, 1), (3, 1), (3, 2), (4, 1), (4, 2))})
ok = True
for n in range(2, 6):
    h = n + 1
    for c in range(1, n + 1):
        kc = Prod(h, [((Fr(h, 2) + c, Fr(1, 2)), 1), ((Fr(h, 2) - c, Fr(-1, 2)), 1),
                      ((-Fr(h, 2) + c, Fr(-1, 2)), 1), ((-Fr(h, 2) - c, Fr(1, 2)), 1)])   # list: blocks may coincide
        ok &= close(lambda t, B: kc(t, B)*kc(-t, B), lambda t, B: 1) and all(abs(kc(t, 0) - 1) < TOL for t in THS)
        ok &= close(lambda t, B: kc(I*PI/2 - t, B), lambda t, B: kc(I*PI/2 + t, B))
report('DG (3.25): particle CDD factors k_c are unitary, satisfy the boundary CDD condition '
       'k(i pi/2 - theta) = k(i pi/2 + theta) (so K_1 k_c obeys crossing-unitarity whenever K_1 does), '
       'and tend to 1 as B -> 0, n = 2..5', ok)

# ---------------------------------------------------------------- 7. a_2 excited boundary states DG (6.7)
def K1nm(nn, mm):
    """DG (6.7); each line is (x(n,m))(h - x(m,n)).  Built as a list: blocks may coincide."""
    return K_uniform(1, 2)*Prod(3, [((Fr(1, 2), -Fr(2*nn + 1, 2)), 1), ((Fr(5, 2), Fr(2*mm + 1, 2)), 1),
                                    ((Fr(1, 2), -Fr(2*nn - 1, 2)), 1), ((Fr(5, 2), Fr(2*mm - 1, 2)), 1),
                                    ((Fr(3, 2), Fr(2*nn - 1, 2)), 1), ((Fr(3, 2), -Fr(2*mm - 1, 2)), 1),
                                    ((-Fr(5, 2), Fr(2*nn + 1, 2)), 1), ((-Fr(1, 2), -Fr(2*mm + 1, 2)), 1)])


def Kex(nn, mm):
    return lambda a: K1nm(nn, mm) if a == 1 else K1nm(mm, nn)          # DG (6.8)


R = range(-3, 4)
ok_ax = all(all(axioms(Kex(nn, mm), 2, 'ex', conj=False)) for nn in R for mm in R)
report('DG (6.7): K^{b_nm} satisfies unitarity, crossing-unitarity (both forms, K_2^{b_nm} = K_1^{b_mn}) and the bulk '
       'bootstrap K_2 = K_1 S_11 K_1, -3 <= n, m <= 3', ok_ax)
report('DG (6.7) at n = m = 0 is DG (6.3)', close(K1nm(0, 0), K_uniform(1, 2)*M_sol(1, 1, 2)))
report('DG (6.7) at (n, m) = (1, 0) is DG (6.6)',
       close(K1nm(1, 0), K_uniform(1, 2)*Prod(3, [((Fr(1, 2), Fr(-1, 2)), 1), ((Fr(-1, 2), Fr(-1, 2)), 1),
                                                  ((Fr(3, 2), Fr(1, 2)), 2), ((Fr(5, 2), Fr(-1, 2)), 1),
                                                  ((Fr(5, 2), Fr(1, 2)), 1), ((Fr(1, 2), Fr(-3, 2)), 1),
                                                  ((Fr(-5, 2), Fr(3, 2)), 1)])))


def bbs(Kold, Knew, part, v):
    """eq-boundary-bootstrap: particle `part` bound at theta = i v; K_b^new = S_{b,part}(th-iv) S(th+iv) K_b^old."""
    return all(close(Knew(b), lambda t, B, b=b: S(b, part, 2)(t - I*v(B), B)*S(b, part, 2)(t + I*v(B), B)*Kold(b)(t, B))
               for b in (1, 2))


ok1 = all(bbs(Kex(nn, mm), Kex(nn + 1, mm), 1, lambda B, nn=nn: PI*(mp.mpf(1)/2 - (2*nn + 1)*B/2)/3)
          for nn in R for mm in R)
ok2 = all(bbs(Kex(nn, mm), Kex(nn, mm + 1), 2, lambda B, mm=mm: PI*(mp.mpf(1)/2 - (2*mm + 1)*B/2)/3)
          for nn in R for mm in R)
ok3 = all(bbs(Kex(-nn, nn), Kex(-nn - 1, nn + 1), 1, lambda B, nn=nn: PI*(mp.mpf(3)/2 - (2*nn + 1)*B/2)/3) for nn in R)
report('boundary bootstrap: particle 1 at i pi(1/2-(2n+1)B/2)/3 takes b_nm to b_{n+1,m}; particle 2 at '
       'i pi(1/2-(2m+1)B/2)/3 takes b_nm to b_{n,m+1}', ok1 and ok2)
report('boundary bootstrap: particle 1 at i pi(3/2-(2n+1)B/2)/3 takes b_{-n,n} to b_{-n-1,n+1} (DG fig. 12a)', ok3)
ok = all(close(K1nm(-nn, nn), K_uniform(1, 2)*Prod(3, {(Fr(1, 2), Fr(2*nn - 1, 2)): 1, (Fr(5, 2), Fr(2*nn + 1, 2)): 1,
                                                      (Fr(3, 2), -Fr(2*nn + 1, 2)): 1, (Fr(3, 2), -Fr(2*nn - 1, 2)): 1}))
         and close(K1nm(nn, -nn), K_uniform(1, 2)*Prod(3, {(Fr(1, 2), -Fr(2*nn + 1, 2)): 1,
                                                          (Fr(5, 2), -Fr(2*nn - 1, 2)): 1,
                                                          (Fr(3, 2), Fr(2*nn - 1, 2)): 1, (Fr(3, 2), Fr(2*nn + 1, 2)): 1}))
         for nn in R)
report('DG (6.9) and (6.10) are (6.7) at m = -n and at (n, -n)', ok)


def K1_a3(n1, n3, printed=False):          # DG (6.11), a_3, b = 2; line (d) as printed or with the crossing partner
    d2 = ((-3, -Fr(2*n1 + 1, 2)) if printed else (-3, Fr(2*n1 + 1, 2)))
    return K_uniform(1, 3)*Prod(4, [((1, -Fr(2*n1 + 1, 2)), 1), ((3, Fr(2*n3 + 1, 2)), 1), ((1, -Fr(2*n1 - 1, 2)), 1),
                                    ((3, Fr(2*n3 - 1, 2)), 1), ((1, Fr(2*n1 - 1, 2)), 1), ((3, -Fr(2*n3 - 1, 2)), 1),
                                    ((-1, -Fr(2*n3 + 1, 2)), 1), (d2, 1)])


def a3_cross(pr):
    return all(close(lambda t, B: K1_a3(n1, n3, pr)(I*PI/2 - t, B),
                     lambda t, B: S(1, 3, 3)(2*t, B)*K1_a3(n3, n1, pr)(I*PI/2 + t, B)) for n1 in (0, 1, -1) for n3 in (0, 2))


report('a_3, b = 2: DG (6.11) with line (d) = (-1-(2n3+1)B/2)(-3+(2n1+1)B/2) reduces to (5.5) at n1 = n3 = 0 and '
       'satisfies crossing-unitarity with K_3^{n1,n3} = K_1^{n3,n1}', close(K1_a3(0, 0), K_uniform(1, 3)*M_sol(1, 2, 3))
       and a3_cross(False))
report('... with line (d) as printed, (-3-(2n1+1)B/2), both fail (correction)',
       not close(K1_a3(0, 0, True), K_uniform(1, 3)*M_sol(1, 2, 3)) and not a3_cross(True))


def F(k, B):                   # energy (units of m_1) of k particles 1 bound in sequence, from the bootstrap
    d = PI*B/6
    return mp.sin(k*d)*mp.cos(PI/6 - k*d)/mp.sin(d)


def E(nn, mm, B):
    return F(nn, B) + F(mm, B)


Bv = mp.mpf('0.45')
ok = all(abs(E(nn + 1, mm, Bv) - E(nn, mm, Bv) - mp.cos(PI*(mp.mpf(1)/2 - (2*nn + 1)*Bv/2)/3)) < TOL for nn in R for mm in R)
ok &= all(abs(E(-nn - 1, nn + 1, Bv) - E(-nn, nn, Bv) - mp.cos(PI*(mp.mpf(3)/2 - (2*nn + 1)*Bv/2)/3)) < TOL for nn in R)
ok &= all(abs(E(nn, -nn, Bv) - mp.sin(nn*PI*Bv/6)**2/mp.sin(PI*Bv/6)) < TOL for nn in R)
report('energies e(b_nm) - e(b_00) = m_1 [F(n) + F(m)], F(k) = sin(k d) cos(pi/6 - k d)/sin d, d = pi B/6, agree with '
       'eq-boundary-bound-state-energy along all three bootstrap steps; e(b_{n,-n}) = m_1 sin^2(n d)/sin d', ok)
BF = Fr(9, 20)
states = sorted({(nn, mm) for nn in R for mm in R if nn + mm > 0 and abs(nn) < 1/(2*BF) + Fr(1, 2)
                 and abs(mm) < 1/(2*BF) + Fr(1, 2)} | {(nn, -nn) for nn in range(-3, 4) if abs(nn) < 3/(2*BF) + Fr(1, 2)},
                key=lambda s: E(*s, Bv))
print('     a_2 boundary states at B = 9/20 (DG sec. 6.2 ranges): state, (e - e_00)/m_1, physical-strip poles of K_1 '
      '(position p: theta = i pi p/3, order)')
for st in states:
    pl = {str(p): o for p, o in K1nm(*st).poles(BF).items() if o > 0}
    print(f'     b_{st[0]},{st[1]}: {mp.nstr(E(*st, Bv), 4)}  {pl}')
ok = True
for st in states:                                     # every state is reached from the ground state by allowed poles
    pass
for (nn, mm) in states:
    if (nn + 1, mm) in states:
        x = Fr(1, 2) - (2*nn + 1)*BF/2
        ok &= 0 < x < Fr(3, 2) and K1nm(nn, mm).poles(BF).get(x, 0) >= 1
    if (nn, mm + 1) in states:
        x = Fr(1, 2) - (2*mm + 1)*BF/2
        ok &= 0 < x < Fr(3, 2) and K1nm(mm, nn).poles(BF).get(x, 0) >= 1
report('at B = 9/20 every creation pole b_nm -> b_{n+1,m}, b_{n,m+1} between listed states is a physical-strip pole', ok)
