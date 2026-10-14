"""B1(g), O.2, O.3: what one loop fixes of the CDD factors (ch. 33 prp-dn-nondiagonal-scalar, sec-an-nondiagonal-k;
ch. 35 eq-dn-neumann-amplitude, sec-neumann-duality).

Book blocks (eq-blocks): (x) = sinh((theta + i pi x/h)/2)/sinh((theta - i pi x/h)/2);
imaginary-coupling blocks (y)_I = sin(pi(mu+y)/2T)/sin(pi(mu-y)/2T), mu = -i T theta/pi, T = h omega/2, so (y)_I = (2y/omega).
P_x = (T/2 - x)_I (T/2 + x)_I; breather amplitude contains sigma(mu_a) sigma(mu_b), mu_{a,b} = mu -+ (T/2 - 1/2).
Checks:
 1. (z)(theta+iv)(z)(theta-iv) = (z + h v/pi)(z - h v/pi) at random points.
 2. P_x(mu_a) P_x(mu_b), continued (omega = -2/B), equals D(q;B) = (h+q+B/2)(q-B/2)(h-q+B/2)(-q-B/2), q = 2x/omega,
    for x = 0, 1/2, omega/2, (omega+1)/2, (T-omega)/2; and D(q;0) = 1.
 3. P_x is unitary and crossing-symmetric for any x: P(mu)P(-mu) = 1, P(T/2-mu) = P(T/2+mu).
 4. d_n, n = 4..7: K^N_1 (eq-dn-neumann-amplitude) times D(-B/2)/D(0) and/or D(1-B/2)/D(1):
    unitarity, crossing-unitarity with the book's S_11, bootstrap closure over the vector-like particles,
    equal value and B-derivative at B = 0 (one loop), different at B = 0.6.
 5. Kim-Cho (29), (40) at O(beta^2), beta^2 = 2 pi B + O(B^2): K^N_1, K^N_2 = K^N_1(th+i pi/6) K^N_1(th-i pi/6) S_11(2th)
    agree, and so do the alternatives; the minimal amplitude (without Xi) does not.
 6. a_3: P_{(omega+1)/2} obeys the soliton-bootstrap condition sigma(mu) = prod_{k=1}^n sigma(mu + T - k omega) (as does
    P_{omega/2}); K^N_1 and the alternative both agree with Kim (33) at O(beta^2); K^(-)_1(theta + i pi) agrees with Kim (33)
    at O(beta^2) but differs from Kim's exact [1/2][3/2][5/2] (42) at B = 0.6.
 7. rho_+- identities: rho_- = A^(-)/A_diag, rho_+ = (-T/2)_I^2 A^(+)/(A_diag P_{omega/2}) from the printed formulas (structural);
    (-T/2)_I^2 unitary and crossing-symmetric, (-T/2)_I alone crossing-antisymmetric.
"""
import os
import mpmath as mp

TC = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'textbook_code')
cwd = os.getcwd(); os.chdir(TC); ns = {}
exec(compile(open('smatrix_simply_laced.py').read().split('\nB = 0.37')[0], 'smatrix_simply_laced.py', 'exec'), ns)
os.chdir(cwd)
mp.mp.dps = 40
I = mp.mpc(0, 1)


def blk(x, th, h):
    return mp.sinh((th + I*mp.pi*x/h)/2)/mp.sinh((th - I*mp.pi*x/h)/2)


def blkI(y, mu, T):
    return mp.sin(mp.pi*(mu + y)/(2*T))/mp.sin(mp.pi*(mu - y)/(2*T))


def S_book(typ, n, a, b, th, B):
    h, bl, mass, _ = ns['setup'](typ, n)
    r = mp.mpc(1)
    for x, k in bl[(a, b)].items():
        r *= (blk(x-1, th, h)*blk(x+1, th, h)/(blk(x-1+B, th, h)*blk(x+1-B, th, h)))**k
    return r


pts = [mp.mpc('0.47', '0.19'), mp.mpc('-0.93', '0.27'), mp.mpc('1.31', '-0.41')]
worst = {}

# 1.
w = 0
for th in pts:
    for z, v, h in ((mp.mpf('0.7'), mp.mpf('1.3'), 4), (mp.mpf('-2.2'), mp.mpf('0.4'), 6)):
        w = max(w, abs(blk(z, th+I*v, h)*blk(z, th-I*v, h)/(blk(z+h*v/mp.pi, th, h)*blk(z-h*v/mp.pi, th, h)) - 1))
print(f'1. block identity: {mp.nstr(w, 3)}')

# 2.
def D(q, B, th, h):
    return blk(h+q+B/2, th, h)*blk(q-B/2, th, h)*blk(h-q+B/2, th, h)*blk(-q-B/2, th, h)


w2 = w0 = 0
for h in (3, 4, 6, 8):
    for B in (mp.mpf('0.6'), mp.mpc('0.37', '0.2')):
        om = -2/B
        T = h*om/2
        for th in pts:
            mu = -I*T*th/mp.pi
            mua, mub = mu - (T/2 - mp.mpf(1)/2), mu + (T/2 - mp.mpf(1)/2)
            P = lambda x, m: blkI(T/2 - x, m, T)*blkI(T/2 + x, m, T)
            for x in (0, mp.mpf(1)/2, om/2, (om+1)/2, (T-om)/2):
                q = 2*x/om
                w2 = max(w2, abs(P(x, mua)*P(x, mub)/D(q, B, th, h) - 1))
                w0 = max(w0, abs(D(q, 0, th, h) - 1))
print(f'2. breather image of P_x equals D(q;B): {mp.nstr(w2, 3)}; D(q;0) = 1: {mp.nstr(w0, 3)}')

# 3.
w3 = 0
for T in (mp.mpf('0.9'), mp.mpf('2.3')):
    for x in (mp.mpf('0.3'), mp.mpf('1.1')):
        for mu in (mp.mpc('0.21', '0.13'), mp.mpc('-0.7', '0.4')):
            P = lambda m: blkI(T/2 - x, m, T)*blkI(T/2 + x, m, T)
            w3 = max(w3, abs(P(mu)*P(-mu) - 1), abs(P(T/2 - mu)/P(T/2 + mu) - 1))
print(f'3. P_x unitary and crossing-symmetric: {mp.nstr(w3, 3)}')


# 4. d_n
def KN1(th, h, B):
    b = lambda x: blk(x, th, h)
    return b(h/2-1)*b(h/2+B/2)*b(h-1+B/2)/(b(h-1)*b(h/2)*b(h/2-1+B/2)*b(B/2))


def KU1(th, h, B):        # eq-dn-uniform-amplitude
    b = lambda x: blk(x, th, h)
    return b(h-B/2)*b(h/2-1)*b(h/2+1-B/2)/(b(h-1)*b(1-B/2)*b(h/2)*b(h/2-B/2))


def variant(th, h, B, flags):
    r = KN1(th, h, B)
    if flags[0]:
        r *= D(-B/2, B, th, h)/D(0, B, th, h)
    if flags[1]:
        r *= D(1-B/2, B, th, h)/D(1, B, th, h)
    return r


def family(first, n, B):
    """vector-like particles a = 1..n-2 by K_{a+1} = K_1(th + i pi a/h) K_a(th - i pi/h) S_{1a}(2 th + i pi(a-1)/h)."""
    h = 2*n - 2
    fam = {1: first}
    for a in range(1, n-2):
        old = fam[a]
        fam[a+1] = (lambda t, a=a, old=old: first(t + I*mp.pi*a/h)*old(t - I*mp.pi/h)
                    * S_book('d', n, 0, a-1, 2*t + I*mp.pi*(a-1)/h, B))
    return fam


w4ax, w4loop, dmin, dmax, w4dual = 0, 0, 10, 0, 0
B0 = mp.mpf('0.6')
for n in (4, 5, 6, 7):
    h = 2*n - 2
    for th in pts:
        w4dual = max(w4dual, abs(KN1(th, h, B0)/KU1(th, h, 2-B0) - 1))
    for flags in ((False, False), (True, False), (False, True), (True, True)):
        first = lambda t, flags=flags: variant(t, h, B0, flags)
        fam = family(first, n, B0)
        for th in pts:
            for a, Ka in fam.items():
                w4ax = max(w4ax, abs(Ka(th)*Ka(-th) - 1),
                           abs(Ka(th)*Ka(th + I*mp.pi)/S_book('d', n, a-1, a-1, 2*th, B0) - 1))
            for a in range(4, n-1):        # second fusing route 2 + (a-2) -> a
                alt = (fam[2](th + I*mp.pi*(a-2)/h)*fam[a-2](th - 2*I*mp.pi/h)
                       * S_book('d', n, 1, a-3, 2*th + I*mp.pi*(a-4)/h, B0))
                w4ax = max(w4ax, abs(alt/fam[a](th) - 1))
            if any(flags):
                ratio = lambda B: variant(th, h, B, flags)/KN1(th, h, B)
                w4loop = max(w4loop, abs(ratio(mp.mpf(0)) - 1), abs(mp.diff(ratio, 0)))
                d = abs(ratio(B0) - 1)
                dmin, dmax = min(dmin, d), max(dmax, d)
print(f'4. d_n K^N_1 = dual of uniform: {mp.nstr(w4dual, 3)}; axioms (unitarity, crossing-unitarity, bootstrap) for '
      f'K^N and the three alternatives: {mp.nstr(w4ax, 3)}; ratio and B-derivative at B=0: {mp.nstr(w4loop, 3)}; '
      f'|ratio - 1| at B = 0.6: {mp.nstr(dmin, 3)} .. {mp.nstr(dmax, 3)}')

# 5. Kim-Cho d_4 one loop: coefficient of B (beta^2 = 2 pi B): i pi/12 (...)
sh, ch, s3 = mp.sinh, mp.cosh, mp.sqrt(3)
kc1 = lambda t: I*mp.pi/12*(sh(t)/ch(t) - sh(t)/(ch(t)-1) + 2*sh(t)/(2*ch(t)+s3) - 2*sh(t)/(2*ch(t)-1))
kc2 = lambda t: I*mp.pi/12*(sh(t)/ch(t) - sh(t)/(ch(t)-1) - 4*sh(t)/(2*ch(t)-s3) + 2*sh(t)/(2*ch(t)+s3) + 2*sh(t)/(2*ch(t)+1))
h = 6
w5, w5min = 0, mp.inf
for th in pts:
    for flags in ((False, False), (True, False), (False, True), (True, True)):
        K1 = lambda B, th=th, flags=flags: variant(th, h, B, flags)
        K2 = lambda B, th=th, flags=flags: (variant(th + I*mp.pi/6, h, B, flags)*variant(th - I*mp.pi/6, h, B, flags)
                                            * S_book('d', 4, 0, 0, 2*th, B))
        w5 = max(w5, abs(K1(0) - 1), abs(K2(0) - 1), abs(mp.diff(K1, 0) - kc1(th)), abs(mp.diff(K2, 0) - kc2(th)))
    # minimal amplitude: remove the images of Xi = P_0 P_{omega/2} P_{(T-omega)/2}, q = 0, 1, h/2 - 1
    Kmin = lambda B, th=th: KN1(th, h, B)/(D(0, B, th, h)*D(1, B, th, h)*D(h/2 - 1, B, th, h))
    w5min = min(w5min, abs(mp.diff(Kmin, 0) - kc1(th)))
print(f'5. Kim-Cho (29), (40) at O(beta^2): K^N and alternatives {mp.nstr(w5, 3)}; minimal amplitude misses by >= {mp.nstr(w5min, 3)}')

# 6. a_3
n, h = 3, 4
w6 = 0
om = mp.mpf('0.61'); T = h*om/2
for x in (om/2, (om+1)/2):
    sig = lambda m: blkI(T/2 - x, m, T)*blkI(T/2 + x, m, T)
    for mu in (mp.mpc('0.19', '0.13'), mp.mpc('-0.4', '0.35')):
        w6 = max(w6, abs(sig(mu)/mp.fprod(sig(mu + T - k*om) for k in range(1, n+1)) - 1))
kim33 = lambda t: I*mp.pi/8*(sh(t)/(ch(t) + 1/mp.sqrt(2)) - sh(t)/(ch(t) - 1))
KNa = lambda th, B: blk(-n, th, h)*blk(n + B/2, th, h)*blk(-B/2, th, h)          # (-n)(n+B/2)(-B/2)
Kminus = lambda th, B: -blk(1, th, h)*blk(h - B/2, th, h)*blk(-1 + B/2, th, h)
kimbr = lambda x, th, B: blk(x - mp.mpf(1)/2, th, h)*blk(x + mp.mpf(1)/2, th, h)/(blk(x - mp.mpf(1)/2 + B/2, th, h)*blk(x + mp.mpf(1)/2 - B/2, th, h))
kim42 = lambda th, B: kimbr(mp.mpf(1)/2, th, B)*kimbr(mp.mpf(3)/2, th, B)*kimbr(mp.mpf(5)/2, th, B)
w6b, fin = 0, []
for th in pts:
    alt = lambda B, th=th: KNa(th, B)*D(1 - B/2, B, th, h)/D(1, B, th, h)
    w6b = max(w6b, abs(mp.diff(lambda B: KNa(th, B), 0) - kim33(th)), abs(mp.diff(alt, 0) - kim33(th)),
              abs(mp.diff(lambda B: Kminus(th + I*mp.pi, B), 0) - kim33(th)),
              abs(mp.diff(lambda B: kim42(th, B), 0) - kim33(th)))
    fin.append(abs(Kminus(th + I*mp.pi, B0)/kim42(th, B0) - 1))
print(f'6. a_3: CDD bootstrap condition for P_(omega/2), P_((omega+1)/2): {mp.nstr(w6, 3)}; one loop vs Kim (33) for K^N_1, '
      f'the alternative, K^(-)_1(theta+i pi) and Kim (42): {mp.nstr(w6b, 3)}; K^(-)_1(theta+i pi)/Kim(42) - 1 at B=0.6: '
      f'{[mp.nstr(f, 3) for f in fin]}')

# 7.
w7 = 0
for T in (mp.mpf('0.9'), mp.mpf('2.3')):
    for mu in (mp.mpc('0.21', '0.13'), mp.mpc('-0.7', '0.4')):
        c = lambda m: blkI(-T/2, m, T)
        w7 = max(w7, abs(c(mu)**2*c(-mu)**2 - 1), abs(c(T/2 - mu)**2/c(T/2 + mu)**2 - 1), abs(c(T/2 - mu)/c(T/2 + mu) + 1))
print(f'7. (-T/2)_I^2 unitary and crossing-symmetric, (-T/2)_I crossing-odd: {mp.nstr(w7, 3)}')
