"""Checks added when the half-line consolidation round was refereed (Part VII: ch. 30 prp-soliton-preserving-half and
the half-solitons after it; ch. 31 prp-solitonic-linear; ch. 33 sec-an-nondiagonal-k (CDD factors) and
prp-dn-nondiagonal-scalar (the factor Xi); ch. 35 sec-neumann-duality, sec-dn-uniform, prp-dn-uniform-classical).
Run from code/.

Book conventions: blocks (x) = sinh((theta + i pi x/h)/2)/sinh((theta - i pi x/h)/2) (eq-blocks); imaginary-coupling blocks
(y)_I = sin(pi(mu+y)/2T)/sin(pi(mu-y)/2T), mu = -i T theta/pi, T = h omega/2, so (y)_I = (2y/omega); continuation
omega -> -2/B; P_x = (T/2 - x)_I (T/2 + x)_I.  Real coupling, units m = beta = 1: V = sum_j n_j (e^{alpha_j.phi} - 1),
B = 2 sum_j C_j sqrt(2 n_j/alpha_j^2)(e^{alpha_j.phi/2} - 1) (eq-boundary-potential).  The d_n S-matrix is the book's
(smatrix_simply_laced.setup, Dorey's block content); only its definitions are loaded.
Checks:
 1. Single-soliton charges of a_n, n = 2..6, by exact rational enumeration of the phase intervals of
    theta_j = 2 pi j b/h + Im xi (ch. 10, sec-topological-charges): h/gcd(b, h) charges for species b, those of species h-b
    are their negatives, and the nontrivial classes of weights mod 2 Lambda_W without a single-soliton charge number
    0, 2, 5, 21 for a_2..a_5; for a_3 they are the classes of e_1 + e_2 and e_1 + e_4.
 2. The block identity (z)(theta+iv)(z)(theta-iv) = (z + hv/pi)(z - hv/pi); the breather image P_x(mu_a) P_x(mu_b),
    mu_{a,b} = mu -+ (T/2 - 1/2), continued, equals D(chi) = (h+chi+B/2)(chi-B/2)(h-chi+B/2)(-chi-B/2), chi = 2x/omega
    (ch. 33); D = 1 at B = 0;
    every P_x is unitary and crossing-symmetric.
 3. d_n, n = 4..7: replacing P_0 by P_{1/2} and/or P_{omega/2} by P_{(omega+1)/2} in Xi multiplies K^N_1
    (eq-dn-neumann-amplitude) by D-ratios; the three alternatives are unitary, crossing-unitary with S_11, close the
    bootstrap of the vector-like particles (two fusing routes for n >= 6), equal K^N_1 with equal B-derivative at B = 0,
    and differ from it at B = 0.6.  For d_4 all four agree with Kim-Cho (hep-th/9505138) (29), (40) at O(beta^2),
    beta^2 = 2 pi B, and the amplitude without Xi does not.
 4. a_3: P_{(omega+1)/2} satisfies the soliton-bootstrap condition sigma(mu) = prod_{k=1}^n sigma(mu + T - k omega), as P_{omega/2}
    does; K^N_1 and the alternative both agree with Kim (hep-th/9506031) (33) at O(beta^2); so does K^(-)_1(theta + i pi),
    which differs from Kim's exact [1/2][3/2][5/2] (his (42)) at B = 0.6.
 5. d_5, particle 3: prod_b S_3b(theta, B = 1) has poles of order +1, -2, +1, -2 at theta = i pi z/8, z = 1..4 (block content,
    confirmed by epsilon scaling), a nonzero residue at 3 i pi/8, p = i(1 + sqrt2) m and Omega = m there; every other pole
    of every particle of d_4, d_5 with 0 < Im theta < pi/2 has p = i m_2/2.
 6. d_5 vacuum phi = f(e_2 - e_4): the spinor directions e_1, e_5 decouple from e_2, e_3, e_4 in the bulk Hessian
    M^2 = sum_j n_j e^{alpha_j.phi} alpha_j alpha_j^T and in Hess B.
 7. Bound states, d_4 and d_5 (prp-dn-uniform-classical): Y = eps' eps^-1 of the decaying solutions obeys
    Y' = M^2 - Omega^2 - Y^2, Y(-L) = U diag(sqrt(m_a^2 - Omega^2)) U^T; at 68 values of Omega^2 in [0, m_min^2), up to
    1e-10 m_min^2 below threshold, Y(0) + Hess B is symmetric with positive smallest eigenvalue, decreasing in Omega^2;
    at Omega = m (d_5) it is 0.897.
 8. Linearized reflection about regular continued vacua (eq-vacuum-pair), by direct ODE integration: a_2 (species 1)
    and a_4 (species 2): diagonal, equal to Delius and Gandenberger's K^(b)_c at B = 0 (prp-solitonic-amplitude (iii)).
Runtime: about half a minute.
"""
import os
from fractions import Fraction as Fr
from itertools import product
from math import gcd

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp

from boundary_common import report

HERE = os.path.dirname(os.path.abspath(__file__))
_cwd = os.getcwd()
os.chdir(HERE)                     # smatrix_simply_laced reads toda_classical.py from the cwd
_ns = {}
exec(compile(open('smatrix_simply_laced.py').read().split('\nB = 0.37')[0], 'smatrix_simply_laced.py', 'exec'), _ns)
os.chdir(_cwd)
setup = _ns['setup']
I = mp.mpc(0, 1)


# ---------------------------------------------------------------- 1. single-soliton charges of a_n
def species_charges(n, b):
    """Charges (0/1 vectors on e_1..e_h, mod the all-ones vector) of a single soliton of species b."""
    h = n + 1
    red = lambda x: (x - (x.numerator // x.denominator)) - (1 if x - (x.numerator // x.denominator) > Fr(1, 2) else 0)
    cuts = sorted({(Fr(1, 2) - Fr(j*b, h)) % 1 for j in range(h)})
    out = set()
    for i, c in enumerate(cuts):
        r = (c + (cuts[i + 1] if i + 1 < len(cuts) else cuts[0] + 1))/2
        t = [red(Fr(j*b, h) + r) for j in range(h)]          # [theta_j]/2 pi in (-1/2, 1/2]
        v = [t[i] - t[(i + 1) % h] for i in range(h)]          # -sum_j alpha_j t_j, alpha_j = e_j - e_{j+1}, alpha_0 = e_h - e_1
        v = [x - min(v) for x in v]
        assert all(x.denominator == 1 for x in v)
        out.add(tuple(int(x) for x in v))
    return out


ok1, info1 = True, {}
for n in range(2, 7):
    h = n + 1
    ch = {b: species_charges(n, b) for b in range(1, h)}
    for b in range(1, h):
        ok1 &= len(ch[b]) == h//gcd(b, h) and all(sorted(set(v)) == [0, 1] and sum(v) == b for v in ch[b])
        ok1 &= {tuple(1 - x for x in v) for v in ch[b]} == ch[h - b]       # species h-b: the negatives
    cls = lambda v: min(tuple(x % 2 for x in v), tuple(1 - x % 2 for x in v))
    reached = {cls(v) for b in ch for v in ch[b]}
    allc = {cls(S) for S in product((0, 1), repeat=h) if 0 < sum(S) < h}
    info1[n] = sorted(allc - reached)
counts = {n: len(m) for n, m in info1.items()}
ok1 &= [counts[n] for n in (2, 3, 4, 5)] == [0, 2, 5, 21]
ok1 &= set(info1[3]) == {cls((1, 1, 0, 0)), cls((1, 0, 0, 1))}
report(f'1. a_n single-soliton charges (exact): h/gcd(b,h) per species, species h-b = negatives; nontrivial classes without '
       f'a single-soliton charge: {counts} (a_2..a_6); a_3: classes of e1+e2 and e1+e4', ok1)


# ---------------------------------------------------------------- blocks
mp.mp.dps = 30


def blk(x, th, h):
    return mp.sinh((th + I*mp.pi*x/h)/2)/mp.sinh((th - I*mp.pi*x/h)/2)


def blkI(y, mu, T):
    return mp.sin(mp.pi*(mu + y)/(2*T))/mp.sin(mp.pi*(mu - y)/(2*T))


def D(q, B, th, h):
    return blk(h + q + B/2, th, h)*blk(q - B/2, th, h)*blk(h - q + B/2, th, h)*blk(-q - B/2, th, h)


def S_book(typ, n, a, b, th, B, _cache={}):
    key = (typ, n)
    if key not in _cache:
        _cache[key] = setup(typ, n)
    h, bl, _, _ = _cache[key]
    r = mp.mpc(1)
    for x, k in bl[(a, b)].items():
        r *= (blk(x - 1, th, h)*blk(x + 1, th, h)/(blk(x - 1 + B, th, h)*blk(x + 1 - B, th, h)))**k
    return r


PTS = [mp.mpc('0.47', '0.19'), mp.mpc('-0.93', '0.27')]
B0 = mp.mpf('0.6')

# ---------------------------------------------------------------- 2. block identity and the breather image of P_x
w2 = 0
for th in PTS:
    for z, v, h in ((mp.mpf('0.7'), mp.mpf('1.3'), 4), (mp.mpf('-2.2'), mp.mpf('0.4'), 6)):
        w2 = max(w2, abs(blk(z, th + I*v, h)*blk(z, th - I*v, h)/(blk(z + h*v/mp.pi, th, h)*blk(z - h*v/mp.pi, th, h)) - 1))
    for h in (3, 4, 6, 8):
        for B in (B0, mp.mpc('0.37', '0.2')):
            om = -2/B; T = h*om/2; mu = -I*T*th/mp.pi
            mua, mub = mu - (T/2 - mp.mpf(1)/2), mu + (T/2 - mp.mpf(1)/2)
            P = lambda x, m: blkI(T/2 - x, m, T)*blkI(T/2 + x, m, T)
            for x in (0, mp.mpf(1)/2, om/2, (om + 1)/2, (T - om)/2):
                q = 2*x/om
                w2 = max(w2, abs(P(x, mua)*P(x, mub)/D(q, B, th, h) - 1), abs(D(q, 0, th, h) - 1))
for T in (mp.sqrt(2)/2, mp.e - 1/mp.pi):
    for x in (mp.mpf('0.3'), mp.sqrt(3)/2):
        for mu in (mp.mpc('0.21', '0.13'), mp.mpc('-0.7', '0.4')):
            P = lambda m: blkI(T/2 - x, m, T)*blkI(T/2 + x, m, T)
            w2 = max(w2, abs(P(mu)*P(-mu) - 1), abs(P(T/2 - mu)/P(T/2 + mu) - 1))
report(f'2. (z)(theta+iv)(z)(theta-iv) = (z+hv/pi)(z-hv/pi); breather image of P_x = D(chi), D = 1 at B = 0; P_x unitary '
       f'and crossing-symmetric [{mp.nstr(w2, 2)}]', w2 < mp.mpf('1e-25'))


# ---------------------------------------------------------------- 3. d_n: Xi is not fixed by one loop
def KN1(th, h, B):
    b = lambda x: blk(x, th, h)
    return b(h/2 - 1)*b(h/2 + B/2)*b(h - 1 + B/2)/(b(h - 1)*b(h/2)*b(h/2 - 1 + B/2)*b(B/2))


def alt(th, h, B, flags):
    r = KN1(th, h, B)
    if flags[0]:
        r *= D(-B/2, B, th, h)/D(0, B, th, h)
    if flags[1]:
        r *= D(1 - B/2, B, th, h)/D(1, B, th, h)
    return r


def family(first, n, B):
    h = 2*n - 2
    fam = {1: first}
    for a in range(1, n - 2):
        old = fam[a]
        fam[a + 1] = (lambda t, a=a, old=old: first(t + I*mp.pi*a/h)*old(t - I*mp.pi/h)
                      * S_book('d', n, 0, a - 1, 2*t + I*mp.pi*(a - 1)/h, B))
    return fam


FLAGS = ((True, False), (False, True), (True, True))
ax, loop, dmin = 0, 0, mp.inf
for n in (4, 5, 6, 7):
    h = 2*n - 2
    for flags in FLAGS:
        fam = family(lambda t, flags=flags: alt(t, h, B0, flags), n, B0)
        for th in PTS:
            for a, Ka in fam.items():
                ax = max(ax, abs(Ka(th)*Ka(-th) - 1), abs(Ka(th)*Ka(th + I*mp.pi)/S_book('d', n, a - 1, a - 1, 2*th, B0) - 1))
            for a in range(4, n - 1):
                route = (fam[2](th + I*mp.pi*(a - 2)/h)*fam[a - 2](th - 2*I*mp.pi/h)
                         * S_book('d', n, 1, a - 3, 2*th + I*mp.pi*(a - 4)/h, B0))
                ax = max(ax, abs(route/fam[a](th) - 1))
            ratio = lambda B: alt(th, h, B, flags)/KN1(th, h, B)
            loop = max(loop, abs(ratio(0) - 1), abs(mp.diff(ratio, 0)))
            dmin = min(dmin, abs(ratio(B0) - 1))
sh, ch, s3 = mp.sinh, mp.cosh, mp.sqrt(3)
kc1 = lambda t: I*mp.pi/12*(sh(t)/ch(t) - sh(t)/(ch(t) - 1) + 2*sh(t)/(2*ch(t) + s3) - 2*sh(t)/(2*ch(t) - 1))
kc2 = lambda t: I*mp.pi/12*(sh(t)/ch(t) - sh(t)/(ch(t) - 1) - 4*sh(t)/(2*ch(t) - s3) + 2*sh(t)/(2*ch(t) + s3)
                            + 2*sh(t)/(2*ch(t) + 1))
kc, nomin = 0, mp.inf
for th in PTS:
    for flags in ((False, False),) + FLAGS:
        K1 = lambda B: alt(th, 6, B, flags)
        K2 = lambda B: alt(th + I*mp.pi/6, 6, B, flags)*alt(th - I*mp.pi/6, 6, B, flags)*S_book('d', 4, 0, 0, 2*th, B)
        kc = max(kc, abs(mp.diff(K1, 0) - kc1(th)), abs(mp.diff(K2, 0) - kc2(th)))
    Kmin = lambda B: KN1(th, 6, B)/(D(0, B, th, 6)*D(1, B, th, 6)*D(2, B, th, 6))
    nomin = min(nomin, abs(mp.diff(Kmin, 0) - kc1(th)))
report(f'3. d_n, n = 4..7: three alternatives to Xi are unitary, crossing-unitary, close the vector-like bootstrap '
       f'[{mp.nstr(ax, 2)}], equal K^N_1 to O(B) [{mp.nstr(loop, 2)}], differ at B = 0.6 [>= {mp.nstr(dmin, 2)}]; d_4: all '
       f'four agree with Kim-Cho (29), (40) at O(beta^2) [{mp.nstr(kc, 2)}], without Xi off by {mp.nstr(nomin, 2)}',
       ax < mp.mpf('1e-25') and loop < mp.mpf('1e-20') and dmin > mp.mpf('1e-4') and kc < mp.mpf('1e-20') and nomin > 0.1)

# ---------------------------------------------------------------- 4. a_3
n, h = 3, 4
om = mp.sqrt(2)/2; T = h*om/2
w4 = 0
for x in (om/2, (om + 1)/2):
    sig = lambda m: blkI(T/2 - x, m, T)*blkI(T/2 + x, m, T)
    for mu in (mp.mpc('0.19', '0.13'), mp.mpc('-0.4', '0.35')):
        w4 = max(w4, abs(sig(mu)/mp.fprod(sig(mu + T - k*om) for k in range(1, n + 1)) - 1))
kim33 = lambda t: I*mp.pi/8*(sh(t)/(ch(t) + 1/mp.sqrt(2)) - sh(t)/(ch(t) - 1))
KNa = lambda th, B: blk(-n, th, h)*blk(n + B/2, th, h)*blk(-B/2, th, h)
Kmi = lambda th, B: -blk(1, th, h)*blk(h - B/2, th, h)*blk(-1 + B/2, th, h)
half = mp.mpf(1)/2
kbr = lambda x, th, B: blk(x - half, th, h)*blk(x + half, th, h)/(blk(x - half + B/2, th, h)*blk(x + half - B/2, th, h))
kim42 = lambda th, B: kbr(half, th, B)*kbr(3*half, th, B)*kbr(5*half, th, B)
fin = mp.inf
for th in PTS:
    w4 = max(w4, abs(mp.diff(lambda B: KNa(th, B), 0) - kim33(th)),
             abs(mp.diff(lambda B: KNa(th, B)*D(1 - B/2, B, th, h)/D(1, B, th, h), 0) - kim33(th)),
             abs(mp.diff(lambda B: Kmi(th + I*mp.pi, B), 0) - kim33(th)), abs(mp.diff(lambda B: kim42(th, B), 0) - kim33(th)))
    fin = min(fin, abs(Kmi(th + I*mp.pi, B0)/kim42(th, B0) - 1))
report(f'4. a_3: P_((omega+1)/2) and P_(omega/2) satisfy the soliton-bootstrap condition; K^N_1, its alternative and '
       f'K^(-)_1(theta + i pi) agree with Kim (33) at O(beta^2) [{mp.nstr(w4, 2)}]; K^(-)_1(theta + i pi) differs from '
       f'Kim (42) at B = 0.6 by >= {mp.nstr(fin, 2)}', w4 < mp.mpf('1e-20') and fin > 0.05)


# ---------------------------------------------------------------- 5. poles of prod_b S_ab(theta, 1)
def prodS(n, a, th):
    return mp.fprod(S_book('d', n, a, b, th, 1) for b in range(n))


ok5, other = True, []
for n in (4, 5):
    h, bl, mass, _ = setup('d', n)
    scale = 2*np.sqrt(2)*np.sin(np.pi/h)/mass[0]          # index 0 is the vector, m_1 = 2 sqrt2 sin(pi/h) m
    m2 = mass[1]*scale
    for a in range(n):
        powers = {}
        for b in range(n):
            for x, k in bl[(a, b)].items():
                for y, p in ((x - 1, k), (x + 1, k), (x, -2*k)):
                    powers[y] = powers.get(y, 0) + p
        for z in range(1, h//2):
            order = powers.get(z, 0) - powers.get(2*h - z, 0)
            t0 = I*mp.pi*z/h
            v = [abs(prodS(n, a, t0 + e)) for e in (mp.mpf('1e-7'), mp.mpf('1e-9'))]
            ok5 &= int(round(float(mp.log(v[1]/v[0])/mp.log(100)))) == order
            if order > 0:
                p = mass[a]*scale*np.sin(np.pi*z/h)
                if abs(p - m2/2) > 1e-12:
                    other.append((n, a + 1, z, order, p, mass[a]*scale*np.cos(np.pi*z/h)))
        if n == 5 and a == 2:
            ok5 &= [powers.get(z, 0) - powers.get(2*h - z, 0) for z in (1, 2, 3, 4)] == [1, -2, 1, -2]
res = [e*prodS(5, 2, 3*I*mp.pi/8 + e) for e in (mp.mpf('1e-9'), mp.mpf('1e-11'))]
ok5 &= len(other) == 1 and other[0][:4] == (5, 3, 3, 1) and abs(other[0][4] - (1 + np.sqrt(2))) < 1e-12 \
    and abs(other[0][5] - 1) < 1e-12 and abs(res[1]) > 1e-3 and abs(res[1] - res[0]) < 1e-8
report(f'5. poles of prod_b S_ab(theta, 1) with 0 < Im theta < pi/2, d_4 and d_5 (orders confirmed by epsilon scaling): all at '
       f'p = i m_2/2 except d_5 particle 3 at 3 i pi/8 (order +1; p = {other[0][4]:.6f} i, Omega = {other[0][5]:.6f}, '
       f'residue {mp.nstr(res[1], 6)})', ok5)


# ---------------------------------------------------------------- 6., 7. d_4, d_5 vacua, Hessians, bound states
def dn_roots(n):
    e = np.eye(n)
    return (np.array([-e[0] - e[1]] + [e[i] - e[i + 1] for i in range(n - 1)] + [e[n - 2] + e[n - 1]]),
            np.array([1, 1] + [2]*(n - 3) + [1, 1]))


M2 = lambda A, nj, phi, x: (A.T*(nj*np.exp(A @ phi(x)))) @ A
HB = lambda A, nj, phi: 0.5*(A.T*(np.sqrt(nj)*np.exp(A @ phi(0.0)/2))) @ A
A5, n5 = dn_roots(5)
c5 = np.tanh(np.log(4 - 2*np.sqrt(2))/4)
phi5 = lambda x: 4*np.arctanh(c5*np.exp(2*min(x, 0.0)))*np.array([0, 1, 0, -1, 0.])   # tanh(f/4) = e^{2(x - x0)}, x <= 0
A4, n4 = dn_roots(4)
sol4 = solve_ivp(lambda s, w: [-w[0]*np.sqrt(2*(1 + w[0])*(3 + w[0]))], [0, 40], [1/8], method='DOP853',
                 rtol=1e-13, atol=1e-300, dense_output=True)                  # w = e^f - 1, s = -x
phi4 = lambda x: np.log1p(sol4.sol(-min(x, 0.0))[0])*np.array([0, 1, -1, 0.])
sp_, mid = [0, 4], [1, 2, 3]
coup = max(np.abs(M2(A5, n5, phi5, x)[np.ix_(sp_, mid)]).max() for x in np.linspace(-6, 0, 13))
coup = max(coup, np.abs(HB(A5, n5, phi5)[np.ix_(sp_, mid)]).max())
report(f'6. d_5 vacuum: the spinor directions e_1, e_5 decouple from e_2, e_3, e_4 in the bulk and boundary Hessians '
       f'[{coup:.0e}]', coup < 1e-14)


def riccati(A, nj, phi, z, L=20.0):
    n = A.shape[1]
    ev, U = np.linalg.eigh(M2(A, nj, lambda x: 0*phi(0.0), 0.0))
    Y0 = (U*np.sqrt(ev - z)) @ U.T
    rhs = lambda x, y: (M2(A, nj, phi, x) - z*np.eye(n) - y.reshape(n, n) @ y.reshape(n, n)).ravel()
    s = solve_ivp(rhs, [-L, 0.0], Y0.ravel(), method='DOP853', rtol=1e-11, atol=1e-13)
    J = s.y[:, -1].reshape(n, n) + HB(A, nj, phi)
    return s.success, np.abs(J - J.T).max(), np.linalg.eigvalsh((J + J.T)/2)[0]


ok7, rows = True, []
for name, A, nj, phi in (('d_4', A4, n4, phi4), ('d_5', A5, n5, phi5)):
    mm = np.linalg.eigvalsh(M2(A, nj, lambda x: 0*phi(0.0), 0.0)).min()
    zs = list(np.linspace(0, mm*(1 - 1e-3), 60)) + [mm*(1 - 10.0**-k) for k in range(4, 11)] + [1.0]
    out = [riccati(A, nj, phi, z) for z in zs]
    mins = [o[2] for o in out]
    ok7 &= all(o[0] and o[1] < 1e-8 for o in out) and min(mins) > 0.3 and all(np.diff(mins[:67]) < 0)
    rows.append(f'{name}: min {min(mins):.4f}, at Omega = m {mins[-1]:.4f}')
ok7 &= abs(mins[-1] - 0.897061) < 1e-5
report(f'7. no classical bound state below threshold (68 frequencies each): smallest eigenvalue of Y(0) + Hess B positive and '
       f'decreasing in Omega^2; ' + '; '.join(rows), ok7)


# ---------------------------------------------------------------- 8. regular continued vacua, direct ODE
def linear_K(n, a, chi, Om, L=22.0):
    h = n + 1
    al = np.array([np.eye(h)[j] - np.eye(h)[(j + 1) % h] for j in range(h)])
    s = np.sin(np.pi*a/h); ma = 2*s
    tau = lambda x: np.array([1 + 2*np.cos(chi + 2*np.pi*j*a/h)/(1 + s)*np.exp(ma*x) + (1 - s)/(1 + s)*np.exp(2*ma*x)
                              for j in range(h)])
    dtau0 = np.array([2*np.cos(chi + 2*np.pi*j*a/h)/(1 + s)*ma + (1 - s)/(1 + s)*2*ma for j in range(h)])
    phi = lambda x: -al.T @ np.log(tau(x))
    dphi0 = -al.T @ (dtau0/tau(0.0))
    C = min(product((1, -1), repeat=h),
            key=lambda C: np.abs(dphi0 + al.T @ (np.array(C)*np.exp(al @ phi(0.0)/2))).max())
    H = 0.5*(al.T*(np.array(C)*np.exp(al @ phi(0.0)/2))) @ al
    Q = np.linalg.qr(np.vstack([np.ones(h), np.eye(h)[:-1]]).T)[0][:, 1:]
    rhs = lambda x, y: np.concatenate([y[h:], ((al.T*np.exp(al @ phi(x))) @ al - Om**2*np.eye(h)) @ y[:h]])
    ys = [solve_ivp(rhs, [0.0, -L], np.concatenate([Q[:, k], -H @ Q[:, k]]).astype(complex), method='DOP853',
                    rtol=1e-12, atol=1e-14).y[:, -1] for k in range(n)]
    w = np.exp(2j*np.pi/h)
    rho = np.array([[w**(j*c) for j in range(h)] for c in range(1, h)]).T
    mc = 2*np.sin(np.pi*np.arange(1, h)/h)
    p = np.sqrt(Om**2 - mc**2)
    Mat = np.block([[rho*np.exp(-1j*p*L), rho*np.exp(1j*p*L)], [rho*1j*p*np.exp(-1j*p*L), -rho*1j*p*np.exp(1j*p*L)]])
    P = np.block([[Q.T, np.zeros((n, h))], [np.zeros((n, h)), Q.T]])
    co = np.array([np.linalg.solve(P @ Mat, P @ y) for y in ys]).T
    K = co[n:] @ np.linalg.inv(co[:n])
    th = np.arccosh(Om/mc)
    dg = [complex(-1/(blk(c, th[c - 1], h)*blk(h - c, th[c - 1], h))
                  * blk(c + a - h/2, th[c - 1], h)*blk(-c - a - h/2, th[c - 1], h)*blk(c - a + h/2, th[c - 1], h)
                  * blk(a - c + h/2, th[c - 1], h)) for c in range(1, h)]
    return np.abs(K - np.diag(np.diag(K))).max(), np.abs(np.diag(K)/np.array(dg) - 1).max()


r8 = [linear_K(2, 1, 0.4, 2.3), linear_K(4, 2, 0.5, 2.3)]
report(f'8. linearized reflection about regular continued vacua by direct ODE (a_2 species 1, a_4 species 2): diagonal '
       f'[{max(r[0] for r in r8):.0e}] and equal to DG K^(b)_c at B = 0 [{max(r[1] for r in r8):.0e}]',
       max(max(r) for r in r8) < 1e-9)
