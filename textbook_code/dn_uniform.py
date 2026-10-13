"""The eps-hat = 0 reflection matrix of the d_n^(1) vector solitons: scalar factor, breather amplitude and the continued
particle amplitude (Part VII, ch. 33, sec-dn-vector-k; ch. 35, sec-dn-uniform; half-line open problem 2).

Book conventions: physical q = -exp(-i pi omega), dn_vector of boundary_common (principal gradation, h = 2n-2),
R-check = rcheck. T = (n-1) omega = h omega/2, mu = -i T theta/pi, blocks (a)_I = sin(pi(mu+a)/2T)/sin(pi(mu-a)/2T);
g(mu, c) = Gamma(mu+c)/Gamma(-mu+c).
Checks:
 1. At eps-hat = 0 the intertwiner V(y) -> V(1/y) is J, antidiagonal with entries +-1, independent of y, J^2 = 1; the
    boundary-state vector sum_ab J_ba v_a (x) v_b = sum +-v_(-i) (x) v_(-i) consists of extremal-weight vectors
    (weights -2 eps_i, multiplicity one in V (x) V), on which R-check is 1; so the matrix part of crossing-unitarity
    has factor 1 (argument, not a check). n = 4, 5.
 2. The vector-vector scalar factor F (normalized to 1 on L(2 omega_1)) is unitary and satisfies
    F(T - mu) = F(mu) c(-mu), c(mu) = sin(pi mu) sin(pi(mu+T-omega))/(sin(pi(mu+omega)) sin(pi(mu+T))).
 3. The boundary scalar factor A solves A(mu)A(-mu) = 1 and A(T/2 - mu)/A(T/2 + mu) = F(2mu), with A(0) = 1 and no
    (-T/2)_I prefactor; the breather amplitude A(mu - T/2 + 1/2) F(2mu) A(mu + T/2 - 1/2) equals the closed form
    (omega/2 - T)_I (1/2 - T)_I (-omega/2 - 1/2)_I (T/2 - omega/2)_I (T/2 + omega/2 + 1/2)_I / ((T/2)_I (T/2 + 1/2)_I).
    n = 4, 5, 6; 30 digits (infinite products: exact head, Stirling tail summed with Hurwitz zeta).
 4. omega -> -2/B turns the breather amplitude into K_1 = (h-B/2)(h/2-1)(h/2+1-B/2)/((h-1)(1-B/2)(h/2)(h/2-B/2))
    (book blocks, eq-blocks), which satisfies unitarity and crossing-unitarity K_1(theta) K_1(theta + i pi) =
    S_11(2 theta) with the book's d_n S-matrix (smatrix_simply_laced), n = 4..7, and tends to the classical factor
    -(h/2-1)(h/2+1)/((1)(h-1)(h/2)^2) as B -> 0.
 5. K_1 equals Castro-Alvaredo and Fring's R^_1 = R~_1 prod_j S_1j(theta, 1 - B/2) (hep-th/0308004, (dn1), (Rhat)),
    with their block (blockw); and their spinor amplitude (dn2) satisfies their crossing-unitarity (cuR) with S_ss for
    n = 4..7, odd n included (a claimed failure for odd n came from swapping S_ss and S_ss').
Runtime: about 30 seconds.
"""
import numpy as np
import mpmath as mp
from boundary_common import dn_vector, q_phys, nodes_of, coideal_gen, intertwiners, report

OM = 0.37

# ---------------------------------------------------------------- 1. matrix part
dev1 = 0
for n in (4, 5):
    q = q_phys(OM)
    N = 2*n
    y = np.exp(OM*(0.41 + 0.23j))
    pairs = [(coideal_gen(dn_vector(n, q, y), j, 1/q, 0), coideal_gen(dn_vector(n, q, 1/y), j, 1/q, 0))
             for j in nodes_of(dn_vector(n, q, y))]
    sols, _ = intertwiners(pairs)
    assert len(sols) == 1
    J = sols[0]/sols[0][0, N - 1]
    anti = np.fliplr(np.diag(np.diag(np.fliplr(J))))
    dev1 = max(dev1, np.linalg.norm(J - anti), np.max(np.abs(np.abs(np.diag(np.fliplr(J))) - 1)),
               np.linalg.norm(J @ J - np.eye(N)))
    pairs2 = [(coideal_gen(dn_vector(n, q, 2*y), j, 1/q, 0), coideal_gen(dn_vector(n, q, 1/(2*y)), j, 1/q, 0))
              for j in nodes_of(dn_vector(n, q, y))]
    J2 = intertwiners(pairs2)[0][0]
    J2 = J2/J2[0, N - 1]
    dev1 = max(dev1, np.linalg.norm(J2 - J))
report(f'd_n, eps-hat = 0: K = J antidiagonal with entries +-1, independent of y, J^2 = 1 (n = 4, 5) [{dev1:.0e}]',
       dev1 < 1e-10)

# ---------------------------------------------------------------- 2./3. scalar factors
mp.mp.dps = 30


def logprod(P, terms, J=60, M=40):
    """log prod_{k>=0} prod_i Gamma(kP + a_i)^{c_i} (balanced): exact head, Stirling tail via Hurwitz zeta."""
    head = sum(c*mp.loggamma(k*P + a) for k in range(J) for c, a in terms)
    tail = sum((-1)**(m + 1)*sum(c*mp.bernpoly(m + 1, a) for c, a in terms)/(m*(m + 1)*P**m)*mp.zeta(m, J)
               for m in range(2, M + 1))
    return head + tail


def g(mu, c, sign=1):
    return [(sign, mu + c), (-sign, -mu + c)]


def Ffac(mu, T, om):
    """F = B_(1+T) B_(1+om)/(B_0 B_(T-om)), B_a = prod_k g(a + 2Tk)/g(a + T + 2Tk)."""
    terms = []
    for a, s in ((1 + T, 1), (1 + om, 1), (0, -1), (T - om, -1)):
        terms += g(mu, a, s) + g(mu, a + T, -s)
    return mp.exp(logprod(2*T, terms))


def Afac(mu, T, om):
    num, den = [1 + 3*T, 1 + om + 2*T, T, 2*T - om], [1 + 2*T, 1 + om + T, 2*T, 3*T - om]
    terms = []
    for a in num:
        terms += g(2*mu, a, 1)
    for a in den:
        terms += g(2*mu, a, -1)
    return mp.exp(logprod(4*T, terms))


def blk(a, mu, T):
    return mp.sin(mp.pi*(mu + a)/(2*T))/mp.sin(mp.pi*(mu - a)/(2*T))


dev2, dev3 = 0, 0
for n, om, mu in ((4, mp.mpf('0.37'), mp.mpc('0.13', '0.21')), (5, mp.mpf('0.53'), mp.mpc('-0.17', '0.11')),
                  (6, mp.mpf('0.29'), mp.mpc('0.21', '-0.07'))):
    T = (n - 1)*om
    s = lambda z: mp.sin(mp.pi*z)
    c = lambda m: s(m)*s(m + T - om)/(s(m + om)*s(m + T))
    dev2 = max(dev2, abs(Ffac(mu, T, om)*Ffac(-mu, T, om) - 1), abs(Ffac(T - mu, T, om)/(Ffac(mu, T, om)*c(-mu)) - 1))
    half = mp.mpf(1)/2
    lhs = Afac(mu - T/2 + half, T, om)*Ffac(2*mu, T, om)*Afac(mu + T/2 - half, T, om)
    rhs = (blk(om/2 - T, mu, T)*blk(half - T, mu, T)*blk(-om/2 - half, mu, T)*blk(T/2 - om/2, mu, T)
           * blk(T/2 + om/2 + half, mu, T)/(blk(T/2, mu, T)*blk(T/2 + half, mu, T)))
    dev3 = max(dev3, abs(Afac(mu, T, om)*Afac(-mu, T, om) - 1),
               abs(Afac(T/2 - mu, T, om)/Afac(T/2 + mu, T, om)/Ffac(2*mu, T, om) - 1), abs(Afac(0, T, om) - 1),
               abs(lhs/rhs - 1))
report(f'd_n vector scalar factor F: unitarity and crossing F(T - mu) = F(mu) c(-mu) (n = 4, 5, 6) '
       f'[{mp.nstr(dev2, 2)}]', dev2 < mp.mpf('1e-20'))
report(f'd_n boundary scalar factor A: unitarity, crossing A(T/2-mu)/A(T/2+mu) = F(2mu), A(0) = 1; breather amplitude '
       f'= closed form (n = 4, 5, 6) [{mp.nstr(dev3, 2)}]', dev3 < mp.mpf('1e-20'))

# ---------------------------------------------------------------- 4. particle level
class _SL:                                                   # the definitions of smatrix_simply_laced, without its checks
    pass


sl = _SL()
_ns = {}
exec(open(__file__.replace('dn_uniform.py', 'smatrix_simply_laced.py')).read().split('\nB = 0.37')[0], _ns)
sl.setup, sl.S = _ns['setup'], _ns['S']


def bth(x, t, h):
    return mp.sinh(t/2 + 1j*mp.pi*x/(2*h))/mp.sinh(t/2 - 1j*mp.pi*x/(2*h))


dev4 = 0
for n in (4, 5, 6, 7):
    h, bl, mass, _ = sl.setup('d', n)
    i1 = int(np.argmin(mass))                               # particle 1, the lightest
    K1 = lambda t, Bv: (bth(h - Bv/2, t, h)*bth(h/2 - 1, t, h)*bth(h/2 + 1 - Bv/2, t, h)
                        / (bth(h - 1, t, h)*bth(1 - Bv/2, t, h)*bth(h/2, t, h)*bth(h/2 - Bv/2, t, h)))
    for Bv in (mp.mpf('0.37'), mp.mpf('1.3')):
        for t in (mp.mpc('0.4', '0.1'), mp.mpc('-1.1', '0.3')):
            om = -2/Bv
            T = (n - 1)*om
            mu = -1j*T*t/mp.pi
            half = mp.mpf(1)/2
            cont = (blk(om/2 - T, mu, T)*blk(half - T, mu, T)*blk(-om/2 - half, mu, T)*blk(T/2 - om/2, mu, T)
                    * blk(T/2 + om/2 + half, mu, T)/(blk(T/2, mu, T)*blk(T/2 + half, mu, T)))
            dev4 = max(dev4, abs(cont/K1(t, Bv) - 1), abs(K1(t, Bv)*K1(-t, Bv) - 1),
                       abs(K1(t, Bv)*K1(t + 1j*mp.pi, Bv)/complex(sl.S(bl[(i1, i1)], h, float(Bv), complex(2*t))) - 1))
    cl = lambda t: -bth(h/2 - 1, t, h)*bth(h/2 + 1, t, h)/(bth(1, t, h)*bth(h - 1, t, h)*bth(h/2, t, h)**2)
    t = mp.mpc('0.6', '0.2')
    dev4 = max(dev4, abs(K1(t, mp.mpf('1e-12'))/cl(t) - 1))
report(f'd_n: continuation gives K_1, which satisfies unitarity and crossing-unitarity with the book d_n S-matrix and '
       f'tends to the classical factor (n = 4..7) [{mp.nstr(dev4, 2)}]', dev4 < 1e-10)

# ---------------------------------------------------------------- 5. Castro-Alvaredo and Fring
# Their block (blockw) ||x||_theta = <x,x>_theta/<x,x>_-theta, theta_h = i pi (2-B)/2h, theta_H = i pi B/2h.
# K_1 = R^_1 = R~_1 prod_j S_1j(theta, 1 - B/2) (their (dn1), (Rhat)); and their spinor amplitude (dn2) satisfies their
# crossing-unitarity R_s(theta) R_sbar(theta + i pi) = S_ss(2 theta) for even and odd n (for odd n, sbar = s').
def _ang(x, th, h, Bv):
    a1, a2 = 1j*mp.pi*(2 - Bv)/(2*h), 1j*mp.pi*Bv/(2*h)
    s = lambda z: mp.sinh(z/2)
    return (s(th + (x - 1)/2*a1 + (x - 1)/2*a2)*s(th + (x + 1)/2*a1 + (x + 1)/2*a2)
            / (s(th + (x - 1)/2*a1 + (x + 1)/2*a2)*s(th + (x + 1)/2*a1 + (x - 1)/2*a2)))


def W(x, th, h, Bv):
    return _ang(x, th, h, Bv)/_ang(x, -th, h, Bv)


dev5 = 0
for n in (4, 5, 6, 7):
    h, bl, mass, _ = sl.setup('d', n)
    i1, s1 = int(np.argmin(mass)), n - 2
    for Bv in (mp.mpf('0.37'), mp.mpf('1.3')):
        th = mp.mpc('0.31', '0.2')
        t = th - 1j*mp.pi
        Rhat = W(h - 1, t, h, Bv)*mp.fprod(W(2*k - 1, t, h, Bv) for k in range(1, h))
        Rhat *= mp.fprod(complex(sl.S(bl[(i1, j)], h, float(1 - Bv/2), complex(th))) for j in range(n))
        K1 = (bth(h - Bv/2, th, h)*bth(h/2 - 1, th, h)*bth(h/2 + 1 - Bv/2, th, h)
              / (bth(h - 1, th, h)*bth(1 - Bv/2, th, h)*bth(h/2, th, h)*bth(h/2 - Bv/2, th, h)))
        Rs = lambda tt: mp.fprod(W(2*k - 1, tt - 1j*mp.pi, h, Bv) for p in range(1, n//2 + 1)
                                 for k in range(2*p - 1, h - 2*p + 2))
        cu = Rs(th)*Rs(th + 1j*mp.pi)/complex(sl.S(bl[(s1, s1)], h, float(Bv), complex(2*th)))
        dev5 = max(dev5, abs(K1/Rhat - 1), abs(cu - 1))
report(f"K_1 = Castro-Alvaredo-Fring's R^_1, and their spinor amplitude satisfies their crossing-unitarity with S_ss "
       f"for n = 4..7 (no error for odd n) [{mp.nstr(dev5, 2)}]", dev5 < 1e-12)
