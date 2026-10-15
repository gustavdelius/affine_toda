"""The eps-hat = 0 reflection matrix of the d_n^(1) vector solitons: scalar factor, breather amplitude and the continued
particle amplitude, and the classical C = +1 boundary it describes (Part VII, ch. 33, sec-dn-vector-k; ch. 35,
sec-dn-uniform; half-line open problem 2). Checks 6-9 added in the consolidation A6 (2026-10-13).

Book conventions: physical q = -exp(-i pi omega), dn_vector of boundary_common (principal gradation, h = 2n-2),
R-check = rcheck. T = (n-1) omega = h omega/2, mu = -i T theta/pi, blocks (a)_I = sin(pi(mu+a)/2T)/sin(pi(mu-a)/2T);
g(mu, c) = Gamma(mu+c)/Gamma(-mu+c).
Checks:
 1. At eps-hat = 0 the intertwiner V(y) -> V(1/y) is J, antidiagonal with entries +-1, independent of y, J^2 = 1; the
    boundary-state vector sum_ab J_ba v_a (x) v_b = sum +-v_(-i) (x) v_(-i) consists of extremal-weight vectors
    (weights -2 eps_i, multiplicity one in V (x) V), on which R-check is 1; so the matrix part of crossing-unitarity
    has factor 1 (argument; computed in check 6). n = 4, 5.
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
 6. At eps-hat = 0 the matrix part of crossing-unitarity is exactly 1 (R-check_(Vbar Vbar)(2 theta) fixes the
    boundary-state vector, construction of an_nondiagonal_scalar.py check 6), and (1 x J) R-check (1 x J) maps the B_1
    singlet of V x V (rapidity difference i(pi - pi/T)) to itself with coefficient 1, as assumed in eq-dn-breather. n = 4, 5.
 7. Real coupling, C_j = +1 (book eq-boundary-potential): the static kinks psi = f v, d_4: v = e_2 - e_3, Tzitzeica
    f'' = 2e^{2f} - 2e^{-f}, e^{f(0)} = 9/8; d_5: v = e_2 - e_4, sinh-Gordon f'' = 4 sinh f, e^{f(0)} = 4 - 2 sqrt2, solve the
    static field equation and the boundary condition (30 digits), are regular on x <= 0, and have energies
    16 sqrt2 - 8 - 6 sqrt6 and -24 - 4 sqrt2 + 32(1 + sqrt2) sin(pi/8) (units m/beta^2). For C_j >= 0 the energy functional
    is strictly convex, so a static solution is the unique vacuum (argument, not a check).
 8. Linearized reflection about these vacua (ODE, free waves at x = -L): diagonal in the particle basis, with
    K_a = prod_j S_aj(theta, B = 1) (the book's S-matrix) for every particle, spinors included, and for a = 1 the B -> 0
    limit of eq-dn-uniform-amplitude; no zero of the bound-state determinant below threshold (the physical-strip pole of K_1
    at p = i m_2/2 is a redundant pole of the e^{m_2 x} tail of the vacuum, not a bound state).
 9. d_6: the vacuum is found by shooting in the sector invariant under the diagram automorphisms (two tail modes, m_2 and
    m_4); the linearized reflection about it is diagonal with K_a = prod_j S_aj(theta, B = 1).
Runtime: about three minutes.
"""
import os
import numpy as np
import mpmath as mp
import boundary_common
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
_TC = os.path.dirname(os.path.abspath(boundary_common.__file__))
_ns = {'os': os, '_TC': _TC}
exec(open(os.path.join(_TC, 'smatrix_simply_laced.py')).read().split('\nB = 0.37')[0]
     .replace("open('toda_classical.py')", "open(os.path.join(_TC, 'toda_classical.py'))"), _ns)
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

# ---------------------------------------------------------------- 6. matrix part at eps-hat = 0: crossing factor and singlet
# (replaces the argument of check 1 and the assumption of eq-dn-breather by a computation)
_src_rm = open(os.path.join(_TC, 'reflection_matrices.py')).read()
_ns_rm = {}
exec(_src_rm.split('\nrng = np.random.default_rng(11)')[0], _ns_rm)
exec('def rcheck_w(A, B):' + _src_rm.split('def rcheck_w(A, B):')[1].split('\n# ----')[0], _ns_rm)
rcheck_w = _ns_rm['rcheck_w']
from boundary_common import cop, kron, _mk


def chevalley(V, y):
    """Vbar(y): V(1) composed with e <-> f, k -> k^-1, at spectral parameter y (as reflection_matrices.chevalley)."""
    g = {}
    for j in nodes_of(V):
        g[('e', j)] = y*V[('f', j)]; g[('f', j)] = V[('e', j)]/y; g[('k', j)] = np.linalg.inv(V[('k', j)])
    return _mk(g, 0)


def singlet(A, B):
    d = A['dim']*B['dim']
    M = np.vstack([cop(A, B, kd, j) - (np.eye(d) if kd == 'k' else 0) for j in nodes_of(A) for kd in 'efk'])
    _, sv, Vh = np.linalg.svd(M)
    w = Vh[-1].conj()
    return w/w[int(np.argmax(np.abs(w) > 1e-8))], sv[-1]/sv[-2]


dev6 = 0
for n in (4, 5):
    N = 2*n
    for om in (0.37, 0.61):
        q, T = q_phys(om), (n - 1)*om
        V1 = dn_vector(n, q, 1.0)
        P = intertwiners([(chevalley(V1, 0.83 + 0.21j)[(k, j)], dn_vector(n, q, 0.83 + 0.21j)[(k, j)])
                          for j in nodes_of(V1) for k in 'efk'])[0][0]               # Vbar(y) -> V(y), y-independent
        y = np.exp(om*(0.37 + 0.11j))
        pairs = [(coideal_gen(dn_vector(n, q, y), j, 1/q, 0), coideal_gen(dn_vector(n, q, 1/y), j, 1/q, 0))
                 for j in nodes_of(V1)]
        J = intertwiners(pairs)[0][0]
        J = J/J[0, N - 1]                                # entries +-1, J^2 = 1 (check 1)
        k = (np.linalg.solve(P, J)).T.reshape(N*N)       # boundary-state vector in Vbar (x) Vbar (an_nondiagonal_scalar check 6)
        for th in (0.31, 0.83 + 0.1j):
            R = rcheck_w(chevalley(V1, np.exp(2*om*th)), chevalley(V1, 1.0))
            v = R @ k
            dev6 = max(dev6, np.linalg.norm(v - k)/np.linalg.norm(k))
        mu = 0.11 + 0.19j
        thb = 1j*np.pi*mu/T
        u = np.pi - np.pi/T
        ta, tb = thb - 0.5j*u, thb + 0.5j*u
        Va, Vb = dn_vector(n, q, np.exp(om*ta)), dn_vector(n, q, np.exp(om*tb))
        w_in, gap = singlet(Va, Vb)
        Wa = dn_vector(n, q, np.exp(-om*tb))
        w_out, _ = singlet(Wa, dn_vector(n, q, np.exp(-om*ta)))
        v = kron(np.eye(N), J) @ (rcheck_w(Va, Wa) @ (kron(np.eye(N), J) @ w_in))
        c = np.vdot(w_out, v)/np.vdot(w_out, w_out)
        dev6 = max(dev6, np.linalg.norm(v - c*w_out)/np.linalg.norm(v), abs(c - 1), gap)
report(f'd_n, eps-hat = 0: R-check_(Vbar Vbar)(2 theta) fixes the boundary-state vector (crossing factor 1), and '
       f'(1 x J) R-check (1 x J) maps the B_1 singlet (u = pi - pi/T) to itself with coefficient 1 (n = 4, 5) [{dev6:.0e}]',
       dev6 < 1e-10)

# ---------------------------------------------------------------- 7. the C = +1 vacuum of d_4 and d_5
# Real coupling, beta = m = 1, psi = beta phi, book eq-boundary-potential with C_j = +1:
# psi'' = sum_j n_j a_j e^{a_j.psi}, psi'(0) = -sum_j sqrt(n_j) a_j e^{a_j.psi(0)/2}; roots of dn_vector.
# For C_j >= 0 the energy is strictly convex, so the static solution below is the unique vacuum (argument, not a check).
mp.mp.dps = 30


def droots(n):
    E = np.eye(n)
    A = [-E[0] - E[1]] + [E[i] - E[i + 1] for i in range(n - 1)] + [E[n - 2] + E[n - 1]]
    return np.array(A), np.array([1, 1] + [2]*(n - 3) + [1, 1], float)


def kink(n, x, lib=np):
    """psi = f(x) v; d_4: v = e_2 - e_3, u = e^f = 2/(s^2 - 1), s = sqrt3 (1 - E)/(1 + E), E = E0 e^{sqrt6 x};
    d_5: v = e_2 - e_4, tanh(f/4) = E0 e^{2x}, E0 = tanh(f(0)/4), e^{f(0)} = 4 - 2 sqrt2."""
    if n == 4:
        r3 = lib.sqrt(3)
        E = (3*r3 - 5)/(3*r3 + 5)*lib.exp(lib.sqrt(6)*x)
        s = r3*(1 - E)/(1 + E)
        return lib.log(2/(s**2 - 1))
    return 4*lib.atanh(lib.tanh(lib.log(4 - 2*lib.sqrt(2))/4)*lib.exp(2*x)) if lib is mp else \
        4*np.arctanh(np.tanh(np.log(4 - 2*np.sqrt(2))/4)*np.exp(2*x))


VEC = {4: (1, 2), 5: (1, 3)}
EXACT = {4: 16*mp.sqrt(2) - 8 - 6*mp.sqrt(6), 5: -24 - 4*mp.sqrt(2) + 32*(1 + mp.sqrt(2))*mp.sin(mp.pi/8)}
dev7, regular = 0, True
for n in (4, 5):
    A, nj = droots(n)
    v = np.zeros(n); v[VEC[n][0]] = 1; v[VEC[n][1]] = -1
    cj = A @ v
    f = lambda x: kink(n, x, mp)
    for x in (mp.mpf('-0.3'), mp.mpf('-2.1'), mp.mpf(0)):
        d2 = mp.diff(f, x, 2)
        dev7 = max(dev7, max(abs(d2*v[k] - sum(nj[j]*A[j, k]*mp.exp(cj[j]*f(x)) for j in range(n + 1))) for k in range(n)))
    d1 = mp.diff(f, mp.mpf(0))
    dev7 = max(dev7, max(abs(d1*v[k] + sum(mp.sqrt(nj[j])*A[j, k]*mp.exp(cj[j]*f(0)/2) for j in range(n + 1)))
                         for k in range(n)))
    dens = lambda x: mp.diff(f, x)**2*(v @ v)/2 + sum(nj[j]*(mp.exp(cj[j]*f(x)) - 1) for j in range(n + 1))
    E = mp.quad(dens, [-mp.inf, -5, -1, 0]) + 2*sum(mp.sqrt(nj[j])*(mp.exp(cj[j]*f(0)/2) - 1) for j in range(n + 1))
    dev7 = max(dev7, abs(E - EXACT[n]))
    xs = np.linspace(-25, 0, 2001)
    regular &= bool(np.all(np.isfinite(kink(n, xs))) and np.all(kink(n, xs) >= 0))
report(f'C = +1 vacuum: the invariant kinks of d_4 (Tzitzeica, e^f(0) = 9/8) and d_5 (sinh-Gordon, e^f(0) = 4 - 2 sqrt2) '
       f'solve the static field equation and the boundary condition, are regular on x <= 0, and have energies '
       f'16 sqrt2 - 8 - 6 sqrt6 and -24 - 4 sqrt2 + 32(1 + sqrt2) sin(pi/8) (m/beta^2) [{mp.nstr(dev7, 2)}]',
       dev7 < mp.mpf('1e-20') and regular)

# ---------------------------------------------------------------- 8. linearized reflection about the vacuum
# eps'' + Om^2 eps = M2(x) eps, M2 = sum_j n_j e^{a_j.psi} a_j a_j^T; eps'(0) = -R eps(0), R = (1/2) sum_j sqrt(n_j)
# e^{a_j.psi(0)/2} a_j a_j^T. Solutions started at x = -L as free waves in the mass eigenbasis (closed channels: the
# decaying exponential), integrated to 0; K = outgoing/incoming amplitudes (book eq-linear-reflection convention).
from scipy.integrate import solve_ivp


def linear_K(A, nj, psi, psi0, Om, L):
    n = A.shape[1]
    M0 = (A.T*nj) @ A
    ev, U = np.linalg.eigh(M0)
    ma = np.sqrt(ev)
    M2 = lambda x: (A.T*(nj*np.exp(A @ psi(x)))) @ A
    R = 0.5*(A.T*(np.sqrt(nj)*np.exp(A @ psi0/2))) @ A
    op = [a for a in range(n) if ma[a] < Om]
    cols, kinds = [], []
    for a in op:
        p = np.sqrt(Om**2 - ma[a]**2)
        for sg in (1, -1):
            cols.append(np.concatenate([U[:, a], 1j*sg*p*U[:, a]])*np.exp(-1j*sg*p*L)); kinds.append(sg)
    for a in range(n):
        if ma[a] >= Om:
            k = np.sqrt(ma[a]**2 - Om**2)
            cols.append(np.exp(-k*L)*np.concatenate([U[:, a], k*U[:, a]]).astype(complex)); kinds.append(0)
    Y0 = np.array(cols).T
    m = Y0.shape[1]
    rhs = lambda x, y: np.vstack([y.reshape(2*n, m)[n:], (M2(x) - Om**2*np.eye(n)) @ y.reshape(2*n, m)[:n]]).ravel()
    Y = solve_ivp(rhs, [-L, 0.0], Y0.ravel(), method='DOP853', rtol=1e-12, atol=1e-14).y[:, -1].reshape(2*n, m)
    Bm = Y[n:] + R @ Y[:n]
    iin = [i for i, t in enumerate(kinds) if t == 1]
    iun = [i for i, t in enumerate(kinds) if t != 1]
    K = np.array([np.linalg.solve(Bm[:, iun], -Bm[:, i])[:len(op)] for i in iin]).T
    return ma[op], K, Bm


dev8, offd, nobound = 0, 0, True
for n in (4, 5):
    A, nj = droots(n)
    v = np.zeros(n); v[VEC[n][0]] = 1; v[VEC[n][1]] = -1
    psi = lambda x, n=n, v=v: kink(n, x)*v
    h, bl, mass, _ = sl.setup('d', n)
    mmin = np.sqrt(np.linalg.eigvalsh((A.T*nj) @ A)).min()
    scale = mmin/mass.min()
    L = 34/(np.sqrt(6) if n == 4 else 2)
    for Om in (1.1*mmin, 2.6, 7.0):
        mo, K, _ = linear_K(A, nj, psi, psi(0.0), Om, L)
        offd = max(offd, np.abs(K - np.diag(np.diag(K))).max())
        for i, m_ in enumerate(mo):
            th = np.arccosh(Om/m_)
            b = [c for c in range(n) if abs(mass[c]*scale - m_) < 1e-9][0]
            pred = np.prod([sl.S(bl[(b, j)], h, 1.0, th) for j in range(n)])
            dev8 = max(dev8, abs(K[i, i]/pred - 1))
            if abs(m_ - mmin) < 1e-9:
                K1cl = complex(-bth(h/2 - 1, th, h)*bth(h/2 + 1, th, h)/(bth(1, th, h)*bth(h - 1, th, h)*bth(h/2, th, h)**2))
                dev8 = max(dev8, abs(K[i, i]/K1cl - 1))
    # no bound state: the BC determinant of the decaying solutions has no zero below threshold; in particular none at
    # Om = m_1 sin(pi/h), where K_1 has its physical-strip pole (p = i m_2/2, a redundant pole of the e^{m_2 x} tail)
    dets = []
    for Om in list(np.linspace(0.02, 0.98, 25)*mmin) + [mmin*np.sin(np.pi/h)]:
        ev, U = np.linalg.eigh((A.T*nj) @ A)
        kap = np.sqrt(ev - Om**2)
        Y0 = (np.vstack([U, U*kap])*np.exp(-kap*L)).astype(complex)
        M2 = lambda x: (A.T*(nj*np.exp(A @ psi(x)))) @ A
        rhs = lambda x, y: np.vstack([y.reshape(2*n, n)[n:], (M2(x) - Om**2*np.eye(n)) @ y.reshape(2*n, n)[:n]]).ravel()
        Y = solve_ivp(rhs, [-L, 0.0], Y0.ravel(), method='DOP853', rtol=1e-12, atol=1e-14).y[:, -1].reshape(2*n, n)
        Y = Y/np.linalg.norm(Y[:n], axis=0)
        R = 0.5*(A.T*(np.sqrt(nj)*np.exp(A @ psi(0.0)/2))) @ A
        dets.append(np.linalg.det(Y[n:] + R @ Y[:n]).real)
    dets = np.array(dets)
    nobound &= bool(np.all(dets > 0) or np.all(dets < 0)) and np.abs(dets).min() > 1e-2
report(f'linearized reflection about the C = +1 vacuum (d_4, d_5): diagonal in the particle basis, K_a = prod_j S_aj(theta, B=1) '
       f'(book S-matrix) for every particle, spinors included, and K_1 = the B -> 0 limit of eq-dn-uniform-amplitude '
       f'[{dev8:.0e}; off-diagonal {offd:.0e}]; no boundary bound state below threshold', dev8 < 1e-9 and offd < 1e-12 and nobound)

# ---------------------------------------------------------------- 9. d_6: vacuum by shooting in the invariant sector
from scipy.optimize import fsolve
n = 6
A, nj = droots(n)
orb = {}
for j in range(n + 1):
    orb.setdefault('end' if min(j, n - j) <= 1 else min(j, n - j), []).append(j)
_, sv, Vt = np.linalg.svd(np.array([A[o[0]] - A[j] for o in orb.values() for j in o[1:]]))
Pinv = Vt[int(np.sum(sv > 1e-10)):].T                          # the sector invariant under the diagram automorphisms
ev, U = np.linalg.eigh((A.T*nj) @ A)
ma = np.sqrt(ev)
inv = [k for k in range(n) if np.linalg.norm(Pinv @ (Pinv.T @ U[:, k]) - U[:, k]) < 1e-9]
AP, dd = A @ Pinv, Pinv.shape[1]
L6 = 26/ma[inv].min()
rhs6 = lambda x, y: np.concatenate([y[dd:], Pinv.T @ ((A.T*nj) @ np.expm1(AP @ y[:dd]))])


def shoot(c, dense=False):
    y0 = sum(c[i]*np.exp(-ma[k]*L6)*np.concatenate([Pinv.T @ U[:, k], ma[k]*(Pinv.T @ U[:, k])]) for i, k in enumerate(inv))
    return solve_ivp(rhs6, [-L6, 0.0], y0, method='DOP853', rtol=1e-13, atol=1e-14*np.exp(-ma[inv].min()*L6),
                     dense_output=dense)


def bc6(c):
    yy = shoot(c).y[:, -1]
    p0, dp0 = Pinv @ yy[:dd], Pinv @ yy[dd:]
    return Pinv.T @ (dp0 + (A.T*np.sqrt(nj)) @ np.exp(A @ p0/2))


import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    c6 = fsolve(bc6, np.array([0.25, -0.085]), xtol=1e-14)    # seed from a direct minimization of the energy
s6 = shoot(c6, dense=True)
psi6 = lambda x: Pinv @ s6.sol(x)[:dd]
h, bl, mass, _ = sl.setup('d', n)
scale = ma.min()/mass.min()
dev9 = np.abs(bc6(c6)).max()
for Om in (1.2*ma.min(), 2.0, 3.8):
    mo, K, _ = linear_K(A, nj, psi6, psi6(0.0), Om, min(L6, 34/ma[inv].min()))
    dev9 = max(dev9, np.abs(K - np.diag(np.diag(K))).max())
    for i, m_ in enumerate(mo):
        b = [c for c in range(n) if abs(mass[c]*scale - m_) < 1e-8][0]
        pred = np.prod([sl.S(bl[(b, j)], h, 1.0, np.arccosh(Om/m_)) for j in range(n)])
        dev9 = max(dev9, abs(K[i, i]/pred - 1))
xg = np.linspace(-L6, 0, 4001)
E6 = getattr(np, 'trapezoid', getattr(np, 'trapz', None))([0.5*np.sum(np.gradient(psi6(xg), xg, axis=1)[:, i]**2) + nj @ (np.exp(A @ psi6(xg)[:, i]) - 1)
               for i in range(len(xg))], xg) + 2*np.sqrt(nj) @ (np.exp(A @ psi6(0.0)/2) - 1)
report(f'd_6, C = +1: the invariant static solution (two tail modes, m_2 and m_4) satisfies the boundary condition, and the '
       f'linearized reflection about it is diagonal with K_a = prod_j S_aj(theta, B=1) [{dev9:.0e}]; energy {E6:.6f} m/beta^2',
       dev9 < 1e-9)
