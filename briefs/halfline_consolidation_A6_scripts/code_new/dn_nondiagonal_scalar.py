"""Scalar factors of the non-diagonal d_n^(1) vector reflection matrices (prp-dn-k with all eps-hat_j^2 = eps_*^2), their
lowest-breather amplitudes and the continued particle amplitudes (Part VII, ch. 33, sec-dn-vector-k; ch. 35, sec-dn-uniform,
sec-boundary-open-problems, problem 2). Consolidation A6, 2026-10-13.

Book conventions: physical q = -exp(-i pi omega), dn_vector and the coideal of boundary_common, R-check from the weight-filtered
solver of reflection_matrices.py (rcheck_w). h = 2n-2, T = (n-1) omega, mu = -i T theta/pi, y = e^{omega theta},
w = (-q)^(1/2) = e^{-i pi omega/2}, blocks (x)_I = sin(pi(mu+x)/2T)/sin(pi(mu-x)/2T), book blocks (x) of eq-blocks.
eps-hat_0 = s eps_*, eps-hat_{j>=1} = eps_*, eps_* = 1/(2 cos(pi omega/2)) (analytic branch); varsigma = (-1)^n s.
K is normalized by its (1,1) entry: K_11 = -sin(pi T/2)/sin(pi mu) (varsigma = +1), -i sin(pi T/2)/cos(pi mu) (varsigma = -1).
Checks:
 1. Unitarity: K(-theta) K(theta) = U 1, U = varsigma (cos 2 pi mu - varsigma cos pi T)(cos 2 pi mu - varsigma cos pi omega)
    / (2 sin^2(pi omega/2) (cos 2 pi mu - varsigma)). n = 4, 5, 6.
 2. Crossing-unitarity of the matrix part (construction of an_nondiagonal_scalar.py check 6, with K: V -> V turned into
    V -> Vbar by the isomorphism P: Vbar(y) -> V(y)): k(theta) = Psi R-check_(Vbar Vbar)(2 theta) k(-theta), with the a_n
    functions Psi_+ = sin pi(omega/2 - mu)/sin pi(omega/2 + mu), Psi_- = cos pi(mu - omega/2)/cos pi(mu + omega/2), for
    varsigma = +-1. n = 4, 5.
 3. The B_1 singlet of V(theta_a) x V(theta_b), theta_b - theta_a = i(pi - pi/T), is mapped by (1 x K) R-check (1 x K) to
    c times itself, c = c^an_varsigma(mu) e_varsigma(mu): c^an_+ = -cos pi mu cos pi(mu + omega/2) cos pi(mu - T)
    /(sin^2(pi omega/2) cos pi(mu - omega/2)), e_+ = cos pi(mu - T/2 + omega/2) cos pi(mu - T/2 - omega/2)/cos^2 pi(mu - T/2);
    c^an_-, e_- the same with cos -> sin. n = 4, 5.
 4. The scalar factor A = A_0 rho g chi solves A(mu)A(-mu) U = 1 and A(T/2 - mu)/A(T/2 + mu) = F(2mu)/Psi(mu) (30 digits;
    n = 4, 5, 6). A_0, F: eq-dn-scalar and prp-dn-r. rho: the a_n factor of prp-an-nondiagonal-scalar without its CDD
    factor, at the d_n value of T; chi = (-T/2)_I for varsigma = +1, 1 for -1;
    g(mu) = v(mu) prod_{k>=1} [v(mu + kT)/v(-mu + kT)]^{(-1)^k}, v_+(x) = G(x + omega/2) G(1 + x - omega/2)/(G(x) G(1 + x)),
    v_-(x) = G(1/2 + x + omega/2) G(1/2 + x - omega/2)/G(1/2 + x)^2  (G = Gamma).
 5. Breather amplitude A(mu - T/2 + 1/2) F(2mu) A(mu + T/2 - 1/2) c(mu) = closed forms in blocks (n = 4, 5, 6; 30 digits).
 6. omega -> -2/B: closed forms K^(+-)_1 in book blocks, unitary, crossing-unitary with the book's S_11, bootstrap closure for
    all vector-like particles (n = 5, 6), and K -> 1 as B -> 0 (n = 4..7).
 7. One loop (B = beta^2/2 pi): d_4 against Kim and Cho, hep-th/9505138, eqs. (K1), (K2): varsigma = +1 agrees after the CDD
    factor (T/2)_I^2 (T/2 - omega/2)_I (T/2 + omega/2)_I (omega/2)_I (T - omega/2)_I, varsigma = -1 (no CDD) agrees at
    theta + i pi; the same against Kim's d_n conjecture (hep-th/9506031, fig. 7) for n = 5, 6, 7, which is
    Castro-Alvaredo-Fring's R~_1(theta + i pi) (hep-th/0308004, (dn1), (KimCF)); the a_3 analogue: K^(-)_1(theta) =
    Kim's (33) at theta + i pi.
 8. Semiclassics (beta^2 = 4 pi/(omega + 1) -> 0): both varsigma give (beta^2/4h) d(-i log A)/d theta -> half that of F(2 theta).
 9. With the CDD factor of check 7, the varsigma = +1 amplitude is the dual (B -> 2-B) of the uniform amplitude
    eq-dn-uniform-amplitude, K^N_1 = (h/2-1)(h/2+B/2)(h-1+B/2)/((h-1)(h/2)(h/2-1+B/2)(B/2)); the varsigma = -1 amplitude is
    K_1/K_1^cl, K_1^cl = -(h/2-1)(h/2+1)/((1)(h-1)(h/2)^2) (n = 4..7). Both relations are those of a_n (sec-neumann-duality).
Runtime: about three minutes.
"""
import os
import itertools
import numpy as np
import mpmath as mp
import boundary_common
from boundary_common import dn_vector, q_phys, nodes_of, coideal_gen, intertwiners, cop, kron, _mk, report

_TC = os.path.dirname(os.path.abspath(boundary_common.__file__))
_src = open(os.path.join(_TC, 'reflection_matrices.py')).read()
_ns = {}
exec(_src.split('\nrng = np.random.default_rng(11)')[0], _ns)
exec('def rcheck_w(A, B):' + _src.split('def rcheck_w(A, B):')[1].split('\n# ----')[0], _ns)
rcheck_w = _ns['rcheck_w']
_ns2 = {'os': os, '_TC': _TC}
exec(open(os.path.join(_TC, 'smatrix_simply_laced.py')).read().split('\nB = 0.37')[0]
     .replace("open('toda_classical.py')", "open(os.path.join(_TC, 'toda_classical.py'))"), _ns2)
setup = _ns2['setup']


# ---------------------------------------------------------------- the matrices
def chevalley(V, y):
    g = {}
    for j in nodes_of(V):
        g[('e', j)] = y*V[('f', j)]; g[('f', j)] = V[('e', j)]/y; g[('k', j)] = np.linalg.inv(V[('k', j)])
    return _mk(g, 0)


def Kmat(n, om, th, s):
    """intertwiner V(y) -> V(1/y) for eps-hat_0 = s eps_*, the others eps_*, normalized by its (1,1) entry."""
    q = q_phys(om)
    T = (n - 1)*om
    eps = np.ones(n + 1)/(2*np.cos(np.pi*om/2))
    eps[0] *= s
    y = np.exp(om*th)
    pairs = [(coideal_gen(dn_vector(n, q, y), j, 1/q, eps[j]), coideal_gen(dn_vector(n, q, 1/y), j, 1/q, eps[j]))
             for j in range(n + 1)]
    sols = intertwiners(pairs)[0]
    assert len(sols) == 1
    K = sols[0]
    mu = -1j*T*th/np.pi
    vs = (-1)**n*s
    k11 = -np.sin(np.pi*T/2)/np.sin(np.pi*mu) if vs == 1 else -1j*np.sin(np.pi*T/2)/np.cos(np.pi*mu)
    return K/K[0, 0]*k11


def U_pred(vs, mu, T, om):
    c2 = np.cos(2*np.pi*mu)
    return vs*(c2 - vs*np.cos(np.pi*T))*(c2 - vs*np.cos(np.pi*om))/(2*np.sin(np.pi*om/2)**2*(c2 - vs))


def Psi(vs, mu, om, lib=np):
    if vs == 1:
        return lib.sin(lib.pi*(om/2 - mu))/lib.sin(lib.pi*(om/2 + mu))
    return lib.cos(lib.pi*(mu - om/2))/lib.cos(lib.pi*(mu + om/2))


def csing(vs, mu, T, om, lib=np):
    tr = (lambda z: lib.cos(lib.pi*z)) if vs == 1 else (lambda z: lib.sin(lib.pi*z))
    can = -vs*tr(mu)*tr(mu + om/2)*tr(mu - T)/(lib.sin(lib.pi*om/2)**2*tr(mu - om/2))
    return can*tr(mu - T/2 + om/2)*tr(mu - T/2 - om/2)/tr(mu - T/2)**2


# 1. unitarity
d1 = 0
for n in (4, 5, 6):
    for om in (0.37, 0.61):
        T = (n - 1)*om
        for s in (1, -1):
            th = 0.37 - 0.21j
            P = Kmat(n, om, -th, s) @ Kmat(n, om, th, s)
            Up = U_pred((-1)**n*s, -1j*T*th/np.pi, T, om)
            d1 = max(d1, np.abs(P - Up*np.eye(2*n)).max()/abs(Up))
report(f'd_n, eps-hat_j^2 = eps_*^2: K(-theta) K(theta) = U_varsigma 1, varsigma = (-1)^n s (n = 4, 5, 6) [{d1:.0e}]', d1 < 1e-10)

# 2. crossing
d2 = 0
for n in (4, 5):
    N = 2*n
    for om in (0.37, 0.61):
        q, T = q_phys(om), (n - 1)*om
        V1 = dn_vector(n, q, 1.0)
        Pm = intertwiners([(chevalley(V1, 0.83 + 0.21j)[(k, j)], dn_vector(n, q, 0.83 + 0.21j)[(k, j)])
                           for j in nodes_of(V1) for k in 'efk'])[0][0]
        for s in (1, -1):
            for th in (0.31, -0.4 + 0.25j):
                km = np.linalg.solve(Pm, Kmat(n, om, 1j*np.pi/2 - th, s)).T.reshape(N*N)
                kp = np.linalg.solve(Pm, Kmat(n, om, 1j*np.pi/2 + th, s)).T.reshape(N*N)
                v = rcheck_w(chevalley(V1, np.exp(2*om*th)), chevalley(V1, 1.0)) @ kp
                f = np.vdot(v, km)/np.vdot(v, v)
                d2 = max(d2, np.linalg.norm(km - f*v)/np.linalg.norm(km),
                         abs(f/Psi((-1)**n*s, -1j*T*th/np.pi, om) - 1))
report(f'crossing-unitarity of the matrix part: k(theta) = Psi_varsigma R-check_(Vbar Vbar)(2 theta) k(-theta), Psi_+-'
       f' as for a_n (n = 4, 5) [{d2:.0e}]', d2 < 1e-10)


# 3. singlet projection
def singlet(A, B):
    d = A['dim']*B['dim']
    M = np.vstack([cop(A, B, kd, j) - (np.eye(d) if kd == 'k' else 0) for j in nodes_of(A) for kd in 'efk'])
    _, sv, Vh = np.linalg.svd(M)
    w = Vh[-1].conj()
    return w/w[int(np.argmax(np.abs(w) > 1e-8))]


d3 = 0
for n in (4, 5):
    N = 2*n
    for om, mu in ((0.53, 0.11 + 0.19j), (1.37, 0.05 - 0.12j)):
        q, T = q_phys(om), (n - 1)*om
        thb = 1j*np.pi*mu/T
        u = np.pi - np.pi/T
        ta, tb = thb - 0.5j*u, thb + 0.5j*u
        Va, Vb = dn_vector(n, q, np.exp(om*ta)), dn_vector(n, q, np.exp(om*tb))
        w_in = singlet(Va, Vb)
        Wa = dn_vector(n, q, np.exp(-om*tb))
        w_out = singlet(Wa, dn_vector(n, q, np.exp(-om*ta)))
        R = rcheck_w(Va, Wa)
        for s in (1, -1):
            v = kron(np.eye(N), Kmat(n, om, ta, s)) @ (R @ (kron(np.eye(N), Kmat(n, om, tb, s)) @ w_in))
            c = np.vdot(w_out, v)/np.vdot(w_out, w_out)
            d3 = max(d3, np.linalg.norm(v - c*w_out)/np.linalg.norm(v), abs(c/csing((-1)**n*s, mu, T, om) - 1))
report(f'B_1 singlet mapped to c_varsigma(mu) = c^an_varsigma(mu) e_varsigma(mu) times itself (n = 4, 5) [{d3:.0e}]', d3 < 1e-10)

# ---------------------------------------------------------------- 4./5. scalar factors, 30 digits
mp.mp.dps = 30


def logprod(P, terms, j0=0, J=60, M=40):
    head = sum(c*mp.loggamma(j*P + a) for j in range(j0, J) for c, a in terms)
    tail = sum((-1)**(k + 1)*sum(c*mp.bernpoly(k + 1, a) for c, a in terms)/(k*(k + 1)*P**k)*mp.zeta(k, J)
               for k in range(2, M + 1))
    return head + tail


def g2(mu, c, sign=1):
    return [(sign, mu + c), (-sign, -mu + c)]


def blk(x, mu, T):
    return mp.sin(mp.pi*(mu + x)/(2*T))/mp.sin(mp.pi*(mu - x)/(2*T))


def F(mu, T, om):
    terms = []
    for a, s in ((1 + T, 1), (1 + om, 1), (0, -1), (T - om, -1)):
        terms += g2(mu, a, s) + g2(mu, a + T, -s)
    return mp.exp(logprod(2*T, terms))


def A0(mu, T, om):
    terms = []
    for a in (1 + 3*T, 1 + om + 2*T, T, 2*T - om):
        terms += g2(2*mu, a, 1)
    for a in (1 + 2*T, 1 + om + T, 2*T, 3*T - om):
        terms += g2(2*mu, a, -1)
    return mp.exp(logprod(4*T, terms))


def rho(vs, mu, T, om):
    h = mp.mpf(1)/2
    if vs == 1:
        t = g2(mu, om/2 - 3*T/2) + g2(mu, 1 - om/2 - T/2) + g2(mu, om/2 - T/2, -1) + g2(mu, 1 - om/2 - 3*T/2, -1)
        return mp.sin(mp.pi*om/2)/mp.sin(mp.pi*(mu - T/2))*mp.exp(logprod(2*T, t, j0=1))*blk(-T/2, mu, T)
    t = g2(mu, h + om/2 - 3*T/2) + g2(mu, h - om/2 - T/2) + g2(mu, h + om/2 - T/2, -1) + g2(mu, h - om/2 - 3*T/2, -1)
    return -1j*mp.sin(mp.pi*om/2)/mp.cos(mp.pi*(mu - T/2))*mp.exp(logprod(2*T, t, j0=1))


def vterms(vs, x, om, sg=1):
    if vs == 1:
        return [(sg, x + om/2), (sg, 1 + x - om/2), (-sg, x), (-sg, 1 + x)]
    h = mp.mpf(1)/2
    return [(sg, h + x + om/2), (sg, h + x - om/2), (-2*sg, h + x)]


def g(vs, mu, T, om):
    t = vterms(vs, mu, om) + vterms(vs, -mu - T, om) + vterms(vs, -mu, om, -1) + vterms(vs, mu - T, om, -1)
    return mp.exp(sum(c*mp.loggamma(a) for c, a in vterms(vs, mu, om)) + logprod(2*T, t, j0=1))


def A(vs, mu, T, om):
    return A0(mu, T, om)*rho(vs, mu, T, om)*g(vs, mu, T, om)


def Ump(vs, mu, T, om):
    c2 = mp.cos(2*mp.pi*mu)
    return vs*(c2 - vs*mp.cos(mp.pi*T))*(c2 - vs*mp.cos(mp.pi*om))/(2*mp.sin(mp.pi*om/2)**2*(c2 - vs))


H = mp.mpf(1)/2


def KB_closed(vs, mu, T, om):
    b = lambda x: blk(x, mu, T)
    if vs == 1:
        return (b(om/2 - T)*b(H - T)**2*b(-om/2 - H)*b(T/2 - om/2)*b(T/2 + om/2 + H)*b(T/2 - H)
                / (b(T/2)*b(H)*b(H - om/2)*b(T + om/2 - H)*b(T/2 - H + om/2)*b(T/2 + H - om/2)))
    return -b(om/2)*b(H - T)*b(-om/2 - H)*b(T/2 + om/2 + H)*b(T/2)/(b(T/2 + H)*b(T/2 + om/2))


d4, d5 = 0, 0
for n, om, mu in ((4, '0.37', ('0.13', '0.21')), (4, '0.61', ('-0.17', '0.11')), (5, '0.53', ('0.21', '-0.07')),
                  (6, '0.29', ('0.08', '0.15'))):
    om, mu = mp.mpf(om), mp.mpc(*mu)
    T = (n - 1)*om
    for vs in (1, -1):
        Af = lambda m: A(vs, m, T, om)
        d4 = max(d4, abs(Af(mu)*Af(-mu)*Ump(vs, mu, T, om) - 1),
                 abs(Af(T/2 - mu)/Af(T/2 + mu)*Psi(vs, mu, om, mp)/F(2*mu, T, om) - 1))
        kb = Af(mu - T/2 + H)*F(2*mu, T, om)*Af(mu + T/2 - H)*csing(vs, mu, T, om, mp)
        d5 = max(d5, abs(kb/KB_closed(vs, mu, T, om) - 1))
report(f'scalar factors A_+-: unitarity A(mu)A(-mu)U = 1 and crossing A(T/2-mu)/A(T/2+mu) = F(2mu)/Psi (n = 4, 5, 6) '
       f'[{mp.nstr(d4, 2)}]', d4 < mp.mpf('1e-20'))
report(f'breather amplitudes: A(mu_a) F(2mu) A(mu_b) c(mu) = closed block forms (n = 4, 5, 6) [{mp.nstr(d5, 2)}]',
       d5 < mp.mpf('1e-20'))


# ---------------------------------------------------------------- 6. particle level
def bb(x, th, h):
    return mp.sinh(th/2 + 1j*mp.pi*x/(2*h))/mp.sinh(th/2 - 1j*mp.pi*x/(2*h))


def K1(vs, th, B, n, cdd=False):
    """continued breather amplitude (blocks of 5 at omega = -2/B), optionally with the CDD factor of check 7."""
    om = -2/B
    T = (n - 1)*om
    mu = -1j*T*th/mp.pi
    k = KB_closed(vs, mu, T, om)
    if cdd:
        c = T/2 - H
        for x in (0, om/2, (T - om)/2):
            for yy in (T/2 - x, T/2 + x):
                k *= blk(yy, mu - c, T)*blk(yy, mu + c, T)
    return k


def K1_closed(vs, th, B, n):
    h = 2*n - 2
    b = lambda x: bb(x, th, h)
    if vs == -1:
        return -b(1)*b(h/2)*b(h/2 + 1 - B/2)/(b(1 - B/2)*b(h + B/2)*b(h/2 - B/2)*b(h/2 + 1))
    return (b(h/2 - 1)*b(h/2 + 1 - B/2)*b(h/2 + B/2)*b(1 + B/2)*b(B/2)*b(h - 1 - B/2)
            / (b(h - 1)*b(h + B/2)**2*b(1 - B/2)*b(h/2)*b(h/2 + 1 + B/2)*b(h/2 - 1 - B/2)))


def Smp(bl, h, B, th):
    r = mp.mpf(1)
    for x, k in bl.items():
        r *= (bb(x - 1, th, h)*bb(x + 1, th, h)/(bb(x - 1 + B, th, h)*bb(x + 1 - B, th, h)))**k
    return r


def boot(n, Kfun, B):
    """K_c(th) = K_a(th + i al) K_b(th - i be) S_ab(2th + i(al - be)) at every forward-channel fusing of vector-like
    particles (poles of the book S-matrix with -i Res > 0 at the mass-triangle angle); returns the worst deviation."""
    h, bl, mass, _ = setup('d', n)
    fus = []
    for a, b, c in itertools.product(range(n - 2), repeat=3):
        ma, mb, mc = mass[a], mass[b], mass[c]
        if not abs(ma - mb) < mc < ma + mb:
            continue
        x = mp.acos((mc**2 - ma**2 - mb**2)/(2*ma*mb))*h/mp.pi
        if abs(x - mp.nint(x)) > 1e-9:
            continue
        u = mp.pi*mp.nint(x)/h
        r = mp.mpf('1e-12')*Smp(bl[(a, b)], h, mp.mpf('0.37'), 1j*u + mp.mpf('1e-12'))
        if not (1e-6 < abs(r) < 1e6 and mp.re(-1j*r) > 0 and abs(mp.im(-1j*r)) < 1e-6*abs(r)):
            continue
        al = mp.atan2(mb*mp.sin(u), ma + mb*mp.cos(u))
        fus.append((a, b, c, al, u - al))
    K = {0: Kfun}
    for _ in range(n):
        for a, b, c, al, be in fus:
            if a in K and b in K and c not in K:
                K[c] = (lambda t, a=a, b=b, al=al, be=be: K[a](t + 1j*al)*K[b](t - 1j*be)
                        * Smp(bl[(a, b)], h, B, 2*t + 1j*(al - be)))
    dev = 0 if len(K) == n - 2 else 1
    for t in (mp.mpc('0.31', '0.17'), mp.mpc('-0.7', '0.05')):
        for a, b, c, al, be in fus:
            dev = max(dev, abs(K[c](t)/(K[a](t + 1j*al)*K[b](t - 1j*be)*Smp(bl[(a, b)], h, B, 2*t + 1j*(al - be))) - 1))
        for a in K:
            dev = max(dev, abs(K[a](t)*K[a](t + 1j*mp.pi)/Smp(bl[(a, a)], h, B, 2*t) - 1), abs(K[a](t)*K[a](-t) - 1))
    return dev


d6 = 0
for n in (4, 5, 6, 7):
    h, bl, _, _ = setup('d', n)
    for vs in (1, -1):
        for Bv in (mp.mpf('0.37'), mp.mpf('1.3')):
            for th in (mp.mpc('0.4', '0.1'), mp.mpc('-1.1', '0.3')):
                d6 = max(d6, abs(K1(vs, th, Bv, n)/K1_closed(vs, th, Bv, n) - 1), abs(K1(vs, th, Bv, n)*K1(vs, -th, Bv, n) - 1),
                         abs(K1(vs, th, Bv, n)*K1(vs, th + 1j*mp.pi, Bv, n)/Smp(bl[(0, 0)], h, Bv, 2*th) - 1))
        d6 = max(d6, abs(K1(vs, mp.mpc('0.6', '0.2'), mp.mpf('1e-15'), n) - 1))
        if n in (5, 6):
            for cdd in (False, True):
                d6 = max(d6, boot(n, lambda t: K1(vs, t, mp.mpf('0.43'), n, cdd), mp.mpf('0.43')))
report(f'omega -> -2/B: K^(+-)_1 = closed forms, unitary, crossing-unitary with the book S_11 (n = 4..7), bootstrap closure '
       f'over the vector-like particles (n = 5, 6), K -> 1 as B -> 0 [{mp.nstr(d6, 2)}]', d6 < 1e-12)

# ---------------------------------------------------------------- 7. one loop
mp.mp.dps = 40
e = mp.mpf('1e-15')
dB = lambda f, t: (f(e, t) - f(-e, t))/(2*e)/(2*mp.pi)          # d/d(beta^2) at 0, B = beta^2/2pi + O(beta^4)
ths = [mp.mpf('0.37'), mp.mpc('1.2', '0.3'), mp.mpc('-0.8', '0.1')]
sh, ch, r3 = mp.sinh, mp.cosh, mp.sqrt(3)
kc1 = lambda t: 1j/24*(sh(t)/ch(t) - sh(t)/(ch(t) - 1) + 2*sh(t)/(2*ch(t) + r3) - 2*sh(t)/(2*ch(t) - 1))
kc2 = lambda t: 1j/24*(sh(t)/ch(t) - sh(t)/(ch(t) - 1) - 4*sh(t)/(2*ch(t) - r3) + 2*sh(t)/(2*ch(t) + r3)
                       + 2*sh(t)/(2*ch(t) + 1))
h4, bl4, _, _ = setup('d', 4)
K2 = lambda vs, cdd: (lambda B, t: K1(vs, t + 1j*mp.pi/6, B, 4, cdd)*K1(vs, t - 1j*mp.pi/6, B, 4, cdd)*Smp(bl4[(0, 0)], 6, B, 2*t))
ok_p = max(max(abs(dB(lambda B, t: K1(1, t, B, 4, True), t) - kc1(t)), abs(dB(K2(1, True), t) - kc2(t))) for t in ths)
ok_m = max(max(abs(dB(lambda B, t: K1(-1, t, B, 4), t) - kc1(t + 1j*mp.pi)), abs(dB(K2(-1, False), t) - kc2(t + 1j*mp.pi)))
           for t in ths)
no_p = min(abs(dB(lambda B, t: K1(1, t, B, 4), t) - kc1(t)) for t in ths[:1])
no_m = min(abs(dB(lambda B, t: K1(-1, t, B, 4), t) - kc1(t)) for t in ths[:1])


def kim(t, B, n):
    h = 2*n - 2
    sq = lambda x: bb(x - H, t, h)*bb(x + H, t, h)/(bb(x - H + B/2, t, h)*bb(x + H - B/2, t, h))
    return mp.fprod(sq(mp.mpf(k)/2) for k in range(1, 2*h - 2, 2))*sq(mp.mpf(h - 1)/2)


def caf(t, B, n):                    # CAF (dn1) with j = 1: R~_1(theta + i pi), blocks (blockw)
    h = 2*n - 2
    a1, a2 = 1j*mp.pi*(2 - B)/(2*h), 1j*mp.pi*B/(2*h)
    s_ = lambda z: mp.sinh(z/2)
    ang = lambda x, th: (s_(th + (x - 1)/2*a1 + (x - 1)/2*a2)*s_(th + (x + 1)/2*a1 + (x + 1)/2*a2)
                         / (s_(th + (x - 1)/2*a1 + (x + 1)/2*a2)*s_(th + (x + 1)/2*a1 + (x - 1)/2*a2)))
    W = lambda x, th: ang(x, th)/ang(x, -th)
    return W(h - 1, t)*mp.fprod(W(2*k - 1, t) for k in range(1, h))


ok_k = 0
for n in (4, 5, 6, 7):
    ok_k = max(ok_k, max(abs(dB(lambda B, t: K1(1, t, B, n, True), t) - dB(lambda B, t: kim(t, B, n), t)) for t in ths),
               max(abs(dB(lambda B, t: K1(-1, t, B, n), t) - dB(lambda B, t: kim(t + 1j*mp.pi, B, n), t)) for t in ths),
               abs(kim(mp.mpc('0.5', '0.2'), mp.mpf('0.6'), n)/caf(mp.mpc('0.5', '0.2'), mp.mpf('0.6'), n) - 1))
# a_3 (h = 4): K^(-)_1 = -(1)(h - B/2)(-1 + B/2) (sec-neumann-duality) against Kim (33) at theta + i pi
k33 = lambda t: 1j/16*(sh(t)/(ch(t) + 1/mp.sqrt(2)) - sh(t)/(ch(t) - 1))
ok_a3 = max(abs(dB(lambda B, t: -bb(1, t, 4)*bb(4 - B/2, t, 4)*bb(-1 + B/2, t, 4), t) - k33(t + 1j*mp.pi)) for t in ths)
report(f'one loop, d_4 (Kim-Cho (K1), (K2)): varsigma = +1 with the CDD factor agrees [{mp.nstr(ok_p, 2)}] (without it: '
       f'{mp.nstr(no_p, 2)}); varsigma = -1 agrees at theta + i pi [{mp.nstr(ok_m, 2)}] (at theta: {mp.nstr(no_m, 2)}); '
       f'the same against Kim\'s d_n conjecture = CAF R~_1(theta + i pi), n = 4..7 [{mp.nstr(ok_k, 2)}]; a_3: '
       f'K^(-)_1(theta) = Kim (33)(theta + i pi) [{mp.nstr(ok_a3, 2)}]',
       ok_p < 1e-12 and ok_m < 1e-12 and ok_k < 1e-12 and ok_a3 < 1e-12 and no_p > 0.01 and no_m > 0.01)

# ---------------------------------------------------------------- 8. semiclassics
mp.mp.dps = 30
th0 = mp.mpf('0.6')
vals = {}
for omv in (80, 160):
    om = mp.mpf(omv)
    T = 3*om
    sc = lambda f: mp.re(-1j*(mp.log(f(-1j*T*(th0 + mp.mpf('1e-8'))/mp.pi)) - mp.log(f(-1j*T*(th0 - mp.mpf('1e-8'))/mp.pi)))
                         / mp.mpf('2e-8'))*4*mp.pi/(om + 1)/(4*6)
    vals[omv] = [sc(lambda m: A(1, m, T, om)), sc(lambda m: A(-1, m, T, om)), sc(lambda m: F(2*m, T, om))/2]
lim = [2*vals[160][i] - vals[80][i] for i in range(3)]          # Richardson, errors O(1/omega)
d8 = max(abs(lim[0] - lim[2]), abs(lim[1] - lim[2]))
report(f'semiclassics (d_4): both varsigma = +-1 give half the time delay of F(2 theta), as a Neumann wall by images '
       f'[{mp.nstr(lim[0], 6)}, {mp.nstr(lim[1], 6)}, {mp.nstr(lim[2], 6)}]', d8 < 2e-3)

# ---------------------------------------------------------------- 9. relations to the uniform amplitude
def Ku(t, B, n):
    h = 2*n - 2
    b = lambda x: bb(x, t, h)
    return b(h - B/2)*b(h/2 - 1)*b(h/2 + 1 - B/2)/(b(h - 1)*b(1 - B/2)*b(h/2)*b(h/2 - B/2))


d9 = 0
for n in (4, 5, 6, 7):
    h = 2*n - 2
    for Bv in (mp.mpf('0.37'), mp.mpf('1.3')):
        for t in (mp.mpc('0.4', '0.1'), mp.mpc('-1.1', '0.3')):
            b = lambda x: bb(x, t, h)
            KN = b(h/2 - 1)*b(h/2 + Bv/2)*b(h - 1 + Bv/2)/(b(h - 1)*b(h/2)*b(h/2 - 1 + Bv/2)*b(Bv/2))
            Kcl = -b(h/2 - 1)*b(h/2 + 1)/(b(1)*b(h - 1)*b(h/2)**2)
            d9 = max(d9, abs(K1(1, t, Bv, n, True)/Ku(t, 2 - Bv, n) - 1), abs(K1(1, t, Bv, n, True)/KN - 1),
                     abs(K1(-1, t, Bv, n)*Kcl/Ku(t, Bv, n) - 1))
report(f'varsigma = +1 with the CDD factor is the dual K_1(theta, 2 - B) of the uniform amplitude, = K^N_1 in closed form; '
       f'varsigma = -1 is K_1/K_1^cl (n = 4..7) [{mp.nstr(d9, 2)}]', d9 < 1e-20)
