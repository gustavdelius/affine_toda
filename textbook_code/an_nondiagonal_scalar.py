"""Scalar factor of the non-diagonal a_n^(1) reflection matrix (s = +1) and the Neumann amplitudes from the soliton side
(Part VII, ch. 33, sec-an-nondiagonal-k, prp-an-nondiagonal-scalar; ch. 35, sec-neumann-duality; half-line open
problem 3). Built from the problem-3 report and its referee scripts.

Book conventions: physical q = -exp(-i pi omega), V/Vbar = an_vector, R-check = rcheck, (eq-an-k-matrix) on the
analytic branch X^p = exp(p omega (theta - i pi/2)), eps_* = 1/(2 cos(pi omega/2)); h = n+1, T = h omega/2,
mu = -i T theta/pi, blocks (x)_I = sin(pi(mu+x)/2T)/sin(pi(mu-x)/2T).
Checks:
 1. Antisoliton and unitarity. J (basis reversal) is charge conjugation: J pi_Vbar(y)(x) J = pi_V(y)(tau x),
    tau: j -> n+1-j. K^T = J K J intertwines Vbar(y) -> V(1/y). K(-theta)^T K(theta) = U 1 with
    U = (s cosh 2T theta - cos pi T)/(2 sin^2(pi omega/2)); K(theta) K(-theta) is not a multiple of 1 (n >= 2).
    n = 2, 3, 4; omega = 0.37, 1.29; s = +-1.
 2. Singlet projection: (1 x K(theta_a)) R-check (1 x K^T(theta_b)) maps the B_1 singlet to c_+(mu) times the singlet,
    c_+ = -cos(pi mu) cos(pi(mu + omega/2)) cos(pi(mu - T))/(sin^2(pi omega/2) cos(pi(mu - omega/2))). n = 2, 3, 4.
 3. The scalar factor A^(+) (closed form of prp-an-nondiagonal-scalar; infinite Gamma products summed with Hurwitz zeta)
    solves A(mu)A(-mu) U = 1 and A(T/2-mu)/A(T/2+mu) = F_11(2mu)/Phi(mu), Phi = sin pi(omega/2-mu)/sin pi(omega/2+mu),
    and the breather identity A(mu_a) F_11(2mu) A(mu_b) c_+(mu) = (-n omega/2)_I (n omega/2 - 1/2)_I (1/2)_I,
    mu_{a,b} = mu -+ (T/2 - 1/2), holds. 30 digits; n = 2, 3, 4.
 4. Continuation omega -> -2/B: (-n omega/2)_I (n omega/2 - 1/2)_I (1/2)_I -> (-n)(n + B/2)(-B/2) = K^N_1
    (eq-neumann-amplitude), n = 1..7.
 5. Kim's one-loop a_3^(1) Neumann amplitudes (kim1995, eqs. (33), (39)) agree with K^N_1, K^N_2 at O(beta^2); without
    the CDD factor (T/2 - omega/2)_I (T/2 + omega/2)_I they become K^N_1/k_1, K^N_2/(k_1 k_1) (DG's k_1) and disagree.
Runtime: about 1 minute.
"""
import numpy as np
import mpmath as mp
from boundary_common import q_phys, an_vector, coideal_gen, rcheck, nodes_of, cop, kron, report


# ---------------------------------------------------------------- the K-matrix
def kmat(n, om, th, s):
    N = n + 1
    L = om*(th - 0.5j*np.pi)
    w = np.exp(-0.5j*np.pi*om)
    q = q_phys(om)
    eps = 1/(2*np.cos(np.pi*om/2))
    K = np.empty((N, N), complex)
    for i in range(N):
        for j in range(N):
            K[i, j] = ((w*np.exp(-N/2*L) - s/w*np.exp(N/2*L))/(eps*(1/q - q)) if i == j
                       else np.exp((j - i - N/2)*L) if j > i else s*np.exp((j - i + N/2)*L))
    return K


def inter_res(K, Vin, Vout, eps, q):
    return max(np.abs(K @ coideal_gen(Vin, j, 1/q, eps[j]) - coideal_gen(Vout, j, 1/q, eps[j]) @ K).max()
               for j in nodes_of(Vin))/np.abs(K).max()


# ---------------------------------------------------------------- 1. antisoliton matrix and unitarity
dJ, dT, dU, minoff = 0, 0, 0, np.inf
for n in (2, 3, 4):
    N = n + 1
    J = np.eye(N)[::-1]
    tau = lambda j: (N - j) % N
    for om in (0.37, 1.29):
        q = q_phys(om)
        es = 1/(2*np.cos(np.pi*om/2))
        th = 0.41 + 0.23j
        y = np.exp(om*th)
        V, Vb = an_vector(n, q, y), an_vector(n, q, y, True)
        dJ = max(dJ, max(np.abs(J @ Vb[(kd, j)] @ J - V[(kd, tau(j))]).max() for j in range(N) for kd in 'efk'))
        T = N*om/2
        for s in (1, -1):
            K = kmat(n, om, th, s)
            eps = es*np.ones(N)
            eps[0] *= s
            dT = max(dT, inter_res(K, V, an_vector(n, q, 1/y, True), eps, q),
                     inter_res(K.T, Vb, an_vector(n, q, 1/y), eps, q), np.abs(J @ K @ J - K.T).max()/np.abs(K).max())
            U = (s*np.cosh(2*T*th) - np.cos(np.pi*T))/(2*np.sin(np.pi*om/2)**2)
            dU = max(dU, np.abs(kmat(n, om, -th, s).T @ K - U*np.eye(N)).max()/abs(U))
            P = K @ kmat(n, om, -th, s)
            minoff = min(minoff, np.abs(P - np.diag(np.diag(P))).max()/abs(P[0, 0]))
report(f'J is charge conjugation [{dJ:.0e}]; K and K^T = JKJ are the soliton and antisoliton intertwiners [{dT:.0e}]; '
       f'K(-theta)^T K(theta) = U 1 [{dU:.0e}]; K(theta)K(-theta) is not scalar [off-diagonal/diagonal >= {minoff:.2f}] '
       f'(n = 2, 3, 4; omega = 0.37, 1.29; s = +-1)', dJ < 1e-14 and dT < 1e-12 and dU < 1e-12 and minoff > 0.1)


# ---------------------------------------------------------------- 2. singlet projection
def singlet(A, B):
    d = A['dim']*B['dim']
    M = np.vstack([cop(A, B, kd, j) - (np.eye(d) if kd == 'k' else 0) for j in nodes_of(A) for kd in 'efk'])
    _, sv, Vh = np.linalg.svd(M)
    w = Vh[-1].conj()
    return w/w[int(np.argmax(np.abs(w) > 1e-8))]            # same normalization for in and out singlets


dP = 0
for n, om, mu in ((2, 0.53, 0.11 + 0.19j), (3, 0.61, 0.17 + 0.23j), (4, 0.47, -0.21 + 0.13j)):
    h, N = n + 1, n + 1
    T, q = h*om/2, q_phys(om)
    th = 1j*np.pi*mu/T
    u = np.pi*(1 - 2/(h*om))
    ta, tb = th - 0.5j*u, th + 0.5j*u
    Va, Vb = an_vector(n, q, np.exp(om*ta)), an_vector(n, q, np.exp(om*tb), True)
    w_in = singlet(Va, Vb)
    Wa = an_vector(n, q, np.exp(-om*tb))
    w_out = singlet(Wa, an_vector(n, q, np.exp(-om*ta), True))
    R = rcheck(Va, Wa)[0]
    v = kron(np.eye(N), kmat(n, om, ta, 1)) @ (R @ (kron(np.eye(N), kmat(n, om, tb, 1).T) @ w_in))
    c = np.vdot(w_out, v)/np.vdot(w_out, w_out)
    cp = (-np.cos(np.pi*mu)*np.cos(np.pi*(mu + om/2))*np.cos(np.pi*(mu - T))
          / (np.sin(np.pi*om/2)**2*np.cos(np.pi*(mu - om/2))))
    dP = max(dP, np.linalg.norm(v - c*w_out)/np.linalg.norm(v), abs(c/cp - 1))
report(f'B_1 singlet is mapped to c_+(mu) times itself (n = 2, 3, 4) [{dP:.0e}]', dP < 1e-12)

# ---------------------------------------------------------------- 3. scalar factor and breather identity
mp.mp.dps = 30


def logprod(P, terms, J=60, M=40):
    """log prod_{j>=1} prod_i Gamma(jP + a_i)^{c_i} for balanced terms: exact head, Stirling tail via Hurwitz zeta."""
    head = sum(c*mp.loggamma(j*P + a) for j in range(1, J) for c, a in terms)
    tail = sum((-1)**(k + 1)*sum(c*mp.bernpoly(k + 1, a) for c, a in terms)/(k*(k + 1)*P**k)*mp.zeta(k, J)
               for k in range(2, M + 1))
    return head + tail


def G(m, a, sign=1):
    return [(sign, m + a), (-sign, -m + a)]


def blk(x, mu, T):
    return mp.sin(mp.pi*(mu + x)/(2*T))/mp.sin(mp.pi*(mu - x)/(2*T))


def scalar_data(n, om):
    h = n + 1
    T, L = h*om/2, h*om

    def Ap(mu):
        a1 = mp.exp(logprod(2*L, G(2*mu, -1.5*L + 1) + G(2*mu, -0.5*L) + G(2*mu, -1.5*n*om - 0.5*om + 1, -1)
                            + G(2*mu, -0.5*n*om - 1.5*om, -1)))           # DG (3.9) without its prefactor
        a0 = mp.sin(mp.pi*om/2)/mp.sin(mp.pi*(mu - T/2))*mp.exp(
            logprod(L, G(mu, om/2 - 1.5*T) + G(mu, 1 - om/2 - T/2) + G(mu, om/2 - T/2, -1) + G(mu, 1 - om/2 - 1.5*T, -1)))
        return a0*a1*blk(T/2 - om/2, mu, T)*blk(T/2 + om/2, mu, T)

    def F11(m):                                                           # prp-an-diagonal-scalar
        t = [(1, m - (h - 1)*om + 1), (1, m - om), (1, -m - L + 1), (1, -m),
             (-1, -m - (h - 1)*om + 1), (-1, -m - om), (-1, m - L + 1), (-1, m)]
        return -mp.exp(logprod(L, t))
    return T, Ap, F11


d3 = 0
for n, om, mu in ((2, mp.mpf('0.53'), mp.mpc('0.11', '0.19')), (3, mp.mpf('0.61'), mp.mpc('0.17', '0.23')),
                  (4, mp.mpf('0.47'), mp.mpc('-0.21', '0.13'))):
    T, Ap, F11 = scalar_data(n, om)
    U = mp.sin(mp.pi*(T/2 + mu))*mp.sin(mp.pi*(T/2 - mu))/mp.sin(mp.pi*om/2)**2
    Phi = mp.sin(mp.pi*(om/2 - mu))/mp.sin(mp.pi*(om/2 + mu))
    cp = (-mp.cos(mp.pi*mu)*mp.cos(mp.pi*(mu + om/2))*mp.cos(mp.pi*(mu - T))
          / (mp.sin(mp.pi*om/2)**2*mp.cos(mp.pi*(mu - om/2))))
    half = mp.mpf(1)/2
    lhs = Ap(mu - T/2 + half)*F11(2*mu)*Ap(mu + T/2 - half)*cp
    rhs = blk(-n*om/2, mu, T)*blk(n*om/2 - half, mu, T)*blk(half, mu, T)
    d3 = max(d3, abs(Ap(mu)*Ap(-mu)*U - 1), abs(Ap(T/2 - mu)/Ap(T/2 + mu)*Phi/F11(2*mu) - 1), abs(lhs/rhs - 1))
report(f'A^(+) solves unitarity and crossing; A(mu_a) F_11(2mu) A(mu_b) c_+ = (-n omega/2)_I (n omega/2 - 1/2)_I (1/2)_I '
       f'(n = 2, 3, 4) [{mp.nstr(d3, 2)}]', d3 < mp.mpf('1e-20'))


# ---------------------------------------------------------------- 4. continuation to K^N_1
def bth(x, t, h):
    return mp.sinh(t/2 + 1j*mp.pi*x/(2*h))/mp.sinh(t/2 - 1j*mp.pi*x/(2*h))


d4 = 0
for n in range(1, 8):
    h = n + 1
    for Bv in (mp.mpf('0.37'), mp.mpf('1.3')):
        for t in (mp.mpc('0.4', '0.1'), mp.mpc('-1.1', '0.3')):
            om = -2/Bv
            T = h*om/2
            mu = -1j*T*t/mp.pi
            cont = blk(-n*om/2, mu, T)*blk(n*om/2 - mp.mpf(1)/2, mu, T)*blk(mp.mpf(1)/2, mu, T)
            KN1 = bth(0, t, h)*bth(h - 1 + Bv/2, t, h)/(bth(h - 1, t, h)*bth(Bv/2, t, h))
            d4 = max(d4, abs(cont/KN1 - 1))
report(f'omega -> -2/B: the breather amplitude becomes K^N_1 of eq-neumann-amplitude (n = 1..7) [{mp.nstr(d4, 2)}]',
       d4 < mp.mpf('1e-25'))

# ---------------------------------------------------------------- 5. Kim's one loop (a_3^(1), h = 4)
mp.mp.dps = 50
h = 4
KN = lambda a, Bv, t: mp.fprod(bth(c - 1, t, h)*bth(h - c + Bv/2, t, h)/(bth(h - c, t, h)*bth(c - 1 + Bv/2, t, h))
                               for c in range(1, a + 1))
k1 = lambda Bv, t: (bth(h/2 + 1 + Bv/2, t, h)*bth(h/2 - 1 - Bv/2, t, h)*bth(-h/2 + 1 - Bv/2, t, h)
                    * bth(-h/2 - 1 + Bv/2, t, h))                                        # DG (3.25), c = 1
dB = lambda f, t, e=mp.mpf('1e-20'): (f(e, t) - f(-e, t))/(2*e)*(1/(2*mp.pi))           # d/d(beta^2) at 0
kim33 = lambda t: 1j/16*(mp.sinh(t)/(mp.cosh(t) + 1/mp.sqrt(2)) - mp.sinh(t)/(mp.cosh(t) - 1))
kim39 = lambda t: 1j/16*(mp.sinh(t)/mp.cosh(t) - mp.sinh(t)/(mp.cosh(t) - 1)
                         - mp.sinh(t)/(mp.cosh(t) - 1/mp.sqrt(2)) + mp.sinh(t)/(mp.cosh(t) + 1/mp.sqrt(2)))
ths = [mp.mpf('0.37'), mp.mpc('1.2', '0.3'), mp.mpc('-0.8', '0.1')]
w_with = max(max(abs(dB(lambda b, t: KN(1, b, t), t) - kim33(t)), abs(dB(lambda b, t: KN(2, b, t), t) - kim39(t)))
             for t in ths)
w_without = max(max(abs(dB(lambda b, t: KN(1, b, t)/k1(b, t), t) - kim33(t)),
                    abs(dB(lambda b, t: KN(2, b, t)/(k1(b, t + 1j*mp.pi/4)*k1(b, t - 1j*mp.pi/4)), t) - kim39(t)))
                for t in ths)
report(f'Kim (33), (39): K^N_1, K^N_2 agree at O(beta^2) [{mp.nstr(w_with, 2)}]; without the CDD factor '
       f'(T/2 - omega/2)_I (T/2 + omega/2)_I they deviate by {mp.nstr(w_without, 3)}', w_with < 1e-25 and w_without > 0.05)

# ---------------------------------------------------------------- 6. the crossing function Phi
# k(theta) = Phi(theta) R-check_{Vbar Vbar}(2 theta) k(-theta), k = sum_ab K_ba(i pi/2 - theta) v_a (x) v_b
# (boundary_coideal.py checks proportionality); Phi = sinh(i pi omega/2 - T theta)/sinh(i pi omega/2 + T theta) for
# s = +1 and cosh(T theta - i pi omega/2)/cosh(T theta + i pi omega/2) for s = -1, independent of n except through T.
d6 = 0
for n in (2, 3, 4):
    N = n + 1
    for om in (0.37, 1.29):
        q, T = q_phys(om), N*om/2
        for s in (1, -1):
            for th in (0.31, 0.83 + 0.1j):
                km = kmat(n, om, 1j*np.pi/2 - th, s).T.reshape(N*N)
                kp = kmat(n, om, 1j*np.pi/2 + th, s).T.reshape(N*N)
                v = rcheck(an_vector(n, q, np.exp(2*om*th), True), an_vector(n, q, 1.0, True))[0] @ kp
                f = np.vdot(v, km)/np.vdot(v, v)
                Phi = (np.sinh(1j*np.pi*om/2 - T*th)/np.sinh(1j*np.pi*om/2 + T*th) if s == 1
                       else np.cosh(T*th - 1j*np.pi*om/2)/np.cosh(T*th + 1j*np.pi*om/2))
                d6 = max(d6, np.linalg.norm(km - f*v)/np.linalg.norm(km), abs(f/Phi - 1))
report(f'crossing-unitarity k(theta) = Phi R-check_(Vbar Vbar)(2theta) k(-theta) with the closed-form Phi, s = +-1, '
       f'n = 2, 3, 4, omega = 0.37, 1.29 [{d6:.0e}]', d6 < 1e-11)

# ---------------------------------------------------------------- 7. semiclassics of A^(+): the Neumann time advance
# Delta t = (1/(m_a sinh theta)) ln(1 - s^2), s = sin(pi/h)/cosh theta (sec-images-neumann); with beta^2 = 4 pi/(omega+1)
# as in reflection_matrices.py check 6: (beta^2/4h) d(-i log A^(+))/d theta -> (1/2) ln(1 - s^2) as omega -> infinity.
mp.mp.dps = 30
ok7, msgs = True, []
for n in (2, 3):
    h = n + 1
    th = mp.mpf('0.6')
    errs = []
    for omv in (40, 80, 160):
        om = mp.mpf(omv)
        T, Ap, _ = scalar_data(n, om)
        f = lambda t: mp.log(Ap(-1j*T*t/mp.pi))
        d = (f(th + mp.mpf('1e-8')) - f(th - mp.mpf('1e-8')))/mp.mpf('2e-8')
        errs.append(abs(mp.re(-1j*d)*4*mp.pi/(om + 1)/(4*h) - mp.log(1 - (mp.sin(mp.pi/h)/mp.cosh(th))**2)/2))
    ok7 &= errs[2] < errs[1] < errs[0] and abs(errs[0]/errs[2] - 4) < 0.5
    msgs.append(f'n = {n}: ' + ', '.join(mp.nstr(e, 2) for e in errs))
report('semiclassical limit of A^(+): the Neumann time advance (1/2) ln(1 - sin^2(pi/h)/cosh^2 theta), errors at '
       'omega = 40, 80, 160 decreasing like 1/omega (' + '; '.join(msgs) + ')', ok7)
