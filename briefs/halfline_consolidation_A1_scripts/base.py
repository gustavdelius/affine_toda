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



def singlet(A, B):
    d = A['dim']*B['dim']
    M = np.vstack([cop(A, B, kd, j) - (np.eye(d) if kd == 'k' else 0) for j in nodes_of(A) for kd in 'efk'])
    _, sv, Vh = np.linalg.svd(M)
    w = Vh[-1].conj()
    return w/w[int(np.argmax(np.abs(w) > 1e-8))]            # same normalization for in and out singlets


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

