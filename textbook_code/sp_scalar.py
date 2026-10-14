"""Scalar factors, soliton bootstrap and breather amplitudes of the soliton-preserving (SP) a_n^(1) reflection matrices
(Part VII, ch. 33, sec-preserving-conjugating, prp-sp-scalar; ch. 35, sec-breather-reflection, half-line open problem 11
in sec-boundary-open-problems).  Consolidation of the problem-11 report; independent of the worker's scripts.

Book conventions: physical q = -exp(-i pi omega), V/Vbar = an_vector, R-check = rcheck (normalized on hw (x) hw),
h = N = n+1, T = h omega/2, mu = -i T theta/pi, blocks (x)_I = sin(pi(mu+x)/2T)/sin(pi(mu-x)/2T).
SP matrices (soliton_preserving.py): K_V = diag(e^{a_j omega theta}) diag(sinh(eta +- T theta)) (+ and a_j = 2j for
j < l, - and a_j = 2j - N otherwise), K_Vbar with etabar = eta + i pi omega (l - N/2) and the opposite signs.
eps = eta/(i pi), epsbar = etabar/(i pi) = eps + omega l - T.  F11 of prp-an-diagonal-scalar, S^T(x) = sin(pi x)/
sin(pi(omega - x)) F11(x).  Ehat_P(z) = Gamma(z) prod_k Gamma(z+kP)/Gamma(1-z+kP) (kP)^(1-2z),
Lam_P(x; c) = Ehat_P(c+x)/Ehat_P(c-x); infinite products summed with a Stirling tail (Hurwitz zeta), 30 digits.
Checks:
 1. Bulk crossing with the pairing v_i <-> vbar_i: Rc_VV(mu) = t(mu) Rc_VbarV(T-mu) and = t(mu) Rc_VVbar(T-mu) in the
    crossed index orders, t = sin(pi mu)/sin(pi(omega-mu)); all transmission entries of Rc_VbarV, Rc_VVbar are 1.
    Hence the soliton-antisoliton S-matrix is S^T(T - mu) Rc.  n = 2, 3; omega = 0.37, 1.29.
 2. SP crossing-unitarity (eq-boundary-crossing-unitarity): with kA = sum_i K_Vbar(i pi/2 - th)_ii v_i (x) vbar_i and
    kB = sum_i K_V(i pi/2 - th)_ii vbar_i (x) v_i, kB(th) = Phi1 Rc_VVbar(2th) kA(-th) and
    kA(th) = Phi2 Rc_VbarV(2th) kB(-th) with Phi1 = -(-q)^(1-l) s, Phi2 = -(-q)^(l-1) s,
    s = sin pi(2mu - T)/sin pi(2mu + T), independent of eta.  n = 2, 3, all l; omega = 0.37, 1.29.
 3. The scalar factors a (on K_V) and bt (on K_Vbar),
      a  = s_ rho(mu) e^{(1-l) omega theta} i Lam_2T(mu; eps+2T) / (Lam_2T(mu; epsbar+T) sin pi(eps+mu)),
      bt = s_ rho(mu) e^{-(1-l) omega theta} i Lam_2T(mu; epsbar+2T) / (Lam_2T(mu; eps+T) sin pi(epsbar+mu)),
      rho = Lam_4T(2mu; 2T+omega)/Lam_4T(2mu; 0) (omega/2)_I,
    solve a(mu)a(-mu)U = 1, bt(mu)bt(-mu)Ubar = 1, a(T/2-mu)Phi1 = F(2mu) bt(T/2+mu), bt(T/2-mu)Phi2 = F(2mu) a(T/2+mu),
    U = sinh(eta + i pi mu) sinh(eta - i pi mu), Ubar likewise with etabar, F(x) = S^T(T - x).  (n, l) = (1,1), (2,0),
    (2,1), (3,0), (3,1), (3,2).  (s_ = +-1; the bootstrap, check 5, fixes s_ = -1 for even n.)  Also: the same
    functions in the book's G_j notation, Lam_P(x;c)/Lam_P(x;c') = calG_{c,c'}(x) = G_0(x,c)/G_0(x,c')
    prod_{j>=1} G_j(x,c)G_j(x,1-c)/(G_j(x,c')G_j(x,1-c')) (period P = 2T in mu; for rho, G_2j in 2mu).
 4. n = 1, l = 1: Ghoshal-Zamolodchikov Dirichlet (Dirichlet limit of eq-bsg-k-full with the book's R_0, sigma;
    xi = pi eps - pi/2, lambda = omega) equals -a (-T/2)_I^2 K_V for every eta.
 5. Soliton bootstrap V^(x n) -> Vbar: Vbar(th) sits in V(th'_1) (x) ... (x) V(th'_n), mu'_k = mu - T + k omega (and not
    in the reversed order).  Reflecting the n solitons (rcheck, K_V) maps the embedding iota to P iota K_Vbar(th), and
    bt(mu) = (-1)^(n+1) P prod_k a(mu'_k) prod_{j<k} F11(mu'_j + mu'_k) for s_ = 1; the exponential solutions
    a e^{m theta}, bt (-1)^m e^{-m theta} (m = +-2) of the same functional equations fail.  n = 1, 2, 3.
 6. CDD factors (sigma_a on the soliton, sigma_b(mu) = sigma_a(T - mu) on the antisoliton) compatible with check 5,
    i.e. sigma_b(mu)/prod_k sigma_a(mu'_k) = +-1: sigma_B = -[(T/2)_I^2 (T/2-omega/2)_I (T/2+omega/2)_I]^-1 and the pair
    (omega)_I, -(T-omega)_I pass for n = 1..5; sigma_A = (-T/2)_I^2 passes for odd n and fails for n = 2, 4.
 7. Breather B_1 -> Bbar_1 (sec-breather-reflection; singlets normalized to 1 on hw (x) lw of their two factors):
    a(mu_a) F(2mu) bt(mu_b) c(mu) = e^{-2 pi i l/h} K0(mu) E(theta), mu_{a,b} = mu -+ (T-1)/2, with
    K0 = (T/2)_I (T/2-1/2)_I (d)_I (1/2-d)_I (1/2+d)_I (T-1/2-d)_I, d = (T - omega)/2,
    E = (sinh th - i sin(pi A/T))/(sinh th + i sin(pi Abar/T)), A = eps + (T-1)/2, Abar = A + omega l - T;
    Bbar_1 -> B_1: e^{2 pi i l/h} K0 Ebar (A <-> Abar).  c from the zero-weight block of Rc_VVbar (closed form checked
    against rcheck) at 30 digits.  n = 1, 2, 3.  For n = 1, l = 1 with sigma_A: Ghoshal's K_B1 (eq-bsg-b1-reflection)
    at Dirichlet.
 8. Continuation omega -> -2/B: K0 -> (h/2)(h/2+B/2)(h/2-1)(h/2-1-B/2)(h/2+1+B/2)/(h/2-1+B/2) (n = 1..7).  With
    L_1 = K0 E, L_n = K0 Ebar and L_{a+1}(th) = L_a(th + i pi/h) S_{a,h-1}(2th + i pi(1-a)/h) L_1(th - i pi a/h)
    (particle_reflection.py's S): the bootstrap gives the soliton-side L_n; unitarity L_a(th) L_{h-a}(-th) = 1;
    crossing L_a(i pi/2 - th) = S_aa(2th) L_a(i pi/2 + th); the bootstrap L_c = L_a S_{a,h-b} L_b at every fusing;
    minimal and sigma_B, n = 2, 3; sigma_A breaks the bootstrap for n = 2.  B -> 0: L_1 -> g_{h/2} g_{h/2-1} E with
    g_x = (sinh th + i sin(pi x/h))/(sinh th - i sin(pi x/h)), not constant for h >= 3.
Runtime: about four and a half minutes.
"""
import os
import numpy as np
import mpmath as mp
import boundary_common
from boundary_common import an_vector, rcheck, q_phys, cop, nodes_of, kron, intertwiners, report

_src = os.path.join(os.path.dirname(os.path.abspath(boundary_common.__file__)), 'particle_reflection.py')
_ns = {}
exec(open(_src).read().split('THS = [')[0], _ns)          # blocks (x) and the a_n S-matrix only
blk, S_ab = _ns['blk'], _ns['S']
mp.mp.dps = 30
I, PI, HALF = mp.mpc(0, 1), mp.pi, mp.mpf(1)/2


# ---------------------------------------------------------------- matrices (double precision, book code)
def K_V(n, om, l, eta, th):
    N, T = n + 1, (n + 1)*om/2
    return np.diag([np.exp(2*j*om*th)*np.sinh(eta + T*th) if j < l
                    else np.exp((2*j - N)*om*th)*np.sinh(eta - T*th) for j in range(N)])


def K_Vbar(n, om, l, eta, th):
    N, T = n + 1, (n + 1)*om/2
    etab = eta + 1j*np.pi*om*(l - N/2)
    return np.diag([np.exp(-2*j*om*th)*np.sinh(etab - T*th) if j < l
                    else np.exp(-(2*j - N)*om*th)*np.sinh(etab + T*th) for j in range(N)])


_rc = {}


def Rc(n, om, a, b, z):
    """rcheck on a(z) (x) b(1), a, b in 'V', 'W' (W = Vbar); falls back to LAPACK gesvd if gesdd fails to converge."""
    key = (n, om, a, b, complex(np.round(z, 13)))
    if key not in _rc:
        q = q_phys(om)
        A, B = an_vector(n, q, z, conj=(a == 'W')), an_vector(n, q, 1.0, conj=(b == 'W'))
        try:
            _rc[key] = rcheck(A, B)[0]
        except np.linalg.LinAlgError:
            import scipy.linalg as sl
            d = A['dim']*B['dim']
            M = np.vstack([np.kron(cop(A, B, kd, j).T, np.eye(d)) - np.kron(np.eye(d), cop(B, A, kd, j))
                           for j in nodes_of(A) for kd in 'efk'])
            R = sl.svd(M, full_matrices=False, lapack_driver='gesvd')[2][-1].conj().reshape(d, d).T
            _rc[key] = R/R[B['hw']*A['dim'] + A['hw'], A['hw']*B['dim'] + B['hw']]
    return _rc[key]


# ---------------------------------------------------------------- 30-digit scalar functions
def _tail(P, terms, J, M):
    return sum((-1)**(k + 1)*sum(c*mp.bernpoly(k + 1, a) for c, a in terms)/(k*(k + 1)*P**k)*mp.zeta(k, J)
               for k in range(2, M + 1))


def logprod(P, terms, J=50, M=36):
    """log prod_{j>=1} prod_i Gamma(jP + a_i)^{c_i} for balanced terms."""
    return sum(c*mp.loggamma(j*P + a) for j in range(1, J) for c, a in terms) + _tail(P, terms, J, M)


def logEhat(z, P, J=50, M=36):
    return (mp.loggamma(z) + sum(mp.loggamma(z + k*P) - mp.loggamma(1 - z + k*P) + (1 - 2*z)*mp.log(k*P)
                                 for k in range(1, J)) + _tail(P, [(1, z), (-1, 1 - z)], J, M))


def Lam(x, c, P):
    return mp.exp(logEhat(c + x, P) - logEhat(c - x, P))


def bI(x, mu, T):
    return mp.sin(PI*(mu + x)/(2*T))/mp.sin(PI*(mu - x)/(2*T))


def F11(m, n, om):
    h, L = n + 1, (n + 1)*om
    t = [(1, m - (h - 1)*om + 1), (1, m - om), (1, -m - L + 1), (1, -m),
         (-1, -m - (h - 1)*om + 1), (-1, -m - om), (-1, m - L + 1), (-1, m)]
    return -mp.exp(logprod(L, t))


class SP:
    def __init__(self, n, om, l, eta):
        self.n, self.l, self.N = n, l, n + 1
        self.om = mp.mpf(om)
        self.T = self.N*self.om/2
        self.eta = mp.mpmathify(eta)
        self.etab = self.eta + I*PI*self.om*(l - mp.mpf(self.N)/2)
        self.eps, self.epsb = self.eta/(I*PI), self.etab/(I*PI)
        self.q = -mp.exp(-I*PI*self.om)
        self.lam = mp.exp(-I*PI*self.om*(1 - l))                     # (-q)^(1-l)

    def s(self, mu):
        return mp.sin(PI*(2*mu - self.T))/mp.sin(PI*(2*mu + self.T))

    def F(self, x):
        y = self.T - x
        return mp.sin(PI*y)/mp.sin(PI*(self.om - y))*F11(y, self.n, self.om)

    def rho(self, mu):
        T = self.T
        return Lam(2*mu, 2*T + self.om, 4*T)/Lam(2*mu, 0, 4*T)*bI(self.om/2, mu, T)

    def a(self, mu, m=0):
        T = self.T
        return (self.rho(mu)*mp.exp((1 - self.l)*self.om*I*PI*mu/T + I*PI*m*mu/T)*I*Lam(mu, self.eps + 2*T, 2*T)
                / (Lam(mu, self.epsb + T, 2*T)*mp.sin(PI*(self.eps + mu))))

    def bt(self, mu, m=0):
        T = self.T
        return ((-1)**m*self.rho(mu)*mp.exp(-(1 - self.l)*self.om*I*PI*mu/T - I*PI*m*mu/T)*I
                * Lam(mu, self.epsb + 2*T, 2*T)/(Lam(mu, self.eps + T, 2*T)*mp.sin(PI*(self.epsb + mu))))


def K0(mu, T, om):
    d = (T - om)/2
    return (bI(T/2, mu, T)*bI(T/2 - HALF, mu, T)*bI(d, mu, T)*bI(HALF - d, mu, T)*bI(HALF + d, mu, T)
            * bI(T - HALF - d, mu, T))


# ---------------------------------------------------------------- 1. bulk crossing, soliton-antisoliton normalization
d1 = dtr = 0
for n in (2, 3):
    N = n + 1
    for om in (0.37, 1.29):
        T = N*om/2
        th = 0.31 + 0.12j
        mu = -1j*T*th/np.pi
        t = np.sin(np.pi*mu)/np.sin(np.pi*(om - mu))
        RVV = Rc(n, om, 'V', 'V', np.exp(om*th))
        RWV, RVW = Rc(n, om, 'W', 'V', np.exp(om*(1j*np.pi - th))), Rc(n, om, 'V', 'W', np.exp(om*(1j*np.pi - th)))
        for a_ in range(N):
            for b_ in range(N):
                for c_ in range(N):
                    for e_ in range(N):
                        lhs = RVV[e_*N + c_, a_*N + b_]
                        d1 = max(d1, abs(lhs - t*RWV[c_*N + b_, e_*N + a_]), abs(lhs - t*RVW[a_*N + e_, b_*N + c_]))
        dtr = max(dtr, max(abs(R[j*N + i, i*N + j] - 1) for R in (RWV, RVW) for i in range(N) for j in range(N) if i != j))
report(f'bulk crossing with v_i <-> vbar_i: Rc_VV(mu) = t(mu) Rc_(VbarV, VVbar)(T - mu) [{d1:.0e}], transmission entries '
       f'= 1 [{dtr:.0e}]; so S_(VbarV) = S^T(T - mu) Rc (n = 2, 3; omega = 0.37, 1.29)', d1 < 1e-12 and dtr < 1e-12)

# ---------------------------------------------------------------- 2. SP crossing functions
d2 = 0
for n in (2, 3):
    N = n + 1
    E = np.eye(N)
    for om in (0.37, 1.29):
        T, q = N*om/2, q_phys(om)
        th = -0.45 + 0.27j
        mu = -1j*T*th/np.pi
        s = np.sin(np.pi*(2*mu - T))/np.sin(np.pi*(2*mu + T))
        R1, R2 = Rc(n, om, 'V', 'W', np.exp(2*om*th)), Rc(n, om, 'W', 'V', np.exp(2*om*th))
        for l in range(N + 1):
            for eta in (0.41 + 0.23j, -0.3 + 1.1j):
                kA = lambda t: sum(K_Vbar(n, om, l, eta, 0.5j*np.pi - t)[i, i]*np.kron(E[i], E[i]) for i in range(N))
                kB = lambda t: sum(K_V(n, om, l, eta, 0.5j*np.pi - t)[i, i]*np.kron(E[i], E[i]) for i in range(N))
                for R, k1, k2, Phi in ((R1, kB, kA, -(-q)**(1 - l)*s), (R2, kA, kB, -(-q)**(l - 1)*s)):
                    v = Phi*(R @ k2(-th))
                    d2 = max(d2, np.linalg.norm(k1(th) - v)/np.linalg.norm(k1(th)))
report(f'SP crossing-unitarity: kB = Phi1 Rc_VVbar(2th) kA(-th), kA = Phi2 Rc_VbarV(2th) kB(-th), Phi1,2 = -(-q)^(+-(1-l)) '
       f'sin pi(2mu-T)/sin pi(2mu+T), independent of eta (n = 2, 3, all l, omega = 0.37, 1.29) [{d2:.0e}]', d2 < 1e-12)

# ---------------------------------------------------------------- 3. scalar factors
d3 = dG = 0
mu = mp.mpc('0.11', '0.19')
for n, l in ((1, 1), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)):
    Sp = SP(n, '0.37', l, mp.mpc('0.41', '0.23'))
    T = Sp.T
    U = mp.sinh(Sp.eta + I*PI*mu)*mp.sinh(Sp.eta - I*PI*mu)
    Ub = mp.sinh(Sp.etab + I*PI*mu)*mp.sinh(Sp.etab - I*PI*mu)
    Phi1, Phi2 = -Sp.lam*Sp.s(mu), -Sp.s(mu)/Sp.lam
    F2 = Sp.F(2*mu)
    d3 = max(d3, abs(Sp.a(mu)*Sp.a(-mu)*U - 1), abs(Sp.bt(mu)*Sp.bt(-mu)*Ub - 1),
             abs(Sp.a(T/2 - mu)*Phi1/(F2*Sp.bt(T/2 + mu)) - 1), abs(Sp.bt(T/2 - mu)*Phi2/(F2*Sp.a(T/2 + mu)) - 1))
    # the same a, bt in the book's G_j notation (G_j(x, c) = Gamma(x + j h omega + c)/Gamma(-x + j h omega + c)):
    # calG_{c,c'}(x; P) = G_0(x,c)/G_0(x,c') prod_j G_j(x,c)G_j(x,1-c)/(G_j(x,c')G_j(x,1-c')), j-th period j P
    Gt = lambda x, c, sg=1: [(sg, x + c), (-sg, -x + c)]
    calG = lambda x, c, c2, P: (mp.gamma(c + x)*mp.gamma(c2 - x)/(mp.gamma(c - x)*mp.gamma(c2 + x))
                                * mp.exp(logprod(P, Gt(x, c) + Gt(x, 1 - c) + Gt(x, c2, -1) + Gt(x, 1 - c2, -1))))
    rhoG = bI(Sp.om/2, mu, T)*calG(2*mu, 2*T + Sp.om, 0, 4*T)
    aG = rhoG*mp.exp(2*PI*I*(1 - l)*mu/(n + 1))*I*calG(mu, Sp.eps + 2*T, Sp.epsb + T, 2*T)/mp.sin(PI*(Sp.eps + mu))
    bG = rhoG*mp.exp(-2*PI*I*(1 - l)*mu/(n + 1))*I*calG(mu, Sp.epsb + 2*T, Sp.eps + T, 2*T)/mp.sin(PI*(Sp.epsb + mu))
    dG = max(dG, abs(aG/Sp.a(mu) - 1), abs(bG/Sp.bt(mu) - 1))
report(f'scalar factors a, bt solve unitarity and crossing-unitarity, (n, l) = (1,1), (2,0), (2,1), (3,0), (3,1), (3,2) '
       f'[{mp.nstr(d3, 2)}]; G_j form = Lam form [{mp.nstr(dG, 2)}]', d3 < mp.mpf('1e-24') and dG < mp.mpf('1e-24'))


# ---------------------------------------------------------------- 4. n = 1: Ghoshal-Zamolodchikov Dirichlet
def R0(u, lam):
    t = lambda s: [(1, -2*lam*s*u/PI), (1, -4*lam + 1 - 2*lam*s*u/PI), (-1, -3*lam - 2*lam*s*u/PI),
                   (-1, -lam + 1 - 2*lam*s*u/PI)]
    return mp.exp(logprod(4*lam, t(1) + [(-c, a) for c, a in t(-1)]))


def sigma(x, u, lam):
    t = lambda s: [(1, HALF + x/PI - lam - lam*s*u/PI), (1, HALF - x/PI - lam - lam*s*u/PI),
                   (-1, HALF - x/PI - 2*lam - lam*s*u/PI), (-1, HALF + x/PI - lam*s*u/PI)]
    return mp.cos(x)/mp.cos(x + lam*u)*mp.exp(logprod(2*lam, t(1) + [(-c, a) for c, a in t(-1)]))


d4 = 0
for eta in (mp.mpc('0.41', '0.23'), mp.mpc('-0.3', '1.1')):
    Sp = SP(1, '0.37', 1, eta)
    lam, T = Sp.om, Sp.T
    xi = PI*Sp.eps - PI/2
    for mu in (mp.mpc('0.11', '0.19'), mp.mpc('-0.23', '0.07')):
        u = PI*mu/lam
        gz = R0(u, lam)*sigma(xi, u, lam)/mp.cos(xi)
        mine = -Sp.a(mu)*bI(-T/2, mu, T)**2
        d4 = max(d4, abs(gz*mp.cos(xi + lam*u)/(mine*mp.sinh(Sp.eta + T*I*u)) - 1),
                 abs(gz*mp.cos(xi - lam*u)/(mine*mp.sinh(Sp.eta - T*I*u)) - 1))
report(f'n = 1: Ghoshal-Zamolodchikov Dirichlet (book R_0, sigma; xi = pi eps - pi/2) = -a (-T/2)_I^2 K_V for all eta '
       f'[{mp.nstr(d4, 2)}]', d4 < mp.mpf('1e-24'))


# ---------------------------------------------------------------- 5. soliton bootstrap V^(x n) -> Vbar
def tensor_rep(reps):
    A = reps[0]
    for B in reps[1:]:
        g = {(kd, j): cop(A, B, kd, j) for j in nodes_of(A) for kd in 'efk'}
        g['dim'], g['hw'] = A['dim']*B['dim'], None
        A = g
    return A


def embedding(n, om, th, order=1):
    T, q = (n + 1)*om/2, q_phys(om)
    mu = -1j*T*th/np.pi
    mus = [mu - T + k*om for k in range(1, n + 1)][::order]
    W = tensor_rep([an_vector(n, q, np.exp(om*1j*np.pi*m/T)) for m in mus])
    Ub = an_vector(n, q, np.exp(om*th), conj=True)
    null, rel = intertwiners([(Ub[(kd, j)], W[(kd, j)]) for j in nodes_of(Ub) for kd in 'efk'], tol=1e-9)
    return null, mus


def boot(n, om, l, eta, th):
    N, T = n + 1, (n + 1)*om/2
    (io_in,), mus = embedding(n, om, th)
    (io_out,), _ = embedding(n, om, -th)
    i0 = np.unravel_index(np.argmax(np.abs(io_in)), io_in.shape)
    io_in, io_out = io_in/io_in[i0], io_out/io_out[i0]
    ths = [1j*np.pi*m/T for m in mus]
    v, pos = io_in.copy(), list(range(n))
    for k in range(n - 1, -1, -1):
        for p in range(k, n - 1):
            M2 = Rc(n, om, 'V', 'V', np.exp(om*(ths[k] + ths[pos[p + 1]])))
            v = np.einsum('ab,ibjm->iajm', M2, v.reshape((N**p, N*N, N**(n - p - 2), N))).reshape(N**n, N)
            pos[p], pos[p + 1] = pos[p + 1], pos[p]
        v = np.einsum('ab,ibjm->iajm', K_V(n, om, l, eta, ths[k]), v.reshape((N**(n - 1), N, 1, N))).reshape(N**n, N)
    target = io_out @ K_Vbar(n, om, l, eta, th)
    P = np.vdot(target.ravel(), v.ravel())/np.vdot(target.ravel(), target.ravel())
    return P, np.linalg.norm(v - P*target)/np.linalg.norm(v), np.abs(io_in - io_out).max(), mus


d5m = d5q = 0
fail_m = mp.inf
rev = []
for n, l in ((1, 1), (2, 1), (3, 2)):
    om, eta, mu = 0.37, 0.41 + 0.23j, 0.11 + 0.19j
    T = (n + 1)*om/2
    P, res, same, _ = boot(n, om, l, eta, 1j*np.pi*mu/T)
    if n > 1:
        rev.append(len(embedding(n, om, 1j*np.pi*mu/T, order=-1)[0]))
    Sp = SP(n, om, l, mp.mpc(eta.real, eta.imag))
    mm = mp.mpc(mu.real, mu.imag)
    mus = [mm - Sp.T + k*Sp.om for k in range(1, n + 1)]
    bulk = mp.mpc(P.real, P.imag)*mp.fprod([F11(mus[j] + mus[k], n, Sp.om) for j in range(n) for k in range(j + 1, n)])
    Q = {m: Sp.bt(mm, m)/(bulk*mp.fprod([Sp.a(x, m) for x in mus])) for m in (0, 2, -2)}
    d5m = max(d5m, res, same)
    d5q = max(d5q, abs(Q[0] - (-1)**(n + 1)))
    fail_m = min(fail_m, abs(Q[2] - (-1)**(n + 1)), abs(Q[-2] - (-1)**(n + 1)))
report(f'soliton bootstrap V^(x n) -> Vbar (n = 1, 2, 3): M iota = P iota K_Vbar [{d5m:.0e}], Vbar not in the reversed '
       f'order {rev}; bt = (-1)^(n+1) P prod a prod F11 [{mp.nstr(d5q, 2)}]; the e^(+-2 theta) solutions fail '
       f'[{mp.nstr(fail_m, 2)}]', d5m < 1e-12 and d5q < 1e-11 and fail_m > 0.1 and rev == [0, 0])

# ---------------------------------------------------------------- 6. CDD factors and the bootstrap
mus6 = [mp.mpc('0.11', '0.19'), mp.mpc('-0.23', '0.07')]
ok6, sA_even = True, mp.inf
for n in range(1, 6):
    om = mp.mpf('0.37')
    T = (n + 1)*om/2
    cands = {'B': lambda m: -1/(bI(T/2, m, T)**2*bI(T/2 - om/2, m, T)*bI(T/2 + om/2, m, T)),
             'om': lambda m: bI(om, m, T), 'A': lambda m: bI(-T/2, m, T)**2}
    for name, sa in cands.items():
        r = [sa(T - m)/mp.fprod([sa(m - T + k*om) for k in range(1, n + 1)]) for m in mus6]
        const = abs(r[0] - r[1]) < 1e-20 and abs(abs(r[0]) - 1) < 1e-20 and abs(r[0].imag) < 1e-20
        if name == 'A' and n % 2 == 0:
            sA_even = min(sA_even, abs(r[0] - r[1]))
        else:
            ok6 &= const
report(f'CDD and bootstrap: sigma_B and the pair ((omega)_I, -(T-omega)_I) satisfy sigma_b/prod sigma_a = +-1 for n = 1..5, '
       f'sigma_A = (-T/2)_I^2 for odd n; sigma_A fails for n = 2, 4 [ratio varies by {mp.nstr(sA_even, 2)}]',
       ok6 and sA_even > 0.1)


# ---------------------------------------------------------------- 7. breather amplitude
def rvv_diag(n, q, z, a, b):
    N = n + 1
    if a == b:
        return 1
    return (q - 1/q)*z**(N - (b - a) % N)/(q*z**N - 1/q)


def t_of(q, Z):
    return (Z - 1)/(q*Z - 1/q)


def rvw0(n, q, om, th, a, b, lib=np):
    """zero-weight block of Rc_VVbar(th): [aN+a, bN+b] = Rc_VV(i pi - th)[aN+b, aN+b]/t(i pi - th)"""
    zc = lib.exp(om*(1j*lib.pi - th))
    return rvv_diag(n, q, zc, a, b)/t_of(q, zc**(n + 1))


d7a = 0
for n in (2, 3):
    N = n + 1
    for om in (0.37, 1.29):
        th = 0.31 + 0.12j
        R = Rc(n, om, 'V', 'W', np.exp(om*th))
        d7a = max(d7a, max(abs(R[a_*N + a_, b_*N + b_] - rvw0(n, q_phys(om), om, th, a_, b_))
                           for a_ in range(N) for b_ in range(N)))


def vec_mp(n, q, y, conj):
    N = n + 1
    g = {}
    for j in range(N):
        a_, b_ = ((j - 1) % N, j % N) if j else (n, 0)
        E = mp.zeros(N, N)
        E[a_, b_] = 1
        H = [0]*N
        H[a_] += 1
        H[b_] -= 1
        g[('e', j)], g[('f', j)] = y*(E.T if conj else E), (E if conj else E.T)/y
        g[('k', j)] = mp.diag([q**(-x if conj else x) for x in H])
    return g


def singlet_mp(A, B, n, i_norm):
    """zero-weight singlet w = sum_k x_k a_k (x) b_k of A (x) B, annihilated by Delta(e_j) = e_j (x) 1 + k_j (x) e_j and
    Delta(f_j) = f_j (x) k_j^-1 + 1 (x) f_j (Part VI coproduct); x[i_norm] = 1."""
    N = n + 1
    cols = []
    for k in range(N):
        ek = mp.zeros(N, 1)
        ek[k] = 1
        col = []
        for j in range(N):
            for terms in ([(A[('e', j)]*ek, ek), (A[('k', j)]*ek, B[('e', j)]*ek)],
                          [(A[('f', j)]*ek, B[('k', j)]**-1*ek), (ek, B[('f', j)]*ek)]):
                v = [0]*(N*N)
                for u_, w_ in terms:
                    for r in range(N):
                        for s_ in range(N):
                            v[r*N + s_] += u_[r]*w_[s_]
                col += v
        cols.append(col)
    nr = len(cols[0])
    big = mp.matrix(nr, N)
    for k in range(N):
        for r in range(nr):
            big[r, k] = cols[k][r]
    others = [k for k in range(N) if k != i_norm]
    Ared, rhs = mp.matrix(nr, N - 1), mp.matrix(nr, 1)
    for r in range(nr):
        rhs[r] = -big[r, i_norm]
        for c_, k in enumerate(others):
            Ared[r, c_] = big[r, k]
    sol = mp.lu_solve(Ared.H*Ared, Ared.H*rhs)
    x = [mp.mpf(1) if k == i_norm else sol[others.index(k)] for k in range(N)]
    return x, mp.norm(big*mp.matrix(x))


def breather(Sp, mu, bar=False):
    """B_1 -> Bbar_1 (bar = False) or Bbar_1 -> B_1 amplitude, minimal scalar factors (s_ cancels)."""
    n, om, N, T, q = Sp.n, Sp.om, Sp.N, Sp.T, Sp.q
    th = I*PI*mu/T
    u = PI*(1 - 1/T)
    ta, tb = th - I*u/2, th + I*u/2
    kV = lambda t, j: (mp.exp(2*j*om*t)*mp.sinh(Sp.eta + T*t) if j < Sp.l
                       else mp.exp((2*j - N)*om*t)*mp.sinh(Sp.eta - T*t))
    kW = lambda t, j: (mp.exp(-2*j*om*t)*mp.sinh(Sp.etab - T*t) if j < Sp.l
                       else mp.exp(-(2*j - N)*om*t)*mp.sinh(Sp.etab + T*t))
    hw = lambda conj: n if conj else 0
    x, r1 = singlet_mp(vec_mp(n, q, mp.exp(om*ta), bar), vec_mp(n, q, mp.exp(om*tb), not bar), n, hw(bar))
    y, r2 = singlet_mp(vec_mp(n, q, mp.exp(-om*tb), not bar), vec_mp(n, q, mp.exp(-om*ta), bar), n, hw(not bar))
    kb = [(kV if bar else kW)(tb, j) for j in range(N)]
    ka = [(kW if bar else kV)(ta, j) for j in range(N)]
    zc = mp.exp(om*(I*PI - (ta + tb)))
    vals = []
    for a_ in range(N):
        if not bar:
            s_ = sum(rvw0(n, q, om, ta + tb, a_, b_, mp)*kb[b_]*x[b_] for b_ in range(N))
        else:
            s_ = sum(rvv_diag(n, q, zc, b_, a_)/t_of(q, zc**N)*kb[b_]*x[b_] for b_ in range(N))
        vals.append(ka[a_]*s_/y[a_])
    spread = max(abs(v/vals[0] - 1) for v in vals) + r1 + r2
    ma, mb = mu - (T - 1)/2, mu + (T - 1)/2
    amp = (Sp.bt(ma) if bar else Sp.a(ma))*Sp.F(2*mu)*(Sp.a(mb) if bar else Sp.bt(mb))*vals[0]
    return amp, spread


def Eeta(Sp, mu, bar=False):
    T = Sp.T
    th = I*PI*mu/T
    A = Sp.eps + (T - 1)/2
    Ab = A + Sp.om*Sp.l - T
    if bar:
        A, Ab = Ab, A
    return (mp.sinh(th) - I*mp.sin(PI*A/T))/(mp.sinh(th) + I*mp.sin(PI*Ab/T))


d7 = d7s = 0
for n, l in ((1, 1), (2, 0), (2, 1), (3, 2)):
    Sp = SP(n, '0.37', l, mp.mpc('-0.3', '1.1'))
    h = n + 1
    for mu in (mp.mpc('0.11', '0.19'), mp.mpc('-0.23', '0.07')):
        for bar in (False, True):
            amp, spread = breather(Sp, mu, bar)
            want = mp.exp((1 if bar else -1)*2*PI*I*l/h)*K0(mu, Sp.T, Sp.om)*Eeta(Sp, mu, bar)
            d7 = max(d7, abs(amp/want - 1))
            d7s = max(d7s, spread)
# n = 1 against Ghoshal at Dirichlet
d7g = 0
for eta in (mp.mpc('0.41', '0.23'), mp.mpc('-0.3', '1.1')):
    Sp = SP(1, '0.37', 1, eta)
    T, lam = Sp.T, Sp.om
    xi = PI*Sp.eps - PI/2
    for mu in (mp.mpc('0.11', '0.19'), mp.mpc('-0.23', '0.07')):
        u = PI*mu/lam
        gh = (mp.cos(u/2 + PI/(4*lam))*mp.cos(u/2 - PI/4 - PI/(4*lam))*mp.sin(u/2 + PI/4)
              / (mp.cos(u/2 - PI/(4*lam))*mp.cos(u/2 + PI/4 + PI/(4*lam))*mp.sin(u/2 - PI/4))
              * (mp.cos(xi/lam) - mp.sin(u))/(mp.cos(xi/lam) + mp.sin(u)))
        sA = lambda m: bI(-T/2, m, T)**2
        d7g = max(d7g, abs(breather(Sp, mu)[0]*sA(mu - (T - 1)/2)*sA(mu + (T - 1)/2)/gh - 1))
report(f'breather: zero-weight block of Rc_VVbar in closed form [{d7a:.0e}]; B_1 -> Bbar_1 = e^(-2 pi i l/h) K0 E, '
       f'Bbar_1 -> B_1 = e^(2 pi i l/h) K0 Ebar ((n, l) = (1,1), (2,0), (2,1), (3,2)) [{mp.nstr(d7, 2)}; singlet and '
       f'projection {mp.nstr(d7s, 2)}]; n = 1 with sigma_A: Ghoshal K_B1 at Dirichlet [{mp.nstr(d7g, 2)}]',
       d7a < 1e-12 and d7 < mp.mpf('1e-24') and d7s < mp.mpf('1e-24') and d7g < mp.mpf('1e-24'))

# ---------------------------------------------------------------- 8. continuation and particle axioms
d8 = 0
THS = [mp.mpc('0.37', '0.11'), mp.mpc('-0.83', '0.29')]
for n in range(1, 8):
    h = n + 1
    H = mp.mpf(h)
    for B in (mp.mpf('0.317'), mp.mpf('1.43')):
        T = h*(-2/B)/2
        for th in THS:
            book = (blk(H/2, th, h)*blk(H/2 + B/2, th, h)*blk(H/2 - 1, th, h)*blk(H/2 - 1 - B/2, th, h)
                    * blk(H/2 + 1 + B/2, th, h)/blk(H/2 - 1 + B/2, th, h))
            d8 = max(d8, abs(K0(-I*T*th/PI, T, -2/B)/book - 1))
SIG = {'none': lambda m, T, om: 1, 'A': lambda m, T, om: bI(-T/2, m, T)**2,
       'B': lambda m, T, om: -1/(bI(T/2, m, T)**2*bI(T/2 - om/2, m, T)*bI(T/2 + om/2, m, T))}


def particle_L(n, B, l, phi, cdd):
    h, om = n + 1, -2/B
    T = h*om/2
    s2 = mp.sin(phi + 2*PI*l/h - PI)
    sg = SIG[cdd]

    def L1(th, bar=False):
        mu = -I*T*th/PI
        if not bar:
            E = (mp.sinh(th) - I*mp.sin(phi))/(mp.sinh(th) + I*s2)
        else:
            E = (mp.sinh(th) - I*s2)/(mp.sinh(th) + I*mp.sin(phi))
        return K0(mu, T, om)*E*sg(mu - (T - 1)/2, T, om)*sg(mu + (T - 1)/2, T, om)

    def L(a):
        if a == 1:
            return L1
        prev = L(a - 1)
        return lambda th: prev(th + I*PI/h)*S_ab(a - 1, h - 1, n)(2*th + I*PI*(2 - a)/h, B)*L1(th - I*PI*(a - 1)/h)
    return L, L1


def axioms(n, B, l, cdd):
    h = n + 1
    L, L1 = particle_L(n, B, l, mp.mpc('0.7', '0.2'), cdd)
    Ls = {a: L(a) for a in range(1, h)}
    w = {'Ln': 0, 'unit': 0, 'cross': 0, 'boot': 0}
    for th in THS:
        w['Ln'] = max(w['Ln'], abs(Ls[n](th)/L1(th, bar=True) - 1))
        for a in range(1, h):
            w['unit'] = max(w['unit'], abs(Ls[a](th)*Ls[h - a](-th) - 1))
            w['cross'] = max(w['cross'], abs(Ls[a](I*PI/2 - th)/(S_ab(a, a, n)(2*th, B)*Ls[a](I*PI/2 + th)) - 1))
            for b in range(1, h):
                if a + b == h:
                    continue
                c, ua, ub = (a + b, PI*b/h, PI*a/h) if a + b < h else (a + b - h, PI*(h - b)/h, PI*(h - a)/h)
                rhs = Ls[a](th + I*ua)*S_ab(a, h - b, n)(2*th + I*ua - I*ub, B)*Ls[b](th - I*ub)
                w['boot'] = max(w['boot'], abs(Ls[c](th)/rhs - 1))
    return w


ok8, worst8 = True, 0
for cdd in ('none', 'B'):
    for n in (2, 3):
        for l in (0, 1):
            w = axioms(n, mp.mpf('0.317'), l, cdd)
            worst8 = max(worst8, max(w.values()))
wA = axioms(2, mp.mpf('0.317'), 0, 'A')
d8lim = 0
for n in (2, 3):
    h = n + 1
    L1 = particle_L(n, mp.mpf('1e-15'), 1, mp.mpc('0.7', '0.2'), 'B')[1]
    g = lambda x, a: (x + I*mp.sin(PI*a/h))/(x - I*mp.sin(PI*a/h))
    for th in THS:
        x = mp.sinh(th)
        E0 = (x - I*mp.sin(mp.mpc('0.7', '0.2')))/(x + I*mp.sin(mp.mpc('0.7', '0.2') + 2*PI/h - PI))
        d8lim = max(d8lim, abs(L1(th)/(g(x, mp.mpf(h)/2)*g(x, mp.mpf(h)/2 - 1)*E0) - 1))
report(f'continuation: K0 -> (h/2)(h/2+B/2)(h/2-1)(h/2-1-B/2)(h/2+1+B/2)/(h/2-1+B/2), n = 1..7 [{mp.nstr(d8, 2)}]; '
       f'particle amplitudes (minimal and sigma_B; n = 2, 3; l = 0, 1): L_n, unitarity, crossing with S_aa, bootstrap '
       f'[{mp.nstr(worst8, 2)}]; sigma_A breaks the bootstrap for n = 2 [{mp.nstr(wA["boot"], 2)}]; '
       f'B -> 0: L_1 -> g_(h/2) g_(h/2-1) E [{mp.nstr(d8lim, 2)}]',
       d8 < mp.mpf('1e-24') and worst8 < mp.mpf('1e-24') and wA['boot'] > 0.1 and wA['cross'] < 1e-24
       and d8lim < mp.mpf('1e-10'))
