"""A2: 30-digit scalar factors for the soliton-preserving a_n^(1) matrices (written from scratch).

Ehat_P(z) = Gamma(z) prod_{k>=1} Gamma(z+kP)/Gamma(1-z+kP) (kP)^(1-2z)     [Ehat(z+P)/Ehat(z) = C_P sin(pi z)]
Lam_P(x; c) = Ehat_P(c+x)/Ehat_P(c-x)                                        [Lam(x)Lam(-x) = 1]
Infinite products: exact head k < J, Stirling tail summed with Hurwitz zeta.
"""
import mpmath as mp

mp.mp.dps = 30
J0, M0 = 60, 40


def _tail(P, terms, J=J0, M=M0, kmin=1):
    """sum_{j>=J} sum_i c_i [logGamma(jP + a_i) - (Stirling leading terms)] for sum c_i = 0."""
    s = 0
    for k in range(kmin, M + 1):
        coef = sum(c*mp.bernpoly(k + 1, a) for c, a in terms)
        s += (-1)**(k + 1)*coef/(k*(k + 1)*P**k)*mp.zeta(k, J)
    return s


def logprod(P, terms, J=J0, M=M0):
    """log prod_{j>=1} prod_i Gamma(jP + a_i)^{c_i}, balanced: sum c_i = sum c_i a_i = sum c_i B_2(a_i) = 0."""
    assert abs(sum(c for c, a in terms)) < 1e-25 and abs(sum(c*a for c, a in terms)) < 1e-20
    assert abs(sum(c*mp.bernpoly(2, a) for c, a in terms)) < 1e-20
    head = sum(c*mp.loggamma(j*P + a) for j in range(1, J) for c, a in terms)
    return head + _tail(P, terms, J, M, kmin=2)


def logEhat(z, P, J=J0, M=M0):
    head = mp.loggamma(z) + sum(mp.loggamma(z + k*P) - mp.loggamma(1 - z + k*P) + (1 - 2*z)*mp.log(k*P)
                                for k in range(1, J))
    return head + _tail(P, [(1, z), (-1, 1 - z)], J, M, kmin=2)


def Lam(x, c, P):
    return mp.exp(logEhat(c + x, P) - logEhat(c - x, P))


def blk(x, mu, T):
    """(x)_I = sin(pi(mu+x)/2T)/sin(pi(mu-x)/2T)."""
    return mp.sin(mp.pi*(mu + x)/(2*T))/mp.sin(mp.pi*(mu - x)/(2*T))


def F11(m, n, om):
    """prp-an-diagonal-scalar (book formula), P = h omega."""
    h = n + 1
    L = h*om
    t = [(1, m - (h - 1)*om + 1), (1, m - om), (1, -m - L + 1), (1, -m),
         (-1, -m - (h - 1)*om + 1), (-1, -m - om), (-1, m - L + 1), (-1, m)]
    return -mp.exp(logprod(L, t))


def ST(x, n, om):
    return mp.sin(mp.pi*x)/mp.sin(mp.pi*(om - x))*F11(x, n, om)


class SP:
    """Soliton-preserving scalar data for (n, omega, l, eta).  mu = -i T theta/pi.
    Physical amplitudes: a(mu) K_V(theta) on the soliton, bt(mu) K_Vbar(theta) on the antisoliton."""

    def __init__(self, n, om, l, eta):
        self.n, self.om, self.l = n, mp.mpf(om), l
        self.N = n + 1
        self.T = self.N*self.om/2
        self.eta = mp.mpmathify(eta)
        self.etab = self.eta + 1j*mp.pi*self.om*(l - mp.mpf(self.N)/2)
        self.eps = self.eta/(1j*mp.pi)
        self.epsb = self.etab/(1j*mp.pi)              # = eps + omega l - T
        self.q = -mp.exp(-1j*mp.pi*self.om)
        self.loglam = -1j*mp.pi*self.om*(1 - l)       # lambda = (-q)^(1-l) on the principal branch
        self.lam = mp.exp(self.loglam)

    # ---- equations
    def U(self, mu):
        return mp.sinh(self.eta + 1j*mp.pi*mu)*mp.sinh(self.eta - 1j*mp.pi*mu)

    def Ub(self, mu):
        return mp.sinh(self.etab + 1j*mp.pi*mu)*mp.sinh(self.etab - 1j*mp.pi*mu)

    def s(self, mu):
        T = self.T
        return mp.sin(mp.pi*(2*mu - T))/mp.sin(mp.pi*(2*mu + T))

    def F(self, x):
        """soliton-antisoliton bulk factor S_{Vbar V}(x) = S^T(T - x)."""
        return ST(self.T - x, self.n, self.om)

    def Phi1(self, mu):
        return -self.lam*self.s(mu)

    def Phi2(self, mu):
        return -self.s(mu)/self.lam

    # ---- solution
    def rho(self, mu):
        T, om = self.T, self.om
        return Lam(2*mu, 2*T + om, 4*T)/Lam(2*mu, 0, 4*T)*blk(om/2, mu, T)

    def alpha0(self, mu):
        T = self.T
        return 1j*Lam(mu, self.eps + 2*T, 2*T)/(Lam(mu, self.epsb + T, 2*T)*mp.sin(mp.pi*(self.eps + mu)))

    def beta0(self, mu):
        T = self.T
        return 1j*Lam(mu, self.epsb + 2*T, 2*T)/(Lam(mu, self.eps + T, 2*T)*mp.sin(mp.pi*(self.epsb + mu)))

    def a(self, mu, m=0):
        return self.rho(mu)*mp.exp(-(self.loglam - 1j*mp.pi*m)*mu/self.T)*self.alpha0(mu)

    def bt(self, mu, m=0):
        return (-1)**m*self.rho(mu)*mp.exp((self.loglam - 1j*mp.pi*m)*mu/self.T)*self.beta0(mu)
