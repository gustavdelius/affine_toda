"""Scalar factors for the d_n vector K-matrices, 30 digits (Stirling/Hurwitz tails as an_nondiagonal_scalar.logprod).
T = (n-1) omega, mu = -i T theta/pi, blocks (x)_I = sin(pi(mu+x)/2T)/sin(pi(mu-x)/2T).
A0: book eq-dn-scalar (eps = 0); F: book prp-dn-r. rho_sig: the a_n non-diagonal/diagonal ratio of prp-an-nondiagonal-scalar
(T general); g_sig: the extra d_n factor, g(mu) g(-mu) = U^an_sig/U^dn_sig, g(T - mu) = g(mu)."""
import mpmath as mp

mp.mp.dps = 30


def logprod(P, terms, j0=0, J=60, M=40):
    """log prod_{j>=j0} prod_i Gamma(jP + a_i)^{c_i} (balanced terms): exact head, Stirling tail (Hurwitz zeta)."""
    head = sum(c*mp.loggamma(j*P + a) for j in range(j0, J) for c, a in terms)
    tail = sum((-1)**(k + 1)*sum(c*mp.bernpoly(k + 1, a) for c, a in terms)/(k*(k + 1)*P**k)*mp.zeta(k, J)
               for k in range(2, M + 1))
    return head + tail


def g2(mu, c, sign=1):
    return [(sign, mu + c), (-sign, -mu + c)]


def blk(x, mu, T):
    return mp.sin(mp.pi*(mu + x)/(2*T))/mp.sin(mp.pi*(mu - x)/(2*T))


class DN:
    def __init__(self, n, om):
        self.n, self.om = n, mp.mpf(om)
        self.T = (n - 1)*self.om

    def F(self, mu):
        T, om = self.T, self.om
        terms = []
        for a, s in ((1 + T, 1), (1 + om, 1), (0, -1), (T - om, -1)):
            terms += g2(mu, a, s) + g2(mu, a + T, -s)
        return mp.exp(logprod(2*T, terms))

    def A0(self, mu):
        T, om = self.T, self.om
        num, den = [1 + 3*T, 1 + om + 2*T, T, 2*T - om], [1 + 2*T, 1 + om + T, 2*T, 3*T - om]
        terms = []
        for a in num:
            terms += g2(2*mu, a, 1)
        for a in den:
            terms += g2(2*mu, a, -1)
        return mp.exp(logprod(4*T, terms))

    def rho(self, mu, sig):
        """a_n ratio, without CDD: sig=+1: sin(pi om/2)/sin pi(mu - T/2) Pi_+(mu); sig=-1: -i sin(pi om/2)/cos pi(mu-T/2) Pi_-."""
        T, om = self.T, self.om
        h = mp.mpf(1)/2
        if sig == 1:
            t = g2(mu, om/2 - 3*T/2) + g2(mu, 1 - om/2 - T/2) + g2(mu, om/2 - T/2, -1) + g2(mu, 1 - om/2 - 3*T/2, -1)
            return mp.sin(mp.pi*om/2)/mp.sin(mp.pi*(mu - T/2))*mp.exp(logprod(2*T, t, j0=1))
        t = (g2(mu, h + om/2 - 3*T/2) + g2(mu, h - om/2 - T/2) + g2(mu, h + om/2 - T/2, -1)
             + g2(mu, h - om/2 - 3*T/2, -1))
        return -1j*mp.sin(mp.pi*om/2)/mp.cos(mp.pi*(mu - T/2))*mp.exp(logprod(2*T, t, j0=1))

    def vterms(self, x, sig, sgn=1):
        om = self.om
        if sig == 1:      # v(x) = G(x+om/2) G(1+x-om/2)/(G(x) G(1+x))
            return [(sgn, x + om/2), (sgn, 1 + x - om/2), (-sgn, x), (-sgn, 1 + x)]
        h = mp.mpf(1)/2   # v(x) = G(1/2+x+om/2) G(1/2+x-om/2)/G(1/2+x)^2
        return [(sgn, h + x + om/2), (sgn, h + x - om/2), (-2*sgn, h + x)]

    def g(self, mu, sig):
        """g(mu) = v(mu) prod_{k>=1} [v(mu+kT)/v(-mu+kT)]^{(-1)^k} = v(mu) prod_{j>=1} v(mu+2jT) v(-mu+(2j-1)T)
        / [v(-mu+2jT) v(mu+(2j-1)T)]."""
        T = self.T
        t = (self.vterms(mu, sig) + self.vterms(-mu - T, sig) + self.vterms(-mu, sig, -1) + self.vterms(mu - T, sig, -1))
        # j-th term: v(mu + 2jT) v(-mu + 2jT - T)/(v(-mu + 2jT) v(mu + 2jT - T)) -> Gamma(2jT + a)
        v0 = mp.exp(sum(c*mp.loggamma(a) for c, a in self.vterms(mu, sig)))
        return v0*mp.exp(logprod(2*T, t, j0=1))

    def chi(self, mu, sig):
        return blk(-self.T/2, mu, self.T) if sig == 1 else 1

    def A(self, mu, sig, cdd=None):
        a = self.A0(mu)*self.rho(mu, sig)*self.g(mu, sig)*self.chi(mu, sig)
        if cdd is not None:
            a *= cdd(mu)
        return a

    # matrix-side data (closed forms verified in t8, t9, t14)
    def Uhat(self, mu, sig):
        T, om = self.T, self.om
        c2 = mp.cos(2*mp.pi*mu)
        return sig*(c2 - sig*mp.cos(mp.pi*T))*(c2 - sig*mp.cos(mp.pi*om))/(2*mp.sin(mp.pi*om/2)**2*(c2 - sig))

    def Psi(self, mu, sig):
        om = self.om
        if sig == 1:
            return mp.sin(mp.pi*(om/2 - mu))/mp.sin(mp.pi*(om/2 + mu))
        return mp.cos(mp.pi*(mu - om/2))/mp.cos(mp.pi*(mu + om/2))

    def csing(self, mu, sig):
        T, om = self.T, self.om
        S = lambda z: mp.sin(mp.pi*z)
        m = mu - mp.mpf(1)/2 if sig == 1 else mu
        cf = S(m - T/2 + om/2)*S(m - T/2 - om/2)*S(m - T)*S(m + om/2)*S(m)/(S(m - T/2)**2*S(m - om/2))
        return -sig*cf/mp.sin(mp.pi*om/2)**2

    def breather(self, mu, sig, cdd=None):
        h = mp.mpf(1)/2
        T = self.T
        return self.A(mu - T/2 + h, sig, cdd)*self.F(2*mu)*self.A(mu + T/2 - h, sig, cdd)*self.csing(mu, sig)
