"""d_n vector soliton scalar factor F (normalized to 1 on L(2w1), identical solitons) and the eps = 0 boundary
scalar factor A.  Variables: T = (n-1) om (= h om/2, h = 2n-2), mu = -i T theta/pi (= t of Part VI)."""
import numpy as np
from scipy.special import loggamma as lg


def rich(f, K):
    return 2*f(2*K) - f(K)


def blk(a, mu, T):
    """DG block (a)_I = sin(pi(mu+a)/2T)/sin(pi(mu-a)/2T)."""
    return np.sin(np.pi*(mu + a)/(2*T))/np.sin(np.pi*(mu - a)/(2*T))


def _L(mu, cs, step, K):
    k = np.arange(K)[:, None]
    c = np.array(cs)[None, :] + step*k
    return lg(mu + c) - lg(-mu + c)


def logF(mu, n, om, K=4000):
    T = (n - 1)*om
    num = [1 + T, 1 + om, T, 2*T - om]
    den = [1 + 2*T, 1 + om + T, 0, T - om]
    def part(KK):
        return np.sum(_L(mu, num, 2*T, KK)) - np.sum(_L(mu, den, 2*T, KK))
    return rich(part, K)


def F(mu, n, om, K=4000):
    """F = B_{1+T} B_{1+om} / (B_0 B_{T-om}),  B_a = prod_k g(a+2Tk)/g(a+T+2Tk),  g(c) = G(mu+c)/G(-mu+c)."""
    return np.exp(logF(mu, n, om, K))


def logA(mu, n, om, K=3000):
    T = (n - 1)*om
    num = [1 + 3*T, 1 + om + 2*T, T, 2*T - om]
    den = [1 + 2*T, 1 + om + T, 2*T, 3*T - om]
    def part(KK):
        return np.sum(_L(2*mu, num, 4*T, KK)) - np.sum(_L(2*mu, den, 4*T, KK))
    return rich(part, K)


def A(mu, n, om, K=3000):
    """eps = 0 boundary scalar factor (minimal; simple zero at mu = T/2, A(0) = 1)."""
    return np.exp(logA(mu, n, om, K))


def r_cross(mu, n, om):
    """F(T - mu)/F(mu) from the crossing factor c_u of the d_n vector R-matrix."""
    T = (n - 1)*om
    s = lambda z: np.sin(np.pi*z)
    return s(mu)*s(mu - T + om)/(s(mu - om)*s(mu - T))


def KB_blocks(mu, n, om):
    """conjectured closed form of A(mu_a) F(2mu) A(mu_b)."""
    T = (n - 1)*om
    xs = [0.5 + T, 0.5 + om/2 + T/2, T/2 - om/2, om/2 + T, 1.5*T - 0.5, 1.5*T, 2*T - om/2 - 0.5]
    return np.prod([blk(x, mu, T) for x in xs], axis=0)


def KB_num(mu, n, om):
    T = (n - 1)*om
    return A(mu - T/2 + 0.5, n, om)*F(2*mu, n, om)*A(mu + T/2 - 0.5, n, om)
