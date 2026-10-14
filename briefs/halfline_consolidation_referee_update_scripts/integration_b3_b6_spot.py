"""Spot checks, written from the book's formulas.
B6: duals of the s = -1 amplitudes, -(a)(h-a) K^N_a (eq-neumann-amplitude), equal CDRS (5.16)
    K_a^0 = prod_{k=1}^a <k-1/2>/<h-k+1/2>, <x> = (x+1/2)/(x-1/2+B/2)   (hep-th/9404108, (5.13), (5.16)),
    and for n = 2 their (5.7) -(1)(2+B/2)/(B/2).
B3: prp-sp-scalar, the first crossing equation of eq-sp-scalar-eqs, n = 2, l = 1:
    -(-q)^{1-l} s(mu) a(T/2-mu) = F(2mu) b~(T/2+mu), F(x) = S^T(T-x), S^T(mu) = sin(pi mu)/sin(pi(omega-mu)) F_11(mu),
    F_11 of prp-an-diagonal-scalar; infinite products by mpmath nprod (Richardson extrapolation).
    Unitarity a(mu)a(-mu)U(mu) = 1 holds factor by factor (checked by hand).
"""
import mpmath as mp
mp.mp.dps = 25
I = mp.mpc(0, 1)


def blk(x, th, h):
    return mp.sinh((th + I*mp.pi*x/h)/2)/mp.sinh((th - I*mp.pi*x/h)/2)


# ---------------------------------------------------------------- B6
worst = 0
for n in (2, 3, 4, 5, 6, 7):
    h = n + 1
    for B in (mp.mpf('0.37'), mp.mpc('1.2', '0.3')):
        for th in (mp.mpc('0.41', '0.17'), mp.mpc('-1.1', '0.6')):
            for a in range(1, n+1):
                KN = mp.fprod(blk(c-1, th, h)*blk(h-c+B/2, th, h)/(blk(h-c, th, h)*blk(c-1+B/2, th, h)) for c in range(1, a+1))
                dual = -blk(a, th, h)*blk(h-a, th, h)*KN
                ang = lambda x: blk(x+mp.mpf(1)/2, th, h)/blk(x-mp.mpf(1)/2+B/2, th, h)
                cdrs = mp.fprod(ang(k-mp.mpf(1)/2)/ang(h-k+mp.mpf(1)/2) for k in range(1, a+1))
                worst = max(worst, abs(dual/cdrs - 1))
                if n == 2 and a == 1:
                    worst = max(worst, abs(dual/(-blk(1, th, h)*blk(2+B/2, th, h)/blk(B/2, th, h)) - 1))
print(f'B6: -(a)(h-a) K^N_a = CDRS (5.16) (and (5.7) for n = 2), n = 2..7, all a: {mp.nstr(worst, 3)}')

# ---------------------------------------------------------------- B3
n, l = 2, 1
h = n + 1
om = mp.mpf('0.41')
T = h*om/2
q = -mp.exp(-I*mp.pi*om)
eta = mp.mpc('0.3', '0.7')
eps = eta/(I*mp.pi)
epsb = eps + om*l - T


def G(j, x, c, step=1):
    return mp.gamma(x + j*step*h*om + c)/mp.gamma(-x + j*step*h*om + c)


def logG(j, x, c):
    return mp.loggamma(x + j*h*om + c) - mp.loggamma(-x + j*h*om + c)


def calG(x, c, cp, step=1):
    """G_{c,c'}(x) = G_0(x,c)/G_0(x,c') prod_{j>=1} G_j(x,c)G_j(x,1-c)/(G_j(x,c')G_j(x,1-c')), G_j with j*step."""
    head = G(0, x, c)/G(0, x, cp)
    tail = mp.nsum(lambda j: (mp.loggamma(x + j*step*h*om + c) - mp.loggamma(-x + j*step*h*om + c)
                              + mp.loggamma(x + j*step*h*om + 1 - c) - mp.loggamma(-x + j*step*h*om + 1 - c)
                              - mp.loggamma(x + j*step*h*om + cp) + mp.loggamma(-x + j*step*h*om + cp)
                              - mp.loggamma(x + j*step*h*om + 1 - cp) + mp.loggamma(-x + j*step*h*om + 1 - cp)),
                   [1, mp.inf])
    return head*mp.exp(tail)


def blkI(y, mu):
    return mp.sin(mp.pi*(mu + y)/(2*T))/mp.sin(mp.pi*(mu - y)/(2*T))


def theta_of(mu):
    return I*mp.pi*mu/T


def rho(mu):
    return blkI(om/2, mu)*calG(2*mu, 2*T + om, 0, step=2)


def a_sp(mu):
    return -I*rho(mu)*mp.exp((1-l)*om*theta_of(mu))*calG(mu, eps + 2*T, epsb + T)/mp.sin(mp.pi*(eps + mu))


def b_sp(mu):
    return -I*rho(mu)*mp.exp(-(1-l)*om*theta_of(mu))*calG(mu, epsb + 2*T, eps + T)/mp.sin(mp.pi*(epsb + mu))


def F11(mu):
    def term(j):
        return (mp.loggamma(mu + j*h*om - (h-1)*om + 1) + mp.loggamma(mu + j*h*om - om) + mp.loggamma(-mu + j*h*om - h*om + 1)
                + mp.loggamma(-mu + j*h*om) - mp.loggamma(-mu + j*h*om - (h-1)*om + 1) - mp.loggamma(-mu + j*h*om - om)
                - mp.loggamma(mu + j*h*om - h*om + 1) - mp.loggamma(mu + j*h*om))
    return -mp.exp(mp.nsum(term, [1, mp.inf]))


def ST(mu):
    return mp.sin(mp.pi*mu)/mp.sin(mp.pi*(om - mu))*F11(mu)


s = lambda mu: mp.sin(mp.pi*(2*mu - T))/mp.sin(mp.pi*(2*mu + T))
Fx = lambda x: ST(T - x)
res = 0
for mu in (mp.mpc('0.13', '0.07'), mp.mpc('-0.21', '0.11')):
    lhs1 = -(-q)**(1-l)*s(mu)*a_sp(T/2 - mu)
    rhs1 = Fx(2*mu)*b_sp(T/2 + mu)
    lhs2 = -(-q)**(l-1)*s(mu)*b_sp(T/2 - mu)
    rhs2 = Fx(2*mu)*a_sp(T/2 + mu)
    res = max(res, abs(lhs1/rhs1 - 1), abs(lhs2/rhs2 - 1))
    U = mp.sinh(eta + I*mp.pi*mu)*mp.sinh(eta - I*mp.pi*mu)
    res = max(res, abs(a_sp(mu)*a_sp(-mu)*U - 1))
print(f'B3: eq-sp-scalar-eqs (unitarity of a, both crossing equations), n = 2, l = 1: {mp.nstr(res, 3)}')
