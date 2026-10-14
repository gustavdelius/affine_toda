"""O.4: leading semiclassical phase derivatives of the d_4 scalar factors (ch. 35, end of sec-dn-uniform), from the
formulas of prp-dn-scalar (eq-dn-scalar, F) and prp-dn-nondiagonal-scalar (rho, g), written independently.
d_4: n = 4, h = 6, T = 3 omega, mu = -i T theta/pi, g(mu, c) = Gamma(mu+c)/Gamma(-mu+c), G_j(mu, a) = g(mu, 2jT + a).
Prints d/dtheta(-i log X)/omega at theta = 0.83 for X = A_+ = A rho_+ g_+, A_- = A rho_- g_-, and F(2 mu) (halved),
at omega = 40.37, 80.37, 160.37, and the first-order (1/omega) extrapolations.
"""
import mpmath as mp
mp.mp.dps = 30
I = mp.mpc(0, 1)


def lg(mu, c):
    return mp.loggamma(mu + c) - mp.loggamma(-mu + c)


def logA(mu, T, om):
    def term(k):
        x = 2*mu
        return (lg(x, 1 + 3*T + 4*T*k) + lg(x, 1 + om + 2*T + 4*T*k) + lg(x, T + 4*T*k) + lg(x, 2*T - om + 4*T*k)
                - lg(x, 1 + 2*T + 4*T*k) - lg(x, 1 + om + T + 4*T*k) - lg(x, 2*T + 4*T*k) - lg(x, 3*T - om + 4*T*k))
    return mp.nsum(term, [0, mp.inf])


def logF(mu, T, om):
    def term(k):
        Bk = lambda a: lg(mu, a + 2*T*k) - lg(mu, a + T + 2*T*k)
        return Bk(1 + T) + Bk(1 + om) - Bk(0) - Bk(T - om)
    return mp.nsum(term, [0, mp.inf])


def blkI(y, mu, T):
    return mp.sin(mp.pi*(mu + y)/(2*T))/mp.sin(mp.pi*(mu - y)/(2*T))


def logrho(s, mu, T, om):
    d = 0 if s == 1 else mp.mpf(1)/2
    if s == 1:
        r = mp.log(blkI(-T/2, mu, T)*mp.sin(mp.pi*om/2)/mp.sin(mp.pi*(mu - T/2)))
    else:
        r = mp.log(-I*mp.sin(mp.pi*om/2)/mp.cos(mp.pi*(mu - T/2)))
    def term(j):
        G = lambda a: lg(mu, 2*j*T + a)
        return G(d + om/2 - 3*T/2) + G(1 - d - om/2 - T/2) - G(d + om/2 - T/2) - G(1 - d - om/2 - 3*T/2)
    return r + mp.nsum(term, [1, mp.inf])


def logv(s, x, om):
    if s == 1:
        return mp.loggamma(x + om/2) + mp.loggamma(1 + x - om/2) - mp.loggamma(x) - mp.loggamma(1 + x)
    h = mp.mpf(1)/2
    return mp.loggamma(h + x + om/2) + mp.loggamma(h + x - om/2) - 2*mp.loggamma(h + x)


def logg(s, mu, T, om):
    return logv(s, mu, om) + mp.nsum(lambda k: (-1)**int(k)*(logv(s, mu + k*T, om) - logv(s, -mu + k*T, om)), [1, mp.inf])


def deriv(f, th, T, eps=mp.mpf('1e-7')):
    mu = lambda t: -I*T*t/mp.pi
    return mp.re(-I*(f(mu(th + eps)) - f(mu(th - eps)))/(2*eps))


th = mp.mpf('0.83')
rows = {}
for om in (mp.mpf('40.37'), mp.mpf('80.37'), mp.mpf('160.37')):
    T = 3*om
    dp = deriv(lambda m: logA(m, T, om) + logrho(1, m, T, om) + logg(1, m, T, om), th, T)/om
    dm = deriv(lambda m: logA(m, T, om) + logrho(-1, m, T, om) + logg(-1, m, T, om), th, T)/om
    dF = deriv(lambda m: logF(2*m, T, om), th, T)/om/2
    rows[om] = (dp, dm, dF)
    print(f'omega = {om}: plus {mp.nstr(dp, 9)}, minus {mp.nstr(dm, 9)}, half bulk {mp.nstr(dF, 9)}', flush=True)
oms = sorted(rows)
ext = [2*rows[oms[2]][i] - rows[oms[1]][i] for i in range(3)]
print('first-order extrapolation from 80.37, 160.37:', [mp.nstr(x, 6) for x in ext],
      ' differences to half bulk:', mp.nstr(abs(ext[0] - ext[2]), 3), mp.nstr(abs(ext[1] - ext[2]), 3))
