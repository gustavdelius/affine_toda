"""S_0 of eq-sg-s0: direct integral, and a Gamma-function product valid on the whole plane.

From F(t) = sinh((pi-xi)t/2)/(sinh(xi t/2)cosh(pi t/2)) = 2(e^{-xi t}-e^{-pi t}) sum_{m,n>=0} (-1)^n e^{-(m xi+n pi)t}
and int_0^inf dt/t e^{-pt} sin(theta t) = arctan(theta/p), summing over m in closed form:
  S_0(theta) = - prod_{n>=0} [ G(a_n+z) G(b_n-z) / (G(a_n-z) G(b_n+z)) ]^{(-1)^n},
  a_n = 1 + n*lam, b_n = (n+1)*lam, z = i*lam*theta/pi.
"""
from mpmath import mp, mpf, mpc, quad, sinh, cosh, sin, exp, pi, inf, loggamma, nsum
mp.dps = 20

def S0_direct(th, xi):
    f = lambda t: sinh(t*(pi - xi)/2)*sin(th*t)/(t*sinh(xi*t/2)*cosh(pi*t/2))
    return -exp(-1j*quad(f, [0, 1, 10, 40, inf]))

def S0(th, lam):
    lam = mpf(lam); z = 1j*lam*th/pi
    def L(n):
        a = 1 + n*lam; b = (n+1)*lam
        return loggamma(a+z) + loggamma(b-z) - loggamma(a-z) - loggamma(b+z)
    # pair consecutive n for absolute convergence
    s = nsum(lambda k: L(2*k) - L(2*k+1), [0, inf])
    return -exp(s)

def ST(th, lam): return sinh(lam*th)/sinh(lam*(1j*pi - th))*S0(th, lam)
def SR(th, lam): return sinh(1j*pi*lam)/sinh(lam*(1j*pi - th))*S0(th, lam)

if __name__ == '__main__':
    for lam in (mpf('0.4'), mpf('0.7'), mpf('1.5'), mpf('2.5')):
        xi = pi/lam
        d = max(abs(S0_direct(th, xi) - S0(th, lam)) for th in (mpc('0.3', '0.2'), mpc('-1.1', '0.4')))
        print(f'lam={lam}: |S0 integral - S0 Gamma product| = {mp.nstr(d,3)}')
    for lam in (mpf('0.7'), mpf('0.4')):
        xi = pi/lam; worst = 0
        for xx in (mpf('0.3'), mpf('1.1'), mpf('-0.6')):
            th = 1j*pi/2 + xx
            lhs = S0_direct(1j*pi/2 - xx, xi)/S0_direct(th, xi)
            rhs = sinh(lam*th)/sinh(lam*(1j*pi - th))
            worst = max(worst, abs(lhs - rhs))
        print(f'crossing S0(i pi - th) = S_T(th), direct integral at th = i pi/2 + x, lam={lam}: max err {mp.nstr(worst,3)}')
    for lam in (mpf('0.4'), mpf('0.7'), mpf('1.2'), mpf('2.5'), mpf('3.7')):
        worst = 0; uni = 0
        for th in (mpc('0.37', '0.21'), mpc('-0.8', '2.3'), mpc('1.3', '-0.9')):
            lhs = S0(1j*pi - th, lam); rhs = ST(th, lam)
            worst = max(worst, abs(lhs - rhs)/abs(rhs))
            uni = max(uni, abs(S0(th, lam)*S0(-th, lam) - 1))
            # S_R crossing to itself, and 2x2 unitarity S_T(th)S_T(-th)+S_R(th)S_R(-th)=1, S_T(th)S_R(-th)+S_R(th)S_T(-th)=0
            u1 = abs(ST(th, lam)*ST(-th, lam) + SR(th, lam)*SR(-th, lam) - 1)
            u2 = abs(ST(th, lam)*SR(-th, lam) + SR(th, lam)*ST(-th, lam))
            c3 = abs(SR(1j*pi - th, lam) - SR(th, lam))
            uni = max(uni, u1, u2, c3)
        print(f'Gamma product, lam={lam}: crossing S0<->S_T rel err {mp.nstr(worst,3)}; unitarity & S_R crossing err {mp.nstr(uni,3)}')
    # convergence strip of the integral for xi > pi: integrand ~ e^{(|Im th| - pi) t}, not e^{(|Im th|-xi)t}
    lam = mpf('0.7'); xi = pi/lam
    for y in (mpf('3.0'), mpf('3.3'), mpf('4.0')):
        F = lambda t: sinh(t*(pi - xi)/2)/(t*sinh(xi*t/2)*cosh(pi*t/2))
        print(f'lam=0.7 (xi={mp.nstr(xi,4)}), Im theta={y}: |integrand| at t=20,40,60:',
              [mp.nstr(abs(F(t)*sinh(y*t)), 3) for t in (20, 40, 60)])
