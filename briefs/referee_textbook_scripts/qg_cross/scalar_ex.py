"""exr-scalar-unitarity and exr-lift-cdd (chapter 21): the product f of sec-scalar-factor, as printed.
Example: N = 1 zero a, T arbitrary.  Factor_k(t) = f1(t+2T(k-1)) f1(-t+2T(k-1/2)) / [f1(t+2T(k-1/2)) f1(-t+2Tk)]."""
from mpmath import mp, mpf, mpc, loggamma, log, pi, sin, sinh, exp, nsum, inf
mp.dps = 30
T, a = mpf('2.37'), mpf('0.83')
def lf1(z, A): return sum(loggamma(z - aj) + loggamma(1 + z + aj) for aj in A) - len(A)*log(pi)
def logfac(t, k, A): return lf1(t + 2*T*(k-1), A) + lf1(-t + 2*T*(k-mpf(1)/2), A) - lf1(t + 2*T*(k-mpf(1)/2), A) - lf1(-t + 2*T*k, A)
t = mpc('0.31', '0.17')
print('log Factor_k(t) for k = 10, 100, 1000:', [mp.nstr(logfac(t, k, [a]).real, 8) for k in (10, 100, 1000)])
print('compare -4 N T log(2Tk)            :', [mp.nstr(-4*T*log(2*T*k), 8) for k in (10, 100, 1000)],
      ' -> printed product diverges/vanishes; t-dependence of each factor -> 0:',
      mp.nstr(abs(logfac(t, 1000, [a]) - logfac(mpc('0.7', '-0.2'), 1000, [a])), 3))
# normalised product: divide each factor by its value at a fixed t0
t0 = mpf('0.1234')
def logf(t, A, K=4000):
    s = sum(logfac(t, k, A) - logfac(t0, k, A) for k in range(1, K))
    return s
def c(t, A): return mp.fprod(sin(pi*(t - aj)) for aj in A)
for tt in (mpc('0.31', '0.17'), mpc('1.1', '-0.4')):
    lhs = c(tt, [a])*c(-tt, [a])*exp(logf(tt, [a]) + logf(-tt, [a]))
    rhs = c(mpc('0.5', '0.3'), [a])*c(-mpc('0.5', '0.3'), [a])*exp(logf(mpc('0.5', '0.3'), [a]) + logf(-mpc('0.5', '0.3'), [a]))
    print(f't={tt}: f(T-t)/f(t) = {mp.nstr(exp(logf(T - tt, [a]) - logf(tt, [a])), 10)};  c(t)c(-t)f(t)f(-t) / same at t=0.5+0.3i = {mp.nstr(lhs/rhs, 10)}')
# reflection formula: f1(t) f1(-t) c(t) c(-t) = 1 exactly
print('f1(t)f1(-t)c(t)c(-t) =', mp.nstr(exp(lf1(t, [a]) + lf1(-t, [a]))*c(t, [a])*c(-t, [a]), 15))
# exr-lift-cdd: F_{a+1}/F_a versus the block product (a+1)(T-a-1), (y) = sinh(th/2 + i pi y/2T)/sinh(th/2 - i pi y/2T), th = i pi t/T
blk = lambda y, th: sinh(th/2 + 1j*pi*y/(2*T))/sinh(th/2 - 1j*pi*y/(2*T))
vals = []
for tt in (mpc('0.31', '0.17'), mpc('1.1', '-0.4'), mpc('-0.6', '0.25')):
    th = 1j*pi*tt/T
    F1 = c(tt, [a])*exp(logf(tt, [a])); F2 = c(tt, [a+1])*exp(logf(tt, [a+1]))
    vals.append(F2/F1/(blk(a+1, th)*blk(T - a - 1, th)))
print('F_{a+1}/F_a / [(a+1)(T-a-1)] at three t:', [mp.nstr(v, 8) for v in vals], '(constant => CDD factor confirmed up to normalisation)')
