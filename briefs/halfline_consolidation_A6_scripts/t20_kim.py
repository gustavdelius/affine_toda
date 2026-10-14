"""Kim's self-dual Neumann conjecture for the d_n vector particle (hep-th/9506031, fig. 7 and the equation for K_{n-2}):
K = prod_{k=1,3,..,2h-3} [k/2] * [(h-1)/2], [x] = (x-1/2)(x+1/2)/((x-1/2+B/2)(x+1/2-B/2)); for d_4 it is Kim-Cho's
(hep-th/9505138) exact conjecture [1/2][3/2][5/2]^2[7/2][9/2]. Compared at O(beta^2) (B = beta^2/2pi) with the s = +-1 amplitudes;
and the CAF identity Kim = R~_1(theta + i pi) (hep-th/0308004 (KimCF), (dn1))."""
import mpmath as mp
from t17_cont import K1, bb

mp.mp.dps = 40


def kim(t, B, n):
    h = 2*n - 2
    sq = lambda x: bb(x - mp.mpf(1)/2, t, h)*bb(x + mp.mpf(1)/2, t, h)/(bb(x - mp.mpf(1)/2 + B/2, t, h)*bb(x + mp.mpf(1)/2 - B/2, t, h))
    r = mp.fprod(sq(mp.mpf(k)/2) for k in range(1, 2*h - 2, 2))
    return r*sq(mp.mpf(h - 1)/2)


def caf_tilde(t, B, n):
    """CAF (dn1) for j = 1: R~_1(theta + i pi) = ||h-1|| prod_{k=1}^{h-1} ||2k-1||, with their block (blockw)."""
    h = 2*n - 2
    a1, a2 = 1j*mp.pi*(2 - B)/(2*h), 1j*mp.pi*B/(2*h)
    s = lambda z: mp.sinh(z/2)
    ang = lambda x, th: (s(th + (x - 1)/2*a1 + (x - 1)/2*a2)*s(th + (x + 1)/2*a1 + (x + 1)/2*a2)
                         / (s(th + (x - 1)/2*a1 + (x + 1)/2*a2)*s(th + (x + 1)/2*a1 + (x - 1)/2*a2)))
    W = lambda x, th: ang(x, th)/ang(x, -th)
    return W(h - 1, t)*mp.fprod(W(2*k - 1, t) for k in range(1, h))


e = mp.mpf('1e-15')
d1 = lambda f, t: (f(e, t) - f(-e, t))/(2*e)/(2*mp.pi)
ths = [mp.mpf('0.37'), mp.mpc('1.2', '0.3'), mp.mpc('-0.8', '0.1')]
for n in (4, 5, 6, 7):
    rp = max(abs(d1(lambda B, t: K1(1, t, B, n, (0, 1, n - 2)), t) - d1(lambda B, t: kim(t, B, n), t)) for t in ths)
    rp0 = max(abs(d1(lambda B, t: K1(1, t, B, n), t) - d1(lambda B, t: kim(t, B, n), t)) for t in ths)
    rm = max(abs(d1(lambda B, t: K1(-1, t, B, n), t) - d1(lambda B, t: kim(t + 1j*mp.pi, B, n), t)) for t in ths)
    rm0 = max(abs(d1(lambda B, t: K1(-1, t, B, n), t) - d1(lambda B, t: kim(t, B, n), t)) for t in ths)
    # beyond one loop: compare exact values
    Bv, t = mp.mpf('0.6'), mp.mpc('0.5', '0.2')
    ex_p = abs(K1(1, t, Bv, n, (0, 1, n - 2))/kim(t, Bv, n) - 1)
    ex_m = abs(K1(-1, t, Bv, n)/kim(t + 1j*mp.pi, Bv, n) - 1)
    cafk = abs(kim(t, Bv, n)/caf_tilde(t, Bv, n) - 1)   # caf_tilde(t) = R~_1(t + i pi)
    print(f'd_{n}: one loop |sig=+1 with CDD (0,1,{n-2}) - Kim| = {mp.nstr(rp, 3)} (without CDD {mp.nstr(rp0, 3)}); '
          f'|sig=-1 - Kim(theta+i pi)| = {mp.nstr(rm, 3)} (vs Kim(theta): {mp.nstr(rm0, 3)}); exact at B=0.6: '
          f'sig=+1/Kim - 1 = {mp.nstr(ex_p, 3)}, sig=-1/Kim(th+i pi) - 1 = {mp.nstr(ex_m, 3)}; Kim/CAF R~_1(th+i pi) - 1 = '
          f'{mp.nstr(cafk, 3)}', flush=True)
