"""Closed form of the CDD-dressed varsigma=+1 breather amplitude and its continuation."""
import mpmath as mp
from t17_cont import KB_imag, bb
from dn_scal import blk
mp.mp.dps = 30
H = mp.mpf(1)/2
d = 0
for n, om, mu in ((4, '0.37', ('0.13', '0.21')), (5, '0.53', ('0.21', '-0.07')), (6, '0.29', ('0.08', '0.15'))):
    om, mu = mp.mpf(om), mp.mpc(*mu); T = (n - 1)*om
    b = lambda x: blk(x, mu, T)
    kbn = b(om/2 - T)*b(T/2 - om/2)*b(T/2 - H)*b(H)*b(T - om/2 - H)/(b(T/2)*b(T/2 - om/2 - H))
    d = max(d, abs(KB_imag(1, mu, T, om, [0, om/2, (T - om)/2])/kbn - 1))
    # continuation of kbn
    for Bv in (mp.mpf('0.37'), mp.mpf('1.3')):
        t = mp.mpc('0.4', '0.1'); h = 2*n - 2
        o2 = -2/Bv; T2 = (n - 1)*o2; m2 = -1j*T2*t/mp.pi
        c = lambda x: blk(x, m2, T2)
        kb2 = c(o2/2 - T2)*c(T2/2 - o2/2)*c(T2/2 - H)*c(H)*c(T2 - o2/2 - H)/(c(T2/2)*c(T2/2 - o2/2 - H))
        bq = lambda x: bb(x, t, h)
        KN = bq(h/2 - 1)*bq(h/2 + Bv/2)*bq(h - 1 + Bv/2)/(bq(h - 1)*bq(h/2)*bq(h/2 - 1 + Bv/2)*bq(Bv/2))
        d = max(d, abs(kb2/KN - 1))
print('max dev', d)
