"""Relations of the continued s=+-1 amplitudes to the uniform d_n amplitude K_1 (eq-dn-uniform-amplitude):
varsigma=+1 with the Kim-fixed CDD factor = K_1 at B -> 2-B; varsigma=-1 = K_1/K_1^cl."""
import mpmath as mp
from t17_cont import K1, bb
mp.mp.dps = 30
def Ku(t, B, n):
    h = 2*n - 2; b = lambda x: bb(x, t, h)
    return b(h - B/2)*b(h/2 - 1)*b(h/2 + 1 - B/2)/(b(h - 1)*b(1 - B/2)*b(h/2)*b(h/2 - B/2))
def Kcl(t, n):
    h = 2*n - 2; b = lambda x: bb(x, t, h)
    return -b(h/2 - 1)*b(h/2 + 1)/(b(1)*b(h - 1)*b(h/2)**2)
def KN(t, B, n):
    h = 2*n - 2; b = lambda x: bb(x, t, h)
    return b(h/2 - 1)*b(h/2 + B/2)*b(h - 1 + B/2)/(b(h - 1)*b(h/2)*b(h/2 - 1 + B/2)*b(B/2))
d = 0
for n in (4, 5, 6, 7):
    for B in (mp.mpf('0.37'), mp.mpf('1.3')):
        for t in (mp.mpc('0.4', '0.1'), mp.mpc('-1.1', '0.3')):
            d = max(d, abs(K1(1, t, B, n, (0, 1, n - 2))/Ku(t, 2 - B, n) - 1), abs(K1(1, t, B, n, (0, 1, n - 2))/KN(t, B, n) - 1),
                    abs(K1(-1, t, B, n)/(Ku(t, B, n)/Kcl(t, n)) - 1))
print('max deviation', d)
