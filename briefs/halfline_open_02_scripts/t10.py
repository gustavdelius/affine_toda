from dn_particle import *
from caf_all import Kd
for n in (5, 7):
    s, sp = n - 1, n
    for th in THS[:2]:
        for B in BS:
            K = Kd(s, n)
            print(n, 'S_ss:', mp.nstr(abs(K(th, B)*K(th + I*PI, B)/S_dn(s, s, n)(2*th, B) - 1), 3),
                  ' S_ss\':', mp.nstr(abs(K(th, B)*K(th + I*PI, B)/S_dn(s, sp, n)(2*th, B) - 1), 3))
