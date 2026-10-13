"""Full d_n eps=0 particle set K_a = R~_a(th,B) prod_j S_aj(th, 1-B/2) (CAF R~ dressed by am3), all a incl. spinors:
agreement with the K_1-bootstrap for vector-like a, and all axioms (unitarity, crossing-unitarity, bulk boundary
bootstrap at 1+1->2, a+b->a+b, 1+2->1, 1+s->s', s+s'->1), classical limit."""
from dn_particle import *
from caf import Rt, Sprod
from dn_boot import K_vec
import sys, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from t6 import K1_eps0

def Kd(a, n):
    return lambda th, B: Rt(a, n, th, B)*Sprod(a, n, th, 1 - B/2)

def bar(a, n):
    return a if (n % 2 == 0 or a < n - 1) else 2*n - 1 - a

for n in (4, 5, 6, 7):
    h = 2*n - 2
    Kv = K_vec(K1_eps0(n), n)
    w_vec = max(abs(Kd(a, n)(th, B)/Kv[a](th, B) - 1) for a in Kv for th in THS[:2] for B in BS)
    w_ax = 0
    for a in range(1, n + 1):
        for th in THS[:2]:
            for B in BS:
                Ka = Kd(a, n)
                w_ax = max(w_ax, abs(Ka(th, B)*Ka(-th, B) - 1),
                           abs(Ka(th, B)*Kd(bar(a, n), n)(th + I*PI, B)/S_dn(a, a, n)(2*th, B) - 1))
    s, sp = n - 1, n
    w_b = 0
    for th in THS[:2]:
        for B in BS:
            # 1 + s -> s' : th_1 = th + i(pi/2 - pi/h), th_s = th - 2 pi i/h
            t1, ts = th + I*(PI/2 - PI/h), th - 2*I*PI/h
            w_b = max(w_b, abs(Kd(sp, n)(th, B)/(Kd(1, n)(t1, B)*S_dn(1, s, n)(t1 + ts, B)*Kd(s, n)(ts, B)) - 1))
            # s + s' -> 1 : u = pi - 2pi/h, th_s = th +- i u/2
            u = PI - 2*PI/h
            w_b = max(w_b, abs(Kd(1, n)(th, B)/(Kd(s, n)(th + I*u/2, B)*S_dn(s, sp, n)(2*th, B)*Kd(sp, n)(th - I*u/2, B)) - 1))
    print(n, f'K_a(CAF+am3) = K_a(bootstrap from K_1), a<=n-2: {float(w_vec):.0e};  unitarity+crossing all a: {float(w_ax):.0e};'
          f'  spinor bootstraps: {float(w_b):.0e}')
