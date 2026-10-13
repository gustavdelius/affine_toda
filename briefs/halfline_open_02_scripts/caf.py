"""Castro-Alvaredo--Fring hep-th/0308004 d_n solutions (dn1)-(dn2) in book blocks, and ambiguity tests."""
from dn_particle import *
from t6 import K1_eps0


def wblk(x, th, B, h):
    """CAF ||x||_theta = ((x-1)/2)((x+1)/2)/(((x-1+B)/2)((x+1-B)/2)) in book blocks (x)_h."""
    return (blk(mp.mpf(x - 1)/2, th, h)*blk(mp.mpf(x + 1)/2, th, h)
            / (blk((x - 1 + B)/2, th, h)*blk((x + 1 - B)/2, th, h)))


def Rt(j, n, th, B):
    """CAF R~_j(theta) (their (dn1), (dn2) give R~_j(theta + i pi))."""
    h = 2*n - 2
    t = th - I*PI
    if j <= n - 2:
        v = mp.mpf(1)
        for k in range(1, j + 1):
            v *= wblk(h - 2*k + 1, t, B, h)
        for p in range(1, j + 1):
            for k in range(p, h - p + 1):
                v *= wblk(2*k - 1, t, B, h)
        return v
    v = mp.mpf(1)
    for p in range(1, n//2 + 1):
        for k in range(2*p - 1, h - 2*p + 2):
            v *= wblk(2*k - 1, t, B, h)
    return v


def Sprod(i, n, th, Bp):
    v = mp.mpf(1)
    for j in range(1, n + 1):
        v *= S_dn(i, j, n)(th, Bp)
    return v


if __name__ == '__main__':
    for n in (4, 5, 6):
        h = 2*n - 2
        K1 = K1_eps0(n)
        th = mp.mpc('0.37', '0.11')
        for B in (mp.mpf('0.317'), mp.mpf('1.43')):
            r = K1(th, B)
            print(n, float(B))
            # self-duality and the CAF crossing check
            print('   CAF R1 unitarity', mp.nstr(abs(Rt(1, n, th, B)*Rt(1, n, -th, B) - 1), 3),
                  ' cross', mp.nstr(abs(Rt(1, n, th, B)*Rt(1, n, th + I*PI, B)/S_dn(1, 1, n)(2*th, B) - 1), 3))
            for name, Bp in (('1-B/2', 1 - B/2), ('B/2', B/2), ('B', B), ('2-B', 2 - B), ('1+B/2', 1 + B/2)):
                for shift in (0, 1):
                    R = Rt(1, n, th + shift*I*PI, B)
                    for e in (1, -1):
                        val = R*Sprod(1, n, th, Bp)**e
                        if abs(val/r - 1) < 1e-12 or abs(val/r + 1) < 1e-12:
                            print('   MATCH: K1 = ', '+-'[abs(val/r - 1) > 1e-12], f'R~(th+{shift} i pi) prod S_1j(th,{name})^{e}')
