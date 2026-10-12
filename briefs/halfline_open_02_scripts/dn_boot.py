"""Boundary bootstrap for the d_n vector-like particles a = 1..n-2 from K_1, in exact block arithmetic."""
from dn_particle import *


def shift(P, s):
    """P(theta + i pi s/h): (x)(th + i pi s/h) = (x + s) / ... -> returns callable"""
    return lambda th, B: P(th + I*PI*s/P.h, B)


def K_vec(K1, n):
    """K_a for a = 1..n-2 by K_{a+1}(th) = K_a(th + i pi/h) S_{a1}(2th + i pi(1-a)/h) K_1(th - i pi a/h), as callables."""
    h = 2*n - 2
    Ks = {1: K1}
    for a in range(1, n - 2):
        Ka, Sa1 = Ks[a], S_dn(a, 1, n)
        Ks[a + 1] = (lambda th, B, Ka=Ka, Sa1=Sa1, a=a:
                     Ka(th + I*PI/h, B)*Sa1(2*th + I*PI*(1 - a)/h, B)*K1(th - I*PI*a/h, B))
    return Ks


def check_vec(K1, n, verbose=False):
    """crossing-unitarity K_a(th) K_a(th + i pi) = S_aa(2th), unitarity, and the general bootstrap a + b -> a + b."""
    h = 2*n - 2
    Ks = K_vec(K1, n)
    res = {}
    for a, Ka in Ks.items():
        u = max(abs(Ka(th, B)*Ka(-th, B) - 1) for th in THS for B in BS)
        c = max(abs(Ka(th, B)*Ka(th + I*PI, B)/S_dn(a, a, n)(2*th, B) - 1) for th in THS for B in BS)
        res[a] = (float(u), float(c))
    boot = 0
    for a in Ks:
        for b in Ks:
            if a + b in Ks:
                for th in THS:
                    for B in BS:
                        lhs = Ks[a + b](th, B)
                        rhs = Ks[a](th + I*PI*b/h, B)*S_dn(a, b, n)(2*th + I*PI*(b - a)/h, B)*Ks[b](th - I*PI*a/h, B)
                        boot = max(boot, float(abs(lhs/rhs - 1)))
    # closure of the vector: 1 + 2 -> 1 at u_12^1 = pi - pi/h: th_1 = th + i(pi - 2pi/h), th_2 = th - i pi/h
    clo = 0
    if 2 in Ks:
        for th in THS:
            for B in BS:
                rhs = Ks[1](th + I*(PI - 2*PI/h), B)*S_dn(1, 2, n)(2*th + I*(PI - 3*PI/h), B)*Ks[2](th - I*PI/h, B)
                clo = max(clo, float(abs(Ks[1](th, B)/rhs - 1)))
    return res, boot, clo
