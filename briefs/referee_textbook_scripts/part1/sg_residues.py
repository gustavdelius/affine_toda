"""Residues of S_T, S_R and of the channel amplitudes S_T +- S_R at the breather poles theta_n = i(pi - n xi)."""
from mpmath import mp, mpf, mpc, pi, sinh, sin, tan, cos
from sg_s0 import S0
mp.dps = 15
for lam in (mpf('1.5'), mpf('3.5'), mpf('4.3')):
    xi = pi/lam
    n = 1
    while n < lam:
        th = 1j*(pi - n*xi)
        s0 = S0(th, lam)
        rT = sinh(lam*th)*s0/lam          # 1/sinh(lam(i pi - th)) has residue 1/lam (times (-1)^n... computed below)
        # exact residue of 1/sinh(lam(i pi - th)) at th_n: d/dth sinh(lam(i pi - th)) = -lam cosh(i pi n) = -lam (-1)^n
        r = 1/(-lam*(-1)**n)
        RT = sinh(lam*th)*s0*r; RR = sinh(1j*pi*lam)*s0*r
        print(f'lam={lam}, n={n}: S0(theta_n)={mp.nstr(s0,8)}, Res S_T={mp.nstr(RT,8)}, Res S_R={mp.nstr(RR,8)}, '
              f'Res(S_T+S_R)={mp.nstr(RT+RR,6)}, Res(S_T-S_R)={mp.nstr(RT-RR,6)}')
        n += 1
    # B1 B1 -> B2 pole of S_11 = f_{xi/pi} at i xi: residue 2 i tan xi
    print(f'   S_11 residue at i xi: 2i tan xi = {mp.nstr(2j*tan(xi),6)}  (B2 exists: {2 < lam})')
    # s B1 -> s pole of S_{sB1} = f_{1/2 - xi/2pi} at i(pi+xi)/2: residue -2i cot(xi/2)
    print(f'   S_sB1 residue at i(pi+xi)/2 (s-channel s+B1->s): {mp.nstr(-2j/tan(xi/2),6)}')
