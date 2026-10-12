"""Calibration of the mixed-representation residue rules on boundary Lee-Yang (h = 6 blocks), K_(1) = (1)(3)/(4).
GZ: Res_{i pi/2} K = (i/2) g^2 ; near theta* = i pi/6:  K ~ -(i/2) Gamma g /(theta - theta*)  (GZ (3.41), (3.45)).
Mixed-rep rules (this work): one-point pole = (i/sqrt2) Gamma g_s / dtheta, where g_s is the coupling defined
from the s-channel alone; at theta = i pi/2 the s- and u-channel trees add with factor (1 - S(i pi)) = 2,
so g_s^2 = g^2/2 and the one-point residue becomes (i/2) Gamma g in modulus: check |Res| = |Gamma g|/2."""
import mpmath as mp
from fractions import Fraction as Fr
from p08_common import Prod, I, PI, blk
mp.mp.dps = 30
h = 6
K = Prod(h, {(1, 0): 1, (3, 0): 1, (4, 0): -1})
S = lambda th: (mp.sinh(th) + I*mp.sin(2*PI/3))/(mp.sinh(th) - I*mp.sin(2*PI/3))
o, R2 = K.laurent(Fr(3), Fr(0))       # theta = i pi 3/6 = i pi/2
o1, R6 = K.laurent(Fr(1), Fr(0))      # theta = i pi/6
eps = mp.mpf('1e-20')
RS = S(2*PI*I/3 + eps)*eps            # Res S at 2 pi i/3
Gam2 = RS/I
g2 = R2/(I/2)
print('Res_{i pi/2} K =', mp.nstr(R2, 12), '  g^2 =', mp.nstr(g2, 12), ' (GZ: g = 2i sqrt(2 sqrt3 - 3): g^2 =',
      mp.nstr(-4*(2*mp.sqrt(3) - 3), 12), ')')
print('Res_{2pi i/3} S = i Gamma^2:  Gamma^2 =', mp.nstr(Gam2, 12))
print('S(i pi) =', mp.nstr(S(I*PI*(1 - mp.mpf('1e-25'))), 6), ' -> s+u channel factor 1 - S(i pi) = 2')
print('|Res_{i pi/6} K| =', mp.nstr(abs(R6), 12), '   |Gamma g|/2 =', mp.nstr(mp.sqrt(abs(Gam2*g2))/2, 12),
      '   |Gamma g|/sqrt2 (s-channel only) =', mp.nstr(mp.sqrt(abs(Gam2*g2))/mp.sqrt(2), 12))
