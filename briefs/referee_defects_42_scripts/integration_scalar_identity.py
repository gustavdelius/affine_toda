"""Integration check (ch. 42, item 3): the scalar identity behind the fused-scalar check of prp-defect-fusion-k (3),
  sigma(i lambda theta_d, -i theta) = (1 + yd^-2)/(1 + y^-2 yd^-2) * Phi(theta - theta_d)/Phi(-theta - theta_d),
and the matrix step k_+ + yd^-2 k_- = 4 yd^-1 P_0^+ (eq-bsg-k-full, u = -i theta, vartheta = lambda theta_d).
sigma(x, u): BPTT form as in eq-bsg-k-full (Gamma product, written here from ch. 29's definition); Phi: eq-sgd-scalar
via defect_common.sg_transmission_scalar (Richardson-summed, ~1e-9 at complex rapidity)."""
import sys
import numpy as np
import mpmath as mp
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'textbook_code'))
from defect_common import sg_transmission_scalar

LAM = 1.37
mp.mp.dps = 25
LL = mp.mpf(LAM)


def log_sigma(x, u):
    x, u = mp.mpmathify(x), mp.mpmathify(u)

    def t(l, uu):
        v, a = LL*uu/mp.pi, x/mp.pi
        return (mp.loggamma(0.5 + a + (2*l - 1)*LL - v) + mp.loggamma(0.5 - a + (2*l - 1)*LL - v)
                - mp.loggamma(0.5 - a + (2*l - 2)*LL - v) - mp.loggamma(0.5 + a + 2*l*LL - v))
    return mp.log(mp.cos(x)/mp.cos(x + LL*u)) + mp.nsum(lambda l: t(l, u) - t(l, -u), [1, mp.inf])


worst_m, worst_s = 0, 0
for th, thd, eta in ((0.33, 0.21, 0.37), (0.21 + 0.17j, 0.21, 0.37), (-0.4 + 0.05j, -0.3, 0.8)):
    y, yd, r = np.exp(LAM*th), np.exp(LAM*thd), np.exp(1j*eta)
    u = -1j*th
    vt = LAM*thd
    kp, km = y*r + 1/(y*r), y/r + r/y
    P0p = np.cos(LAM*u)*np.cos(eta)*np.cosh(vt) - np.sin(LAM*u)*np.sin(eta)*np.sinh(vt)
    P0m = np.cos(LAM*u)*np.cos(eta)*np.cosh(vt) + np.sin(LAM*u)*np.sin(eta)*np.sinh(vt)
    worst_m = max(worst_m, abs(kp + km/yd**2 - 4*P0p/yd), abs(km + kp/yd**2 - 4*P0m/yd))
    lhs = complex(mp.exp(log_sigma(1j*vt, u)))
    rhs = (1 + yd**-2)/(1 + y**-2*yd**-2)*sg_transmission_scalar(th - thd, LAM)/sg_transmission_scalar(-th - thd, LAM)
    worst_s = max(worst_s, abs(lhs/rhs - 1))
print(f'matrix step k_+ + yd^-2 k_- = 4 P_0^+/yd (and the - entry): {worst_m:.1e}')
print(f'scalar identity sigma(i lambda theta_d, -i theta) = ... Phi/Phi: max rel. dev. {worst_s:.1e} (Richardson Phi ~1e-9)')
print('PASS' if worst_m < 1e-12 and worst_s < 1e-7 else 'FAIL')
