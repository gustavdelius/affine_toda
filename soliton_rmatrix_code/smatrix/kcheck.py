import numpy as np
from allS import *
import allS as AS
AS.S0 = lambda theta: S.R(np.exp(mu*theta))
def br(x, l, s): return (x - s*q**l)/(1 - s*x*q**l)
# QSIGN=1: q = e^{-i pi omega}; QSIGN=-1: the physical q = -e^{-i pi omega}. Odd powers of q flip the subscripts:
# k_32 = -<1>_+ becomes +<1>_- (the overall sign of S_32 changes), k_12 = <1>_i <3>_-i becomes <1>_-i <3>_i.
s_ = 1 if float(__import__("os").environ.get("QSIGN", "1")) > 0 else -1
forms = {(3, 3): lambda x: 1.0,
         (3, 2): lambda x: br(x, 1, s_),
         (2, 2): lambda x: br(x, 2, 1),
         (3, 1): lambda x: br(x, 0, 1j)*br(x, 2, -1j),
         (1, 2): lambda x: br(x, 1, s_*1j)*br(x, 3, -s_*1j),
         (1, 1): lambda x: br(x, 2, 1)*br(x, 4, -1)/br(x, 2, -1)}
for (a, b), fm in forms.items():
    vals = []
    for th in (0.13+0.07j, -0.21+0.33j, 0.4-0.1j, 0.05+0.9j):
        x = np.exp(mu*th); vals.append(K_top(a, b, th)[0]/fm(x))
    vals = np.array(vals); C = vals.mean()
    print(f"k_{a}{b}(x) = C * form :  C = {C.real:+.6f}{C.imag:+.6f}i   (spread {np.abs(vals - C).max():.1e})")
# reverse order amplitudes: k_ba(x) should be 1/k_ab(1/x)
for a, b in [(1, 3), (2, 3), (2, 1)]:
    th = 0.17+0.21j; x = np.exp(mu*th)
    kab = K_top(b, a, -th)[0]; kba = K_top(a, b, th)[0]
    print(f"k_{a}{b}(x) * k_{b}{a}(1/x) = {kba*kab:.6f}")
