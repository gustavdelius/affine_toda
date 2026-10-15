import numpy as np
# Truncation scaling of the breather-soliton closed forms: run with JTRUNC=4000 and JTRUNC=16000 and compare.
# The sign of S[B3,3] depends on the sign of q (README, 'Sign of q').
src = open('breathers.py').read()
exec(src.split("results = {}")[0].split("# superset of real pole positions")[0])
exec(open('verify_w.py').read().split("def model(sign, bl, th):")[1].join(["def model(sign, bl, th):", ""]) if False else "")
from fractions import Fraction as F
pred = {(2, 2): (+1, [(F(1,2),F(1,4)), (F(1,2),F(3,4)), (F(3,2),F(1,4)), (F(3,2),F(3,4)), (F(5,2),F(1,4)), (F(5,2),F(3,4)), (F(-1,2),F(1,4)), (F(-5,2),F(-5,4))]),
        (3, 3): (1 if float(__import__('os').environ.get('QSIGN', '1')) > 0 else -1, [(F(1,2),F(1,2)), (F(3,2),0), (F(3,2),1), (F(5,2),F(1,2)), (F(-1,2),0), (F(-5,2),-1)])}
def model(sign, bl, th):
    out = sign
    for A_, B_ in bl:
        y = float(A_)*omega + float(B_); out *= np.sinh(th/2 + 1j*np.pi*y/(2*T))/np.sinh(th/2 - 1j*np.pi*y/(2*T))
    return out
for (a, b), (sg, bl) in pred.items():
    r = np.array([K_top(f'B{a}', b, th)[0]/model(sg, bl, th) for th in (0.37+0.21j, -0.8+1.3j, 1.1-0.4j)])
    print(f"J={J:6d}  S[B{a},{b}]: max |ratio-1| = {np.abs(r - 1).max():.2e}")
