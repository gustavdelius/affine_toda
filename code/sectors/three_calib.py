import numpy as np
from lib import *
# calibration 1: sine-Gordon (sl2) and a_2^(1) vector solitons (sl3): TW: 2-body phases, 3-body sl3 not
for n in (2, 3):
    for mu in (0.4, 1.1, 2.3):
        q = np.exp(1j*mu)
        two = maxdev_scan(lambda x: jimbo(n, q, x), np.linspace(-8, 8, 161))
        three, kr = scan3(lambda x: jimbo(n, q, x), n, L=4, n=300)
        print(f"sl_{n} vector, q=e^(i {mu}): 2-body maxdev {two[0]:.1e};  3-body maxdev {three[0]:.2e} at logx {np.round(three[1],2)}  (Krein pairing err {kr:.0e})")
