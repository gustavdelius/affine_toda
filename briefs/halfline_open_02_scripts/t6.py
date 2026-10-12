from dn_particle import *
from dn_scalar import KB_blocks
import numpy as np

def K1_eps0(n):
    h = 2*n - 2
    return Prod(h, {(h, -Fr(1, 2)): 1, (-1, Fr(1, 2)): 1, (1 - h, 0): 1, (Fr(h, 2) - 1, 0): 1,
                    (Fr(h, 2) + 1, -Fr(1, 2)): 1, (Fr(h, 2), 0): -1, (Fr(h, 2), -Fr(1, 2)): -1})

# 1. continuation check: KB(mu) at omega = -2/B equals K1(theta)  (mu = -i T theta/pi)
for n in (4, 5, 6):
    h = 2*n - 2
    for B in (0.317, 1.43):
        om = -2/B; T = (n - 1)*om
        th = 0.37 + 0.11j
        mu = -1j*T*th/np.pi
        print(n, B, complex(K1_eps0(n)(mp.mpc(th), B)), KB_blocks(mu, n, om))
