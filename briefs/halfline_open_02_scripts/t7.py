from dn_boot import *
from t6 import K1_eps0
for n in (4, 5, 6, 7):
    res, boot, clo = check_vec(K1_eps0(n), n)
    print(n, 'unit/cross per a:', {a: (f'{u:.0e}', f'{c:.0e}') for a, (u, c) in res.items()}, f'boot {boot:.0e}', f'closure 1+2->1 {clo:.1e}')
