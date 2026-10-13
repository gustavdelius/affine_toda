"""Consistency checks of the double-row transfer matrices: commutativity (decides the left-wall convention),
Krein-unitarity T^dag C T = C, Z_N symmetry."""
import numpy as np
from b_transfer import Tdr, Zshift, kron_all
from p07_common import Cmat
n = 2; d = 3
for om in (0.37, 1.43):
    for cfg in ('UU', 'NN', 'UN', 'NU', 'MM', 'NM'):
        for mirror in ('omega', 'swap'):
            th = [0.9, 0.35]; types = ['V', 'V']
            T1 = Tdr(0, types, th, cfg[0], cfg[1], n, om, mirror); T2 = Tdr(1, types, th, cfg[0], cfg[1], n, om, mirror)
            comm = np.abs(T1 @ T2 - T2 @ T1).max()/(np.abs(T1).max()*np.abs(T2).max())
            th3 = [1.2, 0.5, 0.2]; t3 = ['V', 'V', 'V']
            Ts = [Tdr(j, t3, th3, cfg[0], cfg[1], n, om, mirror) for j in range(3)]
            comm3 = max(np.abs(Ts[a] @ Ts[b] - Ts[b] @ Ts[a]).max() for a in range(3) for b in range(a))
            CC = kron_all([Cmat(d)]*3)
            kr = max(np.abs(T.conj().T @ CC @ T - CC).max() for T in Ts)
            s = -1 if 'M' in cfg else 1
            ok_z = all(cfg.count(c) in (0, 2) for c in 'NM') or cfg in ('UN', 'NU')
            Z3 = kron_all([Zshift(n, s if cfg in ('MM',) else 1)]*3)
            zc = max(np.abs(Z3 @ T - T @ Z3).max() for T in Ts)
            print(f'om={om} {cfg} mirror={mirror:5s}: [T1,T2] {comm:.1e}  N=3 max comm {comm3:.1e}  '
                  f'Krein {kr:.1e}  [Z,T] {zc:.1e}')
