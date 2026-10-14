"""Bootstrap at omega = 1.29 (several l off the principal branch of the worker's k)."""
import numpy as np
import mpmath as mp
from t4_bootstrap import boot_matrix
from sp_mp import SP, F11
mp.mp.dps = 30
for n, ls in ((2, (0, 2, 3)), (3, (0, 3, 4))):
    om = 1.29
    T = (n + 1)*om/2
    for l in ls:
        eta = -0.3 + 1.1j
        S = SP(n, mp.mpf(om), l, mp.mpc(eta.real, eta.imag))
        mu = 0.11 + 0.19j
        P, res, same, mus = boot_matrix(n, om, l, eta, 1j*np.pi*mu/T)
        mmu = mp.mpc(mu.real, mu.imag)
        mmus = [mmu - S.T + k*S.om for k in range(1, n + 1)]
        rhs = mp.mpc(P.real, P.imag)
        for j in range(n):
            for k in range(j + 1, n):
                rhs *= F11(mmus[j] + mmus[k], n, S.om)
        Q = {m: S.bt(mmu, m)/(rhs*mp.fprod([S.a(x, m) for x in mmus])) for m in (0, 2, -2)}
        print(f'n={n} om={om} l={l}: matrix res {res:.1e}; Q(m=0) = {mp.nstr(Q[0], 12)}, Q(m=+2) = {mp.nstr(Q[2], 6)}, '
              f'Q(m=-2) = {mp.nstr(Q[-2], 6)}', flush=True)
