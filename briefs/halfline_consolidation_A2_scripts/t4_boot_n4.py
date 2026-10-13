"""n = 4 SP soliton bootstrap (reduced set of points; book intertwiners use the Gram path here)."""
import numpy as np
import mpmath as mp
from t4_bootstrap import boot_matrix
from sp_mp import SP, F11
mp.mp.dps = 30
n, om = 4, 0.37
T = (n + 1)*om/2
for l in (0, 2, 5):
    eta = 0.41 + 0.23j
    S = SP(n, mp.mpf(om), l, mp.mpc(eta.real, eta.imag))
    for mu in (0.11 + 0.19j, -0.23 + 0.07j):
        th = 1j*np.pi*mu/T
        P, res, same, mus = boot_matrix(n, om, l, eta, th)
        mmu = mp.mpc(mu.real, mu.imag)
        mmus = [mmu - S.T + k*S.om for k in range(1, n + 1)]
        rhs = mp.mpc(P.real, P.imag)*mp.fprod([S.a(x) for x in mmus])
        for j in range(n):
            for k in range(j + 1, n):
                rhs *= F11(mmus[j] + mmus[k], n, S.om)
        print(f'n=4 l={l} mu={mu}: matrix res {res:.1e} iota_out=iota_in {same:.1e} Q = {mp.nstr(S.bt(mmu)/rhs, 12)}',
              flush=True)
