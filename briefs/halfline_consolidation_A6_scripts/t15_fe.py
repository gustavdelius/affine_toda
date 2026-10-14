"""Check the functional equations of A^(sig) for d_n, and evaluate the breather vs candidate block form."""
import mpmath as mp
from dn_scal import DN, blk
for n, om, mu in ((4, '0.37', ('0.13', '0.21')), (4, '0.61', ('-0.17', '0.11')), (5, '0.53', ('0.21', '-0.07'))):
    D = DN(n, mp.mpf(om)); mu = mp.mpc(*mu); T = D.T
    for sig in (1, -1):
        A = lambda m: D.A(m, sig)
        uni = A(mu)*A(-mu)*D.Uhat(mu, sig)
        cr = A(T/2 - mu)/A(T/2 + mu)*D.Psi(mu, sig)/D.F(2*mu)
        print(f'n={n} om={om} sig={sig}: unitarity {mp.nstr(uni, 12)}, crossing ratio {mp.nstr(cr, 12)}', flush=True)
