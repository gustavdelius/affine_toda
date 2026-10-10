"""Check: graded transfer matrix of kink 1 = bosonic one times twist exp(i pi (Q-1) m_1), sector by sector."""
from ik_common import *
rng = np.random.default_rng(1)
for xi_pi in (0.55, 0.9, 1.3):
    q = q_of_xi(xi_pi); Rp = setup(q); Rg = lambda x: Qg @ Rp(x)
    for Nn in (2, 3, 4):
        Qs = charges(Nn); m1 = m_first(Nn); worst = 0
        for trial in range(5):
            ls = rng.uniform(-4, 4, Nn - 1)
            Tb = transferN(Rp, 3, ls); Tg = transfer_graded(Rg, ls)
            Om = np.diag(np.exp(1j*np.pi*(Qs - 1)*m1))
            worst = max(worst, np.abs(Tg - Om @ Tb).max())
        print(f"xi={xi_pi}pi N={Nn}: max |T_graded - Omega T_bos| = {worst:.1e}")
