"""Two kinks, Q=0 sector: critical twist alpha_c(psi) above which exp(i a m_1) T_1 has unimodular spectrum for all real theta."""
from ik_common import *
lx = np.concatenate([np.linspace(1e-3, 6, 600), np.linspace(6, 40, 120)])
idx = np.where(charges(2) == 0)[0]; m1 = m_first(2)[idx]
def dev(Rp, a):
    Om = np.exp(1j*np.pi*a*m1)
    return max(np.abs(np.abs(np.linalg.eigvals(Om[:, None]*Rp(s*np.exp(l))[np.ix_(idx, idx)])) - 1).max() for l in lx for s in (1,))
print(" xi/pi   psi/pi   alpha_c/pi   guess 2|psi|-1   guess")
for xi_pi in (0.42, 0.45, 0.5, 0.55, 0.6, 0.64, 0.70, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3):
    xi_pi += 0.0013
    q = q_of_xi(xi_pi); Rp = setup(q); psi = psi_of_xi(xi_pi)
    lo, hi = 0.0, 1.0
    if dev(Rp, 0.0) < 1e-8: print(f"{xi_pi:6.3f} {psi:+.4f}  unbroken at alpha=0"); continue
    for _ in range(30):
        mid = (lo + hi)/2
        if dev(Rp, mid) > 1e-8: lo = mid
        else: hi = mid
    print(f"{xi_pi:6.3f} {psi:+.4f}   {hi:.5f}      {2*abs(psi)-1:.5f}", flush=True)
