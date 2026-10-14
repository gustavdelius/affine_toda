"""Two kinks of a_2^(2), sectors Q=+1 and Q=-1: unimodular for every twist alpha in [-pi, pi]."""
from ik_common import *
Qs = charges(2); m1 = m_first(2); lx = np.linspace(-12, 12, 481); worst = 0
for xi_pi in (0.42, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.1, 1.3, 2.0, 3.0):
    Rp = setup(q_of_xi(xi_pi + 0.0013))
    for Q in (-1, 1):
        idx = np.where(Qs == Q)[0]
        for a in np.linspace(-1, 1, 41):
            Om = np.exp(1j*np.pi*a*m1[idx])
            worst = max(worst, max(np.abs(np.abs(np.linalg.eigvals(Om[:, None]*Rp(np.exp(l))[np.ix_(idx, idx)])) - 1).max() for l in lx))
print('two kinks, Q=+-1, alpha/pi in [-1,1], 12 couplings, log x in [-12,12]: max ||s|-1| =', f'{worst:.1e}')
