"""Soliton-excited-soliton of a_2^(2) (R-part Rcheck(-z)): breaking by total charge Q at alpha = 0 and alpha = pi, two particles."""
from ik_common import *
lx = np.concatenate([np.linspace(0.02, 8, 300), np.linspace(8, 30, 60)])
Qs = charges(2); m1 = m_first(2)
for xi_pi in (0.3, 0.45, 0.6, 0.9, 1.5, 3.0):
    xi_pi += 0.0013; q = q_of_xi(xi_pi); Rp = setup(q)
    Rs = [Rp(-np.exp(l)) for l in lx]; row = []
    for Q in (0, 1, 2):
        idx = np.where(Qs == Q)[0]
        for a in (0.0, 1.0):
            Om = np.exp(1j*np.pi*a*m1[idx])
            d = max(np.abs(np.abs(np.linalg.eigvals(Om[:, None]*R[np.ix_(idx, idx)])) - 1).max() for R in Rs)
            row.append(f"Q={Q},a={a:.0f}pi: {d:.1e}")
    print(f"xi={xi_pi:.4f}pi psi/pi={psi_of_xi(xi_pi):+.3f}: " + "; ".join(row), flush=True)
