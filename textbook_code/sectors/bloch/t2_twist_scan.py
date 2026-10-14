"""Twisted bosonic transfer matrices exp(i a m_1) T_1 of N IK kinks, by total charge Q: max ||s|-1| and Krein pairing."""
from ik_common import *
import sys
Nn = int(sys.argv[1]); ntr = int(sys.argv[2]) if len(sys.argv) > 2 else 60
Qs = charges(Nn); m1 = m_first(Nn)
alphas = np.linspace(0, 1, 9)            # alpha/pi
rng = np.random.default_rng(3)
lss = [rng.uniform(-5, 5, Nn - 1) for _ in range(ntr)]
for xi_pi in (0.30, 0.45, 0.55, 0.62, 0.70, 0.9, 1.1, 1.3, 2.0, 3.0):
    xi_pi += 0.0013                       # avoid roots of unity
    q = q_of_xi(xi_pi); Rp = setup(q); psi = psi_of_xi(xi_pi)
    Ts = [transferN(Rp, 3, ls) for ls in lss]
    out = []
    for Q in range(0, Nn + 1):
        idx = np.where(Qs == Q)[0]
        row = []
        for a in alphas:
            Om = np.exp(1j*np.pi*a*m1[idx]); d = 0; kr = 0
            for T in Ts:
                e = np.linalg.eigvals(Om[:, None]*T[np.ix_(idx, idx)]); d = max(d, np.abs(np.abs(e) - 1).max()); kr = max(kr, krein_err(e))
            row.append((d, kr))
        out.append((Q, row))
    print(f"xi={xi_pi:.2f}pi psi/pi={psi:+.3f}  N={Nn}  (alpha/pi = {', '.join(f'{a:.3g}' for a in alphas)})")
    for Q, row in out:
        print(f"   Q={Q}: dev " + " ".join(f"{d:7.1e}" for d, _ in row) + "   Krein " + " ".join(f"{k:5.0e}" for _, k in row))
    sys.stdout.flush()
