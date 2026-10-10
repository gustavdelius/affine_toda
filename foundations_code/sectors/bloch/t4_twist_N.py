"""N kinks: max ||s|-1| of exp(i a m_1) T_1 by total charge Q >= 0 and twist alpha, random real rapidities; plus Krein error."""
from ik_common import *
import sys
Nn = int(sys.argv[1]); ntr = int(sys.argv[2]); alphas = [float(a) for a in sys.argv[3].split(',')]
xis = [float(a) for a in sys.argv[4].split(',')]
Qs = charges(Nn); m1 = m_first(Nn)
rng = np.random.default_rng(7)
lss = [rng.uniform(-6, 6, Nn - 1) for _ in range(ntr)]
print(f"N={Nn}, {ntr} random rapidity sets, columns alpha/pi = {alphas}")
for xi_pi in xis:
    xi_pi += 0.0013
    q = q_of_xi(xi_pi); Rp = setup(q); psi = psi_of_xi(xi_pi)
    Ts = [transferN(Rp, 3, ls) for ls in lss]
    print(f"xi={xi_pi:.4f}pi psi/pi={psi:+.4f} (|psi|N/pi={abs(psi)*Nn:.2f})")
    for Q in range(0, Nn + 1):
        idx = np.where(Qs == Q)[0]; cells = []; kmax = 0
        for a in alphas:
            Om = np.exp(1j*np.pi*a*m1[idx]); d = 0
            for T in Ts:
                e = np.linalg.eigvals(Om[:, None]*T[np.ix_(idx, idx)]); d = max(d, np.abs(np.abs(e) - 1).max()); kmax = max(kmax, krein_err(e))
            cells.append(d)
        print(f"   Q={Q} (dim {len(idx):2d}): " + " ".join(f"{d:7.1e}" for d in cells) + f"   Krein<= {kmax:.0e}", flush=True)
