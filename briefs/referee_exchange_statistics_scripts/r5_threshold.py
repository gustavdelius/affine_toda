# Two a_2^(2) kinks, Q=0: critical twist alpha_c by bisection (threshold 1e-6), finer rapidity grid incl. |log x|<1e-2, large log x, both signs of log x
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'sectors', 'bloch'))
from ik_common import *
lpos = np.concatenate([np.logspace(-5, -2, 40), np.linspace(1e-2, 6, 900), np.linspace(6, 80, 150)])
lx = np.concatenate([lpos, -lpos])
idx = np.where(charges(2) == 0)[0]; m1 = m_first(2)[idx]
def dev(Rs, a):
    Om = np.exp(1j*np.pi*a*m1)
    return max(np.abs(np.abs(np.linalg.eigvals(Om[:, None]*B)) - 1).max() for B in Rs)
def formula(p): return max(0, min(2*abs(p) - 1, 1 - abs(p)))
for psi_t in [float(v) for v in sys.argv[1].split(',')]:
    k = 0 if psi_t > 0 else 1
    xi_pi = 2/(3*(psi_t + 2*k)) + 0.0013
    psi = psi_of_xi(xi_pi); Rp = setup(q_of_xi(xi_pi)); Rs = [Rp(np.exp(l))[np.ix_(idx, idx)] for l in lx]
    d0 = dev(Rs, 0.0)
    if d0 < 1e-6: ac = 0.0
    else:
        lo, hi = 0.0, 1.0
        for _ in range(22):
            mid = (lo + hi)/2
            if dev(Rs, mid) > 1e-6: lo = mid
            else: hi = mid
        ac = hi
    # also check broken just below and unbroken just above the formula value
    f = formula(psi)
    chk = (dev(Rs, max(f - 0.01, 0)) > 1e-6 if f > 0.01 else None, dev(Rs, f + 0.01) < 1e-6)
    print(f"xi/pi={xi_pi:.4f} psi/pi={psi:+.4f}: dev(alpha=0)={d0:.1e}, alpha_c/pi={ac:.5f}, formula {f:.5f}, diff {abs(ac-f):.1e}; broken at f-0.01: {chk[0]}, unbroken at f+0.01: {chk[1]}", flush=True)
