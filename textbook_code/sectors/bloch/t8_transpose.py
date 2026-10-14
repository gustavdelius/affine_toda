"""Transpose symmetry R_12^T = R_21 (time reversal) and sigma(Omega_a T) = sigma(Omega_{-a} T)."""
from ik_common import *
import spin_common as sc
xs = [0.37, 1.9, 7.3]
q = q_of_xi(0.9013); Rp = setup(q); P = flip(3, 3)
print("IK: |R^T - P R P| =", max(np.abs(Rp(x).T - P @ Rp(x) @ P).max() for x in xs))
for n in (2, 3):
    D, Wn, s2 = sc.make(n, 'c'); Rs = s2(np.exp(-1j*np.pi*1.61)); P = flip(D, D)
    print(f"c_{n} spinor: |R^T - P R P| =", max(np.abs(Rs(x).T - P @ Rs(x) @ P).max() for x in xs))
# spectra at +alpha and -alpha, N=3, IK
Qs = charges(3); m1 = m_first(3); ls = [0.7, -1.9]; T = transferN(Rp, 3, ls)
for a in (0.3, 0.5):
    for Qv in (0, 1):
        idx = np.where(Qs == Qv)[0]
        ep = np.sort_complex(np.linalg.eigvals(np.exp(1j*np.pi*a*m1[idx])[:, None]*T[np.ix_(idx, idx)]))
        em = np.linalg.eigvals(np.exp(-1j*np.pi*a*m1[idx])[:, None]*T[np.ix_(idx, idx)])
        print(f"IK N=3 Q={Qv} alpha={a}pi: max dist sigma(+a) to sigma(-a) = {max(np.min(np.abs(em - z)) for z in ep):.1e}; max||s|-1| = {np.abs(np.abs(ep)-1).max():.2e}")
