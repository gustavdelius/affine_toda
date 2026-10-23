"""Search for a diagonal gauge x^{a m} in which the Izergin-Korepin R-matrix satisfies R^T = P R P (none found)."""
from ik_common import *
q = q_of_xi(0.9013); Rp = setup(q); P = flip(3, 3); m = np.array([1., 0, -1])
for a in (-1, -0.5, -1/3, 0.5, 1/3, 1):
    for b in (0, 1):
        err = 0
        for x in (0.37, 1.9, 7.3):
            Dg = np.kron(np.diag(x**(a*m)), np.diag(x**(-b*a*m)))   # gradation change x^{a m} (x) x^{-b a m}
            R2 = Dg @ Rp(x) @ np.linalg.inv(Dg)
            err = max(err, np.abs(R2.T - P @ R2 @ P).max())
        print(f"a={a:+.3f} b={b}: |R'^T - P R' P| = {err:.1e}")
