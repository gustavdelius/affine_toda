"""Which one-site operators c give C-invariance (c x c) R (c x c) = R or R_21, and Krein relation R^dag = (c x c) R^{-1} (c x c), for IK and spinors."""
from ik_common import *
import spin_common as sc
def tests(name, Rp, D, cands, xs):
    P = flip(D, D)
    for cn, c in cands.items():
        CC = np.kron(c, c); ci = np.linalg.inv(CC)
        e1 = max(np.abs(CC @ Rp(x) @ ci - Rp(x)).max() for x in xs)
        e2 = max(np.abs(CC @ Rp(x) @ ci - P @ Rp(x) @ P).max() for x in xs)
        e3 = max(np.abs(Rp(x).conj().T - CC @ np.linalg.inv(Rp(x)) @ ci).max() for x in xs)
        e4 = max(np.abs(Rp(x).conj().T - CC @ np.linalg.inv(P @ Rp(x) @ P) @ ci).max() for x in xs)
        print(f"  {name} c={cn:10s}: cRc=R {e1:.0e} | cRc=R21 {e2:.0e} | R^dag=cR^-1c {e3:.0e} | R^dag=cR21^-1c {e4:.0e}")
xs = [0.37, 1.9, 7.3]
q = q_of_xi(0.9013); Rp = setup(q)
A = np.fliplr(np.eye(3))
tests("IK", Rp, 3, {"1": np.eye(3), "diag(1,-1,1)": np.diag([1., -1, 1]), "antidiag": A, "antidiag(1,-1,1)": A @ np.diag([1., -1, 1])}, xs)
for n in (2, 3):
    D, Wn, setup2 = sc.make(n, 'c'); Rs = setup2(np.exp(-1j*np.pi*1.61))
    A = np.fliplr(np.eye(D)); g = np.array([(-1.0)**bin(i).count('1') for i in range(D)])
    tests(f"c_{n} spinor", Rs, D, {"1": np.eye(D), "grading": np.diag(g), "antidiag": A, "antidiag*grad": A @ np.diag(g)}, xs)
