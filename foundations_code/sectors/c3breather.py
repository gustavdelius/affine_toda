import numpy as np
from c3fuse import *
sv = add_bound('B3', 3, 3, 3*omega, 1); print("B3 residue sv", np.round(sv[:3], 6))
sv = add_bound('B1', 1, 1, 3*omega + 0.5, 1); print("B1 residue sv", np.round(sv[:3], 6))
for a, b in (('B3', 3), ('B3', 1), ('B1', 3)):
    devs = []; offd = []
    for th in np.linspace(-1.5, 1.5, 13):
        M, res = Spl(a, b, th); s = np.trace(M)/M.shape[0]
        offd.append(np.abs(M - s*np.eye(M.shape[0])).max()); devs.append(abs(abs(s) - 1))
    print(f"S[{a},{b}] at real theta: scalar x identity to {max(offd):.1e}; max ||s|-1| = {max(devs):.1e}")
