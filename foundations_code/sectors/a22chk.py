import numpy as np
from a22 import rep, W
from rsolve import intertwiner
P = np.zeros((9, 9))
for i in range(3):
    for j in range(3): P[j*3+i, i*3+j] = 1
om = 0.13; q = np.exp(-1j*np.pi*om)
for x in [2.0, 1.3+0.4j]:
    R, null, sv, M = intertwiner(rep(x, q), rep(1.0, q), (0, 1), W, W)
    print("null dim", null, sv, M)
    R = R/R[0, 0]
    ev = np.linalg.eigvals(R); print("Rcheck eig", np.round(ev, 6))
    # compare TW forms with Q = q^a
    for name, Q in [("q", q), ("1/q", 1/q), ("q^1/2", q**0.5), ("q^-1/2", q**-0.5)]:
        for X in (x, 1/x):
            cand = [1, (X*Q**4 - 1)/(X - Q**4), (X*Q**6 + 1)/(X + Q**6)]
            if all(np.min(np.abs(ev - c)) < 1e-8 for c in cand): print("   match TW with Q =", name, "X =", "x" if X == x else "1/x")
    print("PR eig", np.round(np.linalg.eigvals(P @ R), 5), np.round(np.abs(np.linalg.eigvals(P @ R)), 5))
