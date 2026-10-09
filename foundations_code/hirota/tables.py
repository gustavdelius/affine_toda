"""Print, for each simply-laced algebra and soliton: zeros/poles of X_ab in the bound-state half-plane,
the net count per energy, isolated vs embedded, multiple zeros, threshold zeros/poles (z=0)."""
import sys
import numpy as np
from lie import Algebra
from compare_dorey import X_dorey
from xtable import soliton_table

for nm in sys.argv[1:]:
    A = Algebra(nm[0], int(nm[1:]))
    Xd = X_dorey(A)
    mmin2 = min(A.mass[b] ** 2 for b in A.mass)
    print(f"=== {nm}  h={A.h}  m^2: " + ", ".join(f"{b}:{A.mass[b]**2:.4f}" for b in sorted(A.mass)))
    for a in sorted(A.mass):
        Ns, groups, thr = soliton_table(A, a, Xd)
        iso, emb, flags = [], [], []
        for lam, net, g in groups:
            det = " ".join(f"X{b}:q{q}^{m:+d}" for (_, b, q, m) in g)
            s = f"{lam:.6f}[net {net:+d}; {det}]"
            if lam < mmin2 - 1e-9:
                iso.append(s)
            else:
                emb.append(s)
            if net < 0:
                flags.append(f"NEGATIVE net at {lam:.5f}")
            if any(abs(m) > 1 for (_, b, q, m) in g):
                flags.append(f"multiple zero/pole at {lam:.5f}: {det}")
            if len(g) > 1 and any(m < 0 for (_, b, q, m) in g) and any(m > 0 for (_, b, q, m) in g):
                flags.append(f"cancellation at {lam:.5f}: {det}")
        print(f" soliton {a}:")
        print("   isolated:", "; ".join(iso))
        print("   embedded:", "; ".join(emb))
        if thr:
            print("   threshold (z=0):", thr)
        for f in flags:
            print("   FLAG", f)
