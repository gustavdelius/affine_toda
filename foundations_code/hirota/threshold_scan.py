"""(i) zeros/poles of X_ab at its own threshold (z=0, Lemma 8.4); (ii) zeros of X_ab in the bound-state
half-plane that sit exactly at the threshold of another channel (threshold resonances, cf. Prop 8.7);
(iii) spectral singularities: zeros at z in i*R\{0} (none possible: all zeros/poles at z = cos(pi q/h))."""
import numpy as np
from lie import Algebra
from compare_dorey import X_dorey
for nm in ['a2','a3','a4','a5','a6','a7','d4','d5','d6','d7','d8','e6','e7','e8']:
    A = Algebra(nm[0], int(nm[1:]))
    Xd = X_dorey(A)
    h = A.h
    thr0, thrc = [], []
    for a in sorted(A.mass):
        for b in sorted(A.mass):
            for q, m in Xd[(a, b)].items():
                if 2 * q == h:
                    thr0.append(f"X_{a}{b}:{'zero' if m>0 else 'pole'}^{abs(m)}")
                elif 2 * q < h and m > 0:
                    lam = A.mass[b] ** 2 * np.sin(np.pi * q / h) ** 2
                    cs = [c for c in A.mass if abs(A.mass[c] ** 2 - lam) < 1e-9]
                    if cs:
                        thrc.append(f"a={a}: X_{a}{b} zero (q={q}) at m_{cs}^2={lam:.4f}")
    print(f"{nm}: own-threshold zeros/poles: {', '.join(thr0) if thr0 else 'none'}")
    print(f"     zeros at another channel's threshold: {'; '.join(thrc) if thrc else 'none'}")
