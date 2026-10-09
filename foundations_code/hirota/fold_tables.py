"""Zero/pole bookkeeping for the folded theories (from data/folds_all.json written by fold.py)."""
import json, numpy as np
D = json.load(open('data/folds_all.json'))
for nm, o in D.items():
    h = o['h']
    m2 = {int(b): v['mass'] ** 2 for b, v in o['species'].items()}
    mmin2 = min(m2.values())
    print(f"=== {nm} (h={h}, central={o['central']}) channel m^2: " + ", ".join(f"{b}:{v:.4f}" for b, v in sorted(m2.items())) + f"  [members {', '.join(f'{b}:{v['members']}' for b, v in o['species'].items())}]")
    sols = sorted(m2)
    for a in sols:
        ent, thr, mult = [], [], []
        for b in sols:
            N = {int(q): m for q, m in o['X'][f"{a},{b}"]['N'].items()}
            for q, m in N.items():
                if 2 * q < h:
                    lam = m2[b] * np.sin(np.pi * q / h) ** 2
                    ent.append((lam, b, q, m))
                    if abs(m) > 1:
                        mult.append((round(lam, 6), b, q, m))
                elif 2 * q == h:
                    thr.append((b, m))
        ent.sort()
        groups = []
        for e in ent:
            if groups and abs(groups[-1][0] - e[0]) < 1e-9:
                groups[-1][1] += e[3]; groups[-1][2].append(e)
            else:
                groups.append([e[0], e[3], [e]])
        iso = [f"{g[0]:.5f}[{g[1]:+d}]" + ("*" if len(g[2]) > 1 or abs(g[2][0][3]) > 1 else "") for g in groups if g[0] < mmin2 - 1e-9]
        emb = [f"{g[0]:.5f}[{g[1]:+d}]" + ("@thr" if any(abs(g[0]-v) < 1e-9 for v in m2.values()) else "") for g in groups if g[0] >= mmin2 - 1e-9]
        print(f"  soliton {a}: isolated {'; '.join(iso)} | embedded {'; '.join(emb)}" + (f" | threshold zeros/poles {thr}" if thr else "") + (f" | multiple: {mult}" if mult else ""))
