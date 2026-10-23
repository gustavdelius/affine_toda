import numpy as np, os
omega = 2.37; T = 3*omega + 1; H = 6 + 2/omega
M = {1: 2*np.cos(np.pi/4 + np.pi/(2*H)), 2: 2*np.cos(np.pi/H), 3: 1.0}
dimname = {1: "singlet", 8: "8-plet", 29: "29-plet"}
entries = []
for a, b in [(3, 3), (3, 1), (3, 2), (1, 1), (1, 2), (2, 2)]:
    f = f"poles_{a}{b}.npy"
    if not os.path.exists(f): continue
    for t0, od, rk, ms, ns, mu_, nu in np.load(f, allow_pickle=True):
        if od == 1 and rk in dimname:
            entries.append((float(ms), int(rk), f"{a}+{b}", float(t0)))
entries.sort()
groups = []
for e in entries:
    for g in groups:
        if abs(g[0] - e[0]) < 1e-5 and g[1] == e[1]: g[2].append((e[2], e[3])); break
    else: groups.append([e[0], e[1], [(e[2], e[3])]])
m_part = {a: 2*M[a]*np.sin(np.pi/(2*T)*(0.5 if a < 3 else 1)) for a in (1, 2, 3)}
print(f"exact soliton masses (M3 = 1): M1 = {M[1]:.5f}, M2 = {M[2]:.5f};  lightest breathers: {m_part}")
print(f"{'mass':>9s}  {'multiplet':9s}  channels (t)")
for m, rk, ch in groups:
    tag = ""
    for a in (1, 2, 3):
        if abs(m - M[a]) < 1e-6 and rk == (8 if a != 2 else 29): tag = f"  <- soliton {a}"
    if rk == 29 and abs(m - 2*M[1]*np.cos(np.pi*omega/(2*T))) < 1e-6: tag = "  <- 2* (excited soliton 2)"
    multi = "  [in %d channels]" % len({c for c, _ in ch}) if len({c for c, _ in ch}) > 1 else ""
    print(f"{m:9.5f}  {dimname[rk]:9s}  " + ", ".join(f"{c}({t:.3f})" for c, t in ch) + tag + multi)
