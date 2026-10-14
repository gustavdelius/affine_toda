import numpy as np, sys, time
from allS import *
H = 6 + 2/omega
Mass = {1: 2*np.cos(np.pi/4 + np.pi/(2*H)), 2: 2*np.cos(np.pi/H), 3: 1.0}
# reference spectrum (units M_3 = 1)
states = {f"M{a}": Mass[a] for a in (1, 2, 3)}
states["M2*"] = 2*Mass[1]*np.cos(np.pi*omega/(2*T))
for a in (1, 2, 3):
    for p2 in range(1, 2*int(T)+2):
        p = p2/2; states[f"B({a})_{p:g}"] = 2*Mass[a]*np.sin(np.pi*p/(2*T))
def name(m):
    hits = [k for k, v in states.items() if abs(v - m) < 1e-6]; return ",".join(hits) if hits else "-"
# superset of real poles of S_33 in t'
cand33 = set()
for j in (1, 2):
    Aj = 2*T*(j-1)
    for a in (omega, 2*omega + 0.5, 3*omega):
        for k in range(0, int(4*T) + 3):
            cand33 |= {a - k - Aj, -1 - a - k - Aj, Aj + T - a + k, 1 + Aj + T + a + k}
cand33 = np.array(sorted(cand33))
rng = np.random.default_rng(1)
def scan(a, b):
    offs = [o for _, o in swaps(a, b)]
    cands = sorted({round(tp - T*o/np.pi, 9) for tp in cand33 for o in offs})
    cands = [c for c in cands if 1e-6 < c < T - 1e-6]
    D = dims[a]*dims[b]; Xr = rng.normal(size=(D, 3)) + 1j*rng.normal(size=(D, 3))
    out = []
    for t0 in cands:
        n = [np.linalg.norm(apply_S(a, b, 1j*np.pi*(t0 + d)/T, Xr)[0]) for d in (1e-4, 1e-6)]
        order = np.log10(n[1]/n[0])/2
        if order > 0.5: out.append((t0, int(round(order))))
    rows = []
    for t0, od in out:
        d = 1e-6
        Xp = np.eye(D, dtype=complex) if D <= 300 else (rng.normal(size=(D, 40)) + 1j*rng.normal(size=(D, 40)))
        Mx, _ = apply_S(a, b, 1j*np.pi*(t0 + d)/T, Xp)
        sv = np.linalg.svd(d**od*Mx, compute_uv=False); rk = int(np.sum(sv > 1e-3*sv[0]))
        u = np.pi*t0/T
        ms = np.sqrt(Mass[a]**2 + Mass[b]**2 + 2*Mass[a]*Mass[b]*np.cos(u))
        mu_ = np.sqrt(max(Mass[a]**2 + Mass[b]**2 - 2*Mass[a]*Mass[b]*np.cos(u), 0))
        rows.append((t0, od, rk, ms, name(ms), mu_, name(mu_)))
    return rows, len(cands)
pairs = [tuple(int(c) for c in p) for p in sys.argv[1:]] or [(3, 3), (3, 1), (3, 2)]
for a, b in pairs:
    t1 = time.time(); rows, nc = scan(a, b)
    np.save(f"poles_{a}{b}.npy", np.array(rows, dtype=object), allow_pickle=True)
    print(f"\n=== S_{a}{b}   (dims {dims[a]}x{dims[b]}; {nc} candidates tested; {len(rows)} poles; {time.time()-t1:.0f}s)")
    print("  (showing poles with multiplet-dimension residues or order > 1)\n      t      order  rank   s-channel mass  [state]          crossed mass  [state]")
    for t0, od, rk, ms, ns, mu_, nu in rows:
        if not (rk in (1, 8, 29) or od > 1): continue
        print(f"  {t0:8.4f}    {od}    {rk:4d}    {ms:8.5f}  {ns:16s}   {mu_:8.5f}  {nu}")
