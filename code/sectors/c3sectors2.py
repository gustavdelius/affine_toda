import numpy as np, time, os, sys
from c3fuse import *
t0 = time.time()
def krein_err(e): return max(np.min(np.abs(e - 1/np.conj(s))) for s in e)
M = {3: 1.0, 1: 2*np.cos(u1/2), 2: 2*np.cos(u2/2)}
print(f"omega={omega}, T={T:.4f}, masses M1={M[1]:.5f} M2={M[2]:.5f} M3=1", flush=True)
# excited solitons
sv = add_bound('2s', 1, 1, omega, 29); print("2* (1+1 at t=omega) residue singular values (29,30th rel):", np.round(sv[27:31], 6), flush=True)
ts = 2*omega + 0.25; us = np.pi*ts/T; Ms = np.sqrt(1 + M[1]**2 + 2*M[1]*np.cos(us))
u3 = np.arctan2(M[1]*np.sin(us), 1 + M[1]*np.cos(us)); u1b = us - u3
sv = add_bound('3s', 3, 1, ts, 8, ubar=(u3, u1b)); print(f"3* (3+1 at t=2w+1/4, M={Ms:.5f}) residue sv:", np.round(sv[6:10], 6), flush=True)
def scan2(a, b, lmax=14, npts=int(os.environ.get("NPTS", "12"))):
    """max over real theta of ||s|-1| for the particle-labelled S_ab; R-only eigenvalues times |prod F33|"""
    sym = sorted(round(o, 12) for _, o in swaps(a, b)) == sorted(round(-o, 12) for _, o in swaps(a, b))
    ths = np.concatenate([-np.linspace(lmax, 1e-3, npts), np.linspace(1e-3, lmax, npts)])/(2*T)
    worst = (0, None, 0); clos = 0; win = []
    for th in ths:
        Mx, res = Spl(a, b, th); clos = max(clos, res)
        e = np.linalg.eigvals(Mx)
        if not sym: e = e*abs(np.exp(scalar_logF(a, b, th)))
        d = np.abs(np.abs(e) - 1).max()
        if d < 1e-6: win.append(2*T*th)
        if d > worst[0]: worst = (d, 2*T*th, krein_err(e))
    lo = [l for l in win]
    return worst, clos, sym, (min(lo), max(lo)) if lo else None
pairs = [(3, 3), (3, 1), (1, 1), (3, 2), (1, 2), (2, 2), ('2s', 3), ('2s', 1), ('3s', 3), ('3s', 1)]
if len(sys.argv) > 1: pairs = [tuple(int(c) if c.isdigit() else c for c in p.split(',')) for p in sys.argv[1:]]
for a, b in pairs:
    (d, lx, kr), clos, sym, win = scan2(a, b)
    print(f"S[{a},{b}] (dims {dims(a)}x{dims(b)}; symmetric offsets {sym}): max ||s|-1| = {d:.2e} at log x = 2T theta = {lx:.2f}; Krein pairing err {kr:.0e}; "
          f"closure {clos:.0e}; unimodular window in log x: {None if win is None else (round(win[0],2), round(win[1],2))}  [{time.time()-t0:.0f}s]", flush=True)
