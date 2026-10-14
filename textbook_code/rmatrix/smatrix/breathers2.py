import numpy as np, time, sys
src = open('breathers.py').read()
src = src.replace("J = 1500", "J = 1000")
exec(src.split("results = {}")[0])
def strip_scan(A_, b):
    offs = [o for _, o in swaps(A_, b)]
    cands = sorted({round(s*(tp - T*o/np.pi), 9) for tp in cand33 for o in offs for s in (1, -1)})
    cands = [c for c in cands if 1e-6 < c < T - 1e-6]
    pz = []
    for t0 in cands:
        n = [abs(scalar(A_, b, 1j*np.pi*(t0 + d)/T)) for d in (1e-3, 1e-5)]
        od = np.log10(n[1]/n[0])/2
        if abs(od) > 0.5: pz.append((t0, int(round(od))))
    return pz
def model2(pz, theta):
    out = 1.0
    for y, od in pz:
        yy = y if od > 0 else -y
        for _ in range(abs(od)): out *= np.sinh(theta/2 + 1j*np.pi*yy/(2*T))/np.sinh(theta/2 - 1j*np.pi*yy/(2*T))
    return out
todo = [tuple(int(c) for c in p) for p in sys.argv[1:]]
for a, b in todo:
    t1 = time.time(); A_ = f'B{a}'
    pz = strip_scan(A_, b)
    r = np.array([scalar(A_, b, th)/model2(pz, th) for th in (0.37+0.21j, -0.8+1.3j, 1.1-0.4j)])
    sgn = r.mean(); spread = np.abs(r - sgn).max()
    cross = scalar(A_, b, 1j*np.pi - (0.3+0.2j))/scalar(A_, b, 0.3+0.2j)
    print(f"S[B{a}, soliton {b}] = {sgn.real:+.3f} x blocks:  poles {[y for y, o in pz if o > 0]}  zeros {[y for y, o in pz if o < 0]}"
          f"  (orders {[o for _, o in pz]}; match {spread:.0e}; crossing {cross.real:+.4f}{cross.imag:+.4f}i) [{time.time()-t1:.0f}s]")
