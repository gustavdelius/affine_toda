import numpy as np, time, os
src = open('f4breather.py').read().split("tB = T - 1; d_ = 1e-7")[0]
exec(src)
t0 = time.time()
def Rx(x): return Rmatrix(solve_fast(x))
found = []
for kph in range(12):
    eps = np.exp(1j*np.pi*kph/6)
    for l in range(0, 31):
        xs = eps*q**(-l)
        try:
            n1 = np.linalg.norm(Rx(xs*(1 + 1e-5))); n2 = np.linalg.norm(Rx(xs*(1 + 1e-7)))
        except Exception: continue
        if n2/n1 > 30:
            sv = np.linalg.svd(1e-7*Rx(xs*(1 + 1e-7)), compute_uv=False); found.append((f"e^(i pi {kph}/6) q^-{l}", int(np.sum(sv > 1e-6*sv[0]))))
print("27 x 27 R-matrix poles (12th roots of unity, l = 0..30):", found, f"[{time.time()-t0:.0f}s]")
# rank drops without poles (zero-type special points) along the same set
drops = []
for kph in range(12):
    eps = np.exp(1j*np.pi*kph/6)
    for l in range(0, 31):
        xs = eps*q**(-l); M = Rx(xs); sv = np.linalg.svd(M, compute_uv=False)
        r = int(np.sum(sv > 1e-8*sv[0]))
        if r < N2: drops.append((f"e^(i pi {kph}/6) q^-{l}", r))
print("rank drops of the top-normalised R (zero-type special points):", drops, f"[{time.time()-t0:.0f}s]")
