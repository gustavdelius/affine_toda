import numpy as np, itertools
n = 3; ht = 2*n
eps = lambda a: 1 if a == n else 2
M = {a: eps(a)*np.sin(a*np.pi/(2*n)) for a in range(1, n+1)}                       # classical (1, sqrt3, 1)
r = {a: (1/(8*n))*(a/(2*n))*(np.cos(a*np.pi/(2*n))/np.sin(a*np.pi/(2*n))) for a in range(1, n+1)}   # dM/M per beta^2 (DG), non-universal part
# fusions: (a, b, c, l, list of allowed phases sigma from the R-matrix tables)
F = [(1,1,2,2,[0,1]), (1,2,1,5,[0,1]), (1,3,3,4,[0.5,-0.5]), (2,2,2,4,[0]), (2,3,3,5,[0.5,-0.5]), (3,3,2,2,[0]), (3,3,1,4,[1])]
rows = []
for a, b, c, l, sig in F:
    u0 = l*np.pi/ht
    assert abs(M[c]**2 - (M[a]**2 + M[b]**2 + 2*M[a]*M[b]*np.cos(u0))) < 1e-12
    du = (M[a]**2*r[a] + M[b]**2*r[b] + np.cos(u0)*M[a]*M[b]*(r[a] + r[b]) - M[c]**2*r[c])/(M[a]*M[b]*np.sin(u0))
    rows.append((a, b, c, l, sig, du))     # du = delta u / beta^2 (radians)
    print(f"{a}+{b}->{c} at u0 = {l}pi/6: one-loop shift delta u = {du:+.6f} beta^2   (delta u/pi = {du/np.pi:+.6f} beta^2)")
# model: delta u/pi = [sigma + 2k - l(sigma_c + 2 k1)/ht] * beta^2 / (ht * 4 pi * kappa) ; sigma_c = 0 for n=3
print("\nSearching kappa > 0, k1 in -3..3, integer k per fusion, sigma from the allowed phases:")
best = []
for k1 in range(-3, 4):
    # for each fusion, B(sigma,k) = sigma + 2k - l*2*k1/ht must equal  ht*4*kappa*du  (du in units of pi*beta^2 -> use du/pi)
    target = lambda du: ht*4*np.pi*(du/np.pi)          # = B / kappa ... B = kappa * target? careful below
    # delta u/pi = B/(ht*4*pi*kappa)  =>  B = ht*4*pi*kappa*(delta u/pi)
    for kappa in np.linspace(0.05, 3, 5901):
        ok, assign = True, []
        for a, b, c, l, sig, du in rows:
            need = ht*4*np.pi*kappa*(du/np.pi) + l*2*k1/ht      # = sigma + 2k
            cands = [(s, (need - s)/2) for s in sig]
            good = [(s, round(kk)) for s, kk in cands if abs(kk - round(kk)) < 2e-3]
            if not good: ok = False; break
            assign.append((a, b, c, good[0]))
        if ok: best.append((k1, kappa, assign))
seen = set()
for k1, kappa, assign in best:
    key = (k1, round(kappa, 2))
    if key in seen: continue
    seen.add(key)
    print(f"  k1 = {k1}, kappa = {kappa:.4f}: " + ", ".join(f"{a}{b}->{c}: sigma={s:+.1f}, k={k}" for a, b, c, (s, k) in assign))
if not best: print("  no consistent assignment found")
