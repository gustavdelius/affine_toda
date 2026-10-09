import numpy as np
from a22cubic import allev
lxs = np.concatenate([np.linspace(1e-3, 6, 600), np.linspace(6, 40, 100)])
res = []
for phi in np.arange(0.0025, 1.0, 0.005):     # q = e^{i pi phi}; spectrum symmetric under q -> conj q (x real) and q -> -q ?
    q = np.exp(1j*np.pi*phi)
    dp = max(np.abs(np.abs(allev(np.exp(l), q)) - 1).max() for l in lxs)
    res.append((phi, dp))
# print intervals of unimodularity
good = [p for p, d in res if d < 1e-6]
print("phi (q=e^{i pi phi}) with unimodular x>0 spectrum:")
iv = []; 
for p in good:
    if iv and abs(p - iv[-1][1] - 0.005) < 1e-9: iv[-1][1] = p
    else: iv.append([p, p])
print([(round(a, 4), round(b, 4)) for a, b in iv])
# symmetry checks
for phi in (0.13, 0.37):
    for qq in (np.exp(1j*np.pi*phi), np.exp(-1j*np.pi*phi), -np.exp(1j*np.pi*phi), 1j*np.exp(1j*np.pi*phi)):
        print(phi, np.round(qq, 3), max(np.abs(np.abs(allev(np.exp(l), qq)) - 1).max() for l in lxs[::5]))
