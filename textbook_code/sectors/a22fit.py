import numpy as np
from a22 import rep, W
from rsolve import intertwiner
om = 0.13; q = np.exp(-1j*np.pi*om)
x0 = 1.37+0.41j
R0 = intertwiner(rep(x0, q), rep(1.0, q), (0, 1), W, W)[0]; R0 /= R0[0, 0]
ev = np.linalg.eigvals(R0); grp = []
for v in ev:
    if not any(abs(v-g[0]) < 1e-7 for g in grp): grp.append([v, 0])
    for g in grp:
        if abs(v - g[0]) < 1e-7: g[1] += 1
print([(np.round(g[0], 4), g[1]) for g in grp])
Pk = []
for g in grp:
    M = np.eye(9, dtype=complex)
    for h in grp:
        if h is not g: M = M @ (R0 - h[0]*np.eye(9))/(g[0] - h[0])
    Pk.append(M)
xs = 1.3*np.exp(2j*np.pi*(np.arange(16)+0.17)/16)
rhos = []
for x in xs:
    R = intertwiner(rep(x, q), rep(1.0, q), (0, 1), W, W)[0]; R /= R[0, 0]
    rhos.append([np.trace(Pm @ R)/np.trace(Pm) for Pm in Pk])
rhos = np.array(rhos)
def label(z):
    s = np.angle(z)/np.pi; best = None
    for L in range(-24, 25):
        for sig, nm in ((0, '+'), (1, '-'), (0.5, '+i'), (-0.5, '-i')):
            d_ = (s - (sig - L*om)) % 2; d_ = min(d_, 2 - d_)
            if best is None or d_ < best[0]: best = (d_, L, nm)
    return f"{best[2]}q^{best[1]:+d}" + ("" if best[0] < 1e-6 and abs(abs(z) - 1) < 1e-6 else f"(?{abs(z):.3f})")
for k in range(len(Pk)):
    vals = rhos[:, k]
    for d in range(0, 5):
        Vm = np.vander(xs, d+1, increasing=True); Am = np.hstack([Vm, -(vals[:, None]*Vm[:, 1:])])
        sol, *_ = np.linalg.lstsq(Am, vals, rcond=None); num, den = sol[:d+1], np.concatenate([[1], sol[d+1:]])
        if np.abs(np.polyval(num[::-1], xs)/np.polyval(den[::-1], xs) - vals).max() < 1e-8: break
    print("channel dim", int(round(np.trace(Pk[k]).real)), "zeros", [label(z) for z in np.roots(num[::-1])], "poles", [label(z) for z in np.roots(den[::-1])])
