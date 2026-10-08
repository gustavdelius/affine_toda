import numpy as np, time, os
from scipy.linalg import null_space, orth
exec(open('f4rep.py').read().split('if __name__ == "__main__":')[0])
E, F = np.load(f"f4rep_w{omega}.npy", allow_pickle=True)
t0 = time.time(); I = np.eye(n)
Kd = {i: K(i) for i in range(5)}; Kinv = {i: np.linalg.inv(Kd[i]) for i in range(5)}
dE = {i: np.kron(E[i], I) + np.kron(Kd[i], E[i]) for i in range(1, 5)}        # Delta(e) = e x 1 + k x e
dF = {i: np.kron(F[i], Kinv[i]) + np.kron(I, F[i]) for i in range(1, 5)}       # Delta(f) = f x k^-1 + 1 x f
WW = (W[:, None, :] + W[None, :, :]).reshape(n*n, 4)
def wkey(w): return tuple(np.round(np.asarray(w)*2).astype(int))
wsp = {}
for idx, w in enumerate(WW): wsp.setdefault(wkey(w), []).append(idx)
# highest-weight vectors of the components
comps = []
for key, idx in wsp.items():
    M = np.vstack([dE[i][:, idx] for i in range(1, 5)])
    ns = null_space(M, rcond=1e-10)
    for c_ in range(ns.shape[1]):
        v = np.zeros(n*n, complex); v[idx] = ns[:, c_]; comps.append((key, v))
print(f"{len(comps)} highest-weight vectors at weights {[k for k, _ in comps]}  [{time.time()-t0:.0f}s]")
# component weight spaces by descending with Delta(f)
simple_keys = {i: wkey(simple[i]) for i in range(1, 5)}
def component_spaces(v, key0, maxdepth=40):
    spaces = {key0: v[:, None]/np.linalg.norm(v)}
    frontier = [key0]
    for _ in range(maxdepth):
        new = []
        for key in frontier:
            B = spaces[key]
            for i in range(1, 5):
                k2 = tuple(a - b for a, b in zip(key, simple_keys[i]))
                if k2 not in wsp: continue
                Y = dF[i] @ B
                if np.linalg.norm(Y) < 1e-12: continue
                cur = spaces.get(k2)
                Z = Y if cur is None else np.hstack([cur, Y])
                O = orth(Z, rcond=1e-10)
                if cur is None or O.shape[1] > cur.shape[1]:
                    spaces[k2] = O; new.append(k2)
        if not new: break
        frontier = list(set(new))
    return spaces
CS = [component_spaces(v, key) for key, v in comps]
dimsC = [sum(B.shape[1] for B in cs.values()) for cs in CS]
print("component dimensions from descent:", dimsC, f" total {sum(dimsC)}  [{time.time()-t0:.0f}s]")
order = np.argsort(dimsC)[::-1]; comps = [comps[i] for i in order]; CS = [CS[i] for i in order]; dimsC = [dimsC[i] for i in order]
th_key = wkey(-theta)
def decompose(w, key):
    """split a vector of weight `key` into its component pieces"""
    blocks = [(j, cs[key]) for j, cs in enumerate(CS) if key in cs]
    Bm = np.hstack([B for _, B in blocks]); coef, res, *_ = np.linalg.lstsq(Bm, w, rcond=None)
    pieces = {}; k = 0
    for j, B in blocks: pieces[j] = B @ coef[k:k+B.shape[1]]; k += B.shape[1]
    return pieces, np.linalg.norm(Bm @ coef - w)
def rho_at(x):
    dE0a = np.kron(x*E[0], I) + np.kron(Kd[0], E[0]); dE0b = np.kron(E[0], I) + np.kron(Kd[0], x*E[0])
    rows = []
    for k_, (key, v) in enumerate(comps):
        w = dE0a @ v
        if np.linalg.norm(w) < 1e-12: continue
        key2 = tuple(a + b for a, b in zip(key, th_key))
        pieces, err = decompose(w, key2)
        rhs = dE0b @ v
        # sum_j rho_j w_j - rho_k rhs = 0
        rowblock = np.zeros((n*n, len(comps)), complex)
        for j, p in pieces.items(): rowblock[:, j] += p
        rowblock[:, k_] -= rhs
        rows.append(rowblock)
    Mx = np.vstack(rows)
    sol, *_ = np.linalg.lstsq(Mx[:, 1:], -Mx[:, 0], rcond=None)
    r = np.concatenate([[1.0], sol]); return r, np.abs(Mx @ r).max()/np.abs(Mx).max()
r, res = rho_at(0.83+0.29j); print(f"eigenvalues from the e_0 intertwining equations: residual {res:.1e}  [{time.time()-t0:.0f}s]")
def label(z):
    s = np.angle(z)/np.pi; best = None
    for L in range(-24, 25):
        for sig, nm in ((0, '+'), (1, '-'), (0.5, '+i'), (-0.5, '-i')):
            d_ = (s - (sig - L*omega)) % 2; d_ = min(d_, 2 - d_)
            if best is None or d_ < best[0]: best = (d_, L, nm)
    return f"{best[2]}q^{best[1]:+d}" + ("" if best[0] < 1e-6 and abs(abs(z) - 1) < 1e-6 else f"(?{abs(z):.3f})")
def ratfit(xs, vals, maxd=10):
    for d in range(0, maxd+1):
        Vm = np.vander(xs, d+1, increasing=True); Am = np.hstack([Vm, -(vals[:, None]*Vm[:, 1:])])
        sol, *_ = np.linalg.lstsq(Am, vals, rcond=None); num, den = sol[:d+1], np.concatenate([[1], sol[d+1:]])
        if np.abs(np.polyval(num[::-1], xs)/np.polyval(den[::-1], xs) - vals).max() < 1e-8*np.abs(vals).max(): return num, den
    return None
xs = 1.3*np.exp(2j*np.pi*(np.arange(24) + 0.17)/24)
rhos = np.array([rho_at(x)[0] for x in xs])
for kk in range(1, len(comps)):
    out = ratfit(xs, rhos[:, kk])
    if out is None: print(f"rho[{dimsC[kk]}]: no rational fit"); continue
    num, den = out; zs, ps = list(np.roots(num[::-1])), list(np.roots(den[::-1]))
    for z in list(zs):
        m = [p_ for p_ in ps if abs(p_ - z) < 1e-6]
        if m: zs.remove(z); ps.remove(m[0])
    print(f"rho[{dimsC[kk]:3d}]/rho[{dimsC[0]}]: zeros {sorted(label(z) for z in zs)}  poles {sorted(label(p_) for p_ in ps)}")
np.save(f"f4_rho_setup_w{omega}.npy", np.array([dimsC], dtype=object), allow_pickle=True)
print(f"[{time.time()-t0:.0f}s]")
