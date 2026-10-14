import numpy as np
from bspin import *
x0 = 1.37+0.41j
def Rnum(x):
    Rc = intertwiner(rep(x), rep(1.0), range(N+1), Wn, Wn)[0]; return Rc/Rc[0, 0]
R0 = Rnum(x0); ev, V = np.linalg.eig(R0)
groups = []
for idx, v in enumerate(ev):
    for g in groups:
        if abs(g[0] - v) < 1e-7*max(1, abs(v)): g[1].append(idx); break
    else: groups.append([v, [idx]])
groups.sort(key=lambda g: -len(g[1]))
lam = [g[0] for g in groups]
Pk = []
for kk in range(len(lam)):
    M = np.eye(D*D, dtype=complex)
    for jj in range(len(lam)):
        if jj != kk: M = M @ (R0 - lam[jj]*np.eye(D*D))/(lam[kk] - lam[jj])
    Pk.append(M)
print("channels (dims):", [int(round(np.trace(P_).real)) for P_ in Pk])
def label(z):
    s = np.angle(z)/np.pi; best = None
    for L in range(-16, 17):
        for sig, nm in ((0, '+'), (1, '-'), (0.5, '+i'), (-0.5, '-i')):
            d_ = (s - (sig - L*omega)) % 2; d_ = min(d_, 2 - d_)
            if best is None or d_ < best[0]: best = (d_, L, nm)
    return f"{best[2]}q^{best[1]:+d}" + ("" if best[0] < 1e-6 and abs(abs(z) - 1) < 1e-6 else f"(?{abs(z):.3f})")
def ratfit(xs, vals, maxd=6):
    for d in range(0, maxd+1):
        Vm = np.vander(xs, d+1, increasing=True); Am = np.hstack([Vm, -(vals[:, None]*Vm[:, 1:])])
        sol, *_ = np.linalg.lstsq(Am, vals, rcond=None); num, den = sol[:d+1], np.concatenate([[1], sol[d+1:]])
        if np.abs(np.polyval(num[::-1], xs)/np.polyval(den[::-1], xs) - vals).max() < 1e-8*np.abs(vals).max(): return num, den
    return None
xs = 1.3*np.exp(2j*np.pi*(np.arange(12) + 0.17)/12)
Rs = [Rnum(x) for x in xs]
rho = np.array([[np.trace(P_ @ R_)/np.trace(P_) for P_ in Pk] for R_ in Rs])
for kk in range(1, len(Pk)):
    r = rho[:, kk]/rho[:, kk-1]; num, den = ratfit(xs, r)
    print(f"rho_{kk}/rho_{kk-1}: zeros {[label(z) for z in np.roots(num[::-1])]}  poles {[label(p) for p in np.roots(den[::-1])]}")
# crossing point and universal crossing factor
def dual(x):
    e, f, k, _ = rep(x); out = []
    for i in range(N+1):
        ki = np.linalg.inv(k[i]); out += [(-ki@e[i]).T, (-f[i]@k[i]).T, ki.T]
    return out
def plain(x):
    e, f, k, _ = rep(x); out = []
    for i in range(N+1): out += [e[i], f[i], k[i]]
    return out
xa = 0.83+0.29j
for ph, nm in ((1, '+'), (-1, '-'), (1j, '+i'), (-1j, '-i')):
    for l in range(-16, 17):
        tt = ph*q**l
        M = np.vstack([np.kron(np.eye(D), A_.T) - np.kron(B_, np.eye(D)) for A_, B_ in zip(dual(xa), plain(xa*tt))])
        s = np.linalg.svd(M, compute_uv=False)
        if s[-1] < 1e-9*s[0]:
            C = np.linalg.svd(M)[2][-1].conj().reshape(D, D); tc = tt; print(f"crossing point: V*(x) = V(x * {nm}q^{l})")
Ci = np.linalg.inv(C)
P = np.zeros((D*D, D*D))
for i in range(D):
    for j in range(D): P[j*D+i, i*D+j] = 1
def pt1(M_): return M_.reshape(D, D, D, D).transpose(2, 1, 0, 3).reshape(D*D, D*D)
cus = []
for x in xs:
    R_ = P @ Rnum(x); lhs = pt1(np.linalg.inv(R_)); rhs = np.kron(Ci, np.eye(D)) @ (P @ Rnum(x*tc)) @ np.kron(C, np.eye(D))
    c_ = np.vdot(rhs.ravel(), lhs.ravel())/np.vdot(rhs.ravel(), rhs.ravel()); cus.append(c_)
num, den = ratfit(xs, np.array(cus))
print(f"crossing factor c_u(x): zeros {[label(z) for z in np.roots(num[::-1])]}  poles {[label(p) for p in np.roots(den[::-1])]}")
np.save(f'bspin_P_n{N}_w{omega}.npy', np.array(Pk)); np.save(f'bspin_C_n{N}_w{omega}.npy', C); np.save(f'bspin_tc_n{N}_w{omega}.npy', np.array([tc]))
