import numpy as np
from common import *
omega = 2.37; q = np.exp(-1j*np.pi*omega)          # physical regime: q a pure phase
S = Spinor33(q)
print(f"q = exp(-i pi {omega}); projector ranks {S.ranks}; closed form reproduces direct solve to {S.check:.1e}")
# antipode for Delta(e)=e(x)1+k(x)e, Delta(f)=f(x)k^-1+1(x)f :  S(e)=-k^-1 e, S(f)=-f k, S(k)=k^-1
def dual(x):
    e, f, k, _ = spin_rep(x, q); out = []
    for i in range(n+1):
        ki = np.linalg.inv(k[i])
        out += [(-ki@e[i]).T, (-f[i]@k[i]).T, ki.T]
    return out
def plain(x):
    e, f, k, _ = spin_rep(x, q); out = []
    for i in range(n+1): out += [e[i], f[i], k[i]]
    return out
x = 0.83+0.29j
found = []
for ph, lab in ((1, '+'), (-1, '-')):
    for l in range(-8, 9):
        t = ph*q**l
        A, B = dual(x), plain(x*t)
        M = np.vstack([np.kron(np.eye(8), Ag.T) - np.kron(Bg, np.eye(8)) for Ag, Bg in zip(A, B)])
        s = np.linalg.svd(M, compute_uv=False)
        if s[-1] < 1e-9*s[0]: found.append((lab, l, t))
print("V*(x) = V(x t) for t =", [f"{lab}q^{l}" for lab, l, _ in found])
lab, l, t = found[0]
A, B = dual(x), plain(x*t)
M = np.vstack([np.kron(np.eye(8), Ag.T) - np.kron(Bg, np.eye(8)) for Ag, Bg in zip(A, B)])
C = np.linalg.svd(M)[2][-1].conj().reshape(8, 8)
def pt1(Mm):
    T = Mm.reshape(8, 8, 8, 8); return T.transpose(2, 1, 0, 3).reshape(64, 64)
Ci = np.linalg.inv(C)
def cu(xx):
    R = P @ S.R(xx); lhs = pt1(np.linalg.inv(R)); rhs = np.kron(Ci, np.eye(8)) @ (P @ S.R(xx*t)) @ np.kron(C, np.eye(8))
    c = np.vdot(rhs.ravel(), lhs.ravel())/np.vdot(rhs.ravel(), rhs.ravel())
    return c, np.abs(lhs - c*rhs).max()/np.abs(lhs).max()
for xx in (0.7+0.2j, 1.9-0.4j):
    c, res = cu(xx); print(f"universal crossing  (R(x)^-1)^t1 = c(x) (C^-1 x 1) R(x t) (C x 1):  residual {res:.1e} at x={xx}")
# fit c(x) as a rational function of x
xs = 1.3*np.exp(2j*np.pi*(np.arange(30)+0.17)/30); cs = np.array([cu(v)[0] for v in xs])
for d in range(1, 6):
    Am = np.hstack([np.vander(xs, d+1, increasing=True), -(cs[:, None]*np.vander(xs, d+1, increasing=True)[:, 1:])])
    sol, *_ = np.linalg.lstsq(Am, cs, rcond=None); a, b = sol[:d+1], np.concatenate([[1], sol[d+1:]])
    res = np.abs(np.polyval(a[::-1], xs)/np.polyval(b[::-1], xs) - cs).max()/np.abs(cs).max()
    if res < 1e-9: break
zs, ps = np.roots(a[::-1]), np.roots(b[::-1])
def lab_q(z):
    l = np.log(z/(np.abs(z)))/(1j*np.pi)            # z = e^{i pi s}  (|z|=1 here) ; express s = sigma - l*omega
    s = (np.angle(z)/np.pi)
    # find integer-ish l and sigma in {0,1}: s = sigma - l*omega mod 2
    best = None
    for L in range(-12, 13):
        for sig in (0, 1, 0.5, -0.5):
            d_ = (s - (sig - L*omega)) % 2
            d_ = min(d_, 2 - d_)
            if best is None or d_ < best[0]: best = (d_, L, sig)
    return f"{'+' if best[2]==0 else ('-' if best[2]==1 else ('+i' if best[2]==0.5 else '-i'))}q^{-best[1]:+d}"
print(f"c(x) rational, degree {d}, fit residual {res:.1e};  zeros at x = {[lab_q(z) for z in zs]};  poles at x = {[lab_q(p) for p in ps]};  |zeros| {np.round(np.abs(zs),6)}")
np.save('C.npy', C); np.save('t.npy', np.array([t]))
