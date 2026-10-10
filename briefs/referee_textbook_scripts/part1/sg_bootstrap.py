"""Bootstrap of the sine-Gordon breathers from eq-sg-smatrix / eq-sg-s0 (book's conventions).

Slot-labelled two-particle S-matrix: S_{ab}^{cd}(theta) maps |a(theta1) b(theta2)> to |c(theta1) d(theta2)>,
as in eq-zf-algebra (c carries theta1).  s=0, sbar=1.
S_ss^ss = S_0, S_{s sbar}^{s sbar} = S_T, S_{s sbar}^{sbar s} = S_R (and C-images).
"""
import itertools, functools
import numpy as np
from mpmath import mp, mpf, mpc, pi, sinh, sin, cos, tan, exp
from sg_s0 import S0 as S0g
mp.dps = 15

@functools.lru_cache(maxsize=None)
def S0c(th, lam):
    return complex(S0g(mpc(th), mpf(lam)))

def Smat(th, lam):
    s0 = S0c(complex(th), lam)
    lamm = mpf(lam); thm = mpc(th)
    st = complex(sinh(lamm*thm)/sinh(lamm*(1j*pi - thm)))*s0
    sr = complex(sinh(1j*pi*lamm)/sinh(lamm*(1j*pi - thm)))*s0
    M = np.zeros((4, 4), complex)       # M[(c,d),(a,b)] = S_ab^cd
    idx = lambda a, b: 2*a + b
    M[idx(0,0), idx(0,0)] = M[idx(1,1), idx(1,1)] = s0
    M[idx(0,1), idx(0,1)] = M[idx(1,0), idx(1,0)] = st
    M[idx(1,0), idx(0,1)] = M[idx(0,1), idx(1,0)] = sr
    return M

def embed(M, i, j, n):
    """act with 4x4 M on slots i,j of (C^2)^{n}"""
    D = 2**n; out = np.zeros((D, D), complex)
    for col in range(D):
        bits = [(col >> (n-1-k)) & 1 for k in range(n)]
        a, b = bits[i], bits[j]
        for c, d in itertools.product(range(2), repeat=2):
            amp = M[2*c+d, 2*a+b]
            if amp == 0: continue
            nb = bits.copy(); nb[i] = c; nb[j] = d
            row = sum(bit << (n-1-k) for k, bit in enumerate(nb))
            out[row, col] += amp
    return out

import sys
lam = float(sys.argv[1]) if len(sys.argv) > 1 else 1.5
xi = float(pi)/lam
print(f'lambda = {lam}, xi = {xi:.6f}')
# 1. Yang-Baxter in the slot-labelled convention
t1, t2, t3 = 0.31+0.1j, -0.42+0.05j, 0.9-0.2j
S12 = embed(Smat(t1-t2, lam), 0, 1, 3); S13 = embed(Smat(t1-t3, lam), 0, 2, 3); S23 = embed(Smat(t2-t3, lam), 1, 2, 3)
print('YBE S12 S13 S23 = S23 S13 S12:', np.abs(S12@S13@S23 - S23@S13@S12).max())

# 2. residue at the first breather pole theta = i(pi - xi), u = pi - xi
u = float(pi) - xi
eps = 1e-6
Res = (Smat(1j*u + eps, lam) - Smat(1j*u - eps, lam))*eps/2      # symmetric: (S(+e)-S(-e))/2 * e  ~ residue
U, sv, Vh = np.linalg.svd(Res)
print('singular values of residue matrix:', np.round(sv, 8))
v = U[:, 0]                                    # image of the residue (bound-state vector in slots 2,3)
print('bound state vector (basis ss, s sbar, sbar s, sbar sbar):', np.round(v/v[1], 8))
print('Res S_T =', Res[1,1], ' Res S_R =', Res[2,1])

# 3. soliton - B1 amplitude by fusion: particle 1 = soliton/antisoliton, (2,3) = B1 at theta_c +- i u/2
def S_sB(th):
    A = embed(Smat(th - 1j*u/2, lam), 0, 1, 3) @ embed(Smat(th + 1j*u/2, lam), 0, 2, 3)   # S12 S13
    vals = []
    for a in range(2):
        e = np.zeros(2); e[a] = 1
        x = np.kron(e, v); y = A @ x
        c = np.vdot(x, y)/np.vdot(x, x)
        vals.append((c, np.linalg.norm(y - c*x)))
    return vals
fx = lambda x, th: (np.sinh(th) + 1j*np.sin(np.pi*x))/(np.sinh(th) - 1j*np.sin(np.pi*x))
for th in (0.37+0.11j, -0.6+0.3j):
    vals = S_sB(th)
    print(f'S_(s,B1)({th}): s: {vals[0][0]:.10f} (resid {vals[0][1]:.1e}),  sbar: {vals[1][0]:.10f} (resid {vals[1][1]:.1e});'
          f' f_(1/2 - xi/2pi) = {fx(0.5 - xi/(2*np.pi), th):.10f}')

# 4. B1 - B1 amplitude by fusing four solitons: (1,2) = B1 at theta_A +- iu/2, (3,4) = B1 at theta_B +- iu/2
def S_BB(th):
    S = lambda i, j, d: embed(Smat(d, lam), i, j, 4)
    # rapidity differences: slot1 - slot3 etc.
    d = {(0,2): th, (0,3): th + 1j*u, (1,2): th - 1j*u, (1,3): th}
    S13, S14, S23, S24 = S(0,2,d[(0,2)]), S(0,3,d[(0,3)]), S(1,2,d[(1,2)]), S(1,3,d[(1,3)])
    x = np.kron(v, v)
    res = {}
    for name, A in [('(S23 S13)(S24 S14)', S23@S13@S24@S14), ('(S13 S14)(S23 S24)', S13@S14@S23@S24),
                    ('S14 S13 S24 S23', S14@S13@S24@S23)]:
        y = A @ x; c = np.vdot(x, y)/np.vdot(x, x)
        res[name] = (c, np.linalg.norm(y - c*x))
    return res
for th in (0.37+0.11j, -0.6+0.3j, 1.2-0.05j):
    print(f'theta = {th}:  f_(xi/pi) = {fx(xi/np.pi, th):.10f}')
    for k, (c, r) in S_BB(th).items():
        print(f'    {k}: S_11 = {c:.10f}  (non-invariance {r:.1e})')
