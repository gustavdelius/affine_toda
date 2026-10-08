import numpy as np
from dtw import rep, weights
from rsolve import intertwiner
q, n = 0.7, 3
W = weights(n); N = 2*n+2; nodes = range(n+1)
def Rhat(x):
    Rc, null, _, _ = intertwiner(rep(n, q, x), rep(n, q, 1.0), nodes, W, W)
    top = 0*N + 0   # w1 (x) w1
    return Rc/Rc[top, top], null
# 1) Yang-Baxter equation (braid form) on V(z1)(x)V(z2)(x)V(z3)
z1, z2, z3 = 1.3+0.4j, 0.6-0.2j, 0.9+0.7j
I = np.eye(N)
R12 = lambda z: np.kron(Rhat(z)[0], I); R23 = lambda z: np.kron(I, Rhat(z)[0])
lhs = R12(z2/z3) @ R23(z1/z3) @ R12(z1/z2); rhs = R23(z1/z2) @ R12(z1/z3) @ R23(z2/z3)
print(f"Yang-Baxter: max |LHS-RHS| / max|LHS| = {np.abs(lhs-rhs).max()/np.abs(lhs).max():.1e}")
# 2) eigenvalue structure at generic x
R, _ = Rhat(1.9+0.3j); ev = np.linalg.eigvals(R)
vals = []
for v in ev:
    for c in vals:
        if abs(c[0]-v) < 1e-6*max(1, abs(v)): c[1] += 1; break
    else: vals.append([v, 1])
print("distinct eigenvalues and multiplicities:", sorted([m for _, m in vals], reverse=True),
      "-> U_q(b_3) content 27+21+7+7+1+1: the 7 and the 1 each occur twice (not multiplicity-free)")
# 3) special points x = +-q^l : classify singular values as delta -> 0
def image_dim(xs):
    dims = []
    for d in (1e-5, 1e-7):
        R, _ = Rhat(xs*(1+d)); R = R/np.linalg.norm(R); s = np.linalg.svd(R, compute_uv=False)
        dims.append(np.sum(s > 1e-3*s[0]))
    return dims
print("special points (image dimension of the normalized R-matrix as x -> x*, at two offsets):")
for l in np.arange(-7, 7.5, 0.5):
    for sgn in (+1, -1):
        dims = image_dim(sgn*q**l)
        if dims[1] < N*N:
            print(f"   x* = {'+' if sgn>0 else '-'}q^{l:+.1f}: image dimension {dims[1]} (offsets 1e-5, 1e-7: {dims})")
