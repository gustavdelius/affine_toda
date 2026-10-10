import numpy as np
from rsolve import intertwiner
omega = 2.37; q = np.exp(-1j*np.pi*omega); D = 2
sp, sm = np.array([[0,1],[0,0]],complex), np.array([[0,0],[1,0]],complex)
W1 = np.array([[0.5], [-0.5]])
def rep(x):   # nodes 0 and 1, both short (q_i = q): U_q(a_1^(1)) doublet = sine-Gordon soliton multiplet
    e = {1: sp.copy(), 0: x*sm}; f = {1: sm.copy(), 0: sp/x}
    k = {1: np.diag(q**(2*W1[:, 0])).astype(complex), 0: np.diag(q**(-2*W1[:, 0])).astype(complex)}
    return e, f, k, {0: q, 1: q}
P = np.zeros((4, 4))
for i in range(2):
    for j in range(2): P[j*2+i, i*2+j] = 1
def Rh(x):
    Rc = intertwiner(rep(x), rep(1.0), range(2), W1, W1)[0]; return Rc/Rc[0, 0]
def dual(x):
    e, f, k, _ = rep(x); out = []
    for i in range(2):
        ki = np.linalg.inv(k[i]); out += [(-ki@e[i]).T, (-f[i]@k[i]).T, ki.T]
    return out
def plain(x):
    e, f, k, _ = rep(x); out = []
    for i in range(2): out += [e[i], f[i], k[i]]
    return out
x0 = 0.83+0.29j
for ph in (1, -1):
    for l in range(-4, 5):
        tt = ph*q**l
        M = np.vstack([np.kron(np.eye(D), A.T) - np.kron(B, np.eye(D)) for A, B in zip(dual(x0), plain(x0*tt))])
        s = np.linalg.svd(M, compute_uv=False)
        if s[-1] < 1e-9*s[0]: C = np.linalg.svd(M)[2][-1].conj().reshape(D, D); tc = tt; print(f"n=1: crossing point x(i pi) = {'+' if ph>0 else '-'}q^{l}")
Ci = np.linalg.inv(C)
def pt1(M): return M.reshape(D, D, D, D).transpose(2, 1, 0, 3).reshape(4, 4)
def cu(x):
    R = P @ Rh(x); lhs = pt1(np.linalg.inv(R)); rhs = np.kron(Ci, np.eye(D)) @ (P @ Rh(x*tc)) @ np.kron(C, np.eye(D))
    return np.vdot(rhs.ravel(), lhs.ravel())/np.vdot(rhs.ravel(), rhs.ravel())
# eigenvalue ratio of the singlet relative to the triplet, to identify the factor <l>
xs = [0.7+0.3j, 1.6-0.2j]
for x in xs:
    ev = np.linalg.eigvals(Rh(x)); print("n=1 eigenvalues of R(x):", np.round(ev, 6), " <2>_+ =", np.round((x - q**2)/(1 - x*q**2), 6))
T = omega; A = [omega]
c_of = lambda t: np.prod([np.sin(np.pi*(t - a)) for a in A])
for t in (0.37+0.21j, 1.3-0.4j):
    x = np.exp(2j*np.pi*t); r = (c_of(T - t)/c_of(t))/cu(1/x)
    print(f"n=1 (sine-Gordon): [c(T-t)/c(t)] / c_u(1/x) = {r.real:+.6f}{r.imag:+.6f}i")
