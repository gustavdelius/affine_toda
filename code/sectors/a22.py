"""U_q(a_2^(2)) on C^3 (Izergin-Korepin); weights m=1,0,-1 ; alpha_1 short (q_1=q), alpha_0=-2 alpha_1 long (q_0=q^4)"""
import numpy as np, sys
from rsolve import intertwiner
def rep(x, q):
    qn = lambda k, Q: (Q**k - Q**-k)/(Q - 1/Q)
    E = lambda a, b: (lambda M: (M.__setitem__((a, b), 1), M)[1])(np.zeros((3, 3), complex))
    # basis index 0:m=1, 1:m=0, 2:m=-1 ; spin-1 of U_q(sl2) with q_1 = q, k_1 = q^{2m}
    c = np.sqrt(qn(2, q) + 0j)
    e1 = c*(E(0, 1) + E(1, 2)); f1 = c*(E(1, 0) + E(2, 1)); k1 = np.diag([q**2, 1, q**-2]).astype(complex)
    e0 = x*E(2, 0); f0 = E(0, 2)/x; k0 = np.diag([q**-4, 1, q**4]).astype(complex)
    return {1: e1, 0: e0}, {1: f1, 0: f0}, {1: k1, 0: k0}, {1: q, 0: q**4}
def check(q, x=0.7+0.2j):
    e, f, k, qi = rep(x, q); err = 0
    for i in (0, 1):
        for j in (0, 1):
            c = e[i]@f[j] - f[j]@e[i]
            tgt = (k[i] - np.linalg.inv(k[i]))/(qi[i] - 1/qi[i]) if i == j else 0
            err = max(err, np.abs(c - tgt).max())
    return err
W = np.array([[1.], [0.], [-1.]])
if __name__ == "__main__":
    P = np.zeros((9, 9))
    for i in range(3):
        for j in range(3): P[j*3+i, i*3+j] = 1
    for om in [0.13, 0.37, 0.61, 0.87, 1.23, 1.71]:
        q = np.exp(-1j*np.pi*om)
        print("omega", om, "rel err", check(q))
        devp, devm = [], []
        for lx in np.linspace(-12, 12, 49):
            for sgn, lst in ((1, devp), (-1, devm)):
                x = sgn*np.exp(lx)
                R = intertwiner(rep(x, q), rep(1.0, q), (0, 1), W, W)[0]; R = R/R[0, 0]
                e = np.linalg.eigvals(P @ R); lst.append(np.abs(np.abs(e) - 1).max())
        print("   x>0 max dev", max(devp), "   x<0 max dev", max(devm))
