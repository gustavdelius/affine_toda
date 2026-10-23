"""A6 task 1: the C = +1 vacuum of d_4, d_5 (closed form) and the linearized reflection matrix about it.
Own code. Units beta = m = 1, book potentials (eq-boundary-potential): psi = beta phi,
V = sum_j n_j (e^{a_j.psi} - 1), B = 2 sum_j C_j sqrt(n_j) (e^{a_j.psi/2} - 1),
psi'' = sum_j n_j a_j e^{a_j.psi}, psi'(0) = -sum_j C_j sqrt(n_j) a_j e^{a_j.psi(0)/2}.
Roots as boundary_common.dn_vector: a_i = e_i - e_{i+1}, a_n = e_{n-1} + e_n, a_0 = -e_1 - e_2.
"""
import sys
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp

sys.path.insert(0, '/home/gustav/Git/affine_toda/code')


def droots(n):
    E = np.eye(n)
    A = [-E[0] - E[1]] + [E[i] - E[i + 1] for i in range(n - 1)] + [E[n - 2] + E[n - 1]]
    nj = [1, 1] + [2]*(n - 3) + [1, 1]
    return np.array(A), np.array(nj, float)


class Kink:
    """Invariant static solution psi = f(x) v for C = +1."""
    def __init__(self, n):
        self.n = n
        E = np.eye(n)
        if n == 4:
            self.v = E[1] - E[2]
            r3 = np.sqrt(3.0)
            self.E0 = (3*r3 - 5)/(3*r3 + 5)
            self.rate = np.sqrt(6.0)
        elif n == 5:
            self.v = E[1] - E[3]
            f0 = np.log(4 - 2*np.sqrt(2))
            self.E0 = np.tanh(f0/4)       # tanh(f/4) = E0 e^{2x}
            self.rate = 2.0
        self.A, self.nj = droots(n)
        self.c = self.A @ self.v

    def f(self, x):
        E = self.E0*np.exp(self.rate*x)
        if self.n == 4:
            s = np.sqrt(3.0)*(1 - E)/(1 + E)
            return np.log(2/(s**2 - 1))
        return 4*np.arctanh(E)

    def fp(self, x):
        f = self.f(x)
        if self.n == 4:
            return np.sqrt(2*(np.exp(2*f) + 2*np.exp(-f) - 3))
        return 4*np.sinh(f/2)

    # mp versions for precise checks
    def f_mp(self, x):
        if self.n == 4:
            r3 = mp.sqrt(3)
            E = (3*r3 - 5)/(3*r3 + 5)*mp.exp(mp.sqrt(6)*x)
            s = r3*(1 - E)/(1 + E)
            return mp.log(2/(s**2 - 1))
        f0 = mp.log(4 - 2*mp.sqrt(2))
        return 4*mp.atanh(mp.tanh(f0/4)*mp.exp(2*x))


def M2(K, x):
    e = np.exp(K.c*K.f(x))
    return (K.A.T*(K.nj*e)) @ K.A


def Rbc(K, C=None):
    C = np.ones(K.n + 1) if C is None else C
    e = np.exp(K.c*K.f(0.0)/2)
    return 0.5*(K.A.T*(C*np.sqrt(K.nj)*e)) @ K.A


def kmatrix(K, Om, L=None, rtol=1e-13):
    """Linearized reflection matrix at frequency Om, mass eigenbasis of M0^2. Returns (masses, U, open idx, Kmat)."""
    n = K.n
    M0 = (K.A.T*K.nj) @ K.A
    ev, U = np.linalg.eigh(M0)
    ma = np.sqrt(ev)
    L = L if L is not None else 34/K.rate
    op = [a for a in range(n) if ma[a] < Om]
    cl = [a for a in range(n) if ma[a] >= Om]
    sols = []
    def rhs(x, y):
        return np.concatenate([y[n:], (M2(K, x) - Om**2*np.eye(n)) @ y[:n]])
    x0 = -L
    starts = []
    for a in op:
        p = np.sqrt(Om**2 - ma[a]**2)
        for sg in (1, -1):
            starts.append(('o', a, sg, np.concatenate([U[:, a]*np.exp(1j*sg*p*x0), 1j*sg*p*U[:, a]*np.exp(1j*sg*p*x0)])))
    for a in cl:
        k = np.sqrt(ma[a]**2 - Om**2)
        starts.append(('c', a, 0, np.exp(-k*L)*np.concatenate([U[:, a], k*U[:, a]]).astype(complex)))
    Y0 = np.array([s[3] for s in starts]).T
    # integrate all solutions together (matrix ODE)
    m = Y0.shape[1]
    def rhsM(x, y):
        Y = y.reshape(2*n, m)
        return np.vstack([Y[n:], (M2(K, x) - Om**2*np.eye(n)) @ Y[:n]]).ravel()
    sol = solve_ivp(rhsM, [x0, 0.0], Y0.ravel(), method='DOP853', rtol=rtol, atol=1e-14)
    Y = sol.y[:, -1].reshape(2*n, m)
    R = Rbc(K)
    Bm = Y[n:] + R @ Y[:n]                   # BC residual of each solution
    # columns: in (sg=+1), out (sg=-1), closed
    iin = [i for i, s in enumerate(starts) if s[0] == 'o' and s[2] == 1]
    iout = [i for i, s in enumerate(starts) if s[0] == 'o' and s[2] == -1]
    icl = [i for i, s in enumerate(starts) if s[0] == 'c']
    Mu = Bm[:, iout + icl]
    Kmat = np.zeros((len(op), len(op)), complex)
    for b, ib in enumerate(iin):
        coef = np.linalg.solve(Mu, -Bm[:, ib])
        Kmat[:, b] = coef[:len(op)]
    return ma, U, op, Kmat


def bstate_det(K, Om, L=None):
    """det of the BC matrix of the decaying solutions for Om below all thresholds (bound state where it vanishes)."""
    n = K.n
    M0 = (K.A.T*K.nj) @ K.A
    ev, U = np.linalg.eigh(M0)
    ma = np.sqrt(ev)
    L = L if L is not None else 34/K.rate
    kap = np.sqrt(ma**2 - Om**2)
    Y0 = np.vstack([U, U*kap])*np.exp(-kap*L)
    def rhsM(x, y):
        Y = y.reshape(2*n, n)
        return np.vstack([Y[n:], (M2(K, x) - Om**2*np.eye(n)) @ Y[:n]]).ravel()
    sol = solve_ivp(rhsM, [-L, 0.0], Y0.ravel(), method='DOP853', rtol=1e-13, atol=1e-14)
    Y = sol.y[:, -1].reshape(2*n, n)
    Y = Y/np.linalg.norm(Y[:n], axis=0)
    return np.linalg.det(Y[n:] + Rbc(K) @ Y[:n])
