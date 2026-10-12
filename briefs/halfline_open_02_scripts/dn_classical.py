"""Real-coupling d_n^(1) Toda on x <= 0 (beta = m = 1): V = sum_j n_j (e^{alpha_j.phi} - 1),
B = 2 sum_j C_j sqrt(n_j) (e^{alpha_j.phi/2} - 1), BC phi'(0) = -sum_j C_j sqrt(n_j) alpha_j e^{alpha_j.phi(0)/2}."""
import numpy as np
from scipy.integrate import solve_ivp


def droots(n):
    E = np.eye(n)
    al = [-E[0] - E[1]] + [E[i] - E[i + 1] for i in range(n - 1)] + [E[n - 2] + E[n - 1]]
    marks = [1, 1] + [2]*(n - 3) + [1, 1]
    return np.array(al), np.array(marks, float)


def mass_matrix(n):
    al, nm = droots(n)
    return (al.T*nm) @ al


def masses(n):
    ev, U = np.linalg.eigh(mass_matrix(n))
    return np.sqrt(np.clip(ev, 0, None)), U


if __name__ == '__main__':
    for n in (4, 5, 6):
        h = 2*n - 2
        m, U = masses(n)
        print(n, 'masses', np.round(m, 6), ' 2 sin(a pi/h)*sqrt2:', np.round([2*np.sqrt(2)*np.sin(a*np.pi/h) for a in range(1, n - 1)], 6),
              ' sqrt2:', np.sqrt(2))


class Static:
    """Static solutions decaying at -inf, parametrized by asymptotic amplitudes c (phi ~ sum_a c_a u_a e^{m_a x})."""
    def __init__(self, n, C, L=None):
        self.n = n
        self.al, self.nm = droots(n)
        self.C = np.asarray(C, float)
        self.m, self.U = masses(n)
        self.L = L if L is not None else 26.0/self.m.min()

    def rhs(self, x, Y):
        n = self.n
        phi, dphi = Y[:n], Y[n:2*n]
        e = np.exp(self.al @ phi)
        return np.concatenate([dphi, (self.al.T*self.nm) @ e, [dphi @ dphi]])

    def solve(self, c, dense=False):
        x0 = -self.L
        ex = np.exp(self.m*x0)
        phi0 = self.U @ (c*ex); dphi0 = self.U @ (c*self.m*ex)
        tail = np.sum((c*ex)**2*self.m/2)            # int_{-inf}^{x0} |phi'|^2 (leading order)
        def blow(x, Y):
            return 60 - np.max(np.abs(Y[:self.n]))
        blow.terminal = True
        s = solve_ivp(self.rhs, [x0, 0.0], np.concatenate([phi0, dphi0, [tail]]), method='DOP853',
                      rtol=1e-12, atol=1e-14, events=blow, dense_output=dense)
        if s.status != 0:
            return None
        return s

    def bpot(self, phi):
        return 2*np.sum(self.C*np.sqrt(self.nm)*(np.exp(self.al @ phi/2) - 1))

    def bc_res(self, phi, dphi):
        return dphi + (self.al.T*(self.C*np.sqrt(self.nm))) @ np.exp(self.al @ phi/2)

    def F(self, c):
        s = self.solve(c)
        if s is None:
            return np.inf
        Y = s.y[:, -1]
        return Y[-1] + self.bpot(Y[:self.n])


def linear_K(n, C, sol, Om, X=None):
    """Linearized reflection about the static solution sol (dense output), x in [-L, 0]: y'' = (M2(x) - Om^2) y,
    y'(0) = -H y(0), H = 1/2 sum_j C_j sqrt(n_j) alpha_j alpha_j^T e^{alpha_j.phi(0)/2}.
    Returns K (open x open) in the mass eigenbasis U (columns), and the open channel indices."""
    al, nm = droots(n)
    C = np.asarray(C, float)
    m, U = masses(n)
    x0 = sol.t[0]
    phi = lambda x: sol.sol(x)[:n] if x >= x0 else np.zeros(n)
    def M2(x):
        e = np.exp(al @ phi(x)); return (al.T*(nm*e)) @ al
    phi0 = sol.sol(0.0)[:n]
    H = 0.5*(al.T*(C*np.sqrt(nm)*np.exp(al @ phi0/2))) @ al
    xs = x0
    open_ = [a for a in range(n) if Om > m[a] + 1e-12]
    closed = [a for a in range(n) if a not in open_]
    def init(a, kind):
        if kind == 'in':
            p = np.sqrt(Om**2 - m[a]**2); f = np.exp(1j*p*xs); fp = 1j*p*f
        elif kind == 'out':
            p = np.sqrt(Om**2 - m[a]**2); f = np.exp(-1j*p*xs); fp = -1j*p*f
        else:
            k = np.sqrt(m[a]**2 - Om**2); f = np.exp(k*xs); fp = k*f
        return U[:, a]*f, U[:, a]*fp
    def rhs(x, Y):
        y = Y[:n]; yp = Y[n:]
        return np.concatenate([yp, (M2(x) - Om**2*np.eye(n)) @ y])
    sols = []
    for kind, lst in (('in', open_), ('out', open_), ('cl', closed)):
        for a in lst:
            y0, yp0 = init(a, kind)
            s = solve_ivp(rhs, [xs, 0.0], np.concatenate([y0, yp0]).astype(complex), method='DOP853', rtol=1e-11, atol=1e-13)
            Y = s.y[:, -1]
            sols.append((kind, a, Y[n:] + H @ Y[:n]))
    nin = len(open_)
    Rin = np.array([r for k, a, r in sols if k == 'in']).T
    Rrest = np.array([r for k, a, r in sols if k != 'in']).T
    solv = np.linalg.solve(Rrest, -Rin)
    return solv[:nin, :], open_
