"""Independent B2 completion from roots and the two closed-form vacua.

Use a matrix logarithmic derivative, not exponentially growing Jost columns.
The check covers sampled subthreshold frequencies, not an interval proof.
"""
import numpy as np
from scipy.integrate import solve_ivp


def data(n):
    e = np.eye(n)
    roots = np.array([-e[0]-e[1]] + [e[j]-e[j+1] for j in range(n-1)]
                     + [e[-2]+e[-1]])
    marks = np.array([1, 1]+[2]*(n-3)+[1, 1])
    if n == 4:
        c = (3*np.sqrt(3)-5)/(3*np.sqrt(3)+5)
        def field(x):
            z = c*np.exp(np.sqrt(6)*x)
            u = np.sqrt(3)*(1-z)/(1+z)
            return np.log(2/(u*u-1))*(e[1]-e[2])
    else:
        c = np.tanh(np.log(4-2*np.sqrt(2))/4)
        def field(x):
            return 4*np.arctanh(c*np.exp(2*x))*(e[1]-e[3])
    def mass(x):
        return (roots.T*(marks*np.exp(roots@field(x))))@roots
    R = (roots.T*(np.sqrt(marks)*np.exp(roots@field(0)/2)))@roots/2
    return mass, R, (roots.T*marks)@roots


def impedance(n, energy, L=20., tolerance=2e-10):
    mass, R, H0 = data(n)
    eig, U = np.linalg.eigh(H0)
    Y0 = (U*np.sqrt(np.maximum(eig-energy, 0)))@U.T
    def rhs(x, y):
        Y = y.reshape(n, n)
        return (mass(x)-energy*np.eye(n)-Y@Y).ravel()
    sol = solve_ivp(rhs, (-L, 0), Y0.ravel(), method='DOP853',
                    first_step=.01, max_step=.15, rtol=tolerance,
                    atol=tolerance/100)
    assert sol.success, sol.message
    Y = sol.y[:, -1].reshape(n, n)
    assert np.linalg.norm(Y-Y.T) < 1e-10
    return np.linalg.eigvalsh((Y+Y.T)/2+R)


def main():
    for n in (4, 5):
        edge = np.linalg.eigvalsh(data(n)[2])[0]
        energies = np.unique(np.r_[np.linspace(-4, edge-.001, 101),
                                   edge-10.**(-np.arange(2, 11)), 1.])
        values = np.array([impedance(n, float(v)) for v in energies])
        i = np.argmin(values[:, 0])
        assert values[i, 0] > .1
        print(f'PASS d{n}: {len(energies)} subthreshold points, '
              f'm_min^2={edge:.12f}; min eig(Y+R)={values[i,0]:.12f} '
              f'at Omega^2={energies[i]:.12f}', flush=True)
        delta = 0
        for en in (1., edge-1e-3, edge-1e-7, edge-1e-10):
            a = impedance(n, en)
            b = impedance(n, en, L=30, tolerance=2e-12)
            delta = max(delta, np.max(np.abs(a-b)))
        assert delta < 1e-8
        print(f'PASS d{n}: domain/tolerance refinement, max eigenvalue change '
              f'{delta:.3e}; at Omega^2=1 min eig={impedance(n,1.)[0]:.12f}',
              flush=True)


if __name__ == '__main__':
    main()
