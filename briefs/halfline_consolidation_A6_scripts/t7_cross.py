"""Crossing-unitarity of the d_n K (exact polynomial normalization from t5_symK): k(theta) = Psi Rcheck_(Vbar Vbar)(2theta)
k(-theta), k(theta) = vec((P^-1 K(i pi/2 - theta))^T), P: Vbar(y) -> V(y). Same construction as an_nondiagonal_scalar.py
check 6, with K: V -> V turned into V -> Vbar by P^-1."""
import sys, pickle
import numpy as np
import sympy as sp
from qk import dn_vector, q_phys, chevalley, rcheck_w, intertwiners, nodes_of, K_raw

w_, y_ = sp.symbols('w y')
n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
D = pickle.load(open(f'symK_d{n}.pkl', 'rb'))
Kf = {name: sp.lambdify((w_, y_), sp.sympify(eval(D[name], sp.__dict__)), 'numpy') for name in D}
N = 2*n


def Kmat(name, om, th):
    w = np.exp(-0.5j*np.pi*om)
    y = np.exp(om*th)
    return np.array(Kf[name](w, y), complex)


def check_intertwiner(name, om, th, s):
    K = Kmat(name, om, th)
    Kn = K_raw(n, om, np.exp(om*th), s)
    r = (K/K[0, 0]) - (Kn/Kn[0, 0])
    return np.abs(r).max()


def Piso(om):
    q = q_phys(om)
    V1 = dn_vector(n, q, 1.0)
    yy = 0.83 + 0.21j
    P = intertwiners([(chevalley(V1, yy)[(k, j)], dn_vector(n, q, yy)[(k, j)]) for j in nodes_of(V1) for k in 'efk'])[0]
    assert len(P) == 1
    P = P[0]
    yy2 = 1.3 - 0.4j
    P2 = intertwiners([(chevalley(V1, yy2)[(k, j)], dn_vector(n, q, yy2)[(k, j)]) for j in nodes_of(V1) for k in 'efk'])[0][0]
    dev = np.abs(P/P[0, -1] - P2/P2[0, -1]).max()
    return P, dev


if __name__ == '__main__':
    np.set_printoptions(linewidth=200, precision=5)
    for om in (0.37, 1.29):
        q = q_phys(om)
        T = (n - 1)*om
        P, devP = Piso(om)
        Pi = np.linalg.inv(P)
        print(f'om={om}: P antidiagonal? {np.abs(P - np.fliplr(np.diag(np.diag(np.fliplr(P))))).max():.1e}, '
              f'y-independent {devP:.1e}')
        for name, s in (('0', 0), ('s+', 1), ('s-', -1)):
            print(f'  {name}: |K_sym - K_book| = {check_intertwiner(name, om, 0.31 + 0.2j, s):.1e}')
            for th in (0.31, 0.83 + 0.1j, -0.4 + 0.25j):
                km = (Pi @ Kmat(name, om, 1j*np.pi/2 - th)).T.reshape(N*N)
                kp = (Pi @ Kmat(name, om, 1j*np.pi/2 + th)).T.reshape(N*N)
                R = rcheck_w(chevalley(dn_vector(n, q, 1.0), np.exp(2*om*th)), chevalley(dn_vector(n, q, 1.0), 1.0))
                v = R @ kp
                f = np.vdot(v, km)/np.vdot(v, v)
                res = np.linalg.norm(km - f*v)/np.linalg.norm(km)
                mu = -1j*T*th/np.pi
                print(f'     th={th}: proportional [{res:.1e}], Psi = {f:.10f}, mu = {mu:.4f}')
