# boundary-state vector k = (C x 1) vec(J^T) of the eps = 0 d_n reflection lies in L(2w1): R-check(x) k = k, so Phi = 1
from dn_common import *
for n in (4, 5):
    for om in (0.37, 1.31):
        q = q_phys(om); N = 2*n; yc = -1/q
        J = Kdn(n, q, 0.9 + 0.2j, np.zeros(n + 1)); J = J/J[N - 1, 0]
        Vd = dual(dn_vector(n, q, 0.7)); Vc = dn_vector(n, q, 0.7*yc)
        C = intertwiners([(Vd[(k, j)], Vc[(k, j)]) for j in range(n + 1) for k in 'efk'])[0][0]
        C = C/np.abs(C).max()
        for kv in (np.kron(C, np.eye(N)) @ J.T.reshape(N*N), np.kron(C, np.eye(N)) @ J.reshape(N*N)):
            worst = 0
            for th in (0.31 + 0.2j, -0.7 + 0.1j):
                R2 = fast_rcheck(dn_vector(n, q, np.exp(2*om*th)), dn_vector(n, q, 1.0))
                worst = max(worst, np.linalg.norm(R2 @ kv - kv)/np.linalg.norm(kv))
            print(n, om, 'R(2th) k = k residual', f'{worst:.1e}')
