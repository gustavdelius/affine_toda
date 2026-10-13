# crossing factor c_u for the d_n vector R-matrix (principal gradation): (R(x)^-1)^{t1} = c_u (C^-1 x 1) R(x x_c) (C x 1)
from dn_common import *
np.set_printoptions(linewidth=200, precision=5, suppress=True)
def pt1(M, N):
    return M.reshape(N, N, N, N).transpose(2, 1, 0, 3).reshape(N*N, N*N)
for n in (4, 5):
    N = 2*n; h = 2*n - 2
    for om in (0.37, 1.31):
        q = q_phys(om); T = h*om/2
        yc = -1/q
        y0 = 0.7 + 0.2j
        # C: V(y0)^* -> V(y0 yc)
        Vd = dual(dn_vector(n, q, y0)); Vc = dn_vector(n, q, y0*yc)
        pairs = [(Vd[(k, j)], Vc[(k, j)]) for j in range(n + 1) for k in 'efk']
        null, _ = intertwiners(pairs); assert len(null) == 1
        C = null[0]
        lab = basis_labels(n)
        # check C is y-independent
        Vd2 = dual(dn_vector(n, q, 1.3)); Vc2 = dn_vector(n, q, 1.3*yc)
        null2, _ = intertwiners([(Vd2[(k, j)], Vc2[(k, j)]) for j in range(n + 1) for k in 'efk'])
        print('n', n, 'om', om, 'C y-indep resid', prop_residual(null2[0], C)[0])
        if om == 0.37:
            nz = np.nonzero(np.abs(C) > 1e-9)
            print('   C entries (row label, col label, val):', [(lab[r], lab[c], np.round(C[r, c]/C[nz[0][0], nz[1][0]], 4)) for r, c in zip(*nz)])
        Ci = np.linalg.inv(C)
        for mu in (0.13 + 0.21j, -0.31 + 0.05j):
            y = np.exp(2j*np.pi*mu/h)   # x = y^h = e^{2 pi i mu}
            Rx = np.eye(N*N)[:, :]  # placeholder
            Rc = fast_rcheck(dn_vector(n, q, y), dn_vector(n, q, 1.0))
            P = np.eye(N*N).reshape(N, N, N, N).transpose(1, 0, 2, 3).reshape(N*N, N*N)
            R = P @ Rc
            Rcc = fast_rcheck(dn_vector(n, q, y*yc), dn_vector(n, q, 1.0))
            Rxc = P @ Rcc
            L = pt1(np.linalg.inv(R), N)
            Rr = kron(Ci, np.eye(N)) @ Rxc @ kron(C, np.eye(N))
            res, cu = prop_residual(L, Rr)
            # candidate closed forms in mu
            s = lambda z: np.sin(np.pi*z)
            cand = {'A': s(mu - om)*s(mu - T)/(s(mu + om)*s(mu + T)) if False else None}
            print('   mu', mu, 'resid', f'{res:.1e}', 'c_u', cu,
                  ' ratios:', cu*s(mu + om)*s(mu + T)/(s(mu)*s(mu + T - om)) )
