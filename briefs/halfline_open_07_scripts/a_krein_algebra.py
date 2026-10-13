"""Problem 7(a): Krein structure for the a_n soliton S- and K-matrices, sine-Gordon and d_n checks."""
import numpy as np
from boundary_common import an_vector, q_phys, k_intertwiner, rcheck, sl2_spin, dn_vector, intertwiners, coideal_gen
from p07_common import Cmat, nu, Rc, k_book, cK, eps_star, A_min
rng = np.random.default_rng(3)

def prop(M, R):
    c = np.vdot(R.ravel(), M.ravel())/np.vdot(R.ravel(), R.ravel())
    return np.linalg.norm(M - c*R)/np.linalg.norm(M), c

# 1. pi(a)^dag = C pi(tau a) C for real y, tau(e_j)=e_nu(j), tau(f_j)=f_nu(j), tau(k_j)=k_nu(j)
for n in (2, 3, 4):
    N = n + 1; C = Cmat(N)
    for om in (0.37, 1.43):
        q = q_phys(om); y = np.exp(0.7)
        err = 0
        for conj in (False, True):
            P = an_vector(n, q, y, conj)
            for j in range(N):
                for k in 'efk':
                    err = max(err, np.abs(P[(k, j)].conj().T - C @ P[(k, nu(j, N))] @ C).max())
        print(f'[1] a_{n} om={om}: pi(a)^dag = C pi(tau a) C on V and Vbar, real y: err {err:.1e}')

# 2. R-checks: R^dag (C x C) R = c (C x C), all four pairs, real rapidities
for n in (2, 3):
    N = n + 1; CC = np.kron(Cmat(N), Cmat(N))
    for om in (0.37, 1.43, 2.6):
        out = []
        for a, b in (('V', 'V'), ('V', 'W'), ('W', 'V'), ('W', 'W')):
            for th in (0.3, -1.1, 2.0):
                R = Rc(n, a, b, np.exp(om*th), om)
                r, c = prop(R.conj().T @ CC @ R, CC)
                out.append((a+b, th, r, c))
        rmax = max(o[2] for o in out); cdev = max(abs(o[3] - 1) for o in out)
        print(f'[2] a_{n} om={om}: R^dag CC R = c CC for VV,VW,WV,WW: prop err {rmax:.1e}, max|c-1| {cdev:.1e}')

# 3. K-matrices: closed form = intertwiner; K^dag C K = cK C; K_Vbar = K^T intertwines; Hermitian analyticity
for n in (2, 3, 4):
    N = n + 1; C = Cmat(N)
    for om in (0.37, 1.43, 2.6):
        q = q_phys(om); es = eps_star(om)
        res = []
        for s in (1, -1):
            for th in (0.4, -0.9, 1.7):
                y = np.exp(om*th)
                e = es*np.ones(N); e[0] *= s
                null, _ = k_intertwiner(an_vector(n, q, y), an_vector(n, q, 1/y, True), [1/q]*N, e)
                K = k_book(n, th, s, om)
                pr = prop(null[0], K)[0]
                rK, c = prop(K.conj().T @ C @ K, C)
                # Vbar -> V: K^T intertwines
                nullb, _ = k_intertwiner(an_vector(n, q, y, True), an_vector(n, q, 1/y), [1/q]*N, e)
                prb = prop(nullb[0], K.T)[0]
                herm = np.abs(K.conj().T - s*k_book(n, -th, s, om)).max()/np.abs(K).max()
                res.append((pr, rK, abs(c/cK(n, th, s, om) - 1), prb, herm))
        r = np.max(np.array(res), axis=0)
        print(f'[3] a_{n} om={om}: K=intertwiner {r[0]:.0e}; K^dag C K prop C {r[1]:.0e}; c=cK(formula) {r[2]:.0e}; '
              f'K_Vbar=K^T {r[3]:.0e}; K(th)^dag = s K(-th) {r[4]:.0e}')
    # unitarity K(-th)^T K(th) = U = s cK, and |A_min|^2 cK = 1, A_min(th)A_min(-th) = 1/U
    om = 0.83
    for s in (1, -1):
        th = 0.61
        U = (k_book(n, -th, s, om).T @ k_book(n, th, s, om))
        print(f'    a_{n} s={s}: K(-th)^T K(th) = s cK 1: {np.abs(U - s*cK(n, th, s, om)*np.eye(N)).max():.1e}; '
              f'|A_min|^2 cK - 1 = {abs(abs(A_min(n, th, s, om))**2*cK(n, th, s, om) - 1):.1e}; '
              f'A_min(th)A_min(-th) s cK - 1 = {abs(A_min(n, th, s, om)*A_min(n, -th, s, om)*s*cK(n, th, s, om) - 1):.1e}')

# 4. coideal stability: tau(b_j(eps)) = b_nu(j)(eps') with eps'_nu(j) = conj(eps_j); test on V(y), random eps
n = 2; N = 3; om = 0.37; q = q_phys(om); C = Cmat(N); y = 1.7
eps = rng.normal(size=N) + 1j*rng.normal(size=N)
V = an_vector(n, q, y)
err = 0
for j in range(N):
    b = coideal_gen(V, j, 1/q, eps[j])
    # tau is antilinear anti-automorphism; on V: pi(tau x) = C pi(x)^dag C
    lhs = C @ b.conj().T @ C
    rhs = coideal_gen(V, nu(j, N), 1/q, np.conj(eps[j]))
    err = max(err, np.abs(lhs - rhs).max())
print(f'[4] tau(b_j(eps)) = b_nu(j)(conj eps_j), c_j = q^-1 preserved: err {err:.1e}')

# 5. sine-Gordon: K for general eps; C-Krein (sigma^x) iff eps real; positive unitarity iff eps0 = conj eps1
def ksg(y, q, e0, e1):
    return np.array([[(q - 1/q)*(e1*y + e0/y), y**2 - y**-2], [y**2 - y**-2, (q - 1/q)*(e0*y + e1/y)]])/(y**2 - y**-2)
om = 0.37; q = q_phys(om); sx = np.array([[0, 1], [1, 0]], complex); y = np.exp(om*0.8)
for name, (e0, e1) in {'real eps': (0.7, -1.3), 'eps0 = conj eps1': (0.4 + 0.9j, 0.4 - 0.9j),
                       'eps0=eps1 real': (0.6, 0.6), 'generic': (0.3 + 0.2j, -0.5 + 0.7j)}.items():
    K = ksg(y, q, e0, e1)
    rC, cC = prop(K.conj().T @ sx @ K, sx); r1, c1 = prop(K.conj().T @ K, np.eye(2))
    print(f'[5] SG {name:18s}: K^dag sx K prop sx: {rC:.1e} (c={cC.real:+.3f}{cC.imag:+.1e}i); K^dag K prop 1: {r1:.1e}')
