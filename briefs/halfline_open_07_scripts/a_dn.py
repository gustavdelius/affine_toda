"""Problem 7(a) for the d_n vector solitons (n = 4): the antilinear anti-automorphism tau_nu (e_j -> e_nu(j), ...),
a Hermitian eta with pi(a)^dag = eta pi(tau_nu a) eta^-1 at real y, Krein property of the R-check (weight-filtered
solver copied from code/reflection_matrices.py) and of the K-matrices of prp-dn-k."""
import numpy as np, scipy.sparse as sps
from boundary_common import dn_vector, q_phys, k_intertwiner, intertwiners, cop, nodes_of
n = 4; N = 2*n; om = 0.37; q = q_phys(om); es = 1/(2*np.cos(np.pi*om/2))
def rcheck_w(A, B):
    dA, dB = A['dim'], B['dim']; js = nodes_of(A)
    win = list(zip(*[np.round(np.diag(cop(A, B, 'k', j)), 9) for j in js]))
    wout = list(zip(*[np.round(np.diag(cop(B, A, 'k', j)), 9) for j in js]))
    cols = [(r, c) for r in range(dA*dB) for c in range(dA*dB) if wout[r] == win[c]]
    blocks = []
    for j in js:
        for kind in 'ef':
            D, Dp = cop(A, B, kind, j), cop(B, A, kind, j); rows, cc, vals = [], [], []
            for i, (r, c) in enumerate(cols):
                for t in np.nonzero(D[c, :])[0]:
                    rows.append(r*dA*dB + t); cc.append(i); vals.append(D[c, t])
                for t in np.nonzero(Dp[:, r])[0]:
                    rows.append(t*dA*dB + c); cc.append(i); vals.append(-Dp[t, r])
            blocks.append(sps.csr_matrix((vals, (rows, cc)), shape=((dA*dB)**2, len(cols))))
    M = sps.vstack(blocks); ev, U = np.linalg.eigh((M.conj().T @ M).toarray())
    R = np.zeros((dA*dB, dA*dB), complex)
    for i, (r, c) in enumerate(cols):
        R[r, c] = U[i, 0]
    return R/R[B['hw']*dA + A['hw'], A['hw']*dB + B['hw']]
nus = {'id': {j: j for j in range(n + 1)}, '(n-1 n)': {**{j: j for j in range(n + 1)}, n - 1: n, n: n - 1}}
y = np.exp(0.6)
for name, nu in nus.items():
    V = dn_vector(n, q, y)
    pairs = [(V[(k, nu[j])], V[(k, j)].conj().T) for j in range(n + 1) for k in 'efk']
    null, _ = intertwiners(pairs)
    print(f'nu={name}: dim of eta-space {len(null)}')
    if len(null) != 1:
        continue
    eta = null[0]; eta = eta/eta[np.unravel_index(np.argmax(np.abs(eta)), eta.shape)]
    print('   eta Hermitian:', np.allclose(eta, eta.conj().T), ' eta^2 prop 1:', np.allclose(eta @ eta, (eta @ eta)[0, 0]*np.eye(N)),
          ' nonzero pattern (row -> col):', [int(np.argmax(np.abs(eta[i]))) for i in range(N)],
          ' entries:', np.round([eta[i, np.argmax(np.abs(eta[i]))] for i in range(N)], 3))
    ee = np.kron(eta, eta)
    for th in (0.3, 1.2):
        R = rcheck_w(dn_vector(n, q, np.exp(om*th)), dn_vector(n, q, 1.0))
        M = R.conj().T @ ee @ R; c = np.vdot(ee.ravel(), M.ravel())/np.vdot(ee.ravel(), ee.ravel())
        print(f'   R-check th={th}: R^dag (eta x eta) R = c (eta x eta): err {np.abs(M - c*ee).max():.1e}, c = {c:.6f}')
    for lab, e in (('eps=0', np.zeros(n + 1)), ('all +eps*', es*np.ones(n + 1)), ('eps_0 = -eps*', es*np.r_[-1, np.ones(n)]),
                   ('eps_n = -eps*', es*np.r_[np.ones(n), -1]), ('eps_{n-1} = -eps*', es*np.r_[np.ones(n - 1), -1, 1])):
        out = []
        for th in (0.3, 1.2):
            yy = np.exp(om*th)
            K = k_intertwiner(dn_vector(n, q, yy), dn_vector(n, q, 1/yy), [1/q]*(n + 1), e)[0][0]
            M = K.conj().T @ eta @ K; c = np.vdot(eta.ravel(), M.ravel())/np.vdot(eta.ravel(), eta.ravel())
            out.append(f'th={th}: err {np.abs(M - c*eta).max()/np.abs(M).max():.1e} c={c.real:+.3f}{c.imag:+.0e}i')
        print(f'   K {lab:16s}: ' + '; '.join(out))
