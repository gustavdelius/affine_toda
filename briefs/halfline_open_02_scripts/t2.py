from dn_common import *
from boundary_common import coideal_gen_tensor
np.set_printoptions(linewidth=200, precision=4, suppress=True)
n = 4; N = 2*n
om = 0.37; q = q_phys(om)
lab = basis_labels(n)
def inv_vecs(ya, yb, eps):
    A, B = dn_vector(n, q, ya), dn_vector(n, q, yb)
    M = np.vstack([coideal_gen_tensor(A, B, j, 1/q, eps[j]) - eps[j]*np.eye(N*N) for j in range(n + 1)])
    _, s, Vh = np.linalg.svd(M)
    return s[-3:]/s[0], Vh[-1].conj()
y = 0.8 + 0.3j
for (ya, yb, name) in ((1/y, y, '1/y,y'), (y, 1/y, 'y,1/y'), (-1/y, y, '-1/y,y'), (1/y/q, y, '1/(yq),y'), (-1/(y*q), y, '-1/(yq), y')):
    s, v = inv_vecs(ya, yb, np.zeros(n + 1))
    print(name, s)
    if s[-1] < 1e-10:
        nz = np.nonzero(np.abs(v) > 1e-8)[0]
        v = v/v[nz[0]]
        print('   support', [(lab[i // N], lab[i % N], np.round(v[i], 4)) for i in nz])
