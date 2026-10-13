from dn_common import *
np.set_printoptions(linewidth=200, precision=4, suppress=True)
n = 4; N = 2*n
om = 0.37; q = q_phys(om)
lab = basis_labels(n)
Ks = []
for y in (0.8 + 0.3j, 1.7 - 0.2j, 0.5j + 0.3):
    K = Kdn(n, q, y, np.zeros(n + 1))
    # normalize by the (0, N-1) entry
    K = K/K[N - 1, 0]
    Ks.append(K)
    print('y', y, 'antidiag', np.round(K[np.arange(N)[::-1], np.arange(N)], 6))
print('y-independent:', max(np.abs(Ks[0] - Ks[i]).max() for i in (1, 2)))
J = Ks[0]
print('J^2 =', np.round(J @ J, 6).diagonal())
print('labels', lab)
