"""Which ordering of V(theta_k), mu_k = mu + T - k omega, contains Vbar(theta) as a subrepresentation?"""
import numpy as np
from sp_lib import an_vector, q_phys, tensor_rep, nodes_of
from boundary_common import intertwiners
for n in (1, 2, 3):
    N = n + 1
    for om in (0.37, 1.29):
        T, q = N*om/2, q_phys(om)
        mu = 0.13 + 0.21j
        th = 1j*np.pi*mu/T
        for order in (+1, -1):
            mus = [mu + T - k*om for k in range(1, n + 1)][::order]
            Vs = [an_vector(n, q, np.exp(om*1j*np.pi*m/T)) for m in mus]
            W = tensor_rep(Vs)
            Ub = an_vector(n, q, np.exp(om*th), conj=True)
            pairs = [(Ub[(kd, j)], W[(kd, j)]) for j in nodes_of(Ub) for kd in 'efk']
            null, rel = intertwiners(pairs, tol=1e-9)
            print(n, om, 'order', order, 'dim Hom(Vbar, V^n) =', len(null), 'smallest rel sv', [f'{r:.1e}' for r in rel[:3]])
