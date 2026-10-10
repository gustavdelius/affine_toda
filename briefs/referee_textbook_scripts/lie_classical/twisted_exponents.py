"""Principal Heisenberg degrees ('exponents') of twisted affine algebras X_N^(k):
eigenvalues of the twisted Coxeter element c = sigma * prod_{sigma-orbits O} s_{i_O} on the Cartan of X_N
(Springer; Kac ch. 14). Print degrees m with eigenvalue exp(2 pi i m/P), P = order of c, and compare P with k*h."""
import numpy as np
from lie import *

def twisted(typ, N, sigma, k, h_book, expect):
    R = RootSystem(typ, N)
    S = np.linalg.solve(R.al, R.al[sigma]).T          # S alpha_i = alpha_sigma(i)
    assert np.allclose(S @ R.al.T, R.al[sigma].T) and np.allclose(S @ S.T, np.eye(N))
    seen = set(); reps = []
    for i in range(N):
        if i in seen: continue
        o = [i]; j = sigma[i]
        while j != i: o.append(j); j = sigma[j]
        seen |= set(o); reps.append(i)
    c = S.copy()
    for i in reps: c = c @ refl(R.al[i])
    P = order_of(c)
    m = sorted(int(round(x)) for x in exponents_of(c, P))
    print(f'{typ}{N} sigma of order {k}: twisted Coxeter order P={P} (k*h_book={k*h_book}, h_book={h_book}); '
          f'degrees mod P: {m}; expected {expect}: {"PASS" if (P == k*h_book and m == expect) else "FAIL"}')

# a_2^(2): A2 flip
twisted('A', 2, [1, 0], 2, 3, [1, 5])
twisted('A', 4, [3, 2, 1, 0], 2, 5, [1, 3, 7, 9])                    # a_4^(2): odd, not divisible by 5, mod 10
twisted('A', 6, [5, 4, 3, 2, 1, 0], 2, 7, [1, 3, 5, 9, 11, 13])
twisted('A', 3, [2, 1, 0], 2, 3, [1, 3, 5])                           # a_3^(2): odd mod 6 (= d_3^(2))
twisted('A', 5, [4, 3, 2, 1, 0], 2, 5, [1, 3, 5, 7, 9])               # a_5^(2): odd mod 10
twisted('A', 7, [6, 5, 4, 3, 2, 1, 0], 2, 7, [1, 3, 5, 7, 9, 11, 13])
twisted('D', 4, [0, 1, 3, 2], 2, 4, [1, 3, 5, 7])                     # d_4^(2): odd mod 8
twisted('D', 5, [0, 1, 2, 4, 3], 2, 5, [1, 3, 5, 5, 7, 9][:0] or [1, 3, 5, 7, 9])  # d_5^(2): odd mod 10
twisted('D', 4, [2, 1, 3, 0], 3, 4, [1, 5, 7, 11])                    # d_4^(3)
twisted('E', 6, [5, 1, 4, 3, 2, 0], 2, 9, [1, 5, 7, 11, 13, 17])     # e_6^(2)
