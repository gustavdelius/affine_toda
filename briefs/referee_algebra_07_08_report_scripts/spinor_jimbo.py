import numpy as np
import scipy.linalg
import itertools

# Build spinor rep of U_q(d_{n+1}^(2)) and solve Jimbo's equations for n=2 and n=3.

Qs = 0.83 * np.exp(0.61j)

def Emat(N, a, b):
    M = np.zeros((N, N), complex)
    M[a, b] = 1.0
    return M

def qnum(n, q):
    return (q**n - q**(-n)) / (q - 1.0/q)

def qfac(n, q):
    prod = 1.0 + 0.0j
    for i in range(1, n+1):
        prod *= qnum(i, q)
    return prod

def qbin(n, k, q):
    return qfac(n, q) / (qfac(k, q) * qfac(n-k, q))

def make_spinor_rep(n, Q=Qs):
    # Basis states: tuples of +1 or -1 (representing +/- 1/2)
    dim = 2**n
    basis = list(itertools.product([1, -1], repeat=n))
    basis_idx = {b: i for i, b in enumerate(basis)}
    
    # Root lengths in Kac:
    # alpha_0^2 = 2 => q_0 = Q
    # alpha_i^2 = 4 for 1 <= i < n => q_i = Q^2
    # alpha_n^2 = 2 => q_n = Q
    qi = {}
    qi[0] = Q
    for i in range(1, n):
        qi[i] = Q**2
    qi[n] = Q
    
    # Cartan matrix of d_{n+1}^(2)
    # Dynkin diagram: 0 <= 1 - 2 - ... - (n-1) => n
    # For n=2: 0 <= 1 => 2 (double bonds at both ends)
    # A_00=2, A_01=-2, A_10=-1, A_11=2, A_12=-1, A_21=-2, A_22=2
    cartan = {i: {j: 0 for j in range(n+1)} for i in range(n+1)}
    for i in range(n+1):
        cartan[i][i] = 2
    for i in range(1, n-1):
        cartan[i][i+1] = -1
        cartan[i+1][i] = -1
    # 0 and 1:
    cartan[0][1] = -2
    cartan[1][0] = -1
    # n-1 and n:
    cartan[n-1][n] = -1
    cartan[n][n-1] = -2

    def rep(x):
        e = {i: np.zeros((dim, dim), complex) for i in range(n+1)}
        f = {i: np.zeros((dim, dim), complex) for i in range(n+1)}
        k = {i: np.zeros((dim, dim), complex) for i in range(n+1)}
        
        # k_i diagonal entries
        # In Kac normalization:
        # For i=1..n-1: alpha_i = e_i - e_{i+1}, length^2 4.
        # h_i = 2 (e_i - e_{i+1}) / 4 * 2 = e_i - e_{i+1}
        # k_i = q_i^{h_i} = (Q^2)^{ (sigma_i - sigma_{i+1})/2 } = Q^{sigma_i - sigma_{i+1}}
        # For i=n: alpha_n = e_n, length^2 2.
        # h_n = 2 e_n / 2 = 2 e_n
        # k_n = q_n^{h_n} = Q^{ 2 sigma_n / 2 } = Q^{sigma_n}
        # For i=0: alpha_0 = -e_1, length^2 2.
        # h_0 = -2 e_1 / 2 = -2 e_1
        # k_0 = q_0^{h_0} = Q^{-sigma_1}
        for idx, b in enumerate(basis):
            # b is tuple of +/- 1
            # k_0:
            k[0][idx, idx] = Q**(-b[0])
            for i in range(1, n):
                k[i][idx, idx] = Q**(b[i-1] - b[i])
            k[n][idx, idx] = Q**(b[n-1])
            
            # e_n, f_n: flips b[n-1]
            if b[n-1] == -1:
                b_new = list(b)
                b_new[n-1] = 1
                new_idx = basis_idx[tuple(b_new)]
                e[n][new_idx, idx] = 1.0
                f[n][idx, new_idx] = 1.0
                
            # e_i, f_i for 1 <= i < n: flips b[i-1] from -1 to 1 and b[i] from 1 to -1
            for i in range(1, n):
                if b[i-1] == -1 and b[i] == 1:
                    b_new = list(b)
                    b_new[i-1] = 1
                    b_new[i] = -1
                    new_idx = basis_idx[tuple(b_new)]
                    e[i][new_idx, idx] = 1.0
                    f[i][idx, new_idx] = 1.0
                    
            # e_0, f_0: flips b[0] from 1 to -1 (carrying x)
            if b[0] == 1:
                b_new = list(b)
                b_new[0] = -1
                new_idx = basis_idx[tuple(b_new)]
                e[0][new_idx, idx] = x * 1.0
                f[0][idx, new_idx] = (1.0 / x) * 1.0
                
        return e, f, k
    
    return rep, cartan, qi

import test_all_three as T

for n in (2, 3):
    rep, cartan, qi = make_spinor_rep(n)
    res = T.check_algebra(rep, cartan, qi)
    print(f"Spinor representation n={n} (dim={2**n}): defining relations residual = {res:.2e}")
