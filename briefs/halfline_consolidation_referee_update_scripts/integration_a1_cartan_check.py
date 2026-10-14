import sympy as sp
from integration_a1_dirichlet import *
for name, rep in [('a2', sl_vector(2)), ('a3', sl_vector(3)), ('d4', so_vector(8)), ('b2', so_vector(5)), ('c2', sp4_vector())]:
    N, Ep, Em, H, E0p, E0m, marks = rep
    Pp = [E0p] + Ep; Pm = [E0m] + Em
    Hs = [comm(a, b) for a, b in zip(Pp, Pm)]
    # affine Cartan matrix A_ij from [h_i, e_j] = A_ij e_j
    A = sp.zeros(len(Pp))
    for i in range(len(Pp)):
        for j in range(len(Pp)):
            C = comm(Hs[i], Pp[j])
            k = [t for t in range(N*N) if Pp[j][t] != 0][0]
            A[i, j] = C[k]/Pp[j][k]
            assert C == A[i, j]*Pp[j]
    # Serre-type check: sum_j n_j h_j proportional to 0 with coroot marks
    print(name, A.tolist())
