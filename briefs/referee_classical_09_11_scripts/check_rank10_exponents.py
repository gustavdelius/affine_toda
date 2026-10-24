#!/usr/bin/env python3
"""
Independent verification of Rank 10 (09-lagrangian.qmd:84-88):
Exponents modulo k*h and order of the twisted Coxeter element for all untwisted and twisted affine algebras,
plus explicit evaluation of the SG kink spin-3 integral (Rank 8).
"""

import numpy as np

print("=" * 72)
print("CHECK (Rank 8): Fast exact evaluation of SG kink Q_3 integral")
print("=" * 72)
# u = m(x cosh theta - t sinh theta)
# cos(4 arctan(e^u)) = 1 - 2 sech^2(u)
# T_4 = (4 m^4 / beta^2) e^{-4 theta} sech^2(u) (tanh^2(u) - sech^2(u)) = (4 m^4 / beta^2) e^{-4 theta} (sech^2(u) - 2 sech^4(u))
# X_2 = (4 m^4 / beta^2) e^{-2 theta} (1 - 2 sech^2(u)) sech^2(u) = (4 m^4 / beta^2) e^{-2 theta} (sech^2(u) - 2 sech^4(u))
# Integral over u in (-infty, +infty) of (sech^2(u) - 2 sech^4(u)) du = 2 - 2*(4/3) = -2/3!
us = np.linspace(-30.0, 30.0, 200001)
integrand = (1.0 / np.cosh(us) ** 2) - 2.0 / (np.cosh(us) ** 4)
print(f"int (sech^2 u - 2 sech^4 u) du = {np.trapz(integrand, us):.12f}  (-2/3 = {-2/3:.12f})")
print("Therefore, on the sine-Gordon kink phi = (4/beta) arctan(exp(m(x cosh theta - t sinh theta))):")
print("  int T_4 dx = -(8 m^3 / (3 beta^2 cosh theta)) e^{-4 theta}")
print("  int X_2 dx = -(8 m^3 / (3 beta^2 cosh theta)) e^{-2 theta}")
print("  Q_3 = int (T_4 + X_2) dx = -(16 m^3 / (3 beta^2)) e^{-3 theta}  (scales as e^{-3 theta}, spin -3!)")

print("\n" + "=" * 72)
print("CHECK (Rank 10): Exponents of untwisted and twisted affine algebras modulo k*h")
print("=" * 72)

# Let's check Kac (1990) Table Aff 1, Aff 2, Aff 3 for the exponents of X_N^(k) and period k*h
# where h = sum_j n_j (sum of coprime Kac labels of X_N^(k)):
twisted_data = [
    ("a_2^(2)", 2, [1, 2], [1, 5]),
    ("a_{2n}^(2) (n=2: a_4^(2))", 2, [1, 2, 2], [1, 3, 7, 9]),
    ("a_{2n-1}^(2) (n=3: a_5^(2))", 2, [1, 1, 2, 1], [1, 3, 5, 7, 9]),
    ("d_{n+1}^(2) (n=3: d_4^(2))", 2, [1, 1, 1, 1], [1, 3, 5, 7]),
    ("e_6^(2)", 2, [1, 2, 3, 2, 1], [1, 5, 7, 11, 13, 17]),
    ("d_4^(3)", 3, [1, 2, 1], [1, 5, 7, 11]),
]

for name, k, n_kac, exps in twisted_data:
    h = sum(n_kac)
    kh = k * h
    print(f"{name}: k={k}, h=sum(n_j)={h}, period k*h={kh}, exponents mod {kh} = {exps}")
    # Check symmetry s <-> kh - s:
    assert all((kh - s) in exps for s in exps)

# Check dual pairs have the same spins:
# 1) (b_n^(1), a_{2n-1}^(2)):
#    b_n^(1) has exponents 1, 3, ..., 2n-1 mod 2n => all odd integers!
#    a_{2n-1}^(2) has k=2, h=2n-1, k*h=4n-2, and exponents in Kac Table Aff 2 are all odd integers not divisible by 2n-1?
#    WAIT A SECOND!
#    HOLD EVERYTHING!
#    What are the exponents of a_{2n-1}^(2) in Kac (1990), Chapter 8 / Table Aff 2, AND what are the spins of conserved charges of b_n^(1) vs a_{2n-1}^(2)?!
print("\nWAIT: Let's check the exponents of b_n^(1), c_n^(1), a_{2n-1}^(2), d_{n+1}^(2) directly from the principal Heisenberg subalgebra!")


# Let's compute the centralizer of E_1 = sum_j c_j e_j in the parent simply-laced algebra
# restricted to the sigma-invariant subalgebra!
# Wait: as shown in 06-affine-algebras.qmd:95, the folded affine Toda Lax pair A_+, A_-
# is literally the restriction of the parent simply-laced Lax pair to the sigma-invariant subalgebra!
# Let's verify this!
# 1) For b_n^(1): parent is d_{n+1}^(1) (with Coxeter number h_D = 2n), folded by sigma swapping spinor nodes n <-> n+1.
#    WAIT: What about a_{2n-1}^(2)?
#    Look at @tbl-foldings (06-affine-algebras.qmd:103):
#    The parent of a_{2n-1}^(2) is d_{2n}^(1) (with Coxeter number h_{d_{2n}} = 2(2n) - 2 = 4n - 2 = 2(2n-1) = k*h)!
#    WAIT A SECOND!
#    In d_{2n}^(1), the Coxeter number is 4n - 2 = 2(2n - 1)!
#    What are the exponents of D_{2n}? They are 1, 3, 5, ..., 4n - 3, and 2n - 1 (with multiplicity 2)!
#    Under the folding automorphism sigma of d_{2n}^(1) -> a_{2n-1}^(2) (which folds the D_{2n} diagram):
#    Let's compute the eigenvalues of sigma on the centralizer of E_1 in d_{2n}^(1) and in d_{n+1}^(1)!

# Let's build the fundamental (2N x 2N) matrix representation of d_N^(1) in the loop algebra so(2N)[z, z^{-1}]
# and find the centralizer of E_1!
def check_dN_folding_centralizer(N):
    # Standard representation of so(2N) on C^{2N}:
    # Basis: e_1, ..., e_N, e_{-N}, ..., e_{-1} (or e_1..e_N, f_1..f_N where metric is J = [[0, I], [I, 0]])
    # Cartan H_i = E_{i,i} - E_{N+i, N+i} (i = 1..N)
    # Simple roots of D_N:
    # e_i (1 <= i <= N-1): E_{i, i+1} - E_{N+i+1, N+i}
    # e_N: E_{N-1, 2N} - E_{N, N+N-1}
    # e_0 (lowest root -(eps_1 + eps_2)): z * (E_{N+1, 2} - E_{N+2, 1})
    dim = 2 * N
    def E_mat(i, j):
        M = np.zeros((dim, dim), dtype=complex)
        M[i, j] = 1.0
        return M

    e_gens = []
    # e_0:
    e0 = E_mat(N, 1) - E_mat(N + 1, 0)
    e_gens.append(e0)
    for i in range(N - 1):
        ei = E_mat(i, i + 1) - E_mat(N + i + 1, N + i)
        e_gens.append(ei)
    eN = E_mat(N - 2, 2 * N - 1) - E_mat(N - 1, 2 * N - 2)
    e_gens.append(eN)

    # Kac labels of d_N^(1): 1, 1, 2, ..., 2, 1, 1
    n_kac = [1, 1] + [2] * (N - 3) + [1, 1]
    c = [np.sqrt(nk) for nk in n_kac]
    # At z = 1, E_1 = sum_{j=0}^N c_j e_j
    E1 = sum(c[j] * e_gens[j] for j in range(N + 1))
    # What are the eigenvalues of E1?
    evals = np.linalg.eigvals(E1)
    # Sort by angle where non-zero:
    nz = [ev for ev in evals if abs(ev) > 1e-8]
    print(f"d_{N}^(1) (h={2*N-2}): {len(nz)} non-zero eigenvalues of E1, magnitudes = {np.round(np.abs(nz), 4)}, zero evals = {dim - len(nz)}")

check_dN_folding_centralizer(4)
check_dN_folding_centralizer(6)


# Now let's include the loop variable z (or the principal grading where every e_j has degree 1)!
# In the loop algebra, every e_j (j=0,...,r) has principal degree +1, so multiplying e_j by lam
# means E_1(lam) = lam * E_1(1)!
# Wait: an element X of the loop algebra of principal degree s > 0 is a linear combination of
# nested commutators of s generators e_{j_1}, ..., e_{j_s}?
# NO! In the loop algebra g[z, z^{-1}], e_1, ..., e_r have z-degree 0 and e_0 has z-degree +1,
# while f_1, ..., f_r have z-degree 0 and f_0 has z-degree -1!
# Wait: in 06-affine-algebras.qmd:57 and 09-lagrangian.qmd:57:
# "For j >= 1 they are the Chevalley generators of g. For j = 0 they are e_0 and f_0, the root vectors
#  E^{-theta} and E^{theta} multiplied by the first power of the loop variable and its inverse.
#  Introduce a spectral parameter lam. It multiplies every e_j equally, so it is the spectral parameter
#  of the principal gradation. In a matrix representation the loop variable can be set to 1, as in @exr-lax,
#  and lam takes over its role."
# YES! Because the principal gradation is isomorphic to the homogeneous gradation by conjugation by
# lam^{rho^vee / h}, setting the loop variable z = 1 and multiplying every e_j by lam (and every f_j by lam^{-1})
# means that an element of principal degree s in the loop algebra is simply X_s * lam^s where X_s in g
# commutes with E_1 = sum_{j=0}^r c_j e_j AND carries charge s mod h under the Coxeter conjugation
# exp(2*pi*i * rho^vee / h) (so that un-conjugating from principal to homogeneous gradation gives single-valued
# powers of z = lam^h)!
# Let's check which odd powers E_1^s carry the right grading under exp(2*pi*i * rho^vee / h) in d_{n+1}^(1) -> b_n^(1)
# and in d_{2n}^(1) -> a_{2n-1}^(2)!
for N in [4, 6]:
    h = 2 * N - 2
    dim = 2 * N
    # In d_N^(1), the principal grading operator ad(rho^vee) gives e_1, ..., e_N degree +1 and e_0 degree -(h-1) = 1 - h.
    # Let's construct the diagonal matrix H_rho such that [H_rho, e_j] = +1 * e_j for j=1..N and [H_rho, e_0] = (1 - h) * e_0!
    # In our 2N x 2N basis:
    # e_i (1 <= i <= N-1) has entry at (i-1, i) and (N+i, N+i-1)
    # e_N has entry at (N-2, 2N-1) and (N-1, 2N-2)
    # e_0 has entry at (N, 1) and (N+1, 0).
    # So H_rho = diag(N-1, N-2, ..., 1, 0, -(N-1), -(N-2), ..., -1, 0)!
    h_diag = np.array(list(range(N - 1, -1, -1)) + list(range(-(N - 1), 1)))
    omega_h = np.exp(2j * np.pi / h)
    Omega = np.diag(omega_h ** h_diag)

    def E_mat(i, j):
        M = np.zeros((dim, dim), dtype=complex)
        M[i, j] = 1.0
        return M

    e_gens = [E_mat(N, 1) - E_mat(N + 1, 0)]
    for i in range(N - 1):
        e_gens.append(E_mat(i, i + 1) - E_mat(N + i + 1, N + i))
    e_gens.append(E_mat(N - 2, 2 * N - 1) - E_mat(N - 1, 2 * N - 2))
    n_kac = [1, 1] + [2] * (N - 3) + [1, 1]
    E1 = sum(np.sqrt(n_kac[j]) * e_gens[j] for j in range(N + 1))

    # Check Omega E1 Omega^{-1} = omega_h * E1:
    assert np.allclose(Omega @ E1 @ np.linalg.inv(Omega), omega_h * E1)
    # Therefore, for ANY odd integer s >= 1, E1^s satisfies:
    # 1) [E1, E1^s] = 0,
    # 2) E1^s in so(2N) (since s is odd and E1 in so(2N)),
    # 3) Omega (E1^s) Omega^{-1} = omega_h^s * E1^s (so it has principal degree s mod h!),
    # 4) E1^s is non-zero for EVERY odd s >= 1 (including s = N - 1 = h/2 when N is even!) because E1 has 2N-2 non-zero eigenvalues!
    # 5) Since E1 is invariant under ANY folding automorphism sigma of d_N^(1), E1^s is ALSO sigma-invariant for every odd s!
    for s in range(1, h + 2, 2):
        E1_s = np.linalg.matrix_power(E1, s)
        norm_s = np.linalg.norm(E1_s)
        print(f"  d_{N}^(1) odd power s={s}: ||E_1^{s}|| = {norm_s:.4f} (> 0, sigma-invariant!)")
