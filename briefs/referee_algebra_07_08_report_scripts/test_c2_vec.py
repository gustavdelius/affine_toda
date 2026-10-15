import sympy as sp
import numpy as np

# Test c_2^(1) vector representation (4-dimensional)
# Finite algebra: c_2 = sp(4).
# Simple roots: alpha_1 (short, length^2 1 in book's long^2=2 normalization), alpha_2 (long, length^2 2).
# Or in Kac: alpha_1^2 = 1, alpha_2^2 = 2.
# Defining rep V = C^4 (highest weight lambda_1).
# Decomposition of V \otimes V for c_2:
# 4 \otimes 4 = 10_sym \oplus 5_anti \oplus 1_sym.
# Highest weights:
# 10: 2 lambda_1 (symmetric, parity +)
# 5: lambda_2 (antisymmetric, parity -)
# 1: 0 (symmetric, parity +)
#
# Let's compute Casimirs C(mu) with long roots of squared length 2:
# In c_2:
# alpha_1 = (1, -1), alpha_2 = (0, 2) in orthonormal basis?
# Wait: standard orthonormal basis e_1, e_2:
# Short roots: +/- e_1 +/- e_2, length^2 2?
# In book's normalization, long roots have length^2 2!
# So short roots have length^2 1!
# Let's use e_1, e_2 with (e_i, e_j) = 1/2 delta_ij:
# Long roots: +/- 2 e_i: squared length = 4 * (1/2) = 2.
# Short roots: +/- e_1 +/- e_2: squared length = (1+1)*(1/2) = 1.
# Simple roots: alpha_1 = e_1 - e_2 (short, length^2 1), alpha_2 = 2 e_2 (long, length^2 2).
# Fundamental weights:
# lambda_1 = e_1:
# <lambda_1, alpha_1^\vee> = 2(e_1, e_1 - e_2)/1 = 2(1/2)/1 = 1.
# <lambda_1, alpha_2^\vee> = 2(e_1, 2e_2)/2 = 0.
# lambda_2 = e_1 + e_2:
# <lambda_2, alpha_1^\vee> = 2(e_1+e_2, e_1-e_2)/1 = 0.
# <lambda_2, alpha_2^\vee> = 2(e_1+e_2, 2e_2)/2 = 2(1)/2 = 1.
# Weyl vector rho = lambda_1 + lambda_2 = 2 e_1 + e_2.
#
# Casimir C(mu) = (mu, mu + 2 rho):
# 1. mu = 2 lambda_1 = 2 e_1:
#    mu + 2 rho = 2 e_1 + 4 e_1 + 2 e_2 = 6 e_1 + 2 e_2.
#    C(2 lambda_1) = (2 e_1, 6 e_1 + 2 e_2) = 2 * 6 * (1/2) = 6.
# 2. mu = lambda_2 = e_1 + e_2:
#    mu + 2 rho = e_1 + e_2 + 4 e_1 + 2 e_2 = 5 e_1 + 3 e_2.
#    C(lambda_2) = (e_1 + e_2, 5 e_1 + 3 e_2) = (5 + 3) * (1/2) = 4.
# 3. mu = 0:
#    C(0) = 0.

# Parities:
# 2 lambda_1 (10): parity +
# lambda_2 (5): parity -
# 0 (1): parity +

# Edges:
# Adjoint of c_2 is 10 = 2 lambda_1.
# Between 2 lambda_1 (+) and lambda_2 (-): opposite parity!
# Edge difference: (C(lambda_2) - C(2 lambda_1)) / 4 = (4 - 6) / 4 = -2/4 = -1/2!
# So rho_{lambda_2} / rho_{2 lambda_1} = <-1/2> = <1/2>!
# Between lambda_2 (-) and 0 (+): opposite parity!
# Edge difference: (C(0) - C(lambda_2)) / 4 = (0 - 4) / 4 = -4/4 = -1!
# Wait: does 0 occur in lambda_2 \otimes ad?
# lambda_2 \otimes 10 = 5 \otimes 10 = 5 \oplus ... (does not contain 0, because 5 is not 10).
# Does lambda_2 occur in 0 \otimes ad?
# 0 \otimes 10 = 10, not 5.
# So there is NO edge between lambda_2 and 0!
# What about between 2 lambda_1 and 0?
# 2 lambda_1 and 0 both have parity + (same parity)! So no untwisted edge!
# Wait! How is rho_0 linked to rho_{2 lambda_1}?
# In c_2^(1), the extended graph has:
# Jimbo's formula for c_n^(1) vector:
# Rcheck = P_{2 lambda_1} + <1/2> P_{lambda_2} + <(n+1)/2> P_0!
# For n=2: <(2+1)/2> = <3/2>!
# So rho_0 = <3/2>!

print("Casimir values for c_2:")
print("C(2 lambda_1) = 6")
print("C(lambda_2) = 4")
print("C(0) = 0")
print("Ratio rho_{lambda_2}/rho_{2 lambda_1} = <1/2> (half-integer argument!)")
