import sympy as sp

# Let's inspect:
# For any affine algebra with simple roots alpha_i and invariant bilinear form (.,.):
# The Cartan matrix is A_ij = 2(alpha_i, alpha_j)/(alpha_i, alpha_i).
# A has rank r, with rank deficiency 1.
# There exist positive integers a_i (Kac labels / marks) such that:
# \sum_j A_ij a_j = 0  (right null vector, \sum a_j alpha_j^\vee? No, \sum a_j alpha_j is null root delta in Kac!)
# Wait, let's verify Kac's convention:
# In Kac, delta = \sum_{i=0}^r a_i alpha_i is the null root: (delta, alpha_j) = 0 for all j.
# (delta, alpha_j) = \sum_i a_i (alpha_i, alpha_j) = 0.
# The Cartan matrix A_ij = 2(alpha_i, alpha_j)/(alpha_i, alpha_i).
# So (alpha_i, alpha_j) = (alpha_i, alpha_i)/2 * A_ij.
# Therefore, \sum_i a_i (alpha_i, alpha_i)/2 * A_ij = 0.
# This means: the vector with components a_i (alpha_i, alpha_i)/2 is a LEFT null vector of A!
# In Kac, the dual Kac labels a_i^\vee (or marks of A^T) satisfy:
# \sum_i a_i^\vee A_ij = 0, with a_0^\vee = 1 (except for A_{2n}^{(2)} where a_0^\vee = 2? or a_0=2?).
# Let's check Kac p. 54-55!

# Let's write out the Cartan matrices and check for:
# 1. c_2^(1) (untwisted non-simply-laced): C_2^(1)
# 2. d_3^(2) (twisted): D_3^(2) = A_3^(2)? (Note D_3^(2) is often called C_2^(1) dual or B_2^(1) dual)
# 3. a_2^(2) (twisted, A_2^(2))
# Let's verify each explicitly.
