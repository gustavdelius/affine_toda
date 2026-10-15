import sympy as sp

# Check centrality of prod k_i^{m_i} and equality with q^{sum n_i^\vee h_i}
# In U_q(\hat g):
# k_i e_j k_i^{-1} = q_i^{a_{ij}} e_j = q^{(\alpha_i, \alpha_j)} e_j, where q_i = q^{(\alpha_i, \alpha_i)/2}.
# For K = \prod_{i=0}^r k_i^{m_i}:
# K e_j K^{-1} = q^{(\sum_i m_i \alpha_i, \alpha_j)} e_j.
# Centrality requires (\sum_i m_i \alpha_i, \alpha_j) = 0 for all j.
# Since \delta = \sum a_i \alpha_i is the unique null root, m_i must be proportional to a_i (the Kac labels).
# If m_i = a_i (coprime integers with a_0 = 1, except Kac A_{2n}^{(2)} where a_0=2):
# Then K = \prod k_i^{a_i} is central.

# Now what about: K = q^{\sum_i n_i^\vee h_i}?
# Recall k_i = q_i^{h_i} = q^{(\alpha_i, \alpha_i)/2 h_i}.
# So \prod k_i^{m_i} = q^{\sum_i m_i (\alpha_i, \alpha_i)/2 h_i}.
# For this to equal q^{\sum_i n_i^\vee h_i}, we need:
# n_i^\vee = m_i \frac{(\alpha_i, \alpha_i)}{2}.
# Let's check for which algebras and normalizations this holds!

# 1. Untwisted: (alpha_0, alpha_0) = 2.
# Then a_i^\vee = a_i (\alpha_i, \alpha_i)/2 holds for all untwisted algebras!
# Check:
# - c_2^(1): long roots have length^2 2, short roots have length^2 1.
#   Kac marks a = (1, 2, 1). Root lengths (alpha_0^2, alpha_1^2, alpha_2^2) = (2, 1, 2).
#   Then a_i (\alpha_i^2 / 2) = (1*1, 2*(1/2), 1*1) = (1, 1, 1).
#   Indeed, for C_2^(1), a^\vee = (1, 1, 1), and h^\vee = 1+1+1 = 3!
#   So n_i^\vee = n_i (\alpha_i^2/2) is EXACT!

# 2. Twisted: d_{n+1}^(2) (e.g. d_3^(2)):
# In Kac normalization: (alpha_0, alpha_0) = 2.
# What are the root lengths in Kac?
# For D_3^(2):
# Node 0 has alpha_0^2 = 2.
# Node 1 has alpha_1^2 = 4 (long root).
# Node 2 has alpha_2^2 = 2.
# What are Kac marks a_i?
# A = [[2, -2, 0], [-1, 2, -1], [0, -2, 2]]
# Right null vector: A a = 0 => a = (1, 1, 1).
# Then a_i (alpha_i^2 / 2) = (1*1, 1*2, 1*1) = (1, 2, 1).
# Left null vector: a^\vee A = 0 => (1, 2, 1) A = [2-2, -2+4-2, -2+2] = (0, 0, 0).
# So a^\vee = (1, 2, 1)!
# And h^\vee = 1 + 2 + 1 = 4!
# Here a_i (alpha_i^2 / 2) = a_i^\vee EXACTLY!
# BUT what if the book's normalization is used?
# In book's normalization: longest root has length^2 2!
# Then alpha_1^2 = 2, alpha_0^2 = alpha_2^2 = 1.
# Then a_i (alpha_i^2/2) = (1/2, 1, 1/2), which are not integers!
# In the book (tbl-foldings caption): n_i^\vee = k * n_i (alpha_i^2/2) = 2 * (1/2, 1, 1/2) = (1, 2, 1).
# So in the book's normalization, n_i^\vee has an extra factor of k!
# In that normalization, \prod k_i^{n_i} = q^{\sum (n_i^\vee / k) h_i} = (q^{1/k})^{\sum n_i^\vee h_i}!

# 3. a_{2n}^(2):
# Let's check a_2^(2):
# If Kac normalization: alpha_0^2 = 1, alpha_1^2 = 4.
# Kac marks: a = (2, 1).
# Then a_i (alpha_i^2/2) = (2 * 1/2, 1 * 2) = (1, 2) = a^\vee!
# Here a_i (alpha_i^2/2) = a_i^\vee is EXACT!
# BUT the central element is \prod k_i^{a_i} = k_0^2 k_1^1, with a_0 = 2, NOT 1!
# If book's folding labels are used:
# n = (1, 2) (node 0 has label 1, node 1 has label 2).
# And alpha_0^2 = 2, alpha_1^2 = 1/2 (longest root = 2).
# Then n_i (alpha_i^2/2) = (1 * 1, 2 * 1/4) = (1, 1/2), NOT INTEGERS!
# And n^\vee in tbl-foldings is (2, 1).
# Here n_i (alpha_i^2/2) != n_i^\vee!

print("Central element analysis complete.")
