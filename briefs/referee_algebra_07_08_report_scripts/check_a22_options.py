# Verify Kac's data for A_{2l}^{(2)} and in particular A_2^(2) (l=1)
# For A_2^(2), there are 2 nodes: 0 and 1.
# The Cartan matrix is 2x2.
# In Kac Table Aff 2:
# For l=1, A_2^(2):
# The diagram has a quadruple bond or double bond?
# Wait! A_2^(2) is rank 1 (finite) + 1 = 2 nodes.
# The Cartan matrix of A_2^(2) has a_00 = 2, a_11 = 2.
# What are a_01 and a_10?
# They must have a_01 * a_10 = 4!
# Either:
# Case A: a_01 = -1, a_10 = -4.
# Case B: a_01 = -4, a_10 = -1.
# Case C: a_01 = -2, a_10 = -2.
# Let's check which one it is!

# In Kac (1990), p. 55:
# For A_2^(2):
# Node 0 has a_0 = 2, a_0^\vee = 1?
# Let's solve:
# If A = [[2, a_01], [a_10, 2]]:
# Null vectors:
# A [a_0, a_1]^T = 0 => 2 a_0 + a_01 a_1 = 0 => a_0/a_1 = -a_01 / 2.
# [a_0^\vee, a_1^\vee] A = 0 => 2 a_0^\vee + a_1^\vee a_10 = 0 => a_0^\vee / a_1^\vee = -a_10 / 2.
# If a_01 = -1, a_10 = -4:
# a_0 / a_1 = 1/2 => a_0 = 1, a_1 = 2.
# a_0^\vee / a_1^\vee = 4/2 = 2 => a_0^\vee = 2, a_1^\vee = 1.
#
# If a_01 = -4, a_10 = -1:
# a_0 / a_1 = 4/2 = 2 => a_0 = 2, a_1 = 1.
# a_0^\vee / a_1^\vee = 1/2 => a_0^\vee = 1, a_1^\vee = 2.

print("Checking A_2^(2) Cartan matrix options")
