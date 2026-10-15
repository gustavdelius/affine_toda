import numpy as np
import itertools

# Let's test the S^2 shift and crossing point numerically and algebraically.
# We will construct representations for:
# 1. U_q(c_2^(1)) vector (dimension 4)
# 2. U_q(d_3^(2)) spinor (dimension 4)
# 3. U_q(a_2^(2)) 3-dimensional representation (fundamental / vector)
#
# Let's carefully define U_q(a_2^(2)):
# a_2^(2) has 2 nodes: 0 and 1.
# What is the Cartan matrix of a_2^(2)?
# In Kac (1990) Table Aff 2:
# A_2^(2) has Cartan matrix:
# [[2, -1],
#  [-4, 2]]  or  [[2, -4], [-1, 2]]?
# Let's check Kac:
# For A_{2l}^{(2)}:
# Node 0 has a_0 = 2, a_0^\vee = 1? Or a_0=2, a_0^\vee = 2?
# Let's check both foldings and definitions.
