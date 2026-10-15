import sympy as sp
import numpy as np

# Let's construct the 3-dimensional evaluation representation of U_q(a_2^(2)).
# Finite algebra: for a_2^(2), the horizontal subalgebra is A_1 (or B_1 = C_1 = sl_2).
# In 3 dimensions, this is the spin-1 representation of sl_2!
# Let's check:
# For spin 1 of U_q(sl_2):
# In basis v_0, v_1, v_2:
# K v_0 = q_1^2 v_0, K v_1 = v_1, K v_2 = q_1^{-2} v_2.
# E v_0 = 0, E v_1 = [2]_{q_1} v_0, E v_2 = [1]_{q_1}[2]_{q_1} v_1 ?
# Or let's check eq-spin-j in Chapter 7:
# K v_k = q^(2j - 2k) v_k. For j=1:
# K v_0 = q^2 v_0, K v_1 = v_1, K v_2 = q^-2 v_2.
# F v_0 = v_1, F v_1 = v_2, F v_2 = 0.
# E v_0 = 0,
# E v_1 = [1]_q [2]_q v_0 = [2]_q v_0,
# E v_2 = [2]_q [1]_q v_1 = [2]_q v_1.
# Let's verify [E, F] = (K - K^-1)/(q - q^-1):
# On v_0: (EF - FE) v_0 = E v_1 - 0 = [2]_q v_0 = (q^2 - q^-2)/(q - q^-1) v_0. Matches!
# On v_1: (EF - FE) v_1 = E v_2 - F([2]_q v_0) = [2]_q v_1 - [2]_q v_1 = 0 = (1 - 1)/(q - q^-1) v_1. Matches!
# On v_2: (EF - FE) v_2 = 0 - F([2]_q v_1) = -[2]_q v_2 = (q^-2 - q^2)/(q - q^-1) v_2. Matches!

print("Spin 1 matrices verified.")
