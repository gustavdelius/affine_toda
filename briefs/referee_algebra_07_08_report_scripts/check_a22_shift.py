# Analytic calculation of S^2 shift for A_{2l}^{(2)}
# Let's consider l=1 (a_2^(2)), where finite algebra is B_1 = sl_2.
# Simple root of sl_2 is alpha_1.
# What is theta_0?
# In A_2, theta_0 is the odd root.
# There are two choices of which node is node 0:
# Choice 1: Node 0 is the LONG root, node 1 is the SHORT root.
# In the book's folding (06:113):
# alpha_0 has squared length 2 (or 4).
# alpha_1 has squared length 1/2 (or 1).
# Ratio is 4:1.
# If alpha_1^2 = 1, then alpha_0^2 = 4.
# Then rho_bar of B_1:
# <rho_bar, alpha_1^\vee> = 1 => 2(rho_bar, alpha_1)/(alpha_1, alpha_1) = 1 => (rho_bar, alpha_1) = 1/2 * 1 = 1/2.
# What is theta_0? alpha_0 = -theta_0.
# Since alpha_0 = - 2 alpha_1 (in weights: 2(alpha_0, alpha_1)/(alpha_0^2) = -1 => 2(alpha_0, alpha_1)/4 = -1 => (alpha_0, alpha_1) = -2).
# Indeed, if alpha_0 = - 2 alpha_1, then (alpha_0, alpha_0) = 4 (alpha_1, alpha_1) = 4 * 1 = 4.
# Then theta_0 = 2 alpha_1.
# Then 2(rho_bar, theta_0) = 2(rho_bar, 2 alpha_1) = 4 (rho_bar, alpha_1) = 4 * (1/2) = 2!
# Then:
# (alpha_0, alpha_0) + 2(rho_bar, theta_0) = 4 + 2 = 6!
# In this scale, alpha_1^2 = 1, alpha_0^2 = 4.
# Then q_0 = q^(alpha_0^2/2) = q^2, q_1 = q^(alpha_1^2/2) = q^(1/2).
# And the shift is q^6!
# In our test_all_three.py:
# We used q_1 = Q, q_0 = Q^4.
# So Q = q^(1/2).
# Then q^6 = (Q^2)^6 = Q^12!
# EXACT MATCH WITH test_all_three.py WHICH RETURNED Q^12!

print("Analytic check done.")
