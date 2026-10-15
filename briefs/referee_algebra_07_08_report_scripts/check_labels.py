# Check Kac labels and central element for untwisted, twisted, and a_{2n}^{(2)}
import numpy as np

# Let's inspect root lengths and Kac labels from Kac Table Aff 1, Aff 2, Aff 3
# In Kac (1990), Chapter 4:
# a_i are the numerical marks (null vector of generalized Cartan matrix A: A a = 0? No, a^T A = 0 or A a = 0?)
# Let's verify standard definitions:
# A is the generalized Cartan matrix: A_ij = 2 (\alpha_i, \alpha_j) / (\alpha_i, \alpha_i) = <\alpha_j, \alpha_i^\vee>.
# Null vector of A: \sum_j A_ij a_j = 0, so \sum_j a_j \alpha_j^\vee has <\alpha_i, \sum a_j \alpha_j^\vee> = 0, which is the canonical central element c = \sum a_i^\vee h_i or \sum a_i h_i?
# Let's be very precise:
# In Kac p. 54:
# The null vector a = (a_0, ..., a_l)^T has \sum_{j=0}^l a_ij a_j = 0.
# The dual null vector a^\vee = (a_0^\vee, ..., a_l^\vee) has \sum_{i=0}^l a_i^\vee a_ij = 0.
# The canonical central element is c = \sum_{i=0}^l a_i^\vee \alpha_i^\vee (or in Chevalley basis, \sum a_i^\vee h_i).
# And \delta = \sum_{i=0}^l a_i \alpha_i is the null root!
# Then (\delta, \alpha_j) = (\sum a_i \alpha_i, \alpha_j) = \sum a_i (\alpha_i, \alpha_i)/2 * (2(\alpha_i, \alpha_j)/(\alpha_i, \alpha_i)) = \sum a_i (\alpha_i, \alpha_i)/2 * a_ij.
# Since \sum a_i^\vee a_ij = 0, for this to be 0 for all j, we must have a_i (\alpha_i, \alpha_i)/2 proportional to a_i^\vee!

print("Testing Kac labels relationship")
