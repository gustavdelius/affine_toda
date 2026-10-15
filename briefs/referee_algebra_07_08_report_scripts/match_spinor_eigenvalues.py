import numpy as np

# Let's inspect the numerical eigenvalues from solve_jimbo for n=2
# z = 1.4 + 0.3j, Q = 0.83 * exp(0.61j)
z = 1.4 + 0.3j
Q = 0.83 * np.exp(0.61j)

# From our Jimbo solver:
# val0 = 1.0 (Lambda^2)
# val1 = 0.91097437 + 0.56795995j (Lambda^1)
# val2 = 0.41972314 + 0.18354013j (Lambda^0)

# Let's check the ratio:
r1 = 0.910974373679805 + 0.5679599525413364j
r2_over_r1 = (0.41972313659223793 + 0.1835401272421379j) / r1
print(f"r1 = {r1}")
print(f"r2 / r1 = {r2_over_r1}")

# Let's test functions of z and Q:
# What if z is replaced by 1/z?
# What if Q is replaced by Q^-1 or Q^2 or Q^(1/2)?
for q_cand, q_name in [(Q, 'Q'), (Q**2, 'Q^2'), (Q**0.5, 'Q^(1/2)'), (1/Q, '1/Q'), (1/Q**2, '1/Q^2'), (1/Q**0.5, '1/Q^(1/2)')]:
    for z_cand, z_name in [(z, 'z'), (1/z, '1/z')]:
        # Test r1 = <1/2>_- or <1/2>_+
        for a in [0.5, 1.0, 1.5, 2.0]:
            f_minus = (1 - z_cand * q_cand**(2*a)) / (z_cand - q_cand**(2*a))
            f_plus  = (1 + z_cand * q_cand**(2*a)) / (z_cand + q_cand**(2*a))
            inv_minus = (z_cand - q_cand**(2*a)) / (1 - z_cand * q_cand**(2*a))
            inv_plus  = (z_cand + q_cand**(2*a)) / (1 + z_cand * q_cand**(2*a))
            for name, val in [('<a >_-', f_minus), ('<a >_+', f_plus), ('1/<a >_-', inv_minus), ('1/<a >_+', inv_plus)]:
                if abs(val - r1) < 1e-4:
                    print(f"r1 matches {name} at a={a} with z={z_name}, q={q_name}! err={abs(val - r1):.2e}")
                if abs(val - r2_over_r1) < 1e-4:
                    print(f"r2/r1 matches {name} at a={a} with z={z_name}, q={q_name}! err={abs(val - r2_over_r1):.2e}")
