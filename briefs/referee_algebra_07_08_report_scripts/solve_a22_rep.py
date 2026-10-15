import numpy as np
import scipy.linalg

# Solve for e0, f0, k0 in 3-dimensional representation of U_q(a_2^(2))
# Let's write a numerical solver for generic q = 0.83 * exp(0.61j)

def test_a22(conv=1):
    q = 0.83 * np.exp(0.61j)
    # conv 1: A = [[2, -1], [-4, 2]], so a01 = -1, a10 = -4.
    # conv 2: A = [[2, -4], [-1, 2]], so a01 = -4, a10 = -1.
    if conv == 1:
        # Long node is 0, short node is 1.
        # In Kac normalization: alpha_0^2 = 4, alpha_1^2 = 1.
        # q0 = q^2, q1 = q^(1/2) (or q0=q, q1=q^(1/4))
        # Let's say q_i = q^(alpha_i^2/2)
        # If alpha_1^2 = 1, alpha_0^2 = 4:
        # q0 = q^2, q1 = q^(1/2).
        # Check: q0^(a01) = (q^2)^(-1) = q^-2.
        # q1^(a10) = (q^(1/2))^(-4) = q^-2.
        # Notice (alpha_0, alpha_1) = -2, so q^(alpha_0, alpha_1) = q^-2!
        A = np.array([[2, -1], [-4, 2]])
        q0 = q**2
        q1 = q**0.5
    else:
        # Node 0 is short, node 1 is long (Kac convention)
        # alpha_0^2 = 1, alpha_1^2 = 4.
        A = np.array([[2, -4], [-1, 2]])
        q0 = q**0.5
        q1 = q**2

    print(f"Testing conv {conv} with A=\n{A}")

test_a22(1)
test_a22(2)
