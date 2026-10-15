import numpy as np

# Test Kac's convention for A_2^(2) on C^3
# In Kac (1990), node 0 is the short root: (alpha_0, alpha_0) = 1, (alpha_1, alpha_1) = 4.
# A = [[2, -4], [-1, 2]]
# So a_01 = -4, a_10 = -1.
# q_0 = q^(1/2), q_1 = q^2 (where q has long root = 2? No, q is the base parameter).
# Let's say q_i = q^(alpha_i^2 / 2).
# If alpha_0^2 = 1 => q_0 = q^(1/2).
# alpha_1^2 = 4 => q_1 = q^2.
# Then q_0^(a_01) = (q^(1/2))^(-4) = q^(-2).
# q_1^(a_10) = (q^2)^(-1) = q^(-2).
# This is consistent!
# Now in this convention, node 0 has (alpha_0, alpha_0) = 1, a_0 = 2, a_0^\vee = 1.
# Node 1 has (alpha_1, alpha_1) = 4, a_1 = 1, a_1^\vee = 2.
# h = a_0 + a_1 = 3, h^\vee = a_0^\vee + a_1^\vee = 3.
# Let's see: node 1 is the finite algebra sl2 (since a_10 = -1, node 1 is long).
# In C^3, the finite algebra on C^3 has:
# e_1: spin 1 of sl_2, with q_1 = q^2.
# e_0: connects the end weights, carries x.
# Let's build this representation and scan for S^2 shift!

Qs = 0.83 * np.exp(0.61j)

def Emat(N, a, b):
    M = np.zeros((N, N), complex)
    M[a, b] = 1.0
    return M

def qnum(n, q):
    return (q**n - q**(-n)) / (q - 1.0/q)

def qfac(n, q):
    prod = 1.0 + 0.0j
    for i in range(1, n+1):
        prod *= qnum(i, q)
    return prod

def qbin(n, k, q):
    return qfac(n, q) / (qfac(k, q) * qfac(n-k, q))

def rep_a22_kac(x, q=Qs):
    # node 1 is the sl2: q_1 = q^2
    q1 = q**2
    c = np.sqrt(qnum(2, q1) + 0j)
    e1 = c * (Emat(3, 0, 1) + Emat(3, 1, 2))
    f1 = c * (Emat(3, 1, 0) + Emat(3, 2, 1))
    k1 = np.diag([q1**2, 1.0, q1**(-2)])
    
    # node 0 carries x: connects 2 to 0
    # [e0, f0] = (k0 - k0^-1)/(q0 - q0^-1), q0 = q^(1/2)
    # k0 must have k0 e1 k0^-1 = q0^(-4) e1 = q^(-2) e1.
    # So k0 = diag(q^-2, 1, q^2).
    # Then [e0, f0] = diag(-1, 0, 1) * (q^-2 - q^2)/(q^(1/2) - q^(-1/2))
    # Note (q^2 - q^-2)/(q^(1/2) - q^(-1/2)) = [4]_{q^(1/2)}
    q0 = q**0.5
    scale0 = np.sqrt(qnum(4, q0) + 0j)
    e0 = scale0 * x * Emat(3, 2, 0)
    f0 = scale0 * Emat(3, 0, 2) / x
    k0 = np.diag([q**(-2), 1.0, q**2])
    return {1: e1, 0: e0}, {1: f1, 0: f0}, {1: k1, 0: k0}

cart_kac = {0: {0: 2, 1: -4}, 1: {0: -1, 1: 2}}
qi_kac = {0: Qs**0.5, 1: Qs**2}

import test_all_three as T
res = T.check_algebra(rep_a22_kac, cart_kac, qi_kac)
print(f"Kac-convention A_2^(2) defining relations residual: {res:.2e}")

for conv in ('7', 'VI'):
    dd = T.scan(rep_a22_kac, conv, power=2, base_q=Qs)
    d1 = T.scan(rep_a22_kac, conv, power=1, base_q=Qs)
    print(f"  Coproduct {conv:2s}: V(x)** = V(x t): t = {T.lab(dd, 'q')}")
    print(f"               V(x)*  = V(x t): t = {T.lab(d1, 'q')}")
