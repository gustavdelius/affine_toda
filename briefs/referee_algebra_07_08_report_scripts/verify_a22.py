import numpy as np

# Verify U_q(a_2^(2)) on C^3
def qnum(n, q):
    return (q**n - q**(-n)) / (q - q**(-1))

def qfac(n, q):
    prod = 1.0 + 0.0j
    for i in range(1, n+1):
        prod *= qnum(i, q)
    return prod

def qbin(n, k, q):
    return qfac(n, q) / (qfac(k, q) * qfac(n-k, q))

def check_rep(e, f, k, cartan, qi):
    err = 0.0
    nodes = list(e.keys())
    for i in nodes:
        for j in nodes:
            c = e[i] @ f[j] - f[j] @ e[i]
            if i == j:
                ki = k[i]
                ki_inv = np.linalg.inv(ki)
                c_target = (ki - ki_inv) / (qi[i] - 1.0 / qi[i])
                err = max(err, np.max(np.abs(c - c_target)))
            else:
                err = max(err, np.max(np.abs(c)))
            
            # Cartan relations
            k_comm_e = k[i] @ e[j] @ np.linalg.inv(k[i]) - (qi[i]**cartan[i][j]) * e[j]
            k_comm_f = k[i] @ f[j] @ np.linalg.inv(k[i]) - (qi[i]**(-cartan[i][j])) * f[j]
            err = max(err, np.max(np.abs(k_comm_e)), np.max(np.abs(k_comm_f)))

            # q-Serre
            if i != j:
                m = 1 - cartan[i][j]
                for X in (e, f):
                    S = np.zeros_like(e[i])
                    for s in range(m + 1):
                        coeff = ((-1)**s) * qbin(m, s, qi[i])
                        term = coeff * (np.linalg.matrix_power(X[i], m - s) @ X[j] @ np.linalg.matrix_power(X[i], s))
                        S += term
                    err = max(err, np.max(np.abs(S)))
    return err

# Let's search for the representation of U_q(a_2^(2)) on C^3
# Let Q = 0.83 * exp(0.61j)
Q = 0.83 * np.exp(0.61j)

# In Kac convention (node 0 is short, alpha_0^2 = 1, alpha_1^2 = 4 in Kac scale, or 1/2 and 2):
# A = [[2, -4], [-1, 2]] (where 0 is short, 1 is long)
# Or in book folding: node 0 is long, node 1 is short:
# A = [[2, -1], [-4, 2]]
print("Script template ready")
