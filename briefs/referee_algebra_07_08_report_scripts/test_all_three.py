import numpy as np

# Comprehensive test of Item 1:
# Representations:
# 1. U_q(a_2^(2)) on C^3
# 2. U_q(d_3^(2)) spinor on C^4
# 3. U_q(c_2^(1)) vector on C^4

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

def check_algebra(rep, cartan, qi, x=1.3+0.2j):
    e, f, k = rep(x)
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
            
            k_comm_e = k[i] @ e[j] @ np.linalg.inv(k[i]) - (qi[i]**cartan[i][j]) * e[j]
            k_comm_f = k[i] @ f[j] @ np.linalg.inv(k[i]) - (qi[i]**(-cartan[i][j])) * f[j]
            err = max(err, np.max(np.abs(k_comm_e)), np.max(np.abs(k_comm_f)))

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

def gens_dual(rep, x, conv, power):
    e, f, k = rep(x)
    out = []
    for i in e:
        ki = np.linalg.inv(k[i])
        if power == 1:
            Se, Sf = (-e[i] @ ki, -k[i] @ f[i]) if conv == '7' else (-ki @ e[i], -f[i] @ k[i])
            out += [Se.T, Sf.T, ki.T]
        else:
            A, Ai = (k[i], ki) if conv == '7' else (ki, k[i])
            out += [A @ e[i] @ Ai, A @ f[i] @ Ai, k[i]]
    return out

def plain(rep, x):
    e, f, k = rep(x)
    out = []
    for i in e:
        out += [e[i], f[i], k[i]]
    return out

def equivalent(A, B):
    D = A[0].shape[0]
    M = np.vstack([np.kron(np.eye(D), a.T) - np.kron(b, np.eye(D)) for a, b in zip(A, B)])
    s = np.linalg.svd(M, compute_uv=False)
    return s[-1] < 1e-8 * s[0]

def scan(rep, conv, power, base_q, roots=(1, -1, 1j, -1j), lmax=24, x0=0.77+0.31j):
    hits = []
    for r in roots:
        for l in range(-lmax, lmax+1):
            t = r * base_q**l
            if equivalent(gens_dual(rep, x0, conv, power), plain(rep, x0 * t)):
                hits.append((r, l))
    return hits

def lab(h, base_name='q'):
    out = []
    for r, l in h:
        ph = np.angle(r) / (2 * np.pi)
        prefix = '' if abs(ph) < 1e-9 else ('-' if abs(abs(ph) - 0.5) < 1e-9 else f'e^(2pi i {ph:.3f})')
        out.append(f"{prefix}{base_name}^{l}")
    return out

# 1. U_q(a_2^(2)) on C^3
# Here q_short = q_1 = Q, q_long = q_0 = Q^4.
# Cartan matrix: A = [[2, -1], [-4, 2]]
def rep_a22(x, Q=Qs):
    c = np.sqrt(qnum(2, Q) + 0j)
    e1 = c * (Emat(3, 0, 1) + Emat(3, 1, 2))
    f1 = c * (Emat(3, 1, 0) + Emat(3, 2, 1))
    k1 = np.diag([Q**2, 1.0, Q**(-2)])
    e0 = x * Emat(3, 2, 0)
    f0 = Emat(3, 0, 2) / x
    k0 = np.diag([Q**(-4), 1.0, Q**4])
    return {1: e1, 0: e0}, {1: f1, 0: f0}, {1: k1, 0: k0}

cart_a22 = {0: {0: 2, 1: -1}, 1: {0: -4, 1: 2}}
qi_a22 = {1: Qs, 0: Qs**4}

# 2. U_q(d_3^(2)) spinor
# Cartan matrix: A = [[2, -2, 0], [-1, 2, -1], [0, -2, 2]]
# q_short = q_0 = q_2 = Q, q_long = q_1 = Q^2
import itertools
def rep_d32(x, Q=Qs):
    n = 2
    W = np.array([[0.5 if b == 0 else -0.5 for b in bits] for bits in itertools.product([0, 1], repeat=n)])
    site = lambda op, i: np.kron(op, np.eye(2)) if i == 0 else np.kron(np.eye(2), op)
    sp_ = Emat(2, 0, 1)
    sm_ = Emat(2, 1, 0)
    e = {1: site(sp_, 0) @ site(sm_, 1), 2: site(sp_, 1), 0: x * site(sm_, 0)}
    f = {1: e[1].T.copy(), 2: e[2].T.copy(), 0: site(sp_, 0) / x}
    k = {1: np.diag(Q**(2*(W[:, 0] - W[:, 1]))), 2: np.diag(Q**(2*W[:, 1])), 0: np.diag(Q**(-2*W[:, 0]))}
    return e, f, k

cart_d32 = {0: {0: 2, 1: -2, 2: 0}, 1: {0: -1, 1: 2, 2: -1}, 2: {0: 0, 1: -2, 2: 2}}
qi_d32 = {0: Qs, 1: Qs**2, 2: Qs}

# 3. U_q(c_2^(1)) vector
# Cartan matrix: A = [[2, -1, 0], [-2, 2, -2], [0, -1, 2]]
# q_long = q_0 = q_2 = Q, q_short = q_1 = Q^(1/2)
def rep_c2(x, Q=Qs):
    P = np.sqrt(Q)
    e = {1: Emat(4, 0, 1) + Emat(4, 2, 3), 2: Emat(4, 1, 2), 0: x * Emat(4, 3, 0)}
    f = {1: Emat(4, 1, 0) + Emat(4, 3, 2), 2: Emat(4, 2, 1), 0: Emat(4, 0, 3) / x}
    k = {1: np.diag([P, 1/P, P, 1/P]), 2: np.diag([1, Q, 1/Q, 1]), 0: np.diag([1/Q, 1, 1, Q])}
    return e, f, k

cart_c2 = {0: {0: 2, 1: -1, 2: 0}, 1: {0: -2, 1: 2, 2: -2}, 2: {0: 0, 1: -1, 2: 2}}
qi_c2 = {0: Qs, 1: np.sqrt(Qs), 2: Qs}

cases = [
    ("U_q(a_2^(2)) fundamental (3-dim)", rep_a22, cart_a22, qi_a22, Qs, "Q"),
    ("U_q(d_3^(2)) spinor (4-dim)", rep_d32, cart_d32, qi_d32, Qs, "Q"),
    ("U_q(c_2^(1)) vector (4-dim)", rep_c2, cart_c2, qi_c2, Qs, "Q"),
]

for name, rep, cart, qi, base_q, base_name in cases:
    res = check_algebra(rep, cart, qi)
    print(f"\n--- {name} ---")
    print(f"Defining relations & Serre residual: {res:.2e}")
    for conv in ('7', 'VI'):
        dd = scan(rep, conv, power=2, base_q=base_q)
        d1 = scan(rep, conv, power=1, base_q=base_q)
        print(f"  Coproduct {conv:2s}: V(x)** = V(x t): t = {lab(dd, base_name)}")
        print(f"               V(x)*  = V(x t): t = {lab(d1, base_name)}")
