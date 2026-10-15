import sympy as sp

# Test if there exists any 2x2 matrix A such that:
# (A \otimes A) \check R(q) (A \otimes A)^{-1} = \check R(-q)
# or (A \otimes A) R(q) (A \otimes A)^{-1} = R(-q)

z, q = sp.symbols('z q')

# Six-vertex R-check(z; q) from eq-six-vertex:
# Basis: 00, 01, 10, 11
Rcheck_q = sp.Matrix([
    [1, 0, 0, 0],
    [0, z*(q**2 - 1)/(q**2 - z), q*(z - 1)/(z - q**2), 0],
    [0, q*(z - 1)/(z - q**2), (q**2 - 1)/(q**2 - z), 0],
    [0, 0, 0, 1]
])

Rcheck_minus_q = sp.Matrix([
    [1, 0, 0, 0],
    [0, z*(q**2 - 1)/(q**2 - z), -q*(z - 1)/(z - q**2), 0],
    [0, -q*(z - 1)/(z - q**2), (q**2 - 1)/(q**2 - z), 0],
    [0, 0, 0, 1]
])

# Let's check eigenvalues of Rcheck_q vs Rcheck_minus_q:
# Trace:
tr_q = sp.simplify(Rcheck_q.trace())
tr_minus_q = sp.simplify(Rcheck_minus_q.trace())
print(f"Trace(Rcheck(q)): {tr_q}")
print(f"Trace(Rcheck(-q)): {tr_minus_q}")
print(f"Are traces equal? {tr_q == tr_minus_q}")

# Eigenvalues:
# Spin 1 eigenvalue is 1 (multiplicity 3)
# Spin 0 eigenvalue is \rho_0(z) = (1 - q^2 z)/(z - q^2).
# Both only depend on q^2!
# So Rcheck(q) and Rcheck(-q) have IDENTICAL eigenvalues!
# But can they be related by A \otimes A?

# Let A = [[a, b], [c, d]]
a, b, c, d = sp.symbols('a b c d')
A = sp.Matrix([[a, b], [c, d]])
A_kron_A = sp.kronecker_product(A, A)

# We want (A \otimes A) Rcheck(q) = Rcheck(-q) (A \otimes A)
diff = sp.simplify(A_kron_A * Rcheck_q - Rcheck_minus_q * A_kron_A)

# Let's solve for a, b, c, d
equations = [sp.simplify(e) for e in diff if e != 0]
sol = sp.solve(equations, [a, b, c, d])
print(f"Solutions for A: {sol}")

# What if A depends on z?
# Let's check for generic numerical z, q:
import numpy as np
z_val = 1.7 + 0.4j
q_val = 0.83 * np.exp(0.61j)

def num_check():
    Rq = np.array(Rcheck_q.subs({z: z_val, q: q_val})).astype(complex)
    Rmq = np.array(Rcheck_minus_q.subs({z: z_val, q: q_val})).astype(complex)
    
    # We want (A \otimes A) Rq - Rmq (A \otimes A) = 0
    # A has 4 entries: a, b, c, d
    # A \otimes A has entries that are quadratic in a, b, c, d.
    # Let's test all diagonal A = diag(a, d):
    # A \otimes A = diag(a^2, ad, da, d^2) = diag(a^2, ad, ad, d^2).
    # Conjugating Rcheck by diag(a^2, ad, ad, d^2):
    # Entry (0,0) -> R00
    # Entry (1,1) -> R11 (since ad / ad = 1)
    # Entry (2,2) -> R22 (since ad / ad = 1)
    # Entry (3,3) -> R33
    # Entry (1,2) -> R12 * (ad / ad) = R12!
    # Entry (2,1) -> R21 * (ad / ad) = R21!
    # SO DIAGONAL A \otimes A DOES NOT CHANGE R12 OR R21 AT ALL!
    print("Diagonal A (x) A leaves Rcheck completely unchanged!")

num_check()
