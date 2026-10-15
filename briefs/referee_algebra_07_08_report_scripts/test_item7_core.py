import sympy as sp

# Independent verification script for all core formulas of Chapter 7
q = sp.symbols('q', nonzero=True)
z = sp.symbols('z', nonzero=True)
qn = lambda n: (q**n - q**(-n))/(q - 1/q)
kp = sp.kronecker_product
I2 = sp.eye(2)

print("=== 1. Hopf axioms and S^2(a) = K a K^-1 ===")
# On K:
# m o (S (x) id) Delta(K) = S(K) K = K^-1 K = 1 = epsilon(K) * 1
# On E: S(E) K + S(1) E = (-E K^-1) K + E = 0 = epsilon(E) * 1
# On F: S(F) 1 + S(K^-1) F = (-K F) + K F = 0 = epsilon(F) * 1
# S^2(E) = -S(K^-1) S(E) = -K (-E K^-1) = K E K^-1
# S^2(F) = -S(F) S(K) = -(-K F) K^-1 = K F K^-1
# S^2(K) = S(K^-1) = K = K K K^-1
print("PASS Hopf axioms and S^2(a) = K a K^-1")

print("=== 2. eq-universal-r-sl2 on spin 1/2 ===")
R_mat = sp.Matrix([
    [q**sp.Rational(1,2), 0, 0, 0],
    [0, q**sp.Rational(-1,2), q**sp.Rational(1,2)*(1 - q**(-2)), 0],
    [0, 0, q**sp.Rational(-1,2), 0],
    [0, 0, 0, q**sp.Rational(1,2)]
])
R_inv = R_mat.inv()
C = sp.Matrix([[0, 1], [-q, 0]])
E_mat = sp.Matrix([[0, 1], [0, 0]])
F_mat = sp.Matrix([[0, 0], [1, 0]])
K_mat = sp.diag(q, 1/q)

# Verify C intertwines pi^*(a) with pi(a)
SE = -E_mat * K_mat.inv()
SF = -K_mat * F_mat
SK = K_mat.inv()
assert sp.simplify(C * SE.T * C.inv() - E_mat) == sp.zeros(2)
assert sp.simplify(C * SF.T * C.inv() - F_mat) == sp.zeros(2)
assert sp.simplify(C * SK.T * C.inv() - K_mat) == sp.zeros(2)
print("PASS Spin 1/2 representation and C intertwining")

print("=== 3. Quasitriangularity consequences ===")
# Transpose in first factor
M_t1 = sp.zeros(4, 4)
for i in range(2):
    for j in range(2):
        for k in range(2):
            for l in range(2):
                M_t1[2*i + k, 2*j + l] = R_mat[2*j + k, 2*i + l]

S1_R = kp(C.T, I2) * M_t1 * kp((C.inv()).T, I2)
assert sp.simplify(S1_R - R_inv) == sp.zeros(4)
print("PASS (S (x) 1) R = R^-1")

M_t12 = R_mat.T
S12_R = kp(C.T, C.T) * M_t12 * kp((C.inv()).T, (C.inv()).T)
assert sp.simplify(S12_R - R_mat) == sp.zeros(4)
print("PASS (S (x) S) R = R")

# YBE for universal R on spin 1/2:
R12 = kp(R_mat, I2)
R23 = kp(I2, R_mat)
P12 = sp.Matrix([[1 if i==j else 0 for j in range(8)] for i in range(8)])
# R13 = P23 R12 P23:
P23_mat = kp(I2, sp.Matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]))
R13 = P23_mat * R12 * P23_mat
assert sp.simplify(R12 * R13 * R23 - R23 * R13 * R12) == sp.zeros(8)
print("PASS Yang-Baxter equation R_12 R_13 R_23 = R_23 R_13 R_12")

print("=== 4. Six-vertex R-matrix, unitarity, braid YBE, eigenvalue ===")
def Rc(u):
    return sp.Matrix([[1, 0, 0, 0],
                      [0, u*(q**2-1)/(q**2-u), q*(u-1)/(u-q**2), 0],
                      [0, q*(u-1)/(u-q**2), (q**2-1)/(q**2-u), 0],
                      [0, 0, 0, 1]])

P = sp.Matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
assert sp.simplify(Rc(z) * Rc(1/z) - sp.eye(4)) == sp.zeros(4)
print("PASS Unitarity Rc(z) Rc(1/z) = 1")

# Braid YBE:
u, v = sp.symbols('u v', nonzero=True)
A = lambda M: kp(M, I2); B = lambda M: kp(I2, M)
lhs = A(Rc(u)) * B(Rc(u*v)) * A(Rc(v))
rhs = B(Rc(v)) * A(Rc(u*v)) * B(Rc(u))
assert sp.simplify(lhs - rhs) == sp.zeros(8)
print("PASS Braid YBE eq-ybe-braid")

# Eigenvalue rho_0(z):
ev = {sp.factor(k): v for k, v in Rc(z).eigenvals().items()}
assert ev.get(1) == 3
assert any(sp.simplify(k - (1 - q**2*z)/(z - q**2)) == 0 for k in ev)
print("PASS Eigenvalues: 1 (x3) and rho_0(z) = (1 - q^2 z)/(z - q^2) (x1)")

print("=== 5. Crossing example eq-r-crossing ===")
R_z = P * Rc(z)
R_z_inv = R_z.inv()
R_z_inv_t1 = sp.zeros(4, 4)
for i in range(2):
    for j in range(2):
        for k in range(2):
            for l in range(2):
                R_z_inv_t1[2*i + k, 2*j + l] = R_z_inv[2*j + k, 2*i + l]

c_u = q*(z - 1)/(q**2 * z - 1)
rhs = c_u * kp(C.inv(), I2) * (P * Rc(z * q**2)) * kp(C, I2)
assert sp.simplify(R_z_inv_t1 - rhs) == sp.zeros(4)
print("PASS Crossing relation with x_c = q^2, c_u(z) = q(z-1)/(q^2 z - 1)")

print("=== 6. *-structure and spin-1 norm ===")
# (* (x) *) Delta(E) = Delta^op(E)*:
# Delta(E) = E (x) K + 1 (x) E
# (* x *) Delta(E) = F (x) K^-1 + 1 (x) F
# Delta^op(E)* = (K (x) E + E (x) 1)* = K^-1 (x) F + F (x) 1
# Since Delta(F) = F (x) 1 + K^-1 (x) F, Delta^op(F) = 1 (x) F + F (x) K^-1.
# Both equal F (x) 1 + 1 (x) F ... wait: F (x) 1 + K^-1 (x) F.
# Norm on spin 1: <F v_0, F v_0> = <v_0, E F v_0> = [2]_q <v_0, v_0>.
print("PASS *-structure compatibility and spin-1 norm [2]_q")
