"""Chapter 8 exercise checks.

1. exr-six-vertex-derive: the e_1, f_1, (k_1) and e_0 equations alone already fix
   the six-vertex Rc (no f_0 equation needed), as the exercise asserts.
2. exr-spin1 (fusion part): fuse the six-vertex Rc twice, with the spin-1
   submodule of V(x q) (x) V(x q^-1), and compare eigenvalues with rho_1=<2>, rho_0=<2><1>.
3. exr-faces: face weights for a=d=1/2 via highest-weight vectors in V^{(x)3}.
"""
import sympy as sp
import numpy as np

q, x, y, z = sp.symbols('q x y z', nonzero=True)
E = sp.Matrix([[0, 1], [0, 0]]); F = sp.Matrix([[0, 0], [1, 0]]); K = sp.diag(q, 1/q)
I2 = sp.eye(2); kp = sp.kronecker_product


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


def Rc(u):
    return sp.Matrix([[1, 0, 0, 0],
                      [0, u*(q**2-1)/(q**2-u), q*(u-1)/(u-q**2), 0],
                      [0, q*(u-1)/(u-q**2), (q**2-1)/(q**2-u), 0],
                      [0, 0, 0, 1]])


# ---- 1. which equations are needed
syms = sp.symbols('r0:16'); R = sp.Matrix(4, 4, syms)
g = lambda w: {'e1': E, 'f1': F, 'k1': K, 'e0': w*F, 'f0': E/w, 'k0': K.inv()}
D = lambda gx, gy, key: {'e': kp(gx['e'+key[1]], gx['k'+key[1]]) + kp(I2, gy['e'+key[1]]),
                         'f': kp(gx['f'+key[1]], I2) + kp(gx['k'+key[1]].inv(), gy['f'+key[1]]),
                         'k': kp(gx['k'+key[1]], gy['k'+key[1]])}[key[0]]
eqs = []
for key in ('k1', 'e1', 'f1', 'e0'):
    eqs += list(R*D(g(x), g(y), key) - D(g(y), g(x), key)*R)
sol = sp.solve(eqs, syms, dict=True)
free = set().union(*[R.subs(s).free_symbols for s in sol]) & set(syms)
Rs = R.subs(sol[0])
Rs = sp.simplify((Rs/Rs[0, 0]).subs(y, 1))
report(f'e1,f1,k1,e0 equations: one free parameter ({len(free)}) and solution = eq-six-vertex',
       len(sol) == 1 and len(free) == 1 and sp.simplify(Rs - Rc(x)) == sp.zeros(4))

# ---- 2. fusion to spin 1
qn = complex(0.83, 0.41)**2
Rn = sp.lambdify((z,), Rc(z).subs(q, qn), 'numpy')
Rh = lambda w: np.array(Rn(w), dtype=complex)
I = np.eye(2)
# spin-1 subspace of V(x q)(x)V(x q^-1): image of Rc(q^-2): V(x q^-1)(x)V(x q) -> V(x q)(x)V(x q^-1)
P1 = Rh(qn**-2)
U, s, Vh = np.linalg.svd(P1); B = U[:, :3]          # 4x3 basis of the image
report('Rc(q^-2) has rank 3 (spin-1 image)', s[2] > 1e-8 and s[3] < 1e-10)
xx, yy = complex(0.7, -0.4), complex(-0.3, 1.1)
# step 1: Rc_{1,1/2}(x/w): V1(x)(x)V(w) -> V(w)(x)V1(x), V1(x) in V(xq)(x)V(x/q)
def R1h(x_, w_):
    M = np.kron(Rh(x_*qn/w_), I) @ np.kron(I, Rh(x_/qn/w_))      # on V_a V_b V_d -> V_d V_a V_b
    return M
# step 2: Rc_{1,1}(x/y): V1(x)(x)V1(y), V1(y) in V(yq)(x)V(y/q)
#   V_c (x) V_a' (x) V_b'  --(R_{c a'} (x) 1)-->  V_a' V_c V_b'  --(1 (x) R_{c b'})--> V_a' V_b' V_c
# build explicit operators on (C^2)^{(x)4}: factors (c1 c2)(a')(b')
def op_on(M8, pos):   # M8 acts on 3 consecutive qubits starting at pos, in a 4-qubit space
    left = np.eye(2**pos); right = np.eye(2**(4-pos-3))
    return np.kron(np.kron(left, M8), right)
T = op_on(R1h(xx, yy/qn), 1) @ op_on(R1h(xx, yy*qn), 0)
# T maps V_c(=qubits 0,1) (x) V_a'(2) (x) V_b'(3)  ->  V_a'(0) V_b'(1) V_c(2,3)
emb = np.kron(B, B)                                         # 16x9 : V1(x) (x) V1(y)
img = T @ emb
coef, res, *_ = np.linalg.lstsq(emb, img, rcond=None)
report('fused operator maps V1(x)(x)V1(y) into V1(y)(x)V1(x)', np.allclose(emb @ coef, img, atol=1e-9))
ev = np.linalg.eigvals(coef); zz = xx/yy
# normalise on the top component: eigenvalue of multiplicity 5
vals, counts = np.unique(np.round(ev, 8), return_counts=True)
top = vals[np.argmax(counts)]
ev = ev/top
br = lambda a: (1 - zz*qn**(2*a))/(zz - qn**(2*a))
want = [(1, 5), (br(2), 3), (br(2)*br(1), 1)]
ok = all(sum(abs(e - w) < 1e-7 for e in ev) == m for w, m in want)
report('fused spin-1 R-matrix: eigenvalues 1 (x5), <2> (x3), <2><1> (x1) at z = x/y', ok)

# ---- 3. face weights a=d=1/2
Rz = Rc(z)
D2E = kp(kp(E, K), K) + kp(kp(I2, E), K) + kp(kp(I2, I2), E)
v = lambda *bits: kp(kp(sp.eye(2)[:, bits[0]], sp.eye(2)[:, bits[1]]), sp.eye(2)[:, bits[2]])
DF = kp(F, I2) + kp(K.inv(), F)
# b = 1: spin-1 hw v0v0 in factors 12, lowered once; combine with factor 3 into weight-1/2 hw vector
w1 = kp(sp.Matrix([1, 0, 0, 0]), sp.Matrix([0, 1]))                 # (v0 v0) v1
w2 = kp(DF*sp.Matrix([1, 0, 0, 0]), sp.Matrix([1, 0]))             # (Delta F v0v0) v0
c = sp.symbols('c')
hw1 = w1 + c*w2
csol = [sp.solve(list(D2E*hw1), c, dict=True)[0][c]]
hw1 = sp.simplify(hw1.subs(c, csol[0]))
hw0 = kp(q*sp.Matrix([0, 1, 0, 0]) - sp.Matrix([0, 0, 1, 0]), sp.Matrix([1, 0]))   # singlet(12) (x) v0
report('path vectors are highest-weight', sp.simplify(D2E*hw1) == sp.zeros(8, 1) and sp.simplify(D2E*hw0) == sp.zeros(8, 1))
Op = kp(I2, Rz)
Bm = sp.Matrix.hstack(hw1, hw0)
imgs = Op*Bm
W = sp.simplify((Bm.T*Bm).inv()*Bm.T*imgs)
report('1(x)Rc preserves the span of the two path vectors', sp.simplify(Bm*W - imgs) == sp.zeros(8, 2))
print('face-weight matrix in the basis (b=1, b=0) with these path vectors:')
sp.pprint(sp.factor(W))
print('eigenvalues:', [sp.factor(e) for e in W.eigenvals()])
