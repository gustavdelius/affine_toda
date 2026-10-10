"""Spin-1 R-matrix of U_q(sl2 hat) and the tensor-product-graph rule (chapter 8).

Checks the q-Serre relation for the spin-1 evaluation representation, solves
Jimbo's equations on spin 1 (x) spin 1, and compares the eigenvalues with
eq-tpg-rule: rho_2=1, rho_1=<2>, rho_0=<2><1>, <a>=(1-z q^{2a})/(z-q^{2a})
(exr-spin1). Takes about a minute.
"""
import sympy as sp

q, x, y = sp.symbols('q x y', nonzero=True)
qn = lambda n: (q**n - q**-n)/(q - 1/q)
E, F = sp.zeros(3), sp.zeros(3)
for k in range(2):
    F[k+1, k] = 1
for k in range(1, 3):
    E[k-1, k] = qn(k)*qn(3-k)
K = sp.diag(q**2, 1, q**-2); I = sp.eye(3); kp = sp.kronecker_product
serre = E**3*F - qn(3)*E**2*F*E + qn(3)*E*F*E**2 - F*E**3
print(('PASS ' if sp.simplify(serre) == sp.zeros(3) else 'FAIL ') + 'q-Serre for spin-1 evaluation rep')
g = lambda z: {'e1': E, 'f1': F, 'k1': K, 'e0': z*F, 'f0': E/z, 'k0': K.inv()}


def cop(gx, gy, i):
    e, f, k = 'e'+i, 'f'+i, 'k'+i
    return {e: kp(gx[e], gx[k]) + kp(I, gy[e]), f: kp(gx[f], I) + kp(gx[k].inv(), gy[f]), k: kp(gx[k], gy[k])}


syms = sp.symbols('r0:81'); R = sp.Matrix(9, 9, syms); eqs = []
for i in '01':
    A, B = cop(g(x), g(y), i), cop(g(y), g(x), i)
    for key in A:
        eqs += list(R*A[key] - B[key]*R)
sol = sp.solve(eqs, syms, dict=True)[0]
Rc = sp.simplify((R.subs(sol)/R.subs(sol)[0, 0]).subs(y, 1))
br = lambda a: (1 - x*q**(2*a))/(x - q**(2*a))
ev = Rc.eigenvals()
want = {1: 5, br(2): 3, br(2)*br(1): 1}
ok = all(any(sp.simplify(k - w) == 0 and ev[k] == m for k in ev) for w, m in want.items())
print(('PASS ' if ok else 'FAIL ') + 'spin-1 eigenvalues rho_2=1 (x5), rho_1=<2> (x3), rho_0=<2><1> (x1)')
