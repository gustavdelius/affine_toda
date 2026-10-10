"""Exact dictionary between the six-vertex conventions of chapter 7 and Part VI.
ch.7 : Delta(E)=E(x)K+1(x)E, Delta(F)=F(x)1+K^-1(x)F ; eval rep e0=xF, f0=E/x, k0=K^-1
VI   : Delta(e)=e(x)1+k(x)e, Delta(f)=f(x)k^-1+1(x)f (rsolve.py, crossing.py) ; same eval rep
"""
import sympy as sp
q, x, y, z, lam, th = sp.symbols('q x y z lambda theta', nonzero=True)
E = sp.Matrix([[0,1],[0,0]]); F = sp.Matrix([[0,0],[1,0]]); K = sp.diag(q, 1/q); I2 = sp.eye(2)
kp = sp.kronecker_product; P = sp.Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
def gens(u, qq=q):
    KK = sp.diag(qq, 1/qq)
    return {'e1': E, 'f1': F, 'k1': KK, 'e0': u*F, 'f0': E/u, 'k0': KK.inv()}
def cop7(a, b, i):
    return [kp(a['e'+i], a['k'+i]) + kp(I2, b['e'+i]), kp(a['f'+i], I2) + kp(a['k'+i].inv(), b['f'+i]), kp(a['k'+i], b['k'+i])]
def copVI(a, b, i):
    return [kp(a['e'+i], I2) + kp(a['k'+i], b['e'+i]), kp(a['f'+i], b['k'+i].inv()) + kp(I2, b['f'+i]), kp(a['k'+i], b['k'+i])]
def solve(cop, qq=q):
    syms = sp.symbols('r0:16'); R = sp.Matrix(4, 4, syms); eqs = []
    for i in '01':
        A, B = cop(gens(x, qq), gens(y, qq), i), cop(gens(y, qq), gens(x, qq), i)
        for M1, M2 in zip(A, B): eqs += list(R*M1 - M2*R)
    sol = sp.solve(eqs, syms, dict=True)[0]
    return sp.simplify((R.subs(sol)/R.subs(sol)[0, 0]).subs(y, 1))
def R7(u, qq):  # eq-six-vertex
    return sp.Matrix([[1,0,0,0],[0,u*(qq**2-1)/(qq**2-u), qq*(u-1)/(u-qq**2),0],[0,qq*(u-1)/(u-qq**2),(qq**2-1)/(qq**2-u),0],[0,0,0,1]])
def ok(name, M): print(('PASS ' if sp.simplify(M) == sp.zeros(*M.shape) else 'FAIL ') + name)
R7s = solve(cop7); RVI = solve(copVI)
ok('ch.7 Jimbo solve = eq-six-vertex', R7s - R7(x, q))
ok('Part VI coproduct is the opposite of ch.7 (Delta_VI = P Delta_7 P), generators e,f,k', sp.Matrix([
    sp.simplify(sum((copVI(gens(x), gens(y), i)[j] - P*cop7(gens(y), gens(x), i)[j]*P for i in '01' for j in range(3)), sp.zeros(4)))]))
ok('R_VI(z;q) = P R_7(1/z;q) P', RVI - P*R7(1/x, q)*P)
ok('R_VI(z;q) = R_7(z;1/q)   [q <-> 1/q]', RVI - R7(x, 1/q))
# q -> -q in eq-six-vertex: only the transmission entry b flips
L = sp.diag(1, 1, -1, 1)
ok('R_7(z;-q) = L R_7(z;q) L^-1, L = diag(1,1,-1,1) (sign on v1(x)v0 only)', R7(x, -q) - L*R7(x, q)*L.inv())
# no factorised (one-particle, rapidity-dependent) diagonal gauge does it
a0, a1, b0, b1 = sp.symbols('a0 a1 b0 b1', nonzero=True)
G = lambda u0, u1: sp.diag(u0, u1)
Rg = kp(G(b0, b1), G(a0, a1))*R7(x, q)*kp(G(a0, a1), G(b0, b1)).inv()  # A(th2)(x)A(th1) R A(th1)^-1(x)A(th2)^-1
print('one-particle gauge leaves b invariant:', sp.simplify(Rg[1, 2] - R7(x, q)[1, 2]) == 0 and sp.simplify(Rg[2, 1] - R7(x, q)[2, 1]) == 0)
a_, b_, cp, cm = 1, R7(x, q)[1, 2], R7(x, q)[1, 1], R7(x, q)[2, 2]
Delta = sp.simplify((a_**2 + b_**2 - cp*cm)/(2*a_*b_))
print('six-vertex invariant Delta = (a^2+b^2-c+c-)/(2ab) =', sp.factor(Delta), ' (odd under q -> -q)')
# sine-Gordon: Zamolodchikov b_Z = sinh(lam th)/sinh(lam(i pi - th)); compare chapter-7 and Part-VI dictionaries
bZ = sp.sinh(lam*th)/sp.sinh(lam*(sp.I*sp.pi - th))
zz = sp.exp(2*lam*th)
for name, R, qq in [('ch.7, q=-e^{i pi lam}', R7, -sp.exp(sp.I*sp.pi*lam)),
                    ('Part VI (R_VI(z;q)=R_7(z;1/q)), q=e^{-i pi w}, w=lam', lambda u, Q: R7(u, 1/Q), sp.exp(-sp.I*sp.pi*lam)),
                    ('Part VI, q=-e^{-i pi w}', lambda u, Q: R7(u, 1/Q), -sp.exp(-sp.I*sp.pi*lam))]:
    b = R(zz, qq)[1, 2]
    r = sp.simplify((b/bZ).rewrite(sp.exp))
    print(f'{name:55s}: b / b_Zamolodchikov = {r}')
# crossing data for q and -q (chapter-7 conventions, already verified for +q in sixvertex.py)
xc = q**2
for qq, C in [(q, sp.Matrix([[0,1],[-q,0]])), (-q, sp.Matrix([[0,1],[q,0]]))]:
    Rm = P*R7(z, qq); lhs = Rm.inv()
    lhs = sp.Matrix(4, 4, lambda r, c: lhs[(c//2)*2 + r % 2, (r//2)*2 + c % 2])   # transpose in factor 1
    rhs = kp(C.inv(), I2)*(P*R7(z*qq**2, qq))*kp(C, I2)
    cu = sp.simplify(lhs[0, 0]/rhs[0, 0])
    print(f'q->{qq}: x_c = q^2, C = {list(C)}, c_u(z) = {sp.factor(cu)},  relation holds: {sp.simplify(lhs - cu*rhs) == sp.zeros(4)}')
