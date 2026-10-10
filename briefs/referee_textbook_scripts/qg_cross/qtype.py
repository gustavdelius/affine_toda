"""q -> -q in the six-vertex model equals using the 'type -1' representation (k_i -> -k_i) of the same U_q."""
import sympy as sp
exec(open('qdict.py').read().split('R7s = solve(cop7)')[0])
def gens_m(u):   # K -> -K for both nodes, same q
    KK = -sp.diag(q, 1/q)
    return {'e1': E, 'f1': F, 'k1': KK, 'e0': u*F, 'f0': E/u, 'k0': KK.inv()}
syms = sp.symbols('r0:16'); Rm = sp.Matrix(4, 4, syms); eqs = []
for i in '01':
    A, B = cop7(gens_m(x), gens_m(y), i), cop7(gens_m(y), gens_m(x), i)
    for M1, M2 in zip(A, B): eqs += list(Rm*M1 - M2*Rm)
sol = sp.solve(eqs, syms, dict=True)[0]
Rt = sp.simplify((Rm.subs(sol)/Rm.subs(sol)[0, 0]).subs(y, 1))
print('type(-1) rep of U_q: R = R_7(z;-q):', sp.simplify(Rt - R7(x, -q)) == sp.zeros(4))
