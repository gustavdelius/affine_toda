"""Unitarity of the exact d_n K (polynomial normalization): K(1/y) K(y) = U(y) 1, factored."""
import sys, pickle
import sympy as sp

w, y = sp.symbols('w y')
for n in (4, 5):
    D = pickle.load(open(f'symK_d{n}.pkl', 'rb'))
    for name in ('s+', 's-'):
        K = sp.sympify(D[name]) if False else eval('sp.' + D[name].replace('Symbol(', 'sp.Symbol(').replace('Integer(', 'sp.Integer(')) if False else sp.parse_expr(D[name]) if False else sp.sympify(eval(D[name], sp.__dict__))
        Ki = K.subs(y, 1/y)
        P = (Ki*K).applyfunc(sp.cancel)
        off = max(abs(sp.Poly(sp.numer(sp.together(P[i, j])), w, y).length()) if i != j and P[i, j] != 0 else 0
                  for i in range(2*n) for j in range(2*n))
        diag = set(sp.factor(P[i, i]) for i in range(2*n))
        print(f'd_{n} {name}: off-diagonal nonzero terms: {off}; diagonal entries distinct: {len(diag)}')
        U = sp.factor(P[0, 0])
        print('   U =', U)
        # also print K[0,0], and corner entries for structure
        print('   K[0,0] =', sp.factor(K[0, 0]), '  K[0,2n-1] =', sp.factor(K[0, 2*n - 1]))
