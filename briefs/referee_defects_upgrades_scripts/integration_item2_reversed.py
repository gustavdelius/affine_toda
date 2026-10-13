"""Item 2 (ch. 41, eq-dt-reversed): the reversed-orientation two-entry matrix intertwines U_q(b^-) on Vbar (x) W for
h = 3, 4, 5, exactly (sympy, symbolic q-powers of the labels, symbolic q, y, c = mu/y_d), and the only coefficient
relation that is not an identity of q-powers is X/y + c q = c/q, i.e. X = -(q - 1/q) mu y/y_d.
Reversed module (orientation -1 of eq-dt-weyl): f_j|n> = c q^{-n_{j+1}}|n - alpha_j>, k_j|n> = q^{alpha_j.n}|n>.
Vbar (boundary_common.an_vector conventions): f_j vb_{j+1} = vb_j / y, k_j vb_i = q^{-(alpha_j)_i} vb_i.
L = flip o M^rev on Vbar (x) W: L|i,n> = q^{n_i}|i,n> + X q^{n_i - n_{i+1}}|i+1, n - alpha_i>.
Equation: L (f_j (x) k_j^-1 + 1 (x) f_j) = (k_j^-1 (x) f_j + f_j (x) 1) L.  No book or referee code."""
import sympy as sp

def check(h, Xval=None):
    q, y, c = sp.symbols('q y c', nonzero=True)
    X = sp.Symbol('X') if Xval is None else Xval(q, y, c)
    zs = list(sp.symbols(f'z0:{h-1}', nonzero=True)); z = zs + [1/sp.prod(zs)]
    E = lambda i: tuple(1 if l == i % h else 0 for l in range(h))
    add = lambda a, b, s=1: tuple(x + s*w for x, w in zip(a, b))
    alpha = lambda j: add(E(j), E(j + 1), -1)
    qn = lambda m, l: z[l % h]*q**m[l % h]
    def apply(op, vec):
        out = {}
        for (i, m), co in vec.items():
            for k, co2 in op(i, m).items():
                out[k] = out.get(k, 0) + co*co2
        return out
    L = lambda i, m: {(i, m): qn(m, i), ((i + 1) % h, add(m, alpha(i), -1)): X*qn(m, i)/qn(m, i + 1)}
    fW = lambda j: (lambda i, m: {(i, add(m, alpha(j), -1)): c/qn(m, j + 1)})
    kWinv = lambda j: (lambda i, m: {(i, m): 1/sp.prod([qn(m, l)**alpha(j)[l] for l in range(h)])})
    fV = lambda j: (lambda i, m: {(j % h, m): 1/y} if i == (j + 1) % h else {})
    kVinv = lambda j: (lambda i, m: {(i, m): q**(alpha(j)[i])})
    def plus(*ops):
        def f(i, m):
            out = {}
            for op in ops:
                for k, co in op(i, m).items():
                    out[k] = out.get(k, 0) + co
            return out
        return f
    comp = lambda a, b: (lambda i, m: apply(a, b(i, m)))
    residuals = []
    for j in range(h):
        Din = plus(comp(fV(j), kWinv(j)), fW(j))          # f_j (x) k_j^-1 + 1 (x) f_j
        Dout = plus(comp(kVinv(j), fW(j)), fV(j))         # k_j^-1 (x) f_j + f_j (x) 1
        for i in range(h):
            v = {(i, (0,)*h): 1}
            l = apply(L, apply(Din, v)); r = apply(Dout, apply(L, v))
            for k in set(l) | set(r):
                d = sp.simplify(l.get(k, 0) - r.get(k, 0))
                if d != 0:
                    residuals.append((j, i, k[0], sp.factor(d)))
    return residuals

ok = True
for h in (3, 4, 5):
    res = check(h)
    rel = {sp.factor(r[3]) for r in res}
    print(f'h = {h}, symbolic X: nonzero components {len(res)}, distinct residual factors: {rel}')
    res2 = check(h, lambda q, y, c: -(q - 1/q)*c*y)
    print(f'h = {h}, X = -(q - 1/q) c y = -(q - 1/q) mu y/y_d: nonzero components {len(res2)}')
    ok &= len(res2) == 0 and len(res) > 0
print('PASS' if ok else 'FAIL', ': item 2 - eq-dt-reversed intertwines exactly iff X/y + c q = c/q (h = 3, 4, 5)')
