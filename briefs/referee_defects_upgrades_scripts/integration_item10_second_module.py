"""Item 10 (ch. 43, prp-dpt-doubled): exact symbolic check, h = 3, 4, 5, of
 (a) the sparse inverse of eq-dt-tvbar: M_Vbar = sum_{d<h} (XU)^d A, U^h = q^-1, (1 - XU) M = (1 - X^h/q) A, so
     H M = (1 - X^h/q) with H = A^-1 (1 - XU), H|i,n> = q^-n_i |i,n> - (X/q)|i+1, n - alpha_i>;
 (b) the e-type equation with psi(e_j) = e_{j-1} for Vbar, directly: M (pi(e_j) (x) 1 + pi(k_j) (x) Y) =
     (pi(e_{j-1}) (x) kappa + 1 (x) Y) M (output flipped to Vbar (x) W), with Y = rho~(e_j), kappa = rho~(k_j) of
     eq-dpt-second-module; and its inverse form (pi(e_j)(x)1 + pi(k_j)(x)Y) H = H (pi(e_{j-1})(x)kappa + 1(x)Y);
 (c) the Cartan-element equation M (pi(k_j) (x) kappa) = (pi(k_{j-1}) (x) kappa) M for V and Vbar;
 (d) the same e-type equation for V (eq-dt-tv), and that Y + c_j q^{n_j} shift still solves V but violates
     kappa Y kappa^-1 = q^2 Y;
 (e) module relations of rho~ (Cartan conjugation q^{a_ij}, prod kappa = 1, q-commutation, q-Serre).
Representations from the chapter (eq-dt-tv, eq-dt-tvbar, eq-dpt-second-module) and boundary_common.an_vector's
conventions (V: e_j = y E_{j,j+1}, k_j = q^{H_j}; Vbar: e_j = y E_{j+1,j}, k_j = q^{-H_j}); states |i, n + m> with n
symbolic (z_l = q^{n_l}, prod z_l = 1) and m an integer shift.  No book or referee code is imported."""
import itertools
import sympy as sp

def run(h):
    q, X, t = sp.symbols('q X t', nonzero=True)
    y = X*t                                           # X = y/t, t = e^{omega(theta_d + i pi/2)}
    zs = list(sp.symbols(f'z0:{h-1}', nonzero=True))
    z = zs + [1/sp.prod(zs)]                          # sum_l n_l = 0 (N in the root lattice, sum s = 0)
    E = lambda i: tuple(1 if l == i % h else 0 for l in range(h))
    add = lambda a, b, c=1: tuple(x + c*w for x, w in zip(a, b))
    alpha = lambda j: add(E(j), E(j + 1), -1)         # alpha_j = e_j - e_{j+1}, indices mod h
    qn = lambda m, l: z[l % h]*q**m[l % h]            # q^{(n + m)_l}
    c = t/(q - 1/q)

    def apply(op, vec):
        out = {}
        for (i, m), co in vec.items():
            for (i2, m2), co2 in op(i, m).items():
                out[(i2, m2)] = out.get((i2, m2), 0) + co*co2
        return out

    def diff(v1, v2):
        keys = set(v1) | set(v2)
        return all(sp.simplify(v1.get(k, 0) - v2.get(k, 0)) == 0 for k in keys)

    # Vbar matrix part (eq-dt-tvbar), as an operator on Vbar (x) W (output flipped): |i, n+m> -> sum_d ...
    def M_Vb(i, m):
        res = {}
        for d in range(h):
            co = X**d*sp.prod([qn(m, i + l) for l in range(d + 1)])
            m2 = add(add(m, E(i), -1), E(i + d))
            res[((i + d) % h, m2)] = res.get(((i + d) % h, m2), 0) + co
        return res
    A = lambda i, m: {(i, m): qn(m, i)}
    U = lambda i, m: {((i + 1) % h, add(m, alpha(i), -1)): qn(m, i + 1)}
    H = lambda i, m: {(i, m): 1/qn(m, i), ((i + 1) % h, add(m, alpha(i), -1)): -X/q}
    ok = True
    for i in range(h):
        v = {(i, (0,)*h): 1}
        # U^h = q^-1
        w = v
        for _ in range(h):
            w = apply(U, w)
        ok &= diff(w, {(i, (0,)*h): 1/q})
        # M = sum (XU)^d A
        acc = {}; w = apply(A, v)
        for d in range(h):
            for k, co in w.items():
                acc[k] = acc.get(k, 0) + co
            w = {k: X*co for k, co in apply(U, w).items()}
        ok &= diff(acc, apply(M_Vb, v))
        # H M = (1 - X^h/q)
        ok &= diff(apply(H, apply(M_Vb, v)), {(i, (0,)*h): 1 - X**h/q})
    print(f'  h = {h}: U^h = 1/q, M_Vbar = sum (XU)^d A, H M_Vbar = 1 - X^h/q: {ok}')

    # generators
    def piVb_e(j):
        return lambda i, m: ({((j + 1) % h, m): y} if i == j % h else {})
    def piVb_k(j):
        return lambda i, m: {(i, m): q**(-alpha(j)[i])}
    def piV_e(j):
        return lambda i, m: ({(j % h, m): y} if i == (j + 1) % h else {})
    def piV_k(j):
        return lambda i, m: {(i, m): q**(alpha(j)[i])}
    def Yop(j, cj=0):
        return lambda i, m: {(i, add(m, alpha(j))): c/qn(m, j) + cj*qn(m, j)}
    def kap(j):
        dj = add(alpha(j), alpha(j - 1), -1)
        return lambda i, m: {(i, add(m, dj)): 1}
    def kapinv(j):
        dj = add(alpha(j), alpha(j - 1), -1)
        return lambda i, m: {(i, add(m, dj, -1)): 1}
    def plus(*ops):
        def f(i, m):
            out = {}
            for op in ops:
                for k, co in op(i, m).items():
                    out[k] = out.get(k, 0) + co
            return out
        return f
    def comp(a, b):                                   # a after b
        return lambda i, m: apply(a, b(i, m))
    def M_V(i, m):                                    # eq-dt-tv: q^{-n_i}|N> v_i + X |N + e_i - e_{i-1}> v_{i-1}
        return {(i, m): 1/qn(m, i), ((i - 1) % h, add(add(m, E(i)), E(i - 1), -1)): X}

    okb = okc = okd = True
    for j in range(h):
        Y, K = Yop(j), kap(j)
        for (Mop, pe, pk, name) in ((M_Vb, piVb_e, piVb_k, 'Vbar'), (M_V, piV_e, piV_k, 'V')):
            Xin = plus(pe(j), comp(pk(j), Y))         # pi(e_j)(x)1 + pi(k_j)(x)Y  (they commute: different factors)
            Xout = plus(comp(pe(j - 1), K), Y)        # pi(e_{j-1}) (x) kappa + 1 (x) Y on the flipped output
            for i in range(h):
                v = {(i, (0,)*h): 1}
                good = diff(apply(Mop, apply(Xin, v)), apply(Xout, apply(Mop, v)))
                cart = diff(apply(Mop, apply(comp(pk(j), K), v)), apply(comp(pk(j - 1), K), apply(Mop, v)))
                if name == 'Vbar':
                    inv = diff(apply(Xin, apply(H, v)), apply(H, apply(Xout, v)))
                    okb &= good and inv
                else:
                    okd &= good
                okc &= cart
        # controls: the opposite rotation (pi(e_{j+1}) on the output) and c -> -c fail for Vbar
        Xin_b = plus(piVb_e(j), comp(piVb_k(j), Y))
        bad1 = plus(comp(piVb_e(j + 1), K), Y)
        Ym = lambda i, m, Y=Y: {k: -co for k, co in Y(i, m).items()}
        Xin_m = plus(piVb_e(j), comp(piVb_k(j), Ym)); bad2 = plus(comp(piVb_e(j - 1), K), Ym)
        ctl1 = all(diff(apply(M_Vb, apply(Xin_b, {(i, (0,)*h): 1})), apply(bad1, apply(M_Vb, {(i, (0,)*h): 1}))) for i in range(h))
        ctl2 = all(diff(apply(M_Vb, apply(Xin_m, {(i, (0,)*h): 1})), apply(bad2, apply(M_Vb, {(i, (0,)*h): 1}))) for i in range(h))
        okb &= (not ctl1) and (not ctl2)
        # free constant: still solves V, breaks kappa Y kappa^-1 = q^2 Y
        Yc = Yop(j, cj=sp.Rational(3, 10))
        Xin = plus(piV_e(j), comp(piV_k(j), Yc)); Xout = plus(comp(piV_e(j - 1), K), Yc)
        for i in range(h):
            v = {(i, (0,)*h): 1}
            okd &= diff(apply(M_V, apply(Xin, v)), apply(Xout, apply(M_V, v)))
        v = {(0, (0,)*h): 1}
        okd &= not diff(apply(comp(comp(K, Yc), kapinv(j)), v), {k: q**2*co for k, co in apply(Yc, v).items()})
    print(f'  h = {h}: Vbar e-type equation (direct and with H) for every node, controls (opposite rotation, c -> -c) fail: {okb}; Cartan-element equation, V and Vbar: '
          f'{okc}; V e-type equation, and the free constant solves V but breaks the Cartan relation: {okd}')
    # (e) module relations
    oke = True
    v = {(0, (0,)*h): 1}
    prodk = v
    for j in range(h):
        prodk = apply(kap(j), prodk)
    oke &= diff(prodk, v)
    for i in range(h):
        for j in range(h):
            aij = 2 if i == j else (-1 if (i - j) % h in (1, h - 1) else 0)
            lhs = apply(comp(comp(kap(i), Yop(j)), kapinv(i)), v)
            oke &= diff(lhs, {k: q**aij*co for k, co in apply(Yop(j), v).items()})
            if (j - i) % h == 1:                        # e_i e_{i+1} = q^-1 e_{i+1} e_i
                oke &= diff(apply(comp(Yop(i), Yop(j)), v), {k: co/q for k, co in apply(comp(Yop(j), Yop(i)), v).items()})
            elif i != j and (i - j) % h != 1:
                oke &= diff(apply(comp(Yop(i), Yop(j)), v), apply(comp(Yop(j), Yop(i)), v))
            if i != j and (i - j) % h in (1, h - 1):   # q-Serre e_i^2 e_j - [2] e_i e_j e_i + e_j e_i^2 = 0
                ei, ej = Yop(i), Yop(j)
                t1 = apply(comp(ei, comp(ei, ej)), v); t2 = apply(comp(ei, comp(ej, ei)), v); t3 = apply(comp(ej, comp(ei, ei)), v)
                keys = set(t1) | set(t2) | set(t3)
                oke &= all(sp.simplify(t1.get(k, 0) - (q + 1/q)*t2.get(k, 0) + t3.get(k, 0)) == 0 for k in keys)
    print(f'  h = {h}: rho~ relations (prod kappa = 1, Cartan conjugation q^a_ij, q-commutation, q-Serre): {oke}')
    return ok and okb and okc and okd and oke

allok = True
for h in (3, 4, 5):
    allok &= run(h)
print('PASS' if allok else 'FAIL', ': item 10 - sparse inverse, Vbar half, Cartan equations and module relations, h = 3, 4, 5, exact')
