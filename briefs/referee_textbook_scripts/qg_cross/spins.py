"""Item C: spins of the non-local charges and the gradation T, in the book's normalisation
(longest roots |alpha|^2 = 2, action 1/2 (d phi)^2 - V, V ~ sum n_j e^{i beta alpha_j.phi}).
Current J_j = exp(-i a_j.phi_L) with a_j.(beta alpha_j) = 8 pi (double pole with the j-th term of V),
a_j || alpha_j  =>  a_j = (4 pi/beta) alpha_j^vee,  h_j = |a_j|^2/8pi = 8 pi/(beta^2 |alpha_j|^2),  s_j = h_j - 1.
u := 4 pi/beta^2, omega := u - 1."""
import sympy as sp, numpy as np
u, w = sp.symbols('u omega', positive=True)
R = sp.Rational
def eps(n, norm2): return [sp.Matrix([sp.sqrt(norm2) if j == i else 0 for j in range(n)]) for i in range(n)]
def data(name, n=None):
    if name == 'a1':
        return {0: sp.Matrix([-sp.sqrt(2)]), 1: sp.Matrix([sp.sqrt(2)])}, [1, 1], 1
    if name == 'c':      # c_n^(1): (e,e)=1/2
        e = eps(n, R(1, 2)); al = {i: e[i-1]-e[i] for i in range(1, n)}; al[n] = 2*e[n-1]; al[0] = -2*e[0]
        return al, [1] + [2]*(n-1) + [1], 1
    if name == 'b':      # b_n^(1): (e,e)=1
        e = eps(n, 1); al = {i: e[i-1]-e[i] for i in range(1, n)}; al[n] = e[n-1]; al[0] = -(e[0]+e[1])
        return al, [1, 1] + [2]*(n-1), 1
    if name == 'a2n-1(2)':   # (e,e)=1/2
        e = eps(n, R(1, 2)); al = {i: e[i-1]-e[i] for i in range(1, n)}; al[n] = 2*e[n-1]; al[0] = -(e[0]+e[1])
        return al, [1, 1] + [2]*(n-2) + [1], 2
    if name == 'd(2)':       # d_{n+1}^(2): (e,e)=1
        e = eps(n, 1); al = {i: e[i-1]-e[i] for i in range(1, n)}; al[n] = e[n-1]; al[0] = -e[0]
        return al, [1]*(n+1), 2
    if name in ('g2', 'd43'):
        a1 = sp.Matrix([sp.sqrt(R(2, 3)), 0])                   # short, |a1|^2 = 2/3
        a2 = sp.Matrix([-sp.sqrt(R(3, 2)), sp.sqrt(R(1, 2))])   # long,  |a2|^2 = 2, a1.a2 = -1
        assert sp.simplify(a1.dot(a2) + 1) == 0
        if name == 'g2': return {1: a1, 2: a2, 0: -(3*a1+2*a2)}, [1, 3, 2], 1
        return {1: a1, 2: a2, 0: -(2*a1+a2)}, [1, 2, 1], 3
    if name in ('f4', 'e6(2)'):   # Bourbaki f4 with long^2 = 2: e_i of norm 1
        E = eps(4, 1)
        al = {1: E[1]-E[2], 2: E[2]-E[3], 3: E[3], 4: R(1, 2)*(E[0]-E[1]-E[2]-E[3])}
        if name == 'f4': al[0] = -(E[0]+E[1]); return al, [1, 2, 3, 4, 2], 1
        al[0] = -E[0]; return al, [1, 1, 2, 3, 2], 2
def analyse(name, n=None):
    al, nj, k = data(name, n); N = len(nj)
    assert all(sp.simplify(x) == 0 for x in sum((nj[j]*al[j] for j in range(N)), sp.zeros(len(al[0]), 1)))
    L2 = [sp.nsimplify(al[j].dot(al[j])) for j in range(N)]
    s = [2*u/L2[j] - 1 for j in range(N)]                         # s_j = 8pi/(beta^2|alpha_j|^2) - 1
    nv = [nj[j]*L2[j] for j in range(N)]; nv = [sp.nsimplify(x/nv[0]) for x in nv]   # labels of g^vee, n_0^vee = 1
    h = sum(nj); hv_book = sum(nj[j]*L2[j]/2 for j in range(N)); hv = sum(nv)
    T2 = sp.expand(sum(nv[j]*s[j] for j in range(N)))
    return dict(nj=nj, L2=L2, s=s, nv=nv, h=h, hv_book=hv_book, hv=hv, k=k,
                T=sp.expand(T2/2), sum_nj_sj=sp.expand(sum(nj[j]*s[j] for j in range(N))/2),
                book=sp.expand((h*u - hv_book)/2), lit=sp.expand((k*h*u - hv)/2))
WU = lambda e: sp.expand(e.subs(u, w + 1))
partVI = {('c', 2): 2*w + R(1, 2), ('c', 3): 3*w + 1, ('c', 4): 4*w + R(3, 2), ('a2n-1(2)', 3): 5*w + 2, ('a2n-1(2)', 4): 7*w + 3,
          ('e6(2)', None): 9*w + 3, ('f4', None): 6*w + R(3, 2), ('g2', None): 3*w + 1, ('d43', None): 6*w + 3, ('a1', None): w,
          ('b', 3): None, ('d(2)', 3): None}
print("u = 4 pi/beta^2, omega = u - 1; T_spin = (1/2) sum_j n^vee_j s_j ; book formula (4 pi h/beta^2 - h^vee_book)/2;"
      " lit formula (4 pi k h/beta^2 - h^vee_Kac)/2")
for (name, n), TVI in partVI.items():
    d = analyse(name, n)
    print(f"{name}{'' if n is None else '_'+str(n)}: |alpha_j|^2={d['L2']}, s_j={[WU(x) for x in d['s']]}, n_j={d['nj']}, n_j^vee={d['nv']}, h={d['h']}, "
          f"h^vee_book={d['hv_book']}, h^vee_Kac={d['hv']}")
    print(f"     T_spin = {WU(d['T'])}   PartVI T = {TVI}   match: {TVI is None or sp.simplify(WU(d['T']) - TVI) == 0};"
          f"  lit formula = {WU(d['lit'])};  book formula (book beta) = {WU(d['book'])};"
          f"  (1/2) sum n_j s_j = {WU(d['sum_nj_sj'])}")
