"""ch.6 line 63: symmetries of the affine diagram vs centre. |Aut(affine)| = |Z| |Aut(finite)|; the special nodes (n_j=1)
form one orbit of size |Z|; number of automorphisms moving node 0 compared with |Z|-1."""
import numpy as np, contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    exec(open('fold.py').read())
for typ, n, Z in [('A', 3, 4), ('A', 4, 5), ('D', 4, 4), ('D', 5, 4), ('D', 6, 4), ('E', 6, 3), ('E', 7, 2), ('E', 8, 1),
                  ('B', 3, 2), ('C', 3, 2), ('F', 4, 1), ('G', 2, 1)]:
    R, A, nj = affine(typ, n); C = cartan(A)
    aut = automorphisms(C); fin = [p for p in aut if p[0] == 0]
    orbit0 = sorted({p[0] for p in aut}); special = [j for j in range(n+1) if nj[j] == 1]
    print(f'{typ.lower()}_{n}^(1): |Aut|={len(aut)}, |Aut fixing 0|={len(fin)}, |Z|={Z}, |Aut|/|Aut_0|={len(aut)//len(fin)}, '
          f'orbit of node 0={orbit0}, nodes with n_j=1={special}, #aut moving node 0={len(aut)-len(fin)}')
