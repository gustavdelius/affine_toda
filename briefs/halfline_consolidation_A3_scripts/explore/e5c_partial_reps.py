"""a_2 partial Dirichlet from diagram reflections in V+Vbar, adjoint and defining reps, lam = +-0.83.
Log: ../logs/e5c_partial_reps.log (run inline on 2026-10-13 with identical code)."""
src = open('e5_partial.py').read()
ns = {'__name__': 'e5defs'}
exec(compile(src[:src.index("print('=== (a)")], 'e5_partial.py[defs]', 'exec'), ns)
E3, nullity_partial, involutions = ns['E3'], ns['nullity_partial'], ns['involutions']
for name, rep, h in [('a2 V+Vbar', E3.v_plus_vbar(2), 3), ('a2 adjoint', E3.adjoint(2), 3), ('a2 defining', E3.defining(2), 3)]:
    for iname, M in involutions(h):
        r = {s: nullity_partial(rep, M, s, 0.83)[0] for s in ('T', 'S')}
        r2 = {s: nullity_partial(rep, M, s, -0.83)[0] for s in ('T', 'S')}
        print(f'{name}: {iname}: lam=0.83 {r}; lam=-0.83 {r2}')
