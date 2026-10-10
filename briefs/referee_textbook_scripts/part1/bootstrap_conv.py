"""Test eq-bootstrap  S_dc(th) = S_da(th - i ubar_ac^b) S_db(th + i ubar_bc^a)  (ubar_ac^b = interior angle
opposite m_b) against: (A) the sign-flipped version, (B) the swapped assignment of angles, on e8 and a2,
using Dorey's formula from textbook_code/smatrix_simply_laced.py."""
import itertools, numpy as np
src = open('smatrix_simply_laced.py').read().split("B = 0.37; ths")[0]
exec(src)
B = 0.37; ths = [0.31 + 0.2j, -0.77 + 1.1j]
for typ, n in [('a', 2), ('a', 4), ('e', 6), ('e', 8)]:
    h, bl, mass, orbits = setup(typ, n)
    cnt = {'book': [0, 0], 'signflip': [0, 0], 'swapped': [0, 0]}
    kf = lambda v: tuple(np.round(v, 6) + 0.0)
    orbset = [set(kf(z) for z in orbits[c]) for c in range(n)]
    for a, b, c in itertools.product(range(n), repeat=3):
        if not any(kf(-(x + y)) in orbset[c] for x in orbits[a] for y in orbits[b]): continue
        ma, mb, mc = mass[a], mass[b], mass[c]
        if not (ma + mb > mc + 1e-9 and mb + mc > ma + 1e-9 and ma + mc > mb + 1e-9): continue
        ub_ac = np.arccos((ma**2 + mc**2 - mb**2)/(2*ma*mc))   # opposite m_b
        ub_bc = np.arccos((mb**2 + mc**2 - ma**2)/(2*mb*mc))   # opposite m_a
        for d in range(n):
            for t in ths:
                cb = [k for k in range(n) if kf(-orbits[c][0]) in orbset[k]][0]   # bound state is cbar
                lhs = S(bl[(d, cb)], h, B, t)
                for key, rhs in [('book', S(bl[(d, a)], h, B, t - 1j*ub_ac)*S(bl[(d, b)], h, B, t + 1j*ub_bc)),
                                 ('signflip', S(bl[(d, a)], h, B, t + 1j*ub_ac)*S(bl[(d, b)], h, B, t - 1j*ub_bc)),
                                 ('swapped', S(bl[(d, a)], h, B, t - 1j*ub_bc)*S(bl[(d, b)], h, B, t + 1j*ub_ac))]:
                    cnt[key][0] += 1; cnt[key][1] += abs(lhs - rhs) < 1e-8*max(1, abs(lhs))
    print(typ+str(n), {k: f'{v[1]}/{v[0]} pass' for k, v in cnt.items()})
