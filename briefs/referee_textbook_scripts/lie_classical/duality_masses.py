"""ch.6 line 111: are untwisted non-simply-laced soliton masses proportional to particle masses of
(i) the Lie dual (g^vee)^(1) (ch.10) or (ii) the affine dual hat g^vee (twisted), as ch.6 says?
Solitons of the folded theory = parent solitons summed over species orbits (species = finite nodes; masses
M_a = 2h m_a/beta^2 with m_a the PF vector)."""
import numpy as np, contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    exec(open('fold.py').read())
def ratios(v): v = np.sort(np.asarray(v)); return np.round(v/v[0], 6)
def particle(name, n, typ, N, pat, f0):
    F = fold(typ, N, pat, f0); return ratios(masses(F['aO'], F['lab']))
def solitons(typ, N, pat):
    F = fold(typ, N, pat, True); R = F['R']
    pf = np.abs(np.linalg.eigh(R.A)[1][:, 0])
    orbs = [[i-1 for i in o] for o in F['orbs'] if 0 not in o]       # species orbits (finite nodes)
    return ratios([pf[o].sum() for o in orbs])
rows = [ # untwisted ghat, parent, pattern, Lie dual (g^vee)^(1) data, affine dual data
 ('c_3^(1)', ('A', 5, [1,1,2,2]), ('b_3^(1)', 'D', 4, [1,1,1,2], True), ('d_4^(2)', 'D', 5, [1,1,2,2], False)),
 ('b_3^(1)', ('D', 4, [1,1,1,2]), ('c_3^(1)', 'A', 5, [1,1,2,2], True), ('a_5^(2)', 'D', 6, [1,2,2,2], False)),
 ('b_4^(1)', ('D', 5, [1,1,1,1,2]), ('c_4^(1)', 'A', 7, [1,1,2,2,2], True), ('a_7^(2)', 'D', 8, [1,2,2,2,2], False)),
 ('g_2^(1)', ('D', 4, [1,1,3]), ('g_2^(1)', 'D', 4, [1,1,3], True), ('d_4^(3)', 'E', 6, [1,3,3], False)),
 ('f_4^(1)', ('E', 6, [1,1,1,2,2]), ('f_4^(1)', 'E', 6, [1,1,1,2,2], True), ('e_6^(2)', 'E', 7, [1,1,2,2,2], False)),
]
for g, par, lie, aff in rows:
    s = solitons(*par)
    pl = particle(lie[0], None, *lie[1:]); pa = particle(aff[0], None, *aff[1:])
    print(f'{g}: soliton ratios {s.tolist()} | {lie[0]} particles {pl.tolist()} match={np.allclose(s, pl)}'
          f' | {aff[0]} particles {pa.tolist()} match={len(s)==len(pa) and np.allclose(s, pa)}')
