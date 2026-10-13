import sys, itertools
import dr
from dr import *
t2s = np.geomspace(1e-4, 3, 12)
def dev2(cfg, types, om, t1=12):
    return max(sector_dev(list(types), [t1, t2], cfg[0], cfg[1], om)[0] for t2 in t2s)
for cfg, types, lo, hi in (('UN', 'VV', 0.69, 0.70), ('UN', 'VW', 0.65, 0.68), ('NN', 'VV', 0.535, 0.54), ('NN', 'VW', 0.50, 0.51)):
    for _ in range(5):
        mid = (lo + hi)/2
        if dev2(cfg, types, mid) > 1e-6: hi = mid
        else: lo = mid
    print(cfg, types, f'threshold in ({lo:.4f}, {hi:.4f})'); sys.stdout.flush()
g3 = np.geomspace(0.003, 8, 7)
for cfg in ('UU', 'UN', 'NN'):
    out = []
    for om in (0.3, 0.5, 0.83, 1.43):
        d = max(sector_dev(['V']*3, [g3[a], g3[b], g3[c]], cfg[0], cfg[1], om)[0]
                for a, b, c in itertools.combinations(range(len(g3))[::-1], 3))
        out.append(f'{om}: {d:.1e}')
    print(cfg, 'VVV zeta=1:', '  '.join(out)); sys.stdout.flush()
# control: the other left-wall convention (matrix on t(-theta) = right-wall matrix on t(theta))
orig = dr.transfer
def transfer_swap(j, types, ths, left, right, om):
    old = dr.wall
    def w2(code, typ, o, th, _old=old):
        return _old(code, typ, o, th)
    return orig(j, types, ths, left, right, om)
src = open('dr.py').read().replace("M = op1(wall(left, bar[t], om, -x), p, N) @ M", "M = op1(wall(left, t, om, -x), p, N) @ M")
ns = {}
exec(src, ns)
print('control, swapped left wall: commutator', ns['sector_dev'](['V', 'V'], [0.7, 0.25], 'N', 'N', 0.37)[1],
      ns['sector_dev'](['V', 'V'], [0.7, 0.25], 'U', 'N', 0.37)[1])
