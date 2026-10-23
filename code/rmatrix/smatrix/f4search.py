import numpy as np, itertools, os
exec(open('f4lift3.py').read().split("print(\"searching integer lift shifts")[0])
# CDS block positions as linear functions of H: evaluate the raw (unnormalised) block list at two H values
def cds_raw(a, b, H):
    Bv = (H - 12)/3
    def brk(x): return [(x - 1, 1), (x + 1, 1), (x - 1 + Bv, -1), (x + 1 - Bv, -1)]
    def sq(x): return brk(x) + brk(H - x)
    spec = {(1, 1): [1, H/3+1], (2, 2): [1, H/3-3, H/3+1, H/3+3], (3, 3): [1, 3, H/3-1, H/3+1, H/3+1, H/3+3],
            (4, 4): [1, 3, 5, H/3-3, H/3-1, H/3-1, H/3+1, H/3+1, H/3+1, H/3+3, H/3+3, H/3+5]}[(a, b)]
    return [z for x in spec for (z, s) in sq(x)]
lin = {}
for a in (1, 2, 3, 4):
    x1, x2 = np.array(cds_raw(a, a, 10.0)), np.array(cds_raw(a, a, 11.0))
    beta = x2 - x1; alpha = x1 - 10*beta; lin[a] = list(zip(alpha, beta))
hits = 0
for sh in itertools.product(range(-2, 3), repeat=3):
    for cdd in (0, 1, -1):
        bb = BB_for(sh, cdd); poles_ = [y for y in bb[0] if y > 0]
        for a in (1, 2, 3, 4):
            Hc = set()
            for y in poles_:
                for al, be in lin[a]:
                    for sgn_ in (1, -1):
                        den = sgn_*y - be*T
                        if abs(den) > 1e-9:
                            H = al*T/den
                            if 3 < H < 40: Hc.add(round(H, 10))
            for H in Hc:
                D_ = cds(a, a, H)
                if len(D_[0]) == len(bb[0]) and D_[1] == bb[1] and np.abs(np.array(D_[0]) - np.array(bb[0])).max() < 1e-7:
                    hits += 1; print(f"   T = {T:.3f}, lift shifts {sh}, sign-CDD {cdd:+d}: IDENTICAL to CDS S_{a}{a} at H = {H:.8f}")
print(f"T = {T:.3f}: {hits} exact matches")
