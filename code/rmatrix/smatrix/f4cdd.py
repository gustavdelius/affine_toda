import numpy as np, itertools
exec(open('f4particle2.py').read().split("Y = lambda U: U*T/H")[0])
ours = Amp.blocks(BBs, BBsign); c11 = cds_blocks(1, 1, H)
need = c11*Amp(list(ours.den.elements()), list(ours.num.elements()), ours.sign)    # CDS / ours
sig = T - 1
def pair(z, p):              # crossing-symmetric unitary pair (z)(T - z), to power p = +/-1
    ys = [z, T - z] if p > 0 else [-z, -(T - z)]
    return Amp.blocks(ys)
def effect(Z): return Z.shift(sig)*Z*Z*Z.shift(-sig)
grid = sorted({round(a*omega + b, 7) for a in range(0, 7) for b in np.arange(-3, 3.01, 0.5)})
grid = [z for z in grid if 0 < z < T]
hits = []
singles = [(z, p) for z in grid for p in (1, -1)]
for (z, p) in singles:
    Z = pair(z, p)
    if effect(Z) == need: hits.append(((z, p),))
for (z1, p1), (z2, p2) in itertools.combinations(singles, 2):
    Z = pair(z1, p1)*pair(z2, p2)
    if effect(Z) == need: hits.append(((z1, p1), (z2, p2)))
print(f"T = {T:.3f}; CDS/ours needs {sum(need.num.values())} sinh factors up, {sum(need.den.values())} down")
print("CDD factors Z (up to two crossing pairs) reproducing CDS:", hits if hits else "none")
for h_ in hits:
    Z = Amp()
    for (z, p) in h_: Z = Z*pair(z, p)
    boot = Z.shift(T/3)*Z.shift(-T/3) == Z
    print(f"   {h_}: self-fusion bootstrap Z(th) = Z(th+i pi/3) Z(th-i pi/3): {boot}")
