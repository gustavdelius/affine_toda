from fractions import Fraction as Fr
from collections import Counter
import numpy as np
# block positions y = A*omega + B (exact); T = 3 omega + 1 ; blocks (y): (y)(-y)=1, (y+2T)=(y), (0)=1, (T)=-1
T = (3, Fr(1)); twoT = (6, Fr(2))
def red(y):
    a, b = Fr(y[0]), Fr(y[1])
    # reduce modulo 2T into (-T, T] using a generic omega value for the ordering
    w = 2.37
    while a*w + b > 3*w + 1 + 1e-9: a, b = a - 6, b - 2
    while a*w + b <= -(3*w + 1) + 1e-9: a, b = a + 6, b + 2
    return (a, b)
def normalize(blocks, sign):
    c = Counter(red(y) for y in blocks); out = Counter()
    for y, k in c.items():
        neg = red((-y[0], -y[1]))
        if y == (0, 0): continue
        if y == (3, 1): sign *= (-1)**k; continue           # (T) = -1
        if y == (-3, -1): sign *= (-1)**k; continue
        out[y] += k
    for y in list(out):
        neg = red((-y[0], -y[1]))
        m = min(out[y], out.get(neg, 0))
        if y != neg and m: out[y] -= m; out[neg] -= m
    return sorted(y for y, k in out.items() for _ in range(k) if k > 0), sign
def shift(blocks, sign, s):
    return normalize([(y[0] + s[0], y[1] + s[1]) for y in blocks] + [(y[0] - s[0], y[1] - s[1]) for y in blocks], sign*sign)
w_ = lambda y: float(y[0]*2.37 + y[1])
def show(y): 
    a, b = y; s = (f"{a}w" if a not in (0, 1, -1) else ("w" if a == 1 else ("-w" if a == -1 else ""))) 
    return (s + (f"{'+' if b >= 0 and s else ''}{b}" if b != 0 or not s else "")).replace("+-", "-")
# numerically identified S[B^a, soliton 3] (block positions in t units, omega = 2.37):
# the sign of S[B^3,3] is + at q = e^{-i pi omega} and - at the physical q = -e^{-i pi omega} (QSIGN=-1); all derived amplitudes contain it squared
SB3 = {3: ([(Fr(-5, 2), -1), (Fr(-1, 2), 0), (Fr(1, 2), Fr(1, 2)), (Fr(3, 2), 0), (Fr(3, 2), 1), (Fr(5, 2), Fr(1, 2))], 1 if float(__import__('os').environ.get('QSIGN', '1')) > 0 else -1),
       1: ([(Fr(1, 2), Fr(1, 2)), (Fr(5, 2), Fr(1, 2))], -1),
       2: ([(0, Fr(1, 4)), (1, Fr(3, 4)), (2, Fr(1, 4)), (3, Fr(3, 4))], 1)}
for a, (bl, sg) in SB3.items():
    print(f"check numeric values S[B{a},3]: {[round(w_(y), 3) for y in bl]}")
s_sol = {1: (1, Fr(1, 4)), 2: (Fr(1, 2), 0)}                      # half fusion angles of 3+3 -> 1, 2 (t units)
s_br = {3: (Fr(3, 2), 0), 1: (Fr(3, 2), Fr(1, 4)), 2: (Fr(3, 2), Fr(1, 4))}   # half fusion angles of lowest breathers
SB = {}
for a in (1, 2, 3):
    SB[(a, 3)] = normalize(*SB3[a])
    for b in (1, 2): SB[(a, b)] = shift(*SB3[a], s_sol[b])
print("\nbreather-soliton amplitudes S[B^a, b] = sign x prod (y), y in t-units (w = omega):")
numeric = {(3,1): [1.435, 1.935, 6.175, 6.675], (3,2): [0.5, 2.87, 3.37, 4.74, 5.24, 7.61], (1,1): [-7.175, -0.935, 3.805, 4.305],
           (1,2): [0.5, 2.87, 5.24, 7.61], (2,1): [0.5, 2.87, 5.24, 7.61]}
for (a, b), (bl, sg) in sorted(SB.items()):
    chk = ""
    if (a, b) in numeric: chk = "  [matches direct numerics]" if sorted(round(w_(y), 3) for y in bl) == sorted(numeric[(a, b)]) else "  [MISMATCH with numerics]"
    print(f"  S[B{a},{b}] = {'+' if sg > 0 else '-'} " + " ".join(f"({show(y)})" for y in bl) + chk)
# breather-breather by two routes
mass = {}
H = 6 + 2/2.37; Tn = 3*2.37 + 1
M = {1: 2*np.cos(np.pi/4 + np.pi/(2*H)), 2: 2*np.cos(np.pi/H), 3: 1.0}
pB = {1: 0.5, 2: 0.5, 3: 1}
mB = {a: 2*M[a]*np.sin(np.pi*pB[a]/(2*Tn)) for a in (1, 2, 3)}
higher = {f"B({a})_{p/2:g}": 2*M[a]*np.sin(np.pi*(p/2)/(2*Tn)) for a in (1, 2, 3) for p in range(1, 12)}
print(f"\nbreather-breather amplitudes (lowest breathers = particles; masses m1={mB[1]:.5f}, m2={mB[2]:.5f}, m3={mB[3]:.5f}):")
for a in (1, 2, 3):
    for b in range(a, 4):
        r1 = shift(*SB[(a, b)], s_br[b]); r2 = shift(*SB[(b, a)], s_br[a])
        ok = (r1 == r2)
        poles = [y for y in r1[0] if 0 < w_(y) < Tn]
        idents = []
        for y in poles:
            u = np.pi*w_(y)/Tn; m = np.sqrt(mB[a]**2 + mB[b]**2 + 2*mB[a]*mB[b]*np.cos(u))
            nm = [k for k, v in higher.items() if abs(v - m) < 1e-5]
            idents.append(f"{show(y)}->{nm[0] if nm else round(m, 4)}")
        zeros = [y for y in r1[0] if w_(y) < 0]
        print(f"  S[B{a},B{b}] = {'+' if r1[1] > 0 else '-'} " + " ".join(f"({show(y)})" for y in r1[0]) +
              f"\n        two fusion routes agree: {ok};  physical-strip poles: {', '.join(idents)}")
