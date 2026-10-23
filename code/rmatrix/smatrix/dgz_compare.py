from fractions import Fraction as Fr
exec(open('bootstrap_exact.py').read().split("# breather-breather by two routes")[0].split("print(\"\\nbreather-soliton amplitudes")[0])
# DGZ c_n^(1) particle S-matrix:  S_ab = prod_{p=1}^{b} {a-b-1+2p}_H {H-a+b+1-2p}_H ,  H = 6 + B   (a >= b)
# {x} = (x-1)(x+1)/((x-1+B)(x+1-B)) ; a block argument x = c + d*B  (H = 6 + B absorbed into c,d)
# identification: x/H = y/T  with  H = 2T/(omega+1/2),  B = -1/(omega+1/2)  =>  y = x(omega+1/2)/2
def to_y(c, d):            # x = c + d*B  ->  y = A*omega + C
    c, d = Fr(c), Fr(d)
    return (c/2, c/4 - d/2)
def brace(c, d):           # {x} with x = c + d B : returns (numerator blocks, denominator blocks) as (c,d) pairs
    return [(c - 1, d), (c + 1, d)], [(c - 1, d + 1), (c + 1, d - 1)]
def dgz(a, b):
    if a < b: a, b = b, a
    num, den = [], []
    for p in range(1, b + 1):
        for (c, d) in [(a - b - 1 + 2*p, 0), (6 - a + b + 1 - 2*p, 1)]:      # H - a + b + 1 - 2p = (6-a+b+1-2p) + 1*B
            n_, d_ = brace(c, d); num += n_; den += d_
    blocks = [to_y(c, d) for c, d in num] + [tuple(-v for v in to_y(c, d)) for c, d in den]   # 1/(x) = (-x)
    return normalize(blocks, 1)
# our breather-breather amplitudes from the exact bootstrap
SB = {}
for a in (1, 2, 3):
    SB[(a, 3)] = normalize(*SB3[a])
    for b in (1, 2): SB[(a, b)] = shift(*SB3[a], s_sol[b])
ok_all = True
for a in (1, 2, 3):
    for b in range(a, 4):
        ours = shift(*SB[(a, b)], s_br[b]); theirs = dgz(a, b)
        same = ours == theirs; ok_all &= same
        print(f"S[B{a},B{b}] (fused solitons)  vs  DGZ c3 particle S_{a}{b} (continued):  {'IDENTICAL' if same else 'DIFFERENT'}"
              + ("" if same else f"\n    ours {ours}\n    DGZ  {theirs}"))
print("all six identical:", ok_all)
