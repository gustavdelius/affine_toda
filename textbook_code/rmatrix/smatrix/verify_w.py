import numpy as np, time
from fractions import Fraction as Fr
src = open('breathers.py').read()
exec(src.split("results = {}")[0].split("# superset of real pole positions")[0])
print(f"omega = {omega}, T = {T:.3f}, J = {J}")
# fundamental soliton fusions at their predicted positions
rng2 = np.random.default_rng(3)
def probe(a, b, t0):
    D = dims[a]*dims[b]; X = rng2.normal(size=(D, 40)) + 1j*rng2.normal(size=(D, 40))
    n = [np.linalg.norm(apply_S(a, b, 1j*np.pi*(t0 + d)/T, X)[0]) for d in (1e-5, 1e-7)]
    Y, _ = apply_S(a, b, 1j*np.pi*(t0 + 1e-7)/T, X); sv = np.linalg.svd(1e-7*Y, compute_uv=False)
    return np.log10(n[1]/n[0])/2, int(np.sum(sv > 1e-3*sv[0]))
for (a, b, t0, lab) in [(3, 3, omega, "3+3->2"), (3, 3, 2*omega + 0.5, "3+3->1"), (1, 1, omega + 0.5, "1+1->2"), (1, 1, omega, "1+1->2*"),
                        (3, 1, 2*omega + 0.75, "3+1->3"), (3, 2, 2.5*omega + 1, "3+2->3"), (1, 2, 2.5*omega + 0.75, "1+2->1")]:
    od, rk = probe(a, b, t0); print(f"  {lab:8s} at t = {t0:.4f}: pole order {od:.2f}, residue rank {rk}")
# closed-form breather amplitudes: predicted block positions (A*omega + B)
F = Fr
pred = {
 (3, 3): (+1, [(F(1,2),F(1,2)), (F(3,2),0), (F(3,2),1), (F(5,2),F(1,2)), (F(-1,2),0), (F(-5,2),-1)]),
 (1, 3): (-1, [(F(1,2),F(1,2)), (F(5,2),F(1,2))]),
 (2, 3): (+1, [(0,F(1,4)), (1,F(3,4)), (2,F(1,4)), (3,F(3,4))]),
 (3, 1): (+1, [(F(1,2),F(1,4)), (F(1,2),F(3,4)), (F(5,2),F(1,4)), (F(5,2),F(3,4))]),
 (3, 2): (-1, [(0,F(1,2)), (1,F(1,2)), (1,1), (2,0), (2,F(1,2)), (3,F(1,2))]),
 (1, 1): (+1, [(F(3,2),F(1,4)), (F(3,2),F(3,4)), (F(-1,2),F(1,4)), (F(-5,2),F(-5,4))]),
 (1, 2): (+1, [(0,F(1,2)), (1,F(1,2)), (2,F(1,2)), (3,F(1,2))]),
 (2, 1): (+1, [(0,F(1,2)), (1,F(1,2)), (2,F(1,2)), (3,F(1,2))]),
 (2, 2): (+1, [(F(1,2),F(1,4)), (F(1,2),F(3,4)), (F(3,2),F(1,4)), (F(3,2),F(3,4)), (F(5,2),F(1,4)), (F(5,2),F(3,4)), (F(-1,2),F(1,4)), (F(-5,2),F(-5,4))])}
def model(sign, bl, th):
    out = sign
    for A_, B_ in bl:
        y = float(A_)*omega + float(B_); out *= np.sinh(th/2 + 1j*np.pi*y/(2*T))/np.sinh(th/2 - 1j*np.pi*y/(2*T))
    return out
print("breather-soliton closed forms vs fused numerics (ratio should be 1):")
for (a, b), (sg, bl) in pred.items():
    t1 = time.time()
    r = np.array([K_top(f'B{a}', b, th)[0]/model(sg, bl, th) for th in (0.37+0.21j, -0.8+1.3j, 1.1-0.4j)])
    kind = "input  " if b == 3 else "derived"
    print(f"  S[B{a},{b}] ({kind}): ratio = {np.round(r, 5)}   max |ratio-1| = {np.abs(r - 1).max():.1e}  [{time.time()-t1:.0f}s]")
