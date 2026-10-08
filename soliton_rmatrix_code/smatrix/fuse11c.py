import numpy as np
src = open('fuse11b.py').read(); exec(src[:src.index("# targeted residue checks")].split("ev = np.linalg.eigvals(M)")[0])
def probe(t0, lab):
    u0 = np.pi*t0/T
    from scipy.optimize import minimize_scalar
    rr = minimize_scalar(lambda u: -np.linalg.norm(S11(1j*u)[0]), bounds=(u0 - 0.006, u0 + 0.006), method='bounded', options={'xatol': 1e-11})
    us_ = rr.x; n1 = np.linalg.norm(S11(1j*(us_ + 1e-5))[0]); n2 = np.linalg.norm(S11(1j*(us_ + 1e-6))[0])
    order = np.log10(n2/n1)                       # ~1 simple pole, ~2 double pole
    A_ = 1e-6*S11(1j*(us_ + 1e-6))[0]; sv = np.linalg.svd(A_, compute_uv=False)
    print(f"   t = {T*us_/np.pi:7.4f}  {lab:34s} pole order ~ {order:4.2f}   residue rank {int(np.sum(sv > 1e-3*sv[0]))}")
print("lowest half-integer breather and the unclassified poles of S_11:")
probe(3*omega + 0.5, "breather p'=1/2 (expected rank 1)")
for t_ in (1.187, 1.604, 1.871, 2.612):
    probe(t_, "unclassified")
