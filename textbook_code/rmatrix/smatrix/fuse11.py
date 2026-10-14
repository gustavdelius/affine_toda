import numpy as np, time
src = open('s33.py').read(); exec(src[:src.index("# ---------------- pole scan")].split("th = 0.31+0.17j\nlhs = P @ Sfull")[0])
# soliton 1 as the bound state 3+3 -> 1 at u1 = pi(2w+1/2)/(3w+1)
u1 = np.pi*(2*omega + 0.5)/(3*omega + 1)
d = 1e-7; Res = d*Sfull(1j*(u1 + d))
U_, sv, _ = np.linalg.svd(Res); r = int(np.sum(sv > 1e-6*sv[0])); B = U_[:, :r]          # W_1 inside V_3 (x) V_3
print(f"soliton-1 multiplet from the 3+3 residue: dimension {r}")
def S11(theta, u=u1):
    v = np.kron(B, B).reshape(8, 8, 8, 8, -1)
    def app(v, M, pos):
        w = np.moveaxis(v, [pos, pos+1], [0, 1]); sh = w.shape
        w = (M @ w.reshape(64, -1)).reshape(sh); return np.moveaxis(w, [0, 1], [pos, pos+1])
    # constituents [th_c - iu/2, th_c + iu/2, th_c' - iu/2, th_c' + iu/2] (image ordering of the residue)
    v = app(v, Sfull(theta + 1j*u), 1); v = app(v, Sfull(theta), 2)
    v = app(v, Sfull(theta), 0);        v = app(v, Sfull(theta - 1j*u), 1)
    img = v.reshape(8**4, -1); BB = np.kron(B, B)
    coef, *_ = np.linalg.lstsq(BB, img, rcond=None)
    return coef, np.abs(BB @ coef - img).max()/max(np.abs(img).max(), 1e-300)
th0 = 0.27 + 0.11j
M, res = S11(th0); print(f"fused S_11 maps W(x)W into W(x)W to {res:.1e}")
U = S11(th0)[0] @ S11(-th0)[0]; print(f"unitarity of fused S_11: {np.abs(U - np.eye(r*r)).max():.1e}")
# eigenvalue structure: should match the vector R-matrix (27+21+7+7+1+1)
ev = np.linalg.eigvals(M); grp = []
for v_ in ev:
    for g_ in grp:
        if abs(g_[0] - v_) < 1e-6*max(1, abs(v_)): g_[1] += 1; break
    else: grp.append([v_, 1])
print("eigenvalue multiplicities of fused S_11:", sorted([m for _, m in grp], reverse=True))
# pole scan in the physical strip
t0 = time.time()
us = np.linspace(2e-4, np.pi - 2e-4, 2401)
nr = np.array([np.linalg.norm(S11(1j*u)[0]) for u in us]); med = np.median(nr)
pk = [us[i] for i in range(1, len(us)-1) if nr[i] > nr[i-1] and nr[i] > nr[i+1] and nr[i] > 20*med]
print(f"scan done [{time.time()-t0:.0f}s]; {len(pk)} poles found")
def rank(u0):
    out = []
    for dd in (1e-6, 1e-8):
        A_ = dd*S11(1j*(u0 + dd))[0]; s_ = np.linalg.svd(A_, compute_uv=False); out.append(int(np.sum(s_ > 1e-3*s_[0])))
    return out
from scipy.optimize import minimize_scalar
print("  u        t=Tu/pi    residue rank")
for p in pk:
    du = np.pi/2400
    rr = minimize_scalar(lambda u: -np.linalg.norm(S11(1j*u)[0]), bounds=(p - du, p + du), method='bounded', options={'xatol': 1e-10})
    print(f"  {rr.x:.6f}  {T*rr.x/np.pi:8.4f}   {rank(rr.x)}")
print("predictions: 1+1->2 at t = w+1/2 =", omega + 0.5, "; half-integer breathers at t = 3w+1/2-k =", [round(3*omega + 0.5 - k, 2) for k in range(8)])
