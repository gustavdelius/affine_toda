import numpy as np, time
src = open('s33.py').read().replace('J = 6000', 'J = 1500'); exec(src[:src.index("# ---------------- pole scan")].split("th = 0.31+0.17j\nlhs = P @ Sfull")[0])
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
        if abs(g_[0] - v_) < 2e-3*max(1, abs(v_)): g_[1] += 1; break
    else: grp.append([v_, 1])
print("eigenvalue multiplicities of fused S_11:", sorted([m for _, m in grp], reverse=True))
# targeted residue checks + coarse scan for unexpected poles
def rank(u0):
    out = []
    for dd in (1e-5, 1e-7):
        A_ = dd*S11(1j*(u0 + dd))[0]; s_ = np.linalg.svd(A_, compute_uv=False); out.append(int(np.sum(s_ > 1e-3*s_[0])))
    return out
cands = [("1+1 -> 2, -q^-2 branch (predicted bound state)", omega + 0.5),
         ("1+1 -> 2, +q^-2 branch k=0", omega), ("1+1 -> 2, +q^-2 branch k=1", omega - 1)]
cands += [(f"breather p'={k+0.5}", 3*omega + 0.5 - k - 1 + 0.5) for k in range(0, 3)]   # t = 3w+1/2 - k'
cands += [(f"breather p'={k+1}", 3*omega - k) for k in range(0, 3)]
print("targeted checks (t = T u / pi):")
for lab, t_ in cands:
    u_ = np.pi*t_/T
    print(f"   t = {t_:7.4f}  {lab:45s} residue rank {rank(u_)}  |S| near: {np.linalg.norm(S11(1j*(u_ + 1e-6))[0]):.2e}")
t1 = time.time()
us = np.linspace(2e-3, np.pi - 2e-3, 700)
nr = np.array([np.linalg.norm(S11(1j*u)[0]) for u in us]); med = np.median(nr)
pk = [us[i] for i in range(1, len(us)-1) if nr[i] > nr[i-1] and nr[i] > nr[i+1] and nr[i] > 20*med]
print(f"coarse scan [{time.time()-t1:.0f}s]: peaks at t =", [round(T*p/np.pi, 3) for p in pk])
