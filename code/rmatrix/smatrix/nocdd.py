import numpy as np, time
src = open('s33.py').read().replace('J = 6000', 'J = 1500')
exec(src[:src.index("# ---------------- pole scan")].split("th = 0.31+0.17j\nlhs = P @ Sfull")[0])
S0 = Sbraid                                  # S_33 = c f R  (no CDD factor)
def fused(B, u):
    r = B.shape[1]; BB = np.kron(B, B)
    def S(theta):
        v = BB.reshape(8, 8, 8, 8, -1)
        def app(v, M, pos):
            w = np.moveaxis(v, [pos, pos+1], [0, 1]); sh = w.shape
            return np.moveaxis((M @ w.reshape(64, -1)).reshape(sh), [0, 1], [pos, pos+1])
        v = app(v, S0(theta + 1j*u), 1); v = app(v, S0(theta), 2); v = app(v, S0(theta), 0); v = app(v, S0(theta - 1j*u), 1)
        coef, *_ = np.linalg.lstsq(BB, v.reshape(8**4, -1), rcond=None); return coef
    return S
def multiplet(u):
    d = 1e-7; U_, sv, _ = np.linalg.svd(d*S0(1j*(u + d))); r = int(np.sum(sv > 1e-6*sv[0])); return U_[:, :r]
def probe(S, t0, lab, D):
    u0 = np.pi*t0/T; ranks = []; norms = []
    for dd in (1e-5, 1e-7):
        M = S(1j*(u0 + dd)); norms.append(np.linalg.norm(M)); sv = np.linalg.svd(dd*M, compute_uv=False); ranks.append(int(np.sum(sv > 1e-3*sv[0])))
    pole = norms[1]/norms[0] > 30
    print(f"   t = {t0:7.4f}  {lab:44s} {'POLE' if pole else 'regular'}  (norm ratio {norms[1]/norms[0]:8.1f})  residue rank {ranks[1] if pole else '-'}")
t_start = time.time()
u1 = np.pi*(2*omega + 0.5)/T; B1 = multiplet(u1)
S11 = fused(B1, u1)
print(f"S_11 from S_33 without CDD (multiplet dim {B1.shape[1]}):")
for t0, lab in [(np.pi/(6 + 2/omega)*T/np.pi, "former spurious pole t = T/H (=1.185)"), (omega + 0.5, "1+1 -> 2 (soliton 2)"), (omega, "1+1 -> 2* (excited soliton 2)"),
                (3*omega + 0.5, "lowest breather p'=1/2 (particle 1)"), (3*omega, "breather p'=1")]:
    probe(S11, t0, lab, 64)
print(f"[{time.time()-t_start:.0f}s]")
u2 = np.pi*omega/T; B2 = multiplet(u2)
S22 = fused(B2, u2)
print(f"S_22 from S_33 without CDD (multiplet dim {B2.shape[1]}):")
for t0, lab in [(2*omega + 0.5, "2+2 -> 2* on the -q^-4 branch (predicted)"), (2*T/3, "exact 2+2 -> 2 at u = 2pi/3 (expected absent)"),
                (2*omega, "2+2 -> (29) on the +q^-4 branch"), (3*omega + 0.5, "lowest breather p'=1/2 (particle 2)"), (np.pi/(6 + 2/omega)*T/np.pi, "t = T/H (check no spurious pole)")]:
    probe(S22, t0, lab, 841)
print(f"[{time.time()-t_start:.0f}s]")
