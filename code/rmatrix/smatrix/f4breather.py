import numpy as np, time, os
src = open('f4sol.py').read().split("# crossing point and crossing-sign test")[0]
exec(src)
t1 = time.time()
# x-independent pieces of the intertwining system:  e0: Xa = x E0x1 + K0xE0 ; Xb = E0x1 + x K0xE0 ; f0: Xa = (1/x) F0xK0^-1 + 1xF0 ; Xb = F0xK0^-1 + (1/x) 1xF0
E01, KE0 = np.kron(E[0], I), np.kron(Kd[0], E[0]); F0K, IF0 = np.kron(F[0], Ki[0]), np.kron(I, F[0])
pre = []
for ti, t in enumerate(types):
    for a, v in enumerate(t['hw']):
        Rv = np.zeros((N2, len(unk)), complex)
        for b in range(t['mult']): Rv[:, uidx[(ti, b, a)]] += t['hw'][b]
        pre.append((R_apply_cols(E01 @ v), R_apply_cols(KE0 @ v), E01 @ Rv, KE0 @ Rv, R_apply_cols(F0K @ v), R_apply_cols(IF0 @ v), F0K @ Rv, IF0 @ Rv))
def solve_fast(x):
    blocks = []
    for (Le1, Le2, Re1, Re2, Lf1, Lf2, Rf1, Rf2) in pre:
        blocks.append(x*Le1 + Le2 - (Re1 + x*Re2)); blocks.append(Lf1/x + Lf2 - (Rf1 + Rf2/x))
    Mx = np.vstack(blocks); sol, *_ = np.linalg.lstsq(Mx[:, 1:], -Mx[:, 0], rcond=None); return np.concatenate([[1.0], sol])
uchk = solve_R(0.7+0.2j)[0]; print(f"fast solve agrees with full solve: {np.abs(solve_fast(0.7+0.2j) - uchk).max():.0e}")
def Rapply_fast(x, V):
    u = solve_fast(x); C_ = Binv @ V; out = np.zeros(V.shape, complex)
    for ti, t in enumerate(types):
        for a in range(t['mult']):
            s, e_ = offs[ti][a]
            for b in range(t['mult']): out += u[uidx[(ti, b, a)]]*(t['B'][b] @ C_[s:e_])
    return out
def Sapp(theta, V): return Fsc(theta)*Rapply_fast(np.exp(mu*theta), V)
tB = T - 1; d_ = 1e-7
Sres = d_*Fsc(1j*np.pi*(tB + d_)/T)*Rmatrix(solve_fast(np.exp(mu*1j*np.pi*(tB + d_)/T)))
U_, sv, _ = np.linalg.svd(Sres); sig = U_[:, 0]; beta = np.pi*tB/T
print(f"breather residue rank-1 check: sv1/sv0 = {sv[1]/sv[0]:.1e}")
hw = np.zeros(n); hw[np.argmax(W @ np.array([8, 4, 2, 1.]))] = 1
def SB(theta):
    v = np.kron(sig, hw).reshape(n, n, n)
    v = np.moveaxis(v, [1, 2], [0, 1]); sh = v.shape
    v = Sapp(theta + 1j*beta/2, v.reshape(N2, -1)).reshape(sh); v = np.moveaxis(v, [0, 1], [1, 2])
    sh = v.shape; v = Sapp(theta - 1j*beta/2, v.reshape(N2, -1)).reshape(sh)
    tgt = np.kron(hw, sig); return np.vdot(tgt, v.ravel())/np.vdot(tgt, tgt)
cand = set()
for j in (1, 2):
    Ajj = 2*T*(j-1)
    for a in A:
        for k in range(int(3*T)): cand |= {a - k - Ajj, -1 - a - k - Ajj, Ajj + T - a + k, 1 + Ajj + T + a + k}
cands = sorted({round(s*(tp + o), 9) for tp in cand for o in (tB/2, -tB/2) for s in (1, -1)})
cands = [c_ for c_ in cands if 1e-6 < c_ < T - 1e-6]; pz = []
for t0 in cands:
    nn_ = [abs(SB(1j*np.pi*(t0 + d)/T)) for d in (1e-3, 1e-5)]; od = np.log10(nn_[1]/nn_[0])/2
    if abs(od) > 0.5: pz.append((t0, int(round(od))))
ys = []; [ys.extend([y if o > 0 else -y]*abs(o)) for y, o in pz]
def model(yl, th_, sg=1):
    out = sg
    for y in yl: out *= np.sinh(th_/2 + 1j*np.pi*y/(2*T))/np.sinh(th_/2 - 1j*np.pi*y/(2*T))
    return out
r = np.array([SB(th_)/model(ys, th_) for th_ in (0.37+0.21j, -0.8+1.3j, 1.1-0.4j)]); sg = int(np.sign(r.mean().real))
print(f"S[B,27]: {len(ys)} blocks, sign {sg:+d}, fit match {np.abs(r/sg - 1).max():.0e} ({len(cands)} candidates) [{time.time()-t1:.0f}s]")
exec(open('e6sol.py').read().split("def red(y):")[1].split("Hc = 2*T")[0].join(["def red(y):", ""]))
BB = norm_blocks([y + tB/2 for y in ys] + [y - tB/2 for y in ys], sg*sg)
Hc = 2*T/(omega + m1)
for a in (1, 2, 3, 4):
    D_ = cds(a, a, Hc); same = len(D_[0]) == len(BB[0]) and np.allclose(D_[0], BB[0], atol=1e-7) and D_[1] == BB[1]
    print(f"   S[B,B] vs CDS S_{a}{a} at H = 2T/(w+{m1}) = {Hc:.5f}: {'IDENTICAL' if same else 'different'}")
print("   ours:", np.round(BB[0], 3), BB[1])
D1 = cds(1, 1, Hc); print("   CDS S_11:", np.round(D1[0], 3), D1[1])
np.save(f"f4BB_{os.environ.get('AOFF','m1')}.npy", np.array([BB[0], BB[1], T, A], dtype=object), allow_pickle=True)
