import numpy as np, time, os
from scipy.linalg import null_space
exec(open('f4breather.py').read().split("tB = T - 1; d_ = 1e-7")[0])
t1 = time.time()
tC = 2*T/3; d_ = 1e-7
Sres = d_*Fsc(1j*np.pi*(tC + d_)/T)*Rmatrix(solve_fast(np.exp(mu*1j*np.pi*(tC + d_)/T)))
U_, sv, _ = np.linalg.svd(Sres); rk = int(np.sum(sv > 1e-6*sv[0])); Bc = U_[:, :rk]
print(f"self-fusion residue at t = 2T/3: rank {rk}")
WW = (W[:, None, :] + W[None, :, :]).reshape(N2, 4)
hwi = int(np.argmax(W @ np.array([8, 4, 2, 1.]))); whw = W[hwi]
mask = np.all(np.isclose(WW, whw), axis=1)
cvec = null_space(Bc[~mask, :]); hwC = Bc @ cvec[:, 0]; hwC /= np.linalg.norm(hwC)
hw = np.zeros(n); hw[hwi] = 1
al = np.pi/3
def K_C27(theta):
    """bound state C (constituents at theta -/+ i pi/3) scattering a 27 at rapidity 0: top amplitude"""
    v = np.kron(hwC, hw).reshape(n, n, n)
    v = np.moveaxis(v, [1, 2], [0, 1]); sh = v.shape
    v = Sapp(theta + 1j*al, v.reshape(N2, -1)).reshape(sh); v = np.moveaxis(v, [0, 1], [1, 2])
    sh = v.shape; v = Sapp(theta - 1j*al, v.reshape(N2, -1)).reshape(sh)
    tgt = np.kron(hw, hwC); return np.vdot(tgt, v.ravel())/np.vdot(tgt, tgt)
ths = [0.37+0.21j, -0.8+1.3j, 1.1-0.4j, 0.2+2.2j, -0.5-0.3j]
rat = np.array([K_C27(th)/Fsc(th) for th in ths])
print("bootstrap test  S_{C,27}(theta) / S_{27,27}(theta)  (top amplitudes):", np.round(rat, 5), f"[{time.time()-t1:.0f}s]")
np.save("f4boot_ratio.npy", np.array([ths, rat], dtype=object), allow_pickle=True)
