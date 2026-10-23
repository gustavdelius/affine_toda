import numpy as np, itertools
ht = 6
# one-loop (DG) angle shifts delta u / beta^2 for c_3 fusions (from oneloop_test.py)
du = {"11>2": 1/72, "12>1": -1/144, "13>3": 1/144, "22>2": 0.0, "23>3": 1/72, "33>2": -1/36, "33>1": -1/72}
l  = {"11>2": 2, "12>1": 5, "13>3": 4, "22>2": 4, "23>3": 5, "33>2": 2, "33>1": 4}
# R-matrix phase data (exponent sigma of x* = e^{i pi sigma} q^{-l}) in our conventions:
#   R_11: +-q^-2 ; R_33: +q^-2, -q^-4 ; R_22: +q^-4 (only + and +i scanned) ; R_13: +-i q^-4 ;
#   R_12 (W built on the + branch): +-q^-5 ; R_23 (same W): +-i q^-5
# unknowns: kappa, k1 (crossing branch; sigma_c = 0 for n=3), sig11 in {0,1} (physical 1+1->2 branch),
#           s2 = phase shift of particle 2 induced by sig11 (0 if sig11=0, +-1/2 if sig11=1),
#           tau = relative phase offset spinor/vector (quarter steps), sig22 in {0,1}
def sigmas(sig11, s2, tau, sig22):
    return {"11>2": [sig11], "33>2": [0], "33>1": [1], "22>2": [sig22],
            "13>3": [tau + 0.5, tau - 0.5],
            "12>1": [0 + s2, 1 + s2],
            "23>3": [tau - s2 + 0.5, tau - s2 - 0.5]}
sols = []
for kappa in np.arange(0.25, 12.01, 0.25):
    for k1 in range(-6, 7):
        for sig11 in (0, 1):
            for s2 in ((0,) if sig11 == 0 else (0.5, -0.5)):
                for tau in np.arange(-1, 1, 0.25):
                    for sig22 in (0, 1):
                        S = sigmas(sig11, s2, tau, sig22); ok = True; ks = {}
                        for f in du:
                            need = 4*ht*kappa*du[f] + l[f]*2*k1/ht        # = sigma + 2k
                            found = None
                            for s in S[f]:
                                kk = (need - s)/2
                                if abs(kk - round(kk)) < 1e-9: found = (s, int(round(kk))); break
                            if found is None: ok = False; break
                            ks[f] = found
                        if ok: sols.append((kappa, k1, sig11, s2, tau, sig22, ks))
kappas = sorted(set(s[0] for s in sols))
print("kappa values admitting a consistent assignment (<=12):", kappas)
mk = min(kappas); shown = 0
print(f"\nall assignments at the minimal kappa = {mk}:")
for s in sols:
    if s[0] == mk and shown < 6:
        kappa, k1, sig11, s2, tau, sig22, ks = s; shown += 1
        print(f"  k1={k1}, 1+1->2 on the {'-' if sig11 else '+'}q^-2 branch, s2={s2:+.1f}, tau={tau:+.2f}, sig22={sig22}: "
              + ", ".join(f"{f}:(sigma={v[0]:+.1f},k={v[1]})" for f, v in ks.items()))
print(f"  ... total assignments at kappa={mk}: {sum(1 for s in sols if s[0]==mk)}")
# n=2 cross-check of the same model against GMW-calibrated data: kappa must be 1
du2 = {"22>1": -1/32, "12>2": None}
print("\nn=2 check (DG one-loop = GMW exact): 2+2->1 needs B = 4*4*kappa*(-1/32) = -kappa/2; model with sigma=0,k=0,sigma_c=1,k1=0 gives B = -2*1/4 = -1/2  ->  kappa = 1")
