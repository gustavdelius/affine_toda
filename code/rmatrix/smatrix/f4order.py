import numpy as np
# exact omega-forms (T = 6w + 3/2): soliton-derived S[B1,B1]
def ours(w):
    T = 6*w + 1.5
    poles = [1, w+1, 2*w+.5, 3*w, 3*w+1.5, 4*w+1, 5*w+.5, 6*w+.5]
    zeros = [4*w+2, 5*w-.5, 3*w+2.5, 3*w-1, w+2, 2*w-.5]
    return T, poles + [-z for z in zeros], -1
def cds11(w):
    T = 6*w + 1.5; H = 2*T/(w + 1); Bv = (H - 12)/3
    def brk(x): return [(x - 1, 1), (x + 1, 1), (x - 1 + Bv, -1), (x + 1 - Bv, -1)]
    ys = []
    for x in (1, H/3 + 1, H - 1, H - H/3 - 1):
        for z, s in brk(x): ys.append(z*T/H if s > 0 else -z*T/H)
    return T, ys, 1
def logS(spec, w, th):
    T, ys, sg = spec(w)
    v = sg*np.prod([np.sinh(th/2 + 1j*np.pi*y/(2*T))/np.sinh(th/2 - 1j*np.pi*y/(2*T)) for y in ys])
    return np.log(v)
ths = [0.3, 0.7, 1.2, 2.0, 3.0]
print("check that the omega-forms reproduce the computed amplitude ordering at w = 2.37: ours T =", ours(2.37)[0])
for w in (200, 400, 800):
    a = np.array([logS(ours, w, th) for th in ths]); b = np.array([logS(cds11, w, th) for th in ths])
    print(f"w = {w:4d}:  w*log S_ours = {np.round(w*a, 5)}\n            w*log S_CDS  = {np.round(w*b, 5)}\n            w*log(Delta) = {np.round(w*(a - b), 5)}   w^2*log(Delta) = {np.round(w*w*(a - b), 4)}")
