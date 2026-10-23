import numpy as np
exec(open('f4R2.py').read().split("def label(z):")[0])
def br(x, l): return (x - q**l)/(1 - x*q**l)
forms = {273: lambda x: br(x, 2), 52: lambda x: br(x, 8), 26: lambda x: br(x, 2)*br(x, 12), 1: lambda x: br(x, 8)*br(x, 18)}
for x in (0.7+0.4j, 1.9-0.6j):
    r, res = rho_at(x)
    print("x =", x, " rho/closed form:", {d: np.round(r[k]/forms[d](x), 8) for k, d in enumerate(dimsC) if d in forms})
# crossing point: V*(x) = V(x t)
def mats(x):
    out = []
    for i in range(5):
        Ei = (x if i == 0 else 1)*E[i]; Fi = F[i]/(x if i == 0 else 1); Ki = Kd[i]
        out.append((Ei, Fi, Ki))
    return out
xa = 0.83+0.29j; found = []
for ph, nm in ((1, '+'), (-1, '-')):
    for l in range(-20, 21):
        tt = ph*q**l; rows = []
        for (Ei, Fi, Ki), (Ej, Fj, Kj) in zip(mats(xa), mats(xa*tt)):
            Kii = np.linalg.inv(Ki)
            for A_, B_ in (((-Kii@Ei).T, Ej), ((-Fi@Ki).T, Fj), (Kii.T, Kj)):
                rows.append(np.kron(np.eye(n), A_.T) - np.kron(B_, np.eye(n)))
        s = np.linalg.svd(np.vstack(rows), compute_uv=False)
        if s[-1] < 1e-9*s[0]: found.append(f"{nm}q^{l}")
print("crossing point: V*(x) = V(x t) for t =", found)
