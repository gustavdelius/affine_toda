"""Independent check of VJS (arXiv:1601.01559) regime-I mass ratios for a_{2n-1}^(2), general n.
Uses the bare thermodynamic Bethe equations (37) of VJS in Fourier space:
   rho + rho_h = s + K rho,  s = (sinh(w g/2)/sinh(w pi/2), 0, ..., 0).
The hole source term is s_phys = (1-K)^{-1} s; masses ~ residues at the pole nearest the real axis."""
import mpmath as mp
mp.mp.dps = 30
def kernel(n, g, w):
    S = mp.sinh
    K = mp.matrix(n, n)
    for j in range(n-1):
        K[j, j] = S(w*(mp.pi/2 - g))/S(w*mp.pi/2)
    K[n-1, n-1] = S(w*(mp.pi/4 - g))/S(w*mp.pi/4)
    for j in range(n-2):
        K[j, j+1] = K[j+1, j] = S(w*g/2)/S(w*mp.pi/2)
    K[n-2, n-1] = K[n-1, n-2] = S(w*g/2)/S(w*mp.pi/4)
    s = mp.matrix(n, 1); s[0] = S(w*g/2)/S(w*mp.pi/2)
    return K, s
def sphys(n, g, w):
    K, s = kernel(n, g, w)
    return mp.lu_solve(mp.eye(n) - K, s)
def check_n2(g, w):
    v = sphys(2, g, w)
    a = mp.cosh(w*(mp.pi/4 - g/2))/mp.cosh(w*(3*mp.pi/4 - g)); b = 1/(2*mp.cosh(w*(3*mp.pi/4 - g)))
    return v[0] - a, v[1] - b
print("n=2 check of dressed source against VJS (58):", [mp.nstr(abs(x), 3) for x in check_n2(mp.mpf('0.7'), mp.mpf('0.37'))])
def nearest_pole(n, g):
    # det(1-K) has zeros on the imaginary axis; scan w = i y
    Hc = 2*n - 1
    y_pred = 2*mp.pi/(Hc*mp.pi - (Hc+1)*g)
    f = lambda y: 1/mp.norm(sphys(n, g, 1j*y))   # vanishes at a pole of s_phys
    # locate by scanning |s_phys| along imaginary axis
    ys = [mp.mpf(k)/400*y_pred*1.5 for k in range(1, 600)]
    vals = [mp.norm(sphys(n, g, 1j*y)) for y in ys]
    k = max(range(len(vals)), key=lambda i: vals[i] if ys[i] < 1.3*y_pred else 0)
    return ys[k], y_pred
def residue_ratios(n, g, y0):
    # residue of each component at w = i y0: lim (w - i y0) s_phys(w)
    eps = mp.mpf('1e-12')
    yp = mp.findroot(lambda y: 1/sphys(n, g, 1j*y)[n-1], y0)
    r = [(sphys(n, g, 1j*(yp + eps))[j]*eps) for j in range(n)]
    return yp, [r[j]/r[n-1] for j in range(n)]
for n in [2, 3, 4, 5, 6]:
    for g in [mp.mpf('0.3'), mp.mpf('0.9'), mp.mpf('1.4')]:
        y0, ypred = nearest_pole(n, g)
        yp, rat = residue_ratios(n, g, ypred)
        H = (2*n - 1) - g/(mp.pi - g)
        pred = [2*mp.sin(a*mp.pi/H) for a in range(1, n)]
        print(f"n={n} gamma={float(g):.2f}: pole y={mp.nstr(yp,12)} (VJS (95): {mp.nstr(ypred,12)}); "
              f"M_a/M_n (a=1..n-1) = {[mp.nstr(mp.re(x),10) for x in rat[:-1]]};  2 sin(a pi/H_VJS) = {[mp.nstr(x,10) for x in pred]}")
