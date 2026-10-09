# Calibration: a_n^(1) (Hollowood) and c_n^(1) (Delius-Grisaru) one-loop soliton masses
import numpy as np
from massformula import dM_channels
pi = np.pi
def Xa(h, a, b):   # roots of the a_{h-1}^(1) factor X_ab(z) = (z-cos(A-B))/(z-cos(A+B))
    A, B = pi*a/h, pi*b/h
    return [(np.cos(A-B)+0j, 1), (np.cos(A+B)+0j, -1)]
def mass(h, b): return 2*np.sin(pi*b/h)   # m = 1
print("a_n^(1):  Delta M_a / m_a   (classical M_a = 2 h m_a / beta^2)")
for h in range(2, 9):
    row = []
    for a in range(1, h):
        ch = [(mass(h, b), Xa(h, a, b)) for b in range(1, h)]
        row.append(dM_channels(ch).real/mass(h, a))
    print(f"  h={h}:", np.round(row, 10), "  -> -h/(2 pi) =", round(-h/(2*pi), 10))
print("\nc_n^(1) from a_{2n-1}^(1): soliton a<n = (a,2n-a), soliton n single; channels b=1..n")
for n in range(2, 7):
    h = 2*n
    Mcl, dM = {}, {}
    for a in range(1, n+1):
        cons = [a, h-a] if a < n else [n]
        ch = []
        for b in range(1, n+1):
            roots = sum((Xa(h, c, b) for c in cons), [])
            ch.append((mass(h, b), roots))
        dM[a] = dM_channels(ch).real
        Mcl[a] = 2*h*sum(mass(h, c) for c in cons)      # times 1/beta^2
    # M_a/M_n = 2 sin(a pi/H), H = 2n + dH beta^2 :  ratio shift r_a = dM_a/Mcl_a - dM_n/Mcl_n (per beta^2)
    out = []
    for a in range(1, n):
        r = dM[a]/Mcl[a] - dM[n]/Mcl[n]
        x = a*pi/h
        dH = -r*h/(x/np.tan(x))    # d log sin(a pi/H) = -x cot x * dH/H
        out.append(dH)
    print(f"  n={n}: dH/beta^2 per a =", np.round(out, 10), "   -1/(4 pi) =", round(-1/(4*pi), 10),
          "  dM_a/M_a^cl(beta^2):", np.round([dM[a]/Mcl[a] for a in range(1, n+1)], 6))
