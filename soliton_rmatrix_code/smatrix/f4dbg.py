import numpy as np, os
src = open('f4breather.py').read().replace("J = 1500", "J = 4000")
src = src.split("ys = []; [ys.extend")[0]
exec(src)
print(f"T = {T:.3f}, zeros {np.round(A, 3)}")
for t0, od in pz: print(f"   t = {t0:8.4f}  order {od:+d}")
ys = []; [ys.extend([y if o > 0 else -y]*abs(o)) for y, o in pz]
def model(yl, th_):
    out = 1.0
    for y in yl: out *= np.sinh(th_/2 + 1j*np.pi*y/(2*T))/np.sinh(th_/2 - 1j*np.pi*y/(2*T))
    return out
for th_ in (0.37+0.21j, -0.8+1.3j, 1.1-0.4j, 0.05+0.3j):
    print(f"   ratio SB/model at {th_}: {SB(th_)/model(ys, th_):.5f}")
