import qsign_patch
from irf_a2 import *
rng=np.random.default_rng(3)
for (p,pp) in [(10,13),(9,11),(11,13),(11,14),(8,11),(13,17)]:
    f=Fast(A2RSOS(p,pp)); out=[]
    for types in [(0,1),(0,0,0),(0,1,0,1)]:
        P,_=f.per(types)
        if len(P)>400: continue
        worst=0; pair=0
        for t in range(8):
            ev=np.linalg.eigvals(f.transfer(types,rng.normal(size=len(types))*1.5))
            worst=max(worst,np.abs(np.abs(ev)-1).max())
            pair=max(pair,max(min(abs(e-1/np.conj(g)) for g in ev) for e in ev))
        out.append(f"{''.join('3' if t==0 else 'b' for t in types)}[{len(P)}]: off {worst:.0e} pairing {pair:.0e}")
    print(f"W3({p},{pp}) r={pp-p} lam={(pp-p)/p:.3f} breathers={'yes' if 3*(pp-p)>p else 'no'}:", "; ".join(out), flush=True)
