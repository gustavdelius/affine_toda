"""Three solitons VVV, sector zeta = 1: is the spectrum broken for every pair of walls at every coupling?
Random and grid rapidity triples at several omega."""
import numpy as np
from b5_boundary_krein import broken
rng=np.random.default_rng(1)
trip=[tuple(sorted(rng.uniform(0.01,6,3))[::-1]) for _ in range(40)]+[(0.9,0.4,0.1),(2.0,0.7,0.05),(5.0,1.1,0.3),(10,5,1e-3),(12,0.5,1e-4)]
for lw,rw in (('U','U'),('U','N'),('N','N')):
    row=[]
    for om in (0.1,0.2,0.35,0.5,0.9,1.3):
        b=max(broken(('V','V','V'),om,lw,rw,list(t)) for t in trip)
        row.append(b)
    print(lw+rw,'max ||Lambda|-1| at omega 0.1,0.2,0.35,0.5,0.9,1.3:',' '.join(f"{x:.2e}" for x in row),flush=True)
