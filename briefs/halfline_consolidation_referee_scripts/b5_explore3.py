import numpy as np
from b5_boundary_krein import broken
grid=np.concatenate([[1e-4,1e-3,1e-2],np.linspace(0.05,3,25),[4,5,6,8,10,12]])
for mus,lw,rw,om in ((('V','W'),'U','N',0.68),(('V','W'),'U','N',0.672),(('V','W'),'N','N',0.52),(('V','W'),'N','N',0.505),(('V','V'),'U','N',0.68),(('V','V'),'N','N',0.52)):
    best=(0,None)
    for t1 in grid:
        for t2 in grid:
            if t1==t2: continue
            b=broken(mus,om,lw,rw,[t1,t2])
            if b>best[0]: best=(b,(t1,t2))
    print(''.join(mus),lw+rw,'omega',om,'max ||Lambda|-1| on the grid',f"{best[0]:.2e}",'at',best[1])
