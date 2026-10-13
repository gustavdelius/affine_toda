import numpy as np, time
from sol_boot import boot_coef
om=0.37
for n in (2,3,4):
    t=time.time()
    th=0.31+0.07j
    cd=boot_coef(n,om,th,0,order=-1); cp=boot_coef(n,om,th,1,order=-1); cm=boot_coef(n,om,th,-1,order=-1)
    print(n,'diag',cd,' s=+1',cp,' s=-1',cm,' P+=',cp[0]/cd[0],' P-=',cm[0]/cd[0], f'{time.time()-t:.1f}s')
