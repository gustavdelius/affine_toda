import sys, time, numpy as np
sys.path.insert(0, sys.argv[1]); from ref0601_a5_lib import *
chi_end=5*np.pi/3
for rt in (1e-8, 1e-10):
    t=time.time(); K,lab=K_out(chi_end+0.1,(1,-1,-1,1,-1,-1),np.cosh(0.6),X=10.0,rtol=rt,atol=rt*1e-3); print(rt,'%.1fs'%(time.time()-t),np.round(K,7).tolist(),flush=True)
