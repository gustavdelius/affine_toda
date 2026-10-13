import sys, time, numpy as np
sys.path.insert(0, sys.argv[1]); from ref0601_a5_lib import *
chi_end=5*np.pi/3; th=0.6; Om=np.cosh(th)
for rt,at in ((1e-8,1e-10),(1e-9,1e-9)):
    for delta in (0.03, 0.01, 0.003, 0.001):
        t=time.time(); K,lab=K_out(chi_end+delta,(1,-1,-1,1,-1,-1),Om,X=10.0,rtol=rt,atol=at)
        off=max(abs(K[0,1]),abs(K[1,0])); dg=max(abs(K[i,i]-K_bp_diag(1,th)) for i in range(2))
        print('rtol %.0e delta %.0e  %.1fs  offdiag %.1e  |diag-BP| %.1e'%(rt,delta,time.time()-t,off,dg),flush=True)
