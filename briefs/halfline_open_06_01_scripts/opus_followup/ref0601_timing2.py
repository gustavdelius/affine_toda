import time, numpy as np, sys
exec(open(sys.argv[1]).read().split("if __name__ == '__main__':")[0])
chi_end=5*np.pi/3
for meth,rt,at,X in (('LSODA',1e-9,1e-12,12.0),('Radau',1e-9,1e-12,12.0),('DOP853',1e-9,1e-12,12.0)):
    for delta in (0.1, 0.03):
        t=time.time()
        try:
            K,lab=K_mine(chi_end+delta,(1,-1,-1,1,-1,-1),np.cosh(0.6),X=X,method=meth,rtol=rt,atol=at)
            print(meth, delta, 'time %.1fs'%(time.time()-t), np.round(K,7).tolist(), flush=True)
        except Exception as e: print(meth, delta, 'ERR', e, flush=True)
