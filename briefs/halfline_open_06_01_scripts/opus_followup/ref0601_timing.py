import time, numpy as np, sys
sys.argv=['x']
exec(open('ref0601_a5_K.py').read().split("if __name__ == '__main__':")[0])
chi_end=5*np.pi/3
for delta in (0.1, 0.03, 0.01):
    t=time.time(); K,lab=K_mine(chi_end+delta,(1,-1,-1,1,-1,-1),np.cosh(0.6)); print(delta, 'mine %.1fs'%(time.time()-t), np.round(K,6).tolist(), flush=True)
