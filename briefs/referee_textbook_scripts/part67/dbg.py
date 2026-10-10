exec(open('an_cross.py').read().split('for N in (3,4,5):')[0])
N=3;q=np.exp(0.7j);x=1.3
Wb,Es,Fs,Ks,inv=fused(N,q,[x*q,x/q])
P=Wb@np.linalg.pinv(Wb)
print(Wb.shape)
for G in Es+Fs: print(np.linalg.norm((np.eye(N**2)-P)@G@Wb))
