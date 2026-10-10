import numpy as np, warnings
warnings.filterwarnings('ignore')
exec(open('floating_check.py').read().split("for n in [3,4]:")[0])
e=np.eye(4)
R=[e[1]-e[2], e[2]-e[3], e[3], 0.5*(e[0]-e[1]-e[2]-e[3]), -(e[0]+e[1])]
mk=[2,3,4,2,1]
assert np.allclose(sum(n*a for n,a in zip(mk,R)),0)
ms,d=shifts(R,mk); o=np.argsort(ms); ms=ms[o]; d=d[o]
def cds(H):
    Hp=1/(1/6-1/H)
    return np.array([np.sin(np.pi/H)*np.sin(2*np.pi/Hp), np.sin(3*np.pi/H)*np.sin(np.pi/Hp),
                     np.sin(2*np.pi/H)*np.sin(2*np.pi/Hp), np.sin(3*np.pi/H)*np.sin(2*np.pi/Hp)])
c=cds(12.0); print("classical", np.round(ms/ms[0],6), "CDS at H=12", np.round(np.sort(c)/np.sort(c)[0],6))
cs=np.sort(c); idx=np.argsort(c)  # map sorted
for i in range(1,4):
    r2=ms[i]**2/ms[0]**2; dr2=r2*(d[i]/ms[i]**2-d[0]/ms[0]**2)
    f=lambda H: (np.sort(cds(H))[i]/np.sort(cds(H))[0])**2
    eps=1e-5; dfdH=(f(12+eps)-f(12-eps))/(2*eps)
    print(f"ratio {i+1}/1: dH = {dr2/dfdH:.6f} beta^2 = {4*np.pi*dr2/dfdH:.6f} beta^2/(4pi)")
