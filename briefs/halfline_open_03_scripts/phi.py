import numpy as np
from boundary_common import q_phys, an_vector
from nk_common import fast_rcheck, k_closed
def Phi(n, s, om, th):
    q=q_phys(om); N=n+1
    ym, yp = np.exp(om*(1j*np.pi/2 - th)), np.exp(om*(1j*np.pi/2 + th))
    km = k_closed(n, ym, s, q).T.reshape(N*N); kp = k_closed(n, yp, s, q).T.reshape(N*N)
    R2 = fast_rcheck(an_vector(n, q, np.exp(2*om*th), True), an_vector(n, q, 1.0, True))
    v = R2 @ kp
    c = np.vdot(v, km)/np.vdot(v, v)
    return c, np.linalg.norm(km - c*v)/np.linalg.norm(km)
def U(n,s,om,th):
    q=q_phys(om); y=np.exp(om*th)
    P=k_closed(n,y,s,q)@k_closed(n,1/y,s,q)
    c=P[0,0]; return c, np.abs(P-c*np.eye(n+1)).max()/abs(c)
if __name__=='__main__':
    for n in (1,2,3,4):
        for s in (1,-1):
            for th in (0.31+0.2j, -0.7+0.1j):
                print(n,s,th,Phi(n,s,0.37,th), U(n,s,0.37,th))
