import numpy as np
from boundary_common import q_phys, an_vector, rcheck
from nk_common import fast_rcheck
np.set_printoptions(precision=4, suppress=True, linewidth=200)
om=0.37; q=q_phys(om)
for n in (2,3):
    N=n+1
    y1,y2=1.3+0.2j,0.7-0.1j
    R=fast_rcheck(an_vector(n,q,y1),an_vector(n,q,y2))
    R2,_=rcheck(an_vector(n,q,y1),an_vector(n,q,y2))
    print('agree',np.abs(R-R2).max())
    z=y1/y2
    # blocks
    for i in range(N):
        for j in range(N):
            if i<j:
                a=i*N+j; b=j*N+i
                M=np.array([[R[a,a],R[a,b]],[R[b,a],R[b,b]]])
                # guess: transmission & reflection
                print(n,i,j, 'R[ij->ij]',R[a,a],'R[ij->ji]',R[b,a],'R[ji->ij]',R[a,b],'R[ji->ji]',R[b,b])
    print('diag', [R[i*N+i,i*N+i] for i in range(N)])
    # test candidate closed forms
    # out index convention: R maps A(x)B -> B(x)A; vec index i*N+j for v_i(x)v_j
