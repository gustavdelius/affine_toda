import numpy as np
from boundary_common import q_phys, an_vector
from nk_common import fast_rcheck
om=0.37; q=q_phys(om)
for n in (2,3,4):
    N=n+1
    for (y1,y2) in [(1.3+0.2j,0.7-0.1j),(0.4+1.1j,0.9+0.3j)]:
        z=y1/y2
        R=fast_rcheck(an_vector(n,q,y1),an_vector(n,q,y2))
        t=R[1*N+0,0*N+1]
        # guess: hom normalisation 1 on hw: Rhom: t = (z^N-1)/(q z^N - q^-1) ; r(i<j) = (q-q^-1) z^N/(q z^N-q^-1)
        Z=z**N
        tg=(Z-1)/(q*Z-1/q)
        print(n, 't ratio', t/tg, ' alt', t/((Z-1)/(Z*q-1/q)*q))
        for d in range(1,N):
            i,j=0,d
            r=R[i*N+j,i*N+j]
            for cand,name in [((q-1/q)*z**(N-d)/(q*Z-1/q),'z^(N-d)'),((q-1/q)*z**(-d)*Z/(q*Z-1/q),'z^-d Z'),((q-1/q)*z**(d)/(q*Z-1/q),'z^d')]:
                print('   d',d,name, r/cand)
