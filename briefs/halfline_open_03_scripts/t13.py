import numpy as np
from nk_common import fast_rcheck
from boundary_common import q_phys, an_vector
om=0.37; n=2; N=3; q=q_phys(om); h=N
for d in (2j*np.pi/h, -2j*np.pi/h):
    y1=np.exp(om*(0.3+d/2)); y2=np.exp(om*(0.3-d/2))
    try:
        R=fast_rcheck(an_vector(n,q,y1),an_vector(n,q,y2)); print(d, np.round(np.linalg.svd(R,compute_uv=False),5))
    except AssertionError as e: print(d,'rcheck fails',e)
