import numpy as np, itertools
from sol_boot import embed, tensor_gen
from boundary_common import q_phys, an_vector, cop
om=0.37; n=2; N=3; q=q_phys(om); E=np.eye(N,dtype=int)
V=an_vector(n,q,1.3)
try:
    X=embed(V,[V],[E[i] for i in range(N)],[E[i] for i in range(N)]); print('identity ok',np.round(X/X[0,0],6))
except AssertionError as e: print('fail',e)
A=an_vector(n,q,1.3); B=an_vector(n,q,0.7)
print('coproduct check', abs(tensor_gen([A,B],'e',1).toarray()-cop(A,B,'e',1)).max(), abs(tensor_gen([A,B],'f',0).toarray()-cop(A,B,'f',0)).max())
