import numpy as np, itertools
from sol_boot import embed
from boundary_common import q_phys, an_vector
om=0.37; n=2; N=3; T=N*om/2; q=q_phys(om)
E=np.eye(N,dtype=int)
th=0.31+0.07j; mu=-1j*T*th/np.pi
for order in (1,-1):
  mus=[mu+T-k*om for k in range(1,n+1)]
  if order<0: mus=mus[::-1]
  ths=[1j*np.pi*m/T for m in mus]
  Vs=[an_vector(n,q,np.exp(om*t)) for t in ths]
  wb=[sum(E[i] for i in tup) for tup in itertools.product(range(N),repeat=n)]
  for c in (1,-1,np.exp(1j*np.pi*om),np.exp(-1j*np.pi*om),-np.exp(1j*np.pi*om),-np.exp(-1j*np.pi*om),np.exp(0.5j*np.pi*om),np.exp(-0.5j*np.pi*om),1j,-1j):
    try:
        X=embed(an_vector(n,q,c*np.exp(om*th),True),Vs,[-E[i] for i in range(N)],wb); print(order,c,'OK')
    except AssertionError as e: pass
