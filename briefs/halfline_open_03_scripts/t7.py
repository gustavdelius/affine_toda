import numpy as np
from proj import cproj
om=0.37
for n in (2,3):
  T=(n+1)*om/2
  for s in (1,-1):
    th=0.31+0.07j
    c0=cproj(n,s,om,th)[0]
    for per,lab in ((1j*np.pi/T,'i pi/T'),(2j*np.pi/T,'2 i pi/T'),(1j*np.pi/om,'i pi/om'),(2j*np.pi/om,'2i pi/om')):
        print(n,s,lab,np.round(cproj(n,s,om,th+per)[0]/c0,6))
