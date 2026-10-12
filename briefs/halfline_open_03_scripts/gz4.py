import numpy as np
from gz3 import X, T, lam
from fitR import fitfun
print('T=',T,'2T=',2*T,'om/2=',lam/2)
print(fitfun(X,T,dmax=10))
