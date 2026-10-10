"""eq-ybe-braid in terms of R = P Rc: R12(x1/x2) R13(x1/x3) R23(x2/x3) = R23 R13 R12 (eq-yang-baxter)."""
import numpy as np
q = (0.83+0.41j)**2
def Rc(u):
    return np.array([[1,0,0,0],[0,u*(q**2-1)/(q**2-u),q*(u-1)/(u-q**2),0],[0,q*(u-1)/(u-q**2),(q**2-1)/(q**2-u),0],[0,0,0,1]])
P = np.eye(4)[[0,2,1,3]]
R = lambda u: P@Rc(u)
I = np.eye(2)
P23 = np.kron(I, P)
R12 = lambda u: np.kron(R(u), I); R23 = lambda u: np.kron(I, R(u)); R13 = lambda u: P23@R12(u)@P23
x1, x2, x3 = 0.7+0.2j, -0.4+1.1j, 1.3-0.5j
lhs = R12(x1/x2)@R13(x1/x3)@R23(x2/x3); rhs = R23(x2/x3)@R13(x1/x3)@R12(x1/x2)
print(('PASS ' if np.allclose(lhs, rhs) else 'FAIL ') + 'R12(x1/x2)R13(x1/x3)R23(x2/x3) = R23 R13 R12 with R = P Rc')
