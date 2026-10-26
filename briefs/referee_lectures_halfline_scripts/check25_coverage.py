import ast
from pathlib import Path
import numpy as np
from scipy.integrate import quad
# No author code executed: parse report labels to inventory the checks.
for p in sorted(Path('code/lectures').glob('halfline_l*_*.py')):
    tree=ast.parse(p.read_text())
    calls=[n for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='report']
    print(p.name,'report call sites',sorted(n.lineno for n in calls))
# Recompute the stationary vacuum's energy using roots and the printed field.
R=np.array([[-1.,0,1],[1,-1,0],[0,1,-1]])
s=np.sqrt(3)/2;ma=np.sqrt(3);chi=.25

def data(x):
    z=np.exp(ma*x);v=chi+2*np.pi*np.arange(3)/3
    tau=((1-z)**2+s*(1-z*z)+4*z*np.cos(v/2)**2)/(1+s)
    tx=ma*(2*np.cos(v)*z+2*(1-s)*z*z)/(1+s)
    return -np.log(tau)@R,-(tx/tau)@R

def en(x):
    phi,dphi=data(x)
    return dphi@dphi/2+sum(np.expm1(R@phi))
phi,_=data(0);E=quad(en,-20,0,epsabs=1e-11)[0]+2*sum(np.array([1,-1,-1])*np.expm1(R@phi/2))
print('25 actual energy',E,'README energy',2*(3-np.sqrt(3)));assert abs(E-2*(4-3*np.sqrt(3)))<1e-10
# Independent algebraic conversion check, which must hold for arbitrary operands.
rng=np.random.default_rng(513);A=rng.normal(size=(4,4));D=rng.normal(size=(4,4));K=rng.normal(size=(2,2));J=rng.normal(size=(2,2))
P=np.eye(4)[[0,2,1,3]];K1=np.kron(K,np.eye(2));K2=np.kron(np.eye(2),J)
fixed=A@K1@(P@D@P)@K2
right=(A@P)@np.kron(np.eye(2),K)@(D@P)@K2
wrong=(P@A)@np.kron(np.eye(2),K)@(P@D)@K2
print('25 corrected braid convention',np.linalg.norm(right-fixed),'printed convention',np.linalg.norm(wrong-fixed))
assert np.linalg.norm(right-fixed)<1e-12 and np.linalg.norm(wrong-fixed)>1
print('PASS 25: scope inventory, independently checked README energy and braid-convention blind spot')
