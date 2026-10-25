"""Secondary integration checks, AFTER independent algebra_checks/loop_check.

Build a new scratch candidate from the frozen author's script and run it.
Never change the original source. Check the numerical BK insertion separately.
"""
import sys,runpy
from pathlib import Path
import numpy as np
root=Path(sys.argv[1]).resolve()
src=(root/'code/rmatrix/smatrix/oneloop_ratio.py').read_text()
src=src.replace('def analyse(name, roots, kac, h, S11, ratio_formula, expected):',
                'def analyse(name, roots, kac, h, S11, ratio_formula, expected, p1=0):')
src=src.replace('tree_over_4m2sinh(m2, C3, C4, 0, th)','tree_over_4m2sinh(m2, C3, C4, p1, th)')
src=src.replace('a5ratio, "1/2 (H = h + B/2)")','a5ratio, "1/2 (H = h + B/2)", p1=1)')
candidate=Path(__file__).resolve().parent/'candidate_oneloop_ratio.py'
candidate.write_text(src)
runpy.run_path(str(candidate))
print('PASS minimal candidate executes all three original examples')

N=15; keep=slice(4,N-4); shift=np.diag(np.ones(N-1),-1)
qvals=(-np.exp(-1j*np.pi*np.sqrt(2)),.83*np.exp(.41j))
def rel(a,b): return np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1e-300)
def dg(x,y,q):
    q3=q*q+1+q**-2
    return x@x@x@y-q3*x@x@y@x+q3*x@y@x@x-y@x@x@x
for q in qvals:
    nums=np.arange(N)-N//2; diag=np.diag(q**(-2.*nums)); inv=np.diag(q**(2.*nums))
    c=.9-.3j; a=1.3-.4j; b=.8+.7j; f=1j*c/(q-1/q)
    P=f*(a*shift+shift.T/a); Q=f*(b*diag+inv/b); rho=c*c*(q+1/q)**2
    err=max(rel((P@Q@Q@P-Q@P@P@Q)[:,keep],rho*(P@P-Q@Q)[:,keep]),
            rel(dg(P,Q,q)[:,keep],rho*(P@Q-Q@P)[:,keep]),
            rel(dg(Q,P,q)[:,keep],rho*(Q@P-P@Q)[:,keep]))
    print('BK finite-band snippet residual',err)
    assert err<1e-12
print('PASS numerical BK snippet at irrational physical parameter and generic complex q')
