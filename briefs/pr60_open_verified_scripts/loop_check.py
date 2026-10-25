"""Fresh Lagrangian tensor expansion and analytic tree block slope, no book imports."""
import itertools
import numpy as np
from scipy.integrate import quad

def tensors(A,n):
    A=np.array(A); n=np.array(n)
    assert np.linalg.norm(n@A)<1e-14
    H=A.T@np.diag(n)@A
    masses,rotation=np.linalg.eigh(H)
    pol=A@rotation
    tensors=[]
    for degree in (3,4):
        arr=np.zeros((len(masses),)*degree)
        for inds in itertools.product(range(len(masses)),repeat=degree):
            arr[inds]=sum(n[j]*np.prod(pol[j,list(inds)]) for j in range(len(n)))
        tensors.append(arr)
    return masses,*tensors

def shifts(ms,C):
    result=[]
    for i,a in enumerate(ms):
        val=0
        for j,b in enumerate(ms):
            for k,c in enumerate(ms):
                if abs(C[i,j,k])<1e-12: continue
                val-=C[i,j,k]**2*quad(lambda z:1/(z*b+(1-z)*c-z*(1-z)*a),0,1,
                                             epsabs=1e-13,epsrel=1e-13)[0]/(8*np.pi)
        result.append(val)
    return np.array(result)

# d/dB Im log {x}_B at B=0. At B=0 the brace is identically 1,
# so any derivative H'(0), or x'(0), cancels exactly.
def slope_brace(x,H,t):
    def slope_block(z):
        p=np.pi*z/(2*H)
        return (1j*np.pi/(2*H)*(1/np.tanh(t/2+1j*p)+1/np.tanh(t/2-1j*p))).imag
    return slope_block(x+1)-slope_block(x-1)

def tree(ms,C,D,i,t):
    a=ms[i]
    invariants=[2*a*(1+np.cosh(t)),0,2*a*(1-np.cosh(t))]
    amp=-D[i,i,i,i]
    for k,b in enumerate(ms): amp-=C[i,i,k]**2*sum(1/(u-b) for u in invariants)
    return amp/(4*a*np.sinh(t))

e=np.eye(3); z=1/np.sqrt(2)
root_c=[z*(e[0]-e[1]),z*(e[1]-e[2]),2*z*e[2],-2*z*e[0]]
root_a=[z*(e[0]-e[1]),z*(e[1]-e[2]),2*z*e[2],-z*(e[0]+e[1])]
e=np.eye(4)
root_f=[e[1]-e[2],e[2]-e[3],e[3],(e[0]-e[1]-e[2]-e[3])/2,-e[0]-e[1]]
def rc(H):
    x=np.sin(np.pi*np.arange(1,4)/H)**2
    return x/x[0]
def ra(H):
    x=np.sort(np.r_[2*np.sin(np.pi*np.arange(1,3)/H),1])**2
    return x/x[0]
def rf(H):
    Hp=6*H/(H-6)
    x=np.sort([np.sin(np.pi/H)*np.sin(2*np.pi/Hp),np.sin(3*np.pi/H)*np.sin(np.pi/Hp),
               np.sin(2*np.pi/H)*np.sin(2*np.pi/Hp),np.sin(3*np.pi/H)*np.sin(2*np.pi/Hp)])**2
    return x/x[0]
for name,A,n,H,ratio,particle,target in [
    ('c3',root_c,[2,2,1,1],6,rc,0,1),('a5tw',root_a,[1,2,1,1],5,ra,1,.5),
    ('f4',root_f,[2,3,4,2,1],12,rf,0,3)]:
    ms,C,D=tensors(A,n); dm=shifts(ms,C)
    dr=dm/ms[0]-ms*dm[0]/ms[0]**2
    step=1e-4
    deriv=(ratio(H-2*step)-8*ratio(H-step)+8*ratio(H+step)-ratio(H+2*step))/(12*step)
    delta=dr[abs(deriv)>1e-7]/deriv[abs(deriv)>1e-7]
    xs=[1,H-1] if name!='f4' else [1,H-1,H/3+1,H-H/3-1]
    kap=[tree(ms,C,D,particle,t)/sum(slope_brace(x,H,t) for x in xs) for t in [.4,.9,1.5,2.3]]
    val=np.mean(delta)/np.mean(kap)
    print(name,'m2=',ms,'dm2=',dm,'deltaH=',delta,'tree/B=',kap,'ratio=',val,flush=True)
    assert np.ptp(delta)<1e-9 and np.ptp(kap)<1e-12 and abs(val-target)<1e-8
    print('PASS',name,'particle-index-corrected ratio',flush=True)
    if name=='a5tw':
        wrong=[tree(ms,C,D,0,t)/sum(slope_brace(x,H,t) for x in xs) for t in [.4,.9,1.5,2.3]]
        print('old particle=0:',wrong,'ratio',np.mean(delta)/np.mean(wrong),flush=True)
        assert np.ptp(wrong)>.01
        print('PASS old a5tw index has wrong rapidity shape',flush=True)

# Reject the old referee replacement before running: its c3 marks do not sum to zero.
bad_marks=np.array([1,2,2,1]); bad_roots=np.array(root_c)/z
print('old report replacement c3 sum(n alpha)=',bad_marks@bad_roots)
assert np.linalg.norm(bad_marks@bad_roots)>1
print('PASS old report wholesale replacement fails affine-root constraint')
