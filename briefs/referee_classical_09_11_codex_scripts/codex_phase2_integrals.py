"""Canonical E,P and independently constructed local Q2,Q4 on A2 tau fields."""
from pathlib import Path
import json,itertools,math
from functools import lru_cache
import mpmath as mp
import sympy as s
mp.mp.dps=32
laws=json.loads(Path('referee_scratch/classical_09_11/codex/codex_phase2_currents_v3.json').read_text())
variables=s.symbols('u1 u2 u3 v1 v2 v3')
funs={spin:(s.lambdify(variables,s.sympify(law['T']),'mpmath'),[s.lambdify(variables,s.sympify(e),'mpmath') for e in law['X']]) for spin,law in laws.items()}

def setup(data,t):
 terms=[];h=3
 for mask in itertools.product([0,1],repeat=len(data)):
  c=mp.mpc(1);k=mp.mpc(0);nu=mp.mpc(0);species=0
  for i,(a,theta,xi) in enumerate(data):
   if mask[i]:
    ma=2*mp.sin(mp.pi*a/h);k+=ma*mp.cosh(theta);nu-=ma*mp.sinh(theta);species+=a;c*=mp.exp(xi)
    for l in range(i):
     if mask[l]:
      b,th,_=data[l];ch=mp.cosh(theta-th)
      c*=(ch-mp.cos(mp.pi*(a-b)/h))/(ch-mp.cos(mp.pi*(a+b)/h))
  terms.append((c*mp.exp(nu*t),k,nu,species))
 @lru_cache(maxsize=18000)
 def densities(x):
  tau=[];dx=[];dt=[];L=[];lxx=[];lxt=[];ltt=[]
  for j in range(3):
   ts=[(c*mp.exp(k*x+2j*mp.pi*j*a/3),k,nu) for c,k,nu,a in terms]
   z=sum(c for c,k,nu in ts);zx=sum(c*k for c,k,nu in ts);zt=sum(c*nu for c,k,nu in ts)
   tau.append(z);dx.append(zx/z);dt.append(zt/z)
   lxx.append(sum(c*k*k for c,k,nu in ts)/z-(zx/z)**2)
   lxt.append(sum(c*k*nu for c,k,nu in ts)/z-zx*zt/z**2)
   ltt.append(sum(c*nu*nu for c,k,nu in ts)/z-(zt/z)**2)
   moments=[mp.mpc(1)]+[sum(c*(nu-k)**n for c,k,nu in ts)/z for n in range(1,4)]
   cumul=[mp.mpc(0)]
   for n in range(1,4):cumul.append(moments[n]-sum(math.comb(n-1,l-1)*cumul[l]*moments[n-l] for l in range(1,n)))
   L.append(cumul)
  # u=(i beta phi)=(log tau0-log tau1,log tau1-log tau2,log tau2-log tau0).
  ux=[dx[j]-dx[(j+1)%3] for j in range(3)];ut=[dt[j]-dt[(j+1)%3] for j in range(3)]
  z=[tau[0]*tau[2]/tau[1]**2,tau[1]*tau[0]/tau[2]**2,tau[2]*tau[1]/tau[0]**2]
  V=sum(1-zz for zz in z)
  kin=-sum(uu*uu+vv*vv for uu,vv in zip(ux,ut))/2
  energy=kin+V;mom=sum(uu*vv for uu,vv in zip(ux,ut));stress=kin-V
  err=max(abs(energy-2*sum(lxx)),abs(-mom-2*sum(lxt)),abs(stress-2*sum(ltt)))
  args=[L[0][n]-L[1][n] for n in range(1,4)]+[L[1][n]-L[2][n] for n in range(1,4)]
  Q=[]
  for spin,(T,X) in funs.items():Q.append(T(*args)+sum(zj*xx(*args) for zj,xx in zip(z,X)))
  return energy,mom,*Q,err
 return densities

cases=[('single1',[(1,mp.mpf(0),.3j)],0),('single2',[(2,mp.mpf(0),.3j)],0),('complex',[(1,mp.mpc('.3','.2'),.4j)],mp.mpf('.1')),('breather',[(1,mp.mpf('.23')+1j*mp.pi/3,1j*mp.pi/2),(2,mp.mpf('.23')-1j*mp.pi/3,1j*mp.pi/2)],mp.mpf('.17')),('breather_shift',[(1,mp.mpf('.23')+1j*mp.pi/3,mp.mpf('.2')+1j*mp.pi/2),(2,mp.mpf('.23')-1j*mp.pi/3,-mp.mpf('.2')+1j*mp.pi/2)],mp.mpf('.43'))]
answers={}
for name,data,t in cases:
 fn=setup(data,t)
 values=[mp.quad(lambda x:fn(x)[j],[-40,-8,-3,-1,0,1,3,8,40]) for j in range(4)]
 answers[name]=values
 E=sum(6*2*mp.sin(mp.pi*a/3)*mp.cosh(theta) for a,theta,xi in data)
 P=sum(6*2*mp.sin(mp.pi*a/3)*mp.sinh(theta) for a,theta,xi in data)
 err=max(abs(values[0]-E),abs(values[1]-P));st=max(fn(mp.mpf(x))[4] for x in [-2,-1,0,1,2])
 print(('PASS' if err<mp.mpf('1e-17') and st<mp.mpf('1e-24') else 'FAIL')+f' {name} canonical E,P and stress: error={mp.nstr(err,5)}, stress={mp.nstr(st,5)}',flush=True)
 print(' E,P,Q2raw,Q4raw = '+', '.join(mp.nstr(v,23) for v in values),flush=True)
 if name in ['complex','breather','breather_shift']:
  for index,spin in [(2,2),(3,4)]:
   q1=answers['single1'][index];q2=answers['single2'][index]
   expected=sum((q1 if a==1 else q2)*mp.exp(spin*theta) for a,theta,xi in data)
   error=abs(values[index]-expected)
   print(('PASS' if error<mp.mpf('1e-16') else 'FAIL')+f' {name} spin{spin} direct integral/additive prediction: residual={mp.nstr(error,6)}, Q/q1={mp.nstr(values[index]/q1,23)}',flush=True)
 # Explicit parity: real coefficient laws in u=i phi are Theta-even. Multiply by 1/q1 (pure imaginary)
 # to normalize one-soliton q1=1; normalized even-spin charges are Theta-odd.
