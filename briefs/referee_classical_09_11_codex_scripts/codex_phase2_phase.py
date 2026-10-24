import mpmath as mp
mp.mp.dps=40
# Euler-Maclaurin leading integral obtained by telescoping Hollowood (5.7).
# beta^2*(-i d log S/dtheta)/(2h), using lambda0=4pi and hbar=1.
def coeff(h,t):
 return mp.pi*mp.quad(lambda x:mp.cot(mp.pi*x-1j*t/2)-mp.cot(-mp.pi*x-1j*t/2),[0,mp.mpf(1)/h])
def A(h,a,b,t):return (mp.cosh(t)-mp.cos(mp.pi*(a-b)/h))/(mp.cosh(t)-mp.cos(mp.pi*(a+b)/h))
for h in [3,4]:
 for a,b in [(1,1),(1,2),(2,2)]:
  t=mp.sqrt(2);pr=mp.mpc(1);total=mp.mpc(0)
  for r in range(1,a+1):
   for q in range(1,b+1):
    shift=mp.pi*(a+1-2*r-b-1+2*q)/h;z=t+1j*shift
    pr*=A(h,1,1,z);total+=coeff(h,z)
  err=max(abs(pr-A(h,a,b,t)),abs(total+mp.log(A(h,a,b,t))))
  print(('PASS' if err<mp.mpf('1e-30') else 'FAIL'),f'h={h} a,b={a,b}: source coefficient={mp.nstr(total,15)}, book logA={mp.nstr(mp.log(A(h,a,b,t)),15)}, residual={mp.nstr(err,3)}')
