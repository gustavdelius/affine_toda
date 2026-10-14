"""Calibrate the breather-singlet construction of b1_dn_matrices.py on a_n, where the book's
c_+ (prp-an-nondiagonal-scalar) was established in an earlier round: singlet in V(th_a) (x) Vbar(th_b),
map (1 x K(th_a)) R (1 x K(th_b)^T) with K from eq-an-k-matrix."""
import sys, numpy as np
sys.path.insert(0,'.'); sys.path.insert(0,'../../textbook_code')
import boundary_common as bc
from b1_dn_matrices import singlet
pi=np.pi
def Kan(n, om, th, s):
    N=n+1; q=bc.q_phys(om); X=np.exp(om*(th-1j*pi/2)); es=1/(2*np.cos(pi*om/2))
    K=np.zeros((N,N),complex)
    for i in range(N):
        for j in range(N):
            if j>i: K[i,j]=X**(j-i-N/2)
            elif j<i: K[i,j]=s*X**(j-i+N/2)
            else: K[i,i]=((-q)**0.5*X**(-N/2)-s*(-q)**(-0.5)*X**(N/2))/(es*(1/q-q))
    return K
def c_an(mu,T,om,s):
    s2=np.sin(pi*om/2)**2
    if s>0: return -np.cos(pi*mu)*np.cos(pi*(mu+om/2))*np.cos(pi*(mu-T))/(s2*np.cos(pi*(mu-om/2)))
    return np.sin(pi*mu)*np.sin(pi*(mu+om/2))*np.sin(pi*(mu-T))/(s2*np.sin(pi*(mu-om/2)))
for n in (2,3):
  om=0.37; q=bc.q_phys(om); T=(n+1)*om/2; N=n+1
  for s in (1,-1):
    for th in (0.31+0.17j,-0.2+0.4j):
      mu=-1j*T*th/pi
      tha=th-1j*pi/2+1j*pi/(2*T); thb=th+1j*pi/2-1j*pi/(2*T)
      ya,yb=np.exp(om*tha),np.exp(om*thb)
      # check K is an intertwiner V(y) -> Vbar(1/y)
      Kt=Kan(n,om,tha,s)
      es=1/(2*np.cos(pi*om/2)); eps=[s*es]+[es]*n; c=[1/q]*N
      V=bc.an_vector(n,q,ya); Vb=bc.an_vector(n,q,1/ya,conj=True)
      r=max(np.abs(Kt@bc.coideal_gen(V,j,c[j],eps[j])-bc.coideal_gen(Vb,j,c[j],eps[j])@Kt).max() for j in range(N))
      A=bc.an_vector(n,q,ya); Bb=bc.an_vector(n,q,yb,conj=True)
      sv=singlet(A,Bb); sv=sv/sv[0]
      R,_=bc.rcheck(bc.an_vector(n,q,ya), bc.an_vector(n,q,1/yb))     # V(th_a) x V(-th_b) -> V(-th_b) x V(th_a)
      w=bc.one_k(Kan(n,om,tha,s),N)@R@bc.one_k(Kan(n,om,thb,s).T,N)@sv
      sv2=singlet(bc.an_vector(n,q,1/yb),bc.an_vector(n,q,1/ya,conj=True))
      i1=0*N+0   # hw of V (v_1) times lowest of Vbar (vbar_1, weight -eps_1)
      sv=sv/sv[i1]; sv2=sv2/sv2[i1]
      cc=w[i1]; prop=np.abs(w-cc*sv2).max()/np.abs(w).max()
      print(f"a_{n} s={s:+d} th={th}: K intertw {r:.1e}, prop {prop:.1e}, cc/c_+- = {cc/c_an(mu,T,om,s):.6f}")
print("K = 1 (eps-hat = 0) in the same construction:")
for n in (2,3):
  om=0.37; q=bc.q_phys(om); T=(n+1)*om/2; N=n+1
  for th in (0.31+0.17j,-0.2+0.4j):
      tha=th-1j*pi/2+1j*pi/(2*T); thb=th+1j*pi/2-1j*pi/(2*T)
      ya,yb=np.exp(om*tha),np.exp(om*thb)
      sv=singlet(bc.an_vector(n,q,ya),bc.an_vector(n,q,yb,conj=True)); sv=sv/sv[0]
      R,_=bc.rcheck(bc.an_vector(n,q,ya), bc.an_vector(n,q,1/yb))
      I=np.eye(N)
      w=bc.one_k(I,N)@R@bc.one_k(I,N)@sv
      sv2=singlet(bc.an_vector(n,q,1/yb),bc.an_vector(n,q,1/ya,conj=True))
      sv=sv/sv[0]; sv2=sv2/sv2[0]; cc=w[0]
      print(f"  a_{n} th={th}: prop {np.abs(w-cc*sv2).max()/np.abs(w).max():.1e}, coefficient {cc:.6f}")
