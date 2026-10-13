"""Independent ch.39 referee checks; no imports from the book's code.

m=beta=1 unless explicitly stated. Roots are rows in ambient R^h.
Run from any directory with PYTHONDONTWRITEBYTECODE=1 (no file access).
"""
import itertools
import math
import numpy as np
from numpy import exp, sin, cos, sinh, cosh, pi


def geometry(h):
    roots = np.eye(h) - np.roll(np.eye(h), 1, axis=1)
    weights = np.zeros((h, h))
    weights[1:] = np.linalg.lstsq(roots[1:], np.eye(h-1), rcond=None)[0].T
    D = sum(2*np.outer(weights[j]-weights[(j+1)%h], weights[j])
            for j in range(1, h))
    P = np.eye(h)-np.ones((h,h))/h
    return roots, weights, D, P-D


def delay(h,a,theta,td,p=0,kind=-1):
    u = (theta-td if kind == -1 else td-theta)-2j*pi*p/h
    x = exp(-u)
    return (x-1j*exp(1j*pi*a/h))/(x-1j*exp(-1j*pi*a/h))


def tau_data(h, species, rapidities, xis, x, t, multipliers=None):
    n = len(species)
    if multipliers is None:
        multipliers = np.ones(n)
    q = exp(2j*pi*np.arange(h)/h)
    ms = 2*sin(pi*np.array(species)/h)
    kx, kt = ms*cosh(rapidities), -ms*sinh(rapidities)
    es = multipliers*exp(kx*x+kt*t+xis)
    tau, tx, tt = (np.zeros(h,complex) for _ in range(3))
    for flags in itertools.product([0,1], repeat=n):
        subset = [i for i in range(n) if flags[i]]
        term = q**sum(species[i] for i in subset)*np.prod(es[subset])
        for i,j in itertools.combinations(subset,2):
            term *= interaction(h,species[i],species[j],rapidities[i]-rapidities[j])
        tau += term
        tx += term*sum(kx[i] for i in subset)
        tt += term*sum(kt[i] for i in subset)
    return tau,tx,tt


def fields(h,data,weight=None):
    roots,_,_,_ = geometry(h)
    tau,tx,tt = data
    f = 1j*roots.T@np.log(tau)
    if weight is not None:
        f += 2*pi*weight
    return f,1j*roots.T@(tx/tau),1j*roots.T@(tt/tau)


def potential(h,phi,psi,td,kind=-1):
    roots,_,D,_ = geometry(h)
    cf = exp(-td if kind==-1 else td)
    cg = 1/cf
    fs = exp(.5j*roots@(D.T@phi+D@psi))
    gs = exp(.5j*roots@D@(phi-psi))
    gp = -.5j*(cf*D@roots.T@fs+cg*D.T@roots.T@gs)
    gq = -.5j*(cf*D.T@roots.T@fs-cg*D.T@roots.T@gs)
    B = -sum(cf*(fs-1)+cg*(gs-1))
    Pd = (-kind)*sum(cf*(fs-1)-cg*(gs-1))
    return gp,gq,B,Pd


def residual(h,left,right,td,kind=-1):
    _,_,D,E = geometry(h)
    phi,px,pt = left
    psi,qx,qt = right
    gp,gq,_,_ = potential(h,phi,psi,td,kind)
    s = -kind
    return max(np.max(abs(px-s*(E@pt+D@qt)+gp)),
               np.max(abs(qx-s*(D.T@pt-E@qt)-gq)))


def interaction(h,a,b,u):
    return (cosh(u)-cos(pi*(a-b)/h))/(cosh(u)-cos(pi*(a+b)/h))


def charge(h,a,phase):
    roots,_,_,_ = geometry(h)
    reduced = (2*pi*np.arange(h)*a/h+phase+pi)%(2*pi)-pi
    return -roots.T@reduced/(2*pi)


def vacuum_EP(h,p,td,kind=-1):
    q = -kind*2j*pi*p/h
    return 2*h*np.array([cosh(td)-cosh(td+q),sinh(td)-sinh(td+q)])


def linear_solve(h,k,theta,td,kind=-1,p=0,source=False,right=False):
    roots,weights,D,E = geometry(h)
    if source:
        D=D.T
        E=(np.eye(h)-np.ones((h,h))/h)-D
    w=exp(2j*pi*k/h)
    rho=roots.T@(w**np.arange(h))
    jump=2*pi*weights[p]
    cf=exp(-td if kind==-1 else td)
    cg=1/cf
    fs=exp(.5j*roots@D@jump)
    gs=1/fs
    F=np.diag(cf*fs)
    G=np.diag(cg*gs)
    A=.25*(D@roots.T@F@roots@D.T+D.T@roots.T@G@roots@D)
    B=.25*(D@roots.T@F@roots@D-D.T@roots.T@G@roots@D)
    C=B.T
    H=.25*(D.T@roots.T@(F+G)@roots@D)
    nu=2*sin(pi*k/h)*cosh(theta)
    mom=2*sin(pi*k/h)*sinh(theta)
    s=-kind
    def eq(R,T):
        if right:
            phi,psi=T*rho,(1+R)*rho
            px,qx=-1j*mom*T*rho,-1j*mom*(1-R)*rho
        else:
            phi,psi=(1+R)*rho,T*rho
            px,qx=1j*mom*(1-R)*rho,1j*mom*T*rho
        pt,qt=-1j*nu*phi,-1j*nu*psi
        return np.r_[px-s*(E@pt+D@qt)+A@phi+B@psi,
                     qx-s*(D.T@pt-E@qt)-C@phi-H@psi]
    zero=eq(0,0)
    mat=np.stack([eq(1,0)-zero,eq(0,1)-zero],axis=1)
    amps=np.linalg.lstsq(mat,-zero,rcond=None)[0]
    return amps,max(abs(eq(*amps)))


def main():
    worst=0.; count=0
    for h in [3,4]:
      for a,p,kind in itertools.product(range(1,h),range(h),[-1,1]):
       for theta in [0.73,0.83+.17j]:
        z=delay(h,a,theta,.31,p,kind)
        for x,t in [(-.21,-.53),(0,.12),(.19,.71)]:
          left=fields(h,tau_data(h,[a],np.array([theta]),np.array([.23+.37j]),x,t))
          right=fields(h,tau_data(h,[a],np.array([theta]),np.array([.23+.37j]),x,t,np.array([z])),geometry(h)[1][p])
          worst=max(worst,residual(h,left,right,.31,kind)); count+=1
    print(f'1 direct exponentials, exact one-soliton profiles: {count} samples, residual {worst:.3e}')
    worst=0.; count=0
    for h in [3,4,6]:
     for a,p,kind,td in itertools.product(range(1,h),range(h),[-1,1],[-.6,.6]):
      theta=td+kind*1j*(pi/2-pi*(a+2*p)/h)
      q=exp(2j*pi*np.arange(h)*a/h)
      vacuum=(q,np.zeros(h,complex),np.zeros(h,complex))
      empty=(np.ones(h,complex),np.zeros(h,complex),np.zeros(h,complex))
      for t in [-.27,.41]:
        sol=tau_data(h,[a],np.array([theta]),np.array([.18+.271j]),0,t)
        first=residual(h,fields(h,sol),fields(h,vacuum,geometry(h)[1][p]),td,kind)
        second=residual(h,fields(h,empty),fields(h,sol,geometry(h)[1][p]),td,kind)
        worst=max(worst,first,second);count+=2
    print(f'2 pole limits: {count} configurations, residual {worst:.3e}')
    for h in [4,6,8]:
      real=[(a,p) for a in range(1,h) for p in range(h) if (2*a+4*p-h)%(4*h)==0]
      print(f'2 real poles h={h}: {real}')
    print('2 zero-envelope counterexample: h=4 a=2 p=1 theta=.6+i*pi/2; cos(Im theta)=0')
    for kind in [-1,1]:
      actual=vacuum_EP(6,2,.6,kind)-vacuum_EP(6,1,.6,kind)
      expected=12*np.array([cosh(.6),sinh(.6)])
      print(f'2 h=6 a=1 p=1 type={kind}: stored EP={actual}, error={max(abs(actual-expected)):.3e}')
    print(f'3 static counterexample h=3 p=1 td=0: EP={vacuum_EP(3,1,0)}')
    shell=0.; ident=0.
    for h in range(2,10):
     for p,kind in itertools.product(range(h),[-1,1]):
      actual=vacuum_EP(h,p,.47,kind)-2*h*np.array([cosh(.47),sinh(.47)])
      shell=max(shell,abs(actual[0]**2-actual[1]**2-(2*h)**2))
      for a in range(1,h):
        th=.47+kind*1j*(pi/2-pi*(a+2*p)/h)
        change=vacuum_EP(h,p+a,.47,kind)-vacuum_EP(h,p,.47,kind)
        ident=max(ident,max(abs(change-4*h*sin(pi*a/h)*np.array([cosh(th),sinh(th)]))))
    print(f'3 mass shell residual {shell:.3e}; 2 absorption EP identity {ident:.3e}')
    for h,a in [(3,1),(3,2),(4,2),(5,2),(6,2),(8,3)]:
      spacing=2*pi*math.gcd(a,h)/h
      boundaries=np.unique(np.round((pi-2*pi*np.arange(h)*a/h)%(2*pi),12))
      phases=(boundaries+.01)%(2*pi)
      charges=np.array([charge(h,a,f) for f in phases])
      hits=set()
      for u in [-7,0,7]:
       arg=np.angle(delay(h,a,u,0))
       for f in np.linspace(0,2*pi,701,endpoint=False)+.00037:
        steps=math.floor((f+arg-boundaries[0])/spacing)-math.floor((f-boundaries[0])/spacing)
        hits.add(abs(steps))
      print(f'3/6 sectors h={h} a={a}: count={len(boundaries)}, spacing/pi={spacing/pi:g}, attainable |steps|={sorted(hits)}')
      if (h,a)==(5,2):
        print('9 exercise charges h=5,a=2 (ambient coordinates):',np.round(charges,3).tolist())
        print('9 exercise example phase=.55, u=3:',np.round(charge(h,a,.55),3), '->', np.round(charge(h,a,.55+np.angle(delay(h,a,3,0))),3))
    mixed=linear=energy=0.
    for h in range(2,11):
     for c,b in itertools.product(range(1,h),repeat=2):
      tc=.39;v=tc+1j*(pi/2+pi*c/h)
      for th in [-.8,.21,1.1]:
       actual=delay(h,b,th,v)*delay(h,b,th,np.conj(v),kind=1)
       mixed=max(mixed,abs(actual-exp(2j*pi*b/h)*interaction(h,c,b,th-tc)))
       actual=delay(h,b,th+1j*pi/2,v)*delay(h,b,th+1j*pi/2,np.conj(v),kind=1)
       z=1j*sinh(th-tc)
       X=(z-cos(pi*(c-b)/h))/(z-cos(pi*(c+b)/h))
       linear=max(linear,abs(actual-exp(2j*pi*b/h)*X))
      ep=-2*h*np.array([cosh(v)+cosh(np.conj(v)),sinh(v)+sinh(np.conj(v))])
      energy=max(energy,max(abs(ep-4*h*sin(pi*c/h)*np.array([cosh(tc),sinh(tc)]))))
    print(f'4 mixed pair h=2..10: soliton {mixed:.3e}, particle {linear:.3e}, EP {energy:.3e}')
    print('4 real-sigma realizations h<=8:',[(h,c,p) for h in range(2,9) for c in range(1,h) for p in range(h) if (4*p-h-2*c)%(4*h)==0])
    err=Rmax=cross=uni=src=0.
    for h in range(2,8):
     for k,kind,p in itertools.product(range(1,h),[-1,1],range(h)):
      th=.83;td=.27
      (R,T),res=linear_solve(h,k,th,td,kind,p)
      Rmax=max(Rmax,abs(R),res)
      err=max(err,abs(T-delay(h,k,th+1j*pi/2,td,p,kind)))
      (Rbar,Tbar),_=linear_solve(h,k,th,td,kind,p,right=True)
      cross=max(cross,abs(Tbar-delay(h,h-k,1j*pi-th+1j*pi/2,td,p,kind)))
      uni=max(uni,abs(delay(h,k,th+1j*pi/2,td,p,kind)*delay(h,h-k,th+1.5j*pi,td,p,kind)-1))
      if p==0 and kind==-1:
       (RR,TT),_=linear_solve(h,k,th,td,source=True)
       mk=2*sin(pi*k/h);w=exp(2j*pi*k/h);sig=exp(-td)
       corr=(1j*exp(-th)*mk-sig*(1-w**-1))/(1j*exp(-th)*mk+sig*(1-w))
       src=max(src,abs(corr-TT))
    print(f'5 Hessian scattering: R/residual {Rmax:.3e}, T {err:.3e}, crossing {cross:.3e}, algebraic unitarity {uni:.3e}, source system {src:.3e}')
    for k in [1,2]:
      th=.83;td=.27;u=th-td;mk=2*sin(pi*k/3);w=exp(2j*pi*k/3);sig=exp(-td)
      corr=(1j*exp(-th)*mk-sig*(1-w**-1))/(1j*exp(-th)*mk+sig*(1-w))
      printed=(1j*exp(-th)*mk-sig*(1-w**-1))/(1j*exp(th)*mk+sig*(1-w))
      b=(exp(-1j*pi/3)*sinh(u/2-1j*pi/6)/sinh(u/2+1j*pi/6) if k==1 else exp(1j*pi/3)*sinh(u/2+2j*pi/3)/sinh(u/2+1j*pi/3))
      print(f'5 source a2 k={k}: corrected-(7.5)={abs(corr-b):.3e}; printed |T|={abs(printed):.9f}')
    rotation=0.
    for h in [3,4,6]:
      roots,_,D,_=geometry(h);P=np.roll(np.eye(h),-1,axis=0)
      # P e_j=e_(j-1), hence rows shift in the opposite direction.
      P=np.column_stack([np.eye(h)[:,(j-1)%h] for j in range(h)])
      rotation=max(rotation,max(abs((D@P+D.T).ravel())))
    print(f'7 leading constraint D*pi=-D^T: {rotation:.3e}; z_limit=omega^a algebraically')
    conj=0.; advances=0
    for h in range(2,10):
     for a,p,kind in itertools.product(range(1,h),range(h),[-1,1]):
      th=.83+.12j;td=.27-.08j
      conj=max(conj,abs(delay(h,h-a,np.conj(th),np.conj(td),-p,kind)-np.conj(delay(h,a,th,td,p,kind))))
      z=delay(h,a,.83,.27,p,kind)
      if cos(2*pi*p/h)<-1e-9:
        assert abs(z)<1; advances+=1
    print(f'8 Theta residual {conj:.3e}; advanced samples {advances}')
    print('8 real negative sigma counterexample: theta_d=i*pi, p=0, h=3, energies',vacuum_EP(3,0,1j*pi))
    multi=0.
    for h in [3,4]:
     for kind,p in itertools.product([-1,1],[0,1]):
      aa=[1,2];th=np.array([.9,-.4]);xi=np.array([.1+.31j,-.23+.19j])
      zz=np.array([delay(h,a,t,.2,p,kind) for a,t in zip(aa,th)])
      L=fields(h,tau_data(h,aa,th,xi,0,.31))
      Q=fields(h,tau_data(h,aa,th,xi,0,.31,zz),geometry(h)[1][p])
      multi=max(multi,residual(h,L,Q,.2,kind))
    bootstrap=0.;reverse_bootstrap=0.
    for h in range(3,9):
     for a in range(1,h):
      for b in range(1,h-a):
       u=.61+.17j
       bootstrap=max(bootstrap,abs(delay(h,a,u+1j*pi*b/h,0)*delay(h,b,u-1j*pi*a/h,0)-delay(h,a+b,u,0)))
       reverse_bootstrap=max(reverse_bootstrap,abs(delay(h,a,u-1j*pi*b/h,0)*delay(h,b,u+1j*pi*a/h,0)-delay(h,a+b,u,0)))
    print(f'9 two-soliton residual {multi:.3e}; bootstrap printed {bootstrap:.3e}; reversed shifts {reverse_bootstrap:.3e}')


if __name__=='__main__':
    main()
