"""Independent checks for corrected Gemini lectures 1, 2 and 6.

Run from any directory. Derivations use the displayed actions and amplitudes,
without importing the book's or referee's checking code.
"""
import sympy as s
import mpmath as mp


def check(name, *residuals):
    ok = all(s.simplify(s.expand(r)) == 0 for r in residuals)
    print(('PASS' if ok else 'FAIL') + ' ' + name, flush=True)
    if not ok:
        raise AssertionError((name, residuals))


# Ising crossing: sinh(i*pi-theta)=sinh(theta), cosh changes sign.
z, c, sh = s.symbols('z c sh', real=True)
T = c*z/(z-s.I*sh)
R = s.I*sh*s.sqrt(1+z*z)/(z-s.I*sh)
check('Ising two-channel unitarity',
      (T*T.subs(z,-z)+R*R.subs(z,-z)-1).subs(c*c,1-sh*sh),
      T*R.subs(z,-z)+R*T.subs(z,-z))
theta, chi = s.symbols('theta chi', real=True)
Ti = s.cos(chi)*s.sinh(theta)/(s.sinh(theta)-s.I*s.sin(chi))
Ri = s.I*s.sin(chi)*s.cosh(theta)/(s.sinh(theta)-s.I*s.sin(chi))
check('Ising same-channel crossing', Ti.subs(theta,s.I*s.pi-theta)-Ti,
      Ri.subs(theta,s.I*s.pi-theta)+Ri)
# At theta=0 and chi=pi/6: R=-1 whereas T(i*pi/2)=sqrt(3).
check('Mixed R=T crossed-channel claim has explicit counterexample',
      Ri.subs({theta:0,chi:s.pi/6})+1,
      Ti.subs({theta:s.I*s.pi/2,chi:s.pi/6})-s.sqrt(3))

fp,gp,phit,psit,V = s.symbols('fp gp phit psit V')
Bphi,Bpsi=(fp+gp)/2,(fp-gp)/2
phix,psix=psit-Bphi,phit+Bpsi
flux=(psit**2+psix**2-phit**2-phix**2)/2+V
check('Momentum functional equation and local primitive',
      (Bphi**2-Bpsi**2)/2-fp*gp/2,
      flux.subs(V,fp*gp/2)-(phit*Bpsi+psit*Bphi))
check('Light-cone sum and difference from sewing',
      phit-psit+phix-psix+fp,phit+psit-phix-psix-gp)

ph,ps,bet,m0=s.symbols('ph ps beta m0', real=True, nonzero=True)
sig=s.symbols('sig', positive=True)
potential=-2*m0/bet**2*(sig*s.cos(bet*(ph+ps)/2)+s.cos(bet*(ph-ps)/2)/sig)
check('Sine-Gordon defect potential gradient normalization',
      s.diff(potential,ph)-m0/bet*(sig*s.sin(bet*(ph+ps)/2)+s.sin(bet*(ph-ps)/2)/sig))

# Set m0=1, p=exp(-theta)>0; solve the two independent plane-wave equations.
p,sigma=s.symbols('p sigma', positive=True)
rr,tt=s.symbols('rr tt')
e1=(sigma-s.I*p)+(sigma-s.I/p)*rr+(sigma+s.I*p)*tt
e2=(1/sigma+s.I/p)+(1/sigma+s.I*p)*rr-(1/sigma-s.I/p)*tt
sol=s.solve([e1,e2],[rr,tt])
trans=(p+s.I*sigma)/(p-s.I*sigma)
check('Linear reflection and transmission derived from sewing',sol[rr],sol[tt]-trans)
x=s.symbols('x', nonzero=True)
t=(x+s.I)/(x-s.I)
sinh_u=(1/x-x)/2
block=(sinh_u-s.I)/(sinh_u+s.I)
check('Classical amplitude square, inverse and crossing', t*t-block,
      1/t**2-1/block,t*t.subs(x,-x)-1)
check('Classical amplitude in sinh/cosh form',t-(s.I-sinh_u)/((1/x+x)/2))
# Negative controls against the report's proposed replacements.
check('Report plane wave has the wrong massive dispersion sign',
      s.cosh(theta)**2-s.sinh(theta)**2-1)
check('Report delay fraction is reciprocal of the actual delay',
      (p-sigma)/(p+sigma)*(p+sigma)/(p-sigma)-1)


# Real kink phi=4 atan(Y); E=iY and theta>0.
Y,delay=s.symbols('Y delay', real=True)
phi=s.Function('phi')
zdelay=(p+sigma)/(p-sigma)
den=(1+Y*Y)*(1+delay*delay*Y*Y)
sin_sum=2*Y*(1+delay)*(1-delay*Y*Y)/den
sin_diff=2*Y*(1-delay)*(1+delay*Y*Y)/den
plus=4*p*Y*(1/(1+Y*Y)-delay/(1+delay*delay*Y*Y))
minus=-4*Y/p*(1/(1+Y*Y)+delay/(1+delay*delay*Y*Y))
check('Soliton delay solves both nonlinear sewing equations',
      (plus+2*sigma*sin_sum).subs(delay,zdelay),
      (minus-2/sigma*sin_diff).subs(delay,zdelay))
check('One-soliton light-cone conversion',
      (p-(-1/p))/2-(p+1/p)/2,
      (p+(-1/p))/2+(1/p-p)/2)
check('Kink logarithm has positive Y derivative',
      s.diff(2*s.I*s.log((1-s.I*Y)/(1+s.I*Y)),Y)-4/(1+Y*Y))
a,b,X,Z=s.symbols('a b X Z')
# tan X=(a+b)/(a-b) tan Z, after removing cos X cos Z.
check('Bianchi algebraic compatibility',
      (a*(Z-X)+b*(Z+X)).subs(X,(a+b)/(a-b)*Z))

# General fused kinetic coefficient: EL(phi)=c psi_t-B_phi;
# EL(psi)=-c phi_t-B_psi-VB_prime. Hence coefficient must be one.
coef,time=s.symbols('coef time')
fld=s.Function('phi')(time); aux=s.Function('psi')(time)
kin=coef*fld*s.diff(aux,time)
el_phi=s.diff(kin,fld)-s.diff(s.diff(kin,s.diff(fld,time)),time)
el_psi=s.diff(kin,aux)-s.diff(s.diff(kin,s.diff(aux,time)),time)
check('Fused symplectic variation has unit coefficient',
      el_phi-coef*s.diff(aux,time),el_psi+coef*s.diff(fld,time))

# Derive fusion directly by multiplication, including independent off-diagonals.
yp,ym,kp,km,g=s.symbols('yp ym kp km g', nonzero=True)
A=s.Matrix([[1,-s.I*yp],[-s.I*yp,1]])
B=s.Matrix([[1,s.I*ym],[s.I*ym,1]])
G=s.diag(g,1/g)
F=B*s.diag(kp,km)*A
expected=s.Matrix([[kp+yp*ym*km,-s.I*(yp*kp-ym*km)],
                   [-s.I*(yp*km-ym*kp),km+yp*ym*kp]])
check('Four entries of fused boundary matrix',*(F-expected))
check('Fused charge-sector factorization',*(G.inv()*F*G-(G.inv()*B)*s.diag(kp,km)*(A*G)))
check('Reverse transmission scalar normalization',*(B*A.subs(yp,ym)-(1+ym*ym)*s.eye(2)))
y,d,e=s.symbols('y d e',nonzero=True)
# e=exp(i eta), d=exp(-lambda theta_d), y=exp(lambda theta).
ks={kp:(y*e+1/(y*e))/2,km:(y/e+e/y)/2,yp:d*y,ym:d/y}
F=F.subs(ks,simultaneous=True)
check('Off-diagonals require opposite constant eta phases',
      F[0,1]+s.I*d*e*(y*y-1/y**2)/2,
      F[1,0]+s.I*d/e*(y*y-1/y**2)/2)
H=s.diag(e,-1) # same similarity as diag(exp(i eta/2),-exp(-i eta/2))
Fsym=H.inv()*F*H
check('Constant diagonal gauge gives book GZ off-diagonals',
      Fsym[0,1]-Fsym[1,0],Fsym[0,1]-s.I*d*(y*y-1/y**2)/2)
lam,td,tdcl=s.symbols('lambda td tdcl',nonzero=True)
check('Conditional UV-IR relation solves vartheta=lambda td=(lambda+1)tdcl',
      (lam*td-(lam+1)*tdcl).subs(td,(lam+1)*tdcl/lam))

# Ordered singlet crossing: state=(s1,s2,n), s_i=+/-1; W charge symbolic q^n=g.
q,y1,y2=s.symbols('q y1 y2',nonzero=True)
def transmit(state,slot,yval):
    spins=list(state[:2]); shift=state[2]
    power=g*q**shift
    diag=1/power if spins[slot]==1 else power
    shifted=shift+spins[slot]
    flipped=spins.copy(); flipped[slot]*=-1
    return [(tuple(spins)+(shift,),diag),(tuple(flipped)+(shifted,),-s.I*yval*diag)]
out={}
for state,amp in [((1,-1,0),1),((-1,1,0),-1)]:
    for mid,c1 in transmit(state,1,y2):
        for end,c2 in transmit(mid,0,y1):
            out[end]=out.get(end,0)+amp*c1*c2
for state in out:
    target=(1+q*y1*y2)*({(1,-1,0):1,(-1,1,0):-1}.get(state,0))
    out[state]=(out[state]-target).subs(y1,q*y2)
check('Ordered singlet fusion cancels all charge shifts',*out.values())

# Particle bootstrap as rational identities; u,v are exp(i gamma_a,b).
u,v=s.symbols('u v',nonzero=True)
def particle(xx,uu): return (xx+uu)/(xx+1/uu)
check('Toda particle crossing and bootstrap a+b<h',
      particle(-x,-1/u)*particle(x,u)-1,
      particle(x/v,u)*particle(x*u,v)-particle(x,u*v))
check('Toda particle bootstrap a+b>h',
      particle(-x*v,u)*particle(-x/u,v)-particle(x,-u*v))
# Derive cancellation in the minimal scalar's logarithmic Fourier kernel.
k,uu,ll=s.symbols('k uu ll', positive=True)
base=s.pi*k/2+s.I*k*(uu+s.I*s.pi/2)
shift=s.pi*k/(2*ll)
res=(s.expand_trig(s.cosh(base+shift)-s.cosh(base-shift)) /
     (2*k*s.sinh(shift)*s.cosh(s.pi*k/2))-s.I*s.sin(k*uu)/(k*s.cosh(s.pi*k/2)))
check('Scalar Fourier kernel loses all lambda dependence',res)
# Evaluate the coupling-cancelled Fourier kernel independently by quadrature.
mp.mp.dps=30
for uv in (-1.2,0.0,0.7):
    integral=mp.quad(lambda k: mp.sin(k*uv)/(k*mp.cosh(mp.pi*k/2)) if k else uv,[0,1,mp.inf])
    expected=(mp.exp(-uv)+1j)/(mp.exp(-uv)-1j)
    assert abs(1j*mp.exp(1j*integral)-expected)<mp.mpf('1e-25')
print('PASS Breather Fourier integral gives the classical amplitude',flush=True)
