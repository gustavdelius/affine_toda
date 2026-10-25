#!/usr/bin/env python3
"""Independent algebra checks for Gemini defects lectures 4 and 5.

Sparse actions are composed on unbounded integer labels; no truncation and no
imports of the book's or referee's implementations are used.
"""
import cmath
import math
import sympy as S


def clean(v):
    return {k: S.simplify(c) for k, c in v.items() if S.simplify(c) != 0}


def act(op, v):
    out = {}
    for k, c in v.items():
        for l, d in op(k).items():
            out[l] = out.get(l, 0) + c*d
    return clean(out)


def equal(a, b):
    return not clean({k: a.get(k, 0)-b.get(k, 0) for k in a.keys() | b.keys()})


def check(name, ok):
    print(('PASS' if ok else 'FAIL') + ' ' + name)
    assert ok, name


q = S.Rational(2)
y, yd = S.Rational(3), S.Rational(5)
kappa = q-1/q
mu = S.I/kappa


def sg_delta(j, state, offset, opposite=False, fock=False):
    a, n = state
    # a=+1 is v_+; k_1=q^a, k_0=q^-a.
    power = 1 if j == 1 else -1
    shift = -1 if j == 1 else 1
    osc = (1-q**(2*n)) if fock and j == 1 else 1
    out = {(a, n+shift): mu/yd*osc*(q**(-power*a) if opposite else 1)}
    if a == power:
        out[(-a, n)] = 1/y*(1 if opposite else q**(-power*(2*n+offset)))
    return clean(out)


def sg_L(state, offset, fock=False):
    a, n = state
    Y = y/yd
    if a == 1:
        diag = q**(-n)
        if fock:
            diag *= 1+Y**2*q**(2*n+2)
        return {(1, n): diag, (-1, n+1): -S.I*Y*q**(-n)}
    down = (1-q**(2*n)) if fock else 1
    return {(1, n-1): -S.I*Y*q**(n+offset)*down, (-1, n): q**(n+offset)}


for fock in [False, True]:
    ok = True
    for offset in [0, 1, 3]:
        for n in range(0 if fock else -3, 4):
            for a in [1, -1]:
                for j in [0, 1]:
                    v = {(a, n): 1}
                    lhs = act(lambda k: sg_L(k, offset, fock), act(lambda k: sg_delta(j, k, offset, fock=fock), v))
                    rhs = act(lambda k: sg_delta(j, k, offset, True, fock), act(lambda k: sg_L(k, offset, fock), v))
                    ok &= equal(lhs, rhs)
    check(('Fock' if fock else 'Weyl')+' L intertwines both Borel generators on unbounded labels', ok)

# The product a_+ a_- depends on n in the Fock module. Z is normal, not central.
for fock in [False, True]:
    u = lambda n: (1-q**(2*n) if fock else 1)/kappa**2
    check(('Fock' if fock else 'Weyl')+' oscillator relation and normal element', all(
        S.simplify(q*u(n)-u(n+1)/q-1/kappa)==0 and
        S.simplify(u(n+1)-1/kappa**2-q**2*(u(n)-1/kappa**2))==0
        for n in range(4)))

# Serre coefficients on each Fock basis vector; no boundary is truncated.
def osc(j, state):
    n = state
    return clean({n+(1 if j == 0 else -1): (1 if j == 0 else 1-q**(2*n))/kappa})
ok = True
for j in [0, 1]:
    for n in range(5):
        out = {}
        for word, c in [([j,j,j,1-j],1),([j,j,1-j,j],-(q*q+1+q**-2)),([j,1-j,j,j],q*q+1+q**-2),([1-j,j,j,j],-1)]:
            v = {n: c}
            for gen in reversed(word):
                v = act(lambda k: osc(gen, k), v)
            for k, d in v.items():
                out[k] = out.get(k, 0)+d
        ok &= not clean(out)
check('Both cubic q-Serre identities on the Fock module', ok)

# Crossing uses a transpose on the particle indices only, not on a charge block.
Q, Y, t = S.symbols('Q Y t', nonzero=True)
ok = True
for n in [-3, 0, 2]:
    for offset in [0, 1]:
        M = S.Matrix([[Q**(-n), -S.I*Y*Q**(n+1+offset)],[-S.I*Y*Q**(-n),Q**(n+1+offset)]])
        # L(-theta)^-1 / scalar, compared to conjugated partial transpose at crossing.
        shiftedY = -Y/Q
        crossed = S.Matrix([[Q**(n+offset),-S.I*shiftedY*Q**(n+1+offset)],[-S.I*shiftedY*Q**(-n),Q**(-n-1)]])
        ok &= S.simplify(M.inv()*Q**offset*(1+Y**2)-crossed)==S.zeros(2)
check('All four crossing entries reduce to Phi Phi_shift (1+Y^2)=q^-s', ok)

qphys = -cmath.exp(-1j*math.pi*math.sqrt(2))
import numpy as np
M = np.array([[qphys**-2,-1j*.7*qphys**4],[-1j*.7*qphys**-2,qphys**4]])
check('Hermitian unitarity for irrational physical coupling', np.max(np.abs(M.conj().T@M-1.49*np.eye(2)))<1e-12)

# Neutral fusion on full spin space, with symbolic q and Y1/Y2.
qsave, ysave, ydsave = q, y, yd
q = Q
Y1,Y2 = S.symbols('Y1 Y2')
def leg(slot, YY, state):
    a,b,n = state
    spin = (a,b)[slot]
    out = {}
    for (newspin, newn), c in ({(1,n):q**-n,(-1,n+1):-S.I*YY*q**-n} if spin==1 else {(1,n-1):-S.I*YY*q**n,(-1,n):q**n}).items():
        spins=[a,b]; spins[slot]=newspin
        out[(*spins,newn)] = c
    return out
okket=okquot=True
for n in [-2,0,3]:
    v={(1,-1,n):1,(-1,1,n):-1}
    fused=act(lambda k:leg(0,Y1,k),act(lambda k:leg(1,Y2,k),v))
    target={k:c*(1+q*Y1*Y2) for k,c in v.items()}
    okket &= equal({k:c.subs(Y1,q*Y2) for k,c in fused.items()}, {k:c.subs(Y1,q*Y2) for k,c in target.items()})
    for a,b in [(1,1),(1,-1),(-1,1),(-1,-1)]:
        fused=act(lambda k:leg(0,Y1,k),act(lambda k:leg(1,Y2,k),{(a,b,n):1}))
        projected={}
        for (aa,bb,nn),c in fused.items():
            if aa != bb:
                projected[nn]=projected.get(nn,0)+aa*c
        expected={n:a*(1+q*Y1*Y2)} if a != b else {}
        okquot &= equal({k:c.subs(Y1,Y2/q) for k,c in projected.items()}, {k:c.subs(Y1,Y2/q) for k,c in expected.items()})
check('Neutral fusion invariant ket at Y1/Y2=q',okket)
check('Neutral fusion quotient covector at Y1/Y2=q^-1',okquot)
q,y,yd=qsave,ysave,ydsave

# Higher-rank lattice module and intertwiners, written directly as state actions.
mu = S.Rational(11,7)
X = -(q-1/q)*mu*y/yd

def shifted(N, i, j):
    out=list(N);out[i]+=1;out[j]-=1
    return tuple(out)

def ar_delta(j, state, conjugate=False, opposite=False):
    i,N=state; h=len(N); jp=(j+1)%h
    sign=-1 if conjugate else 1
    vk=q**(sign*((i==j)-(i==jp)))
    wk=q**(N[j]-N[jp])
    out={(i,shifted(N,jp,j)):mu/yd*q**N[jp]*(1/vk if opposite else 1)}
    if i==(jp if conjugate else j):
        out[(j if conjugate else jp,N)]=1/y*(1 if opposite else 1/wk)
    return clean(out)

def ar_L(state,conjugate=False):
    i,N=state;h=len(N)
    if not conjugate:
        return {(i,N):q**(-N[i]),((i-1)%h,shifted(N,i,(i-1)%h)):X}
    return {((i+j)%h,shifted(N,(i+j)%h,i)):X**j*q**sum(N[(i+k)%h] for k in range(j+1)) for j in range(h)}
for conjugate in [False,True]:
    ok=True
    for h in [3,4,5]:
        for N in [tuple([0]*h),tuple([2,-1,-1]+[0]*(h-3))]:
            for i in range(h):
                for j in range(h):
                    v={(i,N):1}
                    lhs=act(lambda k:ar_L(k,conjugate),act(lambda k:ar_delta(j,k,conjugate),v))
                    rhs=act(lambda k:ar_delta(j,k,conjugate,True),act(lambda k:ar_L(k,conjugate),v))
                    ok &= equal(lhs,rhs)
    check(('Antivector' if conjugate else 'Vector')+' a_r transmission intertwines for h=3,4,5',ok)

# A control catches the original lecture's incompatible Cartan matrix.
old_k = S.diag(1/q,q)
f1 = S.Matrix([[0,0],[1/y,0]])
check('Original bulk Cartan convention fails its declared f1 relation', old_k*f1*old_k.inv() != q**-2*f1)

# The adjacent-generator clock relation is computed from lattice shifts.
ok=True
for h in [3,4,5]:
    N=tuple([2,-1,-1]+[0]*(h-3))
    for j in range(h):
        jp=(j+1)%h
        def f(gen,state):
            end=(gen+1)%h
            return {shifted(state,end,gen):mu/yd*q**state[end]}
        a=act(lambda k:f(j,k),act(lambda k:f(jp,k),{N:1}))
        b=act(lambda k:f(jp,k),act(lambda k:f(j,k),{N:1/q}))
        ok &= equal(a,b)
check('Root-lattice clock relation for h=3,4,5',ok)

# Independent Gamma-product evaluation with the explicit common cutoff.
from scipy.special import loggamma

def log_phi_cutoff(u, lam, count):
    k=np.arange(count,dtype=float)
    args=[k*lam+.25,(k+1)*lam+.75,(k+.5)*lam+.25,(k+.5)*lam+.75]
    z=1j*lam*u/(2*np.pi)
    logratio=sum(sign*np.sum(loggamma(a-z)-loggamma(a+z)) for a,sign in zip(args,[1,1,-1,-1]))
    return .25j*np.pi*(1+lam)-np.log(1+1j*np.exp(lam*u))+logratio

def phi_product(u,lam):
    # Richardson eliminates the leading 1/N tail, not the physical phase.
    return np.exp(2*log_phi_cutoff(u,lam,4096)-log_phi_cutoff(u,lam,2048))

lam=math.sqrt(2)+.5
u=.41
cross=phi_product(u,lam)*phi_product(u+1j*np.pi,lam)*(1+np.exp(2*lam*u))
unit=abs(phi_product(u,lam))**2*(1+np.exp(2*lam*u))
check('Explicit Gamma product phase gives scalar crossing and unitarity',abs(cross-1)<2e-6 and abs(unit-1)<1e-9)
plus=1j*np.pi/(2*lam);minus=-plus
eps=1e-5
regular=phi_product(plus+eps,lam)/phi_product(plus+2*eps,lam)
pole=phi_product(minus+eps,lam)/(2*phi_product(minus+2*eps,lam))
check('Gamma-factor upper-pole cancellation and lower simple-pole scaling',abs(regular-1)<2e-4 and abs(pole-1)<2e-4)
check('Classical logarithmic-delay series',all(abs(2*sum(np.exp(-(2*j+1)*u)/(2*j+1) for j in range(100))-np.log(1/np.tanh(u/2)))<1e-13 for u in [.4,1.1,2.3]))
