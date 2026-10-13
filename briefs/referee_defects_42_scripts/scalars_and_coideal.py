"""Additional independent scalar, exact coideal, and source-convention checks."""
import mpmath as m
import numpy as np
import sympy as s
from independent_checks import ls_test

def scalar_products():
    m.mp.dps=35
    def extrap(vs,ns):
        A=m.matrix([[m.mpf(1)/n**j for j in range(len(ns))] for n in ns])
        return m.lu_solve(A,m.matrix(vs))[0]
    for coupling,th,td in (('1.37',m.mpf('.33'),'.21'),('1.37',m.mpc('.21','.17'),'.21'),('.73',m.mpc('.29','.11'),'-.38')):
        l=m.mpf(coupling);td=m.mpf(td);v=l*td;u=-1j*th
        def logrterm(k,z):
            return (m.loggamma(k*l+m.mpf('.25')-z)+m.loggamma((k+1)*l+m.mpf('.75')-z)
                   -m.loggamma((k+m.mpf('.5'))*l+m.mpf('.25')-z)-m.loggamma((k+m.mpf('.5'))*l+m.mpf('.75')-z))
        def sigterm(k,z):
            a=1j*v/m.pi
            return (m.loggamma(m.mpf('.5')+a+(2*k-1)*l-z)+m.loggamma(m.mpf('.5')-a+(2*k-1)*l-z)
                   -m.loggamma(m.mpf('.5')-a+(2*k-2)*l-z)-m.loggamma(m.mpf('.5')+a+2*k*l-z))
        z1=1j*l*(th-td)/(2*m.pi);z2=1j*l*(-th-td)/(2*m.pi);zu=l*u/m.pi
        # Strip common R0 sigma(eta,u); the ratio is the full fusion scalar identity.
        pref=(1+m.exp(-2*v))/(1+m.exp(-2*l*th-2*v))*(1+1j*m.exp(l*(-th-td)))/(1+1j*m.exp(l*(th-td)))
        sigpref=m.cosh(v)/m.cos(1j*v+l*u)
        total=m.mpc(0);ns=[40,80,160,320,640];values=[]
        for k in range(640):
            total+=logrterm(k,z1)-logrterm(k,-z1)-logrterm(k,z2)+logrterm(k,-z2)-sigterm(k+1,zu)+sigterm(k+1,-zu)
            if k+1 in ns:values.append(pref/sigpref*m.exp(total))
        estimate4=extrap(values[-4:],ns[-4:]); estimate5=extrap(values,ns)
        print('scalar coupling/theta/td:',coupling,m.nstr(th),m.nstr(td),flush=True)
        print('  cutoff residuals:',[m.nstr(abs(z-1),5) for z in values],flush=True)
        print('  extrap4,extrap5,stability:',*[m.nstr(z,7) for z in (abs(estimate4-1),abs(estimate5-1),abs(estimate4-estimate5))],flush=True)

def exact_coideal():
    # PBW letters (kind,Cartan exponent). X has at most ONE e/f in total.
    cart=[('1',a) for a in (-1,0,1)];other=[('e',0),('f',0)]
    B=[(g,h) for g in cart for h in cart]+[(g,h) for g in other for h in cart]+[(g,h) for g in cart for h in other]
    eps=lambda g: int(g[0]=='1')
    def rep(g):
        kind,a=g; k=s.diag(2,s.Rational(1,2))**a
        X=s.eye(2) if kind=='1' else (s.Matrix([[0,1],[0,0]]) if kind=='e' else s.Matrix([[0,0],[1,0]]))
        return X*k
    M=s.Matrix.hstack(*[s.Matrix(rep(g)*eps(h)-eps(g)*rep(h)).reshape(4,1) for g,h in B])
    S=s.Matrix.hstack(*M.nullspace()); dims=[S.cols]
    def delta(g):
        kind,a=g
        if kind=='1':return [(('1',a),('1',a))]
        if kind=='e':return [(('e',a),('1',a)),(('1',a+1),('e',a))]
        return [(('f',a),('1',a-1)),(('1',a),('f',a))]
    maps={}
    for col,(g,h) in enumerate(B):
        for g1,g2 in delta(g):
            for h1,h2 in delta(h):
                key=(g1,h2)
                if key not in maps:maps[key]=s.zeros(21)
                maps[key][B.index((g2,h1)),col]+=1
    while True:
        C=s.Matrix.hstack(*S.T.nullspace()).T
        restrictions=s.Matrix.vstack(*[C*F*S for F in maps.values()])
        T=s.Matrix.hstack(*restrictions.nullspace())
        if T.cols==S.cols:break
        S=S*T;dims.append(S.cols)
    targets=s.zeros(21,5)
    for j,a in enumerate((-1,0,1)):targets[B.index((('1',a),('1',a))),j]=1
    targets[B.index((('e',0),('1',0))),3]=1;targets[B.index((('1',1),('e',0))),3]=1
    targets[B.index((('f',0),('1',-1))),4]=1;targets[B.index((('1',0),('f',0))),4]=1
    print('exact rational coideal dimensions:',dims,'final',S.cols,'span equality',S.row_join(targets).rank()==5,flush=True)

def cz():
    # Equation (2.18) in SOURCE conventions, before book translation.
    x,q,r,eps=s.symbols('x q r eps',nonzero=True)
    Q=s.symbols('Q',nonzero=True); v=s.symbols('v',nonzero=True)
    alpha=4
    T=lambda x:s.Matrix([[Q**alpha/v,eps*x*Q**(alpha+2)],[eps*x*Q**(-alpha),v*Q**(-alpha-2)]])
    KD=s.diag(r*x+1/(r*x),r/x+x/r)
    F=s.simplify(T(1/x).inv()*KD*T(x))
    plus=x*(r*q-eps**2/r/q)+(r**-1*q-eps**2*r/q)/x
    minus=x*(q/r-eps**2*r/q)+(r*q-eps**2/r/q)/x
    # Source similarity removes Q^alpha; diagonal entries and offdiag product are gauge invariant.
    inv=[s.factor((F[0,0]/F[1,1]-plus/minus).subs(q,Q**-2)),
         s.factor((F[0,1]*F[1,0]/F[0,0]/F[1,1]-(x*x-x**-2)**2*eps**2/(plus*minus)).subs(q,Q**-2))]
    print('CZ (2.18) source-convention invariant differences:',inv,flush=True)
    d=s.symbols('d',nonzero=True)
    print('CZ map eps=-i q/d diagonal coefficient identity:',s.simplify(plus.subs(eps,-s.I*q/d)-q*(x*(r+d**-2/r)+(1/r+d**-2*r)/x))==0,flush=True)
    # Generic type-II numerical blocks, constructed directly from the SOURCE displayed T_II.
    qq=np.exp(.71j); QQ=np.exp(-.355j); bp=.37+.19j;bm=-.28+.13j; rr=np.exp(.23j);xx=np.exp(.41)
    def TT(x,alpha):
        A=lambda a:QQ**a-bm*bp.conjugate()*QQ**(-a)*x*x
        D=lambda a:-bp*bm.conjugate()*QQ**(a-4)*x*x+QQ**(-a)
        B=lambda a:x*(bp*QQ**a+bm*QQ**(-a))
        C=lambda a:-x*(bm.conjugate()*QQ**(a-4)+bp.conjugate()*QQ**(-a))
        return np.array([[A(alpha),C(alpha+2)],[B(alpha),D(alpha+2)]])
    kd=np.diag([rr*xx+1/rr/xx,rr/xx+xx/rr])
    ratios=[]
    for a in (0,2,4,6):
        ff=np.linalg.solve(TT(1/xx,a),kd@TT(xx,a))
        ratios.append(ff[0,0]/ff[1,1])
    print('CZ type-II charge-dependent diagonal ratios:',ratios,flush=True)

if __name__=='__main__':
    exact_coideal();cz()
    for kind,tw,direction in (('f',0,-1),('e',1,1)):
        print('I-plus CORRECT limit',kind,tw,direction,[ls_test(5,1,kind,tw,td=direction*t) for t in (2,4,6)],flush=True)
    scalar_products()
