import numpy as np
m=1.;a=np.sin(np.pi/3);th=.6;s=a/np.cosh(th);k=2*m*a*np.cosh(th);U=np.exp(.2+.7j)
for C in [1,-1]:
    D=1+C*s
    j=np.arange(3);z=np.exp(2j*np.pi*j/3)*U
    tau=(2+z+1/z)/D
    direct=k*((z+1/z)/D+2*(1-s*s)/D**2)
    printed=k/D*(z+1/z+2*(1-C*s)/D)
    corrected=k/D*(z+1/z+2*(1-C*s))
    # algebraic square roots continued in j, not principal numerical sqrt.
    psi=np.log(U)/(2j);cos=np.cos(psi+np.pi*j/3)
    rootprod=4*(cos*cos-a*a)/D
    ratios=(direct-m*C*rootprod)/tau
    print('10 C',C,'printed derivative error',np.max(abs(direct-printed)),'correct derivative error',np.max(abs(direct-corrected)),'constant ratios',ratios)
    assert np.max(abs(ratios-(k-m*C)))<1e-12
print('10 correct k=',k,'main printed 2ms=',2*m*s,'exercise printed 2m*a/s=',2*m*a/s)
print('PASS 10: boundary identity with corrected derivative and k*s=2m*sin^2(pi*a/h)')
