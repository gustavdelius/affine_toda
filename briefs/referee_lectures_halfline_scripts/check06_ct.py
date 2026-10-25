import mpmath as m
m.mp.dps=65
B=m.mpf(9)/20;h=3

def block(x,t):return m.sinh(t/2+1j*m.pi*x/(2*h))/m.sinh(t/2-1j*m.pi*x/(2*h))
def K(n,k,t):
    xs=[m.mpf('.5')-(2*n+1)*B/2,m.mpf('2.5')+(2*k+1)*B/2,
        m.mpf('.5')-(2*n-1)*B/2,m.mpf('2.5')+(2*k-1)*B/2,
        m.mpf('1.5')+(2*n-1)*B/2,m.mpf('1.5')-(2*k-1)*B/2,
        -m.mpf('2.5')+(2*n+1)*B/2,-m.mpf('.5')-(2*k+1)*B/2]
    return block(3-B/2,t)/block(2,t)/block(1-B/2,t)*m.fprod(block(x,t) for x in xs)
def S(t):return block(2,t)/block(B,t)/block(2-B,t)
p=m.mpf('.5')+3*B/2; r=3*B/2-m.mpf('.5');v=m.mpf('1.5')-3*B/2
for eps in [m.mpf('1e-10'),m.mpf('1e-20')]:
    R=eps*K(-2,2,1j*m.pi*p/3+eps)
    g2=2*eps*K(-1,1,1j*m.pi*v/3+eps)/1j
    G2=eps*S(2j*m.pi/3+eps)/1j
    der=K(1,-1,1j*m.pi*r/3+eps)/eps
    rhs=g2*G2*der/2
    print('eps',eps,'-i Res',m.nstr(R/1j,20),'identity residual',m.nstr(abs(R-rhs),6),'zero derivative',m.nstr(der,14))
    assert abs(R-rhs)<eps*100
print('PASS 6: simple pole, negative residue, simple reflection zero, residue identity')
