import mpmath as m
m.mp.dps=45
S=lambda t:(m.sinh(t)+1j*m.sin(2*m.pi/3))/(m.sinh(t)-1j*m.sin(2*m.pi/3))
b=lambda x,t:m.sinh(t/2+1j*m.pi*x/12)/m.sinh(t/2-1j*m.pi*x/12)
K=lambda t:b(1,t)*b(3,t)/b(4,t)
e=m.mpf('1e-20'); G2=e*S(2j*m.pi/3+e)/1j;g2=2*e*K(1j*m.pi/2+e)/1j
print('22 Gamma^2',m.nstr(G2,20),'g^2',m.nstr(g2,20));assert abs(G2+2*m.sqrt(3))<m.mpf('1e-18')
F=lambda n,d:m.sin(n*d)*m.cos(m.pi/6-n*d)/m.sin(d)
d=m.pi/12
print('24 degeneracy',m.nstr(F(2,d)+F(-2,d),20),m.nstr(F(1,d),20))
assert abs(F(2,d)+F(-2,d)-F(1,d))<m.mpf('1e-40')
for n in range(1,8):
    d=m.pi*m.sqrt(2)/30
    err=abs(sum(m.cos(m.pi/6-(2*k+1)*d) for k in range(n))-F(n,d));assert err<m.mpf('1e-40')
print('PASS 22,24: independent residues, energy sums, exact B=1/2 degeneracy')
