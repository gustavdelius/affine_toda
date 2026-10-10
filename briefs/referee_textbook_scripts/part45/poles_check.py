import mpmath as mp
def blk(x,H): return lambda th: mp.sinh(th/2+1j*mp.pi*x/(2*H))/mp.sinh(th/2-1j*mp.pi*x/(2*H))
def curly(x,H,B):
    f=[blk(x-1,H),blk(x+1,H),blk(x-1+B,H),blk(x+1-B,H)]
    return lambda th: f[0](th)*f[1](th)/(f[2](th)*f[3](th))
def order_at(S,th0):
    # estimate pole order from |S| scaling
    e1,e2=1e-4,1e-5
    a=abs(S(th0+e1)); b=abs(S(th0+e2))
    return mp.log(b/a)/mp.log(e1/e2)
# c_2: S_22 = {1}{H-1}{3}{H-3}, H=4+B
B=mp.mpf('0.37'); H=4+B
fs=[curly(1,H,B),curly(H-1,H,B),curly(3,H,B),curly(H-3,H,B)]
S=lambda th: fs[0](th)*fs[1](th)*fs[2](th)*fs[3](th)
for x in [B,2,4,H-4,H-2,H-B]:
    print("c_2 S_22 at theta = i pi*%s/H : order ~ %.3f"%(mp.nstr(x,4), float(order_at(S,1j*mp.pi*x/H))))
# a_4: S_34 = {2}{4}{6}, h=5
h=5; B=mp.mpf('0.37')
fs=[curly(p,h,B) for p in (2,4,6)]
S=lambda th: fs[0](th)*fs[1](th)*fs[2](th)
for x in [1,3,5,7-5+0.0]:
    print("a_4 S_34 at theta = i pi*%s/5 : order ~ %.3f"%(x, float(order_at(S,1j*mp.pi*x/h))))
print("c_2 S_22: -i*Res at theta=2 pi i/H, as function of B")
for B in [0.05,0.3,0.6,1.0,1.4,1.7,1.95]:
    B=mp.mpf(B); H=4+B
    fs=[curly(1,H,B),curly(H-1,H,B),curly(3,H,B),curly(H-3,H,B)]
    S=lambda th: fs[0](th)*fs[1](th)*fs[2](th)*fs[3](th)
    th0=1j*mp.pi*2/H; e=mp.mpf('1e-8')
    res=e*S(th0+e)
    print("  B=%.2f  -i Res = %s"%(float(B), mp.nstr(-1j*res,6)))
