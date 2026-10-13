"""Sector determinants of the closed forms in the chapter's projected labels n = N + s (sum 0), symbolic q, X, s."""
import sympy as sp
q, X = sp.symbols('q X', positive=True)
for h in (3, 4, 5):
    s = sp.symbols(f's1:{h}'); s = list(s) + [-sum(s)]
    Q = [sp.Rational(1) - sp.Rational(1, h)] + [-sp.Rational(1, h)]*(h - 1)       # a sector of V (x) W: Q = mu_1 + N + s, N = 0
    Q = [Q[c] + s[c] for c in range(h)]
    mu = lambda i: [sp.Integer(1 if c == i else 0) - sp.Rational(1, h) for c in range(h)]
    MV = sp.zeros(h); MVb = sp.zeros(h); MP = sp.zeros(h); MR = sp.zeros(h)
    Qb = [-x for x in Q]                                                           # a sector of Vbar (x) W
    for i in range(h):
        n = [Q[c] - mu(i)[c] for c in range(h)]          # V: n = Q - mu_i
        MV[i, i] = q**(-n[i]); MV[(i - 1) % h, i] = X
        nb = [Qb[c] + mu(i)[c] for c in range(h)]        # Vbar: n = Q + mu_i
        pw = 1
        for j in range(h):
            b = (i + j) % h; pw = pw*q**nb[b]; MVb[b, i] = X**j*pw
        MP[i, i] = q**(-nb[i]); MP[(i - 1) % h, i] = X*q**(nb[i] - nb[(i - 1) % h] - 1)     # eq-dt-tplus, X = X_+
        MR[i, i] = q**(nb[i]); MR[(i + 1) % h, i] = X*q**(nb[i] - nb[(i + 1) % h])          # reversed orientation, Vbar
    r = h - 1
    for name, M, guess in (('M_V  (I-)', MV, q**r*(1 + (-1/q)**r*X**h)),
                           ('M_Vb (I-)', MVb, None),
                           ('M+_Vb(I+)', MP, q**(-r)*(1 + (-q)**r*X**h)),
                           ('M-_Vb(rev)', MR, None)):
        d = sp.factor(sp.simplify(sp.powsimp(sp.expand(M.det()), force=True)))
        print(f'h={h} {name}: det = {d}' + (f'   matches guess: {sp.simplify(d - guess) == 0}' if guess is not None else ''))
