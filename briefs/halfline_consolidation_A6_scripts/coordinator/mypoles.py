import numpy as np, cmath
exec(open('smatrix_simply_laced.py').read().split('\nB = 0.37')[0])
for n in (4, 5, 6):
    h, blocks, pf, orbits = setup('d', n)
    m1 = pf.min()   # mass units: pf/pf.min(); m_2 from the vector fusion
    vec = [a for a in range(n)]
    for a in range(n):
        f = lambda th: np.prod([S(blocks[(a, b)], h, 1.0, th) for b in range(n)])
        pl = []
        for y in np.arange(1, h//2 + 1) if h % 2 == 0 else []:
            if y >= h/2: continue
            u = np.pi*y/h
            v1, v2 = abs(f(1j*(u + 1e-6))), abs(f(1j*(u + 1e-3)))
            if v1 > 1e3*v2/1e3 and v1 > 1e4: pl.append(y)
        print(f'd_{n} label {a} mass {pf[a]:.4f}: poles at theta = i pi y/h, y =', pl,
              ' p/(m_lightest...) =', [round(pf[a]*np.sin(np.pi*y/h), 6) for y in pl])
    print('   masses', np.round(pf, 4), ' 2 m_v cos(pi/h)/2 =', round(np.cos(np.pi/h)*pf[0], 6))
