from dn_scalar import *
for n in (4, 5, 6):
    for om in (0.37, 1.31):
        T = (n - 1)*om
        for mu in (0.13 + 0.21j, -0.41 + 0.33j):
            u = abs(F(mu, n, om)*F(-mu, n, om) - 1)
            c = abs(F(T - mu, n, om)/F(mu, n, om)/r_cross(mu, n, om) - 1)
            ua = abs(A(mu, n, om)*A(-mu, n, om) - 1)
            ca = abs(A(T/2 - mu, n, om)/A(T/2 + mu, n, om)/F(2*mu, n, om) - 1)
            print(n, om, mu, f'F unit {u:.1e} cross {c:.1e}   A unit {ua:.1e} cross {ca:.1e}')
        print('   F(0+) =', F(1e-7, n, om))
