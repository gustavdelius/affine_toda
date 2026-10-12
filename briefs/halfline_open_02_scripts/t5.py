from dn_scalar import *
for n in (4, 5, 6):
    for om in (0.37, 1.31, 2.2):
        T = (n - 1)*om
        for mu in (0.13 + 0.21j, -0.41 + 0.33j):
            ua = abs(A(mu, n, om)*A(-mu, n, om) - 1)
            ca = abs(A(T/2 - mu, n, om)/A(T/2 + mu, n, om)/F(2*mu, n, om) - 1)
            kb = KB_num(mu, n, om)/KB_blocks(mu, n, om)
            print(n, om, mu, f'A unit {ua:.1e} cross {ca:.1e}  KB/blocks = {kb:.8f}')
