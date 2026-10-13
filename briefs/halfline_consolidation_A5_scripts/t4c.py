import sys
exec(open('t4b.py').read().split('paths = {')[0])
for cfac in ('0.5', '-0.5', '2', '-2'):
    for e in ('1e-2', '1e-3', '1e-4', '1e-6'):
        ee = mp.mpf(e)
        bg = bg_of(d1s, d2s*(1-ee), chi1=9*mp.pi/5 + mp.mpf(cfac)*ee)
        if not regular(bg):
            print(f'cfac {cfac} e={e}: not regular'); continue
        opn, K, Kf, res = kmatrix(bg, P, Om)
        print(f'cfac {cfac} e={e}: tau2(0)={mp.nstr(bg.tau[2](0).real,3)} tau4(0)={mp.nstr(bg.tau[4](0).real,3)} K11={mp.nstr(K[(1,1)],8)} K44={mp.nstr(K[(4,4)],8)} |K41|={mp.nstr(abs(K[(4,1)]),6)} |K14|={mp.nstr(abs(K[(1,4)]),6)} unit={mp.nstr(unitarity_defect(opn,Kf),2)}', flush=True)
