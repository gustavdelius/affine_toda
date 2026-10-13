import sys
exec(open('t4b.py').read().split('paths = {')[0])
# paths moving O(e) along the singular curve but with tau_2(0), tau_4(0) = O(e^2)
for sgn in ('1', '-1'):
    for c in ('1', '10'):
        out = []
        for e in ('1e-2', '1e-3', '1e-4'):
            ee = mp.mpf(e)
            d1 = d1s*(1 + mp.mpf(sgn)*ee); d2c = d2_of(d1)
            bg = bg_of(d1, d2c*(1 - mp.mpf(c)*ee**2))
            if not regular(bg):
                out.append('irregular'); continue
            opn, K, Kf, res = kmatrix(bg, P, Om)
            out.append(f'tau2(0)={mp.nstr(bg.tau[2](0).real,2)} K11={mp.nstr(K[(1,1)],7)} |K41|={mp.nstr(abs(K[(4,1)]),6)}')
        print(f'd1=d1*(1{"+" if sgn=="1" else "-"}e) on curve, d2 lowered by {c}e^2:', ' | '.join(out), flush=True)
