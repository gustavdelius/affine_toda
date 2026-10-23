import mpmath as mp
exec(open('t4b.py').read().split('th = mp.mpf')[0].replace('mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 60', 'mp.mp.dps = 50'))
print('tau_i\'/tau_i(0), i=0,1,3:', [mp.nstr(bgs.tau[i](0, 1).real/bgs.tau[i](0).real, 12) for i in (0, 1, 3)])
print('m1+m2-sqrt5 =', mp.nstr(T.mass(1) + T.mass(2) - mp.sqrt(5), 12))
# residue of tau_2 zero: tau_2'(0), and second-order behaviour
print('tau_2\'(0), tau_4\'(0):', mp.nstr(bgs.tau[2](0, 1).real, 6), mp.nstr(bgs.tau[4](0, 1).real, 6))
src = open('/home/gustav/Git/affine_toda/code/particle_reflection.py').read()
PR = {}; exec(src.split('\n# ---------------------------------------------------------------- 1. blocks')[0], PR)
print('DG b=2 K_1(0.6), a_4:', mp.nstr((PR['K_uniform'](1, 4)*PR['M_sol'](1, 2, 4))(mp.mpf('0.6'), 0), 10))
