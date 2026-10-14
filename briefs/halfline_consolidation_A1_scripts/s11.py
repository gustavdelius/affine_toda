import sys
sys.path.insert(0, '/home/gustav/Git/affine_toda/textbook_code')
src = open('/home/gustav/Git/affine_toda/textbook_code/particle_reflection.py').read()
exec(src.split('# ---------------------------------------------------------------- 1. blocks')[0])
exec(src.split('def cdrs(a, n, tilde):')[0].split('# ---------------------------------------------------------------- 1. blocks')[0])
ns = {}
code = 'def cdrs(a, n, tilde):' + src.split('def cdrs(a, n, tilde):')[1].split('\n\n\nreport')[0]
exec(code)
def Kminus(a, n):
    h = n + 1
    return Prod(h, sign=-1).mul((a, 0)).mul((h - a, 0))*K_uniform(a, n)
for n in range(1, 8):
    print(n, [close(Kminus(a, n).dual(), cdrs(a, n, False)) for a in range(1, n + 1)],
          [close(Kminus(a, n).dual(), lambda t, B, a=a: 1/cdrs(a, n, False)(t, B)) for a in range(1, n + 1)])
# CDRS 5.16 untilded: does it satisfy the axioms?
print([all(axioms(lambda a, n=n: cdrs(a, n, False), n, 'c')) for n in range(1, 6)])
