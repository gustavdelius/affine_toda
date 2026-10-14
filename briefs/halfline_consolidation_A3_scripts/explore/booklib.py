"""Exploration helper: load definitions (not checks) from the book's textbook_code modules.

boundary_classical.py and boundary_images.py run their checks at import; we exec only the definition blocks.
"""
import os

TC = '/home/gustav/Git/affine_toda/textbook_code'


def _between(src, start, stop):
    i = 0 if start is None else src.index(start)
    j = src.index(stop, i + 1)
    return src[i:j]


def load_classical():
    src = open(os.path.join(TC, 'boundary_classical.py')).read()
    ns = {}
    code = _between(src, None, '\n# check the labels')
    code += _between(src, '# ---------------------------------------------------------------- 2. K-matrix compatibility',
                     '\nsigns = lambda r')
    code += _between(src, '# ---------------------------------------------------------------- 8. Sklyanin',
                     '\npairs = [(0.7, 1.3)')
    exec(compile(code, 'boundary_classical.py[defs]', 'exec'), ns)
    return ns


def load_images():
    src = open(os.path.join(TC, 'boundary_images.py')).read()
    ns = {}
    code = _between(src, None, '\n# ---------------------------------------------------------------- 1. boundary conditions')
    exec(compile(code, 'boundary_images.py[defs]', 'exec'), ns)
    return ns


if __name__ == '__main__':
    import time
    t0 = time.time(); C = load_classical(); I = load_images()
    print('loaded in %.2fs' % (time.time() - t0))
    print(sorted(k for k in C if not k.startswith('_'))[:60])
    print(sorted(k for k in I if not k.startswith('_')))
    print('beta, m (classical):', C['beta'], C['m'], ' (images):', I['beta'], I['m'])
