"""Pure-rename step for inserting the chapter on superalgebra Toda theories as
chapter 16 (end of Part IV).  Run once from the repository root on a clean tree:

    python3 briefs/superalgebra_toda_renumber.py

It renames chapters 16-50 to 17-51 (git mv), the chapter-numbered figures and
figure scripts chNN_* (NN >= 16) to ch(NN+1)_*, and fixes only the paths that
would otherwise break: the chapter list in _quarto.yml, the image links and
script paths in the renamed chapters, the output names inside the renamed
figure scripts, the figure-script rows of code/README.md, the path in
code/sectors/bloch/README.md, and the current-file columns of briefs/README.md
(its dated history paragraph keeps the old numbers).  It then stages exactly
these files, so that the index holds the pure-rename commit.  Prose and
hard-coded chapter numbers are left for the following text commit.

The mapping is written to briefs/superalgebra_toda_renames.tsv.
"""
import pathlib
import re
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[1]
DIRS = ['part5-semiclassics', 'part6-soliton-smatrices', 'part7-half-line', 'part8-defects',
        'part9-foundations', 'outlook']


def git(*args):
    subprocess.run(['git', *args], cwd=ROOT, check=True)


chapters = []
for d in DIRS:
    for f in sorted((ROOT / d).glob('[0-9][0-9]-*.qmd')):
        n = int(f.name[:2])
        if 16 <= n <= 50:
            chapters.append((f.relative_to(ROOT).as_posix(), f'{d}/{n + 1:02d}{f.name[2:]}'))
figs = []
for f in sorted(ROOT.glob('part*/figures/ch*_f*.svg')) + sorted(ROOT.glob('code/figures/ch*_f*.py')):
    m = re.match(r'ch(\d+)_(f\d+_.*)', f.name)
    if m and 16 <= int(m.group(1)) <= 50:
        new = f.with_name(f'ch{int(m.group(1)) + 1}_{m.group(2)}')
        figs.append((f.relative_to(ROOT).as_posix(), new.relative_to(ROOT).as_posix()))

renames = chapters + figs
# rename from the highest number down, so that no target exists yet
for old, new in sorted(renames, key=lambda p: p[0], reverse=True):
    git('mv', old, new)
(ROOT / 'briefs' / 'superalgebra_toda_renames.tsv').write_text(
    'old\tnew\n' + ''.join(f'{o}\t{n}\n' for o, n in renames))

# path substitutions: full chapter paths and figure basenames (without extension)
subs = [(o, n) for o, n in chapters]
subs += [(pathlib.Path(o).stem, pathlib.Path(n).stem) for o, n in figs]
# longest first, and via placeholders so that 16->17 and 17->18 do not chain
subs.sort(key=lambda p: -len(p[0]))


def apply(path, skip_lines=()):
    p = ROOT / path
    lines = p.read_text().split('\n')
    out = []
    for i, line in enumerate(lines, start=1):
        if i in skip_lines:
            out.append(line)
            continue
        for k, (o, _) in enumerate(subs):
            line = line.replace(o, f'\x00{k}\x00')
        for k, (_, n) in enumerate(subs):
            line = line.replace(f'\x00{k}\x00', n)
        out.append(line)
    new = '\n'.join(out)
    if new != p.read_text():
        p.write_text(new)
        return True
    return False


touched = []
for path in ['_quarto.yml', 'code/README.md', 'code/sectors/bloch/README.md']:
    if apply(path):
        touched.append(path)
if apply('briefs/README.md', skip_lines=(3,)):
    touched.append('briefs/README.md')
for _, new in chapters:
    if apply(new):
        touched.append(new)
for _, new in figs:
    if new.endswith('.py') and apply(new):
        touched.append(new)
# the gap-map script reads its chapter by path, in two string pieces
gm = ROOT / 'code/figures/ch45_f33_gap_map.py'
gm.write_text(gm.read_text().replace('"part9-foundations" / "44-twelve-gaps.qmd"',
                                     '"part9-foundations" / "45-twelve-gaps.qmd"'))
touched.append('code/figures/ch45_f33_gap_map.py')
touched = sorted(set(touched))
git('add', '--', *touched, 'briefs/superalgebra_toda_renames.tsv', 'briefs/superalgebra_toda_renumber.py')
print(f'{len(renames)} renames; path fixes staged in {len(touched)} files:')
print('\n'.join(touched))
