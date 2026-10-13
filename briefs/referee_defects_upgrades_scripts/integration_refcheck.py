"""Reference check over all book .qmd files (from _quarto.yml chapters, else every .qmd outside referee/brief dirs):
every @sec-/@eq-/@prp-/@fig-/@exr-/@thm-/@cnj-/@tbl-/@lem-/@cor-/@def-/@rem- reference resolves to a defined id,
no id is defined twice, and every citation key [@key ...] / -@key / @key (non-crossref) is in references.bib."""
import re, os, glob, sys
root = '/home/gustav/Git/affine_toda/.claude/worktrees/integrable-defects-framework-a54a84'
files = sorted(f for f in glob.glob(os.path.join(root, '**', '*.qmd'), recursive=True)
               if '/referee_scratch/' not in f and '/_book/' not in f and '/.quarto/' not in f and '/briefs/' not in f)
prefixes = ('sec', 'eq', 'prp', 'fig', 'exr', 'thm', 'cnj', 'tbl', 'lem', 'cor', 'def', 'rem', 'exm', 'nte', 'lst')
defs = {}
refs = []
cites = []
bib = open(os.path.join(root, 'references.bib')).read()
keys = set(re.findall(r'@\w+\s*\{\s*([^,\s]+)\s*,', bib))
for f in files:
    txt = open(f).read()
    # strip code blocks
    txt_nc = re.sub(r'```.*?```', '', txt, flags=re.S)
    for m in re.finditer(r'#((?:%s)-[A-Za-z0-9_\-]+)' % '|'.join(prefixes), txt_nc):
        # definitions appear as {#id ...}
        pass
    for m in re.finditer(r'\{#((?:%s)-[A-Za-z0-9_\-]+)' % '|'.join(prefixes), txt_nc):
        defs.setdefault(m.group(1), []).append(os.path.relpath(f, root))
    for m in re.finditer(r'#\|\s*label:\s*((?:%s)-[A-Za-z0-9_\-]+)' % '|'.join(prefixes), txt):
        defs.setdefault(m.group(1), []).append(os.path.relpath(f, root))
    for m in re.finditer(r'(?<![\w.])-?@((?:%s)-[A-Za-z0-9_\-]*[A-Za-z0-9])' % '|'.join(prefixes), txt_nc):
        refs.append((m.group(1), os.path.relpath(f, root)))
    for m in re.finditer(r'(?<![\w.@])-?@([A-Za-z][A-Za-z0-9_:\-]*[A-Za-z0-9])', txt_nc):
        k = m.group(1)
        if k.split('-')[0] in prefixes:
            continue
        cites.append((k, os.path.relpath(f, root)))
dups = {k: v for k, v in defs.items() if len(v) > 1}
unres = sorted({(k, f) for k, f in refs if k not in defs})
unknown = sorted({(k, f) for k, f in cites if k not in keys})
print(f'{len(files)} files, {len(defs)} ids, {len(refs)} cross-references, {len(cites)} citations')
print('duplicates:', dups)
print('unresolved:', unres[:20])
print('unknown citation keys:', unknown[:20])
print('PASS' if not dups and not unres and not unknown else 'FAIL')
