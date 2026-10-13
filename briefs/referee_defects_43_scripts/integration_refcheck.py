"""Reference check over all book .qmd files: ids defined once, every @sec/eq/prp/fig/exr/thm/cnj/tbl/... reference
resolves, every citation key is in references.bib."""
import re, os, glob, collections
root = '.'
files = [f for f in glob.glob('**/*.qmd', recursive=True)
         if not f.startswith(('referee_scratch', 'briefs', '_book', '.quarto', '.claude', 'node_modules'))]
prefixes = ('sec', 'eq', 'prp', 'fig', 'exr', 'thm', 'cnj', 'tbl', 'lem', 'cor', 'def', 'exm', 'rem', 'alg')
idpat = re.compile(r'\{#((?:' + '|'.join(prefixes) + r')-[A-Za-z0-9_\-]+)')
refpat = re.compile(r'(?<![\w@])-?@((?:' + '|'.join(prefixes) + r')-[A-Za-z0-9_\-]*[A-Za-z0-9_])')
citepat = re.compile(r'(?<![\w@./])-?@([A-Za-z][A-Za-z0-9_:\-]*[A-Za-z0-9_])')
bib = open('references.bib').read()
keys = set(re.findall(r'@\w+\{([^,\s]+),', bib))
defs = collections.Counter(); where = {}
for f in files:
    txt = open(f).read()
    # strip code blocks
    txt2 = re.sub(r'```.*?```', '', txt, flags=re.S)
    for m in idpat.finditer(txt2):
        defs[m.group(1)] += 1; where.setdefault(m.group(1), []).append(f)
dups = {k: where[k] for k, c in defs.items() if c > 1}
unres, badcite = [], []
for f in files:
    txt = re.sub(r'```.*?```', '', open(f).read(), flags=re.S)
    txt = re.sub(r'`[^`\n]*`', '', txt)                     # inline code
    for m in refpat.finditer(txt):
        if m.group(1) not in defs:
            unres.append((f, m.group(1)))
    for m in citepat.finditer(txt):
        k = m.group(1)
        if k.split('-')[0] in prefixes:
            continue
        if k not in keys:
            badcite.append((f, k))
print(f'{len(files)} files, {len(defs)} ids; duplicates: {dups}; unresolved: {unres}; unknown citation keys: {sorted(set(badcite))}')
