#!/usr/bin/env python3
"""Check the cross-references and citations of the book without rendering it.

Run from anywhere: python3 check-refs.py [repo-root]  (default: this file's directory)

Checks, one PASS, FAIL or WARN line each; exit status 1 on any FAIL:
  1. every id ({#sec-...}, or a code cell's "#| label:" / "%%| label:") is defined once;
  2. every cross-reference @sec-/@eq-/@prp-/... resolves to a defined id;
  3. every citation key @key is in references.bib;
  4. (WARN only) cross-references on unnumbered pages (part pages, and pages
     marked .unnumbered): Quarto prints the target's title there instead of its
     number, so the sentence must read correctly with a title;
  5. no status label such as **[Numerical, @sec-x]** is read by Pandoc as a
     citation. Pandoc does so when every ";"-separated part of the bracket holds
     an @key; the brackets are then lost and status-tags.lua cannot tag the
     label. Write **[Numerical]** (@sec-x) instead.

Book files are all .qmd files outside briefs/, referee_scratch/ and build or
tool directories. Fenced code blocks and inline code are ignored, except for
the cell labels that they define.
"""
import collections
import re
import sys
from pathlib import Path

ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).resolve().parent)
SKIP = {'briefs', 'referee_scratch', '_site', '_book', '.quarto', '.claude', '.git', 'node_modules'}
PREFIXES = ('sec', 'eq', 'prp', 'thm', 'lem', 'cor', 'cnj', 'def', 'exm', 'exr', 'sol', 'rem',
            'alg', 'fig', 'tbl', 'lst', 'nte', 'tip', 'wrn', 'imp', 'cau')
STATUS_WORDS = ('Theorem', 'Established', 'Proposition', 'Numerical', 'Conjecture', 'Part')

P = '|'.join(PREFIXES)
ID_DEF = re.compile(r'\{#((?:%s)-[A-Za-z0-9_\-]+)' % P)
CELL_LABEL = re.compile(r'^\s*(?:#|%%%%)\|\s*label:\s*((?:%s)-[A-Za-z0-9_\-]+)' % P, re.M)
XREF = re.compile(r'(?<![\w.@/])-?@((?:%s)-[A-Za-z0-9_\-]*[A-Za-z0-9_])' % P)
CITE = re.compile(r'(?<![\w.@/])-?@([A-Za-z][A-Za-z0-9_:\-]*[A-Za-z0-9_])')
STATUS = re.compile(r'\*\*\[((?:%s)\b[^\]\n]*)\]\*\*' % '|'.join(STATUS_WORDS))
FENCE = re.compile(r'^(```|~~~).*?^\1', re.S | re.M)
INLINE_CODE = re.compile(r'`[^`\n]*`')


def book_files():
    for f in sorted(ROOT.rglob('*.qmd')):
        rel = f.relative_to(ROOT)
        if not SKIP.intersection(rel.parts[:-1]):
            yield rel


def prose(text):
    return INLINE_CODE.sub('', FENCE.sub('', text))


def lineno(text, pos):
    return text.count('\n', 0, pos) + 1


def main():
    files = list(book_files())
    raw = {f: (ROOT / f).read_text() for f in files}
    text = {f: prose(t) for f, t in raw.items()}
    bib = (ROOT / 'references.bib').read_text()
    keys = set(re.findall(r'^@\w+\s*\{\s*([^,\s]+)\s*,', bib, re.M))

    defs = collections.defaultdict(list)
    for f in files:
        for m in ID_DEF.finditer(text[f]):
            defs[m.group(1)].append(f'{f}:{lineno(text[f], m.start())}')
        for m in CELL_LABEL.finditer(raw[f]):
            defs[m.group(1)].append(f'{f}:{lineno(raw[f], m.start())}')

    dups, unres, unknown, unnumbered, status = [], [], [], [], []
    for k, where in sorted(defs.items()):
        if len(where) > 1:
            dups.append(f'{k} defined at {", ".join(where)}')
    for f in files:
        t = text[f]
        first_heading = re.search(r'^#\s.*$', t, re.M)
        is_unnumbered = f.name == 'part.qmd' or bool(
            first_heading and re.search(r'\{[^}]*\.unnumbered', first_heading.group(0)))
        for m in XREF.finditer(t):
            where = f'{f}:{lineno(t, m.start())}'
            if m.group(1) not in defs:
                unres.append(f'{where} @{m.group(1)}')
            if is_unnumbered:
                unnumbered.append(f'{where} @{m.group(1)}')
        for m in CITE.finditer(t):
            k = m.group(1)
            if k.split('-')[0] not in PREFIXES and k not in keys:
                unknown.append(f'{f}:{lineno(t, m.start())} @{k}')
        for m in STATUS.finditer(t):
            parts = m.group(1).split(';')
            if all(re.search(r'@[A-Za-z]', s) for s in parts):
                status.append(f'{f}:{lineno(t, m.start())} {m.group(0)[:80]}')

    results = [
        ('ids defined once', dups, 'FAIL'),
        ('cross-references resolve', unres, 'FAIL'),
        ('citation keys in references.bib', unknown, 'FAIL'),
        ('cross-references on unnumbered pages print titles', unnumbered, 'WARN'),
        ('status labels not parsed as citations', status, 'FAIL'),
    ]
    print(f'{len(files)} files, {len(defs)} ids, {len(keys)} bibliography keys')
    ok = True
    for name, bad, severity in results:
        print(('PASS' if not bad else severity) + f'  {name}' + (f' ({len(bad)})' if bad else ''))
        for b in bad[:30]:
            print('      ' + b)
        if len(bad) > 30:
            print(f'      ... and {len(bad) - 30} more')
        ok = ok and not (bad and severity == 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
