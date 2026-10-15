#!/usr/bin/env python3
"""Check the cross-references and citations of the book without rendering it.

Run from anywhere: python3 check-refs.py [project-root]  (default: this file's directory, the book)

The project root may also be a Quarto project nested in the repository, such as
lecture-notes/defects/; the book's own check skips such nested projects (any
subdirectory with its own _quarto.yml), and their files are checked by running
this script on them. The bibliography is the one named by the project's
_quarto.yml (default references.bib).

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
  6. the list of original results in appendices/E-original-results.qmd (the
     ::: {.original-results} div) matches the chapters. Each item opens with a status
     label and ends with the references of the statements it summarises,
         - **[Theorem]** Both spin-two charges survive ... (@prp-spin-two).
     Every thm-/prp-/lem-/cor-/cnj- id in that final parenthesis must carry a
     label with the same first word and the same colour (status-tags.lua's
     kind) as the item's label. So the item may abbreviate a long label, but
     an upgrade, a downgrade or a refutation in the chapter fails until the
     list is updated. Section ids are allowed; their labels are not checked.
  7. every script path cited as inline code resolves. A code span counts as a
     script citation when its first whitespace-delimited token ends in ".py" or
     starts with textbook_code/, foundations_code/ or soliton_rmatrix_code/.
     The token must name an existing file, or an existing directory if it ends
     in "/", relative to the repository root or to textbook_code/ (the rule of
     script-links.lua, which turns such spans into GitHub links).
  8. (nested projects only) every link [text](book:<id>) names an id that the
     book defines; lecture-notes/book-links.lua turns such links into links to
     the published book.

Book files are all .qmd files outside briefs/, referee_scratch/, build or
tool directories and nested Quarto projects. Fenced code blocks and inline code are ignored, except for
the cell labels that they define.
"""
import collections
import re
import sys
from pathlib import Path

ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).resolve().parent).resolve()


def find_repo(start):
    """The repository root: the nearest directory at or above start that holds textbook_code/."""
    for d in (start, *start.parents):
        if (d / 'textbook_code').is_dir():
            return d
    return start


REPO = find_repo(ROOT)
SKIP = {'briefs', 'referee_scratch', '_site', '_book', '.quarto', '.claude', '.git', 'node_modules'}
PREFIXES = ('sec', 'eq', 'prp', 'thm', 'lem', 'cor', 'cnj', 'def', 'exm', 'exr', 'sol', 'rem',
            'alg', 'fig', 'tbl', 'lst', 'nte', 'tip', 'wrn', 'imp', 'cau', 'hyp')
STATUS_WORDS = ('Theorem', 'Established', 'Proposition', 'Numerical', 'Conjecture', 'Part')

P = '|'.join(PREFIXES)
ID_DEF = re.compile(r'\{#((?:%s)-[A-Za-z0-9_\-]+)' % P)
CELL_LABEL = re.compile(r'^\s*(?:#|%%%%)\|\s*label:\s*((?:%s)-[A-Za-z0-9_\-]+)' % P, re.M)
XREF = re.compile(r'(?<![\w.@/])-?@((?:%s)-[A-Za-z0-9_\-]*[A-Za-z0-9_])' % P)
CITE = re.compile(r'(?<![\w.@/])-?@([A-Za-z][A-Za-z0-9_:\-]*[A-Za-z0-9_])')
STATUS = re.compile(r'\*\*\[((?:%s)\b[^\]\n]*)\]\*\*' % '|'.join(STATUS_WORDS))
FENCE = re.compile(r'^(```|~~~).*?^\1', re.S | re.M)
INLINE_CODE = re.compile(r'`[^`\n]*`')
CODE_SPAN = re.compile(r'(?<!`)`([^`\n]+)`(?!`)')
SCRIPT_DIRS = ('textbook_code/', 'foundations_code/', 'soliton_rmatrix_code/')
ORIGINAL_RESULTS = 'appendices/E-original-results.qmd'


def project_files(root=ROOT):
    """The .qmd files of the project at root, without those of nested projects."""
    for f in sorted(root.rglob('*.qmd')):
        rel = f.relative_to(root)
        if SKIP.intersection(rel.parts[:-1]):
            continue
        if any((root.joinpath(*rel.parts[:i]) / '_quarto.yml').exists() for i in range(1, len(rel.parts))):
            continue
        yield rel


book_files = project_files


def bibliography():
    cfg = ROOT / '_quarto.yml'
    m = re.search(r'^bibliography:\s*(\S+)\s*$', cfg.read_text(), re.M) if cfg.exists() else None
    return (ROOT / (m.group(1) if m else 'references.bib')).resolve()


def prose(text):
    return INLINE_CODE.sub('', FENCE.sub('', text))


def lineno(text, pos):
    return text.count('\n', 0, pos) + 1


def status_kind(label):
    """The colour status-tags.lua gives a label; keep the two in step."""
    if label.startswith('Part'):
        return 'mixed'
    if 'sketch' in label or label.startswith('Numerical'):
        return 'partial'
    if label.startswith(('Theorem', 'Established')):
        return 'proved'
    if label.startswith('Conjecture'):
        return 'conjecture'
    return 'mixed'


def statement_label(t, sid):
    """The first status label inside the statement div {#sid}, or None."""
    m = re.search(r'^:::+\s*\{#%s[\s}]' % re.escape(sid), t, re.M)
    if not m:
        return None
    for line in t[m.end():].split('\n')[:15]:
        if re.match(r'^:::+\s*\{#', line):
            break
        s = STATUS.search(line)
        if s:
            return s.group(1).strip()
    return None


def trailing_refs(item):
    """The last top-level parenthesis of a list item, if the item ends with it."""
    body = item.rstrip().rstrip('.').rstrip()
    if not body.endswith(')'):
        return ''
    depth = 0
    for i in range(len(body) - 1, -1, -1):
        depth += {')': 1, '(': -1}.get(body[i], 0)
        if depth == 0:
            return body[i:]
    return ''


def original_results(text, defs):
    """Check 6. Returns (problems, note)."""
    t = text.get(Path(ORIGINAL_RESULTS), '')
    m = re.search(r'^:::+\s*\{\.original-results\}\s*$(.*?)^:::+\s*$', t, re.S | re.M)
    if not m:
        return [], f'no ::: {{.original-results}} list in {ORIGINAL_RESULTS}'
    start = lineno(t, m.start(1))
    items, cur = [], None
    for k, line in enumerate(m.group(1).split('\n')):
        if line.startswith('- '):
            cur = [start + k, line[2:]]
            items.append(cur)
        elif cur and line.startswith('  ') and line.strip():
            cur[1] += ' ' + line.strip()
        else:
            cur = None
    bad = []
    for n, item in items:
        where = f'{ORIGINAL_RESULTS}:{n}'
        lab = re.match(r'\*\*\[((?:%s)\b[^\]\n]*)\]\*\*' % '|'.join(STATUS_WORDS), item)
        if not lab:
            bad.append(f'{where} item does not open with a status label')
            continue
        lab = lab.group(1).strip()
        refs = XREF.findall(trailing_refs(item))
        if not refs:
            bad.append(f'{where} item does not end with its references, e.g. (@prp-x)')
            continue
        for r in refs:
            if r.split('-')[0] not in ('thm', 'prp', 'lem', 'cor', 'cnj') or r not in defs:
                continue
            f = Path(defs[r][0].rsplit(':', 1)[0])
            target = statement_label(text[f], r)
            if target is None:
                bad.append(f'{where} @{r} has no status label in {f}')
            elif (target.split()[0].rstrip(',;') != lab.split()[0].rstrip(',;')
                  or status_kind(target) != status_kind(lab)):
                bad.append(f'{where} [{lab[:40]}] but @{r} is [{target[:60]}]')
    return bad, f'{len(items)} items'


def script_paths(raw):
    """Check 7. Returns the code spans that cite a script that does not exist."""
    bad = []
    for f, t in raw.items():
        t = FENCE.sub(lambda m: re.sub(r'[^\n]', '', m.group(0)), t)   # keep line numbers
        for m in CODE_SPAN.finditer(t):
            words = m.group(1).split()
            if not words:
                continue
            token = words[0]
            if not (token.endswith('.py') or token.startswith(SCRIPT_DIRS)):
                continue
            want_dir = token.endswith('/')
            rel = token.rstrip('/')
            ok = bool(rel) and not token.startswith('/') and '..' not in token
            ok = ok and any((p.is_dir() if want_dir else p.is_file())
                            for p in (ROOT / rel, REPO / rel, REPO / 'textbook_code' / rel))
            if not ok:
                bad.append(f'{f}:{lineno(t, m.start())} `{m.group(1)[:70]}`')
    return bad


def book_links(text):
    """Check 8. Returns the [text](book:<id>) links whose id the book does not define."""
    book = set()
    for f in project_files(REPO):
        if f.parts and f.parts[0] == 'lecture-notes':
            continue
        book.update(ID_DEF.findall(prose((REPO / f).read_text())))
    bad = []
    for f, t in text.items():
        for m in re.finditer(r'\]\(book:([^)\s]*)\)', t):
            if m.group(1) not in book:
                bad.append(f'{f}:{lineno(t, m.start())} book:{m.group(1)}')
    return bad


def main():
    files = list(book_files())
    raw = {f: (ROOT / f).read_text() for f in files}
    text = {f: prose(t) for f, t in raw.items()}
    bib = bibliography().read_text()
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
    if (ROOT / ORIGINAL_RESULTS).exists():
        orig, note = original_results(text, defs)
        if note.startswith('no '):
            results.append(('original results match their statements', [note], 'WARN'))
        else:
            results.append((f'original results match their statements ({note})', orig, 'FAIL'))
    results.append(('script paths in code spans resolve', script_paths(raw), 'FAIL'))
    if REPO != ROOT:
        results.append(('links into the book resolve', book_links(text), 'FAIL'))
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
