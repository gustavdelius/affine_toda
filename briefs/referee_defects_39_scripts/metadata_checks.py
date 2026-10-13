"""Quarto reference audit with sentence punctuation excluded."""
import hashlib
import re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]  # repository root

chapter=ROOT/'part8-defects/39-defect-solitons.qmd'
text=chapter.read_text()
files=[ROOT/'index.qmd']+list(ROOT.glob('part*/*.qmd'))
ids={v for file in files for v in re.findall(r'\{#([^\s}]+)',file.read_text())}
keys=set(re.findall(r'@\w+\{([^,]+),',(ROOT/'references.bib').read_text()))
refs={v.rstrip('.') for v in re.findall(r'@([A-Za-z][\w:.-]*)',text)}
print(f'distinct references/citations={len(refs)}; unresolved={sorted(refs-ids-keys)}')
print('citation keys='+', '.join(sorted(refs&keys)))
print('chapter sha256='+hashlib.sha256(chapter.read_bytes()).hexdigest())
print('The naive scanner in supplemental_checks.py includes sentence-final periods; this scanner removes them.')
assert not refs-ids-keys
