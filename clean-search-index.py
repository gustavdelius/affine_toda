#!/usr/bin/env python3
"""Strip hidden per-chapter bibliography blocks from the search index.

Run as a project post-render step, after Quarto has generated _site and search.json.
Quarto inserts a hidden <div id="refs"> at the end of each chapter HTML file to
power citation hover tooltips. Because this div falls inside the chapter's final
section (typically 'Notes and references'), Quarto's search indexer extracts its
text and indexes all bibliography entries as part of that section.

This script strips that appended bibliography text from every chapter entry in
_site/search.json, keeping only the dedicated references.html page as the single
search target for bibliography entries.
"""

import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parent
SITE = ROOT / "_site"
SEARCH_JSON = SITE / "search.json"
REFERENCES_PAGE = "references.html"


def main() -> None:
    site_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else SITE
    search_json = site_dir if site_dir.name == "search.json" else site_dir / "search.json"
    if not search_json.exists():
        return

    try:
        data = json.loads(search_json.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as e:
        print(f"clean-search-index: warning: could not read {search_json}: {e}", file=sys.stderr)
        return

    ref_entry = next((d for d in data if d.get("objectID") == REFERENCES_PAGE), None)
    ref_text = ref_entry.get("text", "").strip() if ref_entry else ""
    bib_first_line = ref_text.splitlines()[0] if ref_text else None

    modified = 0
    for item in data:
        if item.get("objectID", "").startswith(REFERENCES_PAGE):
            continue
        text = item.get("text", "")
        idx = text.find(bib_first_line) if bib_first_line else -1
        if idx == -1:
            m = re.search(r"\n{2,}\s*\[1\]\s+[A-Z]", text)
            if m:
                idx = m.start()
        if idx != -1:
            item["text"] = text[:idx].rstrip()
            modified += 1

    if modified:
        search_json.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
