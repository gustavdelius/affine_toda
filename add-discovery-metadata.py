#!/usr/bin/env python3
"""Add canonical links and editor-only citation metadata to the HTML book.

Quarto uses ``book.site-url`` when it creates the sitemap and social-card
metadata, but it does not currently emit canonical links or Highwire citation
tags for an editor-only book. This post-render step supplies them without
misrepresenting the editor as an author.

Running the script more than once is harmless: its previous block is replaced.
"""

from html import escape
from pathlib import Path
import re
from urllib.parse import quote


ROOT = Path(__file__).resolve().parent
SITE = ROOT / "_site"
SITE_URL = "https://gustavdelius.github.io/affine_toda/"
START = "  <!-- discovery-metadata:start -->"
END = "  <!-- discovery-metadata:end -->"
OLD_BLOCK = re.compile(rf"\n?{re.escape(START)}.*?{re.escape(END)}\n?", re.DOTALL)


def canonical_url(path: Path) -> str:
    relative = path.relative_to(SITE).as_posix()
    if relative == "index.html":
        return SITE_URL
    return SITE_URL + quote(relative)


def metadata_block(path: Path) -> str:
    lines = [
        START,
        f'  <link rel="canonical" href="{escape(canonical_url(path), quote=True)}">',
    ]
    if path == SITE / "index.html":
        lines.extend(
            [
                '  <meta name="citation_title" content="Affine Toda field theory: A unified review, with new results">',
                '  <meta name="citation_editor" content="Delius, Gustav">',
                '  <meta name="citation_publication_date" content="2026">',
                '  <meta name="citation_language" content="en">',
                f'  <meta name="citation_public_url" content="{SITE_URL}">',
                f'  <meta name="citation_pdf_url" content="{SITE_URL}affine-toda.pdf">',
                '  <meta name="citation_keywords" content="affine Toda field theory; integrable quantum field theory; solitons; quantum groups; exact S-matrices; integrable boundaries; integrable defects; non-unitary">',
            ]
        )
    lines.append(END)
    return "\n".join(lines)


def main() -> None:
    if not SITE.is_dir():
        return
    for path in SITE.rglob("*.html"):
        text = path.read_text(encoding="utf-8")
        text = OLD_BLOCK.sub("\n", text)
        if "</head>" not in text:
            continue
        text = text.replace("</head>", metadata_block(path) + "\n</head>", 1)
        path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
