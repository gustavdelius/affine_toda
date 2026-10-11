#!/usr/bin/env python3
"""Add citation backlinks to the HTML references page.

Run as the project post-render step, after Quarto has written _site.
Quarto's citeproc renders each chapter on its own, so the references page
has no record of where an entry is cited. This script works on the finished
HTML instead: it gives every citation span an id, then appends to each
bibliography entry links to those spans, in every chapter that cites it.

Only the HTML is changed. The PDF has no backlinks.
Running it again is harmless: old backlinks and ids are replaced.
"""
import os
import re

ROOT = os.path.dirname(os.path.abspath(__file__))
SITE = os.path.join(ROOT, "_site")
CONFIG = os.path.join(ROOT, "_quarto.yml")
REFERENCES = "references.html"

# A citation span, with or without the id this script gave it earlier.
CITATION = re.compile(r'<span(?: id="[^"]*")? class="citation" data-cites="([^"]*)">')
ENTRY = re.compile(r'<div id="ref-([^"]+)" class="csl-entry"')
DIV_TAG = re.compile(r"<div\b|</div>")


def chapter_pages():
    """HTML pages of the book, in the order _quarto.yml lists the .qmd files."""
    with open(CONFIG, encoding="utf-8") as f:
        for line in f:
            match = re.match(r"\s*-\s+([\w./-]+)\.qmd\s*$", line)
            if match:
                yield match.group(1) + ".html"


def entry_end(text, start):
    """Index just past the </div> that closes the div opened at start."""
    depth = 0
    for tag in DIV_TAG.finditer(text, start):
        depth += 1 if tag.group(0) != "</div>" else -1
        if depth == 0:
            return tag.end()
    raise ValueError("unclosed bibliography entry at offset %d" % start)


def backlinks_markup(places):
    links = ", ".join('<a href="%s#%s">%d</a>' % (page, anchor, number)
                      for number, (page, anchor) in enumerate(places, 1))
    return ' <span class="backlinks">cited: %s</span>' % links


def main():
    if not os.path.isdir(SITE):
        return
    cited_in = {}  # citation key -> [(page, anchor), ...] in book order

    for page in chapter_pages():
        if page == REFERENCES or not os.path.exists(os.path.join(SITE, page)):
            continue
        path = os.path.join(SITE, page)
        with open(path, encoding="utf-8") as f:
            text = f.read()
        prefix = "cite-" + page[:-len(".html")].replace("/", "-") + "-"
        count = 0

        def give_id(match):
            nonlocal count
            count += 1
            anchor = prefix + str(count)
            for key in match.group(1).split():
                cited_in.setdefault(key, []).append((page, anchor))
            return '<span id="%s" class="citation" data-cites="%s">' % (anchor, match.group(1))

        text = CITATION.sub(give_id, text)
        with open(path, "w", encoding="utf-8") as f:
            f.write(text)

    path = os.path.join(SITE, REFERENCES)
    if not os.path.exists(path):
        return
    with open(path, encoding="utf-8") as f:
        text = f.read()
    text = re.sub(r' <span class="backlinks">.*?</span>', "", text)

    # Work from the end so that earlier offsets stay valid.
    for match in reversed(list(ENTRY.finditer(text))):
        places = cited_in.get(match.group(1))
        if not places:
            continue
        start = match.start()
        end = entry_end(text, start)
        # Put the links at the end of the entry text, inside the text div.
        text_div = text.find('class="csl-right-inline"', start, end)
        if text_div != -1:
            insert_at = text.find("</div>", text_div, end)
        else:
            insert_at = end - len("</div>")
        text = text[:insert_at] + backlinks_markup(places) + text[insert_at:]

    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


if __name__ == "__main__":
    main()
