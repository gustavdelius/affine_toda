#!/usr/bin/env python3
"""Add citation backlinks to the HTML references page.

Run as the project post-render step, after Quarto has written _site.
Quarto's citeproc renders each chapter on its own, so the references page
has no record of where an entry is cited. This script works on the finished
HTML instead: it finds, for every citation span, the subsection or
subsubsection that contains it, and appends to each bibliography entry one
link per such section, labelled by the section number (e.g. "cited: 3.1,
12.4.2"). Hovering over a link shows the section's title. Several citations
in the same section give one link. A citation before a chapter's first
section links to the top of the chapter and is labelled by its number.

Only the references page is changed, and only its HTML; the PDF has no
backlinks. Running it again is harmless: old backlinks are replaced.
"""
import html
import os
import re

ROOT = os.path.dirname(os.path.abspath(__file__))
SITE = os.path.join(ROOT, "_site")
CONFIG = os.path.join(ROOT, "_quarto.yml")
REFERENCES = "references.html"
DEEPEST_LEVEL = 3  # link to subsubsections (###), not to deeper headings

# A citation span. Earlier versions of this script gave it an id.
CITATION = re.compile(r'<span(?: id="[^"]*")? class="citation" data-cites="([^"]*)"')
# Opening and closing tags of the <section> elements Quarto makes for headings.
SECTION_TAG = re.compile(r'<section\b([^>]*)>|</section>')
ENTRY = re.compile(r'<div id="ref-([^"]+)" class="csl-entry"')
DIV_TAG = re.compile(r"<div\b|</div>")
TAG = re.compile(r"<[^>]+>")
SPACE = re.compile(r"\s+")


def chapter_pages():
    """HTML pages of the book, in the order _quarto.yml lists the .qmd files."""
    with open(CONFIG, encoding="utf-8") as f:
        for line in f:
            match = re.match(r"\s*-\s+([\w./-]+)\.qmd\s*$", line)
            if match:
                yield match.group(1) + ".html"


def attribute(attributes, name):
    match = re.search(r'\b%s="([^"]*)"' % name, attributes)
    return match.group(1) if match else None


def plain_text(fragment):
    """The text of an HTML fragment, without tags and on one line."""
    return SPACE.sub(" ", html.unescape(TAG.sub("", fragment))).strip()


def heading_title(text, start):
    """Title of the heading that opens the section starting at start,
    without its number."""
    match = re.compile(r"<h(\d)\b[^>]*>(.*?)</h\1>", re.S).search(text, start)
    if not match:
        return ""
    inner = re.sub(r'<span class="header-section-number">.*?</span>', "",
                   match.group(2), flags=re.S)
    # A hover title is plain text, so MathJax cannot typeset it: drop the
    # delimiters of inline math and leave its TeX source.
    return re.sub(r"\\[()]", "", plain_text(inner))


def chapter_place(text):
    """(anchor, label, title) for the top of a chapter page."""
    number = re.search(r'<span class="chapter-number">(.*?)</span>', text, re.S)
    title = re.search(r'<span class="chapter-title">(.*?)</span>', text, re.S)
    if not title:
        title = re.search(r"<h1 class=\"title\"[^>]*>(.*?)</h1>", text, re.S)
    title = plain_text(title.group(1)) if title else ""
    label = plain_text(number.group(1)) if number else title
    return "", label, title


def citations_by_place(text):
    """Yield (citation keys, (anchor, label, title)) for every citation in the
    main text of a page, where the place is the innermost enclosing section of
    level DEEPEST_LEVEL or less, or the top of the page."""
    main_start = text.find("<main")
    main_end = text.find("</main>", main_start)
    if main_start == -1 or main_end == -1:
        return
    top = chapter_place(text)
    stack = []  # one entry per open <section>: a place, or None if too deep
    events = [(m.start(), "section", m) for m in SECTION_TAG.finditer(text, main_start, main_end)]
    events += [(m.start(), "cite", m) for m in CITATION.finditer(text, main_start, main_end)]
    for _, kind, match in sorted(events, key=lambda event: event[0]):
        if kind == "section":
            if match.group(0) == "</section>":
                if stack:
                    stack.pop()
                continue
            attributes = match.group(1)
            level = re.search(r'\blevel(\d)\b', attribute(attributes, "class") or "")
            anchor = attribute(attributes, "id")
            if not level or not anchor or int(level.group(1)) > DEEPEST_LEVEL:
                stack.append(None)
                continue
            title = heading_title(text, match.end())
            number = attribute(attributes, "data-number")
            stack.append((anchor, number or title, title))
        else:
            place = next((p for p in reversed(stack) if p is not None), top)
            yield match.group(1).split(), place


def backlinks_markup(places):
    links = []
    for page, (anchor, label, title) in places:
        href = page + ("#" + anchor if anchor else "")
        hover = "%s %s" % (label, title) if label != title else title
        links.append('<a href="%s" title="%s">%s</a>'
                     % (html.escape(href), html.escape(hover), html.escape(label)))
    return ' <span class="backlinks">cited: %s</span>' % ", ".join(links)


def main():
    if not os.path.isdir(SITE):
        return
    cited_in = {}  # citation key -> [(page, place), ...] in book order, no repeats

    for page in chapter_pages():
        if page == REFERENCES or not os.path.exists(os.path.join(SITE, page)):
            continue
        with open(os.path.join(SITE, page), encoding="utf-8") as f:
            text = f.read()
        for keys, place in citations_by_place(text):
            for key in keys:
                places = cited_in.setdefault(key, [])
                if (page, place) not in places:
                    places.append((page, place))

    path = os.path.join(SITE, REFERENCES)
    if not os.path.exists(path):
        return
    with open(path, encoding="utf-8") as f:
        text = f.read()
    text = re.sub(r' <span class="backlinks">.*?</span>(?=</div>)', "", text)

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


def entry_end(text, start):
    """Index just past the </div> that closes the div opened at start."""
    depth = 0
    for tag in DIV_TAG.finditer(text, start):
        depth += 1 if tag.group(0) != "</div>" else -1
        if depth == 0:
            return tag.end()
    raise ValueError("unclosed bibliography entry at offset %d" % start)


if __name__ == "__main__":
    main()
