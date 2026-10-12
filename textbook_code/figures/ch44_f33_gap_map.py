"""Figure F33 (chapter 36, sec-summary-results, exr-gap-map): which chapter addresses which of the twelve gaps H1-H12,
and with what status.  The figure is generated from the text: it parses the table of sec-summary-results in
part9-foundations/44-twelve-gaps.qmd, resolves every cross-reference (@thm-..., @sec-..., ...) in a row to the
chapter file that defines it, and colours the row by the first words of its status column.

Checks: the table has the twelve rows H1..H12 in order; every cross-reference of every row resolves to exactly one
definition in a numbered chapter; every status is classified; each gap is addressed by at least one chapter.
Writes part9-foundations/figures/ch44_f33_gap_map.svg."""
import re
from pathlib import Path

import numpy as np
import style
from style import plt, report

style.setup()
ROOT = style.ROOT
text = (ROOT / "part9-foundations" / "44-twelve-gaps.qmd").read_text()
table = text.split("## Summary of results")[1].split("## Exercises")[0]
rows = [r for r in table.splitlines() if re.match(r"\|\s*H\d+\s*\|", r)]
report("the summary table has the twelve rows H1..H12 in order",
       [re.match(r"\|\s*(H\d+)", r).group(1) for r in rows] == [f"H{i}" for i in range(1, 13)])

# label -> chapter number (from the file name 'NN-....qmd')
chapter_of = {}
for f in sorted(ROOT.glob("part*/[0-9][0-9]-*.qmd")):
    n = int(f.name[:2])
    for m in re.finditer(r"\{#((?:sec|thm|prp|lem|cnj|exr|tbl|eq|fig)-[\w-]+)", f.read_text()):
        chapter_of.setdefault(m.group(1), set()).add(n)

short = {1: "existence of the regularized theory", 2: "sign problem", 3: "continuum limit", 4: "Hilbert space",
         5: "reality of the vacuum energy", 6: "classical reality", 7: "which complex saddles", 8: "semiclassical quantization",
         9: "non-unitarity", 10: "exact S-matrices", 11: "missing classical solutions", 12: "unitary subtheories"}
cells, status = {}, {}
unresolved = []
ambiguous = []
for r in rows:
    cols = [c.strip() for c in r.strip().strip("|").split("|")]
    h = int(cols[0][1:])
    st = cols[2].lower()
    if st.startswith("closed"):
        status[h] = "closed"
    elif st.startswith(("proved here", "reduced", "reclassified", "partly open", "established")):
        status[h] = "partly"
    elif st.startswith("open"):
        status[h] = "open"
    labels = re.findall(r"@((?:sec|thm|prp|lem|cnj|exr|tbl|eq|fig)-[\w-]+)", cols[2])
    chs = set()
    for lab in labels:
        if lab not in chapter_of:
            unresolved.append(lab)
        else:
            if len(chapter_of[lab]) != 1:
                ambiguous.append(lab)
            chs |= chapter_of[lab]
    cells[h] = chs
report("every cross-reference in the table resolves to a definition in a numbered chapter", not unresolved)
report("each cross-reference is defined in exactly one chapter", not ambiguous)
report("every status is classified as closed, partly or open", sorted(status) == list(range(1, 13)))
report("each gap is addressed by at least one chapter", all(cells[h] for h in range(1, 13)))
report("statuses: H11 is the only open gap", [h for h, s in status.items() if s == "open"] == [11])

cols_ch = sorted(set().union(*cells.values()))
col_color = {"closed": style.GREEN, "partly": style.ORANGE, "open": style.RED}
fig, ax = plt.subplots(figsize=(7.0, 3.9))
for i, h in enumerate(range(1, 13)):
    y = 11 - i
    ax.text(-0.75, y, f"H{h}", ha="right", va="center", fontsize=8.5, color=col_color[status[h]], fontweight="bold")
    ax.text(-0.55, y, short[h], ha="left", va="center", fontsize=7.5)
    for j, c in enumerate(cols_ch):
        x = 5.2 + j
        if c in cells[h]:
            ax.add_patch(plt.Rectangle((x - 0.4, y - 0.38), 0.8, 0.76, color=col_color[status[h]], alpha=0.85, lw=0))
        else:
            ax.add_patch(plt.Rectangle((x - 0.4, y - 0.38), 0.8, 0.76, color=style.LIGHTGREY, alpha=0.35, lw=0))
for j, c in enumerate(cols_ch):
    ax.text(5.2 + j, 12.05, str(c), ha="center", va="bottom", fontsize=7.5)
ax.text(5.2 + (len(cols_ch) - 1) / 2, 12.9, "chapter", ha="center", va="bottom", fontsize=8)
for lab, c, x in (("closed", style.GREEN, 0.0), ("partly closed", style.ORANGE, 2.0), ("open", style.RED, 4.6)):
    ax.add_patch(plt.Rectangle((x - 0.9, -1.45), 0.5, 0.5, color=c, alpha=0.85, lw=0))
    ax.text(x - 0.3, -1.2, lab, va="center", fontsize=7.5)
ax.set_xlim(-2.4, 5.2 + len(cols_ch) - 0.5)
ax.set_ylim(-1.8, 13.7)
ax.axis("off")
fig.tight_layout()
style.save(fig, "ch44_f33_gap_map", style.OUT_PART9)
print("chapters:", cols_ch, "status:", status)
raise SystemExit(0 if style.all_ok() else 1)
