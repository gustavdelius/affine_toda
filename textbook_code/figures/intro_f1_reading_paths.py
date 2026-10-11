"""Figure F1b (introduction, sec-reader-guide): the five reading paths.

One row per path; boxes are chapter ranges, coloured by Part. The chapter lists are checked
against the "Reading paths" bullets of 01-introduction.qmd, so that the figure cannot drift
from the text.

Writes figures/intro_f1_reading_paths.svg.
"""
import re

import numpy as np
from matplotlib.patches import FancyBboxPatch

import style
from style import ROOT, plt, report

style.setup()

PATHS = [
    ("Integrable field theory", [(2, 3), (5, 6), (9, 10), (12, 14), (20, 22)]),
    ("Quantum groups in field theory", [(2, 3), (5, 8), (20, 22), (27, 27)]),
    ("Integrable boundaries", [(2, 3), (7, 8), (9, 10), (13, 13), "VII"]),
    ("Non-unitary field theory", [(2, 4), (9, 11), "VIII"]),
    ("For specialists", ["V–VIII"]),
]
PART_OF = lambda c: ("I" if c <= 4 else "II" if c <= 8 else "III" if c <= 11 else "IV" if c <= 15
                     else "V" if c <= 19 else "VI" if c <= 27 else "VII" if c <= 35 else "VIII")
COLORS = {"I": style.BLUE, "II": style.GREEN, "III": style.ORANGE, "IV": style.RED,
          "V": style.PURPLE, "VI": style.TEAL, "VII": "#c2569b", "VIII": style.GREY}

# ---- check against the text of the introduction ---------------------------
text = (ROOT / "01-introduction.qmd").read_text()
bullets = re.findall(r"^- \*(?:A course on|Integrable boundaries|Non-unitary|For specialists)[^\n]*", text, flags=re.M)
report("the introduction has five reading-path bullets", len(bullets) == 5)


def chapters_in(bullet):
    seg = bullet.split("chapters", 1)[1].split(".", 1)[0] if "chapters" in bullet else ""
    seg = seg.split("Part VII")[0]
    out = set()
    for a, b in re.findall(r"(\d+)(?:–(\d+))?", seg):
        out.update(range(int(a), int(b or a) + 1))
    return out


for (name, boxes), bullet in zip(PATHS[:4], bullets[:4]):
    mine = set()
    for b in boxes:
        if isinstance(b, tuple):
            mine.update(range(b[0], b[1] + 1))
    report(f"chapters of '{name}' agree with the text", mine == chapters_in(bullet))
report("'For specialists' is Parts V-VIII in the text", "Parts V–VIII" in bullets[4])

# ---- figure ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.0, 3.3))
ax.set_xlim(0, 11.6)
ax.set_ylim(-0.2, 5.3)
ax.axis("off")
W, H, GAP = 1.05, 0.58, 0.38
x0 = 3.45
for r, (name, boxes) in enumerate(PATHS):
    y = 4.55 - r
    ax.text(0.0, y, name, ha="left", va="center", fontsize=8.5)
    x = x0
    for i, b in enumerate(boxes):
        if isinstance(b, tuple):
            lab = f"ch. {b[0]}–{b[1]}" if b[0] != b[1] else f"ch. {b[0]}"
            part = PART_OF(b[0])
        else:
            lab = f"Part {b}" if b in ("VII", "VIII") else f"Parts {b}"
            part = b if b in ("VII", "VIII") else "V"
        w = W if len(lab) < 11 else 1.35
        ax.add_patch(FancyBboxPatch((x, y - H / 2), w, H, boxstyle="round,pad=0.02,rounding_size=0.08",
                                    fc=COLORS[part] + "33", ec=COLORS[part], lw=1.4))
        ax.text(x + w / 2, y, lab, ha="center", va="center", fontsize=8)
        if i < len(boxes) - 1:
            ax.annotate("", xy=(x + w + GAP, y), xytext=(x + w + 0.03, y),
                        arrowprops=dict(arrowstyle="->", color=style.GREY, lw=1.1))
        x += w + GAP
# legend of parts
lx = 0.0
for part in ("I", "II", "III", "IV", "V", "VI", "VII", "VIII"):
    ax.add_patch(FancyBboxPatch((lx, -0.12), 0.28, 0.2, boxstyle="round,pad=0.01,rounding_size=0.04",
                                fc=COLORS[part] + "33", ec=COLORS[part], lw=1.2))
    ax.text(lx + 0.38, -0.02, part, ha="left", va="center", fontsize=7.5)
    lx += 0.38 + 0.2 + 0.1 * len(part) + 0.35
ax.text(lx + 0.0, -0.02, "Part (colour)", ha="left", va="center", fontsize=7.5)
fig.tight_layout()
style.save(fig, "intro_f1_reading_paths", style.OUT_INTRO)
raise SystemExit(0 if style.all_ok() else 1)
