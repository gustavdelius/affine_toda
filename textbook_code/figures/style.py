"""Shared style for the figures of the introduction and Parts I to VII.

Every figure script imports this module. Figures are designed for a canvas about 7 inches wide: the book
shows them at text width (about 5.4 inches in the PDF), so text of 8-9 pt stays readable. It fixes the matplotlib backend, one
palette that reads on both the light and the dark theme of the book (mid-tone
colours, grey text, transparent background), and the helpers `save` and `report`.

Figures are written as SVG into the book tree and committed, so that the HTML
build needs no Python:

  * `figures/`                    for the introduction (index.qmd),
  * `part1-background/figures/`   for chapters 2-4,
  * `part2-algebra/figures/`      for chapters 5-8,
  * `part3-classical/figures/`    for chapters 9-11,
  * `part4-real-coupling/figures/` for chapters 12-15,
  * `part5-semiclassics/figures/`  for chapters 16-19,
  * `part6-soliton-smatrices/figures/` for chapters 20-27,
  * `part7-half-line/figures/` for chapters 28-35.

Each script also prints one PASS or FAIL line per formula it plots, like the
other scripts in textbook_code/.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT_INTRO = ROOT / "figures"
OUT_PART1 = ROOT / "part1-background" / "figures"
OUT_PART2 = ROOT / "part2-algebra" / "figures"
OUT_PART3 = ROOT / "part3-classical" / "figures"
OUT_PART4 = ROOT / "part4-real-coupling" / "figures"
OUT_PART5 = ROOT / "part5-semiclassics" / "figures"
OUT_PART6 = ROOT / "part6-soliton-smatrices" / "figures"
OUT_PART7 = ROOT / "part7-half-line" / "figures"

GREY = "#868b94"          # text, axes: readable on white and on the dark theme
LIGHTGREY = "#868b9455"   # grid lines, light fills
BLUE = "#4c8bf5"
RED = "#e4572e"
GREEN = "#2fa86b"
ORANGE = "#f2a33a"
PURPLE = "#9b6bd6"
TEAL = "#2bb3b1"

_ok = True


def setup():
    plt.rcParams.update({
        "mathtext.fontset": "cm",
        "font.family": "serif",
        "font.size": 8.5,
        "svg.fonttype": "path",          # no font needed to view the SVG
        "axes.edgecolor": GREY,
        "axes.labelcolor": GREY,
        "axes.linewidth": 0.8,
        "xtick.color": GREY,
        "ytick.color": GREY,
        "text.color": GREY,
        "axes.facecolor": "none",
        "figure.facecolor": "none",
        "savefig.transparent": True,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 1.8,
        "legend.frameon": False,
    })


def save(fig, name, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"{name}.svg"
    fig.savefig(path, bbox_inches="tight", pad_inches=0.08, dpi=250)  # dpi only affects rasterized parts
    plt.close(fig)
    print(f"wrote {path.relative_to(ROOT)}")


def report(name, ok):
    global _ok
    _ok = _ok and bool(ok)
    print(("PASS " if ok else "FAIL ") + name)


def all_ok():
    return _ok
