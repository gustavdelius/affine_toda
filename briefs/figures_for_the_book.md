# Brief: adding figures to the book (lessons from the introduction and Part I)

You are adding diagrams to chapters that have none. The introduction and Part I (chapters 2–4) received nine figures on 2026-10-11 (merged to `main` as `4999958`). This brief records the decisions, the recipe and the pitfalls, so that you can do Parts II–VII without rediscovering them. Read it fully before drawing anything. The process described here worked; the pitfalls cost real time.

## What exists

| Item | Location |
|---|---|
| Shared style and helpers | `code/figures/style.py` |
| One script per figure, named `<chapter>_<fN>_<topic>.py` | `code/figures/` |
| Generated SVGs, committed | `figures/` (introduction) and `part1-background/figures/` (Part I) |
| One README row per script | `code/README.md`, section "Figures" |
| Embedded in | `01-introduction.qmd`, `part1-background/0[234]*.qmd` |

Figure ids used so far: `fig-dependencies`, `fig-reading-paths`, `fig-two-regimes`, `fig-theta-strip`, `fig-factorization`, `fig-bound-state-poles`, `fig-sg-solutions`, `fig-sg-coupling-axis`, `fig-dimer-krein`, `fig-krein-geometry`. Reuse the naming pattern: `fig-` plus a topic, never a number.

## Decisions already made (do not reopen them without a reason)

1. **Pre-generated SVG, committed.** The GitHub HTML build installs no Python, so figures must not be generated at render time. Put the script in `code/figures/`, run it, commit both the script and the SVG. For a new Part, create `partN-.../figures/`.
2. **matplotlib for anything quantitative or geometric.** Mermaid only for plain dependency graphs, and only when the layout is a simple left-to-right flow. The user said PDFs are not built on GitHub, so mermaid is allowed, but see the pitfalls: its layout is hard to control.
3. **Every script checks what it plots.** Each prints one `PASS` or `FAIL` line per formula, in the convention of the other `code/` scripts, and exits non-zero on failure. A figure that is not backed by a check can silently disagree with the text. Examples: the residue of $f_x$ equals $2i\tan\pi x$; a breather solves the sine-Gordon equation by finite differences; the reading-paths figure parses `01-introduction.qmd` and compares its chapter lists with the bullets.
4. **One transparent SVG for both themes.** Mid-tone colours and grey text (`style.GREY = #868b94`) read on the light and on the dark theme. No light/dark pairs.
5. **Design for a 7 inch canvas with 8–9 pt text.** See the first pitfall.

## The recipe

1. **Read the chapter.** Pick a figure only where a reader must otherwise assemble geometry in their head, or where a statement is easy to misread. Tie it to a section id and, where possible, an exercise (`exr-strip` got the strip map). Aim for one or two per chapter. Skip decoration.
2. **Check the current text first.** Another agent may have renamed symbols (the reading fixes renamed $\Theta_{s-1}\to\mathcal X_{s-1}$ and Fisher's $h,g\to\nu,\kappa$ while the figures were being planned). Use the symbols as printed now, and re-read `briefs/part<N>_reading_fixes.md` (section C, notation clashes) so that you do not bake a symbol into a figure that is about to be renamed.
3. **Write the script** from this skeleton:

```python
"""Figure F<N> (chapter <k>, sec-<id>): what it shows, and which formulas it checks.
Writes part<k>-.../figures/<name>.svg."""
import numpy as np
import style
from style import plt, report
style.setup()
# ---- checks: one report(name, bool) per formula that the figure plots ----
# ---- figure ----
fig, ax = plt.subplots(figsize=(7.0, 3.0))
...
fig.tight_layout()
style.save(fig, "<name>", style.OUT_PART1)   # add OUT_PARTn to style.py for other Parts
raise SystemExit(0 if style.all_ok() else 1)
```

4. **Look at it.** Convert and open the image:
   `rsvg-convert -b white -w 900 fig.svg -o /path/scratch.png` (also `-b '#212529'` for the dark theme), then Read the PNG. Make a contact sheet scaled to 520 px per figure (PIL) to judge readability at PDF scale. Do not skip this: it caught overlapping labels, literal `\n` strings and tiny text.
5. **Embed it** with a figure attribute list, alt text and a pointer from the running text:
   `![Caption.](figures/<name>.svg){#fig-<topic> fig-alt="One sentence a screen reader can use."}`
   Add `(@fig-<topic>)` in a nearby sentence so the figure is referenced. Insert by Python replacement with an assertion that the marker string occurs exactly once; it stopped two misplaced insertions.
6. **Add the README row** (script, figure id, book location, what is checked, runtime).
7. **Verify** (next section), then commit the script, the SVG, the `.qmd` and the README together.

## Verification (all four, in this order)

- Run every figure script: each must exit 0 with no `FAIL`. (52 checks in the first batch.)
- `quarto render <chapter>.qmd --to html` for each edited chapter (about 3 s each). Grep the log for `warning|unable to resolve`, and the HTML for `?@fig` (unresolved reference) and for the figure id.
- Build the PDF: `quarto render <any one chapter>.qmd --to pdf`. In this book project that builds the whole book (about 60 s, 293 pages) into `_site/affine-toda.pdf`, which is git-ignored. Find the page with `pdftotext -f N -l N` plus a caption phrase (join lines with `tr '\n' ' '`, because captions wrap), render it with `pdftoppm -r 80 -png`, and look at it.
- Confirm the SVG files reached `_site/`, and that they are small (target tens of kB; see the contour pitfall).

## Pitfalls (each of these happened)

- **Figures are shown at text width.** The PDF text block is about 5.4 inches wide. A 10–11 inch figure with 10 pt text arrives at 5 pt and is unreadable. Design at about 7 inches and keep text at 8–9 pt (7–7.5 pt for annotations). Panels with long text labels must be shortened or split; move formulas into the caption.
- **Mermaid layout is not under your control.** A `flowchart LR` with `direction TB` subgraphs rendered as four huge stacked blocks that spilled over two PDF pages, in the default yellow subgraph style. The plain part-level graph worked. Use mermaid only for simple flows; otherwise draw with matplotlib (the reading-paths figure) so that you control proportions.
- **matplotlib 3.6 mathtext is limited.** `\mathcal` needs braces (`\mathcal{E}`); there is no `\#`; there is no `\operatorname` for new functions (use `\mathrm{sech}`). Strings that mix `r"..."` and `"\n"` pieces are fine, but see the next point.
- **Escaping when you rewrite scripts with Python.** Replacement strings written as `"\\\\n"` or `"\\\\int"` land in the file as two backslashes, which prints a literal `\n` or breaks mathtext. After any scripted edit, grep the scripts for `\\\\`. For anything longer than a line, rewrite the file whole with the Write tool. Chained `str.replace` calls that overlap (a font-size remap followed by a pattern containing the old size) corrupt each other: apply one-pass regex maps, or rewrite.
- **Filled contours explode the SVG.** `contourf` on a 600×600 grid made a 2.4 MB SVG. In matplotlib 3.6 `rasterized=True` is ignored by `contourf`; use `for c in cs.collections: c.set_rasterized(True)` and `style.save(..., dpi=250)` (already set). Result: 139 kB.
- **Sorting complex numbers is not a robust comparison.** `np.sort_complex` reorders conjugate eigenvalues with equal real parts up to rounding, and a correct formula failed its check. Compare sets by nearest distance instead.
- **The Write tool refuses a file changed since you last read it,** including by your own `sed`. Read it again first.
- **Do not `pkill -f` a name that appears in your own command line.** It matched the shell running the command (exit 144). Match the process by its exact command via `/proc/<pid>/cmdline`, and note that a script may already have finished.
- **Unbuffered output for long scripts.** `python3 -u` and a redirect to a scratch file, or you see nothing until the end.

## Git and coordination

- `origin/main` moves while you work (PRs #8 and #9 merged during the Part I figures). Before merging: `git fetch`, check `git log HEAD..origin/main` and the overlapping files, fast-forward local `main` to `origin/main`, then merge your branch. Re-render the touched chapters on the merged tree before pushing. Commit and push only when the user asks, and end commit messages with the attribution line given in the session.
- Text fixes from the reading briefs and figures touch the same chapters. Keep them in separate commits.
- Do not touch `_quarto.yml`, `_site/` (ignored) or other agents' briefs.

## Where figures would help in Parts II–VII

These are suggestions drawn from reading Parts I–VI. Verify each against the current chapter text before drawing, and give each one a check. I have not read Part VII closely, so its entries are guesses from its section titles.

| Part | Figure idea | What a check could be |
|---|---|---|
| II | Root systems of $A_2$ and $G_2$ with the Weyl chambers, and the Coxeter-plane polygons of Dorey's orbits | orbit sizes $=h$; roots are the stated lengths |
| II | Affine Dynkin diagrams with Kac labels, and the folding of $a_3^{(1)}\to c_2^{(1)}$ and $d_4^{(1)}\to g_2^{(1)}$ | labels sum to $h$; folded labels $1,1,2$ and $1,2,3$ |
| II | The six-vertex weights and the tensor-product graph for a small case | eigenvalues from `eq-tpg-rule` against `sixvertex.py` |
| III | Time delay in a two-soliton collision with the interaction coefficient $A_{ab}$; the weights of $\Lambda^2\mathbb C^4$ that a single soliton does and does not reach (missing charges) | $0<A_{ab}<1$; $h/\gcd(a,h)$ of $\binom ha$ weights |
| IV | Poles of $S_{ab}$ in the physical strip with bound-state labels ($a_2$, $e_8$); $B(\beta)$ and its duality; $H$ against $B$ for a dual pair | $B(4\pi/\beta)=2-B$; fusing angles give the stated masses |
| V | Thimbles of $S=\frac{x^2}2+\frac{ix^3}3$ in the complex plane; zeros and poles of transmission factors in the $z$-plane next to the spectrum | critical points $0$ and $i$; $S(i)=-\frac16$; the $a_3^{(1)}$ heat trace closed form |
| VI | Pole towers in the strip variable $t$ for sine-Gordon, $c_n^{(1)}$ and $a_{2n-1}^{(2)}$; $H_s$ against $H_b$ for $g_2^{(1)}$ and $f_4^{(1)}$; a map of the $(p,r)$ plane showing where all quantum dimensions are positive | positivity iff $r\equiv\pm1\pmod{2p}$ (a script exists: `rsos_qdims.py`, and the sign criterion was checked for $p\le30$) |
| VII | The collapse thresholds of the Coulomb gas on the coupling axis; the pairing $(s,1/\bar s)$ of Krein-unitary eigenvalues; a map of the sector status | thresholds from the text; pairing from the transfer-matrix scripts in `foundations_code/` |

## Not done

- Only F2 was viewed on the dark theme; the others use the same palette but were not looked at there.
- The "full set" ideas for Part I were left out: fluctuation-operator potentials, the $B$–$\xi$ continuation line, RSOS height paths, Lee–Yang zeros, the Bethe–Yang circle, a map of the validity criteria C1–C6 and a Part I concept map.
- The figures use Part I's symbols as they stood on 2026-10-11. If notation changes again, rerun the scripts and look at the results.
