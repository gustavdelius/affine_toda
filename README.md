# Affine Toda field theory

A unified review of affine Toda field theories, with new results: two-dimensional quantum field theories with one theory for each affine Kac–Moody algebra. They are integrable, so many quantities can be computed exactly, and they display particles and bound states, solitons and breathers, strong–weak duality, quantum-group symmetry, integrable boundaries and, at imaginary coupling, non-unitarity.

The book collects results scattered over five decades of papers and presents them in a logical order, in one notation and with one set of conventions. It is a reference, not a course; companion lecture notes teach parts of the subject step by step. One course accompanies Part VII on integrable boundaries: *Lectures on affine Toda field theory on the half-line* (`lecture-notes/half-line/`, drafted by Gemini). Two courses accompany Part VIII on integrable defects: *Integrable defects, from Bäcklund to Baxter* (`lecture-notes/defects/`, drafted by Claude) and *Lectures on integrable defects in affine Toda field theory* (`lecture-notes/defects-gemini/`, drafted by Gemini). The home page of the site says what the site holds and compares the courses.

**Read the book and the lecture notes:** [gustavdelius.github.io/affine_toda](https://gustavdelius.github.io/affine_toda/) · [PDF of the book](https://gustavdelius.github.io/affine_toda/affine-toda.pdf)

**Coordinating editor:** Gustav Delius

**Thanks** to Niall MacKay, Gerard Watts and Ed Corrigan for discussions and feedback. The list will grow as more readers respond. Being thanked does not imply having reviewed or endorsed the book.

## Status

> **Most of the text and the supporting code were drafted by Claude, an AI model made by Anthropic**, working from outlines, source papers and briefs assembled by the coordinating editor. Gustav Delius coordinates the project and its editorial process. Further Claude, Codex and Gemini instances have been used for automated critique and consistency checks, not as human referees. Validation is by AI agents only; human readers contribute by giving feedback. The prose, many of the proofs and the numerical checks are machine-drafted, and errors may remain. Where a result matters to you, check the proof, and run the scripts.

Each main statement from Part V on carries a status label:

| Label | Meaning |
|---|---|
| **[Theorem]** | Complete proof given in the book but without outside check |
| **[Proposition, sketch]** | Proof outline that a specialist could complete |
| **[Established]** / **[Established; reformulated]** | Proved or computed in the cited literature, possibly restated |
| **[Numerical]** | Checked numerically by a script in this repository, not proved |
| **[Conjecture]** | Stated with its supporting evidence |

The record of what has been checked, and how, is in [`briefs/`](briefs/) (for example `briefs/referee_textbook.md`, and `briefs/referee_half_line_*.md` for Part VII).

## Contents

| Part | Chapters | Subject |
|---|---|---|
| I. Physical background | 1–4 | Factorized scattering, sine-Gordon theory, non-Hermitian quantum mechanics |
| II. Algebraic background | 5–8 | Root systems, affine Kac–Moody algebras, quantum groups, R-matrices |
| III. Classical theory | 9–11 | Lagrangians, integrability, solitons, reality of their energies |
| IV. Real coupling | 12–16 | Quantum theory at real coupling: exact particle S-matrices for every algebra, and the theories based on Lie superalgebras |
| V. Semiclassics | 17–20 | Semiclassical quantization of complex solitons (Picard–Lefschetz, one-loop masses) |
| VI. Soliton S-matrices | 21–28 | Exact soliton S-matrices, including $c_n^{(1)}$, $a_{2n-1}^{(2)}$, $e_6^{(2)}$ and the open case $f_4^{(1)}$ |
| VII. Half-line | 29–36 | Integrable boundaries, boundary solitons, boundary quantum groups, quantum symmetric pairs, reflection amplitudes |
| VIII. Defects | 37–44 | Integrable defects: scattering axioms, defect conditions as Bäcklund transformations, solitons through defects, Borel quantum-group symmetry, transmission amplitudes, relation to boundaries |
| IX. Foundations | 45–50 | Whether the imaginary-coupling theories exist as quantum field theories (in progress) |

An unnumbered preface comes first. Chapter 1 is the introduction and reader's guide, and the outlook is chapter 51. The guide sets out six reading paths: integrable field theory, quantum groups in field theory, integrable boundaries, integrable defects, non-unitary quantum field theory, and a path for specialists through Parts V–IX.

The book's own list of its original results is [Appendix E](appendices/E-original-results.qmd), and [Appendix D](appendices/D-corrections.qmd) lists every correction the book makes to published papers.

## Repository layout

| Path | What it contains |
|---|---|
| `index.qmd`, `preface.qmd`, `01-introduction.qmd`, `part*/`, `outlook/`, `appendices/` | The book, as Quarto source (`.qmd`) with bibliography `references.bib`. `index.qmd` is the home page of the site, in HTML only: it says what the site holds and compares the lecture courses. `preface.qmd` is the preface |
| `figures/`, `part*/figures/` | The committed SVG figures |
| `lecture-notes/` | Lecture notes that accompany the book, each a Quarto book of its own: `lecture-notes/half-line/` (on Part VII, drafted by Gemini), `lecture-notes/defects/` and `lecture-notes/defects-gemini/` (both on Part VIII; the second was drafted by Gemini), with `lecture-notes/book-links.lua`, which turns `[text](book:<id>)` into links to the published book, and `lecture-notes/_navbar.yml`, the navbar the courses share |
| `code/` | All the book's check scripts. Flat scripts check Parts I–IV, VII and VIII, appendix A and parts of Part IX; subdirectories hold the larger computations |
| `code/figures/` | Scripts that draw the figures and check the formulas they plot |
| `code/lectures/` | Check and figure scripts of the lecture notes; see its `README.md` |
| `code/hirota/`, `oneloop/`, `thimbles/` | Hirota fluctuations, Jordan types, one-loop soliton masses and thimbles (Part V) |
| `code/rmatrix/` | Construction of quantum-affine representations and R-matrices, soliton S-matrices and one-loop tests (Part VI, appendix B) |
| `code/sectors/`, `irf/`, `vernier/`, `review614/` | Soliton sectors and Krein unitarity, RSOS kinks, the Vernier–Jacobsen–Saleur lattice model and the continuum-limit lemmas (Parts VI and IX) |
| `briefs/` | Referee reports, reading-fix briefs and open research tasks |
| `add-backlinks.py`, `add-discovery-metadata.py`, `home-page.lua`, `part-preamble.lua`, `status-tags.lua`, `theme*.scss` | Quarto build helpers, discovery metadata and styling |
| `publish-pdf.sh` | Builds the PDF locally and attaches it to the `pdf-latest` release |
| `check-refs.py`, `render-scratch.sh` | Checks cross-references and citations without rendering; renders the book in a scratch copy and reports unresolved references |
| `AGENTS.md`, `CLAUDE.md`, `.claude/skills/` | Working rules for the AI agents that draft and referee the book, and the Claude Code skills for referee briefs and their integration |
| `requirements.txt` | Python packages needed by the scripts |
| `LICENSE` | MIT licence for the code |

[Appendix C](appendices/C-code.qmd) describes the code for the reader and says which directory serves which part of the book. `code/README.md` is the index of every script, and `code/rmatrix/README.md` describes the scripts of Part VI.

## Building the book

You need [Quarto](https://quarto.org). The HTML build needs nothing else.

```bash
quarto render              # HTML site, written to _site/
quarto preview             # live preview while editing
```

The PDF needs a LaTeX installation:

```bash
quarto render --to pdf
```

Two helpers check the book while it is being edited:

```bash
python3 check-refs.py      # ids, cross-references, citation keys, status labels and the list of original results, in under a second
./render-scratch.sh        # render in a scratch copy, so a running preview is undisturbed (--pdf for the PDF too)
```

The lecture notes in `lecture-notes/` are separate Quarto projects, which the book's render skips. Check and render them with

```bash
python3 check-refs.py lecture-notes/defects
```

```bash
./render-scratch.sh --project lecture-notes/defects
```

and likewise for `lecture-notes/defects-gemini` and `lecture-notes/half-line`.

The GitHub Actions workflow (`.github/workflows/publish.yml`) renders the HTML on each push to `main` and deploys it to GitHub Pages, with the lecture notes of each `lecture-notes/<name>/` under `lectures/<name>/`. It does not render the PDF, because Quarto hangs on the Mermaid figure in chapter 1 on the runner. The PDF is built locally with `./publish-pdf.sh`, which uploads it to the `pdf-latest` release and redeploys the site. That script needs `gh` logged in to GitHub.

## Running the supporting code

The scripts are Python 3 and use numpy, scipy, sympy, mpmath and matplotlib (for the figures). Install them with:

```bash
python3 -m pip install -r requirements.txt
```

- **`code/`** — each script runs on its own and prints one `PASS` or `FAIL` line per check. For example:
  ```bash
  cd code && python3 sine_gordon.py
  ```
  Runtimes range from seconds to a few minutes; the table in [`code/README.md`](code/README.md) maps each script to the equations and results it checks.
- **`code/figures/`** — each script writes an SVG into the book tree and checks the formulas it plots, printing `PASS` or `FAIL` per check. The SVGs are committed, so building the book does not require running them. To regenerate one, run `python3 code/figures/<name>.py` from the repository root.
- **Subdirectories of `code/`** — run each script from its own directory, as [`code/README.md`](code/README.md) and the subdirectory READMEs describe. Many print tables rather than `PASS`/`FAIL` lines. Some need `rsolve.py` or `dtw.py` copied alongside them, and the READMEs say which.

Some scripts take a long time, and a few depend on a sign convention for $q$. The `QSIGN` environment variable selects the physical sign, and `code/README.md` and `code/rmatrix/README.md` explain the default.

## Contributing and feedback

The book was produced by an iterative process of drafting, checking and fixing, recorded in `briefs/`. If you find an error, the most useful report is a specific one: the chapter and equation label (`eq-…`, `prp-…`, `sec-…` ids in the source), what is wrong, and, where possible, a script that shows it. Errors in the published literature that the book has found are recorded in Appendix D, and that list should be kept in step with the chapter notes.

## Online publication and discovery plan

The book is a living online publication and is not intended for print publication. The HTML website is the canonical reading edition, while occasional numbered releases will provide stable versions for citation.

The repository is prepared for the following publication process:

1. **Search and sharing.** The build supplies canonical URLs, a description, subject keywords, honest editor-only citation metadata and social-card metadata. The published site includes a sitemap and a `robots.txt` file so that search engines can discover every chapter.
2. **Stable releases.** Substantial milestones will be tagged with numbered versions. Each release will preserve the PDF and source corresponding to that version.
3. **DOI and preservation.** The first stable release will be deposited in Zenodo as a *Publication — Book*, with Gustav Delius identified as editor. Zenodo's concept DOI will identify the evolving book; each archived release will also have its own version DOI. The reserved DOI must be added to `_quarto.yml`, the citation page and the PDF before the first deposit is published.
4. **Scholarly indexes.** After the DOI has been registered, the book will be submitted to INSPIRE as a book. The DOI and book record can also be added to relevant ORCID and institutional-repository records.
5. **Feedback.** Readers who give feedback are thanked in the introduction (`sec-about`) and in this README.

The [citation page](citation.qmd) gives the current editor-only citation in several formats. Google Search Console and Bing Webmaster Tools should be connected after deployment and given the generated sitemap; those steps require the site owner's accounts and cannot be configured in this repository alone.

## License

- **Code** (the Python scripts, the Lua filters and the build scripts) is released under the [MIT License](LICENSE).
- **Book text and figures** (the `.qmd` sources, the rendered book and the SVG figures) are released under the [Creative Commons Attribution 4.0 International licence](https://creativecommons.org/licenses/by/4.0/legalcode). You may share and adapt them, including commercially, provided you give appropriate credit to the project and its contributors.

Copyright in individual contributions remains with the applicable rights holders. Editorial credits and acknowledgements do not by themselves imply copyright ownership.
