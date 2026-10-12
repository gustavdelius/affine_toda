# Affine Toda field theory

A graduate textbook and research monograph on affine Toda field theories: two-dimensional quantum field theories with one theory for each affine Kac–Moody algebra. They are integrable, so many quantities can be computed exactly, and they display particles and bound states, solitons and breathers, strong–weak duality, quantum-group symmetry, integrable boundaries and, at imaginary coupling, non-unitarity.

The book is written for a reader who knows quantum field theory at the level of a first graduate course, and it leads from sine-Gordon theory to current research.

**Read the book:** [gustavdelius.github.io/affine_toda](https://gustavdelius.github.io/affine_toda/) · [PDF](https://gustavdelius.github.io/affine_toda/affine-toda.pdf)

Authors: Gustav Delius, Niall MacKay and Gerard Watts.

## Status

> **Most of the text and the supporting code were written by Claude, an AI model made by Anthropic**, working from outlines, source papers and briefs supplied by the authors. The authors directed the mathematics and physics and are responsible for it, and the text was checked in repeated review passes by the authors and by further Claude instances acting as referees. The prose, many of the proofs and the numerical checks are machine-drafted, and errors may remain. Where a result matters to you, check the proof, and run the scripts.

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
| IV. Real coupling | 12–15 | Quantum theory at real coupling: exact particle S-matrices for every algebra |
| V. Semiclassics | 16–19 | Semiclassical quantization of complex solitons (Picard–Lefschetz, one-loop masses) |
| VI. Soliton S-matrices | 20–27 | Exact soliton S-matrices, including $c_n^{(1)}$, $a_{2n-1}^{(2)}$, $e_6^{(2)}$ and the open case $f_4^{(1)}$ |
| VII. Half-line | 28–35 | Integrable boundaries, boundary solitons, boundary quantum groups, quantum symmetric pairs, reflection amplitudes |
| VIII. Defects | 36–43 | Integrable defects: transmission of particles and solitons, defect conditions, the underlying quantum-group structure, relation to boundaries |
| IX. Foundations | 44–49 | Whether the imaginary-coupling theories exist as quantum field theories (in progress) |

Chapter 1 is the introduction and reader's guide, and the outlook is chapter 50. The guide sets out five reading paths: a course on integrable field theory, a course on quantum groups in field theory, integrable boundaries, non-unitary quantum field theory, and a path for specialists through Parts V–IX.

The book's own summary of new material is in the introduction, and [Appendix D](appendices/D-corrections.qmd) lists every correction the book makes to published papers.

## Repository layout

| Path | What it contains |
|---|---|
| `index.qmd`, `part*/`, `outlook/`, `appendices/` | The book, as Quarto source (`.qmd`) with bibliography `references.bib` |
| `figures/`, `part*/figures/` | The committed SVG figures |
| `textbook_code/` | Symbolic and numerical checks of the textbook material, Parts I–VII and appendix A |
| `textbook_code/figures/` | Scripts that draw the figures and check the formulas they plot |
| `foundations_code/` | Research computations behind Parts V and IX |
| `soliton_rmatrix_code/` | Construction of quantum-affine representations and R-matrices, soliton S-matrices and one-loop tests (Part VI) |
| `briefs/` | Referee reports, reading-fix briefs and open research tasks |
| `add-backlinks.py`, `part-preamble.lua`, `status-tags.lua`, `theme*.scss` | Quarto build helpers and styling |
| `publish-pdf.sh` | Builds the PDF locally and attaches it to the `pdf-latest` release |
| `requirements.txt` | Python packages needed by the scripts |
| `LICENSE` | MIT licence for the code |

[Appendix C](appendices/C-code.qmd) says which directory serves which part of the book. Each code directory has its own README with a file-by-file description.

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

The GitHub Actions workflow (`.github/workflows/publish.yml`) renders the HTML on each push to `main` and deploys it to GitHub Pages. It does not render the PDF, because Quarto hangs on the Mermaid figure in chapter 1 on the runner. The PDF is built locally with `./publish-pdf.sh`, which uploads it to the `pdf-latest` release and redeploys the site. That script needs `gh` logged in to GitHub.

## Running the supporting code

The scripts are Python 3 and use numpy, scipy, sympy, mpmath and matplotlib (for the figures). Install them with:

```bash
python3 -m pip install -r requirements.txt
```

- **`textbook_code/`** — each script runs on its own and prints one `PASS` or `FAIL` line per check. For example:
  ```bash
  cd textbook_code && python3 sine_gordon.py
  ```
  Runtimes range from seconds to a few minutes; the table in [`textbook_code/README.md`](textbook_code/README.md) maps each script to the equations and results it checks.
- **`textbook_code/figures/`** — each script writes an SVG into the book tree and checks the formulas it plots, printing `PASS` or `FAIL` per check. The SVGs are committed, so building the book does not require running them. To regenerate one, run `python3 textbook_code/figures/<name>.py` from the repository root.
- **`foundations_code/`** and **`soliton_rmatrix_code/`** — run each script from its own directory, as the directory READMEs describe. Some scripts need `rsolve.py` or `dtw.py` copied alongside them, and the READMEs say which.

Some scripts take a long time, and a few depend on a sign convention for $q$. The `QSIGN` environment variable selects the physical sign, and the soliton and foundations READMEs explain the default.

## Contributing and feedback

The book was produced by an iterative process of drafting, refereeing and fixing, recorded in `briefs/`. If you find an error, the most useful report is a specific one: the chapter and equation label (`eq-…`, `prp-…`, `sec-…` ids in the source), what is wrong, and, where possible, a script that shows it. Errors in the published literature that the book has found are recorded in Appendix D, and that list should be kept in step with the chapter notes.

## License

- **Code** (the Python scripts, the Lua filters and the build scripts) is released under the [MIT License](LICENSE).
- **Book text and figures** (the `.qmd` sources, the rendered book and the SVG figures) are released under the [Creative Commons Attribution 4.0 International licence](https://creativecommons.org/licenses/by/4.0/legalcode). You may share and adapt them, including commercially, provided you give appropriate credit to the authors.

Copyright is held by Gustav Delius, Niall MacKay, Gerard Watts and contributors to this repository.
