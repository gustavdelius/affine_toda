# Brief: referee pass on Part III (chapters 9–11), phase 1 of 2

You are an adversarial mathematical-physics referee. You are working in the repository of a Quarto book on affine Toda field theory (worktree `.claude/worktrees/referee-briefs-part-iii`, branch `claude/referee-briefs-part-iii`).

## How this pass works

This pass has two phases.

- **Phase 1 (this brief).** You decide what to check. Read the chapters cold and write your own ranked list of the claims most at risk. Then check as many as you can, most important first.
- **Phase 2.** When you report phase 1 complete, the editor sends you a second list of items, written by the coordinator who drafted these chapters. You then check what your own list missed, and compare the two lists.

The book was drafted by Claude, and earlier referee passes on these chapters were also run by Claude. Errors that Claude's own checklist would not find are the ones this pass is for. So phase 1 must rest on your judgment alone, and this brief deliberately does not tell you where to look.

**Rules.**
- Other agents may be editing files in this worktree at the same time.
- Do not edit any existing file.
- Do not run git commands that change the index or the working tree (no commit, checkout, stash, reset, rebase or restore).
- Do not run `quarto render`.
- Do not spawn subagents.
- Write scratch scripts only under `referee_scratch/classical_09_11/`.
- Write the phase-1 report to `briefs/referee_classical_09_11_phase1_report.md`. Do not change it after you have reported phase 1 complete.
- **Until phase 1 is complete,** do not read the earlier referee brief and report that cover these chapters, or the reading pass's brief: `briefs/referee_textbook.md`, `briefs/referee_textbook_report.md`, `briefs/referee_textbook_scripts/` and `briefs/part3_reading_fixes.md`. They contain earlier reviewers' judgments. Do not look for the phase-2 list. You may read everything else, including the git history.

## Context

Part III, "Classical affine Toda field theory" (`part3-classical/`), defines the theories the book studies and their classical properties:

- **Chapter 9**, `part3-classical/09-lagrangian.qmd` (`sec-theories`): Lagrangians, Lax pair, conserved charges, masses and couplings, Dorey's fusing rule.
- **Chapter 10**, `part3-classical/10-solitons.qmd` (`sec-classical-solitons`): Hirota solitons, vertex operators, topological charges, soliton masses and Lie duality, multi-solitons and time delays, fusion, breathers.
- **Chapter 11**, `part3-classical/11-reality.qmd` (`sec-classical`): the antilinear symmetry $\Theta$, reality of energy and charges, singular solutions, the missing-charge problem.

The part page is `part3-classical/part.qmd`. The figures and their captions are part of the text under review.

**Statements and labels.** Two statements carry labels: `prp-dorey-rule` (**[Established]**) and `prp-hirota-energy` (**[Established; reformulated]**). The labels are defined in `README.md`, section "Status labels". The book lists its claims of originality in Appendix E (`appendices/E-original-results.qmd`); one item there concerns Part III. Appendix D (`appendices/D-corrections.qmd`) lists the book's corrections to published papers; it has none for Part III.

**History.** These are facts from the git log; use them as you see fit.
- The chapters were drafted on 2026-10-09 (commit de48bf7).
- One referee pass, by Claude agents, was run on 2026-10-10 as part of a book-wide pass. Its replacement text was integrated in commit 0e57166.
- A reading pass then edited the chapters (790bedf). Figures were added (4321e19). Small edits followed (e6a0dbc, 1834767, cf4597d).
- `git diff 0e57166 HEAD -- part3-classical/` shows the text that no referee has seen.
- This is the first pass on Part III by a model other than Claude.

**Check scripts.** Each prints one PASS or FAIL line per check. They are the drafter's checks, not evidence.
- `code/toda_classical.py` (11 checks, about 3 minutes). It is cited in chapters 9 and 10; its row in `code/README.md` lists what it checks.
- Four figure scripts, each 3–14 s: `code/figures/ch9_f13_dorey_d4.py` (10 checks), `ch10_f14_topological_charges.py` (4), `ch10_f15_two_soliton_shift.py` (7) and `ch11_f16_breather_theta.py` (4).

  The figure scripts overwrite the committed SVGs in `part3-classical/figures/` whatever `OUT_DIR` says. If you run them, report it, and do not restore the files yourself; the editor will.

## Notation

- **Coordinates.** Metric signature $(+,-)$, with $\partial_\pm=\partial_t\pm\partial_x$, so that $\partial_+\partial_-=\partial_t^2-\partial_x^2$. These are twice $\partial/\partial x^\pm$ for $x^\pm=t\pm x$.
- **Roots.** The longest roots of the finite algebra have squared length 2, and $\vec\alpha^\vee=2\vec\alpha/\vec\alpha^2$. The Kac labels $n_j$ satisfy $\sum_jn_j\vec\alpha_j=0$ with $n_0=1$, and $h=\sum_jn_j$. The Cartan matrix is $a_{ij}=\vec\alpha_i^\vee\cdot\vec\alpha_j$. Some of the literature gives the longest root of a twisted algebra $X^{(k)}$ squared length $2k$; then $\beta_{\rm lit}^2=\beta^2/k$.
- **Potentials** (`eq-imaginary-potential`): $V_{\rm real}=\frac{m^2}{\beta^2}\sum_jn_j(e^{\beta\vec\alpha_j\cdot\vec\phi}-1)$ and $V_{\rm imag}=-\frac{m^2}{\beta^2}\sum_jn_j(e^{i\beta\vec\alpha_j\cdot\vec\phi}-1)$, with $\beta>0$. For sine-Gordon, $\beta_{\rm SG}^2=2\beta^2$ and the particle mass is $m_0=2m$ (`sec-sg-dictionary`).
- **Symbols that clash with other parts of the book.**
  - In chapters 10–11, $\omega=e^{2\pi i/h}$ is a root of unity. It is *not* the coupling $\omega=4\pi/\beta^2-1$ used elsewhere.
  - $\xi$ is a soliton's position constant, not the sine-Gordon $\xi$.
  - $E$ is the exponential of `eq-an-soliton`.
- **Conjugate species and $\Theta$.** For $a_n^{(1)}$, $\bar a=h-a$. $\Theta$ acts on fields by $\vec\phi\mapsto-\bar{\vec\phi}$.
- **Breathers.** Chapter 10 writes the constituent rapidities as $\theta_0\pm iu$; chapter 3 writes $\pm i(\frac\pi2-u)$.
- **The quantum $q$,** mentioned in passing in chapter 10: $q=e^{-4\pi^2i/\beta^2}$ with the coproduct of Part VI, and $q^{-1}$ with that of chapter 7 (`eq-q-physical`).
- **Time delays.** $\Delta t=t_{\rm out}-t_{\rm in}$, as defined in `sec-multi-solitons`.

## Read first

- `01-introduction.qmd`, section `sec-conventions`; `README.md`, section "Status labels".
- The background that Part III uses:
  - chapter 2, `eq-charge-eigenvalue`;
  - chapter 3, `sec-sg-classical`, `eq-sg-breather`, `eq-sg-s0`, `eq-sg-classical-mass`;
  - chapter 4, `sec-antilinear`;
  - chapter 5, `sec-weyl-coxeter`, `sec-weights`, `sec-perron-frobenius`;
  - chapter 6, `sec-affine-dynkin`, `sec-affine-folding`, `sec-gradations`, `sec-vertex-operators`.
- Read the later parts that Part III points to only as far as you need them.
- To find an object: `grep -rn '{#<id>}' part*/ appendices/ outlook/ *.qmd`. Figure ids are followed by attributes: `{#fig-x fig-alt=…}`.

**Sources.** Texts already in `~/referee-copy/papers/` (or `papers/`, if you run through `run_referee_agy.sh`) are named `<arXiv id>.txt`, with `/` replaced by `_`, and the PDFs are next to them. Fetch any other source from arXiv. These are the arXiv ids of the sources cited in Part III:
- Hollowood, hep-th/9110010 (`hollowood1992`);
- Olive–Turok–Underwood, hep-th/9305160 (`olive1993b`);
- Kneipp–Olive, hep-th/9404030 (`kneipp1996`);
- Aratyn–Constantinidis–Ferreira–Gomes–Zimerman, hep-th/9212086 (`aratyn1993`);
- MacKay–McGhee, hep-th/9208057 (`mackay1993`);
- McGhee, hep-th/9307035 (`mcghee1994`);
- Fring–Johnson–Kneipp–Olive, hep-th/9405034 (`fring1994`);
- Harder–Iskandar–McGhee, hep-th/9409035 (`harder1995`);
- Freeman, hep-th/9408092 (`freeman1995`);
- Zhu–Caldi, hep-th/9404097 (`zhu1995`);
- Khastgir–Sasaki, hep-th/9507001 (`khastgir1996`);
- Takács–Watts, hep-th/9810006 (`tw1999`);
- Delius–Grisaru–Zanon, hep-th/9201067 (`dgz1992`);
- Feigin–Frenkel, hep-th/9310022 (`feigin1993`);
- Niedermaier, hep-th/9401078 (`niedermaier1994`).

The others are pre-arXiv papers and books. Olive–Turok–Underwood's `olive1993a` (Nucl. Phys. B 401 (1993) 663) is not on arXiv. If you cannot read a source, say so and mark that citation "unclear"; do not reconstruct what it probably says.

## How to work

For each claim you check, do **one independent computation**: a derivation written out, or a script of your own. You may run the book's scripts for comparison, but they are not evidence. A verdict without its own evidence counts as "unclear". `briefs/referee_defects_38_40_report.md` shows the expected depth.

Hirota solutions are finite sums of exponentials, so they can be evaluated exactly or to high precision; use `mpmath` where cancellations are severe.

## Phase 1: your own list

1. **Read** the three chapters and the part page in full before checking anything: text, captions, exercises and notes.
2. **Write the risk list** as section 1 of the phase-1 report. It holds the claims you judge most likely to be wrong, overstated, misattributed or mislabelled, ranked, with roughly 15–30 entries.
   - For each entry give its location (id, or file and line), the claim quoted, why it is at risk, and how you will check it.
   - Write this list before you start checking, and do not reorder or delete entries afterwards. Append claims you discover later at the end, marked "added during checking".
   - What counts as risky is for you to decide. Formulas, signs and normalizations, status labels, claims of originality (search the literature for them), whether citations say what their sources say, consistency with other chapters, and the exercises are all fair game.
3. **Check** the entries in order of rank, at least the top ten, and more if you can. For each, give the verdict (correct / wrong / unclear), your own computation and its output, and, for errors, replacement text ready to paste, with file and line.
4. **List** the entries you did not check, and why.
5. **Stop.** End your phase-1 session with a final message that says "Phase 1 complete" and gives the report path. Do not start phase 2 on your own; the editor will send it.

## Phase-1 report format (`briefs/referee_classical_09_11_phase1_report.md`)

1. A header with the date, your model, and the commit you reviewed (`git rev-parse HEAD`).
2. Section 1: the risk list, as written before checking, with any later additions marked.
3. Section 2: one subsection per checked entry, in rank order.
4. Section 3: the entries you did not check, and why.
5. Section 4: other observations, labelled as unverified.

In phase 2 you will write a second report, `briefs/referee_classical_09_11_report.md`. Your phase-1 report stays as it is.
