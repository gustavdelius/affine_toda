# Brief: referee pass on the half-line lecture notes, phase 1 of 2

You are an adversarial mathematical-physics referee. You are working in the repository of a Quarto book on affine Toda field theory and its companion lecture notes (worktree `.claude/worktrees/referee-halfline-lectures`, branch `claude/referee-halfline-lectures`). If you were started through `briefs/run_referee_agy.sh`, you are in a throwaway copy of that worktree.

## How this pass works

This pass has two phases.

- **Phase 1 (this brief).** You decide what to check. Read the lectures cold and write your own ranked list of the claims most at risk. Then check as many as you can, most important first.
- **Phase 2.** When you report phase 1 complete, the editor sends you a second list of items, written by the coordinator of this pass. You then check what your own list missed, and compare the two lists.

The lectures were drafted by Gemini. The book chapters they accompany were drafted by Claude, and Claude coordinates this pass. Errors that neither model's own checklist would find are the ones this pass is for. So phase 1 must rest on your judgment alone, and this brief deliberately does not tell you where to look.

**Rules.**
- Other agents may be editing files in this worktree at the same time.
- Do not edit any existing file, and do not run figure scripts in a shared worktree (they overwrite tracked SVGs in place).
- Do not run git commands that change the index or the working tree (no commit, checkout, stash, reset, rebase or restore).
- Do not run `quarto render`.
- Do not spawn subagents.
- Write scratch scripts only under `referee_scratch/lectures_halfline/`.
- Write the phase-1 report to `briefs/referee_lectures_halfline_phase1_report.md`. Do not change it after you have reported phase 1 complete.
- **Until phase 1 is complete,** do not read the earlier referee briefs and reports on Part VII of the book, which these lectures restate: `briefs/referee_half_line_*.md`, `briefs/halfline_*` (files and directories), `briefs/referee_bsg_dirichlet_sign*` and `briefs/reader_review_part7.md`. They contain earlier reviewers' judgments. No earlier pass has reviewed the lectures themselves. Do not look for the phase-2 list. You may read everything else, including the git history.

## Context

**What is under review.** `lecture-notes/half-line/` is a standalone Quarto book, "Lectures on affine Toda field theory on the half-line" (`_quarto.yml`), of about 2300 lines:

- `index.qmd`: preface, the "five big ideas", syllabus;
- Lecture 1, `01-wall-and-euclidean-trick.qmd` (`sec-lec1-wall-euclidean`): half-line charges, the reflection equations, boundary crossing-unitarity, boundary bootstraps, boundary scaling Lee–Yang;
- Lecture 2, `02-sklyanin-and-rigidity.qmd` (`sec-lec2-sklyanin-rigidity`): the boundary potential, Sklyanin's double row, the classical $K$-matrices of $a_1^{(1)}$ and $a_2^{(1)}$, Dirichlet and soliton-preserving boundaries;
- Lecture 3, `03-mirror-images-and-vacua.qmd` (`sec-lec3-mirror-images`): mirror solitons, time delays, real-coupling vacua, linearized reflection;
- Lecture 4, `04-coideals-and-sine-gordon.qmd` (`sec-lec4-coideals`): boundary non-local charges, the coideal, the sine-Gordon $K$-matrix, the sign of $q$, the UV–IR relation, boundary bound states;
- Lecture 5, `05-a2-intertwiners-and-q-onsager.qmd` (`sec-lec5-a2-qsp`): the $a_2^{(1)}$ intertwiners, the quantum Serre relation of the coideal, quantum symmetric pairs and $q$-Onsager algebras, the Neumann and uniform boundaries;
- Lecture 6, `06-breathers-particles-and-duality.qmd` (`sec-lec6-particles-duality`): soliton scalar factors, breather fusion, particle reflection amplitudes, duality, solitonic boundaries, the $a_2^{(1)}$ bound states, a boundary Coleman–Thun diagram;
- `references.qmd`, with `references-notes.bib` for entries used only by the course (the others come from the book's `references.bib`).

The six figures in `lecture-notes/half-line/figures/` and their captions are part of the text under review.

**Relation to the book.** Each lecture accompanies chapters of Part VII of the book, "Affine Toda theory on the half-line" (`part7-half-line/`, chapters 29–36); the syllabus table in `index.qmd` gives the pairing. A link `[text](book:<id>)` points to the book object with that id (`lecture-notes/book-links.lua`). The book is a review with new results; step-by-step teaching belongs in the lecture notes (`AGENTS.md`). When a lecture and the book disagree, decide which is right. An error in the book is a finding too.

**Statements and labels.** The lectures contain no labelled statements. Their claims are made in prose, displayed equations, "blackboard calculations", figure captions and 18 exercises (`exr-lec1-…` to `exr-lec6-…`). Many of them restate statements of the book, which carry status labels; the labels are defined in `README.md`, section "Status labels". The book lists its original results in Appendix E (`appendices/E-original-results.qmd`, paragraph "Half-line"), and its corrections to published papers in Appendix D (`appendices/D-corrections.qmd`).

**History.** These are facts from the git log; use them as you see fit.
- Gemini drafted the lectures and their check scripts on 2026-10-24 (ffb9a5c).
- LaTeX fixes for the PDF followed (c2bab20).
- On 2026-10-25 the sign of $\phi_0$ in the Dirichlet relation of boundary sine-Gordon theory was corrected in the book and in Lectures 4 and 5 (ed0585a).
- A navbar was added (f76fdb9). `git diff ffb9a5c HEAD -- lecture-notes/half-line/` shows every edit since the draft.
- This is the first referee pass on the lectures.
- The lecture text under review is that of `main` at 035351d. Record `git rev-parse HEAD` in your report header.

**Check scripts.** Six scripts, written by the drafter from the lectures' definitions, check the lectures and draw their figures: `code/lectures/halfline_l1_f1_reflection_crossing.py`, `halfline_l2_f1_sklyanin_rigidity.py`, `halfline_l3_f1_mirror_images.py`, `halfline_l4_f1_coideal_sg.py`, `halfline_l5_f1_a2_serre.py` and `halfline_l6_f1_particles_duality.py`. They print 5, 5, 4, 4, 4 and 4 PASS lines and take 1–10 s each. Their rows in `code/lectures/README.md` say what they check.
- They are the drafter's checks, not evidence.
- They overwrite the committed SVGs in `lecture-notes/half-line/figures/` whatever `OUT_DIR` says. Run them only in a throwaway copy, and say so in your report if you did.
- The book's own scripts for Part VII are in `code/` (index: `code/README.md`), for example `code/boundary_sine_gordon.py`, `code/boundary_coideal.py`, `code/boundary_vacua.py`, `code/particle_reflection.py` and `code/qsp_checks.py`. They are not evidence either.

## Notation

- **The half-line.** The field lives on $x\le0$, with the wall at $x=0$ (`sec-conventions`). In multi-particle states and tensor products the factors are in spatial order, so the wall is at the right end.
- **Coordinates and chirality.** Metric signature $(+,-)$. The book writes $\partial_\pm=\partial_t\pm\partial_x$, which is twice $\partial/\partial x^\pm$ for $x^\pm=t\pm x$; Lecture 1 writes $x^\pm=\frac12(t\pm x)$. Left-movers are functions of $x^+$ (`AGENTS.md`, "Chirality and spin"). Bernard and LeClair label chiralities the other way.
- **Roots and potentials.** The longest roots have squared length 2. The Kac labels $n_j$ satisfy $\sum_jn_j\vec\alpha_j=0$ with $n_0=1$, and $h=\sum_jn_j$. For $a_n^{(1)}$, $h=n+1$, every $n_j=1$ and $\vec\alpha_j^2=2$, and $\bar a=h-a$. The potentials are $V_{\rm real}=\frac{m^2}{\beta^2}\sum_jn_j(e^{\beta\vec\alpha_j\cdot\vec\phi}-1)$ and $V_{\rm imag}=-\frac{m^2}{\beta^2}\sum_jn_j(e^{i\beta\vec\alpha_j\cdot\vec\phi}-1)$, with $\beta>0$ (`eq-imaginary-potential`). For sine-Gordon, $\beta_{\rm SG}^2=2\beta^2$ and the particle mass is $m_0=2m$.
- **Boundary potential** (`eq-boundary-potential`). At real coupling $\mathcal B=\frac{2m}{\beta^2}\sum_jC_j\sqrt{2n_j/\vec\alpha_j^2}\,\big(e^{\beta\vec\alpha_j\cdot\vec\phi/2}-1\big)$. $C_j\equiv0$ is the Neumann boundary and $C_j\equiv+1$ the uniform boundary. The book's $C_j$ is Delius and Gandenberger's, and $-\varepsilon_j$ of Delius (1998). The book calls a sign pattern solitonic if $\prod_jC_j=1$ and not all $C_j=+1$ (`sec-real-coupling-vacua`).
- **Quantum group.** $\omega=\frac{4\pi}{\beta^2}-1$. For sine-Gordon, Lecture 4's $\lambda=\frac{8\pi}{\beta_{\rm SG}^2}-1$ is the same number. The physical deformation parameter is $q=-e^{-i\pi\omega}=e^{-4\pi^2i/\beta^2}$ (`eq-q-physical`), with the coproduct of Part VI: $\Delta(e_j)=e_j\otimes1+k_j\otimes e_j$, $\Delta(f_j)=f_j\otimes k_j^{-1}+1\otimes f_j$, $\Delta(k_j)=k_j\otimes k_j$. Chapter 7 uses the opposite coproduct, so $R_{\rm VI}(x;q)=R_7(x;1/q)$. The map $q\to-q$ is a Klein-factor twist, not a gauge.
- **The coideal.** $b_j=e_j+q_j^{-1}f_jk_j+\hat\epsilon_jk_j$, with $\Delta(B_{\hat\epsilon})\subset U_q(\hat{\mathfrak g})\otimes B_{\hat\epsilon}$, and $\epsilon_*^2=1/\big((1-q)(1-q^{-1})\big)$.
- **Spectral parameters and blocks.** $y=e^{\omega\theta}$ and $u=-i\theta$; $T=\frac{h\omega}2$ and $\mu=-\frac{iT\theta}\pi$. The real-coupling blocks are $(x)=\sinh(\frac\theta2+\frac{i\pi x}{2h})/\sinh(\frac\theta2-\frac{i\pi x}{2h})$ (`eq-blocks`), and $(x)_I$ are the imaginary-coupling blocks. $f_x(\theta)=\frac{\sinh\theta+i\sin\pi x}{\sinh\theta-i\sin\pi x}$. The continuation to real coupling is $\omega\to-2/B$, with $B=\frac{2\beta^2}{4\pi+\beta^2}\in(0,2)$; weak–strong duality is $B\to2-B$.
- **Boundary sine-Gordon.** The book uses the parameters $(\eta,\vartheta)$ of Bajnok, Palla, Takács and Tóth (`eq-bsg-eta-theta`), not Ghoshal and Zamolodchikov's $(\xi,k)$, whose $\xi$ clashes with the book's $\xi=\pi/\lambda$ and with a soliton's position $\xi$. On $x\le0$ the Dirichlet relation is $\eta=-4\pi\phi_0/\beta_{\rm SG}$ (`eq-bsg-dirichlet`). Ghoshal and Zamolodchikov's (5.27) has the sign for a bulk on $x\ge0$.
- **Symbols used for more than one thing in the lectures.**
  - $s$: a spin (Lecture 1); $\sin\frac{\pi a}h/\cosh\theta$ (Lecture 3); $C_0C_1C_2$ (Lecture 2); $\prod_j\hat\epsilon_j/\epsilon_*$ (Lectures 5 and 6).
  - $\lambda$: Sklyanin's spectral parameter (Lecture 2); the sine-Gordon $\omega$ (Lecture 4); with a subscript, the boundary coupling $\lambda_b$.
  - $\zeta=e^{2\pi i/h}$ (Lecture 3); chapter 10 of the book calls it $\omega$.
  - $\sigma$: an anti-automorphism, a singular value, the scalar factor $\sigma(\eta,u)$.
  - $B$: the boundary operator, the coupling, the coideal $B_{\hat\epsilon}$.
  - $T$: a monodromy matrix, the topological charges $T_j$, $\frac{h\omega}2$.
  - $\epsilon_k$ (signs, Lecture 3), $\epsilon_j$ (coefficients of the boundary perturbation, Lecture 4), $\hat\epsilon_j$ and $\epsilon_*$.
- **Roots of unity.** $q$ is a root of unity exactly when $\omega$ is rational, and representations then degenerate. Test at generic couplings such as $\omega=\sqrt2$ or $e/2$, and avoid $B=\frac1{2k}$ in residue checks.

## Read first

- `AGENTS.md`: under "Writing book text", the bullet "Lecture notes"; and the section "Conventions and known traps".
- `01-introduction.qmd`, section `sec-conventions`; `README.md`, section "Status labels".
- The book chapters that each lecture accompanies, as far as you need them:
  - Lecture 1: `part7-half-line/29-boundary-scattering.qmd` (`sec-boundary-scattering`);
  - Lecture 2: `31-boundary-conditions.qmd` (`sec-boundary-conditions`);
  - Lecture 3: `32-half-line-solitons.qmd` (`sec-half-line-solitons`);
  - Lecture 4: `30-boundary-sine-gordon.qmd` (`sec-boundary-sine-gordon`) and `33-boundary-quantum-groups.qmd` (`sec-boundary-qg`);
  - Lecture 5: chapter 33, `34-soliton-reflection.qmd` (`sec-soliton-reflection`) and `35-quantum-symmetric-pairs.qmd` (`sec-qsp`);
  - Lecture 6: chapters 30 and 34, and `36-particle-reflection.qmd` (`sec-particle-reflection`).
- Background: chapter 21 (`sec-qg-symmetry`, `sec-qg-coproduct`, `eq-q-physical`) and, for Hirota solitons, chapter 10 (`sec-hirota`).
- To find an object: `grep -rn '{#<id>}' part*/ appendices/ outlook/ *.qmd lecture-notes/half-line/`.

**Sources.** Texts already in `~/referee-copy/papers/` (or `papers/`, if you run through `run_referee_agy.sh`) are named `<arXiv id>.txt`, with `/` replaced by `_`, and the PDFs are next to them. Fetch any other source from arXiv. These are the arXiv ids of the sources the lectures cite:
- Ghoshal–Zamolodchikov, hep-th/9306002 (`gz1994`);
- Corrigan–Dorey–Rietdijk–Sasaki, hep-th/9404108 (`cdrs1994`);
- Corrigan–Dorey–Rietdijk, hep-th/9407148 (`cdr1995`);
- Bowcock–Corrigan–Dorey–Rietdijk, hep-th/9501098 (`bcdr1995`);
- Saleur–Skorik–Warner, hep-th/9408004 (`ssw1995`);
- de Vega–González-Ruiz, hep-th/9306089 (`devega1994`);
- Sasaki, hep-th/9311027 (`sasaki1993`);
- Fring–Köberle, hep-th/9304141 (`fring_koberle1994a`);
- Kim, hep-th/9506031 (`kim1995`);
- Bowcock, hep-th/9609233 (`bowcock1998`);
- Corrigan, hep-th/9707235 (`corrigan1998`);
- Mezincescu–Nepomechie, hep-th/9709078 (`mezincescu_nepomechie1998`);
- Dorey–Pocklington–Tateo–Watts, hep-th/9712197 (`dptw1998`);
- Gandenberger, hep-th/9806003 (`gandenberger1999a`);
- Delius, hep-th/9807189 (`delius1998a`) and hep-th/9809140 (`delius1998b`);
- Dorey–Tateo–Watts, hep-th/9810098 (`dtw1999`);
- Delius–Gandenberger, hep-th/9904002 (`delius_gandenberger1999`);
- Bajnok–Palla–Takács–Tóth, hep-th/0106070 (`bptt2002`);
- Delius–MacKay, hep-th/0112023 (`delius_mackay2003`);
- Baseilhac–Koizumi, hep-th/0507053 (`baseilhac_koizumi2005`);
- Baseilhac–Belliard, arXiv:0906.1215 (`baseilhac_belliard2010`);
- Kolb, arXiv:1207.6036 (`kolb2014`);
- Balagović–Kolb, arXiv:1507.06276 (`balagovic_kolb2019`);
- Appel–Vlaar, arXiv:2007.09218 (`appel_vlaar2022`).

Four are not on arXiv: Cherednik, Theor. Math. Phys. 61 (1984) 977 (`cherednik1984`); Sklyanin, J. Phys. A 21 (1988) 2375 (`sklyanin1988`); Bernard–LeClair, Commun. Math. Phys. 142 (1991) 99 (`bernard1991`); Letzter, J. Algebra 220 (1999) 729 (`letzter1999`). If you cannot read a source, say so and mark that citation "unclear"; do not reconstruct what it probably says.

## How to work

For each claim you check, do **one independent computation**: a derivation written out, or a script of your own. You may run the drafter's or the book's scripts for comparison, but they are not evidence. A verdict without its own evidence counts as "unclear". `briefs/referee_defects_38_40_report.md` shows the expected depth.

The lectures teach by calculation. A blackboard calculation claims that each displayed step follows from the one before, so check the steps, not only the result. An exercise claims that its statement is true and that a student can solve it with what the lectures provide.

## Phase 1: your own list

1. **Read** `index.qmd` and the six lectures in full before checking anything: text, captions and exercises.
2. **Write the risk list** as section 1 of the phase-1 report. It holds the claims you judge most likely to be wrong, overstated, misattributed or misleading for a student, ranked, with roughly 15–30 entries.
   - For each entry give its location (id, or file and line), the claim quoted, why it is at risk, and how you will check it.
   - Write this list before you start checking, and do not reorder or delete entries afterwards. Append claims you discover later at the end, marked "added during checking".
   - What counts as risky is for you to decide. Formulas, signs and normalizations, derivations, agreement with the book and with its status labels, attributions and whether citations say what their sources say, claims that a script checks something, figure captions and the exercises are all fair game.
3. **Check** the entries in order of rank, at least the top ten, and more if you can. For each, give the verdict (correct / wrong / unclear), your own computation and its output, and, for errors, replacement text ready to paste, with file and line.
4. **List** the entries you did not check, and why.
5. **Stop.** End your phase-1 session with a final message that says "Phase 1 complete" and gives the report path. Do not start phase 2 on your own; the editor will send it.

## Phase-1 report format (`briefs/referee_lectures_halfline_phase1_report.md`)

1. A header with the date, your model, and the commit you reviewed (`git rev-parse HEAD`).
2. Section 1: the risk list, as written before checking, with any later additions marked.
3. Section 2: one subsection per checked entry, in rank order.
4. Section 3: the entries you did not check, and why.
5. Section 4: other observations, labelled as unverified.

In phase 2 you will write a second report, `briefs/referee_lectures_halfline_report.md`. Your phase-1 report stays as it is.
