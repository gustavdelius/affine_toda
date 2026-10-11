# Brief: referee chapter 35 (breather and particle reflection amplitudes; open problems of Part VII)

You are a mathematical-physics referee. The repository is `/home/gustav/Git/affine_toda` (branch `claude/affine-toda-boundary-dae69e`, PR #6), a Quarto book on affine Toda field theory. Do not edit any file in it. Work in your own scratch directory, and deliver your full report in your final message, because subagents may be unable to write report files.

## Context

Part VII, "Affine Toda theory on the half-line" (`part7-half-line/`, chapters 28–35), was drafted on 2026-10-10. This brief covers chapter 35 (`35-particle-reflection.qmd`, `sec-particle-reflection`). Its sections:
- the boundary bootstrap from the soliton $K$ to the lowest breather (`prp-breather-amplitude`);
- the continuation in the coupling from breathers to real-coupling particles (`cnj-reflection-continuation`, with Delius–Gandenberger's three assumptions);
- the $a_n$ particle reflection amplitudes for the uniform $(+\dots+)$ boundary (`prp-uniform-amplitude`, `eq-uniform-amplitude`), for Neumann (`cnj-neumann-amplitude`, `eq-neumann-amplitude`) and for the "solitonic" boundary conditions (`prp-solitonic-amplitude`, `eq-solitonic-factor`);
- the $a_2$ boundary spectrum (`prp-a2-excited`, `eq-a2-excited`, `eq-a2-energies`, `tbl-a2-boundary-states`);
- the numbered open problems of the whole part (`sec-boundary-open-problems`).

Conventions:
- Blocks $(x)$ and $\{x\}$ are those of `eq-blocks`, with $h=n+1$ and $B$ from `eq-b-of-beta`.
- Diagonal crossing-unitarity is `eq-boundary-crossing-diagonal` (ch. 28); the bulk boundary bootstrap is `eq-boundary-bootstrap-bulk`.
- The identification of the soliton $K$ belonging to each boundary is in `sec-epsilon-constraint` (ch. 32). The diagonal $K$ belongs to $(+\dots+)$, and Neumann is the non-diagonal $s=+1$ matrix.

Find any object with `grep -rn '{#<id>}' part*/`.

Sources (fetch with `curl -sL https://arxiv.org/pdf/<id>`):
- Delius–Gandenberger hep-th/9904002 (DG), the main source;
- Gandenberger hep-th/9806003 (G98);
- Kim hep-th/9506031;
- Perkins–Bowcock hep-th/9807146;
- Corrigan–Dorey–Rietdijk–Sasaki hep-th/9404108, Corrigan–Dorey–Rietdijk hep-th/9407148;
- Bowcock hep-th/9609233;
- Dorey–Tateo–Watts hep-th/9810098;
- Corrigan hep-th/9707235;
- Ghoshal hep-th/9310188, Corrigan–Delius hep-th/9909145, Chenaghlou–Corrigan hep-th/0002065.

Experience in this repository is that every adversarial referee pass finds real errors.

## Already verified (do not redo unless you find a reason)

**By script** (`textbook_code/particle_reflection.py`, 43 checks, all PASS):
- **Bulk input.** Block identities; $a_n$ S-matrix unitarity, crossing and C-invariance.
- **Uniform boundary, DG (3.16), for $n=1..7$:**
  - it follows from DG (3.12) by iterating (3.15);
  - unitarity, and crossing-unitarity in the book form and in DG (3.17);
  - $K_a=K_{h-a}$, and the bulk bootstrap at every fusing;
  - the $B\to0$ limit is $-1/((a)(h-a))$ = `eq-linear-reflection` at $C=+1$;
  - it equals the CDRS conjecture (5.16) with $\langle x\rangle\to\langle\tilde x\rangle$ (control included);
  - its physical-strip poles are simple and coupling-independent (DG (6.1)).
- **Continuation.** The closed form (3.11) continues to (3.12).
- **$n=1$:**
  - DG (3.11) is Ghoshal's $K_{B_1}$ at $\eta=\pi/2$, $\vartheta=0$ (the $\hat\epsilon=0$ point);
  - (3.16) is `eq-shg-reflection` at $E=B/2$, $F=0$;
  - (3.20) is the Neumann point $E=1-B/2$;
  - Chenaghlou–Corrigan's $E=0$ differs from $E=B/2$ at $O(B^2)$.
- **Neumann, DG (3.20):**
  - it satisfies the axioms and equals (3.16) with $B\to2-B$;
  - its $B\to0$ limit is 1;
  - neither amplitude is self-dual;
  - its $O(\beta^2)$ term equals Kim's one-loop $a_3$ result (his (33), (39)), with a control, and G98 (6.11) for $a_2$.
- **G98 $K^{(-)}$ ($s=-1$).** It satisfies the axioms, and its dual is the CDRS minimal $C=-1$ conjecture (5.7).
- **Solitonic amplitudes, DG (5.5), for $n=2..6$, all $b$:**
  - the axioms hold, and $M_a^{(b)}=M_a^{(h-b)}$;
  - the $B\to0$ limit is (5.7);
  - DG (6.3), (6.18), (6.19), (6.25) are special cases, and so is (6.26) once a missing block is restored. As printed, (6.26) violates crossing-unitarity.
- **CDD factors.** DG (3.25) particle CDD factors are consistent.
- **$a_2$ excited states, DG (6.7), for $-3\le n,m\le3$:**
  - the axioms hold;
  - the boundary bootstrap at the three kinds of creation pole;
  - (6.6), (6.9), (6.10) are special cases.
- **DG (6.11) for $a_3$, $b=2$.** With line (d) corrected to $(-3+(2n_1+1)B/2)$ it reduces correctly; as printed it fails.
- **Energies.** The closed form $e_{n,m}-e_{0,0}=m_1[F(n)+F(m)]$ (`eq-a2-energies`) agrees with every bootstrap step.
- **Poles.** At $B=9/20$ all creation poles between the listed states are in the physical strip.

**By hand:**
- The statement of DG's assumptions, and the graded evidence for `cnj-reflection-continuation`.
- The rotor interpretation of the $b_{n,-n}$ tower as the remnant of the classical vacuum modulus. This is the drafting agent's own interpretation and is labelled as such.

## To be checked

Report on each item: correct / wrong (with the correction) / unclear (with what would settle it). Items are ordered roughly by risk.

1. **Corrections to Delius–Gandenberger.** Confirm or refute each against the journal version (Nucl. Phys. B 554 (1999) 325), not only the arXiv text:
   - (6.26): a block is missing, and the label should be $K_3^{(2)}$;
   - (6.11) line (d): sign;
   - (3.11): it is not the product of the three scalar factors alone, because the matrix part of the soliton S-matrix on the reflected singlet contributes too.
2. **The breather bootstrap for $n\ge2$** (`prp-breather-amplitude`). The singlet projection at matrix level, DG §3.4 / G98 appendix, is quoted, not redone. Redo it for $n=2$ with the diagonal $K$ and the book's $\check R$, or say what fails.
3. **Neumann.**
   - Is `cnj-neumann-amplitude` correctly labelled?
   - Does the text say clearly that the Neumann soliton $K$ is the non-diagonal $s=+1$ matrix (ch. 32–33), and that the soliton-route derivation of the Neumann particle amplitude exists only for $a_2$ (G98)?
   - Is Kim's one-loop comparison quoted correctly? Kim's results are for $a_3$, in his appendix.
4. **The $n=1$ discrepancy** (open problem 12). DG (3.16) at $n=1$ gives $E=B/2$, while Chenaghlou–Corrigan extended to $C_0=C_1=1$ gives $E=0$, and the two differ at $O(B^2)$. Which is right, or is the extension illegitimate? Relate it to `cnj-shg-e` (ch. 29) and Corrigan–Taormina.
5. **Unchecked quoted results:**
   - Bowcock (4.21) against DG (5.7);
   - Perkins–Bowcock's $O(\beta^2)$ uniform result for $a_2$ (quoted via G98);
   - the $a_4$ excited amplitudes (6.20), (6.21), (6.28), (6.29), checked only at the vacuum;
   - DG's Coleman–Thun explanations of the non-bound-state poles, which are quoted. The script checks only pole positions and the bootstrap. Check at least the $a_2$ diagrams of DG §6.1, against Dorey–Tateo–Watts's boundary Coleman–Thun rules.
6. **`eq-a2-energies` and `tbl-a2-boundary-states`.**
   - Are the closed form and the table right?
   - Do the existence ranges of the states (DG §6.2) match the poles listed?
   - Is the "quantized rotor" interpretation of the $b_{n,-n}$ tower ($\approx m_1(\pi B/6)n^2$) plausible, or misleading?
7. **The open problems** (`sec-boundary-open-problems`). Is each one crisp, correct and actually open? Check later literature for each, especially:
   - the Neumann amplitudes for general $n$;
   - boundaries with $\prod C_j=-1$;
   - particle reflection amplitudes beyond $a_n$;
   - duality on the half-line;
   - the non-unitary half-line theory and its relation to the Dirichlet walls of Part VIII (`sec-collapse-thresholds`, `thm-window`). That item contains the drafting agent's observation that a Part VIII Dirichlet wall fixes Re φ, which is half of the soliton-preserving condition. Is this correctly framed as a question?
8. **Status labels.** Should `prp-uniform-amplitude` be [Numerical] or [Established]? How should a result labelled [Numerical] be made conditional on the [Conjecture] `cnj-reflection-continuation`?
9. **Exercises.** For every `exr-`, including the open-ended `exr-modulus-rotor`, check that it is well posed and solvable with the text, and that no hint is wrong.
10. **Bibliography.** Check `perkins_bowcock1999`, `corrigan1998` and `bcz2004`, and the use of `delius_gandenberger1999`, `gandenberger1999a`, `kim1995`, `cdrs1994`, `cdr1995`, `bowcock1998`, `dtw1999` and `gandenberger1998`. Check authors, titles, journal data and DOIs in `references.bib`, and that each citation supports its sentence.
11. **Prior art and later work.** Search the INSPIRE citations of DG. Are the (6.11) and (6.26) misprints, or the closed form of the $a_2$ boundary energies, already in the literature? Has the Neumann or general-$n$ problem been solved since?

## Deliverable

One entry per item: a verdict, the justification, and, for errors, the corrected text ready to paste, with file and line. Include any scripts you wrote, either inline or as a list of what they check with their output. Do not mark anything verified on the strength of a citation alone.
