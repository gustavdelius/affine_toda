# Brief: residues of boundary Coleman–Thun diagrams and masked boundary bound states (half-line open problem 8)

You are a mathematical physicist with expertise in the analytic S-matrix, Landau singularities and the boundary bootstrap. The repository is `/home/gustav/Git/affine_toda`. Do not edit any file in it. Work in your own scratch directory and deliver your full report in your final message: subagents may be unable to write report files. Do NOT spawn subagents. This brief was written from the book, from abstracts and from machine summaries of papers. Treat every claim in it as unverified until you have rechecked it.

## The problem

Open problem 8 of `sec-boundary-open-problems` (`part7-half-line/35-particle-reflection.qmd`):

> **Masked bound states.** Delius and Gandenberger explain many poles by anomalous thresholds. In the bulk of non-simply-laced theories, Coleman–Thun poles can mask genuine bound states [DG sec. 7]. Deciding this on the half-line requires computing the residues of boundary Coleman–Thun diagrams, which has not been done.

Two concrete targets in the book:

1. **The unexplained pole** (`sec-boundary-spectrum-a2`). On $b_{-n,n}$ the pole $p=\frac12+(2n-1)\frac B2$, which DG attribute to creating $b_{-n+1,n}$, is in the physical strip for $B<\frac2{2n-1}$, but $b_{-n+1,n}$ is in the spectrum only for $B<\frac1{2n-1}$. At $B=\frac9{20}$ this is the pole $\frac{47}{40}$ of $K_1$ on $b_{-2,2}$ (`tbl-a2-boundary-states`). "Either a Coleman–Thun diagram explains it, or the spectrum is larger than DG's; this is open."
2. **Coincident diagrams.** At $B=\frac1{2k}$ levels cross, and two diagrams contribute to the same pole. For $a_4$ at $n_2=-n_3$, a pole gets two diagrams and fails the naive bootstrap. For $a_3$, $b=2$, the analysis is incomplete: DG (6.11)–(6.14) need states with negative labels whose amplitudes were never computed.

## What the book already has

Read `sec-boundary-bound-states`, `sec-boundary-coleman-thun` (ch. 28, including the Mattsson–Dorey lemmas); `sec-coleman-thun` and `sec-dorey-formula` (bulk); `sec-boundary-spectrum-a2` (ch. 35). Scripts: `textbook_code/particle_reflection.py` (functions `K1nm`, `Kex`, `bbs`, `K1_a3`, `F`, `E`: excited-state amplitudes, bootstrap steps, energies), `boundary_axioms.py`, and the bulk pole tools in `foundations_code/bulk_poles.py`. The counting rule used so far: a diagram with $P$ on-shell propagators and $L$ loops gives a pole of order $P-2L$, lowered by one for every reflection at a zero (Dorey–Tateo–Watts).

## Literature

- Dorey–Tateo–Watts (`dtw1999`), "Generalisations of the Coleman–Thun mechanism and boundary reflection factors": the boundary counting rule, boundary Lee–Yang.
- Mattsson–Dorey (`mattsson_dorey2000`) and Bajnok–Palla–Takács–Tóth (`bptt2002`): Dirichlet and general boundary sine-Gordon, where every pole is accounted for.
- Bajnok–Böhm–Takács, hep-th/0309119 (NPB 682 (2004) 585) and "Boundary reduction formula", hep-th/0207079 (not in `references.bib`): Landau equations, Cutkosky rules and the Coleman–Norton interpretation **with a boundary**. This is the framework in which boundary Coleman–Thun residues can be defined.
- Dorey–Polvara, "From tree- to loop-simplicity in affine Toda theories I", arXiv:2206.09368 (JHEP 09 (2022) 220; not in `references.bib`): the Laurent coefficients at $2N$th-order bulk poles of simply-laced ATFT from one loop, by cutting loop diagrams into products of tree amplitudes. This is the bulk method to transplant.
- Homrich et al., "Nonperturbative anomalous thresholds", arXiv:2212.06157 (not in `references.bib`): explicit non-perturbative formulas for anomalous-threshold singularities, tested against the Coleman–Thun poles of the $E_8$ model. On-shell diagrams give residues as products of on-shell couplings with a geometric factor.
- Riva, hep-th/0102163 (NPB 604 (2001) 511) and thesis hep-th/0106268 (not in `references.bib`): boundary bound-state structure of the $E_n$ ATFTs, with "a more complicated relation between the spectra of bound states" than in the bulk. This is the closest existing pole analysis beyond $a_n$.

## Suggested plan of attack

1. **Residue formula.** Write the residue of a boundary on-shell diagram as a product of bulk three-point couplings $\Gamma_{ab}^c$ (from residues of $S$), boundary couplings $g^\beta_{a\alpha}$ (from residues of $K$ at bound-state poles, `sec-boundary-bound-states`), reflection factors evaluated on the internal lines, and a Jacobian from the diagram's geometry. Derive it by the boundary Landau analysis of Bajnok–Böhm–Takács, or postulate it by analogy with Homrich et al. and calibrate it on boundary Lee–Yang (DTW) and on the sine-Gordon Dirichlet spectrum (Mattsson–Dorey).
2. **Reflection at a zero.** Make "lowered by one at a zero" quantitative: the residue picks up $K'_a$ at the zero. Check the formula on a case where DTW's rule is known to work.
3. **The $a_2$ test.** For the pole $p=\frac12+(2n-1)\frac B2$ on $b_{-n,n}$ in the window $\frac1{2n-1}<B<\frac2{2n-1}$, enumerate the candidate diagrams with the existing states and fusings. Compute their total residue and compare it with the residue of the amplitude (`eq-a2-excited`). A match explains the pole; a mismatch means a missing state, whose energy and amplitude then follow from the bootstrap.
4. **Coincident diagrams.** At $B=\frac1{2k}$, and for the $a_4$ case $n_2=-n_3$, check that the sum of the two diagram residues equals the amplitude's residue (or its leading Laurent coefficient for double poles).
5. **Masking (the stated problem).** Pick a case where a pole is explained by a Coleman–Thun diagram but the residue sum does **not** match the amplitude: the remainder must come from a bound state. Natural places to look: the $a_3$, $b=2$ amplitudes DG (6.11)–(6.14) (with the corrected (6.11), see the notes of ch. 35), and Riva's $E_n$ cases. Report every masked state found.

## Known traps

- The physical strip on the half-line is $0<\operatorname{Im}\theta<\frac\pi2$; poles at $\frac{i\pi p}h$ with $p$ in book blocks.
- Signs of residues: for boundary bound states, the sign of $-i\operatorname{Res}K$ relates to the positivity of $(g^\beta_{a\alpha})^2$. In a non-unitary context (imaginary coupling) signs can be either way; at real coupling they should be fixed.
- The boundary energies (`eq-a2-energies`) use the renormalized $m_1$.
- DG (6.11) line (d) has a sign misprint (book's notes), and DG (6.26) lacks a block.

## Deliverable

- The residue formula, its derivation or calibration, and its tests.
- The verdict on the $a_2$ unexplained pole, with numbers.
- Any masked boundary bound states found, with their energies.
- Proposed replacement text for the paragraph "This range leaves one pole unexplained" and for open problem 8.
- Scripts.

## Feasibility (assessment when the brief was written, 2026-10-11)

**Medium.** The diagram enumeration and the numerical residues are concrete, and the bulk analogue is now well developed (Dorey–Polvara, Homrich et al.). The risky part is a trustworthy boundary residue formula including reflections at zeros; calibration on boundary Lee–Yang and sine-Gordon makes it testable.
