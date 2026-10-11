# Brief: completeness of the boundary coideals and their classical boundaries (half-line open problem 5)

You are a mathematician or mathematical physicist with expertise in quantum symmetric pairs, the reflection equation and classical integrable boundaries. The repository is `/home/gustav/Git/affine_toda`. Do not edit any file in it. Work in your own scratch directory and deliver your full report in your final message: subagents may be unable to write report files. Do NOT spawn subagents. This brief was written from the book, from abstracts and from machine summaries of papers. Treat every claim in it as unverified until you have rechecked it.

## The problem

Open problem 5 of `sec-boundary-open-problems` (`part7-half-line/35-particle-reflection.qmd`):

> **Completeness of the coideals.** Do the coideals of `sec-boundary-qg` and their generalizations (`sec-generalized-satake`) exhaust the integrable boundaries? Which of them give soliton-preserving reflection matrices, and what are the corresponding classical boundaries beyond the one of `sec-soliton-preserving`? The $s=-1$ solution is a coideal without a known classical boundary condition. Baseilhac and Belliard propose a classification of integrable fixed and dynamical boundaries through generalized q-Onsager algebras, which is a partial answer (`sec-qsp`).

Restrict to $a_n^{(1)}$ and the vector solitons first, where the reflection equation has been solved completely in the literature.

## What the book already has

Read all of ch. 34 (`sec-qsp`), especially `tbl-physics-qsp`, `prp-trig-k`, `prp-hypothesis-ii` and `sec-generalized-satake`; and `sec-sklyanin`, `sec-durham-classification`, `sec-soliton-preserving`, `sec-background-field` (ch. 30). Script: `textbook_code/qsp_checks.py`. Key facts:

- Classical side: with $\sigma$ = transposition and constant $K$, the boundaries are those of `tbl-durham`; with $\sigma$ = Hermitian conjugation and $K=1$, the soliton-preserving one. Other $\sigma$ (diagram automorphisms composed with transposition) are not studied in the book. Background-field boundaries need an inner automorphism of odd order.
- Quantum side: three soliton-conjugating points for $a_n$ ($\hat\epsilon=0$, $s=\pm1$), torus conjugates, and the generalized Satake families of Regelskis–Vlaar.

## Literature

- **Regelskis–Vlaar, arXiv:1803.06491 (J. Phys. A 51 (2018); not in `references.bib`):** the complete set of invertible solutions of the untwisted and twisted reflection equations for the $A^{(1)}_{N-1}$ Bazhanov–Jimbo R-matrix, all obtained by affinization from constant solutions. For the vector representation this settles the quantum side of "do the coideals exhaust the solutions".
- Regelskis–Vlaar (`regelskis_vlaar2016`): generalized Satake diagrams and the K-matrices of their quantum pair algebras (Thm. 9.3, a conjecture for general $n$; App. D).
- Appel–Vlaar (`appel_vlaar2022`, `appel_vlaar2025`) and Appel–Vlaar, "Tensor K-matrices for quantum symmetric pairs", arXiv:2402.16676 (not in `references.bib`): universal and trigonometric K-matrices for QSPs of affine type.
- Baseilhac–Belliard (`baseilhac_belliard2010`): the generalized q-Onsager algebras and a proposed classification of fixed and dynamical boundaries.
- Avan–Doikou, arXiv:0809.2734 (not in `references.bib`): boundary Lax pairs for $A_n^{(1)}$ Toda, **two distinct boundary types already for $A_2^{(1)}$**; a candidate classical home for coideals not covered by `tbl-durham`.
- Arnaudon–Avan–Crampé–Doikou–Frappat–Ragoucy, math-ph/0406021 (not in `references.bib`): classification of SP and SNP solutions for $sl(N)$ spin chains.
- Lima-Santos and Malara, nlin/0412058 (not in `references.bib`): regular solutions of the reflection equation for non-exceptional vertex models ($b$, $c$, $d$ and twisted). This is the starting point beyond $a_n$.
- W. Wang's survey of ıquantum groups, arXiv:2112.10911, for the current state of affine QSP theory.

## Suggested plan of attack

1. **Quantum census for $a_n$, vector representation.** From Regelskis–Vlaar 2018, list all invertible solutions (twisted = soliton-conjugating, untwisted = soliton-preserving) up to torus gauge and up to the equivalences that do not change physics. For each, identify the coideal: a Kolb QSP, a generalized Satake family, or a generalized q-Onsager point such as $s=\pm1$. Check a sample numerically in book conventions with `qsp_checks.py`. Answer: does every solution come from a coideal of the type in `tbl-physics-qsp`?
2. **Physical filter.** Keep the solutions that satisfy boundary unitarity and crossing-unitarity with some scalar factor (the matrix part of crossing must factor). These are the candidates for integrable boundaries of the field theory.
3. **Classical census.** Generalize Sklyanin's construction in `sec-sklyanin` to $\sigma=\gamma\circ{\rm transposition}$ and $\sigma=\gamma\circ\dagger$ with $\gamma$ a diagram automorphism of $a_n^{(1)}$ (rotations and reflections of the cycle). Solve the analogue of (`eq-boundary-k-equation`) for constant $K$ and find the boundary conditions. Compare with Avan–Doikou's two $A_2$ types. Every new classical boundary is a candidate partner for an unassigned quantum solution, the $s=-1$ point first.
4. **Matching.** Use the $\Theta$ criterion (`sec-boundary-theta`), the soliton-preserving/conjugating type, the number of free parameters, and the first-order coideal parameters (`prp-dm-first-order`) to pair classical and quantum boundaries. List the unpaired ones on each side.
5. **Beyond $a_n$ (stretch).** Repeat step 1 for $d_n$ with Lima-Santos–Malara and Regelskis–Vlaar 2016's D-families. Compare with `prp-dn-k` and `tbl-durham`.

## Known traps

- Kolb's parameter $c^{\rm Kolb}=-q^{-1}$ at the physical point, not $-q$; the $s=\pm1$ coideals are not Kolb QSPs (their $s\notin\mathcal S$), though Kolb's Serre obstruction vanishes at $\hat\epsilon^2=\epsilon_*^2$ (`sec-letzter-kolb`).
- Twisted versus untwisted reflection equation (`tbl-k-types`) and the homogeneous/principal gradations.
- Torus dressing changes $c_j$ but not physics (`exr-torus-dressing`).

## Deliverable

- The quantum census table for $a_n$ vector with coideal identifications, and the physical filter.
- The classical census for the automorphism-twisted $\sigma$, with the new boundary conditions if any.
- The proposed pairing, and the remaining unpaired items; in particular, a candidate classical boundary for $s=-1$ or an argument that none exists among constant-$K$ boundaries.
- Proposed text for open problem 5 and an update of `tbl-physics-qsp`.

## Feasibility (assessment when the brief was written, 2026-10-11)

**Medium** for $a_n$ vector: Regelskis–Vlaar 2018 already did the hard quantum classification, so steps 1–2 are translation and checking, and step 3 is a finite computation. **Low** for full completeness ("exhaust the integrable boundaries") in general representations and algebras, which is a research programme.
