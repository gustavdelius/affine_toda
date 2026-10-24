# Brief: referee pass on chapter 16 (affine Toda theories based on Lie superalgebras)

You are an adversarial mathematical-physics referee. You are working in the repository of a Quarto book on affine Toda field theory (worktree `.claude/worktrees/superalgebra-toda`, branch `claude/superalgebra-toda`).

**Rules.**
- Other agents may be editing files in this worktree at the same time.
- Do not edit any existing file.
- Do not run git commands that change the index or the working tree (no commit, checkout, stash, reset, rebase or restore).
- Do not run `quarto render`.
- Do not spawn subagents.
- Write scratch scripts only under `referee_scratch/superalgebra_toda/`.
- Write the report to `briefs/referee_superalgebra_toda_report.md`.

## Context

Chapter 16 (`part4-real-coupling/16-superalgebra-toda.qmd`, `sec-super-toda`) is new. It was drafted on 2026-10-23 by a Claude agent from the handover `briefs/draft_superalgebra_toda.md`. Reviewed commit: the commit "Add chapter 16, affine Toda theories based on Lie superalgebras" on branch `claude/superalgebra-toda`, which also adds this brief; find its hash with `git log --format=%h -1 --grep="Add chapter 16, affine Toda"` and give it in your report header. The drafting report is `briefs/draft_superalgebra_toda_report.md`; read its sections on the verification record and on limitations first.

Also under review are the integration edits of the same commit: `part2-algebra/06-affine-algebras.qmd` (paragraph "Superalgebras" in `sec-lie-algebraic-data`), the Notes of chapters 12 and 14, `sec-ode-im` (chapter 15), the Notes of chapter 23 (`sec-known-smatrices-notes`, paragraph "Supersymmetric scattering"), the new sections `sec-bsg-super` (chapter 30) and `sec-sgd-super` (chapter 38), and three new entries of Appendix D (`appendices/D-corrections.qmd`: Delius–Grisaru–Penati–Zanon 1991, Gandenberger 1998, Evans–Madsen 1999). Also under review are the two sentences pointing to the new chapter at the ends of the relevant paragraphs in `sec-lagrangians` and `sec-lax-pair` of chapter 9, added in the next commit on the branch, "Chapter 9: point to the superalgebra chapter". The renumbering of chapters 16–50 to 17–51 in the same change is mechanical and not under review.

Check scripts (all in `code/super_toda/`, run from the repository root; each prints PASS/FAIL lines and exits nonzero on a failure; see `code/super_toda/README.md`):

| Script | Checks | Runtime |
|---|---|---|
| `superspace.py` | 9 | seconds |
| `components.py` | 31 | about 25 s |
| `lax_currents.py` | 11 | about 10 s |
| `dgpz_smatrix.py` | 27 | about 1 min |
| `n2_solitons.py` | 7 | about 80 s |
| `gandenberger_scalar.py` | 2 | seconds |

The drafter flags for the referee: the sign $\gamma=-\sin\pi\Omega$ (fixed by the drafter's own tree amplitude, not by the source); the sign of $\lambda$ in the first-order soliton equation of Evans–Hollowood; the identification of the super sinh-Gordon algebra as $A(1,0)^{(2)}$ rather than $\mathfrak{osp}(1|2)^{(1)}$; and the three Appendix D entries.

## Notation

- Light cone $x^\pm=\frac12(t\pm x)$, $\partial_\pm=\partial_t\pm\partial_x$ (the book's $\partial_\pm$).
- Superspace: $D_\pm=\partial/\partial\theta^\pm-i\theta^\pm\partial_\pm$, $D_\pm^2=-i\partial_\pm$; $\Phi=\phi+i\theta^+\psi_++i\theta^-\psi_-+i\theta^+\theta^-\sigma$; $\int d^2\theta\,\theta^+\theta^-=i$. These are the conventions of Evans–Hollowood 1992 and Evans–Madsen 1997; Evans–Hollowood 1991, Delius–Grisaru–Penati–Zanon and Ito–Zhu use others.
- The coupling enters as $\Phi\to\beta\Phi$, $\mathcal L\to\mathcal L/\beta^2$. The book's roots have longest squared length 2 where the bosonic sector is one Toda theory; the super sinh-Gordon odd root has $\vec\alpha^2=\frac12$. Delius–Grisaru–Penati–Zanon normalize even simple roots to squared length 2, so $\beta_{\rm book}^2=2\beta_{\rm DGPZ}^2$; their $\{x\}$ is the book's block (`eq-blocks`) with $B=2\Omega$.
- Scattering operator $\check S$ in braid form: in-basis $a\otimes b$ mapped to out-states in spatial order; at zero coupling it is the graded permutation. Tensor products of even operators are plain Kronecker products.

Read first: `briefs/superalgebra_toda_conventions.md`; `sec-lagrangians`, `sec-lax-pair`, `sec-blocks`, `sec-non-simply-laced`, `sec-classical-solitons` and `eq-soliton-mass` of the book. Source texts in `~/referee-copy/papers/*.txt` if present, otherwise on arXiv: hep-th/9703065 (Evans–Madsen 1997, v2), hep-th/9712227 (Evans–Madsen 1999), hep-th/9207073 (Evans–Hollowood 1992), hep-th/9205045 (Penati–Zanon), hep-th/9608190 (Evans–Madsen 1996), hep-th/9703158 (Gandenberger), 2206.08024 (Ito–Zhu), hep-th/0103029 (Nepomechie), hep-th/0601014 and 1412.2579 (super sinh-Gordon defects), hep-th/9201067 (DGZ 1992). The two 1991 papers of Delius, Grisaru, Penati and Zanon (Nucl. Phys. B359 125, Phys. Lett. B256 164) and Evans–Hollowood 1991 (Nucl. Phys. B352 723) are not on arXiv; ask the editor for copies if they are not in `~/referee-copy/papers/`.

To find an object: `grep -rn '{#<id>}' part*/ appendices/ outlook/ *.qmd`.

For each item, do **one independent computation**: a derivation written out, or a script of your own. You may run the book's scripts for comparison, but they are not evidence. A verdict without its own evidence counts as "unclear". `briefs/referee_defects_38_40_report.md` shows the expected depth.

## Items (most important first)

1. **Component Lagrangian and signs (`eq-super-action`, `eq-super-eom`, `eq-super-components`, `eq-ssg-components`).**
   - Derive the component Lagrangian from the superspace action in the stated conventions, with your own Grassmann bookkeeping. Check the signs of the fermion kinetic term, the Yukawa term and the auxiliary-field term, and the $\theta^+\theta^-$ term for even roots.
   - Check that the stated field equation is the Euler–Lagrange equation of the action, in components.
   - For super sinh-Gordon, check the masses ($2m$ for both), the bosonic truncation as the book's $a_1^{(1)}$ theory with the same $m$ and $\beta$, and the continuation $\beta\to i\beta$.
   - Check `exr-super-masses` ($M_B^2=M_F^2$ for purely odd systems).
2. **The graded Lax pair (`eq-super-lax`, `eq-super-zero-curvature`, `eq-ssg-realization`).**
   - Verify the order-by-order argument in $\lambda$, including the sign of the curvature and $\eta_j=-i\mu_j$, with an explicit Grassmann-envelope convention of your own (state it).
   - Check the $\mathfrak{sl}(2|1)$ realization of the four relations, and the claim that a realization in the loop algebra of $\mathfrak{osp}(1|2)$ with $e_0,f_1\propto F_-$ fails.
   - Is the identification of the super sinh-Gordon algebra as $A(1,0)^{(2)}\simeq C(2)^{(2)}=\mathfrak{osp}(2|2)^{(2)}$, and not $\mathfrak{osp}(1|2)^{(1)}$, correct and correctly attributed (Evans–Madsen 1997 Table 4, Ito–Zhu Sec. 3)? Is the sentence about $B(0,1)^{(1)}$ right?
3. **Supercurrents and the Evans–Madsen correction (`eq-ssg-currents`, Appendix D).**
   - Verify $D_-J=D_+\bar J$ on shell for both currents with the stated $\bar J$, by your own on-shell reduction.
   - Verify from Evans–Madsen 1999 (5.47) that their component equations (5.48)–(5.51) follow from $D_-J=D_+\bar J$ and not from their printed (5.44).
   - Check the claimed normalization of $\bar J$ against their Lagrangian (1.5). Is this a misprint, as Appendix D says, or a convention we have misread? Quote the source lines.
   - Check the statement about the spin-4 density and `eq-spin3-current`.
4. **The boson–fermion S-matrix (`eq-dgpz-operator` to `eq-dgpz-san`, `cnj-dgpz-smatrix`).**
   - From DGPZ (2.4) and App. A, check that `eq-dgpz-operator` represents their amplitudes correctly, including the sign of the $\psi\psi\to\phi\phi$ entry and the claim that the free limit is the graded permutation.
   - Check that `eq-dgpz-ybe` is equivalent to their (2.7) and partners, and that it fails for the ungraded assignment.
   - Evaluate $T=\Pi(\theta)\Pi(i\pi-\theta)$ yourself and check crossing, `eq-dgpz-unitarity`, unitarity of $\check S$, and that $T$ has no poles in the physical strip.
   - Check crossing $D(i\pi-\theta)=iR(\theta)$ and its attribution to the Majorana crossing phase. Is the phase derived from the physical charge-conjugation matrix, or only consistent with it?
   - **Sign of $\gamma$.** Compute the tree amplitude of $\phi_n\phi_n\to\phi_n\phi_n$ from the Lagrangian yourself (cubic and quartic couplings in the mass basis, contact and three exchange channels, `eq-s-from-m`) and decide whether $\gamma=-\sin\pi\Omega$.
   - Check the bound-state poles, the residue signs and the rank-one claim, and the bootstrap amplitude $S_{an}$.
   - Is `cnj-dgpz-smatrix` labelled and scoped correctly, and is its evidence list faithful to DGPZ Sec. 4 and the letter?
5. **The free-fermion condition (Appendix D, DGPZ (B.1)).** Decide whether (B.1) as printed contradicts DGPZ's own (2.4) and (2.34), or whether a convention for $D$ in their Appendix B makes it right. Quote the lines.
6. **$C^{(2)}(n+1)$ and the masses (`eq-dgpz-masses`, `eq-dgpz-cn`).**
   - Reconstruct the simple roots of $A^{(2)}(0,2n-1)$ and $C^{(2)}(n+1)$ from Evans–Madsen's diagram rules (Table 2(b)) and check the masses for $n=2,3$.
   - Check unitarity and crossing of the fermion amplitudes for $n$ even and odd, and the statements about the folded $c_n^{(1)}$ S-matrix and the fermion box diagram.
   - Is "Guided by the tree amplitudes, Delius et al. assume ..." a faithful account?
7. **The $N=2$ family (`eq-n2-condition`, `eq-n2-masses`, `eq-n2-central`, `cnj-n2-sine-gordon`).**
   - Check that `eq-n2-condition` forces null roots and a bipartite graph, and the claims for $A(1,1)^{(1)}$ and $D(2,1;\alpha)^{(1)}$.
   - Check the realization of the $\mathfrak{sl}(m|m)^{(1)}$ roots, and that at zero fermions the potential is $V_{\rm real}+V_{\rm imag}$ with the book's $m=\mu$.
   - Check that the one-soliton solves the field equation, its mass, and the central-charge formula on the asymptotic vacua.
   - Derive the first-order equation from Evans–Hollowood's (5.5) yourself, and decide the sign of $\lambda$.
   - Is the notation paragraph ($\mathfrak{sl}$ vs $\mathfrak{psl}$ vs $\mathfrak{gl}$) correct?
   - Is the scope of the BPS/central-charge statements (classical only) stated accurately?
8. **The $D(2,1;\alpha)^{(1)}$ family.** Check the root data, masses, kinetic signs under the three twisted reality conditions, the Yukawa coupling and its normalization, and the limits (foldings to $C(2)^{(2)}$ and $B(1,1)^{(1)}$, $\alpha\to0$ to $A(1,1)^{(1)}$) against Evans–Madsen 1999.
9. **Scope of quantum-integrability claims (`sec-super-quantum-currents`, `sec-super-fixed-masses`).** Check against Penati–Zanon and Evans–Madsen 1999 Sec. 6 what exactly is proved: to which order, in which scheme, for which reality conditions. Check the one-loop mass claims and their attribution (Grisaru–Penati–Zanon, PLB 253 (1991) 357; note that DGPZ's NPB reference list gives PLB 234 (1990) 88 for the same work).
10. **Gandenberger (`sec-known-smatrices-notes`, Appendix D).** Check the sign in his (5.10) by your own computation from (3.28) and (5.6), and his (5.17). Is the paragraph's statement of what the identities do and do not imply accurate?
11. **Boundary and defect sections (`sec-bsg-super`, `sec-sgd-super`).** Check each sentence against Nepomechie 2001, Gomes–Ymai–Zimerman 2006 and Aguirre–Gomes–Spano–Zimerman: the role of the boundary fermion, the two parameters, what was conserved, the status of the boundary S-matrix; the defect fields, the type-II character, the Lax and defect-matrix statements. Are classical and quantum results kept apart?
12. **Fidelity of citations and prior art.** For every citation in chapter 16 and the integration edits, check that the source says what the text attributes to it (quote the line), including `prp-super-bosonic`, `prp-super-susy` (its list of algebras), the claims attributed to Evans–Madsen 1996 and to Ito–Zhu, and the Notes. The chapter claims no original results and adds no item to Appendix E; check that this is right, and that nothing in the chapter is presented as new that is in the literature.
13. **Bibliography.** Spot-check the 23 new entries in `references.bib` against INSPIRE or the publishers (journal, volume, page, year, title, DOI); the drafter checked all but `vandeleur1989` against INSPIRE.

## Report format

- A header with the date, your model and the commit reviewed.
- One section per item: the verdict (correct / wrong / unclear); your own computation (script inline or path, and its output); and, for errors, replacement text ready to paste, with file and line.
- End with other observations, labelled as unverified, and a list of the claims you could not check.
