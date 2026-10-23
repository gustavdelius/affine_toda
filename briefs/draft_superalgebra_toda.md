# Handover: affine Toda theories based on Lie superalgebras

Date: 2026-10-23. Prepared by Codex for a Claude drafting session.

Repository inspected: `/home/gustav/Git/affine_toda`, HEAD `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.
Paths and chapter numbers below describe that state; resolve existing material by section id if the tree has moved.

## Task and authorization

You are Claude, drafting for the Quarto book *Affine Toda field theory*, coordinated by Gustav Delius. **Do not spawn subagents.**

The editor asked whether supersymmetric affine Toda theories and Toda theories based on affine superalgebras deserve a chapter. After a literature survey, Codex recommended one focused chapter at the end of Part IV, with short connections in existing chapters. The editor approved implementation in the style of the existing chapters, then requested this handover so that Claude, rather than Codex, does the drafting.

Implement the new chapter, its verification scripts, bibliography and associated integration edits. This is a drafting assignment, not a referee assignment. Follow `AGENTS.md`; this brief does not override it. Do not commit, push, merge or publish without an explicit request. Work in an isolated worktree under `.claude/worktrees/`, check the current coordination records, and avoid files under active referee review. Record any deferred integration edits by path and id instead of silently omitting them.

The main checkout contained unrelated untracked files when inspected. Preserve all unrelated work. This handover is the only file written by Codex for this task.

## Editorial outcome

Write a chapter titled **Affine Toda theories based on Lie superalgebras**, placed after the current chapter 15, *Beyond the S-matrix*, at the end of Part IV. Suggested chapter id: `sec-super-toda`; suggested filename after renumbering: `part4-real-coupling/16-superalgebra-toda.qmd`.

The title deliberately includes theories without spacetime supersymmetry. Explain the distinction at the start. The chapter should give a coherent construction and worked examples, then a selective account of quantum results. It should not become a catalogue of superalgebras or an unexplained list of papers. Aim for roughly 4,000–6,000 words of connected exposition, allowing the worked derivations to determine the final length.

Part IV's present title is *The real-coupling quantum theory*. Keep the new chapter's core discussion on massive particle theories and explain the excursions into complex reality conditions and solitons. Update the part introduction to make this scope clear. Do not present all superalgebra theories as real-coupling, unitary or purely massive.

## Read first and match the style

- `AGENTS.md`, `README.md`, `01-introduction.qmd`, especially `sec-conventions`.
- `part2-algebra/06-affine-algebras.qmd`: the algebraic knowledge available to the reader.
- `part3-classical/09-lagrangian.qmd`: `sec-lagrangians`, `sec-lax-pair`, `sec-masses-couplings`.
- `part3-classical/11-reality.qmd`: the existing treatment of reality and antilinear symmetry.
- `part4-real-coupling/12-perturbation-theory.qmd`, `13-simply-laced.qmd`, `14-non-simply-laced.qmd`, `15-beyond-smatrix.qmd`.
- `part6-soliton-smatrices/21-construction.qmd`, `22-known-smatrices.qmd`, `27-rsos.qmd` for later connections.
- `part7-half-line/29-boundary-sine-gordon.qmd` and `part8-defects/37-defect-sine-gordon.qmd` for the short extension sections.
- `briefs/figures_for_the_book.md` if a figure would help.

Use chapters 9 and 13 as models for the exposition: a short learning-goals callout; a concrete physical question; definitions immediately before use; numbered displayed formulas; an explicit low-rank example; derivations; exercises; Notes and references. Use short declarative sentences, British spelling and the established notation. Explain superspace and fermion parity rather than assuming a prior course in supersymmetry. Do not reproduce the more rhetorical style of the preface.

Use the book's statement divs and status tags. Published results should normally be labelled **[Established]** with precise attribution. A conjectured exact S-matrix remains conjectural even if its algebraic identities are checked. Numerical evidence and a finite-order perturbative check do not establish an all-orders result. Do not introduce original research claims merely to give the chapter theorem boxes.

## Source reading: leads, not verified formulas

The survey established a substantial literature, but it was not a formula-by-formula referee pass. Some evidence came from abstracts and bibliographic records. Independently read the primary papers, including errata, before importing their formulas or scope claims. Build a source-to-book dictionary with page/equation locators. Do not treat this brief or the earlier chat as evidence.

1. **J. Evans and T. Hollowood, *Supersymmetric Toda field theories*, Nucl. Phys. B352 (1991) 723–768.** DOI: https://doi.org/10.1016/0550-3213(91)90105-7 . Foundational source for the superalgebra zero-curvature construction, supersymmetry and integrable deformations. Evans–Madsen's bibliography lists an erratum in Nucl. Phys. B382 (1992) 662: obtain and account for it before using the original formulas.
2. **G. W. Delius, M. T. Grisaru, S. Penati and D. Zanon, *Exact S-matrix and perturbative calculations in affine Toda theories based on Lie superalgebras*, Nucl. Phys. B359 (1991) 125–167.** Bibliographic record: https://pure.york.ac.uk/portal/en/publications/exact-s-matrix-and-perturbative-calculations-in-affine-toda-theor/ . Read the full paper to select the models and scattering formulas. Also consult the companion *The exact S-matrices of affine Toda theories based on Lie superalgebras*, Phys. Lett. B256 (1991) 164–172, https://doi.org/10.1016/0370-2693(91)90668-G . These are central sources, not optional historical footnotes. Do not confuse them with `dgz1992` or `dgz1996`.
3. **J. Evans and T. Hollowood, *Integrable N=2 Supersymmetric Field Theories* (1992).** https://arxiv.org/abs/hep-th/9207073 ; Phys. Lett. B293 (1992) 100–110. Source for massive families based on the authors' `sl(m,m)^(1)` notation, supersymmetric solitons and conjectured scattering. Explain the precise meaning of the algebra notation and any quotient/degeneracy issues before translating it into modern notation.
4. **S. Penati and D. Zanon, *Supersymmetric, Integrable Toda Field Theories: The B(1,1) Model* (1992).** https://arxiv.org/abs/hep-th/9205045 ; Phys. Lett. B288 (1992) 297–305. Quantum W-currents; a higher-spin current conserved to all perturbative orders; one-loop masses and charges of the affine model. Distinguish the finite-algebra conformal theory from its affinization.
5. **J. M. Evans and J. O. Madsen, *Integrability vs. Supersymmetry* (1996).** https://arxiv.org/abs/hep-th/9608190 . A caution about particular `(1,0)` supersymmetrizations of ordinary conformal Toda theories. Do not turn their results into a universal no-go theorem for every supersymmetric extension of every affine Toda model.
6. **J. M. Evans and J. O. Madsen, *Dynkin Diagrams and Integrable Models Based on Lie Superalgebras* (1997), revised arXiv version 2000.** https://arxiv.org/abs/hep-th/9703065 ; Nucl. Phys. B503 (1997) 715–746. Main organizing source for root systems, bosonic subsectors, supersymmetry and reality. Use the corrected version. The introduction explicitly says: “These integrable models contain both bosonic and fermionic fields, but they need not be supersymmetric in general.” It also explains why the massive/massless distinction is subtler than for ordinary affine Toda theory.
7. **J. M. Evans and J. O. Madsen, *Quantum Integrability of Coupled N=1 Super Sine/Sinh-Gordon Theories and the Lie Superalgebra D(2,1;alpha)*.** https://arxiv.org/abs/hep-th/9712227 ; Int. J. Mod. Phys. A14 (1999) 2551–2580. Exceptional family, reality conditions, positive kinetic energy, bounded bosonic potentials and quantum conserved currents. Separate the algebra parameter from the Toda coupling.
8. **G. M. Gandenberger, *Trigonometric S-Matrices, Affine Toda Solitons and Supersymmetry*.** https://arxiv.org/abs/hep-th/9703158 ; Int. J. Mod. Phys. A13 (1998) 4553–4590. Already in `references.bib` as `gandenberger1998`, and already cited in chapter 22. Read the final discussion of minimal N=1/N=2 scattering and folding; make its precise relationship to Toda soliton amplitudes explicit without asserting a field-theory equivalence from an amplitude identity alone.
9. **K. Ito and M. Zhu, *ODE/IM correspondence and supersymmetric affine Toda field equations* (2022).** https://arxiv.org/abs/2206.08024 ; Nucl. Phys. B985 (2022) 116004. Linear problems for affine superalgebras with purely odd simple roots and their bosonic reductions. Distinguish results obtained from the conformal limit/bosonic reduction from claims about a fully established massive quantum correspondence.
10. **R. I. Nepomechie, *The boundary supersymmetric sine-Gordon model revisited* (2001).** https://arxiv.org/abs/hep-th/0103029 . Two-parameter boundary interactions and proposed breather reflection amplitudes. Read the boundary fermion construction and the exact scope of the conservation tests. For an optional finite-size connection: C. Ahn and R. I. Nepomechie, https://arxiv.org/abs/hep-th/0005170 .
11. **J. F. Gomes, L. H. Ymai and A. H. Zimerman, *Classical Integrable Super sinh-Gordon equation with Defects* (2006).** https://arxiv.org/abs/hep-th/0601014 . Also **A. R. Aguirre, J. F. Gomes, N. I. Spano and A. H. Zimerman, *N=1 super sinh-Gordon model with defects revisited* (2015)**, https://arxiv.org/abs/1412.2579 . Sources for a short classical defect extension, with explicit defect matrices and modified charges in the later paper. Do not describe these classical results as established quantum transmission matrices.

Check the primary papers' bibliographies for earlier construction and scattering sources. This list is a starting point, not an assertion of priority or completeness. Verify all new BibTeX metadata and reuse existing entries where possible.

## Required chapter content

### 1. What changes when the algebra is graded

Introduce even and odd generators, graded commutators, root parity and isotropic roots. Explain the role of the chosen simple-root system; do not transfer ordinary Lie-algebra Weyl-group arguments without justification. Explain why an auxiliary Lie superalgebra is a different statement from a spacetime supersymmetry algebra acting on physical states.

Give a small explicitly specified root/Gram-matrix example that will be used again. State the invariant bilinear form and reality conditions. Ordinary longest-root-squared-length-two conventions cannot simply be imposed on odd isotropic roots.

### 2. Construction in superspace and in components

Define light-cone coordinates, Grassmann derivatives, integration convention and scalar superfield. State how the source's N=1/N=2 notation translates into `(p,q)` notation in two dimensions. Give a graded zero-curvature representation for the chosen class and explain its hypotheses, including the role of purely odd simple roots for the manifestly supersymmetric construction.

Work through auxiliary-field elimination. Display the resulting bosonic potential and Yukawa couplings for a concrete example and check the relation between the root data and the bosonic subsector. Explain why this potential is not obtained merely by adding a fermion to the ordinary Toda Lagrangian. State the limits of any classification presented.

### 3. A worked super sine/sinh-Gordon example

Choose and verify the precise affine superalgebra/root presentation from the sources; do not guess its identification from the ordinary bosonic theory. Use one convention consistently from action to equations, masses and scattering. Show the equal boson/fermion masses around a supersymmetric vacuum, and explain the bosonic truncation and analytic continuation. The bosonic action may require a coupling rescaling to match chapter 3: derive the dictionary.

### 4. Quantum integrability and scattering

Give a substantive account of the Delius–Grisaru–Penati–Zanon results, identifying exactly which models and amplitudes are treated. Distinguish their hypotheses and supersymmetry properties model by model. Include a nontrivial formula or low-rank worked scattering example, with a reproducible check.

Explain multiplet mixing and the graded tensor product. Define whether the scattering operator includes the graded permutation before writing Yang–Baxter, unitarity or crossing equations. Explain the status and limitations of any scalar-factor/matrix-factor decomposition. Verify the crossing convention against the physical charge-conjugation matrix, not an ungraded template.

Use the B(1,1) current calculation as a concrete example of the evidence for quantum integrability. Summarize precisely what the paper proves; reproducing the full all-loop computation is not required. A classical or finite-order script does not verify the all-loop claim.

### 5. N=2 theories, solitons and reality

Present the main massive N=2 family and one explicit, checked soliton or BPS calculation if the source and length permit. Define supercharges, central charge and boundary conditions before using a BPS bound. Distinguish classical statements from quantum statements, including possible anomalies. Do not import bosonic soliton mass or charge formulas unchanged.

Give a concise account of the exceptional D(2,1;alpha) family and explain why real forms and positive kinetic terms matter. No broad claim of unitarity follows from a real classical spectrum. Cross-reference chapter 11 without identifying its antilinear symmetry with supersymmetry.

### 6. Connections, exercises and notes

End the main exposition with a short connection to ODE/IM and pointers to the boundary and defect extensions elsewhere in the book. Include about six exercises, progressing from a graded-bracket or auxiliary-field calculation to a mass/scattering check. Supply enough data that each exercise is self-contained. In Notes and references, separate established constructions, conjectured scattering data and issues not treated here. Attribute the results; no claims of novelty are requested.

## Convention dictionary and verification

Create `briefs/superalgebra_toda_conventions.md` and `briefs/draft_superalgebra_toda_report.md`. Keep a checkpoint in `briefs/draft_superalgebra_toda_checkpoint.md` for long work. In the conventions document, record the source equation, book equation and explicit transformation for every imported convention: metric and light-cone normalization; left/right Grassmann derivatives and integration order; kinetic and Yukawa signs; root parity and invariant form; mass and coupling scales; algebra names; supercharge algebra; graded permutation; tensor order; crossing and charge conjugation.

If a quantum-group comparison is included, read `eq-q-physical` and the chapter 7/Part VI coproduct dictionary. Do not identify a superalgebra deformation parameter with the bosonic physical q, or dismiss a sign change as a gauge transformation, without a derivation. Use generic irrational parameter choices where generic q is needed.

Add focused scripts under `textbook_code/`, preferably a small `super_toda/` subdirectory with its own README if several scripts are needed. Index the directory/scripts in `textbook_code/README.md`. Cite each script by full repository path at the formula it checks. Print one PASS or FAIL line per check and make a failed check return a nonzero exit status.

Required meaningful checks, adjusted to the actual examples chosen:

- Derive the component potential after auxiliary-field elimination and compare with the displayed formula.
- Compute a bosonic Hessian and fermion mass matrix in the same convention; check the claimed multiplet masses and bosonic truncation.
- Check the selected graded Lax representation or an explicitly displayed classical conserved current, including odd signs. A custom Grassmann representation must document its ordering and check its algebra.
- Test the displayed scattering example against the particular unitarity, symmetry and graded Yang–Baxter identities claimed. Test crossing if a complete scalar factor and crossing prescription are actually supplied. Do not claim a complete physical S-matrix check from its matrix part alone.
- Verify any explicit soliton/BPS formulas against equations and asymptotic boundary data.

Prefer exact symbolic identities where practical. For numerical tests, give parameter domains, tolerances and points away from poles; count poles structurally if making pole claims. Do not extend shared bosonic helpers into graded code unless necessary; new local helpers reduce the risk of changing existing results.

## Integration into the existing book

Use section ids below; filenames describe the baseline before renumbering.

| Existing location | Required change |
|---|---|
| Chapter 6, `sec-affine-algebras` / `sec-lie-algebraic-data` | Add a short introduction or forward pointer to affine superalgebras, parity and isotropic roots. Keep the detailed construction in the new chapter. If still under review, record the exact deferred patch. |
| Chapter 9, `sec-lagrangians` / `sec-lax-pair` | Clarify the scope of the ordinary bosonic formulas and point to the new construction. Do not rewrite the current derivation as if it already applied to superalgebras. |
| Chapters 12–14, appropriate Notes | Connect the quantum superalgebra literature to perturbative mass tests and particle scattering where useful; avoid implying that all their diagonal bosonic formulas carry over. |
| Chapter 15, `sec-ode-im` | Add a concise, accurately scoped paragraph citing Ito–Zhu. |
| Chapter 22, `sec-known-smatrices` and Notes | Add the verified content of Gandenberger's supersymmetry connection. State precisely the relationship between amplitudes; avoid an unsupported equivalence of theories. |
| Chapter 29, `sec-boundary-sine-gordon` | Add a short subsection on supersymmetric boundaries, conserved boundary supercharges and any necessary boundary degrees of freedom, citing Nepomechie. It should make sense to a reader following only the boundary reading path. |
| Chapter 37, `sec-defect-sine-gordon` | Add a short subsection on super sinh-Gordon defects, extra fermionic variables where present, frozen Bäcklund relations and modified classical charges. Cite the 2006/2015 sources and distinguish classical from quantum results. |
| Introduction, README contents, Part IV preamble, outlook | Update scope, chapter ranges and reading paths. Remove or qualify universal wording such as “every affine algebra” where it refers only to the ordinary bosonic class. Do not blindly replace every occurrence. |
| `references.bib`, `briefs/README.md`, code READMEs | Add verified bibliography entries and index the new drafting record and scripts. |

If adding an explicit new boundary/defect formula, check it under the same script-at-claim rule. Short accurately attributed descriptions may be preferable to introducing several more full calculations.

## Numbering and file migration

Insertion makes the new chapter 16. Existing chapters 16–50 become 17–51. Parts V–IX retain their part numbers. Preserve all existing section/equation/statement ids. Dated briefs retain their historical numbers.

Audit filenames, `_quarto.yml`, README tables, introduction paths, figure script names and paths, figure captions, script imports, lecture-note links, `textbook_code/figures/style.py`, and other references to renamed paths. Search the actual repository; do not assume this list is exhaustive. Keep a mapping in the drafting report and update the numbering history in `briefs/README.md`. Avoid gratuitous regeneration of figures.

`AGENTS.md` requests a pure-rename commit before text changes, while commits require an explicit editor request. Do not commit on the strength of this handover. If the editor has not separately authorized commits, prepare the rename mapping and preserve a reproducible separation of mechanical renames from content changes (for example, a saved rename manifest and separate patches). Record the pending two-commit sequence for the editor. Continue the authorized drafting work; do not silently make an exception to the no-commit rule.

## Validation and handback

Run every new or affected check script and `python3 check-refs.py`. Render through `./render-scratch.sh`, never `quarto render` in the working checkout. Check the chapter's HTML, especially superscripts, graded brackets, parity labels, root tables and cross-references. Use the scratch PDF workflow if introducing notation likely to expose LaTeX issues. Report what was actually run, with outputs and unresolved failures.

The drafting report must contain:

- Final paths and stable ids, rename mapping and any integration edits deferred because of active review.
- The source-to-book convention dictionary, source equation/page locators, and what was checked by derivation, script or source reading.
- Exact script commands, outputs, check counts and approximate runtimes.
- An `## Integration` section recording which suggestions from this handover were verified, modified or rejected, and why.
- Every conjectural formula, unverified source claim and limitation; any proposed correction to a paper must quote the exact source line and meet the chapter-plus-Appendix-D rule.
- Appendix E impact. A literature chapter normally adds no original-result item; update the appendix if existing scope/status changes or a genuinely original result is introduced.
- A proposed independent referee brief against the finished state. Prioritize graded component/Lax signs, normalization, scattering/crossing conventions, the N=2 family, the scope of quantum-integrability claims and fidelity to sources. Include evidence requirements and a fixed reviewed commit once the editor has committed the draft. Do not launch a referee or claim that this handover is an independent verification of the resulting chapter.

Deliver the implemented draft and integration edits, rather than only a new outline. If source access prevents a particular advanced calculation, keep the chapter coherent with verified material, state the limitation, and identify the missing work precisely in the report.
