# Drafting report: chapter 16, affine Toda theories based on Lie superalgebras

Date: 2026-10-23. Model: Claude Opus 5.5 (drafting, not refereeing). Worktree `.claude/worktrees/superalgebra-toda`, branch `claude/superalgebra-toda`, based on origin/main `b7087eb` (the handover was written at `745586e`; nothing it relies on changed in between). Handover: `briefs/draft_superalgebra_toda.md` (Codex). Nothing has been committed or pushed.

## Summary

- New chapter `part4-real-coupling/16-superalgebra-toda.qmd` (`sec-super-toda`), about 6,100 words: graded algebras and their roots; superspace, the component Lagrangian and the bosonic sector; the graded Lax pair; the super sinh-Gordon theory as worked example (realization of the Lax pair, supercurrents); masses that do not float, the boson–fermion S-matrix of Delius–Grisaru–Penati–Zanon (DGPZ) and the $C^{(2)}(n+1)$ S-matrix, quantum conserved currents; the $N=2$ family, its solitons and reality conditions, and $D(2,1;\alpha)^{(1)}$; connections; six exercises; Notes.
- Six check scripts in `code/super_toda/` (87 checks, all passing), with a README.
- Integration edits in chapters 6, 12, 14, 15, 23, 30, 38, the introduction, preface, Part IV preamble, outlook, README, Appendices D and E, the code READMEs, `references.bib` (23 entries) and `briefs/README.md`.
- Chapters 16–50 renumbered to 17–51. The pure-rename commit is **staged** in the index; all other edits are unstaged (see "Commit sequence").
- Three misprints in the sources, entered in Appendix D; two signs fixed by our own checks.

## Final paths and ids

**New chapter** `part4-real-coupling/16-superalgebra-toda.qmd`. Ids: `sec-super-toda`, `sec-super-roots`, `eq-super-bracket`, `eq-super-chevalley`, `tbl-super-examples`, `sec-super-superspace`, `eq-super-derivatives`, `eq-super-field`, `eq-super-action`, `eq-super-eom`, `eq-super-components`, `prp-super-bosonic` [Established], `prp-super-susy` [Established], `sec-super-lax`, `eq-super-lax`, `eq-super-zero-curvature`, `sec-super-shg`, `eq-ssg-action`, `eq-ssg-components`, `eq-ssg-realization`, `eq-ssg-currents`, `sec-super-scattering`, `sec-super-fixed-masses`, `eq-dgpz-masses`, `sec-super-dgpz`, `eq-dgpz-operator`, `eq-dgpz-ybe`, `eq-dgpz-amplitudes`, `eq-dgpz-pi`, `eq-dgpz-unitarity`, `eq-dgpz-san`, `eq-dgpz-cn`, `cnj-dgpz-smatrix` [Conjecture], `sec-super-quantum-currents`, `sec-super-n2`, `eq-n2-condition`, `eq-n2-masses`, `eq-n2-central`, `cnj-n2-sine-gordon` [Conjecture], `sec-super-connections`, `sec-super-toda-exercises`, `exr-super-gl11`, `exr-super-auxiliary`, `exr-super-masses`, `exr-super-b11`, `exr-super-n2`, `exr-super-dgpz`, `sec-super-notes`.

**New sections elsewhere:** `sec-bsg-super` (chapter 30, `part7-half-line/30-boundary-sine-gordon.qmd`), `sec-sgd-super` (chapter 38, `part8-defects/38-defect-sine-gordon.qmd`).

**Scripts:** `code/super_toda/{superspace,components,lax_currents,dgpz_smatrix,n2_solitons,gandenberger_scalar}.py`, `code/super_toda/README.md`.

**Records:** this report; `briefs/superalgebra_toda_conventions.md` (source-to-book dictionary); `briefs/draft_superalgebra_toda_checkpoint.md`; `briefs/superalgebra_toda_renumber.py` and `briefs/superalgebra_toda_renames.tsv`; `briefs/referee_superalgebra_toda.md` (proposed referee brief).

## Renumbering and commit sequence

`briefs/superalgebra_toda_renumber.py` performed 77 renames with `git mv`: the 35 chapter files 16–50 → 17–51 (Parts V–IX keep their part numbers; the outlook becomes `outlook/51-outlook.qmd`), 21 chapter-numbered SVGs `part*/figures/chNN_*` and 21 figure scripts `code/figures/chNN_*` with NN ≥ 16. It fixed only paths that would otherwise break, in 43 files: the chapter list of `_quarto.yml`, image links and script paths in the renamed chapters, output names in the renamed figure scripts (and the chapter path read by the gap-map script), the figure rows of `code/README.md`, the path in `code/sectors/bloch/README.md`, and the current-file columns of `briefs/README.md` (its dated history paragraph keeps old numbers). It staged exactly these files. The full mapping is in `briefs/superalgebra_toda_renames.tsv`.

AGENTS.md asks for a pure-rename commit followed by a text commit, and the handover does not authorize commits. The two commits are prepared but not made:

1. **Commit 1 (pure rename)**: the index as it stands. `git diff --cached -M --stat` shows 83 files (77 renames, path fixes). Suggested message: "Insert chapter 16 (superalgebra Toda): renumber chapters 16–50 to 17–51, file renames".
2. **Commit 2 (text)**: everything else. Stage by name, never with `git add -A`: the new chapter, `code/super_toda/`, the edited chapters and appendices, `references.bib`, the READMEs, the regenerated `figures/intro_f1_reading_paths.svg` and `part9-foundations/figures/ch45_f33_gap_map.svg`, `code/figures/intro_f1_reading_paths.py`, and the briefs files (including the handover, which was untracked in the main checkout).

Hard-coded chapter numbers updated in the text commit: the dependency diagram and reading paths of the introduction (and the reading-paths figure, regenerated from the text; its script's part ranges updated), the parts table of `README.md`, the part headings of Appendix E, two "ch. NN" references in chapter 41 (replaced by `@sec-` references), the gap-map caption's alt text (and the figure, whose column labels are chapter numbers), and the half-line range in `code/README.md`. `briefs/README.md` records the renumbering. Lecture notes contain no chapter paths.

**Coordination.** The open Part III integration branches `part3-int-ch9-10` and `part3-int-ch11` (no commits since 2026-10-15) were told to record edits to Part IX and chapter 21 (now 22) as pending; those files are renamed here, so their patches will need the new paths. The rename-only commit keeps git's rename detection working for a later merge.

## Chapter 9 edit, applied while Part III is frozen

Chapter 9 is frozen (`briefs/README.md`, "Part III: two-phase referee pass": files under review until the reports are integrated). The integration of those reports stalled on 2026-10-15 with no book edits and was restarted on 2026-10-23 in a separate session. On 2026-10-23 the editor decided to apply this two-sentence edit now, in its own commit on `claude/superalgebra-toda` ("Chapter 9: point to the superalgebra chapter"), rather than wait. Edits to `part3-classical/09-lagrangian.qmd`:

- In `sec-lagrangians`, appended to the paragraph defining $\vec\phi$, the sums and the $n_j$ (ending "… because $\sum_jn_j\vec\alpha_j=0$."): "Throughout this chapter $\hat g$ is an affine Lie algebra. Replacing it by an affine Lie superalgebra gives theories with fermions, which are the subject of @sec-super-toda; the formulas of this chapter do not apply to them without change."
- In `sec-lax-pair`, appended after "`code/toda_classical.py` checks the equivalence symbolically for $a_2^{(1)}$.": "The same construction works in superspace for an affine Lie superalgebra, with odd generators for the odd simple roots (@sec-super-lax)."

**The Part III integration must merge this commit.** Its branch `part3-referee-integration` was cut before it. Neither sentence changes a claim under review. In the four Part III reports, the only replacement text for these two sections is Codex's for the sentences that specify the generators (`09-lagrangian.qmd:57` at 2d36c9c, phase 1, item 5), a different paragraph. Gemini's phase 2 verifies the paragraph that ends at the second insertion (09:59–73) as correct, and no report targets the paragraph that ends at the first (09:28). Whichever of the two branches reaches main second must keep both sentences when it resolves any conflict in these two paragraphs.

**Resolved 2026-10-24.** The integration was merged into main as PR #44 (with PR #45), and main was then merged into this branch. The `sec-lax-pair` paragraph conflicted, because #44 extended the same sentence ("… and, with $3\times3$ matrices, for $a_2^{(2)}$."); the resolution keeps #44's text and appends the superalgebra sentence. The `sec-lagrangians` sentence merged without conflict. The other conflicts (`index.qmd`, `references.bib`, `briefs/README.md`, the gap-map SVG) were resolved by keeping both sides; the gap map was regenerated. Two stale details from main were corrected in the merge: the gap-map alt text now lists H6 among the orange rows (#44 made H6 partly open), and the docstring of `code/an_soliton_scalar.py` (#45) says chapter 22 for `sec-scalar-factor`.

Chapters 10 and 11 needed no edit; chapter 16 cross-references `sec-classical-solitons`, `eq-soliton-mass` and `sec-theta-classical` without changing them.

## Sources and locators

Read in full or in the relevant sections (PDFs and text in the session scratchpad only):

| Source | Read | Used for |
|---|---|---|
| Evans–Madsen 1997, hep-th/9703065 v2 (2000) | all | §2 node colours, Figure 4 sign rules; (2.1); (3.1)–(3.5); Props. 1–7; Tables 1, 2, 4, 5; (5.1)–(5.12); (7.10)–(7.12) |
| Evans–Madsen 1999, hep-th/9712227 | all | (1.1)–(1.6); (2.7)–(2.12); (3.13)–(3.32) and the reality table; §4 limits; §5 (5.35)–(5.60); §6 quantum currents |
| Evans–Hollowood 1992, hep-th/9207073 | all | (1.1)–(1.9); (2.1)–(2.14); §3; (5.1)–(5.10) |
| Evans–Hollowood 1991, NPB 352 (scan from the editor) | §§3.1–3.3, 4, 7 | superspace Lax construction (3.9)–(3.18); centre of $\mathfrak{sl}(m|m)$; deformations. Erratum NPB 382 (1992) 662 **not obtained**; no formula is taken from this paper alone |
| DGPZ, NPB 359 (1991) 125 (scan) | §§1–3, 5, App. A, App. B, refs | (1.1)–(1.5), (2.1)–(2.51), (3.1)–(3.17), (A.1)–(A.9), (B.1)–(B.6) |
| DGPZ, PLB 256 (1991) 164 (scan) | all | (1)–(33), the $C^{(2)}(n+1)$ one-loop double pole |
| Penati–Zanon 1992, hep-th/9205045 | all | (4), (16)–(18), (24), (28), (33)–(35), conclusion |
| Evans–Madsen 1996, hep-th/9608190 | §§1–4 | scope of the $(1,0)$ result |
| Gandenberger 1997, hep-th/9703158 | §§3.2 (3.28), 4 (4.4), 5 | scalar-factor identities (5.6), (5.10), (5.14), (5.17), (5.21) |
| Ito–Zhu 2022, 2206.08024 | abstract, §§1–3, 5 | scope of the ODE/IM result |
| Nepomechie 2001, hep-th/0103029 | all | (1)–(35) |
| Gomes–Ymai–Zimerman 2006, hep-th/0601014 | all | (2.1)–(2.55), §4 |
| Aguirre–Gomes–Spano–Zimerman, 1412.2579 | abstract, introduction | scope |
| DGZ 1992, hep-th/9201067 | §§1, 2.3 | relation of $a_{2n-1}^{(2)}$ to $A^{(2)}(0,2n-1)$ |

Not read: Ahn–Nepomechie hep-th/0005170 (optional in the handover; not used); Kobayashi–Uematsu, Fendley–Intriligator, Schoutens, Ahn 1991, Hollowood–Mavrikis (cited only for attributions relayed from Evans–Hollowood 1992, DGPZ and Gandenberger, labelled as such).

The source-to-book dictionary, with equation locators, is `briefs/superalgebra_toda_conventions.md`.

## Verification record

**By script** (`code/super_toda/`, outputs below):

- Component Lagrangian (EM97 (3.4)) from the superspace action for generic rank-2 data with a symbolic Gram matrix, with all signs; auxiliary field; super sinh-Gordon Lagrangian, truncation to the book's $a_1^{(1)}$ (same $m$, $\beta$), masses $2m$; the field equation equals the Euler–Lagrange equations (super sinh-Gordon).
- $M_B^2=M_F^2$ for four purely odd systems; $B(1,1)^{(1)}$ masses and opposite kinetic signs; $D(2,1;\alpha)^{(1)}$ Cartan matrix, masses, kinetic signs under the three twists, Yukawa coupling (EM normalization); $A(1,1)^{(1)}$ with twisted reality; DGPZ masses for $n=2,3,4$.
- The $N=2$ condition: solutions for $A(1,1)^{(1)}$, none for $C(2)^{(2)}$, $B(1,1)^{(1)}$, $D(2,1;\alpha)^{(1)}$.
- Superspace zero curvature with the $\mathfrak{sl}(2|1)$ realization; the $\mathfrak{osp}(1|2)$ obstruction; the two supercurrents on shell, by reduction in jet space.
- DGPZ S-matrix: $T$ (with a Stirling-series tail, accurate to about $10^{-12}$), crossing, unitarity relation, no poles in the strip; (2.7) and the braid equation, failure of the ungraded version; unitarity including $\check S\check S^\dagger=1$; crossing and the spinor phase; $B-F=2T$; the free-fermion condition; poles, residue signs, rank one; tree level against DGPZ (2.5) and against our own tree amplitude from the potential; the bootstrap (2.45)$=$(2.47); $C^{(2)}(n+1)$ fermion amplitudes for $n=2,\dots,5$.
- $N=2$ solitons for $m=3,4$: Cartan matrix, $V_{\rm real}+V_{\rm imag}$, field equation, energy by quadrature, central-charge formula, first-order equation.
- Gandenberger (5.10) up to sign, and (5.17).

**By derivation in the text:** the zero-curvature argument (the script confirms it in one realization); $M_B^2=M_F^2$ (exercise); the $N=2$ condition implies null roots and a bipartite graph; the supercharge conservation law's components ($D_-J=D_+\bar J$ gives EM98 (5.48)–(5.51), checked by hand); the crossing of tree amplitudes (2.5) under $\theta\to i\pi-\theta$.

**By source reading only (attributed, not checked by us):** the classifications in `prp-super-bosonic`, `prp-super-susy` and the purely-massive statement (EM97 Props. 1, 2, 4, 6); the one-loop mass results (GPZ, as stated by DGPZ); the one-loop double-pole analyses (DGPZ Sec. 4, letter); Penati–Zanon's all-order spin-3 current and one-loop results; Evans–Madsen's quantum spin-7/2 currents; DGZ 1992 Sec. 2.3; Ito–Zhu's ODEs; Nepomechie's boundary results; the defect papers; the $N=2$ sine-Gordon S-matrix conjecture; DGPZ's remark on $n=1$; the existence of an S-matrix for super sine-Gordon (Ahn).

## Script commands and outputs

Run from the worktree root (times on this machine, with a Quarto render running in parallel):

```
python3 code/super_toda/superspace.py
python3 code/super_toda/components.py
python3 code/super_toda/lax_currents.py
python3 code/super_toda/dgpz_smatrix.py
python3 code/super_toda/n2_solitons.py
python3 code/super_toda/gandenberger_scalar.py
```

Exit codes, counts and times of the final run (all 87 checks pass):

```
superspace exit=0 time=2s PASS=9 FAIL=0
components exit=0 time=26s PASS=31 FAIL=0
lax_currents exit=0 time=10s PASS=11 FAIL=0
dgpz_smatrix exit=0 time=61s PASS=27 FAIL=0
n2_solitons exit=0 time=77s PASS=7 FAIL=0
gandenberger_scalar exit=0 time=6s PASS=2 FAIL=0
```

Full output:

`superspace.py`:
```
PASS Grassmann generators anticommute and square to zero
PASS d/dtheta is an odd derivation acting from the left
PASS D_+^2 = -i d_+, D_-^2 = -i d_-, {D_+, D_-} = 0 on a superfield
PASS {Q_+, D_+} = {Q_+, D_-} = 0 and Q_+^2 = i d_+
PASS berezin(th+ th-) = i and berezin(1) = berezin(th+) = 0
PASS sexp agrees with the Taylor series on a nilpotent element
PASS gl(1|1): {E12, E21} = E11 + E22 for odd matrix units
PASS envelope: [a E12, b E21] = -a b {E12, E21} with odd a, b
PASS gl(2|1): graded Jacobi identity on matrix units
```

`components.py`:
```
PASS generic rank-2 data: int d^2theta and sigma elimination give eq-super-components (Evans-Madsen (3.4)), including the fermion kinetic and Yukawa signs
PASS the auxiliary field is sigma = sum_odd mu_j alpha_j exp(alpha_j.phi) (as a vector)
PASS super sinh-Gordon: L = 1/2 d+phi d-phi + (i/2)(psi+ d-psi+ + psi- d+psi-) - (mu^2/4beta^2)(e^{sqrt2 beta phi} + e^{-sqrt2 beta phi} - 2) + i mu psi+ psi- cosh(beta phi/sqrt2) (eq-ssg-components)
PASS bosonic truncation psi = 0 is the a_1^(1) theory of eq-imaginary-potential with m = mu/2 and the same beta
PASS super sinh-Gordon: boson mass^2 = mu^2 = (2m)^2, the a_1^(1) particle mass of sec-tree-level
PASS super sinh-Gordon: fermion mass = mu (Yukawa coefficient at phi = 0), equal to the boson mass
PASS C(2)^(2): phi = 0 is stationary and M_B^2 = M_F^2 (exr-super-masses)
PASS B(1,1)^(1): the realization reproduces the Evans-Madsen Cartan matrix and 2a1+a2+a3 = 0
PASS B(1,1)^(1): phi = 0 is stationary and M_B^2 = M_F^2 (exr-super-masses)
PASS B(1,1)^(1): boson masses^2 in the ratio 1 : 4 (masses m and 2m; Penati-Zanon m_1^2 = 2, m_2^2 = 8)
PASS B(1,1)^(1): the heavy and light sinh-Gordon sectors have kinetic terms of opposite sign (Evans-Madsen (5.4))
PASS D(2,1;alpha): alpha_i^2 = 0, a_12 = alpha, a_13 = 1, a_23 = -1-alpha (Evans-Madsen (2.7)) and the affine Cartan matrix (2.12)
PASS D(2,1;alpha)^(1): phi = 0 is stationary and M_B^2 = M_F^2 (exr-super-masses)
PASS D(2,1;alpha)^(1): boson masses m(1+alpha), m, m alpha with m = 2 mu (Evans-Madsen (3.23))
PASS D(2,1;alpha)^(1): fermion mass matrix is diagonal in the eps basis with entries -m(1+alpha)... up to sign: |M_F| = m(1+alpha), m, m alpha
PASS D(2,1;alpha)^(1), alpha > 0: only the twisted reality condition (1) gives a definite (negative) kinetic term (Evans-Madsen table after (3.28))
PASS D(2,1;alpha)^(1): in Evans-Madsen field normalization the Yukawa coupling of psibar_1 psi_2 phi_3 has modulus 4 mu sqrt(alpha(1+alpha)) (Evans-Madsen after (3.32))
PASS A(1,1)^(1): realization of the Evans-Madsen Cartan matrix with two null eigenvectors
PASS A(1,1)^(1), gamma = 1: phi = 0 is stationary and M_B^2 = M_F^2 (exr-super-masses)
PASS A(1,1)^(1), gamma = 1: the two bosons have equal mass, independent of gamma (Evans-Hollowood 1992 sec. 2, Evans-Madsen (5.7))
PASS A(1,1)^(1), generic gamma: phi = 0 is stationary and M_B^2 = M_F^2 (exr-super-masses)
PASS A(1,1)^(1), generic gamma: the two bosons have equal mass, independent of gamma (Evans-Hollowood 1992 sec. 2, Evans-Madsen (5.7))
PASS A(1,1)^(1) with twisted reality: U = sinh-Gordon(phi^1) + sine-Gordon(w) and the kinetic term becomes positive (Evans-Madsen (7.10))
PASS A^(2)(0,2n-1), n = 2: phi = 0 stationary; boson and fermion masses^2 are those of Delius-Grisaru-Penati-Zanon (2.2)/(3.2) with h = 3
PASS C^(2)(n+1), n = 2: phi = 0 stationary; boson and fermion masses^2 are those of Delius-Grisaru-Penati-Zanon (2.2)/(3.2) with h = 4
PASS A^(2)(0,2n-1), n = 3: phi = 0 stationary; boson and fermion masses^2 are those of Delius-Grisaru-Penati-Zanon (2.2)/(3.2) with h = 5
PASS C^(2)(n+1), n = 3: phi = 0 stationary; boson and fermion masses^2 are those of Delius-Grisaru-Penati-Zanon (2.2)/(3.2) with h = 6
PASS A^(2)(0,2n-1), n = 4: phi = 0 stationary; boson and fermion masses^2 are those of Delius-Grisaru-Penati-Zanon (2.2)/(3.2) with h = 7
PASS C^(2)(n+1), n = 4: phi = 0 stationary; boson and fermion masses^2 are those of Delius-Grisaru-Penati-Zanon (2.2)/(3.2) with h = 8
PASS A(1,1)^(1): exactly two solutions lambda_j = +-i(-1)^j of (lambda_i+lambda_j)a_ij = 0 (Evans-Hollowood 1992 (2.9))
PASS no second supersymmetry for C(2)^(2), B(1,1)^(1) or D(2,1;alpha)^(1) (alpha = 0.7)
```

`lax_currents.py`:
```
PASS sl(2|1) realization of the C(2)^(2) relations with alpha^2 = 1/2 (eq-ssg-realization); h has zero supertrace
PASS osp(1|2): {F+, F-} is diagonal, but {F-, F-} = 2F-^2 and {F+, F+} are nonzero, so {e_0, f_1} cannot vanish when e_0 and f_1 are both built from F-
PASS superfield equation, theta^0 component: sigma = (mu/beta) a (e^{beta a phi} - e^{-beta a phi}), the same as the Euler-Lagrange equation of sigma
PASS superfield equation, theta^+ theta^- component: d_+d_-phi equals the Euler-Lagrange equation of the component Lagrangian (Evans-Madsen (3.5) in these conventions)
PASS superfield equation, theta^- component: d_-psi_+ = i yuk psi_-, the Euler-Lagrange equation of the component Lagrangian
PASS superfield equation, theta^+ component: d_+psi_- = -i yuk psi_+
PASS zero curvature: D_+A_- + D_-A_+ + {A_+, A_-} = (i beta/|alpha|) (field equation) h for all lambda (eq-super-zero-curvature), with eta_0 = eta_1 = -i mu
PASS super sinh-Gordon: D_-T = D_+Tbar on shell with T = D_+Phi D_+^2Phi - D_+^3Phi and Tbar = -i m D_+ e^{-Phi} (eq-ssg-currents, spin 3/2)
PASS super sinh-Gordon: D_-J = D_+Jbar on shell for J = T D_+T and Jbar = i m (D_+Phi (D_+^2Phi)^2 - D_+^2Phi D_+^3Phi) e^{-Phi} (eq-ssg-currents, spin 7/2)
PASS the printed Evans-Madsen 1999 (5.44)-(5.46) do not hold, either as printed or with the D_+- exchanged (a misprint; sec-app-corrections)
PASS at psi = 0 the theta^+ component of T D_+T is ((d_+phi)^2 - d_+^2 phi)^2, which differs from (d_+^2 phi)^2 + (d_+phi)^4 by a total d_+ derivative (the sinh-Gordon density of eq-spin3-current)
```

`dgpz_smatrix.py`:
```
PASS n = 2: T(theta) = Pi(theta)Pi(i pi - theta) obeys crossing and the unitarity relation (2.37), and is bounded in the physical strip (max |T| on the grid 0.946)
PASS n = 3: T(theta) = Pi(theta)Pi(i pi - theta) obeys crossing and the unitarity relation (2.37), and is bounded in the physical strip (max |T| on the grid 0.976)
PASS n = 2: the Yang-Baxter equations (2.7) of Delius et al. and their B <-> F partners
PASS n = 2: braid Yang-Baxter equation for the graded operator Sc (eq-dgpz-ybe)
PASS n = 2: with the ungraded assignment (free limit = ordinary permutation) the braid Yang-Baxter equation fails
PASS n = 2: unitarity Sc(theta)Sc(-theta) = 1, and Sc Sc^dagger = 1 for real theta
PASS n = 2: crossing B, F, T symmetric and D(i pi - theta) = i R(theta)
PASS n = 2: B - F = 2T and the free-fermion condition a1 a2 + b1 b2 = c1 c2 + d1 d2 of the eight-vertex weights of Sc, i.e. BF + T^2 + D^2 - R^2 = 0 since d1 d2 = D(-D); the sign printed in (B.1) of Delius et al., BF + T^2 - D^2 - R^2 = 0, fails
PASS n = 3: the Yang-Baxter equations (2.7) of Delius et al. and their B <-> F partners
PASS n = 3: braid Yang-Baxter equation for the graded operator Sc (eq-dgpz-ybe)
PASS n = 3: with the ungraded assignment (free limit = ordinary permutation) the braid Yang-Baxter equation fails
PASS n = 3: unitarity Sc(theta)Sc(-theta) = 1, and Sc Sc^dagger = 1 for real theta
PASS n = 3: crossing B, F, T symmetric and D(i pi - theta) = i R(theta)
PASS n = 3: B - F = 2T and the free-fermion condition a1 a2 + b1 b2 = c1 c2 + d1 d2 of the eight-vertex weights of Sc, i.e. BF + T^2 + D^2 - R^2 = 0 since d1 d2 = D(-D); the sign printed in (B.1) of Delius et al., BF + T^2 - D^2 - R^2 = 0, fails
PASS crossing phase of the Majorana spinors: u(theta + i pi) = i v(theta)
PASS n = 2: simple poles of B, F, D at theta = i pi(2(n-a)-1)/h, i.e. s = m_a^2 = 8 sin^2(pi a/h), with residues i x positive in B and F and a residue matrix of rank one (Omega = 0.1, 0.31, 0.6); poles of R, B, F at theta = 2 pi i a/h, a = 1..n-1
PASS n = 3: simple poles of B, F, D at theta = i pi(2(n-a)-1)/h, i.e. s = m_a^2 = 8 sin^2(pi a/h), with residues i x positive in B and F and a residue matrix of rank one (Omega = 0.1, 0.31, 0.6); poles of R, B, F at theta = 2 pi i a/h, a = 1..n-1
PASS n = 2: order-Omega terms of T, R, D, B equal the tree amplitudes (2.5) of Delius et al. with beta^2 = 4 pi Omega
PASS n = 3: order-Omega terms of T, R, D, B equal the tree amplitudes (2.5) of Delius et al. with beta^2 = 4 pi Omega
PASS n = 2: the order-Omega term of B equals the tree amplitude computed here from the bosonic potential sum q_j e^(alpha_j.phi) + e^(2 alpha.phi) (contact, s-, t- and u-channel exchange)
PASS n = 3: the order-Omega term of B equals the tree amplitude computed here from the bosonic potential sum q_j e^(alpha_j.phi) + e^(2 alpha.phi) (contact, s-, t- and u-channel exchange)
PASS n = 2: the bootstrap amplitude (2.45) equals the block product prod_p {h/2 - a + 2p - 1} (2.47) with B = 2 Omega
PASS n = 3: the bootstrap amplitude (2.45) equals the block product prod_p {h/2 - a + 2p - 1} (2.47) with B = 2 Omega
PASS C^(2)(2+1): the fermion amplitudes (3.10) are unitary; crossing exchanges S_chichi and S_chichibar (h = 4)
PASS C^(2)(3+1): the fermion amplitudes (3.13) are unitary; crossing maps S_++ and S_+- to themselves (h = 6)
PASS C^(2)(4+1): the fermion amplitudes (3.10) are unitary; crossing exchanges S_chichi and S_chichibar (h = 8)
PASS C^(2)(5+1): the fermion amplitudes (3.13) are unitary; crossing maps S_++ and S_+- to themselves (h = 10)
```

`n2_solitons.py`:
```
PASS m = 3: null roots on a cycle with alternating products +-1, rank 4 = 2m-2, null vectors (1,0,1,0,..), (0,1,0,1,..)
PASS m = 4: null roots on a cycle with alternating products +-1, rank 6 = 2m-2, null vectors (1,0,1,0,..), (0,1,0,1,..)
PASS m = 3: at psi = 0 the potential is V_real(phi^(1)) + V_imag(phi^(2)) of a_2^(1) (eq-imaginary-potential) with the book m equal to mu
PASS the one-soliton (5.3) solves the complex Toda equation (5.2) (m = 3, 4; a = 1, 2)
PASS its energy at rest is M_a = 4 m mu sin(pi a/m)/beta^2 (5.4) = 2 h m_a/beta^2 with m_a = 2 mu sin(pi a/h)
PASS the central-charge formula (5.9) from the asymptotic values gives the same M_a
PASS for a != m/2 the soliton satisfies the first-order equation (5.7) with J alpha_j = i(-1)^j alpha_j, sigma = tan(pi a/2m + pi/4) and lambda = -sqrt((1-v)/(1+v)); with lambda = +sqrt((1-v)/(1+v)) it fails for either sign of J
```

`gandenberger_scalar.py`:
```
PASS N = 2: the partial products of (5.6) and of (3.28) at lambda = 1/h differ by the factor -(M - mu)/(M + mu), so G^(N=2)_11 = -F_11 (Gandenberger (5.10) holds up to the sign -1)
PASS N = 1: at lambda = 1/h the product (4.4) for a_2n^(2) equals the product in (5.14) factor by factor, so (5.17) holds as printed (n = 1, 2)
```

Also rerun: `python3 code/figures/intro_f1_reading_paths.py` (7 PASS, after updating its path table) and `python3 code/figures/ch45_f33_gap_map.py` (6 PASS); `python3 check-refs.py` (all PASS; the one WARN about `citation.qmd` predates this work).

**Render.** `./render-scratch.sh <scratch>` (HTML): `quarto render` exited with 0 after 612 s; "No unresolved cross-references in the HTML"; one warning, "Could not convert TeX math U_q(\widehat{sl}_2)", from chapter 3 and present before this work. In the rendered chapter, statements are numbered Proposition 16.1–16.2 and Conjecture 16.1–16.2, and all `@sec-`, `@eq-`, `@tbl-` references resolve. `./render-scratch.sh --pdf <scratch>`: `quarto render` exited with 0 after 796 s; no unresolved cross-references; the same single chapter-3 warning; PDF of 669 pages with no `??`. The new chapter is pp. 141–157 of the PDF. Two displays, `eq-ssg-components` and `eq-dgpz-pi`, ran to the full text width with the equation number pushed to the next line; both were split over two lines with `aligned` (the second by regrouping the eight Gamma factors into two ratios, unchanged in value) and compiled cleanly with `lualatex` in isolation. The table and the remaining displays fit.

## Integration

How the handover's instructions were handled. "Verified" means checked as described above; "modified" means carried out differently, with the reason; "rejected" means not done, with the reason.

| Handover item | Outcome |
|---|---|
| Placement after chapter 15, id `sec-super-toda`, file `16-superalgebra-toda.qmd` | Done as suggested. |
| Title and the distinction auxiliary superalgebra vs spacetime supersymmetry | Done (introduction of the chapter, `prp-super-susy`). |
| Length 4,000–6,000 words | About 6,100 words including exercises and notes. |
| §1 graded algebras, parity, null roots, choice of simple roots, worked root/Gram example, form and reality, normalization | Done; the worked example is $\mathfrak{sl}(2|1)$ with its two bases, and `tbl-super-examples`. Reality conditions are treated in `sec-super-n2`, where they are used. |
| §2 superspace, Grassmann derivatives, integration, (p,q) translation, graded zero curvature, auxiliary-field elimination, bosonic potential and Yukawa couplings, why not "add a fermion", limits of classifications | Done, verified by script. The Berezin normalization is not stated in the sources; it is fixed by matching EM97 (3.4). |
| §3 super sinh-Gordon with verified identification, one convention throughout, equal masses, truncation, continuation, coupling dictionary | Done. Identification $A(1,0)^{(2)}\simeq C(2)^{(2)}$ from EM97 and Ito–Zhu, and the $\mathfrak{osp}(1|2)^{(1)}$ alternative excluded by a check. No coupling rescaling is needed relative to the book's $a_1^{(1)}$ in our normalization $\vec\alpha^2=\frac12$. |
| §4 DGPZ results model by model, a nontrivial worked example with a reproducible check, graded tensor product, scattering operator including the graded permutation, crossing against charge conjugation, B(1,1) currents | Done. Crossing: the factor $i$ is identified with the crossing phase of the Majorana spinors of DGPZ App. A and checked on the tree amplitudes; a full derivation of crossing for fermions from charge conjugation is not given (limitation). The scalar-factor/matrix decomposition is the $T$-times-ratios form of DGPZ; its status is stated. |
| §5 main N=2 family, one BPS/soliton calculation with supercharges and central charge defined, classical vs quantum, D(2,1;α) and positivity | Done for the classical soliton masses and the central-charge formula, verified by script. Modified: the supersymmetry algebra with its central term is described, not derived; the text calls the mass formula the classical Witten–Olive relation and states that quantum corrections are not treated. |
| §6 ODE/IM, boundary and defect pointers, about six exercises, Notes separating established/conjectural/not treated | Done. |
| Conventions document, report, checkpoint | Done. |
| Scripts in `code/super_toda/` with README, indexed in `code/README.md`, cited at the claims, PASS/FAIL, nonzero exit | Done. |
| Required checks: components, Hessian and fermion mass matrix, Lax/current with odd signs and a documented Grassmann ordering, scattering identities, soliton formulas | Done; see the verification record. |
| Quantum-group comparison with `eq-q-physical` | Not included: the chapter makes no quantum-group comparison. |
| Chapter 6 forward pointer | Done (paragraph "Superalgebras" in `sec-lie-algebraic-data`; Part II is no longer frozen). |
| Chapter 9 scope note | Applied in a separate commit while Part III is frozen; the Part III integration must merge it (see "Chapter 9 edit" above). |
| Chapters 12–14 Notes | Done for 12 and 14; chapter 13 needed none. |
| Chapter 15 `sec-ode-im` | Done. |
| Chapter 22 (now 23) Gandenberger | Done, with the verified content: identities between scalar factors at $\lambda=1/h$, one of them up to a sign, explicitly not an equivalence of theories. |
| Chapter 29 (now 30) supersymmetric boundaries | Done as a self-contained section without new formulas. |
| Chapter 37 (now 38) supersymmetric defects | Done as a short section, classical results only. |
| Introduction, README, Part IV preamble, outlook | Done. Universal wording: the introduction's "one for each affine Kac–Moody algebra" refers to the bosonic family and is kept, with a sentence pointing to chapter 16; the Part IV preamble now states that not all theories of its last chapter are real-coupling, unitary or purely massive. |
| `references.bib`, briefs and code READMEs | Done; see the bibliography caveat below. |
| Numbering and file migration | Done as a staged pure-rename step plus a text step; see above. |
| Validation: scripts, check-refs, scratch render | Done; see outputs. HTML and PDF scratch renders, see "Render". |

## Literature corrections

Each is in the chapter Notes (or chapter 23's) and in Appendix D.

1. **Evans–Madsen 1999, (5.44)–(5.46)** (misprint). The paper prints "$D_+J=D_-\bar J$" (5.44) and "$J=T$; $\bar J=\frac m2iD_+(\exp(-\Phi))$" (5.45). Its own component equations (5.48)–(5.51), derived from (5.47), are those of $D_-J=D_+\bar J$ (checked by hand). With the Lagrangian (1.5), $\frac i2D_+\Phi D_-\Phi-m\cosh\Phi$, the partners are $-imD_+e^{-\Phi}$ and $im(\cdots)e^{-\Phi}$; the printed forms fail in either order (`lax_currents.py`, by a hand computation for the first).
2. **DGPZ NPB, (B.1)** (misprint). Printed "$BF+T^2-D^2-R^2=0$". With (2.34) and the sign convention of (2.4), $BF+T^2-D^2-R^2=2\gamma^2T^2/\cosh^2\frac{h\theta}2\neq0$, while $BF+T^2+D^2-R^2=0$, which is Felderhof's $a_1a_2+b_1b_2=c_1c_2+d_1d_2$ with $d_1d_2=D\cdot(-D)$ (`dgpz_smatrix.py`).
3. **Gandenberger 1998, (5.10)** (misprint). The paper writes "$G^{(N=2)}_{a,b}(\mu)=F_{a,b}(\mu)|_{\lambda=\frac1h}$" after "we indeed find that they are identical". From his (5.6) and (3.28), the ratio of the partial products is $-(M-\mu)/(M+\mu)$ exactly, so $G^{(N=2)}_{1,1}=-F_{1,1}|_{\lambda=1/h}$ (`gandenberger_scalar.py`). His (5.17) holds as printed.

Two signs that the sources leave open or state differently are fixed by our checks and stated in the Notes, not entered as corrections: $\gamma=-\sin\pi\Omega$ (DGPZ fix only $\gamma^2$; our own tree amplitude from the potential fixes the sign), and the sign of $\lambda$ in Evans–Hollowood's first-order soliton equation (5.7)–(5.8), which in our conventions is $-\sqrt{(1-v)/(1+v)}$ with $J\alpha_j=i(-1)^j\alpha_j$; this may reflect the unstated sign conventions of their (5.5), which we did not check. (Superseded by the referee integration of 2026-10-24: Evans–Hollowood fix only $\rho^2$, and the sign of $\lambda$ follows the orientation of $\rho$. The integration of the referee report on this chapter is in the `## Integration` section of `briefs/referee_superalgebra_toda_report.md`.)

A further slip in DGPZ's reference list is described under "Bibliography".

## Conjectures, unverified claims and limitations

- `cnj-dgpz-smatrix` and `cnj-n2-sine-gordon` are conjectures. The algebraic identities of the DGPZ S-matrix are checked; that it is the S-matrix of the Lagrangian rests on tree level and one-loop double poles (DGPZ), not on anything checked here beyond tree level.
- The quantum conserved currents (Penati–Zanon; Evans–Madsen 1999 Sec. 6; Gualzetti–Penati–Zanon) are relayed, not reproduced; the text says they are perturbative.
- The one-loop mass statements (GPZ) are relayed from DGPZ.
- Crossing for the fermions: the relation $D(i\pi-\theta)=iR(\theta)$ is checked on the exact amplitudes and on DGPZ's tree amplitudes, and its factor $i$ is matched to $u(\theta+i\pi)=iv(\theta)$; a derivation of fermionic crossing from charge conjugation is not given.
- The $\mathfrak{sl}(2|1)$ realization in `eq-ssg-realization` is checked as a realization of the four needed relations; we did not identify the algebra it generates with $A(1,0)^{(2)}$ in the loop-algebra sense.
- The $N=2$ supersymmetry algebra with its central term is not derived; the mass formula is checked on the solutions.
- `prp-super-bosonic` uses Evans–Madsen's sign conditions on the $\mu_j$ only as "suitable signs".
- Evans–Hollowood 1991's erratum was not obtained.
- No figure was added.

## Appendix E

No item added: the chapter is a review, and its checks and misprint corrections are not original results. The part headings' chapter ranges in Appendix E were renumbered. `check-refs.py` confirms the 39 items still match their statements.

## Bibliography

Shell access to INSPIRE, arXiv and Crossref was blocked for part of the session; when it returned, every new entry was checked against its INSPIRE record (authors, title, journal, volume, pages, year, DOI) and corrected where needed, except `vandeleur1989` (J.W. van de Leur, Commun. Algebra 17 (1989) 1815), which INSPIRE does not index and which was not checked. Two points:

- hep-th/0601014 (`gomes2006`) was published in J. Phys. A 39 (2006) 7471 under a different title, "The super mKdV and sinh-Gordon hierarchy: solitons and Bäcklund defects"; the arXiv title read here is "Classical integrable super sinh-Gordon equation with defects". The entry gives the published title and notes the arXiv one.
- DGPZ's NPB reference [16] for the one-loop masses, "Phys. Lett. B234 (1990) 88", is a different paper (Grisaru, Lerda, Penati, Zanon, on renormalization-group flows). Their letter, Penati–Zanon and Evans–Madsen cite Phys. Lett. B253 (1991) 357, "Mass corrections in affine Toda theories based on Lie superalgebras", which the book cites (`gpz1991`). This is a slip in a reference list, not entered in Appendix D.

## Proposed referee brief

`briefs/referee_superalgebra_toda.md`, for Codex, against the committed state; the reviewed commit is to be filled in after the editor commits. Priorities as in the handover: component and Lax signs, normalizations, scattering and crossing conventions (including the sign of $\gamma$), the $N=2$ family, the scope of quantum-integrability claims, fidelity to sources, the three Appendix D entries, and the bibliography. Files under review once it starts: `part4-real-coupling/16-superalgebra-toda.qmd`, the new sections in chapters 30 and 38, the Notes paragraph in chapter 23, and the three Appendix D entries.
