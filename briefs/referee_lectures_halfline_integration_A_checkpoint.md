# Integration agent A checkpoint: index, Lectures 1 and 2

Scope: items I1–I6, L1-1 to L1-19, L2-1 to L2-10 of `referee_lectures_halfline_combined.md`, and G1–G4 in `lecture-notes/half-line/index.qmd`, `01-wall-and-euclidean-trick.qmd`, `02-sklyanin-and-rigidity.qmd`, `code/lectures/halfline_l1_f1_reflection_crossing.py`, `code/lectures/halfline_l2_f1_sklyanin_rigidity.py` and their two SVGs.

Scratch: `/tmp/claude-1000/-home-gustav-Git-affine-toda/af9e27ec-122e-43f3-a547-2dbab370a65e/scratchpad/integ_A/` (verification scripts `integA_l1_checks.py`, `integA_l2_checks.py`; copied to `briefs/referee_lectures_halfline_scripts/integration_A_*` at the end).

## Log (one line per finished item)

- (progress) Independent verification done for all L1/L2 mathematical items: `integA_l1_checks.py` (49 PASS) and `integA_l2_checks.py` (31 PASS) in scratch. No lecture file edited yet at this point. Next: edit Lecture 1, then Lecture 2, then index, then scripts/figures.
- L1-1 modified (L1 text): repair (a) of P1 2.20 (density T−X, flux T+X), which is the book's convention (`sec-factorization`: (T_{s+1},X_{s-1}) with ∂_-T=∂_+X is spin −s, Q_{-s}=∫(T−X)). Also relabelled spins to the book's (first current is spin −s) and switched Θ→𝒳 with barred parity partner; s=1 normalization Q_{±1}=E±P gives Q̂_1=2(E+𝓑), Σ_1=2𝓑. Evidence: integA_l1_checks.py (sympy identities + explicit s=1 current). Note: the referee's counterexample T=t+x, X=t−x satisfies both forms; T=t−x, X=t+x refutes the printed one.
- L1-2 verified (L1 text): boundary term acts as constant H^{(s)} on asymptotic states; cancels only if the boundary returns to the same state. Evidence: book `sec-half-line-kinematics` display, by hand.
- L1-3 verified (L1 text + Exercise 1.1): equal-coefficient hypothesis stated, γ^{(−s)}=γ^{(s)} assumption stated, general case referred to book `prp-boundary-elastic`; Exercise 1.1(2) now: RHS→r, LHS→0/1/∞, so θ'_max=θ_1 then r=1; stem says single species or equal γ^{(s)}. By hand (+ numeric illustration in integA_l1_checks.py).
- L1-4 verified (L1 text): Š=S_12P=PS_21, K in second slot, term-by-term conversion displayed; P S_12 convention exchanges 12↔21. Evidence: integA_l1_checks.py random-matrix test.
- L1-5 verified (L1 eq + Exercise 1.2 display): Ř_{μ̄,ν}→Ř_{ν,μ̄}. Evidence: domain tracer in integA_l1_checks.py; matches book eq-reflection-twisted.
- L1-6 modified (L1 text + caption): rotation (t_E,−x) as in the book figure; replaced the inconsistent E/p swap by a derivation: the rotation is the continued boost with rapidity −iπ/2, shifts rapidities by +iπ/2; reflected leg → θ̃=iπ/2−θ, incoming leg crossed → −θ̃. Evidence: integA_l1_checks.py (determinants, image of x<0, tanh identity, inconsistency of printed relation), boost formula by hand.
- L1-7 verified (L1 text): pair-only statement qualified (zero-rapidity single particle allowed; GZ (3.42) term ∝ g^c A_c(0)); measure tied to 2πδ normalization (book: dθ/4π). Source read: GZ TeX (3.29)–(3.45), book eq-boundary-state.
- L1-8 verified (L1 Blackboard 2): lower indices S_{cd}^{ba}; derivation rewritten (in/out expansions, ZF at θ1=t, θ2=−t). Evidence: by hand + integA_l1_checks.py (SG coideal K via code/boundary_common.py: corrected residual 7e-16, printed 0.505).
- L1-9 verified (L1 text): product form attributed to Corrigan–Delius (2.7); CDRS credited with the equivalent shifted a_2 form K_1(θ)K_2(θ−iπ)=S_11(2θ). Source read: 9909145 `cu`, 9404108 `\crossing`.
- L1-10 verified (both bootstraps): S_{bā}→S_{ba} (excited) and bulk-fusion bootstrap reordered to the book's K_a S_ab K_b. Evidence: integA_l1_checks.py a_2 tests (S_11 holds 2e-41, S_12 fails 0.12; dressing S_baS_ba unitary+crossing, printed violates).
- L1-11 verified (L1 text): "signals" → "can signal", every pole must be explained; fixed/moving diagnostic replaced by the qualified sentence (Coleman–Thun alternative). By hand vs book `sec-boundary-bound-states`.
- L1-12 verified (eq-lec1-one-point-residue): written as K_ā^b ≃ −(i/2)Γ_ab g/(θ−i(π/2−u/2)) from GZ (3.41) with the Jacobian; same Jacobian gives +(i/2)g² at iπ/2 from GZ (3.45); added the iπ/6 residue check i(3−√3). Evidence: integA_l1_checks.py residues (1/ε scaling), GZ source.
- L1-13 verified (Step 2 rewritten): complete derivation (a) shifted product −ε²(3)(5)/(2), (b) doubling identity (x)|_{2θ}=−(x/2)(6−x/2)^{-1} derived via sinh 2u and cosh z=−i sinh(z+iπ/2), (c) S_LY(2θ)=−(1)(2)/((4)(5)), (d) product ε²(1)(3)/(4) ⇒ ε=1. New label eq-lec1-doubling. Evidence: integA_l1_checks.py.
- L1-14 verified: Γ²=−2√3 with residue computed in the text; "Cardy". Evidence: integA_l1_checks.py.
- L1-15 verified (Steps 2, 4): status of DPTW identifications (TBA/TCSA, book prp-ly-reflection Numerical), b–h fitted numerically, b=0 equal amplitudes ≠ equal boundaries (DPTW conclusion read); shift K_{b+2}=K_{4−b} with the missing cancellation line; third CDD condition shown; fixed the crossing parenthesis ((x)|_{iπ/2−θ} = −(6−x)|_{iπ/2+θ}, not (6−x)^{-1}); "miracle" removed. Evidence: integA_l1_checks.py.
- L1-16 verified (Exercise 1.3(4)): f_x=−(6x)(6−6x), with a hint; heading links both book exercises (exr-boundary-cdd, exr-crossing-forms). Evidence: integA_l1_checks.py.
- L1-19 verified: "inits"→"in its", "therotated"→"the rotated" (passage rewritten).
- G3 (L1): all typed book numbers removed (Chapter/Section/Figure/Proposition/Exercise), exercise title prefixes "Exercise 1.x:" removed. check-refs lecture-notes/half-line PASS.
- Next: L1-17 figure + L1-18 script updates (script still to be edited), then Lecture 2, index.
- L1-17 verified (script figure + SVG regenerated): annotation "b=0: K_b=K_(1) (equal amplitudes); zero of (1−b)^{-1} cancels the moving pole"; added b=−2 marker (zero of (1+b) cancels the iπ/6 pole; K_{−2} = GZ's second minimal solution); y-axis title "Im θ" matching tick labels; panel (a) uses θ̃, t̃_E=0 and the diagonal crossing form S_{aā}(2θ̃)K_ā. Caption rewritten ("candidate pole positions", schematic). Visual check of rendered PNG done.
- L1-18 verified (L1 text + script): script claims now say what is tested (term-by-term conversion on arbitrary matrices + six-vertex example; three crossing forms for Lee–Yang at sample rapidities; doubling identity, S(2θ) and final bootstrap; figure "drawn by"). Script rewritten: 10 checks, all PASS, ~1 s.
- Next: Lecture 2 (L2-1..L2-10), then index (I1–I6), then G checks and final.
- L2-1 verified (L2 intro, :73 vacuum, classification paragraph after eq-lec2-k-a2): "proves the restriction in the defining a_2 representation; book summarizes the rest with its evidence (prp-durham Numerical, BCDR necessity, non-simply-laced continuous nodes)"; ground state "if it exists ... Lecture 3 constructs backgrounds for some patterns". By hand vs book tbl-durham/prp-durham.
- L2-2 verified (L2 text after eq-lec2-dt-tau): K_-=1 justified by vacuum asymptotics for σ=T (a_t(λ)+a_t(1/λ)^T=0 at a vacuum); involution needs the ultralocal bracket and both reflection equations (K_-=1 satisfies its own). Evidence: integA_l2_checks.py (vacuum identity); book prp-sklyanin (b).
- L2-3 verified (L2 text under eq-lec2-k-a1): solves node equations for arbitrary C; invertible for generic λ (λ^4 det K = λ^8 − 4C0C1λ^6 + (4C0²+4C1²−2)λ^4 − 4C0C1λ² + 1); vanishes at C=0, λ=1. Evidence: integA_l2_checks.py (sympy).
- L2-4 modified (Branch I rewritten): propagation is simpler than P1's: with k13=k23=k31=k32=0 from node 1, node 2's spectator equations λk12=0, k21/λ=0 kill k12,k21 for every C_2 (no need to pass through C_2=0); then active equations give C_j k11=0. Mixed patterns excluded explicitly. Evidence: integA_l2_checks.py (sympy + 125-pattern scan).
- L2-5 verified (caption): orange slice (C,+1,−1) vanishes only at C=±1; green at 0,±1; slices sampled at one λ; schematic left panel. Evidence: integA_l2_checks.py. Figure annotation to be adjusted in the script.
- L2-6 verified (Dirichlet paragraph): Cartan condition is the anticommutator {K,H}=0; in sl_{n+1}, n≥2, forces K=0 (explicit choice of traceless H); a_1 keeps the antidiagonal K; scope of book prp-no-dirichlet stated (untwisted, evaluation reps w/o trivial summands, stated σ; field-dependent K not excluded). Evidence: integA_l2_checks.py (nullspace dims sl_2,3,4).
- L2-7 verified (coweight paragraph): added the root-phase argument ((−1)^{α_j·ν} w_j − w̄_j) and narrowed the script claim to "illustrates for one coweight of a_2". Evidence: by hand + integA_l2_checks.py (field-based check).
- L2-8 verified (Exercise 2.1): d_4 vector representation, roots, doublets/spectators named; node-local with C_j/√n_j; determinant 1−C_j²/n_j; C_2²=2 branch. Book finding G5(b) (book `exr-node-local` and ch. 31 "couple the nodes") confirmed: pending for coordinator. Evidence: integA_l2_checks.py (32 patterns: normalized all exist; unnormalized |C|=1 none, min σ 0.039; C_2²=2 all exist).
- L2-9 verified (Exercise 2.2 heading): links both book exercises (exr-sg-k-classical, exr-gz-dictionary); Table 31.2 typed number removed.
- L2-10 verified (Exercise 2.3 rewritten): D flips C_j→t_jC_j at the same λ, t_j=d_ad_b, Πt_j=1, two classes of four; inversion flips s; together one orbit of eight; s=+1 class = uniform + three solitonic. Evidence: integA_l2_checks.py.
- G3 (L2): all typed book numbers and "Exercise 2.x:" prefixes removed; factorial-looking "!" removed from displayed equations (e.g. "1 − C_1^2!").
- Next: L2 script (checks + figure annotation), then index.
- L2-5 (figure part) done: annotation "Neumann C=0 (green slice only)", green legend "equal couplings", panel (a) box now "τ̇=0 ⟸ compatibility" (sufficient, not iff). SVG regenerated, PNG inspected.
- L2 script rewritten: 9 checks (a_1 formula + det polynomial + zero at C=0,λ=1; a_2 8 patterns; 125-pattern rigidity scan + orange slice at C=0; classical RE for both K + K_-=1 left-end RE + vacuum condition [2 lines]; Dirichlet anticommutator dims; d_4 normalization 32 patterns; coweight flips/inversion; soliton-preserving coweight wall), all PASS, ~1.6 s.
- Next: index (I1–I6, G1–G3 there).
- I1 verified (index §1): rotation (x̃,t̃_E)=(t_E,−x) (det +1, x≤0 → t̃_E≥0, as in the book's figure; the printed (−t_E,x) maps x≤0 to negative new time); garbled clause fixed ("spatial translations give conservation of momentum"). Evidence: integA_l1_checks.py.
- I2 verified (index §1): pairs at nonzero rapidity, zero-momentum single particle allowed with a one-point coupling. Source: GZ (3.29) footnote 5, (3.42).
- I3 verified/modified (index lines 9–12, §2, §3, §5, syllabus): rigidity scoped to simply-laced a_n (other types keep continuous nodes); couplings discrete "within Sklyanin's construction with field-independent K"; Dirichlet scope = book prp-no-dirichlet; antisoliton return "for the boundaries considered" (parity); mirror collision behind/at/in front (not universal); backgrounds "explicit regular ... on allowed modulus intervals; believed not proved ground states; not every pattern"; vacuum puzzle restated (product +1 only, from a_3 not all, book tbl-missing-vacua); §5 "examine the proposed dictionary", no WKB, "exchanged exactly" for the proposed amplitudes, conjectural identification; a_n intertwiner dichotomy qualified (a_2 proved, book n≤7 + sketch); syllabus L1 "reflection factors", L5 "conjectural argument", L6 "CDD freedom", "proposed duality", "rotor-like", CT "double pole reduced by a zero", WKB removed; "every matrix multiplication ... without hiding a single step" → "decisive ... by hand"; "zero skipped algebra" → "step by step; where a long computation is taken from the book, the text says so"; "Miracle"/"Reconciliation" headings toned down. Wording consistent with B (mirror, vacua) and C (rotor-like, proposed amplitudes, CDD). By hand vs book labels.
- I4 verified (index §3, §5): duality "proposed by Corrigan in rank one and conjectured by Gandenberger (a_2) and Delius–Gandenberger (a_n)" with @gandenberger1999a added; images: SSW (sine-Gordon) and Delius (affine Toda) [@delius1998a], backgrounds [@bowcock1998; @delius_gandenberger1999]. Sources: P2 §9.1/§9.6 quotes; agent C read Corrigan/Gandenberger/DG sources (C checkpoint).
- I5 verified (index study notes): P1 2.25 replacement, plus "each lecture says what its script checks".
- I6 verified (index line 11): "soliton multiplet is self-conjugate (V≃V̄: the antisoliton lies in the same doublet)"; solitons are not their own antiparticles. By hand.
- G3 (index): "[Ch. NN]" links in the syllabus replaced by chapter-title links; "Part VII" removed from link text and prose. Also coideal convention note (Kolb's right coideals related by a coproduct-reversing automorphism), consistent with agent B's L4 text.
- Checks: check-refs lecture-notes/half-line PASS (7/7); check-refs (book) PASS except a pre-existing WARN at citation.qmd:11 (not my file).
- Final: scratch scripts copied to `briefs/referee_lectures_halfline_scripts/integration_A_l1_checks.py` (51 PASS, 0 FAIL; imports only `code/boundary_common.py`) and `integration_A_l2_checks.py` (27 PASS, 0 FAIL). (The progress line above miscounted them as 49/31.) Lecture scripts: L1 10 PASS in 1.4 s, L2 9 PASS in 1.7 s. Both SVGs have real content changes (kept). check-refs PASS. No commits made.

## Table

| Item | Verdict | How checked | Files changed |
|---|---|---|---|
| I1 | verified | determinants/image of x<0 (integration_A_l1_checks.py); book fig-boundary-state-rotation | index.qmd |
| I2 | verified | GZ TeX (3.29) fn. 5, (3.42) | index.qmd |
| I3 | verified (wording own) | book labels (prp-durham, prp-no-dirichlet, prp-solitonic-vacua, tbl-missing-vacua, cnj-neumann-*); B/C wording | index.qmd |
| I4 | verified | P2 §9.1/§9.6 source quotes; agent C's source reads; bib keys | index.qmd |
| I5 | verified | P1 2.25 text, own audit of the two scripts | index.qmd |
| I6 | verified | by hand (SG doublet) | index.qmd |
| L1-1 | modified (repair (a) + book spin labels, 𝒳 notation) | sympy (integration_A_l1_checks.py); book sec-factorization | 01 |
| L1-2 | verified | book sec-half-line-kinematics display | 01 |
| L1-3 | verified | by hand; numeric illustration | 01 (text + Exercise 1.1) |
| L1-4 | verified | random-matrix identity (scratch + L1 script check 1) | 01, L1 script |
| L1-5 | verified | domain tracer; book eq-reflection-twisted | 01 (eq + exercise) |
| L1-6 | modified (derivation via continued boost −iπ/2 instead of P1's prose) | determinants, tanh identity, boost by hand | 01 (text + caption), L1 script (figure) |
| L1-7 | verified | GZ TeX; book eq-boundary-state (dθ/4π) | 01 |
| L1-8 | verified | ZF by hand; SG coideal K proportionality (scratch + L1 script check 3) | 01, L1 script |
| L1-9 | verified | CDRS `\crossing`, Corrigan–Delius `cu` read | 01 |
| L1-10 | verified (+ bulk bootstrap reordered to K_a S_ab K_b) | a_2 numerics (scratch + L1 script check 4) | 01, L1 script |
| L1-11 | verified | book sec-boundary-bound-states | 01 |
| L1-12 | verified (+ iπ/6 residue check added) | Jacobian by hand; residues by 1/ε scaling (scratch + L1 script check 8) | 01, L1 script |
| L1-13 | verified (complete derivation written) | block identities (scratch + L1 script checks 6–7) | 01, L1 script |
| L1-14 | verified | residue by hand + numerics | 01, L1 script |
| L1-15 | verified (+ fixed crossing parenthesis, third CDD condition shown, shift step filled) | numerics; DPTW conclusion read | 01, L1 script |
| L1-16 | verified | numerics; heading links both book exercises | 01, L1 script |
| L1-17 | verified (+ b=−2 cancellation marked) | numerics; rendered PNG | L1 script, l1 SVG, 01 caption |
| L1-18 | verified | own reading of script | 01, L1 script |
| L1-19 | verified | — | 01 |
| L2-1 | verified | book labels | 02 |
| L2-2 | verified | vacuum identity, K_- RE numerics (L2 script check 4) | 02, L2 script |
| L2-3 | verified | sympy determinant | 02, L2 script |
| L2-4 | modified (shorter propagation: node-2 spectator equations kill k12,k21 for any C_2) | sympy + 125-pattern scan | 02, L2 script |
| L2-5 | verified | numerics; rendered PNG | 02 caption, L2 script, l2 SVG |
| L2-6 | verified | nullspace dims sl_2..sl_4; book prp-no-dirichlet scope | 02, L2 script |
| L2-7 | verified (+ root-phase argument added) | by hand + numerics | 02, L2 script |
| L2-8 | verified | d_4 vector rep, 32 patterns | 02 (Exercise 2.1), L2 script |
| L2-9 | verified | book exercise ids | 02 |
| L2-10 | verified | numerics (flips, classes, inversion) | 02 (Exercise 2.3), L2 script |
| G1 | verified (rows in index, L1, L2) | book labels | index, 01, 02 |
| G2 | verified (Index:69, L1:162/264/266/404, L2:172/238/354) | own script audit | index, 01, 02 |
| G3 | verified | check-refs; grep for typed numbers | index, 01, 02 |
| G4 | verified | rendered PNGs | both scripts, both SVGs |

## Pending

- Book finding G5(b) confirmed: chapter 31 (`sec-durham-classification`, bullet "It is node-local") and `exr-node-local` say that dropping n_j makes the factors "couple the nodes" / "no longer decouple node by node". Independently, the altered equations stay node-local with C_j/√n_j (integration_A_l2_checks.py: only the central coefficient changes; K exists iff (C_2/√2)²=1). Suggested book wording: "The equations remain node-local, but the allowed values of the parameters acquire n_j dependence (C_2^2=2 at the central node of d_4^{(1)})."
- Minor book wording (optional): `sec-boundary-lee-yang` says "Unitarity and crossing-unitarity are linear in K"; unitarity K(θ)K(−θ)=1 is quadratic. The point (both are invariant under K→−K) stands.
- Record note (not a book finding): P1 2.20's counterexample T=t+x, Θ=t−x satisfies both the printed and the corrected equation; T=t−x, Θ=t+x refutes the printed one. The finding itself is correct.
- Nothing left unfinished in my items.
- Script states for `code/lectures/README.md` (coordinator):
  - `halfline_l1_f1_reflection_crossing.py`: 10 checks, ~1.4 s: (1) braid ↔ S/K term-by-term conversion with Š=S_12P on arbitrary matrices (fails with PS_12); (2) six-vertex example in both forms, K in slot 2; (3) crossing-unitarity index order S_cd^ba for the sine-Gordon coideal K (SG reps/R/K built from the coproduct in the script; S_dc^ba fails); (4) a_2 diagonal bootstrap channels (bulk fusion with S_11, excited dressing S_baS_ba unitary/crossing; S_bā variants fail); (5) three diagonal crossing forms for Lee–Yang K_b; (6) Lee–Yang doubling identity, S(2θ)=−(1)(2)/((4)(5)), shifted product; (7) Lee–Yang unitarity and bootstrap sign (ε=+1 holds, −1 fails); (8) residues Γ²=−2√3, g²=−4(2√3−3), Res_{iπ/6}K_(1)=i(3−√3) with square −Γ²g²/4; (9) CDD factor Φ_b, f_x=−(6x)(6−6x), ratio; (10) K_b^(1)=K_{b+2}=K_{4−b}, K_0=K_(1), K_{−2} = GZ second minimal. Figure: panel (a) schematic in θ̃, t̃_E; panel (b) candidate pole positions with b=0 and b=−2 cancellations.
  - `halfline_l2_f1_sklyanin_rigidity.py`: 9 checks, ~1.7 s: (1) a_1 K formula, λ^4 det K polynomial, K=0 at C=0,λ=1; (2) a_2 K for 8 sign patterns; (3) 125-pattern rigidity scan incl. mixed patterns, orange slice at C=0; (4,5) classical RE for a_1 and a_2 K, K_-=1 left-end RE, vacuum condition; (6) Dirichlet anticommutator nullspace (a_1: 2, a_2/a_3: 0); (7) d_4 vector normalization (32 patterns, unnormalized none, C_2²=2 exists); (8) a_2 coweight flips/inversion classes; (9) soliton-preserving wall with K_+=e^{iπν·H}. Figure: panel (a) "⟸", panel (b) slices labelled, "Neumann C=0 (green slice only)".
