# Ch. 43 drafting report: "Breather and particle transmission amplitudes"

Files written (nothing else touched, nothing committed):
- `part8-defects/43-particle-transmission.qmd` (266 lines, replaces the skeleton; all skeleton ids kept)
- `code/defect_particle_transmission.py` (new; imports only from `boundary_common`, `defect_common`)
- `briefs/draft_ch43_checkpoint.md`, `briefs/draft_ch43_report.md`
- Scratch (not in repo): `/tmp/claude-1000/-home-gustav-Git-affine-toda--claude-worktrees-integrable-defects-framework-a54a84/cfddf2dd-7e99-4c1c-b8f5-c799ea33aced/scratchpad/ch43/` (`breather.py`, `breather_all.py`, `doubled.py`, `b6b.py`, `b6c.py`)

Notation: ω = 4π/β² − 1 is always the coupling (as in ch. 41); ch. 39's root of unity is written e^{2πi/h}; γ_a = πa/h. In problem 11 the ch. 40 twist is written ω_m = q^{−m/2} and flagged as "a twist, not the coupling".

## 1. Chapter summary

**Main finding.** The lowest-breather transmission factors computed from the ch. 41 soliton matrices do not depend on the coupling, for every species. They equal the classical linear factors T^cl_a of ch. 39. So continuing to real coupling is trivial, and the only quantum content of particle transmission is the renormalized defect rapidity θ_d(σ, β).

**Sections** (all skeleton ids kept):
- `sec-breather-transmission`:
  - sine-Gordon recap;
  - the a_r breathers B^(a)_1 and their constituents. ū = π/2 − π/(hω); species h−a sits at θ−iū and species a at θ+iū. Built from h vector solitons, B^(a)_1 is the Λ^h V string with its two halves pulled together by 2π/(hω);
  - `prp-lightest-breather-universal` and what it means;
  - the block form 𝕋^cl_a = −e^{iπa/h}/(h−a) and 𝕋^cl_{+,a} = e^{iπa/h}/(a);
  - class-p defects; higher breathers.
- `sec-an-particle-transmission`:
  - `eq-dpt-particle`; the axioms;
  - Dorey fusings of a_r, of both kinds (c = a+b and c = a+b−h), with their ū's;
  - `prp-an-particle-transmission`;
  - unitarity, CDD ambiguities, no bound states.
- `sec-transmission-continuation`:
  - four assumptions: minimality, vacuum, parameter, integrability;
  - `prp-dpt-vacuum`, `cnj-transmission-continuation`, and the evidence for it;
  - the one-loop prediction for a_r;
  - comparison with reflection (ch. 35).
- `sec-transmission-limits`:
  - the classical limit;
  - one loop for sinh-Gordon (Bajnok–Simon Sec. 4.2);
  - the free limit;
  - weak–strong duality (the defect family maps to itself; conditional remark);
  - θ_d → ±∞, the transparent line and the π-line.
- `sec-defect-spectrum`: table `tbl-dpt-spectrum` covering Lee–Yang, sinh-Gordon, SG types I and II, a_r solitons, a₂^(2), and a_r particles for real and imaginary θ_d; a reading of the table; defects in front of boundaries.
- `sec-defect-finite-volume` (brief): the defect TBA (BS), DTCSA (BHW), Kormos–Runkel–Watts, and Jiang (abstract fetched from arXiv:1703.03562 and paraphrased).
- `sec-defect-open-problems`: opens with `prp-dpt-doubled`, the test of the ch. 42 a_r expectation, then 15 numbered problems.

**Results and status labels:**

| Id | Label | Content |
|---|---|---|
| `prp-lightest-breather-universal` | [Numerical] for 2 ≤ r ≤ 4 | B^(a)_1 crosses as T^cl_a(θ−θ_d) times the identity on 𝒲. Independent of ω, of the defect label and of the offset s (classes p = 0, 1, 2 for a_3). Ordering: V̄⊗V gives T^cl_1, V⊗V̄ gives T^cl_r; no singlet with the shifts reversed (a = 1, r). cz2007 (7.5) = (T^cl_2, T^cl_1) in book labels. a_5: separate run, all five species, not in the script. |
| `prp-an-particle-transmission` | [Numerical] (as solutions of the axioms) | For h = 2..8, both types: unitarity and crossing, bootstrap at every fusing, phases and Hermitian analyticity, no physical-strip poles, equal to the linearized classical factors. The Lagrangian identification is deferred to the conjecture (Penati–Polvara caveat). |
| `prp-dpt-vacuum` | [Theorem] | Real coupling, σ > 0: B_d ≥ 0 with equality only at φ = ψ = 0 (convexity, Σα_j = 0). σ < 0: unbounded below. |
| `cnj-transmission-continuation` | [Conjecture] | If the real-coupling type-I defect with σ > 0 is quantum integrable, its particles transmit with `eq-dpt-particle` (𝕋⁺_a for I₊), with no reflection, at a renormalized θ_d(σ, β) → −ln σ. |
| `prp-dpt-doubled` | [Numerical] | The ch. 42 a_r expectation holds for type I₋ (a_2, a_3): Δ(U_q(𝔟⁻)) and (1⊗π)Δ(U_q(𝔟⁺)), with ψ(e_j) = e_{j−1}. Every other twist fails. |

**Exercises (8):** `exr-dpt-telescope`, `exr-dpt-cz`, `exr-dpt-vacuum`, `exr-dpt-bs`, `exr-dpt-limits`, `exr-dpt-imaginary`, `exr-dpt-cdd`, `exr-dpt-doubled`.

**New ids** (no clashes; `dpt-` was unused):
- `eq-dpt-particle`, `prp-dpt-vacuum`, `prp-dpt-doubled`, `tbl-dpt-spectrum`;
- the 8 `exr-dpt-*` ids above.

A python check confirmed that every `@`-reference in the chapter resolves to an existing id and every citation key is in `references.bib`.

## 2. Script output and runtime

`python3 code/defect_particle_transmission.py` (from `code/`): **8 PASS, 0 FAIL**: 62 s wall time when run alone (final run, after the last edits), 120 s on a first run under load from parallel scratch jobs; check (f) takes most of it.

```
PASS (a) sine-Gordon: ... = T^cl_1 at h = 2, for lambda = 1.37, 2.6 [2e-13; proportionality 8e-16]
PASS (b1) ... a_2, a_3, a_4, omega = 1.37, generic offset [max dev 4e-09, proportionality 1e-15]; with the shifts reversed there is no singlet [smallest rel. sv 0.09]; cz2007 (7.5) = (T^cl_2, T^cl_1) in book labels [8e-16]
PASS (b2) ... h = 3, 4, 5, all a, omega = 0.6, 1.37, 3.1, real and complex theta [dev 2e-08, proportionality 3e-14]; the same for the classes p = 0, 1, 2 (offset s = lambda_p) and other defect labels (a_3) [2e-09]
PASS (c) ... h = 2..8: [1e-15], Hermitian analyticity [8e-16], bootstrap at every fusing [2e-15] (fusing angles are poles of eq-an-smatrix, masses close); CDD S_ae(theta - theta_0) [2e-14]; no pole in the physical strip for real theta_d
PASS (d) classical limit: linearized defect conditions give R = 0 and T = T^cl_k (I_-), T^cl_+,k (I_+), r = 2, 3, 4 [3e-15]; theta_d -> +-inf limits [2e-16]
PASS (e) sinh-Gordon: BS classical factor and (16) = T^cl with mu = -theta_d^cl [3e-16]; BS one-loop = O(b^2) term of T^cl(theta - (1 - B/2) theta_d^cl) [remainder/b^2 = 9e-08 at b^2 = 1e-4]
PASS (f) ch. 42 expectation for a_r, type I_-: Delta(U_q(b^-)) and (1 (x) pi) Delta(U_q(b^+)) hold [6e-15]; untwisted and pi^-1-twisted e, and pi^+-1-twisted f, fail [min 0.67]
PASS (g) real-coupling vacuum: B_d/|fields|^2 >= 2.34e-01 on random samples, B_d(0, 0) = 0; sigma < 0 unbounded below
```

## 3. Verification record

**By script** (`defect_particle_transmission.py`):
- (a) the SG B_1 re-check;
- (b1) a_2, a_3, a_4 breathers from the closed forms `eq-dt-tv` and `eq-dt-tvbar` with g and ḡ, on the ch. 41 windows. This includes the singlet orderings and cz2007 (7.5);
- (b2) all species B^(a)_1 for h = 3, 4, 5 and ω = 0.6, 1.37, 3.1:
  - the singlet is unique in the weight-zero subspace;
  - exact dict propagation through `eq-dt-tv`, with no window;
  - label and class independence for a_3;
- (c) the particle-factor axioms for h ≤ 8:
  - fusing angles located as poles of the book's `eq-an-smatrix`, and mass closure;
  - the CDD example S_ae(θ−θ_0);
  - the pole positions;
- (d) independent linearization of the defect conditions with `defect_common.defect_potential`/`dmatrix` (complex-step Hessian), r = 2, 3, 4, both types. Also D = 1 + i cot γ_k on the π-eigenvector of eigenvalue e^{2πik/h}, and the θ_d → ±∞ limits;
- (e) Bajnok–Simon's classical formula, their (16), and their one-loop correction against 𝕋^cl(θ − (1 − B/2)θ_d^cl);
- (f) the doubled-boundary test: least squares for θ-independent Y and κ, with V and V̄ jointly; a_2 all nodes, 3 rapidities; a_3 nodes 0, 1, 2 rapidities;
- (g) vacuum positivity.

**By scratch run (not in the repo script):**
- a_5 (h = 6): all five species at ω = 1.37, θ = 0.4 + 0.1i give T^cl_a to ≤ 5e-9. Runtime about 10 min, which is why it is not in the script.
- (f) for all four nodes of a_3: same pattern (holds 3e-15, fails 0.67–0.93).

**By hand:**
- the bootstrap telescoping identity for both kinds of fusing;
- the block forms −e^{iπa/h}/(h−a) and e^{iπa/h}/(a) (also checked numerically in scratch);
- Hermitian analyticity;
- the proof of `prp-dpt-vacuum`;
- ∂_μ𝕋^cl = −1/(1 − i sinh(θ+μ));
- the duality remark θ'^cl = (2−B)/B · θ_d^cl;
- the pole window for imaginary θ_d (`exr-dpt-imaginary`);
- the a_r TBA driving terms for a and ā are complex conjugate;
- the derivation of the two doubled-boundary equations from Δ_A (`eq-dbd-coproduct`) and the weights of Y and κ;
- the string picture of B^(a)_1: the gap between the halves is 2π/h − 2π/(hω).

**Unchecked:**
- **Breather factors not computed:** the type-I₊ breather factors (no I₊ scalar factor exists); breathers of class-p defects through the scalar g(θ − 2πip/h) (the scalar of a class-p defect is not derived; only the offset independence is checked); the higher-breather product formula for a_r.
- **Doubled-boundary module structure:** the module relations of the second structure ρ̃ for a_r. The e/π solution has free constants, and it was not checked that some choice of them gives a U_q(𝔟⁺) module.
- **Not attempted:** a proof of `prp-lightest-breather-universal` for all r.
- **Quoted only (from sources or earlier chapters):** BS defect TBA and ε_def (from ch. 37); BHW DTCSA; Kormos–Runkel–Watts; Jiang (abstract only); cz2010 (3.19); cz2011 (3.14)–(3.15); the a₂^(2) unstable windows; KL (5.30)–(5.37); PP's anomaly statements; the class-p pole positions for general a (by hand in ch. 41).
- **Background statement:** whether a_r particles are exactly the continued B^(a)_1 (the breather test, ch. 21).

## 4. Corrections and clarifications to the sources (in the Notes)

- **cz2007 (7.3)–(7.5).** Their defect matrix is the transpose of the book's, so:
  - their T^{b1}_1 is the book's particle 2 (= r) and their T^{b1}_2 is particle 1;
  - their rapidity assignment (their species 2 at θ − iū) is the book's V(θ−iū)⊗V̄(θ+iū);
  - their observation that the factors "appear to depend on β only via η" holds exactly, and for every species of a_2, a_3 and a_4.
- **Clarification of cz2007 Sec. 7.** η tends to the classical parameter as β → 0, but for sinh-Gordon theory it differs at O(β²) (BS (23)).
- **Clarification of BS Sec. 4.2.** Their one-loop factor is exactly the O(b²) term of the exact factor with θ_d = (1 − B/2)θ_d^cl, B/2 = b²/(8π + b²).
- **Book-internal, ch. 42 `sec-jump-coideal`.** The "rotation-twisted" Borel symmetry is (1⊗π)Δ(U_q(𝔟⁺)) with π(e_j) = e_{j−1}, the direction of πα_j = α_{j−1}. The opposite rotation fails.

## 5. Ch. 43 obligations from the status file

| From | Obligation | Where addressed |
|---|---|---|
| ch. 37 | `prp-lightest-breather-universal` | result: `prp-lightest-breather-universal` [Numerical] r ≤ 4 (a_5 in text); problem 9 (proof) |
| ch. 37 | Penati–Polvara vs Bajnok–Simon | problem 1; `sec-transmission-continuation` (assumption 4, one-loop prediction) |
| ch. 37 | SG continuation of `cnj-defect-parameter` | problem 3; `sec-transmission-continuation` "Comparison with reflection" |
| ch. 37 | stable defects at imaginary θ_d | `tbl-dpt-spectrum`, `exr-dpt-imaginary`, problem 10 |
| ch. 38 | e₆–e₈ defects | problem 5 |
| ch. 38 | `cnj-momentum-integrability` | problem 4 |
| ch. 38 | time-like Liouville integrability for a_r | problem 4 |
| ch. 38 | quantum fate of the second chirality | problem 2 |
| ch. 40 | proof of `cnj-defect-tq` for all N | problem 11; prior art (Rossi–Weston, BLZ, Weston §2.4.1) stated as known, per referee items 10 and 17 |
| ch. 40 | meaning of ω = q^{−m/2} | problem 11 |
| ch. 40 | KL same-chirality terms beyond first order | problem 2 |
| ch. 40 | q̂ vs q | problem 11 |
| ch. 40 | convergence of Runkel's defect operators | problem 11 |
| ch. 39 | I₋I₊ pair = soliton S-matrix × π-line? type II? | problem 6 |
| ch. 39 | lightest-breather factors = T^cl for all a_r | result (r ≤ 4, a_5 in text); problem 9 |
| ch. 41 | bootstrap the lightest breather (cz2007 §7) | result: `prp-lightest-breather-universal` items 2, 4 |
| ch. 41 | a₂^(2) breather unchecked | problem 5 (quoted cz2011 (3.14)–(3.15)) |
| ch. 41 | (h−1)-dimensional image of the residue | problem 10 |
| ch. 41 | candidate stable states at θ_d = 0 | `tbl-dpt-spectrum`; problem 10 |
| ch. 41 | I₊ species-1 closed form and class shift | problem 7 |
| ch. 41 | fused I₋I₊ pair vs soliton S-matrix | problem 6 |
| ch. 41 | remaining a₂^(2) work | problem 5 |
| ch. 41 | Robertson's fused matrices | problem 6 |
| ch. 41 | Penati–Polvara vs axiomatic matrices | problem 1 |
| ch. 41 | q-sign for h odd | problem 8 |
| ch. 41 | hypothesis (ii) for r ≥ 3 | problem 7 |
| ch. 41 | U_q(a_r^(1)) universal R on 𝒲 | problem 7 |
| ch. 41/42 | a_r doubled-boundary expectation, testable with the ch. 41 matrices | TESTED: `prp-dpt-doubled` (holds, type I₋); remaining parts in problem 12 |
| ch. 42 | doubled boundaries of type-II defects and cz2012 (2.19) | problem 12 |
| ch. 42 | quantizing fused a_n dynamical boundaries, involution | problem 13 |
| ch. 42 | non-invertible lines (doubled) | problem 12 |
| ch. 42 | Borel-leg tensor K-matrices as boundary Q-operators | problem 13 |
| ch. 42 | affine diagonal-type QSP | problem 12 |
| ch. 42 | algebra of the two Borel halves | problem 12 |
| ch. 42 | boundary crossing-unitarity of K_f | problem 13 |
| brief | anomalies and quantum integrability | problem 1 |
| brief | type I beyond a_r | problem 5 |
| brief | Krein/Θ with defects | problem 14 |
| brief | real-coupling vacua | `prp-dpt-vacuum`; problem 14 |
| brief | finite volume | `sec-defect-finite-volume`; problem 15 |
| ch. 35 | Part VII problem 10 (defects) | intro to the problem list: "is the subject of this part" (the global rewrite of ch. 35 problem 10 is still pending, not done here) |

## 6. Requests

**Edits elsewhere** (for the integration agent):
1. **Ch. 42, `sec-jump-coideal`, "Beyond sine-Gordon".** Replace "We have not tested this; it needs the a_r transmission matrices of @sec-soliton-transmission." with a pointer: "@prp-dpt-doubled confirms it for type I₋ of a_2 and a_3, with π(e_j) = e_{j−1}; the module relations of the second structure are not checked."
2. **Ch. 37, `sec-sg-breather-transmission`.** "@sec-breather-transmission asks whether ..." can become "it does, for every species of a_2–a_4 (@prp-lightest-breather-universal)".
3. **Ch. 41, end of `sec-transmission-semiclassics` ("Lightest breathers").** Optionally state the result: "they are, for every species (@prp-lightest-breather-universal)".
4. **Ch. 40, `tbl-boundary-defect-algebra`, third column.** It could mention the a_r case via `prp-dpt-doubled`.
5. **Ch. 35, problem 10.** Still on the global pending list; it could now point to `sec-defect-open-problems`.
6. **Status file.** Add the ch. 43 row (`defect_particle_transmission.py`, 8 PASS, 62 s) and the verification record of §3.
7. **`code/README.md`.** Add a row for `defect_particle_transmission.py` (ch. 43, 8 checks, about 1 min). No changes to `defect_common.py`; all helpers are local to the script.

**Missing bib keys:** none needed; every cited key exists.
- Jiang 2017 is `jiang2017` (arXiv:1703.03562, verified).
- For the a_n breather–soliton amplitude, problem 6 uses the existing `gandenberger1998` (hep-th/9703158), as ch. 35 does.

**Conventions-brief additions (proposed §11, ch. 43):**
- **Breathers.** The lowest breather B^(a)_1 is the U_q singlet with species h−a at θ − iū and species a at θ + iū, ū = π/2 − π/(hω). It exists only in this order (checked for a = 1, r). V̄⊗V gives particle 1 and V⊗V̄ gives particle r. CZ's labels are swapped (transposed D).
- **Particle factors.** 𝕋_a = (e^{−D} + e^{iπa/h})/(e^{−D} + e^{−iπa/h}) = −e^{iπa/h}/(h−a) and 𝕋⁺_a = e^{iπa/h}/(a), in blocks (`eq-blocks`) with D = θ − θ_d. 𝕋̄_a(θ) = 𝕋_ā(iπ − θ). Hermitian analyticity holds for real θ_d.
- **Transmission bootstrap (no S-factor).** For a+b < h: c = a+b, ū = (πb/h, πa/h). For a+b > h: c = a+b−h, ū = (π(h−b)/h, π(h−a)/h).
- **Doubled-boundary equations (from Δ_A).**
  - e-type: 𝕋[π(e_j)⊗1 + π(k_j)⊗Y] = [Y⊗1 + κ⊗π(e_{ψj})]𝕋.
  - f-type: 𝕋[π(f_j)⊗κ + 1⊗Y] = [Y⊗π(k_{ψj})^{-1} + 1⊗π(f_{ψj})]𝕋.
  - For a_r type I₋, the twist ψ = π acts by e_j ↦ e_{j−1}.
- **Real-coupling vacuum.** B_d ≥ 0 for σ > 0; unbounded below for σ < 0.

**Flags for the referee:**
- The duality remark and the a_r guess θ_d = (1 − B/2)θ_d^cl are explicitly conditional.
- "Breathers equal T^cl for every r" is stated as an expectation, labelled [Numerical] only for r ≤ 4.
- The CDD paragraph's claim, that a soliton CDD factor generically spoils the coupling independence, is argued, not computed.
