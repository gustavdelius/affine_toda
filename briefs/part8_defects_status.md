# Part VIII (defects): drafting status, verification records, pending edits

Working record for Part VIII, kept in the repository because the session scratchpad is ephemeral. Each drafting agent's verification record is copied here when its chapter is integrated. These lists become the "Already verified" and "To be checked" sections of the referee briefs `referee_defects_*.md`.

## Chapter status

| Ch | Drafted | Script(s) | Checks | Refereed |
|---|---|---|---|---|
| 36 | 2026-10-11 | `defect_axioms.py` | 20 PASS (2 s) | Gemini (light) + Codex deep, integrated 2026-10-12 |
| 37 | 2026-10-11 | `defect_classical_sg.py`, `defect_sg_transmission.py` | 27 + 24 PASS (18 s, 103 s) | Gemini (light) + Codex deep, integrated 2026-10-12 |
| 38 | 2026-10-12 | `defect_classical_an.py`, `defect_lax_rmatrix.py` | 24 + 6 PASS (114 s, 2 s) | Codex deep, integrated 2026-10-12 |
| 39 | 2026-10-12 | `defect_solitons_an.py` | 9 PASS (10 s) | no |
| 40 | 2026-10-12 | `defect_universal_r.py`, `defect_transfer.py` | 16 + 9 PASS (92 s, 29 s) | Codex deep, integrated 2026-10-12 |
| 41 | 2026-10-12 | `defect_an_transmission.py` | 21 PASS (100–176 s) | no |
| 42 | 2026-10-12 | `defect_doubling.py` | 17 PASS (30 s) | no |
| 43 | 2026-10-12 | `defect_particle_transmission.py` | 8 PASS (62 s) | no |

## Chapter 36 (`sec-defect-scattering`)

**By script** (`defect_axioms.py`):
- the three triangle forms are equivalent, with positive and negative controls;
- the no-go reductions (symbolic ZF engine), and S² = 1;
- the DMS Ising amplitudes;
- all Lee–Yang statements of `prp-ly-transmission` and the text: energies, bound state, pure-phase lines at sample points, fusion onto the boundary;
- the sinh-Gordon CDD statements of `exr-defect-cdd`;
- `prp-doubling` items 1–3 (item 1 for the sine-Gordon Ř with 𝒲 = spin 1; items 2–3 for Lee–Yang);
- `prp-twist-defect`, including the seam point;
- the Copetti–Córdova–Komatsu crossing rules and the relating factor.

**By hand (drafting agent):**
- proofs of `prp-rt-nogo`, `prp-twist-defect` and `prp-doubling` item 3;
- unitarity with reflection from the ZF algebra;
- Hermitian analyticity ⇔ Re b ≡ 3 (mod 6), i.e. Re b = ±3 modulo 12 (rechecked by the Codex deep pass);
- the kink-rescaling reading of the CCK factor;
- Bajnok–Holló–Watts μ̂ + μ̄̂ = h(b) of ch. 28;
- the signs and arguments of the doubled K; 𝕋 depending on θ−θ_d and 𝕋̄ on θ+θ_d.

**Independently checked by the Codex deep referee pass (2026-10-12; scripts in `briefs/referee_defects_36_37_deep_scripts/`):**
- `prp-doubling` item 3 with the non-diagonal sine-Gordon 𝕋 and its scalar (full doubled K, |n| ≤ 5), with a negative control;
- `prp-rt-nogo`: both reductions, given label-preserving R, 𝕋 and independent out-states (hypotheses now stated); the direct-sum counterexample to cancelling 𝕋 in CAFG's argument;
- the CCK (18)/(21) crossing rules on all 14 M₅ faces, YBE and unitarity of both forms, the direction of the factor, and fusion-vertex covariance of the kink rephasing;
- the Lee–Yang analyticity condition, the |1⁺⟩, |2⁺⟩ bootstrap and Coleman–Thun pole orders at b = 1.6, the TBA driving term and the Lüscher sign; the crossed-channel factor at b = −3 is positive and < 1, not a phase.

**Unchecked:**
- CAFG's non-diagonal argument beyond the cancellation caveat (their (22)–(24));
- the BHW TCSA conclusions;
- Mintchev–Ragoucy–Sorba (cited generally);
- Jiang and Kormos–Runkel–Watts (cited by abstract);
- Runkel's T-system statements;
- "no interacting example known" for the CAFG loophole (two natural sine-Gordon attempts fail; no general argument);
- that the defect operator, as a trace over 𝒲, commutes with the transfer matrix (left to `defect_transfer.py`, ch. 40);
- the remarks on unstable states (the sine-Gordon resonance reading now restricted to λ > 1);
- the doubled theory's lack of parity invariance (parity combined with copy exchange not examined; Codex 'unverified' list);
- the general equivalences of pure transmission, momentum conservation, deformability and topological lines (Codex 'unverified' list).

**Book-internal corrections (2026-10-12):** the TBA reliability statement (BHW's phase condition concerns real-rapidity 𝕋 at b = −3 + iα, not the crossed-channel factor); the CCK paragraph no longer claims a comparison with amplitudes printed in ch. 27; the Lee–Yang "fundamental range" is qualified.

**Corrections to sources recorded in the Notes:**
- DMS (2.12) needs no parity, and their (2.10) is empty only if S_ab = S_ba.
- The DMS χ branch.
- CAFG (32) omits a factor S_ba.
- CAFG (22)–(23): cancelling 𝕋 needs nondegeneracy (direct sum of reflecting and transmitting sectors) — added 2026-10-12.
- Book-internal: the twist of `prp-bloch-sectors` acts at a common seam.

**Pending edits elsewhere (integration phase):**
- `textbook_code/README.md` row for `defect_axioms.py`.
- `part8-foundations/40-scattering.qmd`: pointer from `prp-exchange-twist` and `prp-bloch-sectors` to `prp-twist-defect`, plus one sentence on the common seam.
- `part7-half-line/28-boundary-scattering.qmd` `sec-boundary-lee-yang`: pointer to `sec-defect-lee-yang` (K_b = 𝕋̄_b K_𝟙 𝕋_b).
- `part6-soliton-smatrices/27-rsos.qmd`: pointer to the modified crossing in `sec-topological-lines`.
- Conventions brief: record the ZF ordering, the 𝕋/𝕋̄ index rule (𝕋̄ parametrized by θ > 0, depending on θ+θ_d), `eq-defect-crossing`, 𝕃 = P𝕋 with R = PŘ, and that the reflected copy scatters with PŠP.

**Obligations for later chapters:**
- Ch. 40 must deliver `prp-defect-twist-limit`, which ch. 36 cites.
- Ch. 42 must decide the symmetry of the doubled boundary; ch. 36 leaves it open, as the pre-registered `cnj-twisted-diagonal` requires.

## Chapter 37 (`sec-defect-sine-gordon`)

**By script.**
- `defect_classical_sg.py` (27 checks):
  - Bäcklund compatibility (sympy); vacuum → soliton; Bianchi → two-soliton with A₁₁.
  - The defect conditions are the Bäcklund transformation with a = −σ, and z is the Bianchi coefficient.
  - I₊ = I₋ ∘ (ψ → −ψ) for a₁, with 𝕋^cl_{I₊} = −𝕋^cl.
  - 𝕋̄^cl = 𝕋^cl(iπ−θ); there is no parity.
  - Absorption and emission (θ_d of both signs); two solitons through the defect.
  - Odd-sector delay and the A₁₁ identities; the M/2 mass shell; 𝕋_odd² = f_{1/2}.
  - Θ holds iff σ is real.
  - Dirichlet fusion gives Ghoshal–Zamolodchikov (GZ) with C_{0,1} = cosh(θ_d ± iβφ_D/2).
  - BCZ 2003 (2.6)–(2.9).
- `defect_sg_transmission.py` (24 checks):
  - The scalar factor Φ: unitarity, crossing, its poles; the BCZ (3.47) phase drift.
  - 𝕋 is unitary at both ±q.
  - **Crossing with C = σ_x holds at q_phys and fails at −q_phys; at −q it is restored only with C′ = [[0,1],[−1,0]].**
  - The resonance at θ_d − iπ/(2λ) (rank-one residue, mass shell, Im E < 0); the 𝕋̄ pole at −θ_d − iπ/(2λ).
  - The odd sector by bootstrap (the s = 1 intertwiner at θ_d − iπ; not unitary).
  - B₁ = 𝕋^cl; B₂ and B_n agree with BCZ (5.10), (5.12)–(5.13).
  - Semiclassical approach of the phase to the classical delay.
  - S₀ Gamma product; S_{sB₁} = the defect pair; Bajnok–Simon fusion identity; BS (17) ⇔ 1 − B/2.
  - Windows: |n| ≤ 5, equations imposed on |n| ≤ 3.

**By hand (drafting agent).**
- Compatibility and permutability proofs; dE/dt and dP/dt.
- The crossing condition ΦΦ′(1+Y²) = q^{−s} for general s (only s = 0 is checked numerically); q^{−s/2}Φ for s ≠ 0.
- No poles of Φ in the physical strip (Codex deep pass: none in 0 < Im θ < π).
- The rational form of 𝕋^cl.
- The momentum criterion forcing exponential potentials (F‴/F′ = G‴/G′); the boost argument.
- The KL dictionary ν ↔ −λθ_d.
- Penati–Polvara normalization γ/2 = b²/8π for real sinh-Gordon coupling b, which continues to −β_SG²/8π (corrected 2026-10-12; source text 1902.10690 rechecked: γ ≡ β²/2π, action 1/(2β²), V = e^{√2φ}+e^{−√2φ}).
- The classical limit of Bajnok–Simon's ε_def; E_D as the continuation of `eq-bsg-dirichlet`.
- 𝕋̄𝕋 = (cosh θ_d + i sinh θ)/(cosh θ_d − i sinh θ).
- The stable-pole window for imaginary θ_d.

**Unchecked.**
- **Type II, quoted only:** all CZ 2009b/2010 type-II statements. Their resonance at η_p + iπ/2γ is a misprint for η_p − iπ/2γ (Gemini pass, 2026-10-11).
- **Attributions to sources:**
  - CZ 2004 gates.
  - Weston's identification, including his Borel subalgebra in book conventions.
  - BCZ Sec. 6 (moving defects) and (3.52)–(3.56).
  - Bajnok–Simon (18), (23); (17) and its scale (14) were rederived in the Codex deep pass.
  - Penati–Polvara's anomaly statements.
  - KL (5.30)–(5.37).
  - The necessity direction of the BCZ classification.
- **Unverified checks:** Bianchi permutability for a general seed; the interpretation of the 𝕋_{B₂} poles.
- **Citations:** `habibullin_kundu2008` is cited by title only.

**Corrections to sources recorded in the Notes.**
- BCZ 2005 (3.42), (3.47) need an endpoint normalization; the paired product r of (3.44) diverges like K^{(γ+1)/4} on its own, and only the common-cutoff ratio r/r̄ converges (corrected 2026-10-12).
- BCZ 2005 Sec. 2 two-defect remark: a directional clarification — correct as a temporal comparison for the overtaken soliton; the left-to-right spatial factors are inverse (reworded 2026-10-12).
- Bajnok–Simon B_BS = B/2; their (17) ⇔ θ_d = (1 − B/2)θ_d^cl.
- BCZ 2005 (2.2): +δ(x).

**Book-internal correction applied to ch. 36 (2026-10-11).** The excited defect is the image of the residue in 𝒲⊗V(iu), with transmission Ř₂₃𝕋₁₂ (the particle sits to the right). Ch. 36 had "a subspace of V(iu)⊗𝒲".

**Pending edits elsewhere.**
- `textbook_code/README.md` rows: `defect_classical_sg.py` (ch. 36–37, 27 checks, ~20 s) and `defect_sg_transmission.py` (ch. 37/40, 24 checks, ~100 s); the new helpers in `defect_common.py`.
- Pointers:
  - ch. 3 `sec-sg-classical` → `sec-backlund`;
  - ch. 10 `sec-multi-solitons` → `prp-bianchi`;
  - ch. 29 `sec-bsg-classical` and `cnj-shg-e` → `sec-sg-breather-transmission`, `cnj-defect-parameter`.
- Optional bib: Rogers–Schief, *Bäcklund and Darboux Transformations* (2002); Habibullin, sine-Gordon on the semi-axis. Confirm what `habibullin_kundu2008` contains.

**Obligations for later chapters.**
- **Ch. 38:**
  - dictionary rows for `tbl-defect-dictionary`;
  - the cz2009a momentum argument in `sec-defect-momentum`;
  - "delay = Bianchi coefficient" for a_r in `sec-defect-backlund-an`;
  - type-II conditions and delay in book conventions;
  - the a_r version of Θ iff σ real.
- **Ch. 39:** the absorption/emission direction rule, and the odd-sector σ → −σ correspondence for a_r.
- **Ch. 40:**
  - reconcile Weston's Borel subalgebra with I₋ ↔ b⁻;
  - relate I₊ = I₋ ∘ (ψ → −ψ) to b⁺ via BCZ's (a,b) labels (3.57)–(3.59);
  - treat the s ≠ 0 scalar.
- **Ch. 41:**
  - reuse the C vs C′ crossing test and the Wigner semiclassical test;
  - reconcile the sign of the CZ 2010 (3.17) pole.
- **Ch. 42:** `prp-defect-fusion-k` with the sinh-Gordon identity and the cz2012 soliton version.
- **Ch. 43:**
  - `prp-lightest-breather-universal`;
  - open problems: Penati–Polvara vs Bajnok–Simon; the sine-Gordon continuation of `cnj-defect-parameter`; stable defects at imaginary θ_d.

## Chapter 38 (`sec-defect-conditions`)

**By script.**
- `defect_classical_an.py` (24 checks):
  - **Momentum-conservation scan:** 23 affine diagrams of rank ≤ 4 (untwisted and twisted) with 78 automorphisms. Only π^{±1} of a_r⁽¹⁾ pass. The failures: 64 automorphisms fix a node, 5 leave 1−P singular, and π^{±2} of a₄⁽¹⁾ fails by cross terms.
  - **Further checks:**
    - cz2009a's gradient form of the momentum condition;
    - the Bäcklund transformation reproduces both field equations, and a non-integrable 𝖣 fails;
    - the a₁, a₂, a₃ delays (formula, an independent fit, and validity at several x) and their identities;
    - Θ maps σ to σ̄;
    - the type-II conditions for sine-Gordon and Tzitzéica;
    - BCZ (5.5) with B = 1 − 𝖣 is solved by A = 2Σ(μ·u)³ (their density (5.2) has ⅓A; the drafted "A = ⅔Σ" was a factor-3 slip, corrected 2026-10-12);
    - the cubic spin-2 obstruction vanishes for r = 2–8;
    - explicit B₊₂ and B₋₂ for a₂, both types, σ symbolic: flux − dB/dt = 0 exactly (sympy; since 2026-10-12);
    - the mirror defect (𝖣ᵀ, −E) has a one-dimensional space of L̃ (since 2026-10-12).
- `defect_lax_rmatrix.py` (6 checks):
  - the Dirac bracket K⁻¹ on the defect values;
  - {L_d ⊗ L_d} = +[r, L_d⊗L_d] for a₁–a₄, both types and couplings, with coefficient exactly 1;
  - the bulk sign c = −1;
  - involution of τ_d, with controls;
  - regularity of L̃ at λ → 0 and λ → ∞.

**By hand (drafting agent).**
- The automorphism/orbit argument behind `prp-ar-only`; the cot identity.
- The spin assignment at the two ends of the λ-plane.
- Dirac bracket = K⁻¹.
- Θ for type II.
- The type-II rescaling prefactors (sine-Gordon, Tzitzéica).
- The Penati–Polvara dictionary row (their Euclidean orientation not checked).
- The time-like canonical transformation, derived in the text and not scripted.

**Not in the repo.** The a₃ spin-±2 run (about 6 min per case, from the phase-0 session); not independently reproduced (Codex deep referee).

**Referee (Codex deep, 2026-10-12).** Report `briefs/referee_defects_38_40_report.md` (items 1–6, 15–17 for ch. 38), scripts in `briefs/referee_defects_38_40_scripts/` (rerun, outputs identical); integrated the same day, see its "Integration (2026-10-12)" section. All items accepted or correct; none rejected. Main changes: the tensor normalization A = 2Σμ³; `prp-spin-two` label (Theorem for (1) and for a₂ in (3), exact; Numerical a₃; sketch (2)), "every spin of the bulk hierarchy", vacuum-plus-analyticity argument for nonzero leading diagonal entries; `prp-defect-lax` mirror; `prp-defect-involution` is now **[Numerical for a₁–a₄]** and explicitly about independent defect values (the step to involution exact); ch. 30 `prp-sklyanin` transport sign fixed; `prp-defect-theta` scope (massive, β ≠ 0, σ ∈ ℂ^×); PP dictionary γ/2 = b_SG²/8π → −β_SG²/8π; corrected CZ 2009b (4.5); no priority claimed for the spin-two cancellation.

**Unchecked, quoted only.**
- Bristow–Bowcock's b_n, c_n, d_n potentials and d₄ delays; Bristow's Lax matrices.
- Robertson's folded defects and fusing rules.
- cz2011's a_n type II.
- Habibullin–Kundu's lattice construction; Avan–Doikou's discrepancy claim; Doikou's gl_N defects.
- Caudrelier's Liouville-integrability proof.
- The Fordy–Gibbons attribution.
- cz2009a's "partner must be a neighbour" step ((2.30)–(2.35)).

**Corrections to sources recorded in the Notes (script-checked).**
1. BCZ 2004 §5 and Penati–Polvara App. B: both spin-±2 charges survive. The obstruction (5.11), and the cubic terms of PP (B.4), vanish identically by the cot identity.
2. BCZ 2004 (5.43) is the transpose of their (3.28).
3. cz2007 (2.21): the self-conjugate delay is −coth, not +coth.
4. cz2009b (4.4): the corrected e^{λ₁} form, with the consequent change to (4.5).
5. cz2009b (3.7): 1/σ multiplies all of f.
6. cz2009b (4.11): the printed z₂ is z₃.
7. cz2009b (4.5) (added 2026-10-12, from the referee): with the corrected (4.4), e^{2i Im λ₁} has its last two factors inverted, so the phase of e^{λ₁} winds by 2π when ρ, ρ̃ > 0, not when both are negative as CZ state.

**Pending edits elsewhere.**
- **Ch. 30 `prp-sklyanin`: DONE 2026-10-12** (re-derived: local {a_x⊗a_x} = +[r, a₁+a₂]δ, four-site lattice transport c = −1). Still open: `boundary_classical.py`'s `tau_bracket` docstring says {T⊗T} = +[r, TT]; harmless (homogeneous), but the docstring could note the sign (file outside the integration list).
- **Appendix D (literature corrections; not on this branch):** when the Part VIII corrections are added, include corrections 1–7 above (7: CZ 2009b (4.5) winding reversal) and state BCZ's tensor as A = 2Σ(μ·u)³ (density ⅓A), not ⅔Σ.
- **Conventions brief:** BCZ's A_abc = 2Σμμμ (density ⅓A); {T⊗T} = −[r, T⊗T] for ∂_xT = −a_xT while {a_x⊗a_x} = +[r, a₁+a₂]δ; PP γ/2 = b_SG²/8π (sinh-Gordon); the defect bracket is for independent defect values.
- **Pending (ch. 38):** retain a reproducible a₃ spin-±2 run (e.g. a separate slow script).
- `textbook_code/README.md` rows: `defect_classical_an.py` (24 checks, ~55 s), `defect_lax_rmatrix.py` (6 checks, <1 s).
- **Conventions brief:**
  - the monodromy T₊L_d⁻¹T₋ and the bracket results;
  - λ → ∞ goes with negative spins; B₊ₛ carries σ⁻¹g_j and B₋ₛ carries σf_j;
  - the a_r delay `eq-dc-ar-delay` with its identities;
  - the type-II sectors and normalizations;
  - corrections 2–6.
- **Optional bib:** Fordy–Gibbons; Aguirre–Araujo–Gomes–Zimerman.

**Obligations for later chapters.**
- **Ch. 39:**
  - use `eq-dc-ar-delay`;
  - the poles of z_a at θ = θ_d − i(π/2 − πa/h);
  - the selection rule;
  - the type-I₊ a_r delays (not yet derived);
  - charged and odd sectors.
- **Ch. 40:** the classical limit of the Borel L-operator should reproduce L_d and its Weyl-type Poisson algebra.
- **Ch. 41:**
  - `eq-dc-ar-delay` for the semiclassical test;
  - the corrected Tzitzéica delay labels.
- **Ch. 42:** the K-matrix dressed by L_d.
- **Ch. 43 (open problems):**
  - e₆–e₈ defects;
  - `cnj-momentum-integrability`;
  - time-like Liouville integrability for a_r;
  - the quantum fate of the second chirality (Penati–Polvara assumed Q₋₂ lost);
  - the equal-time reduction: can the defect values be treated as independent coordinates Poisson-commuting with the bulk (for a_r; Caudrelier's time-like route for sine-Gordon)?
  - the all-rank Drinfeld–Sokolov argument of `prp-spin-two` (2): infinite-line regularization, exceptional defect values;
  - priority of the spin-two correction is not established (bounded literature search).

## Chapter 40 (`sec-defect-qg`)

**By script.**
- `defect_universal_r.py` (16 checks):
  - q-oscillator relations on Weyl and Fock modules; Borel relations and q-Serre.
  - The shift-module classification.
  - Khoroshkin–Tolstoy ℛ on V⊗V against PŘ, with error O(x⁴); controls without the imaginary-root factor or the height-3 roots fail.
  - ℛ on V⊗(Weyl or Fock, both Borels) against `eq-l-operator`.
  - ℛ₂₁⁻¹ on V⊗𝒲⁻ equals c·(ch. 37 matrix) with c → q^{−s/2}.
  - Fock intertwiners are unique; RTT holds for both Borels and both modules.
  - The Hom dimensions of hypothesis (ii) of `thm-t-triangle` on windows, and (since 2026-10-12) exactly without truncation: End(V⊗V⊗𝒲) = scalars for the Weyl modules of 𝔟± (transported loop commutant) and for ℱ⁺ (lowest-block cyclicity), sympy at a rational point.
  - `prp-sg-t-unique` for s = 0, 0.3, 1; (since 2026-10-12) at q = e^{iπ/3} the shift S³ commutes with the 𝒲⁻ generators, so uniqueness needs q not a root of unity.
  - Weston (2.9) intertwines, and its transpose has the 𝔟⁻ invariants.
  - The ψ → −ψ two-label map gives 𝕋_{I₊} with s₊ = −1−s.
  - With q^{−s/2}Φ, unitarity (real s) and crossing hold for s ≠ 0.
  - No extension of 𝔟⁻ by e_j acts on 𝒲.
  - Classical-limit structure: Weyl relations, factorization, RTT.
- `defect_transfer.py` (9 checks; N ≤ 6 for the Weyl TQ, λ = 1.37):
  - Fock trace stable under window growth at |ω| = 0.8.
  - Fock commutators vanish, with controls.
  - Fock TQ with no fitted constant.
  - Weyl factorization, q^{mn} for 𝒲⁻ and q^{−mn} for 𝒲⁺ (since 2026-10-12), and the Weyl window trace does not converge.
  - Weyl commutators need the sector twist ω = q^{−m/2}.
  - Weyl TQ for N = 2, …, 6 (N = 5, 6 added 2026-10-12).
  - Twist limit; T₁ tends to the Bloch-twisted transfer matrix.
  - The charge-conjugation limit, and 𝒟 → q^{ms/2}.

**By hand (drafting agent).**
- The proof of `thm-t-triangle`, and the 𝔟⁻ uniqueness proof.
- q-Serre via Z.
- `eq-dq-massless-defect`, and the chirality of the two defect-potential terms.
- The τ anti-isomorphism argument.
- The Weyl factorization.
- The Poisson brackets {u, w±} and the dictionary to L_d.
- The Ising and Lee–Yang eigenvalues.
- The train argument.

**Exploratory only (scratchpad).** The fit that produced the TQ coefficients; that ±iy_F both work; that the TQ coefficients do not depend on ω.

**Unchecked.**
- The Fock exact sequence in book conventions.
- The Lee–Yang eigenvalue assignment against BHW.
- Hernandez–Jimbo and Frenkel–Hernandez (cited qualitatively).
- Runkel and Petkova–Zuber (quoted).
- KL's same-chirality terms not surviving exactly: inferred from check 13, not derived. (2026-10-12: the exact statement that the generic Weyl module of 𝔟⁻ does not extend to U_q(ŝl₂) is now proved in the text; the all-orders massive charges are still not constructed.)
- The duality q̂ vs q (a structural remark).
- Cyclicity of the averaged Weyl trace (sketch; the text now calls it a formal twisted trace).

**Referee (Codex deep, 2026-10-12).** Report `briefs/referee_defects_38_40_report.md` (items 7–17 for ch. 40), integrated the same day (see its Integration section). All accepted or correct; none rejected. Main changes: global scope "q not a root of unity, q^s fixed"; `thm-t-triangle` (iii) needs a nonzero matrix element, and Hopf closure is called sufficient, not necessary; hypothesis (ii) now also proved exactly on the infinite modules; `prp-sg-t-unique` hypotheses corrected (q⁴ ≠ 1 was insufficient: q = e^{iπ/3} counterexample), level-zero argument written out; `prp-two-borels` (4) unitarity only for real s; `prp-defect-transfer` (1) domain |q| = 1, (2) restricted to 𝒲⁻ (𝒲⁺ has exponent −mn, support ω² = q^m), formal twisted trace, ω = ±q^{−m/2}; `cnj-defect-tq` numerical for N ≤ 6; bootstrap reading qualified (image/kernel identification not established); `prp-defect-twist-limit` scope (λ > 0, real rapidities, real s) and the Gamma ratio tends to e^{−iπ(1+λ)/4}, not 1; Schur needs the full chiral algebras; compact-boson lattice caveat; exact no-extension argument; Ising disorder field sits in the ε-twisted sector (`exr-dq-verlinde` fixed); BLZ trace convergence domain and analytic continuation; Rossi–Weston credited for the delta factor, stripped Baxter matrix, neutral-sector TQ (to N = 10) and commutation.

**Results to note.**
- The universal-R direction is the opposite of the plan's guess: 𝕋_{I₊} ∝ P(π⊗ρ⁺)(ℛ) and 𝕋_{I₋} ∝ P(π⊗ρ⁻)(ℛ₂₁⁻¹), with ℛ ∈ U_q(𝔟⁻)⊗U_q(𝔟⁺).
- Defect operators need inverse L-operators.
- Weyl traces converge for no twist; they factorize by charge sector.

**Corrections to sources recorded in the Notes.**
- Weston's 𝔟⁺ is the book's 𝔟⁻ via his transposition.
- Weston's "regularized trace" for r₁ = r₂ = 0 at |q| = 1: now worded as a clarification (no ordinary convergent twist; formal delta factor), credited to Rossi–Weston 2002 §4.3, (5.5)–(5.10), not as an error.
- BCZ (3.57)–(3.59): 𝕋_{I₊} ≠ (1⊗C)𝕋_{I₋} in the same labels.
- BGKNR's ℛ is the book's ℛ₂₁⁻¹; BLZ III's B₋ is the book's U_q(𝔟⁺) with q → 1/q.

**Pending edits elsewhere.**
- Ch. 36 `sec-defect-circle`: DONE 2026-10-12 (inverse L-operators, convergent twisted trace, Fock = Q-operators, Weyl → `prp-defect-transfer`).
- `textbook_code/README.md` rows for `defect_universal_r.py` (16 checks, ~90 s) and `defect_transfer.py` (9 checks, ~30 s), plus the new `defect_common` helpers (`q_oscillator_module`, `fock_module_b`, `borel_l_operator`, `flip_VW`); `defect_classical_an.py` is now 24 checks, ~110 s; `defect_lax_rmatrix.py` 6 checks, ~2 s.
- Docstring trap in `defect_common.sg_tmatrix`: its Weyl module uses μ = 1, not the module the matrix intertwines.
- Conventions brief §6:
  - q generic (not a root of unity) unless stated; at q = e^{iπ/3} extra intertwiners (S³⊗1)𝕋;
  - Weyl factorization exponents: +mn for 𝒲⁻, −mn for 𝒲⁺ (supports ω² = q^{∓m});
  - the Gamma ratio of `eq-sgd-scalar` tends to e^{−iπ(1+λ)/4} as θ → −∞ (cancels the prefactor);
  - Fock modules are done;
  - the 𝕃± closed forms and the universal-R direction;
  - the universal-R normalization q^{±s/2};
  - the defect operator with inverses, and the sector twist;
  - the BGKNR/BLZ translation rules;
  - the massless defect.
- Optional bib: Verlinde 1988.
- The ch. 37 items now resolved: q^{−s/2}Φ is checked numerically, and Weston's Borel is identified.

**Obligations for later chapters.**
- **Ch. 41:**
  - the a_r Borel modules (CZ 2011) and hypothesis (ii) of `thm-t-triangle` for them;
  - fix the s-dependence of the scalars through the universal-R normalization;
  - keep q-sign statements convention-qualified.
- **Ch. 42:**
  - decide `cnj-twisted-diagonal`, including how the RTT algebra sits in the doubled reflection-equation algebra;
  - 𝒲 as a dynamical boundary module;
  - the third column of `tbl-boundary-defect-algebra`.
- **Ch. 43 (open problems):**
  - a proof of `cnj-defect-tq` for all N;
  - the meaning of ω = q^{−m/2};
  - KL's same-chirality terms beyond first order;
  - q̂ vs q;
  - convergence of Runkel's defect operators;
  - a proof of `cnj-defect-tq` via Weston's exact sequence (Sec. 2.4.1) in book conventions, with the scalar factors, and the image/kernel identification of the bootstrap reading; a detailed comparison with Rossi–Weston's stripped Q-matrix and their q^n = 1 condition in charged sectors;
  - analytic control of the Fock Q-operators at physical q (auxiliary–auxiliary intertwiners, mutual commutation) and a canonical regularization of the bilateral Weyl trace;
  - defects at roots of unity (q rational in λ): reducible Weyl modules and extra intertwiners;
  - the compact-boson realization of the massless defect gluing (T-duality interface vs a line of one theory);
  - whether every transmitting-defect symmetry must be a Hopf subalgebra (only sufficiency is shown).

## Chapter 39 (`sec-defect-solitons`)

The full drafting report is in `briefs/draft_ch39_report.md`: summary, script output, verification record, corrections and requests.

**Key results.**
- **Jump classes.** A jump in class p ∈ Λ_W/Λ_R acts like a jump-free defect at θ_d ± 2πip/h.
- **The general delay** `eq-an-defect-delay` covers both types and all classes. The type-I₊ delays are new, derived by time reversal and confirmed by a fit. As θ_d → −∞, type I₋ becomes the π-rotation line.
- **`prp-stored-charge` and the selection rule.** The outgoing charge moves at most a/gcd(a,h) steps along the π-orbit, and the defect never produces the missing charges.
- **`prp-ds-absorption` [Theorem].**
  - Absorption occurs at θ_d ∓ i(π/2 − π(a+2p)/h).
  - The absorbed or emitted envelope moves with velocity tanh θ_d.
  - Absorption at real rapidity happens iff a + 2p ≡ h/2 mod 2h.
- **Linear transmission.** T_k = z_k(θ−θ_d+iπ/2).
- **`prp-defect-frozen-soliton` [Theorem], in a twisted form.** A type-I₋/I₊ pair at conjugate complex parameters gives ω^b A_cb and ω^b X_cb, with the soliton's energy-momentum.

**Corrections to sources recorded in the Notes.**
- cz2007 (7.11): the denominator, and the implicit m = 1.
- cz2007 §2: sector boundaries are 2π·gcd(a,h)/h apart.
- "Only self-conjugate solitons are absorbed at real rapidity" (cz2007, cz2009a) holds only for root-lattice jumps.
- cz2009a after (4.8): static energies are also real for h even with classes differing by h/2.

**Unchecked.**
- Multi-soliton passage beyond two solitons (argued from the Bäcklund property).
- The charge change through the singular soliton (argued).
- The quantum statements of cz2007 §6 and cz2009a (5.8).

**Pending edits elsewhere.**
- `textbook_code/README.md` row for `defect_solitons_an.py`.
- Ch. 38, end of `sec-defect-backlund-an`: restrict "no absorption at real rapidity except self-conjugate" to class 0, and point to `prp-ds-absorption`.
- Optional pointers: ch. 37 `prp-half-soliton`, ch. 10 `sec-topological-charges`, ch. 16 `eq-an-transmission`, ch. 36 `sec-topological-lines`.
- Conventions brief:
  - jump classes;
  - the I₊ delays and linear factors;
  - T = z continued;
  - the strip convention for complex rapidities;
  - notation: γ_a = πa/h, M̂ = 2hm/β², ν for frequency.

**Obligations for later chapters.**
- **Ch. 41:**
  - the selection rule in the quantum a_r matrices;
  - the Wigner test with |z| and arg z;
  - charged sectors at θ_d ± 2πip/h;
  - absorption poles in charged classes against cz2009a's unstable-state windows.
- **Ch. 43:**
  - is the I₋I₊ pair a soliton S-matrix times the π-line, and is it a type-II defect?
  - lightest-breather factors equal to T^cl for all a_r.

## Chapter 41 (`sec-soliton-transmission`)

The full drafting report is in `briefs/draft_ch41_report.md`. The exploration checkpoint is `briefs/draft_ch41_checkpoint.md`. The first agent was cut off by the quota, and a second agent finished from the checkpoint.

**Re-run on integration.**
- `defect_an_transmission.py`: 21 PASS in 176 s, run in parallel with five other scripts.
- `defect_common.py` gained 251 appended lines and no other change. The other scripts that import it still pass: axioms 20, doubling 17, transfer 9, universal_r 14, sg_transmission 24.

**Convention.** In this chapter ω = 4π/β² − 1, as in Part VI (Corrigan–Zambon's γ), and X = e^{ω(θ−θ_d−iπ/2)}, which is CZ's x̂. The chapter writes ch. 39's root of unity as e^{2πi/h}. See the clash noted under the global pending edits.

**Key results.**
- **`prp-an-transmission` [Established; numerical re-derivation].**
  - The a_r Weyl module `eq-dt-weyl` (offset s = λ_p for class p) has unique intertwiners.
  - Closed forms `eq-dt-tv`, `eq-dt-tvbar`, `eq-dt-tplus` agree with cz2009a (5.8), (5.11) in the same gauge.
  - The triangle relation holds.
  - **Selection rule.** M_V obeys the classical rule exactly. M_V̄ tunnels to the forbidden weights with relative size min(|X|^j, |X|^{j−h+1}).
- **`prp-transmission-crossing` [Numerical].**
  - The matrices are not unitary for r ≥ 2.
  - Crossing holds with C = 1, c_I = 1 + (−q)^r X^h and c_II = 1 − qX^h, independently of s.
  - The bootstrap gives ḡ = ∏g. The scalar factor is given in the convergent form `eq-dt-g`.
  - **The q-sign, in the book's conventions.** At −q crossing fails for h odd with both C and C′, and for h even it holds with C′ = diag((−1)^i).
- **`prp-a22-transmission` [Numerical].** Algebra, unique intertwiners, the triangle relation, and crossing at x_c = −q^{−6} with an antidiagonal C. What remains is listed in the chapter.
- **`prp-a2-unstable` [Numerical; windows Established].**
  - The residue has rank h−1 per sector (`eq-dt-pole`).
  - The window is `eq-dt-window`.
  - In book labels, species 1 is unstable for ω > 1, which includes the classical limit; species 2 for 1/2 < ω ≤ 2 (cz2007's window).
  - **Kernel bootstrap.** It gives the class-1 module with θ_d → θ_d + 2πi/h, the quantum form of `eq-ds-class-shift`.
  - Classical real absorptions become resonances at θ_d − iπ/(hω).
- **`prp-transmission-semiclassical` [Numerical].** Wigner test against `eq-dc-ar-delay`: deviations 0.092/0.025/0.0063 (a₂) at ω = 10/40/160.

**Verification record.**
- **By script:** all items above (21 checks); window sizes are in the report.
- **By hand:**
  - q-Serre;
  - the gauge match with cz2009a;
  - (6.8) from the crossing relations;
  - the universal-R reading of the diagonal entries (the U_q(a_r^(1)) universal R was not evaluated on 𝒲);
  - the class-p poles and windows;
  - the resonance claim;
  - the physical-strip statement.
- **Unchecked:**
  - Robertson's fused a₃ matrix;
  - cz2007's a₂ matrices entry by entry;
  - uniqueness for r ≥ 4, and hypothesis (ii) for r ≥ 3;
  - the scalar factors of species 1 < a < r;
  - the excited-defect scalar;
  - the type-I₊ species-1 closed form and the I₊ class shift;
  - the θ_d = 0 candidate poles;
  - most of a₂^(2): gauge invariants, scalar, dictionary, pole positions, semiclassics;
  - a Klein-type ±q relation for h odd.

**Corrections to sources recorded in the Notes.**
- cz2009a (6.10): the Gamma functions are paired across different k, so the partial products drift like K^{−z}. The equal-k pairing converges.
- cz2009a §6: for a < h/2 the window is 1/γ < a, not "h/2 − 1 or (h−1)/2".
- cz2007 §6: their species 1 is the book's species 2, because of the transposed defect matrix.
- cz2007 (6.14): (−q)^{+1/2}.
- cz2009a §5: the tunnelling in (5.11) is made quantitative.

**Pending edits elsewhere.**
- `textbook_code/README.md` row for `defect_an_transmission.py` and the ch. 41 helpers.
- Ch. 36 `sec-defect-bound-states` and ch. 37 `sec-sg-defect-spectrum`: "the a₂ unstable state exists only in a window of couplings" holds only for species 2; species 1 is unstable for ω > 1.
- Optional: ch. 39, after `prp-defect-frozen-soliton`, a pointer to `sec-defect-open-problems`.
- Conventions brief:
  - the a_r Weyl module, μ and X;
  - orientation −1 = the charge-conjugate defect;
  - C and C′ for the q-sign test;
  - c_I and c_II are independent of s;
  - the quantum class shift;
  - `defect_crossing_residual` uses Cᵀ.

**Obligations for later chapters.**
- **Ch. 43:**
  - `prp-lightest-breather-universal`: bootstrap the lightest breather from these matrices (cz2007 §7 for a₂); the a₂^(2) breather is unchecked.
- **Ch. 43 open problems:**
  - the (h−1)-dimensional image of the residue;
  - candidate stable states at θ_d = 0;
  - the I₊ species-1 closed form and the class shift;
  - the fused I₋I₊ pair versus the soliton S-matrix;
  - the remaining a₂^(2) work;
  - Robertson's fused matrices;
  - Penati–Polvara versus the axiomatic matrices;
  - the q-sign for h odd;
  - hypothesis (ii) for r ≥ 3;
  - the U_q(a_r^(1)) universal R on 𝒲;
  - the a_r doubled-boundary expectation of ch. 42, testable with these matrices.

## Chapter 42 (`sec-defect-boundary`)

The full drafting report is in `briefs/draft_ch42_report.md`, including the proposed row texts for the third column of `tbl-boundary-defect-algebra` and a draft conventions §10. The drafting checkpoint is `briefs/draft_ch42_checkpoint.md`. Re-run on integration: 17 PASS, 30 s.

**Key results.**
- **`prp-reflected-copy`.** The reflected copy is acted on through Δ^op, so the doubled bulk algebra is U⊗U^cop. Through ω_q the reflected copy is the conjugate multiplet: parity maps solitons to antisolitons.
- **`prp-diagonal-coideal`.** Every invertible line has a twisted-diagonal symmetry (1⊗ψ)Δ(U): the identity for the transparent defect, a torus twist for the twist line, the node swap τ for the charge-conjugation line. Numerically the largest left coideal inside the commutant of the doubled K is exactly this twisted diagonal.
- **`cnj-twisted-diagonal` is refuted as stated.** For type I₋ and type I₊ alike, no twisted diagonal is a symmetry; the residuals are 0.93–0.97, stable under window growth.
- **`prp-diagonal-borel`. What holds instead:**
  - Δ(U_q(𝔟^∓)), with 𝒲 the Weyl module as a dynamical boundary module. This is the expected outcome.
  - (1⊗τ)Δ(U_q(𝔟^±)) on a second module structure ρ̃ of 𝒲 (`eq-dbd-second-module`). This was not predicted.
  - The two halves do not combine: any coideal containing both contains e_j⊗1, which fails.
  - As θ_d → +∞ the untwisted diagonal is restored, and as θ_d → −∞ the τ-twisted one.
  - The RTT algebra is the off-diagonal part of the doubled reflection-equation algebra, with 𝕃⁻¹ adjoined.
- **`prp-defect-fusion-k`.**
  - A type-I₋ defect fused onto Dirichlet equals Ghoshal–Zamolodchikov with (η, ϑ) = (η_D, λθ_d), sector by sector and with the scalar factor, so no CDD factor appears. With the UV–IR conjecture this gives λθ_d = (λ+1)θ_d^cl, consistent with `cnj-defect-parameter`.
  - The sinh-Gordon identity (Bajnok–Simon) is re-checked.
  - The classical dressing is K ↦ σ(L_d(1/λ)⁻¹) K_ψ L_d(λ)⁻¹.
- **`prp-fused-bdg`.**
  - The dressing gives conserved generating functions for every a_r boundary in the Durham table, including Neumann, and hence new classical dynamical a_n boundaries. These partly answer Part VII problem 9; involution and quantization are not examined.
  - The sine-Gordon defect potential equals BDG's H_osc at complex (p̂, q̂).
  - The fused system is not the BDG boundary.
- **`sec-tensor-k`.** Appel–Vlaar, with one numerical check. Jump defects fall outside their categories.
- **`sec-noninvertible-boundaries`.** Shimamori–Yamaguchi; their modified boundary crossing becomes the usual Ghoshal–Zamolodchikov form after kink renormalization (`exr-dbd-boundary-crossing`).

**Verification record.**
- **By script:** every [Numerical] item (`defect_doubling.py`, 17 checks).
  - Window sizes: N = 5, 7, 9 for the twisted-diagonal tests and fusion; N = 6, 8 for the diagonal Borel.
  - The scalar identity agrees only to 3e-9, the accuracy of the Richardson-summed Φ at complex rapidity.
- **By hand:**
  - the [Theorem] items;
  - the ω_q⊗id map;
  - ρ̃ and the e_j⊗1 obstruction;
  - the dressing argument;
  - the Shimamori–Yamaguchi renormalization.
- **Unchecked or quoted:**
  - the theorems of Bao–Wang, Wang, Appel–Vlaar and Shimamori–Yamaguchi;
  - the reading ρ̃ ↔ (ψ → −ψ), flagged as an interpretation;
  - boundary crossing-unitarity of K_f;
  - the classical dressing for type I₊;
  - the a_r expectation;
  - two uncited CFT remarks.
- **Scope flag for the referee:** `prp-fused-bdg` (new classical a_n dynamical boundaries) goes slightly beyond "no new research". Check that its label and wording are modest.

**Corrections to sources recorded in the Notes.**
- Bao–Wang Rem. 4.10(2) is a finite-type remark. Kolb 2014 and Appel–Vlaar 2022 assume an indecomposable Cartan matrix, so affine ĝ⊕ĝ is not covered as stated. Kolb–Yakimov's "diagonal type" is unrelated.
- cz2012 (2.18) in book conventions is GZ with (η_D, λθ_d); the identity extends to the scalar factor.
- Shimamori–Yamaguchi's modified crossing is equivalent to conventional crossing after kink renormalization. The Ward–Takahashi identities, not crossing, distinguish the normalizations.

**Pending edits elsewhere.**
- Ch. 40 `tbl-boundary-defect-algebra`:
  - the third column and its status entry (row texts in the report §6);
  - "where it is decided" → "where it is settled".
- Ch. 40 `sec-perturbed-defects`: "twisted by charge conjugation they survive (@prp-diagonal-borel)".
- Pointers to `prp-fused-bdg` from:
  - ch. 35 problem 9;
  - ch. 29 `sec-bsg-dynamical` or ch. 30 `sec-dynamical-boundary`.
- Pointers:
  - ch. 34 `sec-universal-k` → `sec-diagonal-qsp`;
  - optionally, ch. 37 `sec-sg-breather-transmission` → `prp-defect-fusion-k`.
- Repository files:
  - `textbook_code/README.md` row for `defect_doubling.py`;
  - conventions brief §10 (text in the report);
  - optional bib entries: Ito–Terwilliger, Kolb 2020, Lusztig's book.
- The papers 2402.16676 and 2504.08375 have been copied to `~/referee-copy/papers/`.

**Obligations for later chapters.**
- **Ch. 43 (open problems):**
  - the a_r expectation Δ(U_q(𝔟⁻)) + (1⊗π)Δ(U_q(𝔟⁺)), tested with the ch. 41 matrices;
  - the doubled boundaries of type-II defects and of cz2012 (2.19);
  - quantizing the fused a_n dynamical boundaries, and their involution;
  - the doubled boundaries of non-invertible lines;
  - tensor K-matrices with a Borel leg as boundary Q-operators;
  - the affine diagonal-type QSP, not covered by the published universal-K theorems;
  - the algebra generated by the two Borel halves;
  - boundary crossing-unitarity of K_f.

## Chapter 43 (`sec-particle-transmission`)

The full drafting report is in `briefs/draft_ch43_report.md`; the checkpoint is `briefs/draft_ch43_checkpoint.md`. Re-run on integration: 8 PASS, 62 s.

**Notation.** ω = 4π/β² − 1 is always the coupling. Ch. 39's root of unity is written e^{2πi/h}. The ch. 40 twist is written ω_m = q^{−m/2}.

**Key results.**
- **`prp-lightest-breather-universal` [Numerical, 2 ≤ r ≤ 4].**
  - The lowest breather B^(a)_1 of every species crosses a type-I₋ defect as 𝕋^cl_a(θ−θ_d) times the identity on 𝒲, so it is independent of the coupling, the defect label and the offset s.
  - A scratch run covers all five species of a₅.
  - cz2007 (7.5) is the pair (𝕋^cl_2, 𝕋^cl_1) in book labels.
- **`prp-an-particle-transmission` [Numerical, as solutions of the axioms].** For h = 2..8 and both types:
  - unitarity, crossing and Hermitian analyticity;
  - the bootstrap at every fusing;
  - no physical-strip poles for real θ_d;
  - equality with the linearized classical factors.
- **`prp-dpt-vacuum` [Theorem].** At real coupling with σ > 0, ℬ_d ≥ 0, with equality only at φ = ψ = 0. For σ < 0 it is unbounded below.
- **`cnj-transmission-continuation` [Conjecture].** If the real-coupling defect is quantum integrable, the only quantum content of particle transmission is the renormalized θ_d(σ, β). The chapter makes a one-loop prediction for a_r: one common δθ_d for all species.
- **Bajnok–Simon.** Their one-loop factor is exactly the O(b²) term of the exact factor with θ_d = (1 − B/2)θ_d^cl.
- **`prp-dpt-doubled` [Numerical].** Ch. 42's a_r expectation holds for type I₋ of a₂ and a₃: Δ(U_q(𝔟⁻)) and (1⊗π)Δ(U_q(𝔟⁺)), with π(e_j) = e_{j−1}. The opposite rotation fails.
- `tbl-dpt-spectrum` and 15 open problems in `sec-defect-open-problems`. They cover every ch. 43 obligation from ch. 37–42; the mapping is a table in the report, §5.

**Verification record.**
- **By script:** checks (a)–(g).
- **By scratch run:** a₅ every species (about 10 min); (f) on all nodes of a₃.
- **By hand:** the telescoping identities, block forms, `prp-dpt-vacuum`, and the derivation of the duality relation.
- **Unchecked:**
  - the type-I₊ breathers (no I₊ scalar factor yet);
  - class-p breathers;
  - higher breathers for a_r;
  - the module relations of ρ̃ for a_r;
  - a proof for all r.
- **Quoted only:** the defect TBA, DTCSA, Kormos–Runkel–Watts, and Jiang (from the abstract only).
- **Flags for the referee:**
  - the duality remark and θ_d = (1 − B/2)θ_d^cl for a_r are conditional;
  - "for every r" is an expectation;
  - "a soliton CDD factor destroys the coupling independence" is argued, not computed.

**Corrections to sources recorded in the Notes.**
- cz2007 (7.3)–(7.5): their particle labels are swapped relative to the book (transposed defect matrix). Their remark that the factors "depend on β only via η" holds exactly, for every species of a₂–a₄.
- cz2007 §7: η differs from the classical parameter at O(β²) (Bajnok–Simon (23)).
- Bajnok–Simon §4.2: the one-loop factor is the O(b²) term with θ_d = (1 − B/2)θ_d^cl.

**Pending edits elsewhere.**
- Ch. 42 `sec-jump-coideal`, "Beyond sine-Gordon": replace "We have not tested this…" with a pointer to `prp-dpt-doubled`.
- Ch. 37 `sec-sg-breather-transmission`: a pointer to `prp-lightest-breather-universal`.
- Optional:
  - ch. 41, the lightest-breathers paragraph: state the result;
  - ch. 40, `tbl-boundary-defect-algebra`: mention the a_r case.
- Ch. 35 problem 10: point to `sec-defect-open-problems` (part of the global rewrite).
- Repository files:
  - `textbook_code/README.md` row;
  - conventions brief §11 (text in the report).

## Global pending edits

- **Notation clash ω.** Ch. 39 uses ω = e^{2πi/h}; ch. 41 uses ω = 4π/β² − 1, Part VI's coupling. Decide globally, e.g. by renaming ch. 39's root of unity.
- Harmonize theorem titles: Part VII and ch. 40, 42 put a `## Title` inside every `::: {#prp-…}` block; ch. 36, 38, 39 have none and ch. 37 has one.
- Rewrite Part VII open problem 10 (`sec-boundary-open-problems`, ch. 35): integrable interacting defects are purely transmitting.
- Naming: the conventions scripts call the two signs of 𝖣 "types I, II"; the book's names are type I₊ and type I₋ (type II is the Corrigan–Zambon defect with an auxiliary field).
