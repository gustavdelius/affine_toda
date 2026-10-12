# Part VIII (defects): drafting status, verification records, pending edits

Working record for Part VIII, kept in the repository because the session scratchpad is ephemeral. Each drafting agent's verification record is copied here when its chapter is integrated. These lists become the "Already verified" and "To be checked" sections of the referee briefs `referee_defects_*.md`.

## Chapter status

| Ch | Drafted | Script(s) | Checks | Refereed |
|---|---|---|---|---|
| 36 | 2026-10-11 | `defect_axioms.py` | 20 PASS (2 s) | Gemini (light) + Codex deep, integrated 2026-10-12 |
| 37 | 2026-10-11 | `defect_classical_sg.py`, `defect_sg_transmission.py` | 27 + 24 PASS (18 s, 103 s) | Gemini (light) + Codex deep, integrated 2026-10-12 |
| 38 | 2026-10-12 | `defect_classical_an.py`, `defect_lax_rmatrix.py` | 24 + 6 PASS (57 s, <1 s) | no |
| 39 | 2026-10-12 | `defect_solitons_an.py` | 9 PASS (10 s) | no |
| 40 | 2026-10-12 | `defect_universal_r.py`, `defect_transfer.py` | 14 + 9 PASS (46 s, 15 s) | no |
| 41 | | `defect_an_transmission.py` | | |
| 42 | | `defect_doubling.py` | | |
| 43 | | `defect_particle_transmission.py` | | |

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
    - the cubic spin-2 obstruction vanishes for r = 2–8;
    - explicit B₊₂ and B₋₂ for a₂ (dB/dt equals the flux to 1e-15).
- `defect_lax_rmatrix.py` (6 checks):
  - the Dirac bracket K⁻¹ on the defect values;
  - {L_d ⊗ L_d} = +[r, L_d⊗L_d] for a₁–a₃, both types and couplings, with coefficient exactly 1;
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

**Not in the repo.** The a₃ spin-±2 run (about 6 min per case, from the phase-0 session).

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

**Pending edits elsewhere.**
- **Ch. 30 `prp-sklyanin` and its proof:** with ∂_{x₊}T = −a_xT, the transport satisfies {T⊗T} = −[r, T⊗T]; the text implies +. To verify before editing. Involution is unaffected, and `boundary_classical.py`'s `tau_bracket` uses +, which is harmless.
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
  - the quantum fate of the second chirality (Penati–Polvara assumed Q₋₂ lost).

## Chapter 40 (`sec-defect-qg`)

**By script.**
- `defect_universal_r.py` (14 checks):
  - q-oscillator relations on Weyl and Fock modules; Borel relations and q-Serre.
  - The shift-module classification.
  - Khoroshkin–Tolstoy ℛ on V⊗V against PŘ, with error O(x⁴); controls without the imaginary-root factor or the height-3 roots fail.
  - ℛ on V⊗(Weyl or Fock, both Borels) against `eq-l-operator`.
  - ℛ₂₁⁻¹ on V⊗𝒲⁻ equals c·(ch. 37 matrix) with c → q^{−s/2}.
  - Fock intertwiners are unique; RTT holds for both Borels and both modules.
  - The Hom dimensions of hypothesis (ii) of `thm-t-triangle`.
  - `prp-sg-t-unique` for s = 0, 0.3, 1.
  - Weston (2.9) intertwines, and its transpose has the 𝔟⁻ invariants.
  - The ψ → −ψ two-label map gives 𝕋_{I₊} with s₊ = −1−s.
  - With q^{−s/2}Φ, unitarity and crossing hold for s ≠ 0.
  - No extension of 𝔟⁻ by e_j acts on 𝒲.
  - Classical-limit structure: Weyl relations, factorization, RTT.
- `defect_transfer.py` (9 checks; N ≤ 4, λ = 1.37):
  - Fock trace stable under window growth at |ω| = 0.8.
  - Fock commutators vanish, with controls.
  - Fock TQ with no fitted constant.
  - Weyl factorization, and the Weyl window trace does not converge.
  - Weyl commutators need the sector twist ω = q^{−m/2}.
  - Weyl TQ for N = 2, 3, 4.
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
- KL's same-chirality terms not surviving exactly: inferred from check 13, not derived.
- The duality q̂ vs q (a structural remark).
- Cyclicity of the averaged Weyl trace (sketch).

**Results to note.**
- The universal-R direction is the opposite of the plan's guess: 𝕋_{I₊} ∝ P(π⊗ρ⁺)(ℛ) and 𝕋_{I₋} ∝ P(π⊗ρ⁻)(ℛ₂₁⁻¹), with ℛ ∈ U_q(𝔟⁻)⊗U_q(𝔟⁺).
- Defect operators need inverse L-operators.
- Weyl traces converge for no twist; they factorize by charge sector.

**Corrections to sources recorded in the Notes.**
- Weston's 𝔟⁺ is the book's 𝔟⁻ via his transposition.
- Weston's "regularized trace" fails for r₁ = r₂ = 0 at |q| = 1.
- BCZ (3.57)–(3.59): 𝕋_{I₊} ≠ (1⊗C)𝕋_{I₋} in the same labels.
- BGKNR's ℛ is the book's ℛ₂₁⁻¹; BLZ III's B₋ is the book's U_q(𝔟⁺) with q → 1/q.

**Pending edits elsewhere.**
- Ch. 36 `sec-defect-circle`: the product is of 𝕃⁻¹ for t_ω; for Weyl modules the trace is a distribution in the twist, supported at ω² = q^{−m}; "Baxter's Q-operators" applies to Fock modules.
- `textbook_code/README.md` rows for `defect_universal_r.py` (14 checks, ~50 s) and `defect_transfer.py` (9 checks, ~15 s), plus the new `defect_common` helpers (`q_oscillator_module`, `fock_module_b`, `borel_l_operator`, `flip_VW`).
- Docstring trap in `defect_common.sg_tmatrix`: its Weyl module uses μ = 1, not the module the matrix intertwines.
- Conventions brief §6:
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
  - convergence of Runkel's defect operators.

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

## Global pending edits

- Rewrite Part VII open problem 10 (`sec-boundary-open-problems`, ch. 35): integrable interacting defects are purely transmitting.
- Naming: the conventions scripts call the two signs of 𝖣 "types I, II"; the book's names are type I₊ and type I₋ (type II is the Corrigan–Zambon defect with an auxiliary field).
