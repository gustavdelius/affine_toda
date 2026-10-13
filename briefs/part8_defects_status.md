# Part VIII (defects): drafting status, verification records, pending edits

Working record for Part VIII, kept in the repository because the session scratchpad is ephemeral. Each drafting agent's verification record is copied here when its chapter is integrated. These lists become the "Already verified" and "To be checked" sections of the referee briefs `referee_defects_*.md`.

## Chapter status

| Ch | Drafted | Script(s) | Checks | Refereed |
|---|---|---|---|---|
| 36 | 2026-10-11 | `defect_axioms.py` | 20 PASS (2 s) | Gemini (light) + Codex deep, integrated 2026-10-12 |
| 37 | 2026-10-11 | `defect_classical_sg.py`, `defect_sg_transmission.py` | 27 + 24 PASS (18 s, 103 s) | Gemini (light) + Codex deep, integrated 2026-10-12 |
| 38 | 2026-10-12 | `defect_classical_an.py`, `defect_lax_rmatrix.py`, `defect_spin_two_a3.py` | 24 + 6 + 9 PASS (114 s, 2 s, 151 s) | Codex deep, integrated 2026-10-12; upgrades pass (item 7) integrated 2026-10-13 |
| 39 | 2026-10-12 | `defect_solitons_an.py` | 9 PASS (48 s) | Codex deep, integrated 2026-10-13 |
| 40 | 2026-10-12 | `defect_universal_r.py`, `defect_transfer.py` | 16 + 9 PASS (92 s, 29 s) | Codex deep, integrated 2026-10-12 |
| 41 | 2026-10-12 | `defect_an_transmission.py` | 28 PASS (192 s) | Codex deep, integrated 2026-10-13; upgrades pass (items 1–3, 5, 6) integrated 2026-10-13 |
| 42 | 2026-10-12 | `defect_doubling.py` | 20 PASS (23 s) | Codex deep, integrated 2026-10-13; upgrades pass (item 8) integrated 2026-10-13 |
| 43 | 2026-10-12 | `defect_particle_transmission.py` | 11 PASS (about 1 min; optional (b3) $a_5$ about 2.5 min) | Codex (items 1–4) + Claude (items 5–8), integrated 2026-10-13; upgrades pass (items 9–11, Codex) integrated 2026-10-13 |

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
- `textbook_code/README.md` row for `defect_axioms.py` — DONE 2026-10-12 (global integration).
- `part9-foundations/48-scattering.qmd`: pointer from `prp-exchange-twist` and `prp-bloch-sectors` to `prp-twist-defect`, plus one sentence on the common seam — DONE 2026-10-12 (global integration). New paragraph after the proof of `prp-bloch-sectors`.
- `part7-half-line/28-boundary-scattering.qmd` `sec-boundary-lee-yang`: pointer to `sec-defect-lee-yang` (K_b = 𝕋̄_b K_𝟙 𝕋_b) — DONE 2026-10-12 (global integration).
- `part6-soliton-smatrices/27-rsos.qmd`: pointer to the modified crossing in `sec-topological-lines` — DONE 2026-10-12 (global integration). In the Notes, citing cck2024.
- Conventions brief: record the ZF ordering, the 𝕋/𝕋̄ index rule (𝕋̄ parametrized by θ > 0, depending on θ+θ_d), `eq-defect-crossing`, 𝕃 = P𝕋 with R = PŘ, and that the reflected copy scatters with PŠP — DONE earlier (conventions brief §7).

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
- `textbook_code/README.md` rows: `defect_classical_sg.py` (ch. 36–37, 27 checks, ~20 s) and `defect_sg_transmission.py` (ch. 37/40, 24 checks, ~100 s); the new helpers in `defect_common.py` — DONE 2026-10-12 (global integration).
- Pointers (all DONE 2026-10-12):
  - ch. 3 `sec-sg-classical` → `sec-backlund` — DONE 2026-10-12;
  - ch. 10 `sec-multi-solitons` → `prp-bianchi` — DONE 2026-10-12;
  - ch. 29 `sec-bsg-classical` and `cnj-shg-e` → `sec-sg-breather-transmission`, `cnj-defect-parameter` — DONE 2026-10-12 (Dirichlet bullet; new status bullet after `cnj-shg-e`).
- Optional bib: Rogers–Schief, *Bäcklund and Darboux Transformations* (2002); Habibullin, sine-Gordon on the semi-axis. Confirm what `habibullin_kundu2008` contains — not done (optional): no text needs the books; `habibullin_kundu2008` is still cited by title only.

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
- `defect_spin_two_a3.py` (9 checks; new 2026-10-13): generalizes check 5/11's a₂ construction to general r via an
  exact (sympy) Gram–Schmidt root-data basis, matching `defect_common.an_data`/`dmatrix` for r = 1..4; the closed-form
  A_abc = 2Σμ_kμ_kμ_k solves the BCZ (5.5) matching equation exactly for a₂ and a₃; the cubic obstruction (5.11) and
  the on-shell bulk conservation are exact (sympy) for a₃; explicit B₊₂, B₋₂ for a₃, both types, σ symbolic: flux −
  dB/dt = 0 EXACTLY (sympy), with the same velocity structure as a₂; controls (spin ±1 passes, a non-integrable 𝖣
  fails). Runtime about 162 s (2.7 min), dominated by the non-integrable-𝖣 control.
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

**Not in the repo.** (Resolved 2026-10-13: see `defect_spin_two_a3.py` above. The old a₃ spin-±2 run from the phase-0 session, about 6 min per case, numerical only and not reproducible, is superseded by the new script's exact construction.)

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
- **Ch. 30 `prp-sklyanin`: DONE 2026-10-12** (re-derived: local {a_x⊗a_x} = +[r, a₁+a₂]δ, four-site lattice transport c = −1). Still open: `boundary_classical.py`'s `tau_bracket` docstring says {T⊗T} = +[r, TT]; harmless (homogeneous), but the docstring could note the sign (file outside the integration list). — Docstring note DONE 2026-10-12 (global integration; `boundary_classical.py` rerun, 43 PASS).
- **Appendix D (literature corrections; not on this branch):** when the Part VIII corrections are added, include corrections 1–7 above (7: CZ 2009b (4.5) winding reversal) and state BCZ's tensor as A = 2Σ(μ·u)³ (density ⅓A), not ⅔Σ — DONE 2026-10-12 (global integration).
- **Conventions brief:** BCZ's A_abc = 2Σμμμ (density ⅓A); {T⊗T} = −[r, T⊗T] for ∂_xT = −a_xT while {a_x⊗a_x} = +[r, a₁+a₂]δ; PP γ/2 = b_SG²/8π (sinh-Gordon); the defect bracket is for independent defect values — DONE 2026-10-12 (global integration). (§5, and the PP sign corrected in §9.)
- **Pending (ch. 38):** retain a reproducible a₃ spin-±2 run (e.g. a separate slow script) — DONE 2026-10-13: `textbook_code/defect_spin_two_a3.py` (9 checks, about 162 s), generalizing check 5/11's a₂ construction to general r; the result is now EXACT (sympy) for a₃, not just numerical, so `prp-spin-two` item (3) and its label were upgraded to [Theorem] for a₃ too (`part8-defects/38-defect-conditions.qmd`). Verification: ran the new script to completion (all 9 PASS) and reran `defect_classical_an.py` (24 PASS, docstring-only change) in this session (label upgrade not independently refereed). Refereed 2026-10-13 (upgrades report item 7): correct, with an independent Laurent-polynomial Euler test; the `boundary_G` docstring had the Θ sign wrong (the flux is G = (T+Θ)(φ) − (T+Θ)(ψ), dQ_bulk/dt = εG, defect contribution −εF for D_tF = G, charge spin −ε·spin) — docstring fixed, script rerun 9 PASS (151 s).
- `textbook_code/README.md` rows: `defect_classical_an.py` (24 checks, ~55 s), `defect_lax_rmatrix.py` (6 checks, <1 s) — DONE 2026-10-12 (global integration).
- **Conventions brief** (all items DONE 2026-10-12, §5 and §8):
  - the monodromy T₊L_d⁻¹T₋ and the bracket results;
  - λ → ∞ goes with negative spins; B₊ₛ carries σ⁻¹g_j and B₋ₛ carries σf_j;
  - the a_r delay `eq-dc-ar-delay` with its identities;
  - the type-II sectors and normalizations;
  - corrections 2–6 — DONE 2026-10-12 (§5 'Added from chapter 38' and §8 traps).
- **Optional bib:** Fordy–Gibbons; Aguirre–Araujo–Gomes–Zimerman — not done (optional; the text does not cite them).

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
- `textbook_code/README.md` rows for `defect_universal_r.py` (16 checks, ~90 s) and `defect_transfer.py` (9 checks, ~30 s), plus the new `defect_common` helpers (`q_oscillator_module`, `fock_module_b`, `borel_l_operator`, `flip_VW`); `defect_classical_an.py` is now 24 checks, ~110 s; `defect_lax_rmatrix.py` 6 checks, ~2 s — DONE 2026-10-12 (global integration).
- Docstring trap in `defect_common.sg_tmatrix`: its Weyl module uses μ = 1, not the module the matrix intertwines — DONE 2026-10-12: re-derived (residual 1e-15 with μ = i/(q−q⁻¹), 3.7 with μ = 1) and stated in the docstring, the README and the conventions brief §6; behaviour unchanged.
- Conventions brief §6 (all items DONE 2026-10-12):
  - q generic (not a root of unity) unless stated; at q = e^{iπ/3} extra intertwiners (S³⊗1)𝕋;
  - Weyl factorization exponents: +mn for 𝒲⁻, −mn for 𝒲⁺ (supports ω² = q^{∓m});
  - the Gamma ratio of `eq-sgd-scalar` tends to e^{−iπ(1+λ)/4} as θ → −∞ (cancels the prefactor);
  - Fock modules are done;
  - the 𝕃± closed forms and the universal-R direction;
  - the universal-R normalization q^{±s/2};
  - the defect operator with inverses, and the sector twist;
  - the BGKNR/BLZ translation rules;
  - the massless defect — DONE 2026-10-12 (all §6 items above, 'Added from chapter 40').
- Optional bib: Verlinde 1988 — not done (optional).
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
- **`prp-stored-charge` and the selection rule.** For θ > 0 (incidence from the left) the outgoing charge moves at most min(a, h−a)/gcd(a,h) steps along the π-orbit; for a = h/2 one step iff θ > θ_d (I₋). One-soliton transmission never produces the missing charges.
- **`prp-ds-absorption` [Theorem].**
  - Absorption occurs at θ_d ∓ i(π/2 − π(a+2p)/h), real θ_d; meromorphic configurations (zeros of τ_j on real spacetime at complex θ_abs).
  - The absorbed or emitted envelope moves with velocity tanh θ_d, provided a + 2p ≢ 0 mod h; for a + 2p ≡ 0 mod h there is no localized envelope (since 2026-10-13).
  - Absorption at real rapidity happens iff a + 2p ≡ h/2 mod 2h (one class per species).
- **Linear transmission.** T_k = z_k(θ−θ_d+iπ/2); `prp-an-linear-transmission` is a [Theorem] for every r since 2026-10-13 (exact elimination).
- **`prp-defect-frozen-soliton` [Theorem], in a twisted form.** A type-I₋/I₊ pair at conjugate complex parameters gives ω^b A_cb and ω^b X_cb, with the soliton's energy-momentum.

**Corrections to sources recorded in the Notes.**
- cz2007 (7.11): the denominator, and the implicit m = 1.
- cz2007 §2: sector boundaries are 2π·gcd(a,h)/h apart.
- "Only self-conjugate solitons are absorbed at real rapidity" (cz2007, cz2009a) holds only for root-lattice jumps, the setting of their delay — now a *Clarification*, not *Incomplete* (2026-10-13).
- cz2009a after (4.8): removed from the corrections on 2026-10-13. CZ state a sufficient condition (same representation ⇒ real), which is correct; the h/2 classes and θ_d = 0 are now an addition in the Notes paragraph, and the Appendix D entry is deleted.
- cz2007 §2 bullet: the step bound is min(a, h−a)/gcd(a,h), not a/gcd(a,h) (Notes and Appendix D, 2026-10-13).

**Unchecked.**
- Multi-soliton passage beyond two solitons (expected from the Bäcklund property and Bianchi permutability, not written out; the text now says it is unproved beyond two).
- The charge change through the singular soliton (argued).
- The quantum statements of cz2007 §6 and cz2009a (5.8).

**Pending edits elsewhere.**
- `textbook_code/README.md` row for `defect_solitons_an.py` — DONE 2026-10-12 (global integration).
- Ch. 38, end of `sec-defect-backlund-an`: restrict "no absorption at real rapidity except self-conjugate" to class 0, and point to `prp-ds-absorption` — DONE 2026-10-12 (global integration). Also 'neighbouring weight' → 'along its cycle of charges (@prp-stored-charge)'.
- Optional pointers: ch. 37 `prp-half-soliton`, ch. 10 `sec-topological-charges`, ch. 16 `eq-an-transmission`, ch. 36 `sec-topological-lines` — DONE 2026-10-12 (all four: → `prp-defect-frozen-soliton`, `prp-stored-charge`, `eq-an-linear-transmission`, the π-line limits).
- Conventions brief (all items DONE 2026-10-12, §10):
  - jump classes;
  - the I₊ delays and linear factors;
  - T = z continued;
  - the strip convention for complex rapidities;
  - notation: γ_a = πa/h, M̂ = 2hm/β², ν for frequency — DONE 2026-10-12 (conventions brief §10).

**Referee (Codex deep, 2026-10-12; integrated 2026-10-13).** Report `briefs/referee_defects_39_report.md`, scripts and outputs in `briefs/referee_defects_39_scripts/` (rerun, all numbers reproduced). See its "Integration (2026-10-13)" section; checkpoint `briefs/integration_39_checkpoint.md`. Items 1–9 accepted or partially accepted; none rejected.

**Independently checked (referee scripts, rerun; integration scratch with own code).** Exact rational sewing identity for both types at all orders in E (arbitrary rapidity, class phase and v = e^{iπa/h}); weight congruence α_j·Dλ_p ≡ −2p/h (exact, a₂, a₃); direct-exponential sewing with negative control; pole limits (768 configurations); class energies, algebraic mass shell, absorption EP identity (h ≤ 9); real-pole lists h = 4, 6, 8; zero-envelope and τ-zero examples; sign(θ) charge bookkeeping; selection-rule enumeration incl. (8,3) three steps; mixed-pair factorization, particle factor and EP (h ≤ 10), real-σ realizations; Hessian scattering both directions and the CZ (7.11) source system; D π = −Dᵀ; Θ on delays; two-soliton passage (h = 3, 4, classes 0, 1); bootstrap with either shift order. Integration: degree-argument residual at arbitrary complex E, exact sympy elimination of the linear system, step sets at fixed θ_d, θ_d → −∞ class energy, class-1 breather control.

Main changes:
- `prp-an-defect-delay`: meromorphic scope at complex θ; the all-orders argument is now a degree argument (residual × ∏τ_jτ'_j is a polynomial of degree ≤ 2h in E), not "as they must for a Bäcklund pair".
- `prp-ds-absorption`: real θ_d; the case a + 2p ≡ 0 mod h (cos χ = 0, no envelope) excluded from absorption; θ_d = 0 stationary; meromorphic configurations with τ-zeros (example a₂, a = 1: τ₀ = 0 at (−2.31, −4.30), m = 1); one class per species for real absorption; stored-energy reality at θ_d = 0.
- Class energies: "energy real iff 2p ≡ 0" was wrong at θ_d = 0; now energy and momentum both real iff 2p ≡ 0, energy alone also at θ_d = 0.
- `prp-stored-charge`: θ > 0, incidence from the left, phases off the boundaries; a = h/2 case; attainability over the full family, fixed-θ_d restriction (a₄, a = 2, I₊, θ_d = 0: only k ≤ 1); θ < 0 paragraph; "one-soliton transmission does not create the missing charges".
- `prp-defect-frozen-soliton`: hypotheses (1 ≤ b,c < h, real θ_c, finite nonzero effective parameters); proof of item 1 by the limits 1 and ω^b of each factor (covers coincident zeros/poles); two I₋ factors give A for b = h/2.
- `prp-an-linear-transmission`: [Numerical, r ≤ 4] → [Theorem] for every r (exact sympy elimination, both types and sides).
- Rotation limit: only on the class-0 branch; class-p energy diverges.
- Reality section: Θ-symmetry iff σ real (θ_d ∈ ℝ + iπℤ), class-0/h/2 breather, |z| = 1 when cos(2πp/h) = 0, energy-only reality at θ_d = 0.
- Exercises `exr-ds-class`, `exr-ds-breather`; multi-soliton text ("unproved beyond two").
- Notes: cz2007 spacing bullet (min(a, h−a)/gcd), absorption bullet (setting of their delay), cz2009a (4.8) bullet moved to the Notes paragraph as an addition.
- Script `defect_solitons_an.py`: check 3's energy test had used (a₂, a = 1, p = 1), the zero-envelope case, and passed only because cos χ = 6×10⁻¹⁷ in floating point (T = 5×10¹⁷); now excluded and checked separately. New: degree-argument check at arbitrary complex E; θ_d = 0 reality; sign(θ) bookkeeping; fixed-θ_d step sets; a = h/2; cos χ = 0 cases; τ-zero; one class per species; θ_d → −∞ class energy; |z| = 1 at cos = 0; class-1 breather control; θ_d = iπ Θ-symmetry; sympy elimination; asymptotic obstruction. 9 PASS, 48 s.
- Minimal consistency edits: ch. 38 (imaginary θ_d ∉ iπℤ breaks Θ), ch. 43 (soliton-pole energy complex for real θ_d ≠ 0 except real absorptions; imaginary θ_d ∉ iπℤ). Appendix D: cz2007 (2.23) entry relabelled *Clarification*; cz2009a (4.8) entry removed; known-corrections cz2007 (2.22) entry min(a, h−a)/gcd. Conventions brief §10 updated.

**Open problems and pending items from the referee pass (2026-10-13).**
- Not verified: the identification of the charged complex classes with quantum sectors excluded by Θ (classical conjugation only).
- Not verified: global regularity, integration contours and conserved integrated energies through the τ-zeros of complex-rapidity absorption configurations; the book now calls them meromorphic and states no prescription.
- Not verified: whether the I₋I₊ pair has the auxiliary degrees of freedom of a fused type-II defect (already an open question in ch. 43).
- Not verified: priority of the class-shifted and time-reversed observations beyond the cited papers.
- Noticed, not pursued: for θ_d ∈ ℝ + i(π/h)ℤ a coweight translation composed with Θ, (φ, ψ) ↦ (−φ̄, −ψ̄ + 2πλ_q/β), maps the σ-theory to itself (σ̄ω^{±q} = σ). Whether this twisted antilinear map is a usable reality structure (the constant in `eq-ds-class-shift`, the bulk) is open; relevant to stable excited defects at imaginary θ_d (ch. 37, 43).
- Multi-soliton passage for N ≥ 3: write out Bianchi permutability for a_r (open).

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
- **`prp-an-transmission` [Established (1); Theorem (2)–(5)]** (relabelled 2026-10-13; items renumbered, the selection rule is now item 5; item 4 upgraded from [Numerical] the same day, see below).
  - (1) Closed forms `eq-dt-tv`, `eq-dt-tvbar` agree with cz2009a (5.8), (5.11) in the same gauge.
  - (2) Type I₊: `eq-dt-tplus` for V̄; for V the inverse of the reversed-orientation two-entry matrix `eq-dt-reversed` (sector −Q, offset −s, X → −X₊), via the Chevalley involution ϖ (e↔f, k→k⁻¹, (ϖ⊗ϖ)Δ = Δ^op ϖ). New, proved in the integration; since 2026-10-13 the proof also contains the all-rank coefficient check of `eq-dt-reversed` (X/y + cq = c/q, c = μ/y_d). Refereed: upgrades report item 2 (correct).
  - (3) **Uniqueness for every r ≥ 2** among charge-conserving maps, θ off the (sector- and s-independent) determinant zeros; for q not a root of unity (ω irrational) every intertwiner conserves charge. Proof: Ω_j = Δ(f_j)Δ(f_{j+1}) − q⁻¹Δ(f_{j+1})Δ(f_j) has no pure-𝒲 term; its image (V) or kernel (V̄) is a coordinate line/hyperplane, so End is diagonal, then scalar. New proof; refereed in the upgrades report, item 1 (correct). Hypothesis since 2026-10-13: the matrices defined and invertible on their nonempty sectors (the I₊ V matrix, an inverse, has det q^{−r}/(1 − qX₊^h): poles, not zeros). At q⁵ = 1 (ω = 3/5) the shift by 5α₁ gives non-charge-conserving extra intertwiners. Note ω = 1.37 is rational too (q²⁰⁰ = 1).
  - (4) **The triangle relation for every r ≥ 2**, both types, pairs VV, V̄V̄, VV̄, V̄V (2026-10-13, checkpoint `briefs/triangle_all_r_checkpoint.md`). Proof: ℛ ∈ U_q(𝔟⁻)⊗̂U_q(𝔟⁺) = q^t(1 + (q−q⁻¹)Σ f_j⊗e_j + …) (type I₊) and ℛ₂₁⁻¹ = (1 − (q−q⁻¹)Σ e_j⊗f_j + …)q^{−t} (type I₋), evaluated on V⊗𝒲, are formal power series in y_d/y resp. y/y_d with finite coefficients, charge-conserving intertwiners over ℂ((x)); the uniqueness proof of (3) is field-theoretic, so P·(evaluated ℛ) = φ·M with φ = 1 + O(x²) (orders 0 and 1 agree exactly, all four cases); on V⊗V it is ψ·Ř (generic uniqueness of Ř, rank argument); the YBE on V(y₁)⊗V(y₂)⊗𝒲 is a series in u, uv, v and gives the triangle relation times ψφφ, invertible in ℂ((u))((v)); rational ⇒ identical; one q throughout. Hypothesis (ii) is not used; it remains proved only for a₂ (exact rank certificate).
  - (5) **Selection rule.** M_V obeys the classical rule exactly. M_V̄ has entries to every weight; the ones forbidden for that incoming soliton have relative size min(|X|^j, |X|^{j−h+1}).
- **`prp-transmission-crossing` [Numerical].**
  - The matrices are not unitary for r ≥ 2.
  - Crossing holds with C = 1, c_I = 1 + (−q)^r X^h and c_II = 1 − qX^h, independently of s.
  - The bootstrap gives ḡ = ∏g. The scalar factor is given in the convergent form `eq-dt-g`.
  - **The q-sign, in the book's conventions.** At −q crossing fails for h odd with both C and C′, and for h even it holds with C′ = diag((−1)^i).
- **`prp-a22-transmission` [Numerical].** Algebra, intertwiners unique up to a scalar on the windows tested, the triangle relation, and crossing at x_c = −q^{−6} with an antidiagonal C. What remains is listed in the chapter. cz2011 (3.12): checked that, grouped at equal k, it solves (3.8), (3.10) with z = i(θ−ν)/ξ, provided (3.2) is read as x = e^{2πθ/ξ} (printed q^{2πθ/ξ}, a misprint; 2026-10-13); its outside Gamma poles are at ν − i(ξ/2 + π/3), ν − i(ξ/2 + 2π/3), so the printed positions (ξ) and windows are wrong; from (3.12) the decaying-energy windows are π/3 < ξ < 4π/3 and −π/3 < ξ < 2π/3 (physical status unchecked).
- **`prp-a2-unstable` [Established (1); Numerical (2), (3)]** (relabelled 2026-10-13).
  - The residue has rank h−1 per sector (`eq-dt-pole`), for generic ω.
  - The window `eq-dt-window` is open at both ends, max(0, a−h/2) < 1/ω < a, for the pole in −π/2 < Im θ_a < 0; at ϑ_a = π/2 the pole is real, Im E = 0 (zero width; cz2009a (6.15) includes it).
  - In book labels, species 1 is unstable for ω > 1 (β² < 2π), with the classical absorption as the limit ω → ∞; species 2 for 1/2 < ω < 2 (cz2007's window, strict).
  - **Kernel of the residue.** The one-dimensional kernel carries the class-1 matrix (offset s+λ₁, X^h → q²X^h), i.e. θ_d → θ_d + 2πi/h + 2πik/(hω) for any integer k; every k has the same classical limit. The excited defect is the (h−1)-dimensional image (ch. 36), not the kernel.
  - **[Conjecture]** class-p scalar g(θ − 2πip/h) (k = 0); the charged-class windows and the resonances at θ_d − iπ/(hω) from real absorptions (a₅ species 1 class 1) rest on it; k = 1 gives the opposite sign of Im E.
- **`prp-transmission-semiclassical` [Numerical].** Wigner test against `eq-dc-ar-delay`: deviations 0.092/0.025/0.0063 (a₂) at ω = 10/40/160.

**Verification record.**
- **By script:** all items above (21 checks); window sizes are in the report.
- **By hand:**
  - q-Serre;
  - the gauge match with cz2009a;
  - (6.8) from the crossing relations;
  - the universal-R reading of the diagonal entries; since 2026-10-13 the evaluated universal R is φ·M with φ = 1 + O(x²) (proof of item 4), φ not computed; the first-order terms are also checked by script;
  - the class-p poles and windows;
  - the resonance claim;
  - the physical-strip statement.
- **Unchecked:**
  - Robertson's fused a₃ matrix;
  - cz2007's a₂ matrices entry by entry;
  - hypothesis (ii) for r ≥ 3 (no longer needed: uniqueness and the triangle relation are proved for every r since 2026-10-13);
  - the scalar factors of species 1 < a < r;
  - the excited-defect scalar (bootstrap on the image of the residue) and the class-p scalar (conjecture);
  - the I₊ scalar factor and class shift (the I₊ V matrix is now known as an inverse, item 2);
  - the θ_d = 0 candidate poles (note: det M_V̄ = (q − X^h)^r vanishes to order r at the species-r pole X^h = q; its residue rank is unchecked);
  - most of a₂^(2): gauge invariants, the physical scalar (only the functional equations and poles of cz2011 (3.12) are checked), dictionary, semiclassics;
  - a Klein-type ±q relation for h odd.

**Corrections to sources recorded in the Notes.**
- cz2009a (6.10): the Gamma functions are paired across different k; normalized at z = 0 the partial products are the equal-k ones times Γ(A_K−z)/Γ(A_K), a drift like (hγK)^{−z}. The equal-k pairing converges.
- cz2009a §6: for a < h/2 the window is 1/γ < a, not "h/2 − 1 or (h−1)/2"; and (6.15) must be strict at ϑ_a = π/2 (2026-10-13).
- cz2011 (3.2): the spectral variable is x_i = e^{2πθ_i/ξ}, not the printed q^{2πθ_i/ξ}; only then τ²x = e^{−2πiz} and (3.12) solves (3.8), (3.10) (PDF checked; misprint; upgrades integration 2026-10-13; Appendix D).
- cz2011 (3.12)–(3.13): pole positions with ξ instead of ξ/2, hence wrong windows (in the corrected (3.2) convention); as printed the products converge only after normalization at z = 0, then drift like K^{−2z} (2026-10-13; Appendix D).
- cz2007 §6: their species 1 is the book's species 2, because of the transposed defect matrix.
- cz2007 (6.14): (−q)^{+1/2}.
- cz2009a §5: the tunnelling in (5.11) is made quantitative.

**Pending edits elsewhere.**
- `textbook_code/README.md` row for `defect_an_transmission.py` and the ch. 41 helpers — DONE 2026-10-12 (global integration).
- Ch. 36 `sec-defect-bound-states` and ch. 37 `sec-sg-defect-spectrum`: "the a₂ unstable state exists only in a window of couplings" holds only for species 2; species 1 is unstable for ω > 1 — DONE 2026-10-12 (global integration).
- Optional: ch. 39, after `prp-defect-frozen-soliton`, a pointer to `sec-defect-open-problems` — DONE 2026-10-12 (global integration).
- Conventions brief (all items DONE 2026-10-12, §11):
  - the a_r Weyl module, μ and X;
  - orientation −1 = the charge-conjugate defect;
  - C and C′ for the q-sign test;
  - c_I and c_II are independent of s;
  - the quantum class shift;
  - `defect_crossing_residual` uses Cᵀ — DONE 2026-10-12 (conventions brief §11).

**Obligations for later chapters.**
- **Ch. 43:**
  - `prp-lightest-breather-universal`: bootstrap the lightest breather from these matrices (cz2007 §7 for a₂); the a₂^(2) breather is unchecked.
- **Ch. 43 open problems:**
  - the excited-defect bootstrap on the (h−1)-dimensional image of the residue, and its scalar factor (open problem 10, updated 2026-10-13);
  - the class-p scalar factor g(θ − 2πip/h) and the resonance reading of the real absorptions (conjectural since 2026-10-13; open problem 10);
  - the a₂^(2) physical scalar factor and windows, given the cz2011 (3.12) pole slip (open problem 5, updated 2026-10-13);
  - candidate stable states at θ_d = 0;
  - the I₊ scalar factor and class shift;
  - the fused I₋I₊ pair versus the soliton S-matrix;
  - the remaining a₂^(2) work;
  - Robertson's fused matrices;
  - Penati–Polvara versus the axiomatic matrices;
  - the q-sign for h odd;
  - hypothesis (ii) for r ≥ 3 (not needed for the triangle relation since 2026-10-13);
  - the scalar series φ of the universal R evaluated on V⊗𝒲, and its relation to g (open problem 7, updated 2026-10-13);
  - the a_r doubled-boundary expectation of ch. 42, testable with these matrices.

**Referee (Codex deep, 2026-10-12/13; integrated 2026-10-13).** Report `briefs/referee_defects_41_report.md`, scripts and outputs in `briefs/referee_defects_41_scripts/` (rerun, outputs reproduced byte for byte), integration scripts `integration_*.py` there; checkpoint `briefs/integration_41_checkpoint.md`; record in the report's "Integration (2026-10-13)" section. Items 1–8 accepted or accepted with modification; none rejected. Main changes:
- `prp-an-transmission`: relabelled per item; uniqueness now a theorem for every r (Ω argument, charge conservation, roots of unity); the I₊ V matrix via the Chevalley involution; new `eq-dt-reversed` with the normalization μ → μ/q (cz2007 x̂ = qX otherwise); hypothesis (ii) for a₂ by an exact certificate.
- Windows strict, nearest strip; endpoint zero width; species 2 of a₂ 1/2 < ω < 2; classical limit as a limit (also ch. 36 l. 215, ch. 37 `sec-sg-defect-spectrum`); figure F40 open circles.
- Kernel vs image: the kernel restriction is a matrix equivalence, not the excited defect; the shift is fixed only mod 2πi/(hω); class-p scalar, charged-class windows and resonances are conjectures (ch. 43 table, reading paragraph, open problem 10).
- cz2009a (6.10): exact finite identity with normalization at z = 0 (text, exr-dt-gamma, Appendix D).
- cz2011 (3.12)–(3.13): pole positions ξ → ξ/2, windows doubled; flagged in ch. 41, Appendix D (new entry), ch. 43 table and open problem 5.
- Penati–Polvara: spin-3 *current* (spin-±2 charges; Q₊₂ tested), r ≤ 5, renormalizations of the potential, derivative coupling and improvement terms (ch. 41 caveat, ch. 43 condition 4 and open problem 1; ch. 37 l. 240 already correct: spin-3 charge of sinh-Gordon = spin-4 current).
- Lagrangian caveat rewritten; Notes; exercises serre (cubic for r = 1), intertwiner, crossing (φφ̄ constant), gamma, window ((5,5) representative p = −1), semiclassical (regular soliton, threshold formula).
- Script: 21 → 27 checks, about 200 s.

**Status upgrades without an independent second referee pass (2026-10-13).** `prp-an-transmission` (2) existence of the I₊ V intertwiner for every r, and (3) uniqueness for every r, both [Theorem] with proofs written in the integration; (5) selection rule labelled [Theorem] (read off the closed forms). Hypothesis (ii) for a₂ was certified independently by the referee as well. Item (4), the triangle relation for every r, upgraded [Numerical] → [Theorem] the same day by the universal-R argument (proof written by a Claude agent, no referee pass; supporting checks: a₄ on complete sectors and the first-order terms in `defect_an_transmission.py`, a₂–a₅ in a scratch check).

**Referee pass of these upgrades (Codex, 2026-10-13; integrated the same day).** Report `briefs/referee_defects_upgrades_report.md`, items 1, 2, 3, 5, 6: every upgrade correct, no label downgraded. Integrated: item 1 (uniqueness hypothesis 'defined and invertible'), item 2 (all-rank coefficient check of `eq-dt-reversed` in the proof), item 3 (type-I₊ selection rule: same relative bound with X₊, equal but not unit moduli at θ = θ_d), item 5 (cz2011 (3.2) misprint; 'so the slip is in the pole positions' qualified; normalization at z = 0), item 6 (Cartan lift with a fixed log q; W need not lie in Frenkel–Hernandez's category O; no pinpoint citations, the sources are not available locally). Script: check (j) gains the control with the printed (3.2). Record: the report's Integration section and `briefs/integration_upgrades_checkpoint.md`.

**Pending edits elsewhere (2026-10-13).**
- `appendices/C-code.qmd` row for `defect_an_transmission.py`: 27 checks, about 200 s (was 21, 100–180 s) — DONE 2026-10-13 (coordinator, after the integration).
- Ch. 42: nothing found that repeats a corrected claim.

## Chapter 42 (`sec-defect-boundary`)

The full drafting report is in `briefs/draft_ch42_report.md`, including the proposed row texts for the third column of `tbl-boundary-defect-algebra` and a draft conventions §10. The drafting checkpoint is `briefs/draft_ch42_checkpoint.md`. Re-run on integration: 17 PASS, 30 s.

**Key results.**
- **`prp-reflected-copy`.** The reflected copy is acted on through Δ^op, so the doubled bulk algebra is U⊗U^cop. Through ω_q the reflected copy is the conjugate multiplet: parity maps solitons to antisolitons.
- **`prp-diagonal-coideal`.** Every invertible line has a twisted-diagonal symmetry (1⊗ψ)Δ(U): the identity for the transparent defect, a torus twist for the twist line, the node swap τ for the charge-conjugation line. Numerically the largest left coideal inside the commutant of the doubled K is exactly this twisted diagonal.
- **`cnj-twisted-diagonal` is refuted as stated.** For type I₋ and type I₊ alike, no twisted diagonal is a symmetry; the residuals are 0.93–0.97, stable under window growth. *Since the referee integration (2026-10-13): refuted for every Hopf twist by an exact y-power argument, for arbitrary θ-independent operators on 𝒲; Hopf automorphisms = torus ⋊ diagram automorphisms by Kolb 2014 Thm. 3.2 (Twietmeyer).*
- **`prp-diagonal-borel`. What holds instead:**
  - Δ(U_q(𝔟^∓)), with 𝒲 the Weyl module as a dynamical boundary module. This is the expected outcome.
  - (1⊗τ)Δ(U_q(𝔟^±)) on a second module structure ρ̃ of 𝒲 (`eq-dbd-second-module`). This was not predicted.
  - The two halves do not combine: any coideal containing both contains e_j⊗1, which fails.
  - ~~As θ_d → +∞ the untwisted diagonal is restored, and as θ_d → −∞ the τ-twisted one.~~ Wrong (referee, verified 2026-10-13): these are matrix limits only. Type I₋: diagonal (twist line t_n = q^(−2n−s) on each charge sector) as θ_d → +∞, conjugation with a charge shift as θ_d → −∞; type I₊ reversed. On all of 𝒲 no twisted diagonal acts in either limit (χ(Δk₁) and χ(Δk₁⁻¹) forced to be the same operator).
  - ~~The RTT algebra is the off-diagonal part of the doubled reflection-equation algebra, with 𝕃⁻¹ adjoined.~~ Since 2026-10-13: with vanishing diagonal blocks and boundary unitarity (a quotient) the doubled RE presentation is equivalent to RTT with 𝕃⁻¹ adjoined; no embedding; Δ(U_q(𝔟)) is a left coideal, not a Hopf subalgebra under Δ_A.
- **`prp-defect-fusion-k`.**
  - A type-I₋ defect fused onto Dirichlet equals Ghoshal–Zamolodchikov with (η, ϑ) = (η_D, λθ_d), sector by sector and with the scalar factor, so no CDD factor appears. With the UV–IR conjecture this gives λθ_d = (λ+1)θ_d^cl, consistent with `cnj-defect-parameter`. *Since 2026-10-13 (4) is conditional on the quantum UV–IR matching at finite coupling; the classical limit fixes only θ_d/θ_d^cl → 1.*
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
- Ch. 40 `tbl-boundary-defect-algebra` (DONE 2026-10-12):
  - the third column and its status entry (row texts in the report §6);
  - "where it is decided" → "where it is settled" — DONE 2026-10-12, with the row texts of the report §6 merged into the cells as changed by the ch. 38–40 integration (status cell: '[Numerical] (sine-Gordon; a₂, a₃)').
- Ch. 40 `sec-perturbed-defects`: "twisted by charge conjugation they survive (@prp-diagonal-borel)" — DONE 2026-10-12 (global integration).
- Pointers to `prp-fused-bdg` from (DONE 2026-10-12):
  - ch. 35 problem 9 — DONE 2026-10-12;
  - ch. 29 `sec-bsg-dynamical` or ch. 30 `sec-dynamical-boundary` — DONE 2026-10-12 (ch. 29 `sec-bsg-dynamical`).
- Pointers (DONE 2026-10-12):
  - ch. 34 `sec-universal-k` → `sec-diagonal-qsp` — DONE 2026-10-12;
  - optionally, ch. 37 `sec-sg-breather-transmission` → `prp-defect-fusion-k` — DONE 2026-10-12.
- Repository files (DONE 2026-10-12, optional bib not done):
  - `textbook_code/README.md` row for `defect_doubling.py` — DONE 2026-10-12;
  - conventions brief §10 (text in the report) — DONE 2026-10-12 (now §12, in chapter order);
  - optional bib entries: Ito–Terwilliger, Kolb 2020, Lusztig's book — not done (optional; cited only through other papers).
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

**Integration of the Codex deep referee report (2026-10-13).** Report `briefs/referee_defects_42_report.md` (Integration section at its end); checkpoint `briefs/integration_42_checkpoint.md`; own scripts `briefs/referee_defects_42_scripts/integration_*.py`. All eight items verified (some modified); none rejected.
- Labels changed: `cnj-twisted-diagonal` → [Conjecture; refuted for Hopf twists, exact argument below]; `prp-diagonal-borel` → [Theorem (1), (3)–(5), and (2) for type I₋; Numerical (2) for type I₊] (upgrades of (3), (5) and (2)-for-I₋, by hand + exact sympy; refereed 2026-10-13 in the upgrades report, item 8: correct, except that (5) and the proof wrote 𝕃 = P𝕋 = 𝒟 + y𝒪, which holds only for the matrix part — fixed, see below); `prp-defect-fusion-k` (4) → "by hand, given the quantum UV–IR matching"; `tbl-dbd-symmetries` status cells; ch. 40 table status cell.
- Book script 20 PASS (was 17): exact y-power obstruction, ρ̃ exact for symbolic n, (4) for both types, limits for both types, CZ (2.17)/(2.19).
- Upgrades referee pass (2026-10-13, item 8): 𝕋 = Φ𝕋̂ with Φ the scalar factor; the proof works with the matrix part 𝕃̂ = P𝕋̂ = 𝒟 + y𝒪 (Φ cancels from every intertwining equation); (5) states the limits for 𝕃̂ (the full 𝕃 has the same θ_d → +∞ limit since Φ → 1; 𝕃/(y/y_d) → 0 since |Φ| = (1+Y²)^{−1/2}; 𝕃/(ΦY) has the conjugation limit, |ΦY| → 1); type-I₊ twist matrices diag(q^n, q^{−n−s}); caption 'becomes proportional to'; script and figure-script docstrings. No label change. Script 20 PASS (23 s).
- New correction to the literature: CZ 2012 (2.19), two misprints (h and the sign of s₋), chapter Notes and Appendix D.

## Chapter 43 (`sec-particle-transmission`)

The full drafting report is in `briefs/draft_ch43_report.md`; the checkpoint is `briefs/draft_ch43_checkpoint.md`. Re-run on integration: 8 PASS, 62 s. After the referee integration of 2026-10-13: 11 PASS, 35–50 s (see the integration subsection below).

**Notation.** ω = 4π/β² − 1 is always the coupling. Ch. 39's root of unity is written e^{2πi/h}. The ch. 40 twist is written ω_m = q^{−m/2}.

**Key results.**
- **`prp-lightest-breather-universal` [Numerical, 2 ≤ r ≤ 4].**
  - The lowest breather B^(a)_1 of every species crosses a type-I₋ defect as 𝕋^cl_a(θ−θ_d) times the identity on 𝒲, so it is independent of the coupling, the defect label and the offset s.
  - All five species of a₅: optional check (b3) of the script (`DPT_A5=1`, about 2.5 min, 4e-9); before 2026-10-13 an unrecorded scratch run.
  - The lowest breather is physical only for hω > 2 (Gandenberger 1998 §3.3.1); h = 3, ω = 0.6 in (b2) is an algebraic continuation. Only B₁ of sine-Gordon is coupling-independent; B_n, n ≥ 2, are not.
  - cz2007 (7.5) is the pair (𝕋^cl_2, 𝕋^cl_1) in book labels.
- **`prp-an-particle-transmission` [Theorem since 2026-10-13, was Numerical for h = 2..8; as solutions of the axioms].** For every h ≥ 2 and both types (full hand proof in the proof block; exact sympy in check (c)):
  - unitarity, crossing and Hermitian analyticity;
  - the bootstrap at every fusing;
  - no physical-strip poles for real θ_d;
  - equality with the linearized classical factors.
- **`prp-dpt-vacuum` [Theorem].** At real coupling with σ > 0, ℬ_d ≥ 0, with equality only at φ = ψ = 0. For σ < 0 it is unbounded below.
- **`cnj-transmission-continuation` [Conjecture].** If the real-coupling defect is quantum integrable, the only quantum content of particle transmission is the renormalized θ_d(σ, β). The chapter makes a one-loop prediction for a_r: one common δθ_d for all species.
- **Bajnok–Simon.** Their one-loop factor is exactly the O(b²) term of the exact factor with θ_d = (1 − B/2)θ_d^cl.
- **`prp-dpt-doubled` [Theorem since 2026-10-13 (upgrades integration), type I₋, every r ≥ 2, generic q (q³, q⁴ ≠ 1 suffice), finite θ_d, any s; was Numerical for a₂, a₃].** Ch. 42's a_r expectation holds for type I₋ (originally checked for a₂ and a₃): Δ(U_q(𝔟⁻)) and (1⊗π)Δ(U_q(𝔟⁺)), with π(e_j) = e_{j−1}. The opposite rotation fails.
  - Since 2026-10-13: the second module ρ̃ is explicit (`eq-dpt-second-module`: ρ̃(e_j)|N⟩ = e^{ω(θ_d+iπ/2)}/(q−q⁻¹) q^{−(N+s)_j}|N+α_j⟩, ρ̃(k_j) = plain shift by α_j−α_{j−1}); derived by hand for V (free constant c_j q^{(N+s)_j} killed by the Cartan relation, q⁴ ≠ 1), V̄ numerical (exact labels, h = 3, 4, 5) until the upgrades integration, which added the all-rank proof by the sparse inverse M_V̄⁻¹ = H/(1 − X^h/q), H|i,N⟩ = q^{−n_i}|i,N⟩ − (X/q)|i+1, N−α_i⟩ (from the upgrades report, item 10; check (h) now also checks H M_V̄ = 1 − X^h/q); h = 2 is `eq-dbd-second-module` in the basis τ^{−n}q^{n(n−1)/2+ns/2}|n⟩, τ = −ie^{iπω/2}, with M_V times q^{s/2}. The failing twists fail by an exact y-power argument (e-type ψ(j) ≠ j−1 at y², f-type ψ(j) ≠ j at y⁻¹). The two Borel halves do not form one coideal (q³ ≠ 1, as `prp-diagonal-borel` (4)).
- **Soliton CDD condition** (`eq-dpt-soliton-cdd`, 2026-10-13): ∏_{k<h} f(D − iπ + 2πik/h) = 1. The simpler f(D)f(D+iπ) = 1 is equivalent only for h = 2; for h > 2 it is not necessary, it is sufficient for even h (pair k with k + h/2) and not sufficient in general for odd h (corrected 2026-10-13, upgrades report item 11; do not write 'neither necessary nor sufficient'). The h = 3 case is in cz2007 right after (6.9), cited; the novelty claim is the all-rank form and derivation; T^cl_1 satisfies it for every h; its breather multiplier depends on ω (computed).
- `tbl-dpt-spectrum` and 15 open problems in `sec-defect-open-problems`. Residue (on each fixed-charge sector; infinite rank on the whole direct sum) at the species-1 pole of 𝕋_V: rank h−1; at the species-r pole of 𝕋_V̄ (X^h = q): rank 1; other species not computed. Lee–Yang bound-state poles are in 𝕋̄_b. For θ_d = iv, Θ is broken only for v ∉ πℤ; at v = π (in the window for a < h/2) σ = −1, unbounded below. They cover every ch. 43 obligation from ch. 37–42; the mapping is a table in the report, §5.

**Verification record.**
- **By script:** checks (a)–(g).
- **By scratch run:** (f) on all nodes of a₃. (The a₅ run is now the optional check (b3).)
- **By hand:** the telescoping identities, block forms, `prp-dpt-vacuum`, and the derivation of the duality relation.
- **Unchecked:**
  - the type-I₊ breathers (no I₊ scalar factor yet);
  - class-p breathers;
  - higher breathers for a_r;
  - a proof for all r.
- **Quoted only:** the defect TBA, DTCSA, Kormos–Runkel–Watts. Jiang 2017 was read in full at the referee integration (abstract, intro, (5.5)–(5.6), §7); the "(from the abstract)" qualifier was removed and the bullet rewritten from the paper.
- **Flags for the referee:**
  - the duality remark and θ_d = (1 − B/2)θ_d^cl for a_r are conditional;
  - "for every r" is an expectation;
  - "a soliton CDD factor destroys the coupling independence" is argued, not computed. (Computed since 2026-10-13, with the corrected condition `eq-dpt-soliton-cdd`.)

**Corrections to sources recorded in the Notes.**
- cz2007 (7.3)–(7.5): their particle labels are swapped relative to the book (transposed defect matrix). Their remark that the factors "depend on β only via η" holds exactly, for every species of a₂–a₄.
- cz2007 §7: η differs from the classical parameter at O(β²) (Bajnok–Simon (23)).
- Bajnok–Simon §4.2: the one-loop factor is the O(b²) term with θ_d = (1 − B/2)θ_d^cl.

**Pending edits elsewhere.**
- Ch. 42 `sec-jump-coideal`, "Beyond sine-Gordon": replace "We have not tested this…" with a pointer to `prp-dpt-doubled` — DONE 2026-10-12 (global integration).
- Ch. 37 `sec-sg-breather-transmission`: a pointer to `prp-lightest-breather-universal` — DONE 2026-10-12 (global integration).
- Optional (both DONE 2026-10-12):
  - ch. 41, the lightest-breathers paragraph: state the result — DONE 2026-10-12;
  - ch. 40, `tbl-boundary-defect-algebra`: mention the a_r case — DONE 2026-10-12.
- Ch. 35 problem 10: point to `sec-defect-open-problems` (part of the global rewrite) — DONE 2026-10-12 (global integration).
- Repository files (DONE 2026-10-12):
  - `textbook_code/README.md` row — DONE 2026-10-12;
  - conventions brief §11 (text in the report) — DONE 2026-10-12 (now §13, in chapter order).

**Integration of the referee report (2026-10-13).** Report `briefs/referee_defects_43_report.md` (Codex items 1–4, Claude items 5–8; Integration section at its end); checkpoint `briefs/integration_43_checkpoint.md`; own scripts `briefs/referee_defects_43_scripts/integration_*.py`. All eight items verified (items 6 and 8 with small modifications); none rejected.
- Labels changed: `prp-an-particle-transmission` [Numerical, h = 2..8] → **[Theorem]** for every h ≥ 2 (full hand proof added to the proof block; exact sympy in check (c)); `tbl-dpt-spectrum` row "a_r particles, real θ_d" status → "[Theorem] for the factors". **No independent referee pass of the upgraded proof** at the time; refereed 2026-10-13 in the upgrades report, item 9 (correct; the type-I₊ case is now written out in the same variable x). All other labels unchanged at that integration (`prp-dpt-doubled` stayed [Numerical]; upgraded to [Theorem] in the upgrades integration, see above).
- New ids: `eq-dpt-second-module` (proof of `prp-dpt-doubled`), `eq-dpt-soliton-cdd` (CDD paragraph).
- Book script 11 PASS (was 8): exact identities in (c), explicit σ < 0 direction in (g), (h) second module, (i) soliton CDD condition, (j) residue ranks, Lee–Yang direction, v = π; optional (b3) a₅.
- Minimal consistency edits outside ch. 43: ch. 37 `sec-sg-defect-spectrum` ("Θ is broken unless v ∈ πℤ"; v = π in the window for 1/2 < λ < 1); ch. 42 "Beyond sine-Gordon" (the module relations are no longer open). Ch. 41 not touched (its rank h−1 claim is for 𝕋_V only, correct). No Part VII or index edits needed. No Appendix D entries (all corrections book-internal).
- Open, recorded in the report's Integration section: a proof of the V̄ half of the second module (sketch: the e-type and Cartan equations for V imply them for the fused (T⊗1)(1⊗T), hence for V̄ through the bootstrap embedding; not written out) — DONE 2026-10-13 by a different route (sparse inverse, upgrades integration); the species-r excited defect (rank-1 residue) as the simplest case of problem 10.

## Global pending edits

- **Notation clash ω.** Ch. 39 uses ω = e^{2πi/h}; ch. 41 uses ω = 4π/β² − 1, Part VI's coupling. Decide globally, e.g. by renaming ch. 39's root of unity — DONE 2026-10-12: not renamed. Part III (ch. 10, `eq-an-soliton`) uses ω = e^{2πi/h} and Part VI (ch. 20–21) ω = 4π/β² − 1, and ch. 20 already warns about it, so Part VIII inherits a book-wide double use; ch. 39 now says so in one sentence (and its script docstring), as does the conventions brief §10. (Part VII ch. 31 writes the root of unity ζ instead; a rename to ζ in ch. 39 remains possible.)
- Harmonize theorem titles: Part VII and ch. 40, 42 put a `## Title` inside every `::: {#prp-…}` block; ch. 36, 38, 39 have none and ch. 37 has one — DONE 2026-10-12 (global integration). Titles added to all 28 untitled blocks of ch. 36–39.
- Rewrite Part VII open problem 10 (`sec-boundary-open-problems`, ch. 35): integrable interacting defects are purely transmitting — DONE 2026-10-12 (global integration).
- Naming: the conventions scripts call the two signs of 𝖣 "types I, II"; the book's names are type I₊ and type I₋ (type II is the Corrigan–Zambon defect with an auxiliary field) — not needed: the scripts were already renamed in the conventions commit bc915ba (they say I_-/I_+ throughout); the remaining "type II" strings are the genuine Corrigan–Zambon defects (`defect_classical_an.py` check 10) and Weston's type-II module, and 'I'/'II' in `defect_an_transmission.py` label the crossing relations (I), (II).
- **Part VII edits from the ch. 42 referee integration (2026-10-13).** Do not apply them on the defects branch; apply them on the branch that edits Part VII. The first two were verified with `boundary_common.an_vector`: V(1/y)∘ω_q = V̄(y) to 2e-16, while V(y)∘ω_q = V̄(1/y).
  - `part7-half-line/33-soliton-reflection.qmd` line 72: replace "In the principal gradation, at the physical $q$, $\bar V(y)$ is $V(y)$ composed with $e_j\leftrightarrow f_j$, $k_j\to k_j^{-1}$. It is the antisoliton at the same rapidity," by "In the principal gradation $\bar V(y)$ has $e_j=yF_j$, $f_j=E_j/y$ and $k_j=q^{-H_j}$: it is $V(1/y)$ composed with $e_j\leftrightarrow f_j$, $k_j\to k_j^{-1}$. At the physical $q$ it is the antisoliton at the same rapidity,". — PENDING
  - `part7-half-line/34-quantum-symmetric-pairs.qmd` line 72: replace "the conjugate multiplet: $\bar V=V\circ\omega_q$ (@sec-an-diagonal-k)." by "the conjugate multiplet: $\bar V(y)=V(1/y)\circ\omega_q$ in the principal gradation (@sec-an-diagonal-k)." — PENDING
  - `part7-half-line/35-particle-reflection.qmd` problem 9 (line ~317): replace "Classically, a type-I defect fused onto a boundary of @tbl-durham gives integrable dynamical $a_n^{(1)}$ boundaries with $n$ boundary variables (@prp-fused-bdg); their quantization is open." by "Classically, a type-I defect fused onto a boundary of @tbl-durham gives dynamical $a_n^{(1)}$ boundaries with $n$ boundary variables and a conserved Lax generating function (@prp-fused-bdg); whether their charges are in involution, and their quantization, are open." (`prp-fused-bdg` proves conservation only, not involution.) — PENDING

- **From the upgrades integration (2026-10-13), outside the files that integration could edit.** `prp-dpt-doubled` is now a [Theorem] for type I₋ and every r ≥ 2. (a) Ch. 40 `tbl-boundary-defect-algebra`, status cell of the doubled-boundary column: "$a_2$, $a_3$ [Numerical]" → "$a_r$, type I$_-$ [Theorem]" — PENDING. (b) `index.qmd` l. 51, "or for $a_2^{(1)}$ and $a_3^{(1)}$ by the rotation of the Dynkin diagram": not wrong, but now weaker than the theorem; optionally "or for $a_r^{(1)}$ by the rotation …" — PENDING (optional).

## Global integration (2026-10-12)

Done by the global-integration agent; checkpoint `briefs/global_integration_part8_checkpoint.md`.
- Every pending edit above is marked DONE or "not done/not needed" with a reason. Still open: a reproducible $a_3$ spin-±2 script (ch. 38), and the optional bibliography entries.
- `index.qmd`: Part VIII in the part list, status note, presentation of known material and original results ("Defects", ch. 36–43); node VIII in the dependency map (edges from II, III, VI, VII; none to IX); reading path "Integrable defects"; conventions sentence; figure `intro_f1_reading_paths` regenerated (six paths).
- `outlook/50-outlook.qmd`: three defect connections (spin chains and circuits, AdS/CFT defects, Bäcklund transformations and Q-operators) and pointers to `sec-defect-open-problems`.
- Appendix D: the Part VIII corrections to published papers, with a fifth kind, *clarification*. Not listed there, by the appendix's rules: convention-only items (Bajnok–Simon $B_{\rm BS}=B/2$; BGKNR/BLZ translation; cz2007 (7.3)–(7.5) and Sec. 6 species labels), items that confirm or extend without correcting (cz2007 Sec. 7 remarks, Bajnok–Simon Sec. 4.2, cz2009a Sec. 5, cz2012 (2.18)), and book-internal items. cz2007 Sec. 2 after (2.22) is under "Corrections made earlier by others" (cz2009a noted it).
- Appendix C: section `sec-app-code-defects` with all `defect_*.py` scripts. `textbook_code/README.md`, `briefs/README.md` and the root `README.md` cover Part VIII.
- Conventions brief: ch. 38 and ch. 40 items folded into §5, §6 and §8 (PP sign in §9 corrected); new §10 (ch. 39), §11 (ch. 41), §12 (ch. 42, the report's "§10") and §13 (ch. 43, the report's "§11").

## Integration of the upgrades referee report (2026-10-13)

Report `briefs/referee_defects_upgrades_report.md` (Codex, commit inspected 4c78d75; Integration section at its end); checkpoint `briefs/integration_upgrades_checkpoint.md`; own scripts `briefs/referee_defects_upgrades_scripts/integration_*.py` with outputs. Every item re-derived before editing; all eleven verified, none rejected; no label downgraded.
- Ch. 38 (item 7): `defect_spin_two_a3.py` `boundary_G` docstring (flux G = (T+Θ)(φ) − (T+Θ)(ψ), dQ_bulk/dt = εG, defect contribution −εF, spin −ε·spin); 9 PASS, 151 s. The a₃ [Theorem] of `prp-spin-two` (3) is now refereed.
- Ch. 41 (items 1, 2, 3, 5, 6): uniqueness hypothesis "defined and invertible"; all-rank coefficient check of `eq-dt-reversed`; type-I₊ selection-rule sentence; cz2011 (3.2) misprint (x = e^{2πθ/ξ}), qualified pole-slip statement and normalization at z = 0 (also Notes, Appendix D); Cartan lift and category-O remark (no pinpoint citations: TK1992, FH2015 not available locally). Script check (j) gains the printed-(3.2) control; 28 PASS, 192 s.
- Ch. 42 (item 8): the limits of `prp-diagonal-borel` (5) and the proof are for the matrix part 𝕃̂ = P𝕋̂, 𝕋 = Φ𝕋̂; caption; script and figure-script docstrings; 20 PASS, 23 s (figure script PASS, SVG restored: ids only).
- Ch. 43 (items 9, 10, 11): type-I₊ particle proof in the same variable; **`prp-dpt-doubled` [Numerical] → [Theorem]** (type I₋, every r ≥ 2, generic q), proof complete (V half, V̄ half by the sparse inverse, module relations, failed twists, coideal obstruction); CDD sentence (even h sufficient, not necessary; odd h not sufficient), CZ 2007 after (6.9) cited, novelty qualified; "on each fixed-charge sector" in `tbl-dpt-spectrum` and problem 10; problem 12 and ch. 42 "Beyond sine-Gordon" follow the theorem. Script (h) gains H M_V̄ = 1 − X^h/q, (i) the even/odd-h sinh-type test; 11 PASS, about 1 min.
- Labels changed: only `prp-dpt-doubled`. Its proof is the referee's (item 10), re-derived by hand and by an exact script here; the integration text itself has had no referee pass.
- Reference check: 65 files, 1164 ids (no new ids), no duplicates or unresolved references, all citation keys in `references.bib`.
- Pending: see "Global pending edits" (ch. 40 table cell, optional `index.qmd` l. 51).
