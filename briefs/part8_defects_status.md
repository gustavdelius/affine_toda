# Part VIII (defects): drafting status, verification records, pending edits

Working record for Part VIII, kept in the repository because the session scratchpad is ephemeral. Each drafting agent's verification record is copied here when its chapter is integrated. These lists become the "Already verified" and "To be checked" sections of the referee briefs `referee_defects_*.md`.

## Chapter status

| Ch | Drafted | Script(s) | Checks | Refereed |
|---|---|---|---|---|
| 36 | 2026-10-11 | `defect_axioms.py` | 20 PASS (2 s) | no |
| 37 | 2026-10-11 | `defect_classical_sg.py`, `defect_sg_transmission.py` | 27 + 24 PASS (18 s, 103 s) | no |
| 38 | | `defect_classical_an.py` (conventions part), `defect_lax_rmatrix.py` | | |
| 39 | | `defect_solitons_an.py` | | |
| 40 | | `defect_universal_r.py`, `defect_transfer.py` | | |
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
- Hermitian analyticity ⇔ Re b = ±3;
- the kink-rescaling reading of the CCK factor;
- Bajnok–Holló–Watts μ̂ + μ̄̂ = h(b) of ch. 28;
- the signs and arguments of the doubled K; 𝕋 depending on θ−θ_d and 𝕋̄ on θ+θ_d.

**Unchecked:**
- CAFG's non-diagonal argument (their (22)–(24));
- the defect TBA and Lüscher formulas (Bajnok–Simon (25), BHW (94)–(96));
- Bajnok–Simon's |2⁺⟩ state and its Coleman–Thun closure;
- the BHW TCSA conclusions;
- Mintchev–Ragoucy–Sorba (cited generally);
- Jiang and Kormos–Runkel–Watts (cited by abstract);
- Runkel's T-system statements;
- "no interacting example known" for the CAFG loophole;
- that the defect operator, as a trace over 𝒲, commutes with the transfer matrix (left to `defect_transfer.py`, ch. 40);
- the remarks on unstable states.

**Corrections to sources recorded in the Notes:**
- DMS (2.12) needs no parity, and their (2.10) is empty only if S_ab = S_ba.
- The DMS χ branch.
- CAFG (32) omits a factor S_ba.
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
- No poles of Φ in the physical strip.
- The rational form of 𝕋^cl.
- The momentum criterion forcing exponential potentials (F‴/F′ = G‴/G′); the boost argument.
- The KL dictionary ν ↔ −λθ_d.
- Penati–Polvara normalization γ/2 = β_SG²/8π.
- The classical limit of Bajnok–Simon's ε_def; E_D as the continuation of `eq-bsg-dirichlet`.
- 𝕋̄𝕋 = (cosh θ_d + i sinh θ)/(cosh θ_d − i sinh θ).
- The stable-pole window for imaginary θ_d.

**Unchecked.**
- **Type II, quoted only:** all CZ 2009b/2010 type-II statements. In particular their resonance at η_p + iπ/2γ has the opposite sign to BCZ (4.1) and is not reconciled.
- **Attributions to sources:**
  - CZ 2004 gates.
  - Weston's identification, including his Borel subalgebra in book conventions.
  - BCZ Sec. 6 (moving defects) and (3.52)–(3.56).
  - Bajnok–Simon (17), (18), (23).
  - Penati–Polvara's anomaly statements.
  - KL (5.30)–(5.37).
  - The necessity direction of the BCZ classification.
- **Unverified checks:** Bianchi permutability for a general seed; the interpretation of the 𝕋_{B₂} poles.
- **Citations:** `habibullin_kundu2008` is cited by title only.

**Corrections to sources recorded in the Notes.**
- BCZ 2005 (3.42), (3.47) diverge; the KL pairing (3.44) converges.
- BCZ 2005 Sec. 2 two-defect remark: the left-to-right factors are inverse.
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

## Global pending edits

- Rewrite Part VII open problem 10 (`sec-boundary-open-problems`, ch. 35): integrable interacting defects are purely transmitting.
- Naming: the conventions scripts call the two signs of 𝖣 "types I, II"; the book's names are type I₊ and type I₋ (type II is the Corrigan–Zambon defect with an auxiliary field).
