# Part VIII (defects): drafting status, verification records, pending edits

Working record for Part VIII, kept in the repository because the session scratchpad is ephemeral. Each drafting agent's verification record is copied here when its chapter is integrated. These lists become the "Already verified" and "To be checked" sections of the referee briefs `referee_defects_*.md`.

## Chapter status

| Ch | Drafted | Script(s) | Checks | Refereed |
|---|---|---|---|---|
| 36 | 2026-10-11 | `defect_axioms.py` | 20 PASS (2 s) | no |
| 37 | | `defect_classical_sg.py` (conventions part), `defect_sg_transmission.py` | | |
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

## Global pending edits

- Rewrite Part VII open problem 10 (`sec-boundary-open-problems`, ch. 35): integrable interacting defects are purely transmitting.
- Naming: the conventions scripts call the two signs of 𝖣 "types I, II"; the book's names are type I₊ and type I₋ (type II is the Corrigan–Zambon defect with an auxiliary field).
