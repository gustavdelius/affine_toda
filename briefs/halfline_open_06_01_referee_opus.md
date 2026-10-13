<!-- Opus referee follow-up, 2026-10-12. Model: Claude Opus 5.5, effort max. Prompt: briefs/halfline_open_06_01_opus_prompt.md. Combined session with 12b (see halfline_open_12b_referee_opus.md). Saved verbatim by the coordinator (section 06_01 of the hand-back). Work directory: /tmp/claude-1000/-home-gustav-Git-affine-toda/0e38a99d-702a-40fb-99b7-74dbe0413cf6/scratchpad/ref_12b_0601/ (W). -->

Headline (from the hand-back): The TeX of BP's Table 2 supports the worker's 12/9 reading, not Gemini's. But computing K as the limit of regular configurations gives `++----` a non-diagonal K, against BP's table, so the count that survives computation is 6/15.

# 06_01

## Item 1 — a_5 diagonal count
**Gemini's reading of the TeX is REFUTED; BP's diagonal entry for `--++++` is CORRECTED by computation.**

**Reading BP's tables.**
- Table 2 (bb1.tex l.1355): `-++-++` and `--++++` share one row with diagonal K: K_1 = (2)(4)/((1)(5)), K_2 = −(1)(3)²(5)/((2)(4)), K_3 = −(2)²(4)²/(3)².
- l.1356: the "not diagonal" row has exactly three entries: `--++-+`, `--+-++`, `--+--+`.
- The five-pattern block Gemini describes exists only in Table 1 (l.1127: same 2 2̄ background, C_BP = √3−1).
- So BP's counts for the 21 missing a_5 patterns are 12 diagonal and 9 non-diagonal, as the worker said:
  - diagonal: the `+-+---` orbit (6) and the `++----` orbit (6);
  - non-diagonal: `++-++-` (3) and the exceptional `+++-+-` (6).

**Computation as a limit of regular configurations.**
- Background: the species-2 stationary pair τ_j = 1 + 2d cos(χ+2πj/3)e^{√3x} + (d²/4)e^{2√3x}, with d = 4−2√3.
- At χ_end = 5π/3, τ_2(0) = τ_5(0) = 0 (simple zeros in x). For χ = χ_end ± δ the members are regular, with exact patterns `+--+--` and `-+--+-` (BC residual 1e-12).
- Method (`ref0601_a5_tauK.py`): BP's own §5 τ-function method. Each linear wave is a species-c "soliton" with σ = ±ip or κ, dressed by the background through Hirota factors. The linearized boundary condition with the target pattern's boundary Hessian is imposed at x=0 in 50-digit arithmetic, so no ODE is solved.
- Checks:
  - each dressed wave solves the linearized Toda equation to 1e-50;
  - the regular patterns reproduce BP/DG's diagonal K to ≤1e-20 for δ ≥ 1e-5;
  - for target patterns at δ=0.1 it agrees to 6 digits in every entry with the worker's independent ODE solver `linrefl2.py` (`ref0601_xcheck.py`).
- Result at θ=0.6 (only the mass-1 pair open):
  - `++----`, converged limit: K = [[0.155633+0.677156i, −0.48997−0.52647i], [−0.21095+0.68756i, same diagonal]]; |K_11| = 0.695, |K_15| = 0.719;
  - same from both sides to 4e-4 at δ=1e-3 and 5e-6 at δ=1e-5; unitary to 1e-20.
  - BP's diagonal value is 0.856552+0.516061i, so `++----` is **not diagonal**.
- Off-diagonal limits of the other three singular patterns: `++--+-` 0.091 or 0.272 and `++-+--` 0.272 or 0.091 depending on side; `++-++-` 0.322.
- At Ω=2.3 (all channels open) the `++----` off-diagonals are 0.0174.
- Precision: below δ≈1e-8, 50 digits are not enough (the control fails), so those rows are discarded.
- Middle-soliton orbit (`ref0601_a5_midK.py`, τ_j = 1+(−1)^j d e^{2x}, d→1):
  - for `+-+---`, `+-----` and `+-+-+-`, the off-diagonals fall linearly in 1−d (1.4e-6 at 1−d = 1e-6);
  - the `+-+---` K_11 → −(1)(5) to 6e-7, the C≡−1 value = DG's b=3 amplitude;
  - this agrees with BP's row 1.
- Why it differs from BP: BP do not take this limit. They perturb the singular background directly and require the perturbation of a vanishing τ_i to vanish at x=0 (l.1198–1216). Their criterion (l.1269–1292) is that K is diagonal iff A_{i−1} = −A_{i+1} at each vanishing τ_i, which `--++++` satisfies. The limit of regular configurations depends on the absolute neighbour signs, not only on their product. So the classical K about a boundary-singular background depends on how the linearization is defined.
- My finite-δ ODE solvers were abandoned: both integration directions are ill-conditioned (error ~ rtol·τ(0)⁻²; the regular control drifts by 0.3 at δ ≤ 0.01). `linrefl2` is reliable only for δ ≳ 0.03.

**Corrected rows for `tbl-missing-vacua`.**
- BP boundary-singular infimum available: a_3 2/2, a_4 5/5, a_5 21/21.
- Classical K diagonal, computed as a limit of regular configurations:
  - a_3: 2 (= DG b=2; established by the worker and Gemini);
  - a_4: 0 (BP; worker's numerics indicative);
  - a_5: 6 (the `+-+---` orbit, = DG b=3). The `++----` orbit (6) is non-diagonal, although BP's Table 2 lists it as diagonal. `++-++-` (3) is non-diagonal (computed). The exceptional orbit (6) is non-diagonal per BP, not recomputed.
- So 6/15. That is Gemini's number, but reached for a different reason.

## Item 2 — honest wording: CORRECTED
**What BP say about these backgrounds.**
- They call them "vacuum solutions"/"background configurations" and mark them "Singular at x=0: yes" in Table 1.
- They accept a singularity at x=0 but reject one at x<0 as having "infinite energy" (l.812–820, 1055).
- Finiteness of the energy is by cancellation: E_bound has a 1/x pole with coefficient ∝ (1−A_i), which vanishes when A_i = 1 (C_i = −1), the "marginally stable" case (l.788–802). In the BC both sides go as 1/x and their coefficients must match (l.984–994).
- They never construct the backgrounds or K as limits of regular configurations. They concede that linearization is unclear at a vanishing τ (l.1198–1200) and that "a large subset" of boundary conditions lack consistent classical scattering (l.1286–1288).

**Proposed text, `sec-boundary-missing-charges`, last paragraph.**
> "Bowcock and Perkins supply the missing vacua as energy infima [@bowcock_perkins2003, sec. 4, Table 1]. For $a_3^{(1)}$ the alternating patterns are the endpoint $d=1$ of @prp-self-conjugate-vacuum, where $\tau_1(0)=\tau_3(0)=0$ and $\vec\phi$ diverges logarithmically at the wall. The infimum $E=-8\,m/\beta^2$ saturates the bound $E_P\ge-4\#\{j:C_j=+1\}$ [Established] and is approached by regular configurations but not attained [Numerical]. The reflection matrix obtained as the limit of those about the regular members $d<1$ is diagonal and equals Delius and Gandenberger's $b=2$ amplitude at $B\to0$ [Numerical]; this is consistency with a limit, not a test against a vacuum, since perturbation theory about a field that diverges at the wall is not defined here. For $a_4^{(1)}$ and $a_5^{(1)}$ the infima are boundary-singular stationary pairs or 'exceptional' multi-soliton backgrounds. Bowcock and Perkins's classical reflection factors about them [Table 2] rest on a prescription for perturbations of a vanishing $\tau_j$; computed instead as limits of regular configurations, the $a_5^{(1)}$ pattern $(+,+,-,-,-,-)$ is not diagonal, although their Table 2 lists it as diagonal. The complex bound states of @cnj-missing-charge are not needed for the energy infima."

**`tbl-boundary-types`, last row.**
- Vacuum: "[Numerical] energy infimum, not attained: boundary-singular (τ_j(0)=0 on nodes with C_j=−1) [BP]; φ=0 for C≡−1, h odd".
- Classical K: "no unique definition; as a limit of regular configurations, diagonal (C≡−1 value) only at the middle-soliton endpoint; non-diagonal at the singular pair endpoints computed".

**`sec-real-coupling-vacua`, closing paragraph.**
> "For non-solitonic patterns no regular static solution is known except $\vec\phi=0$ for $C\equiv-1$ with $h$ odd. Bowcock and Perkins give for every pattern of $a_2^{(1)}$–$a_5^{(1)}$ a static Hirota configuration with $\tau_j(0)=0$ on some nodes with $C_j=-1$; its field diverges at the wall, its finite energy comes from a cancellation between bulk and boundary, and it is an infimum approached by regular configurations, not attained [Numerical]. It is a limit of regular static solutions of other (solitonic) patterns, not of its own."

## Item 3 — claim 8: CORRECTED
**What holds.**
- Within BP's static Hirota families the conclusion holds. At all 26 endpoints of the pair families of a_2–a_5, the regular patterns on either side have opposite signs on the two neighbours of each vanishing node, and both neighbours flip across the endpoint (`ref0601_claim8.py`).
- So a ∏C=−1 pattern meeting the endpoint, which has an odd number of flips, is never realized by a regular member.
- The cleaner reason: regular continued-soliton solutions always satisfy solitonic patterns, ∏C=+1 (Delius's ∏A_i = (−1)^{n+1}).

**Why the worker's argument is incomplete.**
- The property is not specific to ∏C=−1. The missing ∏C=+1 patterns that meet an endpoint are not limits of regular solutions of their own pattern either: a_5 `++----` and `++-++-`, and a_3 `+-+-`.
- The "opposite neighbour signs" statement holds only for pair families: the middle-soliton family's regular members (C≡−1) have equal neighbour signs.
- The result also assumes the Hirota family exhausts the regular static solutions (BP's assumption).
- These vacua are limits of regular configurations, just not of solutions.

**Book wording:** "these vacua are limits of regular static solutions of other (solitonic) patterns, never of their own".

## Scripts (06_01 part)
- W/ref0601_a5_tauK.py and .out (exact τ-method, pair endpoint)
- W/ref0601_a5_midK.py and .out (middle soliton)
- W/ref0601_xcheck.py (cross-check against linrefl2)
- W/ref0601_converged.out
- W/ref0601_claim8.py and .out
- Abandoned ODE attempts: W/ref0601_a5_K.py, W/ref0601_a5_Kout.py, W/ref0601_a5_lib.py, W/ref0601_timing.py, W/ref0601_timing2.py, W/ref0601_t3.py, W/ref0601_t4.py, W/ref0601_a5_K.out
- Sources: W/dl_bp/x/bb1.tex
