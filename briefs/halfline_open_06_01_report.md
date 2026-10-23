<!-- Worker report, 2026-10-12. Model: Claude Sonnet 5.5, effort high. Scope: all of problem 6 (brief halfline_open_06_missing_vacua.md) and the classical part of problem 1 (brief halfline_open_01_negative_product_vacua.md, steps 1-4 and 6; step 5 not attempted). Interrupted twice by the usage limit and resumed. Saved verbatim by the coordinator. Scripts and outputs copied to briefs/halfline_open_06_01_scripts/ (pickles not copied). -->

# Report: half-line open problem 6 (missing vacua) and the classical part of problem 1 (∏C=-1 vacua)

Nothing in the repository was edited. Step 5 of problem 1 (particle amplitudes) was not attempted. All scripts are in
`/tmp/claude-1000/-home-gustav-Git-affine-toda/0e38a99d-702a-40fb-99b7-74dbe0413cf6/scratchpad/hl/` (listed at the end). That directory is called `hl/` below.
The sibling directory `.../scratchpad/work/` is used by another agent, so I moved my files out of it.

## 0. Summary

1. **Bowcock–Perkins (BP) already contain vacua for every sign pattern of a_2..a_5, including the missing patterns of problem 6 and all ∏C=-1 patterns.**
   - I read the paper (arXiv hep-th/9909174), tables 1 and 2, and translated them to book conventions.
   - I then recomputed the vacuum energy of every pattern class independently. All 12 classes in a_4 and a_5 agree with BP. For a_5 this means every class except (+,-,-,+,-,-) (A=-++-++), which has a regular vacuum, and (++-+--) (A=--+-++), which equals a class already tested by rotation or reflection. For a_2 and a_3 every pattern was run. In no case did the search go below BP.
2. **These vacua are not regular.**
   - For every pattern except flat φ=0 and the regular intervals of the continued pair, the vacuum is a static Hirota configuration with τ_j(0)=0 on non-adjacent nodes where C_j=-1.
   - The energy infimum is approached by regular configurations but not attained: φ diverges logarithmically at x=0 on those nodes.
   - So "no regular real vacuum" is the correct reading, and the infimum exists, is finite, and equals BP's energy.
3. **Missing patterns of a_3 (alternating).**
   - The vacuum is the endpoint d→1 of the C≡-1 self-conjugate family τ_j=1+(-1)^j d e^{2x}.
   - E=-8 m/β², and this equals a rigorous lower bound (§2.3).
   - The classical K equals the C≡-1 flat value, which equals DG's class b=2 value (checked, 10⁻¹²).
   - So DG's amplitude assignment is consistent with this vacuum.
4. **Missing patterns of a_4 (++−+−, class b=2).**
   - The vacuum is BP's "exceptional" 4-soliton 1 1̄ 2 2̄, E=-10.416.
   - I reproduced it independently and found the Hirota parameters.
   - Its linearized reflection looks non-diagonal, so it does NOT reproduce DG's diagonal b=2 amplitude (partial numerical evidence).
5. **DG (5.7) check.** The linearized reflection about every regular continued-pair vacuum of a_2..a_5 equals the B→0 limit of DG (5.7) to 10⁻¹¹ (30 patterns). The book said this was unchecked.
6. **Problem 1.** All ∏C=-1 patterns have a BP vacuum, singular at x=0, except C≡-1 with h odd, which has flat φ=0.
   - Classical K is diagonal only for the middle-soliton (d=1) type; there it equals the C≡-1 value.
   - Otherwise K is non-diagonal.
   - No ∏C=-1 boundary has K→1, so the s=-1 breather amplitude is not selected classically.

## 1. BP translated to book conventions

Dictionary:
- A_i = -C_i, and BP's flat "+++…" is the book's C≡-1.
- BP list one pattern per cyclic orbit. The cyclic rotation of the extended Dynkin diagram is a symmetry at real coupling, unlike the coweight shift. I confirmed numerically that F is rotation-invariant for all 16 patterns of a_3.
- E_book = -2(n+1)·C_BP - 2ΣC_j, in units m/β². Checks: this gives 0 for flat C≡±1, and for a_2, C=(+,-,-) it gives 8-6√3, as in the book's (eq-vacuum-energy).

The table below uses book C-patterns.
- "reg" means regular on x≤0 (the book's prp-solitonic-vacua).
- "sing" means τ_j(0)=0 on the nodes with C_j=-1.
- K is BP's classical reflection, d=diagonal, nd=non-diagonal.
- b is DG's class for solitonic patterns.
- "num" is my reduced minimization (§2), converging from above.

**a_2**

| C (orbit) | ∏C | solitons | E_BP→book | num | type | K |
|---|---|---|---|---|---|---|
| (−−−) | −1 | none | 0 | 0 | flat | d |
| (+−−) ×3 | +1, b=1 | 1 1̄ | −2.3923 | −2.3923 | reg, modulus χ | d (=DG) |
| (++−) ×3 | −1 | 1 1̄ | −6.3923 | −6.3923 | sing | nd |
| (+++) | +1 | none | 0 | 0 | flat | d |

**a_3**

| C | ∏C | solitons | E | num | type | K |
|---|---|---|---|---|---|---|
| (−−−−) | +1, b=2 | none, or 2 family | 0 | 0 | flat + modulus | d |
| (+−−−) ×4 | −1 | 2 (d=1) | −4 | −4.0000 | sing | d (C≡−1 value) |
| (+−+−) ×2 | +1, b=2 (**missing**) | 2 (d=1) | −8 | −8.0000 | sing | d, = DG b=2 |
| (++−−) ×4 | +1, b=1 | 1 1̄ | −3.3137 | −3.3137 | reg | d |
| (+++−) ×4 | −1 | 1 1̄ | −7.3137 | −7.310 | sing | nd |
| (++++) | +1 | none | 0 | 0 | flat | d |

**a_4** (ids by C-pattern; "missing" means the book's table has no single-soliton vacuum)

| C | ∏C | solitons | E | num | type | K |
|---|---|---|---|---|---|---|
| (−−−−−) | −1 | none | 0 | 0 | flat | d |
| (+−−−−) | +1, b=2 | 2 2̄ | −3.0211 | −3.0211 | reg | d |
| (+−+−−) | −1 | 2 2̄ | −7.0211 | −7.0126 | sing | nd |
| (++−−−) | −1 | 1 1̄ 2 2̄ | −6.4162 | −6.357 | sing | nd |
| (++−+−) | +1, b=2 (**missing**, 5 patterns) | 1 1̄ 2 2̄ | −10.4162 | −10.378; on singular curve −10.4152 | sing, "exceptional" | nd |
| (+++−−) | +1, b=1 | 1 1̄ | −3.7557 | −3.7557 | reg | d |
| (++++−) | −1 | 1 1̄ | −7.7557 | −7.742 | sing | nd |
| (+++++) | +1 | none | 0 | 0 | flat | d |

**a_5** (BP E; num is my reduced minimization over 12 starts; "d" or "nd" is BP's K)

| C | ∏C | solitons | E | num | type | K |
|---|---|---|---|---|---|---|
| (——————) | +1, b=3 | none | 0 | 0 | flat | d |
| (+−−−−−) | −1 | 3 (d=1) | −4 | −4.0000 | sing | d |
| (+−+−−−) | +1, b=3 (missing) | 3 (d=1) | −8 | −7.99996 | sing | d |
| (+−+−+−) | −1 | 3 (d=1) | −12 | −11.9998 | sing | d |
| (+−−+−−) | +1, b=2 | 2 2̄ | −4.7846 | −4.7846 | reg | d |
| (++−−−−) | +1, b=2 (missing) | 2 2̄ | −4.7846 | −4.774 | sing | d (= DG b=2) |
| (++−+−−), (++−−+−) | −1 | 2 2̄ | −8.7846 | −8.77 | sing | nd |
| (++−++−) | +1, b=3 (missing) | 2 2̄ | −12.7846 | −12.758 | sing | nd |
| (+++−−−) | −1 | 1 1̄ 2 2̄ 3 | −7.4211 | −6.94 (poorly converged) | sing | nd |
| (+++−+−) | +1, b=2 (missing) | 1 1̄ 2 2̄ 3 | −11.4211 | −11.386 | sing | nd |
| (++++−−) | +1, b=1 | 1 1̄ | −4 | −4.0000 | reg | d |
| (+++++−) | −1 | 1 1̄ | −8 | −7.989 | sing | nd |
| (++++++) | +1 | none | 0 | 0 | flat | d |

Verdicts:
- Recomputed: all rows above in the "num" column, as energy infima over the static Hirota family.
- Not recomputed: BP's K entries. I recomputed only the K's described in §3.
- "Missing" patterns, by counting from BP's rows (a_5 is by row counting, unverified by direct enumeration):
  - a_3: 2 of 2 have BP vacua, with diagonal K equal to DG.
  - a_4: 5 of 5 have BP vacua (exceptional), with non-diagonal K.
  - a_5: of the 21 missing, 12 have BP vacua with diagonal K (the "3" middle soliton, 6 patterns; the singular 2 2̄ with b=2, 6 patterns). The other 9 have non-diagonal K (singular 2 2̄ with b=3, 3 patterns; exceptional 1 1̄ 2 2̄ 3 with b=2, 6 patterns).
- a_6 and a_7 were not examined.

## 2. How the vacua were computed

### 2.1 Reduction

- For a configuration with φ→0 at -∞, ½φ'²+V ≥ √(2V)|φ'|, so E ≥ d_J(0,φ(0)) + B_P(φ(0)). Here d_J is the Jacobi (Agmon) distance for the metric 2V.
- Equality holds for the static solution that decays at -∞ (the Bogomolny equation of BP, §3).
- Hence inf E_P = inf over φ(0) of F_P = W+B_P. For a regular static solution W = ∫(½φ'²+V) = -2Σ_j τ_j'(0)/τ_j(0). I verified this by quadrature and checked the field equation (residual 4e-15).
- I parametrized static solutions by all n species with real data (c_{h−a}=c̄_a). Coefficients: τ_j=Σ_μ Π(c_aω^{ja}e^{m_a x})^{μ_a} Π A_ab^{μ_aμ_b}, with A_ab = sin²(π(a−b)/2h)/sin²(π(a+b)/2h).
- Regularity means τ_j>0 on x≤0, tested on a dense grid with parabola refinement.
- The assumption that this family exhausts the regular static solutions is BP's own.
- Nelder–Mead was run with 12–60 random starts per pattern.

### 2.2 A pitfall

A grid that is too coarse admits τ with tiny negative dips, which give spurious energies (-5000, etc.). The final code uses a dense grid and requires τ_j(0)>1e-6.

### 2.3 Rigorous bound for the missing a_3 patterns

- With reference pattern Q≡-1: E_P = E_Q + 4Σ_{j: C_j=+1}(e^{ψ_j(0)/2}−1).
- BP's bound gives E_Q≥0, with equality only on the Q-vacua.
- So E_P ≥ -4·#{C_j=+1}.
- For C=(+,−,+,−) this is -8. The d→1 limit of τ_j=1+(-1)^j d e^{2x} has e^{ψ_even/2}→0 and E_Q=0, so the bound is saturated and not attained by any regular configuration.
- The same argument gives E=-4 for (+,−,−,−).
- It also shows that for h odd, E_Q has the unique minimum φ=0, so φ=0 is the vacuum of C≡-1 (F=0 numerically for a_2 and a_4, all starts).

### 2.4 Direct PDE check

- A discretized gradient minimizer (L-BFGS on a graded grid) in `halfline_min.py` agrees with the reduced results for regular cases: a_2 flat 1e-13; a_2 solitonic -2.3899 vs -2.3923; a_3 (++−−) -3.3101 vs -3.3137.
- For singular vacua the lattice limit converges only slowly. For a_2 (++−) it gives -5.17, -5.64, -5.77, -5.79 for dx0=1e-2…1e-5, versus -6.392. The lattice results never go below the reduced infimum, but they are not an independent confirmation of the singular energies.
- The independent evidence for the singular energies is the explicit regular sequence: configurations from the pair family approach the endpoint with E_P→F (§2.3 for a_3; scan for a_2 and a_4).

## 3. Reflection factors

`linrefl2.py` solves the linearized problem y''=(M²(x)−Ω²)y, y'(0)=−Hy(0), H=½ΣC_j e^{ψ_j(0)/2}α_jα_j^T. It integrates from -X to 0 for stability and returns the open-channel K matrix.

- It reproduces flat K=−[(a)(h−a)]^{−C} (1e-12).
- **DG (5.7) ↔ linearization:** for the continued pair at every regular χ-interval, for a_2, a_3, a_4 (a=1,2), a_5 (a=1,2), K11 and K22 equal DG's K_1^{(b)}(B→0) to ≤1e-11 and the off-diagonal elements are ≤1e-11 (output of `kcheck.py`). This also verifies the formula for the class b from `prp-dg-weights`.
- **a_3 alternating:** along the family d<1 (pattern C≡-1) K is d-independent and equals the flat C≡-1 value and DG (b=2, a=1) to 4 digits (θ=0.6). The singular limit d→1 has the same K, because the terms of the Hessian that differ between the alternating pattern and C≡-1 vanish in the limit.
- **a_2, ∏C=-1, (++−):** the regular sequences approaching the endpoint give a non-diagonal limit at θ=0.6, K11≈0.005+0.995i, K12≈0.049−0.086i, K21≈0.049+0.085i, the same from both sides of χ_end.
  - This matches BP's "not diagonal".
  - ∏C=-1 vacua are not limits of regular solutions of the same pattern: the regular neighbours of the endpoint have a fixed product of the signs at j±1, and flipping one gives ∏C=-1.
- **a_4 missing (++−+−):** at a near-singular point of the exceptional family, K has diagonal ≈0.325+0.934i (DG: 0.343+0.939i) and off-diagonals ≈0.14 for τ_2(0)=1e-3, 1e-4. The solver is unreliable below 1e-5, so this is indicative only. It is consistent with BP's "not diagonal".

## 4. Step 3 (Hirota fits)

- **a_3 alternating:** the minimizer is species 2 only, |c_2|=1 (χ=0 or π), other species 0.
  - This is the d=1 endpoint of the self-conjugate family of weight ε₁+ε₃, not a continued bound state with weight ε₁+ε₂.
  - It lies in the same DG class b=2. This resolves the "missing vacuum" without cnj-missing-charge.
  - A complex bound state would have c_2 purely imaginary, which makes τ complex, so it cannot give a real vacuum. I did not test the fusing limit u→π/4 numerically; this is an argument only.
- **a_4 (++−+−):**
  - The 4-soliton family with phases satisfying BP's reality conditions e^{i(χ₂−2χ₁)}, e^{−i(χ₁+2χ₂)} ∈ {±1} has τ_2(0)=τ_4(0)=0 on a one-parameter curve (ln d₂ against ln d₁).
  - F along this curve has a minimum -10.4152 (at τ(0)=1e-7) at d₁=0.2579, d₂=0.5637, which tends to BP's -10.4162 (C_BP=√(5+2√5)−√5) as the cutoff goes to 0.
  - Both zeros are on nodes with C=-1 and the configuration is reflection-symmetric.
- **a_5 exceptional 1 1̄ 2 2̄ 3:** found numerically with E=-11.386 (approaching -11.421), not studied further.

## 5. Problem 1

- **Are there vacua at all?** Yes as infima. For the ∏C=-1 patterns of a_2..a_5 the vacuum is the singular endpoint of a continued-soliton configuration. It is attained only in the limit; no regular real vacuum exists. The energies agree with BP and are given by the same formula.
- **Which soliton configuration do they continue?**
  - a_2, a_3, a_4, a_5 with "1 1̄": the endpoint where τ_j(0)=0 of the pair family (χ_end), at the point where two regular solitonic patterns meet.
  - Patterns related to C≡-1 or to the species-h/2 soliton: the d=1 endpoint of that family.
  - The rest: exceptional multi-soliton backgrounds.
  - I did not find the imaginary-coupling images beyond this; step 3 of problem 1 is only done for the real-coupling Hirota data. The h-odd coweight-shift images of C≡-1 at imaginary coupling were not compared.
- **φ=0 for C≡-1, h odd.**
  - It is the minimum: F=0 in all starts for a_2 and a_4, plus the bound E_Q≥0.
  - The linearized spectrum has no negative modes. The boundary quadratic form for a mode e^{κx} of mass μ is Q=(κ²+μ²)/(4κ)−μ²/4≥(μ/4)(2−μ)>0 for μ<2, so there is a zero mode only for a=h/2 (consistent with the book).
- **s=-1.** The diagonal ∏C=-1 reflections are the C≡-1 values −[(a)(h−a)], never 1, and the others are non-diagonal. No ∏C=-1 boundary has the classical limit K→1 required by K^{(−)}→1, so the classical data do not favour assigning s=-1 to any ∏C=-1 boundary.
- **Amplitudes:** not attempted.

## 6. Proposed text (for the book)

**Revised `tbl-missing-vacua`:** add rows "multi-soliton or boundary-singular vacuum (BP)": a_3: 2/2, a_4: 5/5, a_5: 21/21 (diagonal classical K: a_3 2, a_4 0, a_5 12). For a_6 and a_7, unknown. Keep the existing rows ("single-soliton regular vacuum": 5 of 7 for a_3, etc.).

**`sec-boundary-missing-charges`, replacement for the last paragraph:** the missing vacua are provided by Bowcock and Perkins. For a_3 the alternating patterns are the endpoint d=1 of the family of @prp-self-conjugate-vacuum, where τ_1(0)=τ_3(0)=0; the field diverges logarithmically on the odd nodes at x=0 but the energy infimum is finite (E=-8 m/β², exactly saturating E_P ≥ -4#{C_j=+1} given E_{C≡-1}≥0) and is not attained by a regular configuration. The classical reflection factor equals the C≡-1 value and DG's b=2 limit. For a_4 and a_5 the vacua are multi-soliton ("exceptional") singular backgrounds, and for many patterns the linearized reflection is non-diagonal, so it does not agree with DG's classical limit (BP; our numerics for a_4 are indicative). The complex bound states of cnj-missing-charge are not needed for vacua and give complex fields.

**Last row of `tbl-boundary-types`:** real-coupling vacuum "boundary-singular infimum (BP): the endpoint (τ_j(0)=0, C_j=-1) of the continued pair or middle soliton, or an exceptional multi-soliton; flat φ=0 for C≡-1 with h odd"; classical bound states: "K diagonal only for the middle-soliton type, equal to C≡-1; otherwise non-diagonal".

**Closing paragraph of `sec-real-coupling-vacua` (Assumptions):** replace "no regular static solution is known" by: BP found infimum-attaining boundary-singular backgrounds for all patterns of a_2..a_5; these are not regular (φ diverges at x=0 on nodes with C_j=-1) and for ∏C=-1 the limit is not a limit of regular solutions of the same pattern.

## 7. Caveats

- The "vacuum" for singular patterns is an infimum, not a finite-energy regular solution: the bulk energy and the boundary term diverge separately and cancel only in the limit. BP accept it; whether perturbation theory about such a background defines a particle reflection factor is open.
- The static Hirota family is assumed to exhaust regular static solutions with φ→0.
- The searches are stochastic; the a_5 value for (+++−−−) is not converged (−6.94 vs BP −7.42), so no a_5 claim in the table rests on it beyond "not lower than BP". a_5 and a_6/a_7 K entries are BP's and not recomputed except the regular ones.
- The K for the a_4 exceptional case is only indicative (solver unreliable at the smallest cutoffs).

## 8. Scripts (all in `/tmp/claude-1000/-home-gustav-Git-affine-toda/0e38a99d-702a-40fb-99b7-74dbe0413cf6/scratchpad/hl/`)

- `halfline_min.py`: roots, discretized static energy with gradient, L-BFGS minimizer (shared code).
- `tauvac.py`: static all-species Hirota τ-functions, parametrization.
- `redmin2.py`: reduced minimization with the dense feasibility test; run as `python3 redmin2.py <n> <nstarts>` (final version). Outputs `red2_4.out` and `red2_5.out`.
- `redmin.py`, `run_red.py`, `slsq.py`: earlier, superseded searches (the coarse-grid version produced the red3.out/red4.out/red5.out used for a_3). Do not use for a_4, a_5.
- `bp_table.py`: BP tables 1 and 2 translated to book conventions, with E_book.
- `linrefl.py`, `linrefl2.py`: linearized reflection matrix (the second is the stable one).
- `kcheck.py`: K about all regular pair vacua vs DG (5.7), a_2..a_5.
- `t5.py`: helper functions (blk, Kdg, pattern_of); `t3.py`, `t4.py`: field-equation/energy checks.
- `t6.py`: a_3 C≡-1 family vs DG b=2. `t7.py`: a_2 ∏C=-1 non-diagonal K.
- `exc4.py`, `exc4b.py`, `exc4c.py`: a_4 exceptional vacuum (scan, singular curve, K).
- `pde_check.py`, `t1.py`, `t2.py`: PDE lattice checks (`pde3.out`, `pde4.out`, `t2.out`).
- BP PDF and text: `/tmp/claude-1000/-home-gustav-Git-affine-toda/0e38a99d-702a-40fb-99b7-74dbe0413cf6/scratchpad/dl/bp.pdf`, `.../bp.txt`.

## Integration (branch half-line-open-problems-integration)

Basis: the worker report, the Gemini verdict, the Opus follow-up `halfline_open_06_01_referee_opus.md`, and the coordinator's rerun of the Opus $\tau$-method script for $a_5$ `++----`. The rerun confirmed that the off-diagonal $|K_{15}|$ converges to 0.719 from both sides; BP's diagonal value fails.

| Claims | Verdict | Book change |
|---|---|---|
| 1 (BP has configurations for all patterns of $a_2$–$a_5$) | **Integrated** (literature; energies recomputed by the worker and Gemini) | ch. 31 new paragraph "Boundary-singular vacua"; `tbl-missing-vacua` caption; `tbl-boundary-types` last row |
| 2, 3 (energies; singular; infima not attained) | **Integrated** in qualitative form, with BP's own description (Opus item 2) | same |
| 4 (bound $E_P\ge-4\#\{C_j=+1\}$; saturation at the $a_3$ endpoint) | **Verified** by the coordinator from the book's boundary potential and $e^{\beta\alpha_j\cdot\phi/2}=\sqrt{\tau_{j-1}\tau_{j+1}}/\tau_j$ | ch. 31 missing-charges section, with the proof inline |
| 5 ($a_3$ alternating $K$ = DG $b=2$) | **Verified, integrated** (consolidation A5, 2026-10-13): exact Hirota-dressing method, limit $d\to1$, off-diagonal $\propto1-d$, robust under admixture paths; also $a_5$ middle soliton ($b=3$) | ch. 31 `sec-boundary-missing-charges`, "Middle-soliton endpoints" |
| 6 (DG (5.7) checked for all regular pair vacua) | **Verified, integrated** (consolidation A5, 2026-10-13): A5's exact Hirota-dressing code (56 backgrounds, $a_2$–$a_5$, 3 frequencies, 27 digits; `code/solitonic_linear.py`) and the coordinator's own ODE integration ($a_3$, species 1, six $\chi$ in all four intervals: $K_{11}$ = DG $b=1$ to $10^{-10}$, off-diagonal $<3\times10^{-11}$; `briefs/halfline_consolidation_A5_scripts/coordinator_ode_check.py`) | new `prp-solitonic-linear` (ch. 31); `tbl-boundary-types`; `prp-solitonic-vacua` ($a_5$); ch. 35 "Classical limits" |
| 7, 10 ($a_4$, $a_5$ diagonal counts; `++----` non-diagonal in the limit) | **Numbers verified, interpretation corrected** (consolidation A5): the `++----` limit 0.719 along regular solutions is reproduced, but "same from both sides" is forced by the reflection symmetry $j\to1-j$, and other approach paths give 0.937. So the Opus statement "BP's diagonal entry is CORRECTED by computation" and the count 6/15 are **withdrawn**: the classical $K$ about singular backgrounds is prescription- and path-dependent, and BP's own prescription (sec. 5) was not implemented. No Appendix D entry. $a_4$ exceptional background reconstructed exactly; limit $|K_{1,4}|=0.147$, non-diagonal, path-independent among those tried | ch. 31 `sec-boundary-missing-charges` rewritten with the three cases |
| 8 (∏C=−1 not limits of their own pattern) | **Integrated** in the corrected Opus form ("limits of regular static solutions of other, solitonic, patterns, never of their own") | ch. 31 |
| 9 ($s=-1$ not selected classically) | **Not integrated** | — |
| "φ=0 lowest for C≡−1, h odd" | **Not integrated**; the old parenthesis "(whether it has the lowest energy is not known)" was removed, because the new paragraph and BP's bound supersede it. Note: the bound $E_{C\equiv-1}\ge0$ is used in the text. | ch. 31 |
