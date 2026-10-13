# Ch. 43 drafting checkpoint

Scratch dir: `/tmp/claude-1000/-home-gustav-Git-affine-toda--claude-worktrees-integrable-defects-framework-a54a84/cfddf2dd-7e99-4c1c-b8f5-c799ea33aced/scratchpad/ch43/`

## Status (update after each step)

- [x] Step 0: read conventions, status (all ch. 43 obligations), ch. 35 model, ch. 36/37/39/40/41/42 relevant sections, referee items 10, 17, sources digest (cz2007 §7 (7.1)–(7.5), bajnok_simon2008).
- [x] Step 1: scratch script `breather.py`: a_r lightest breather from ch. 41 matrices (singlet in V⊗V̄, both orders), compare with T^cl_k of eq-an-linear-transmission and cz2007 (7.5).
- [x] Step 2: particle-factor checks (unitarity, crossing, Dorey bootstrap, CDD), sinh-Gordon/BS O(b²) check.
- [x] Step 3: optional (f): a_r doubled-boundary expectation — TESTED, holds (see finding 4).
- [x] Step 4: `textbook_code/defect_particle_transmission.py` written: 8 PASS lines (a, b1, b2, c, d, e, f, g), 120 s wall (f takes ~97 s). Scratch `b6b.py` running a_5 (h = 6) breathers for the record.
- [x] Step 5: chapter written (266 lines), all refs and bib keys exist (python check), moved prp-dpt-doubled before the open-problem list.
- [x] Step 6: report written; final script run 8 PASS in 62 s. DONE.
- a_5 (h = 6) scratch: all five species at ω = 1.37, θ = 0.4+0.1i give T^cl_a to ≤ 5e-9 (`b6b.py`, `b6c.py`; gesdd SVD failed once for a = 4, gesvd fallback now in the script).

## Decisions

- ω notation: write the coupling as ω = 4π/β² − 1 (as ch. 41 and Part VI), and e^{2πi/h} explicitly for ch. 39's root of unity; γ_k = πk/h.
- Bulk bootstrap for transmission has no S-factor (ch. 36 "Bulk bound states"): 𝕋_c(θ) = 𝕋_a(θ+iū_a)𝕋_b(θ−iū_b). a_r fusings: c = a+b (a+b<h), ū_a = πb/h, ū_b = πa/h; c = a+b−h (a+b>h), ū_a = π(h−b)/h, ū_b = π(h−a)/h. By hand: T^cl satisfies both exactly (telescoping cosh identity).
- Breather B_1: pole of S_{1r} at θ_12 = iπ − 2πi/(hω); constituents θ ± i(π/2 − π/(hω)); cz2007 (7.3) puts antisoliton 2 at θ − iū, soliton 1 at θ + iū (their labels).

## Findings (with numbers)

1. **Lightest breathers, a_r, from the ch. 41 matrices** (scratch `breather.py`, dense with windows; `breather_all.py`, exact dict propagation, no window). Breather constituent shift ū = π/2 − π/(hω). The singlet exists only for "left constituent at θ−iū, right at θ+iū" (other order: no singlet). Results:
   - V̄(θ−iū)⊗V(θ+iū) (species r left, species 1 right) → factor = T^cl_1(θ−θ_d) (eq-an-linear-transmission, γ_k = πk/h); V(θ−iū)⊗V̄(θ+iū) → T^cl_{h−1}. Residual of proportionality to |N⟩⊗singlet 1e-16; deviation from T^cl ~1e-9 (scalar g numerics). a_2, a_3, ω = 1.37, R2 = 14.
   - All species: B^(a)_1 = (species h−a at θ−iū) ⊗ (species a at θ+iū), each built from vector solitons at θ_c + iπ(2m−n+1)/h: the U_q singlet is unique in the weight-zero subspace of V^{⊗h}; transmitted factor = T^cl_a(θ−θ_d) for every a, h = 3, 4, 5, ω = 0.6, 1.37, 3.1, θ real and complex, generic offset s; deviations ≤ 2e-8, proportionality residual ≤ 3e-14. Runtime 16 s total (h=5 dominates via scalar factors).
   - cz2007 (7.5) in book labels: their T^{b1}_1 = T^cl_2 = T^cl_{h−1}, their T^{b1}_2 = T^cl_1 (by hand; consistent with ch. 39 item 5 CZ's T_k = T^cl_{h−k}); their (7.3) ordering (their species 2 left at θ−iū, species 1 right) = book V(θ−iū)⊗V̄(θ+iū) since their species 1 = book species 2.
2. **Particle factors** (by hand): T^cl_a satisfies unitarity+crossing T_a(θ)T_ā(θ+iπ)=1, |T_a|=1 on real θ, Hermitian analyticity 𝕋̄_a(θ) = conj 𝕋_a(−θ*) for real θ_d, and the bootstrap at both a_r fusings (c = a+b, ū = (πb/h, πa/h); c = a+b−h, ū = (π(h−b)/h, π(h−a)/h)) exactly (cosh telescoping). Bootstrap for transmission has no S-factor.
3. **Bajnok–Simon O(b²)** (0712.4292 Sec. 4.2): T−(θ, Bϑ(μ)) = T−(θ, μ) + (μb²/8π)/(1 − i sinh(θ+μ)); this equals the expansion of the exact T with Bϑ = μ(1 − b²/8π) since ∂_x T−(θ,x) = −1/(1 − i sinh(θ+x)) (by hand). Their classical T∓(k) = −(m sinh μ ∓ iω)/(m cosh μ − ik). Dictionary: μ = −θ_d^cl, B_BS ϑ = −θ_d.
4. **(f) a_r doubled-boundary expectation, a_2 type I₋** (scratch `doubled.py`): with θ-independent unknown Y (weight α_j), κ on 𝒲, solved jointly for V and V̄ at 3 rapidities (R2=18): e-generators twisted by π (ψ(e_j) = e_{j−1}) hold [2.6e-15]; untwisted e and π^{-1}-twisted e fail [0.93]; f untwisted holds [6.5e-15] (theorem: Δ(U_q(b⁻))); f twisted by π^{±1} fails [0.85, 0.88]. Equations: e: T[π(e_j)⊗1 + π(k_j)⊗Y] = [Y⊗1 + κ⊗π(e_{ψ(j)})]T; f: T[π(f_j)⊗κ + 1⊗Y] = [Y⊗π(k_{ψ(j)})^{-1} + 1⊗π(f_{ψ(j)})]T (from Δ_A of eq-dbd-coproduct, block P𝕋). Weights: e: Y α_j, κ α_j−α_{ψ(j)}; f: Y −α_{ψ(j)}, κ α_j−α_{ψ(j)}. All j for a_2 and a_3 (scratch, 150 s): same pattern, holds 3e-15, fails 0.67–0.93. The e/π solution has a null space (free constants, as for sine-Gordon in ch. 42); the module relations of the second structure ρ̃ for a_r were NOT checked.
5. **Real-coupling vacuum** (by hand + script g): B_d^real ≥ 0 for σ > 0 by convexity (Σ_j α_j = 0 ⇒ Σ_j log f_j = 0 ⇒ Σ(f_j − 1) ≥ 0, same for g_j), equality iff φ = ψ = 0; unbounded below for σ < 0.
6. **Linearization** (script d): with an_data's π eigenvector of eigenvalue e^{2πik/h}, D = 1 + i cot γ_k ✓, R = 0, T = T^cl_k (I₋) and T^cl_{+,k} (I₊), r = 2, 3, 4.

## Chapter plan (sections, new ids)

New ids: `eq-dpt-particle`, `prp-dpt-vacuum` [Theorem], `prp-dpt-doubled` [Numerical], `tbl-dpt-spectrum`, exercises `exr-dpt-telescope`, `exr-dpt-cz`, `exr-dpt-vacuum`, `exr-dpt-bs`, `exr-dpt-limits`, `exr-dpt-imaginary`, `exr-dpt-cdd`, `exr-dpt-doubled`. Open problems: ~14 grouped items covering every ch. 43 obligation.

## Dead ends

(none yet)
