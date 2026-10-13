<!-- Worker report, 2026-10-12. Model: Claude Opus 5.5, effort xhigh. Scope: full brief halfline_open_11_soliton_preserving_scalar.md, given problem 3's report as a method template. Interrupted three times by the usage limit and resumed. Saved verbatim by the coordinator. Scripts copied to briefs/halfline_open_11_scripts/; they also import unchanged copies of textbook_code/*.py and of gz.py, nk_common.py, rhonp.py from briefs/halfline_open_03_scripts/, which must be copied alongside. -->

## Report: problem 11, scalar factors and breather amplitudes for soliton-preserving reflection in a_n^(1)

**Status.** Steps 1–4 are done for n = 1, 2, 3 (some parts up to n = 4). The checks are numerical, at ω = 0.37, with closed forms fitted and stated below. Step 5 has a verdict. Not done:
- the soliton bootstrap that would fix the CDD factor for n ≥ 2;
- crossing and bulk-bootstrap checks of the continued breather amplitude;
- a check against the text of Regelskis–Vlaar 2018 (I did not open the paper);
- numerical verification of the coideal.

Nothing in the repo was edited.

**Verdict in one paragraph.**
- The soliton-preserving (SP) diagonal matrices are the de Vega–González-Ruiz (DVGR) family. They are consistent: the braid reflection equation, the mixed V/V̄ reflection equation, unitarity and the matrix part of crossing all hold.
- They have Gamma-product scalar factors. For n = 1 the result reproduces Ghoshal–Zamolodchikov Dirichlet exactly, with a CDD factor (−T/2)_I².
- But for n ≥ 2, **no SP matrix is Z_N-symmetric**, while Delius's boundary at a vacuum is.
- The SP boundary reflects the breather **B_1 → B̄_1**, i.e. particle 1 → particle n, never B_1 → B_1.
- So the continued amplitudes are neither Dirichlet (−1) nor Neumann (+1) for n ≥ 2. Their B → 0 limit depends on θ.
- This agrees with problem 3: the Neumann boundary is soliton-conjugating with s = +1. The SP route does **not** give a second derivation of the Neumann amplitudes.

### 1. Which matrix (book conventions, principal gradation)

**Braid-form Ř_VV in closed form** (checked against `rcheck`):
- Ř v_i⊗v_i = v_i⊗v_i.
- Ř v_i⊗v_j = t v_j⊗v_i + r_d v_i⊗v_j, with d = (j−i) mod N.
- t = (Z−1)/(qZ−q⁻¹) and r_d = (q−q⁻¹) z^{N−d}/(qZ−q⁻¹), where Z = z^N.
- Ř depends on j−i only through d, so it is Z_N-symmetric.

**Diagonal solutions** (solved analytically by 2×2 blocks):
- On each block the reflection equation is q-independent. Its general solution is k_i/k_j = y^{N−2d}(ξY−1)/(ξ−Y), with Y = y^N.
- Transitivity across blocks leaves exactly one family: a cyclic block S = {s,…,s+l−1} and one parameter η (ξ = e^{2η}).
- With T = Nω/2 and s = 0, the matrix is
  K_V(θ) = diag(e^{a_j ωθ}) · diag(sinh(η+Tθ) for j∈S, sinh(η−Tθ) for j∉S),
  where a_j = 2j on S and a_j = 2j−N off S.
- Checked: the braid reflection equation holds to 1e-15 for n = 2, 3 and all l.

**Identification.**
- The gradation change K_hom = diag(y^{−2j}) K_prin turns this exactly into DVGR's K_a = e^{u} sinh(ξ+u) for a ≤ l and e^{−u} sinh(ξ−u) for a > l, with u = Tθ and ξ = η. The book's unchecked statement is therefore correct.
- Its Z_N rotations (cyclic intervals S) are the other DVGR orientations. l = 0 gives diag(y^{2j}).
- **Coideal (by hand, not verified numerically).** K intertwines e_j and f_j for every node j except 0 and l. So the Satake data are X = {all nodes except 0, l} and τ: j → l−j, which swaps nodes 0 and l. The finite part is AIII-type, gl(l)⊕gl(N−l).
- This is presumably the diagonal family of Regelskis–Vlaar 2018, but I did not check it against the paper.

**Z_N obstruction (proved by hand).**
- A Z_N-covariant diagonal K needs k_i/k_{i+1} to be constant. The block formula then requires a Möbius function of y^N to equal y^{2−N}, which is impossible for N ≥ 3.
- For N = 2 the covariant points are η = iπ/2 (K ∝ 1, GZ Dirichlet at φ₀ = 0) and K ∝ diag(1,−1).
- Delius's condition at a vacuum v = 0 is Z_N-invariant, because the rotation α_j → α_{j+1} is linear and preserves Λ_W.
- Topological charge conservation (Re φ(0) fixed) forces K to be diagonal.
- **Conclusion:** for n ≥ 2 Delius's boundary at a vacuum has no Z_N-symmetric quantum reflection matrix in V. Every SP matrix breaks Z_N down to the reflection that swaps nodes 0 and l.

**Classical side (by hand).**
- With σ = † and a constant diagonal K, the compatibility condition gives K_{i+1}/K_i = e^{±iβα_j·Re φ(0)}. Hence Re φ(0) ∈ (π/β)Λ_W^∨ and K_+ = e^{iπμ·H}.
- This **contradicts the book's statement in `sec-soliton-preserving`** that coweights outside 2Λ_W^∨ fail. They fail with K_+ = 1, but they work with K_+ = e^{iπμ·H}. This is the a_n analogue of sine-Gordon Dirichlet at φ₀ = π/β. It should be rechecked in `boundary_classical.py`.
- Real-coupling Dirichlet φ(0) = φ₀ is **not classically integrable for n ≥ 2** with a field-independent K, any linear σ, in the fundamental representation. The Cartan term forces σ = −Ad_K on h; the step terms then need ψ_j = λ² at every node; and the product around the cycle gives λ^{2N} = 1, which fails. Sine-Gordon escapes because there E_{α₀} = E_{−α₁}.

### 2. Antisoliton matrix and crossing-unitarity

**Antisoliton matrix.** The mixed reflection equation (Ř_{V̄V}, Ř_{VV̄}) fixes it uniquely (residual 1e-15, n ≤ 4):
- K_V̄(θ) = diag(e^{−a_j ωθ}) · diag(sinh(η̄−Tθ) on S, sinh(η̄+Tθ) off S),
- with η̄ = η + iπω(l − N/2).
- So the antisoliton sees a shifted parameter, and η̄ = η only for l = N/2.

**Unitarity.** It is scalar: K(θ)K(−θ) = U·1 with U = sinh(η+Tθ) sinh(η−Tθ).

**Crossing.**
- The boundary state pairs a soliton with its antisoliton, so the crossed channel uses Ř_{VV̄}. Write k_A = Σ (K_V̄)_ii(iπ/2−θ) v_i⊗v̄_i and k_B = Σ (K_V)_ii(iπ/2−θ) v̄_i⊗v_i.
- Then k_B(θ) = Φ₁ Ř_{VV̄}(2θ) k_A(−θ) and k_A(θ) = Φ₂ Ř_{V̄V}(2θ) k_B(−θ), with residual 1e-15 for n = 1..3 and all l.
- Closed form (n ≤ 4, two values of ω): Φ_{1,2} = −(−q)^{±(1−l)} · sin π(2μ−T)/sin π(2μ+T).
  - It does not depend on η.
  - Rescaling K_V̄ by λ = (−q)^{1−l} makes both equal to −s(μ), where s = sin π(2μ−T)/sin π(2μ+T).
- The transmission entries of Ř_{V̄V} all equal 1. So the soliton–antisoliton bulk factor is F(2μ) = S^T(T−2μ), the crossed transmission amplitude.

**Scalar equations** (a multiplies K_V, b multiplies λK_V̄):
- a(μ)a(−μ) = 1/U and b(μ)b(−μ) = 1/Ū, with Ū = λ² sinh(η̄ + iπμ) sinh(η̄ − iπμ).
- a(T/2−μ) = Y b(T/2+μ) and b(T/2−μ) = Y a(T/2+μ), with Y = −S^T(T−2μ)/s(μ).

### 3. Minimal scalar factor (closed form)

**Notation.**
- Ê_P(z) = Γ(z) Π_{k≥1} Γ(z+kP)/Γ(1−z+kP) · (kP)^{1−2z}. Then Ê_P(z+P)/Ê_P(z) ∝ sin πz and Ê_P(z)/Ê_P(z−1) = (P/π) sin(π(z−1)/P).
- Λ_P(x;c) = Ê_P(c+x)/Ê_P(c−x).
- ε = η/(iπ) and ε̄ = ε + ωl − T.

**Identity derived and checked:** F₁₁(x−2T) = F₁₁(x) · sin πx sin π(x−2T)/(sin π(x−2T+ω) sin π(x−ω)).

**Result:**
- a(μ) = i R(μ) (ω/2)_I γ(μ)/sin π(ε+μ)
- b(μ) = (i/λ) R(μ) (ω/2)_I γ̄(μ)/sin π(ε̄+μ)
- R(μ) = Λ_{4T}(2μ; 2T+ω)/Λ_{4T}(2μ; 0)
- γ = λ^{−2μ/T} e^{kμ} Λ_{2T}(μ; ε+2T)/Λ_{2T}(μ; ε̄+T)
- γ̄ = λ^{2μ/T} e^{−kμ} Λ_{2T}(μ; ε̄+2T)/Λ_{2T}(μ; ε+T)
- k is fixed numerically on the minimal branch; it is 0 for l = 1.

**Checks.**
- All four equations hold to about 1e-7 for (n,l) = (1,1), (2,0), (2,1), (3,0), (3,1), (3,2).
- **n = 1:** a/[GZ R₀ σ(ξ,u)/cos ξ], with ξ = π/2 − πε, equals −i(T/2)_I² for every η. So GZ Dirichlet is our minimal solution times the CDD factor (−T/2)_I². That factor is unitary and crossing-symmetric.
- **CDD for n ≥ 2 is not fixed.** Two natural candidates:
  - σ_A = (−T/2)_I²;
  - σ_B = −[(T/2)_I² (T/2−ω/2)_I (T/2+ω/2)_I]⁻¹, minus the inverse of problem 3's CDD, which satisfies DG's bootstrap condition.
  - Both reduce to GZ at n = 1. I did not run the SP soliton bootstrap V^{⊗n} → V̄ that would decide between them.

### 4. Breather amplitudes

**Structure.**
- The singlet in V(θ_a)⊗V̄(θ_b) is mapped exactly (residual 1e-15) onto the singlet in **V̄(−θ_b)⊗V(−θ_a)**. That is B_1 → B̄_1.
- Projection coefficient: c = U(μ_a) c₀(μ), with c₀ independent of η and l. For B̄ → B the analogue involves Ū and c̄₀, and c̄₀ depends on l.

**Closed form** (singlets normalised on v₀⊗v̄₀; overall sign fixed by n = 1):
- K_{1→n}(μ) = −K₀(μ) E_η(θ) × CDD,
- K₀ = (T/2)_I (T/2−1/2)_I (δ)_I (1/2−δ)_I (1/2+δ)_I (T−1/2−δ)_I, with δ = (T−ω)/2. Fitted for n = 1, 2, 3 (residual ≤ 5e-6).
- E_η = (sinh θ − i sin(πA/T))/(sinh θ + i sin(πĀ/T)), with A = ε + (T−1)/2 and Ā = A + ωl − T. Derived analytically from the Ê shift rules and confirmed by η-ratio fits for n ≤ 3, l = 0, 1, 2.
  - For l = 0 and l = N, E_η = 1, so η is not a genuine parameter there.
- K_{n→1} is the same with A ↔ Ā. Then K_{1→n}(θ) K_{n→1}(−θ) = 1.

**n = 1 check.** With (−T/2)_I² the amplitude equals Ghoshal's Dirichlet K_{B₁}, up to the singlet sign convention, for all η. This was checked directly to 1e-6. Continued at η = iπ/2 it gives exactly K_D = (2−B/2)(1+B/2)/(1).

**Continuation** ω → −2/B, with (a)_I → (−aB):
- K₀ → (h/2)(h/2−1)(h/2+B/2)(h/2+1+B/2)(h/2−1−B/2)/(h/2−1+B/2).
- σ_A contributes (B/2)²/(h+B/2)².
- σ_B contributes (B/2)²(B/2+1)(B/2−1)/[(h+B/2)²(h−1+B/2)(h+1+B/2)].
- Both CDD contributions tend to 1 as B → 0.

**Axioms.** Unitarity holds by construction. Crossing and the bulk boundary bootstrap of the continued amplitudes were **not checked**.

### 5. Which real-coupling boundary — verdict

- **B → 0 limit.** At η = iπ/2, the analogue of φ₀ = 0, where K_V ∝ diag(e^{a_jωθ}):
  K_{1→n} → −(s+i)(s+i cos(π/h)) / [(s−i cos(π/h))(s−i cos(2πl/h))], with s = sinh θ.
  - This is −1 for h = 2.
  - For h ≥ 3 it depends on θ for every l, and no η makes it constant: E_η supplies one zero–pole pair, but K₀ has two to cancel.
  - So it is **neither −1 (Dirichlet) nor +1 (Neumann)**. In addition, particles reflect a → h−a, which neither Dirichlet nor Neumann on all components can do. (Classically, Neumann on τ-even components with Dirichlet on τ-odd ones gives a → h−a with amplitude 1.)
- **Which statement generalises.**
  - The n = 1 statement (Dirichlet, K_D) is a sine-Gordon accident: there V̄ ≅ V, B̄_1 = B_1, and real Dirichlet is integrable.
  - "σ = † with K_+ = 1 at real coupling is Neumann" is true classically, and quantum mechanically that boundary is the soliton-conjugating s = +1 matrix of problem 3, not an SP matrix.
  - Algebraically, at real λ, σ = † has the r-matrix structure of σ = T (twisted equation). The braid (SP) equation corresponds to τ ≠ id (Appel–Vlaar), as the book's ch. 34 already says.
- **Comparison with problem 3: no conflict; the two results support each other.**
  - Neumann reflects a → a and is Z_N-symmetric, as the s = +1 soliton-conjugating matrix is.
  - SP matrices reflect B_1 → B̄_1 and break Z_N.
  - The brief's hoped-for "second soliton route to Neumann" does not exist.
  - Problem 3's warning about `k_closed` is irrelevant here; I did not use it.
  - Problem 3's method carried over, with two changes: SP unitarity is plain K(θ)K(−θ) ∝ 1, and the antisoliton matrix is not JKJ but K_V̄ with the shifted η̄.

### Proposed book text (short)

- **`sec-preserving-conjugating`.** Replace "(unchecked in book conventions)" with the K_V formula above, the gradation map to DVGR, K_V̄ with η̄ = η + iπω(l − N/2), and the crossing function −s(μ). Add as **[Theorem, by hand]**: "For n ≥ 2 no invertible diagonal solution of the braid reflection equation is Z_N-covariant. Every SP matrix singles out two nodes, 0 and l."
- **Ch. 35, new paragraph "Soliton-preserving boundaries" [Numerical, n ≤ 3]:**
  > The scalar factor of the SP matrices is a Gamma product (formula above); for n = 1 with CDD factor (−T/2)_I² it is Ghoshal–Zamolodchikov Dirichlet. The boundary maps B_1 to B̄_1, with amplitude −K₀E_η. Continued, its B → 0 limit depends on θ, so it describes neither Dirichlet nor Neumann for n ≥ 2. The CDD factor for n ≥ 2 and the identification with a classical boundary are open.
- **Open problem 11:** mark as partially solved. What remains: the SP soliton bootstrap and CDD, crossing and bootstrap of the particle amplitudes, a classical boundary for the DVGR family, and the RV 2018 match.
- **`sec-soliton-preserving` correction:** coweights not in 2Λ_W^∨ are allowed with K_+ = e^{iπμ·H} (by hand; to be checked).

### Scripts

All in `/tmp/claude-1000/-home-gustav-Git-affine-toda/0e38a99d-702a-40fb-99b7-74dbe0413cf6/scratchpad/p11/`. This is a temp directory, so copy it soon. It also holds unchanged copies of `textbook_code/*.py` and of problem 3's `nk_common.py`, `gz.py` and `rhonp.py`.

| Script | What it does |
|---|---|
| `sp_core.py` | K_V, K_V̄, the crossing function Φ |
| `e1.py`–`e4.py` | Ř closed form; solving the block reflection equation |
| `t1.py` | braid reflection equation; mixed equation solved for K_V̄ |
| `t2.py`–`t6.py` | antisoliton block constants; fit of η̄ |
| `t7.py` | crossing-unitarity, matrix part |
| `t8.py`, `t9.py` | closed form of Φ |
| `scal.py` | Ê, Λ, F₁₁, S^T, R, γ, γ̄, `amps` |
| `t10.py`, `hfun.py` | GZ σ relations |
| `t11.py`–`t13.py`, `t15.py` | shift relations, unitarity and crossing of the scalar factors |
| `t14.py` | rational fits; crossing residual of R |
| `t16.py`, `t17.py` | n = 1 against GZ and Ghoshal |
| `bproj.py` | breather singlet projection |
| `breather.py` | breather amplitude, block fitting |
| `t18.py`, `t19.py` | η-ratio fits and fit of K₀ |

Not done, ranked: SP soliton bootstrap to fix the CDD (adapt problem 3's `sol_boot.py`); crossing and bootstrap of the continued amplitudes; n = 4 breather; RV 2018 match; numerical coideal check.

## Integration (branch half-line-open-problems-integration)

| Claims | Verdict | How it was checked | Book change |
|---|---|---|---|
| 1, 2 (SP matrices in book conventions; DVGR map) | **Verified, integrated** | The coordinator tested the formula in the book's own braid reflection equation (`boundary_common.re_preserving`, `rcheck`, physical $q$): $\le10^{-15}$ for $n=2,3$, $l=0..N$; flipped exponents fail. The gradation map to DVGR was checked numerically up to a scalar. New `textbook_code/soliton_preserving.py` (3 checks, PASS). | ch. 33 `sec-preserving-conjugating`: formula replaces "(unchecked in book conventions)". |
| 3 ($\mathbb Z_N$ obstruction) | **Verified**, including the loopholes the coordinator's Opus prompt asked about | A diagonal gauge commutes with diagonal $K$; a scalar or twisted shift still forces equal neighbouring ratios $k_j/k_{j+1}$, which the family violates (also checked numerically for $n=2,3,4$). | Stated in ch. 33 (no separate proposition). |
| 6 ($K_{\bar V}$ with $\bar\eta=\eta+i\pi\omega(l-N/2)$; scalar unitarity) | **Verified, integrated** | The mixed braid equation was derived from factorization in book conventions and solved to $2\times10^{-15}$; the control $\bar\eta=\eta$ fails where it differs. Uniqueness was not checked, and the text does not claim it. | ch. 33. |
| 9 ($B_1\to\bar B_1$; particles $a\to h-a$) | **Verified by structural argument** (SP $K$ keeps $V\to V$, $\bar V\to\bar V$, so the constituent order of the singlet is reversed relative to the conjugating case), consistent with the worker's and Gemini's numerical projection | | ch. 33 bullet; open problem 11. |
| 10 (neither Dirichlet nor Neumann for $n\ge2$) | **Integrated in the structural form** (from $a\to h-a$). The $\theta$-dependent $B\to0$ formula, which comes from a fitted $K_0$, is **not** integrated. | | ch. 33; open problem 11. |
| 7, 8 (crossing function; Gamma-product scalar factors; CDD candidates; breather amplitude $-K_0E_\eta$) | **Not integrated**: no independent recomputation yet (Gemini reran the worker's scripts) | | Open problem 11 lists them as open. |
| 4 (book's coweight exclusion wrong; $K_+=e^{i\pi\mu\cdot H}$) | **Pending**: by hand only; Gemini restated without a check | | — |
| 5 (real Dirichlet not integrable for $n\ge2$) | **Pending**: Gemini's re-derivation drops the constants $e^{\beta\alpha_j\cdot\phi_0/2}$ and is incomplete; prior art (BCDR) not checked | | — |
