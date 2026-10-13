<!-- Opus referee follow-up, 2026-10-12. Model: Claude Opus 5.5, effort max. Prompt: briefs/halfline_open_03_opus_prompt.md. Combined session with 12a (see halfline_open_12a_referee_opus.md). Saved verbatim by the coordinator. Scripts: /tmp/claude-1000/-home-gustav-Git-affine-toda/0e38a99d-702a-40fb-99b7-74dbe0413cf6/scratchpad/ref_03_12a/ (s03/). -->

## 03: Neumann amplitudes from the s=+1 soliton matrix

### Item 1, claim 1 against the book: CONFIRMED, with a refinement. The book is incomplete in ch. 28 and wrong as written in ch. 33.

**What the book says.**
- **ch. 28, `eq-boundary-unitarity`:** "K(θ)K(−θ)=1 … and the obvious modification for conjugating boundaries". The modification is never specified, so this is incomplete.
- **ch. 33, `sec-k-procedure` step 4:** "Boundary unitarity, K(θ)K(−θ)=1". This is wrong for soliton-conjugating non-diagonal K, and it is even ill-typed, since K(−θ) maps V(−θ)→V̄(θ).
- **No damage to results.** Nothing in ch. 28–35 applies it to the non-diagonal matrix. It is used only for sine-Gordon and K=1 (ε̂=0), where K_V̄=K.
- **The ch. 28 crossing-unitarity statement for twisted aₙ** (k=Σ K_ba ν̄_a⊗ν̄_b with Ř_V̄V̄) is consistent with the worker.
- **ch. 34–35** contain no use of K(θ)K(−θ) for the non-diagonal K.

**Analytic derivation in book conventions** (`boundary_common.an_vector`; any q):
- **(a) Charge conjugation.** J π_V̄(y)(x) J = π_V(y)(τx) for every eⱼ, fⱼ, kⱼ, where τ: j→n+1−j fixes node 0. Checked to be exact (residual 0).
- **(b) The antisoliton matrix.** If ε̂ is τ-invariant, conjugating K π_V(y)(bⱼ) = π_V̄(1/y)(bⱼ) K by J shows that JK(y)J intertwines V̄(y)→V(1/y). By uniqueness K_V̄ ∝ JKJ.
- **(c) Bulk C-invariance.** τ is a Hopf automorphism of the book's coproduct and J maps highest weight to highest weight. Hence Ř_V̄V̄ = (J⊗J)Ř_VV(J⊗J) and Ř_V̄V = (J⊗J)Ř_VV̄(J⊗J), numerically to 3e-15.
- **(d) Same scalar.** With a C-invariant boundary the antisoliton amplitude is A·JKJ with the same A. The sign of q plays no role in this step.
- **Numerical check.** n=2..5, ω=0.37, 1.29, 2.6, s=±1: residuals 1.9e-15. K(−θ)ᵀK(θ) is scalar to 1.8e-15; K(θ)K(−θ) has off-diagonal entries ≥0.38 of the diagonal.

**Is JKJ=Kᵀ an identity?** Yes, for (`eq-an-k-matrix`) as printed. The matrix is Toeplitz (Kᵢⱼ depends only on j−i), hence persymmetric. This holds for all θ, ω, s and either branch. But (`eq-an-k-matrix`) is itself the τ-symmetric gauge representative. For the other sign patterns of ε̂ in the same class s (gauge transforms K′=DKD, d_i=±1, D²=1):
- K′ᵀ is still the antisoliton intertwiner (1.9e-15);
- JK′J is not (residual ≥2.0).

In a general torus gauge (not parity-symmetric), K′_V̄ = D²K′ᵀD², so neither statement holds. So "antisoliton = JKJ" needs τ-invariant ε̂. "Antisoliton = Kᵀ" holds throughout the parity-symmetric normalization.

**Proposed text.**
- **ch. 28, after (`eq-boundary-unitarity`):** replace "and the obvious modification for conjugating boundaries" with: "For a soliton-conjugating boundary the reflected particle is in the conjugate multiplet, and the second reflection is that of the antiparticle: $K_{\bar\mu}(-\theta)K_\mu(\theta)=1$, with $K_{\bar\mu}(-\theta):V_{\bar\mu}(-\theta)\to V_\mu(\theta)$. For the $a_n^{(1)}$ matrices of @sec-an-nondiagonal-k, $K_{\bar V}=K^{\rm T}$, not $K$."
- **ch. 33 step 4:** "Boundary unitarity (-@eq-boundary-unitarity), $K_{\bar\mu}(-\theta)K_\mu(\theta)=1$, and …".
- **sec-an-nondiagonal-k, refining the worker's "Unitarity" paragraph:** "Let $J$ reverse the basis. Then $J\pi_{\bar V(y)}(x)J=\pi_{V(y)}(\tau x)$, where $\tau:j\mapsto n+1-j$ fixes node 0. So $J$ is charge conjugation, and $\check R_{\bar V\bar V}=(J\otimes J)\check R_{VV}(J\otimes J)$. Since the $\hat\epsilon$ of (-@eq-an-k-matrix) is $\tau$-invariant, $JKJ$ is the intertwiner $\bar V(\theta)\to V(-\theta)$. As (-@eq-an-k-matrix) is a Toeplitz matrix, $JKJ=K^{\rm T}$. The other sign patterns are $DKD$ with $D^2=1$; for them $K^{\rm T}$ is still the antisoliton matrix but $JKJ$ is not. Boundary unitarity is $A(\theta)A(-\theta)K(-\theta)^{\rm T}K(\theta)=1$; $K(\theta)K(-\theta)$ is not a multiple of $\mathbf 1$ for $n\ge2$." Then add the worker's U formula and its proof.

### Item 2, Kim's one loop: CONFIRMED

- **Kim's equations.** In hep-th/9506031 (TeX `9506031.tex`), (33) is lines 909–912: K₁=1+(iβ²/16)(sinhθ/(coshθ+1/√2) − sinhθ/(coshθ−1)). (39) is lines 958–962. My equation count matches the book's citation. B(β)=β²/(2π(1+β²/4π)) is at line 176, and his blocks (x) at line 480 are the book's.
- **With the factor.** At h=4 the O(β²) terms of K^N₁ and K^N₂ equal (33) and (39) identically (residual 2e-32 and 5e-32 at three θ, 50 digits).
- **What the factor continues to.** σ(μ_a)σ(μ_b), with σ=(T/2−ω/2)_I(T/2+ω/2)_I at n=3 and ω→−2/B, equals DG's k₁ (c=h/2−1=1) to 6e-51.
- **k₁ at O(β²).** k₁ = 1 + (iβ²/8)(sinhθ/(coshθ+1/√2) − sinhθ/(coshθ−1/√2)) + O(β⁴), which is nonzero (verified to 7e-32).
- **Without the factor.** K₁ = K^N₁/k₁ and K₂ = K^N₂/[k₁(θ+iπ/4)k₁(θ−iπ/4)]. Their O(β²) coefficients deviate from (33) by up to 0.123 and from (39) by up to 0.66.
- **Conclusion.** Kim's one loop selects the factor at n=3.

Proposed addition to the proviso: "For $n=3$ without this factor the amplitudes become $K^N_1/k_1$ and $K^N_2/[k_1(\theta+\frac{i\pi}4)k_1(\theta-\frac{i\pi}4)]$, which contradict Kim's (33) and (39) at $O(\beta^2)$."

### Item 3, the breather identity at high precision: CONFIRMED (the 37–40 digits are mine)

**What was evaluated.** The worker's closed form A⁽⁺⁾ = a₀⁽⁺⁾ a₁ (T/2−ω/2)_I(T/2+ω/2)_I and the book's F₁₁. Every infinite Γ-product was computed two independent ways at 40 digits:
- **(I)** a Malmsten integral for the tail, with E₃-regularised numerator;
- **(II)** a Stirling series summed with Hurwitz ζ (J=60, 45 terms).

None of the worker's truncations were used.

**|LHS/RHS − 1| for the breather identity:**

| n | ω | μ | method (I) | method (II) |
|---|---|---|---|---|
| 3 | 0.61 | 0.17+0.23i | 9.3e-41 | 1.6e-37 |
| 4 | 0.47 | −0.21+0.13i | 2.0e-40 | 4.9e-37 |
| 3 | 1.37 | 0.29−0.11i | 2.7e-40 | 1.3e-36 |

- **Point values.** LHS = RHS = 0.0191110522994533768285 − 2.3458278182419274313658 i at the first point, and −0.1353858512005747452013 − 0.2310583560841050676635 i at the second.
- **Functional equations of the closed form.** A(μ)A(−μ)U = 1 and A(T/2−μ)/A(T/2+μ)·Φ/F₁₁(2μ) = 1 both hold to ≤4e-40 (I). F₁₁ unitarity holds to ≤2e-43.
- **Matrix side at the same points.** With my own K and the book's `rcheck`, the singlet maps to c₊·singlet with |c/c₊−1| ≤ 1.6e-15.
- **Pitfall.** The in and out singlets must be normalised identically; they are the same vector (to 1e-15), because the singlet depends only on y₁/y₂. With arbitrary SVD normalisation I first got c = −c₊ for n=3. The worker normalised correctly.
- **Proposed book wording:** "…equals the expression above to $10^{-36}$ at generic points for $n=3,4$ (infinite products by integral representation), and to $10^{-7}$ for $n\le5$; the matrix part $c_+$ to $10^{-15}$."

### Item 4, prior art: CONFIRMED that none is found

**What I searched.**
- INSPIRE citations of hep-th/9806003 (29), hep-th/9911178 (34), hep-th/9904002 (28) and hep-th/0110218 (6), screened by abstract.
- INSPIRE full text: "scalar factor", "crossing" and "Neumann" within those citation sets, plus "soliton non-preserving" with "crossing-unitarity" or "scalar factor".
- TeX sources grepped: Delius–MacKay hep-th/0112023 (line 908, "undetermined overall factor k(θ)"), Regelskis–Vlaar 1803.06491 ("up to an overall scalar factor"), Baseilhac–Belliard 0906.1215, Belliard–Fomin 1106.1317 ("scalar" means c-number), and Doikou hep-th/0006197 (rational SU(3) chain, no boundary S-matrix).
- Doikou–Nepomechie hep-th/9906069 treats the critical-chain boundary S-matrices, not the non-diagonal SNP K.
- Also checked: Doikou 0803.0943, Avan–Doikou 0809.2734 (classical), Fateev–Onofri hep-th/0203131 and Castro-Alvaredo–Fring hep-th/0308004 (particle amplitudes only). Kuniba–Okado–Yoneyama, Lima-Santos, Yang–Zhang and Nepomechie hep-th/0204181 give matrices without scalar factors.
- Bowcock–Umezawa returned no INSPIRE hits; Bowcock's half-line papers are classical.

**Result.** No published scalar factor of the non-diagonal aₙ⁽¹⁾ K, and no soliton-side Neumann derivation for n≥3.

## Scripts (03 part)

All in `/tmp/claude-1000/-home-gustav-Git-affine-toda/0e38a99d-702a-40fb-99b7-74dbe0413cf6/scratchpad/ref_03_12a/`. Downloaded sources are in `dl/` and the copies of the book and worker code in `code03/`; run the 03 scripts with `PYTHONPATH=code03:s03`.
- `s03/item1_antisoliton.py`: charge conjugation J, the antisoliton intertwiner, persymmetry, C-invariance, non-τ-symmetric patterns.
- `s03/item2_kim.py`: Kim (33)/(39) against K^N, k₁ continuation, the version without the factor.
- `s03/hp_breather.py` (output `s03/hp_breather.out`): 40-digit breather identity and functional equations, integral and Hurwitz methods.
- `s03/item3_projection.py`: c₊ from the matrix projection at the same points.
- `s03/parse_inspire.py`, `s03/grep_abs.py`: INSPIRE screening (data in `dl/inspire/`).
- `count_eq.py`: Kim's equation numbering.
