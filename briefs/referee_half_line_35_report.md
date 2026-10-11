## Referee report: chapter 35 (`part7-half-line/35-particle-reflection.qmd`)

I read everything from the worktree and edited nothing in the repository. Scratch scripts and outputs are in `/tmp/claude-1000/-home-gustav-Git-affine-toda--claude-worktrees-affine-toda-boundary-dae69e/e8814860-aa66-4a52-8bab-4efe4e4077d9/scratchpad/referee_35/`. I re-ran `textbook_code/particle_reflection.py` and all 43 checks pass.

**Headline findings:**
- **(A) The chapter's "correction" of DG (3.11) is wrong.** The real error is a misquoted prefactor in ch. 33's `eq-an-diagonal-scalar`. It affects ch. 33 and 35 and two scripts.
- **(B) The matrix-level breather bootstrap now works for n = 2, 3.** It confirms `prp-breather-amplitude`, but only with Delius–Gandenberger's own prefactor.
- **(C) The "rotor" reading is now backed by a calculation.** The classical collective-coordinate inertia of the modulus χ equals exactly the value the b_{n,−n} tower requires.
- **(D) Delius–Gandenberger's a_2 spectrum (and so the table) leaves one listed pole unaccounted for.**
- **(E) Later literature is missing:**
  - Ahn–Kim–Rim give a non-perturbative TBA check of the uniform and Neumann amplitudes.
  - Castro-Alvaredo–Fring and Baseilhac–Belliard bear on open problems 2, 5 and 9.

---

### 1. Corrections to Delius–Gandenberger
I could not get the journal version. ScienceDirect returns 403 to scripted access. The arXiv v2 (2 Jun 1999, comment "minor corrections") is the latest version, and INSPIRE confirms NPB 554 (1999) 325, DOI 10.1016/S0550-3213(99)00304-1. All verdicts below are against v2, read from the rendered PDF.

**(6.26): correct.**
- As printed (p. 34), the left side is K₂^{(2)} = K₃^{(1)}, and the second line has only three blocks, (5/2−B/2)(5/2+B/2)(−1/2−B/2).
- The fourth block must be (−9/2+B/2). Crossing-unitarity needs each extra factor built from pairs (x)(h−x) mod 2h. With h = 5, (5/2−B/2)(5/2+B/2) is such a pair; (−1/2−B/2) needs the partner 11/2+B/2 ≡ −9/2+B/2.
- K₃^{(1)} is a label misprint: charge conjugation gives K_a^{(b)} = K_{h−a}^{(b)}, so it should read K₃^{(2)}.
- Unchecked: whether the journal fixed these.

**(6.11) line (d): correct.**
- Printed (p. 26) as (−1−(2n₃+1)B/2)(−3−(2n₁+1)B/2).
- The pairing x(n₁,n₃), h−x(n₃,n₁) with h = 4 gives the partner 5+(2n₁+1)B/2 ≡ −3+(2n₁+1)B/2.
- At n₁ = n₃ = 0 the corrected version gives K₁(1−B/2)²(3+B/2)², which is eq-solitonic-factor at b = 2. The printed version does not.

**(3.11): the chapter's correction is WRONG.** DG (3.11) is correct as printed: the breather amplitude is exactly A·F₁₁·A.
- `r1_scalar_product.py` and `r1b_prefactor.py`, n = 1–4, two values of λ, two values of μ: A(μ−hλ/4+½) F₁₁(2μ) A(μ+hλ/4−½) divided by the closed form (3.11) is 1 to 1e-10. That is the accuracy of the infinite-product extrapolation.
- This holds with DG's (3.9) exactly as printed (PDF p. 11): prefactor (−(n+1)λ/4)_I and **no overall minus sign**.
- The book's `eq-an-diagonal-scalar` (ch. 33 line 91) prints −(h/4·ω)_I instead. Then A_book = A_DG·f(θ) with f = −(h/4·ω)_I² = −(1−i sinh θ)/(1+i sinh θ), verified to 1e-26.
- With the book's A, the product picks up an extra factor f(θ_a)f(θ_b), verified exactly. That is the "further scalar function of μ" the chapter attributes to the matrix part. It does not come from the matrix part (item 2).
- `r1c_f11_sg.py`: F₁₁ at h = 2 equals the sine-Gordon S₀ (eq-sg-s0), ratio 1 to 1e-10. So at n = 1 the DG prefactor reproduces Ghoshal's K_{B₁}(η = π/2, ϑ = 0), and the book's prefactor contradicts its own `prp-bsg-breather-k`.

Corrected text:
- **ch. 33 line 91:** `A(\mu)=\Big(-\frac h4\omega\Big)_I\prod_{k=1}^\infty …` (drop the leading minus, move the sign inside the block).
- **ch. 33 line 105, replace the last sentence with:** "The equations (-@eq-an-diagonal-scalar-eqs) do not fix the prefactor: $(-\frac h4\omega)_I$ and $-(\frac h4\omega)_I$ differ by $-(1-i\sinh\theta)/(1+i\sinh\theta)$, which is unitary and crossing-symmetric and drops out of the semiclassical limit. Delius and Gandenberger's $(-\frac h4\omega)_I$ is the choice for which the breather bootstrap gives DG (3.11), and for $n=1$ Ghoshal's $K_{B_1}$ at $\hat\epsilon=0$."
- **ch. 33 line 260:** same content as the line-105 replacement.
- **Code:** in `reflection_matrices.py`, `A1` should be `blk(-h*lam/4, …)` without the leading minus; its checks still pass. The docstring of `particle_reflection.py`, lines 15–18, must be fixed.
- **ch. 35 line 286:** replace the DG (3.11) bullet with: "- **DG (3.11).** Correct as printed, with DG's prefactor in (3.9): the matrix part of $S_{\bar V\bar V}$ acts trivially on the reflected singlet (@sec-breather-reflection)."

### 2. Breather bootstrap for n ≥ 2: done, it works
`r2_matrix_bootstrap.py` uses `boundary_common.an_vector` and `rcheck`, i.e. the book's Ř in the principal gradation, normalized to 1 on v₁⊗v₁. Setup: n = 2, 3; ω = 0.37, 1.83; two values of θ; both orderings.
- K = 1 intertwines the ε = 0 coideal (V→V̄ and V̄→V).
- The singlet exists in V(θ−iu/2)⊗V̄(θ+iu/2), with u = π(1−2/(hω)).
- (1⊗K)Ř(2θ)(1⊗K)w = 1·w_out with deviation below 1e-15. The incoming and outgoing singlets are the same vector (|w_in−w_out| ~ 1e-16), because in the principal gradation the singlet depends only on y_a/y_b.

Why: after the first reflection the state is Σcᵢ vᵢ⊗vᵢ. Ř acts on every vᵢ⊗vᵢ by the identical-soliton eigenvalue, which is 1. So K_{B₁} = A·F₁₁·A exactly. `prp-breather-amplitude` [Established] is confirmed, but only after the ch. 33 prefactor fix.

**ch. 35 line 37, replace the paragraph with:** "For $n=1$ the identity with Ghoshal's formula is checked in `textbook_code/particle_reflection.py`, and the matrix-level projection is in `textbook_code/boundary_sine_gordon.py`. For $n\ge2$ the projection is simple. After the first reflection the singlet becomes $\sum_ic_i\,v_i\otimes v_i$, a sum of identical-soliton states, on which $\check R$ acts by its identical-soliton eigenvalue, normalized to 1. So $K_{B_1}(\mu)=A(\mu-\frac{h\omega}4+\frac12)F_{11}(2\mu)A(\mu+\frac{h\omega}4-\frac12)$, and with (-@eq-an-diagonal-scalar) this equals (-@eq-breather-k-imaginary). (We checked the matrix projection for $n=2,3$ and the scalar identity for $n=1,\dots,4$.) The prefactor of (-@eq-an-diagonal-scalar) matters: with $-(\frac h4\omega)_I$ instead, the product acquires $f(\theta_a)f(\theta_b)$ with $f=-(1-i\sinh\theta)/(1+i\sinh\theta)$, and for $n=1$ it would contradict Ghoshal."

### 3. Neumann
- **Label: correct.** `cnj-neumann-amplitude` [Conjecture] is right.
- **Text: correct.** Lines 125 and 287 state clearly that Neumann is the s = +1 non-diagonal matrix and that the soliton-route derivation exists only for a_2 (G98 (6.4)).
- **Kim citation: correct.**
  - Kim's (33) and (39) are in his appendix "a₃ affine Toda theory" (rendered PDF pp. 16–17). The script's formulas match.
  - Kim's identity (41) is quoted exactly at line 120. I also derived it.
- **"for a₂ with Kim's result as quoted by Gandenberger" (line 66): slightly inaccurate.** Kim's paper gives no explicit a₂ one-loop result, and G98 quotes only its own expansion while asserting a match with Kim (footnote e). Suggested: "…and for $a_2^{(1)}$ Gandenberger reports agreement with Kim's one-loop formulas [@gandenberger1999a, sec. 6]."
- **Kim's exact conjecture, checked** (`r3_kim.py`). Kim's (17), K₁ = ∏[k/2], agrees with DG (3.20) at O(B) for n = 2–6 (and K₂ for a₃) and differs at O(B²), with coefficient about 0.35–0.48. This confirms line 123.

### 4. The n = 1 discrepancy
**Verdict: both results are right in their own terms; the extension is legitimate, given its assumptions.**
- `cnj-shg-e` is stated and derived only for ½ < a < 1 (ε < 0). Going to C₀ = C₁ = 1 (a = 0) is an analytic continuation in a.
- That continuation is exactly the continued Zamolodchikov/BPT relation (`cnj-bsg-uv-ir`), as Corrigan–Taormina and BPTT observed.
- Through the same relation, E = B/2 (the ε̂ = 0 point) corresponds to C₀ = C₁ = cos(πB/(2(2−B))) = 1 − O(B²). Check: 2a(1−B/2) = B/2 gives a = B/(2(2−B)). This equals the line-101 sine-Gordon statement continued with λ → −2/B.
- So if `cnj-bsg-uv-ir` holds, the sinh-Gordon uniform boundary C = 1 has E = 0, and DG (3.16) at n = 1 describes a renormalized C slightly below 1. This is a boundary-parameter renormalization, possible only because the C_j are continuous for n = 1.
- Chenaghlou hep-th/0108077 is second order in C₀−C₁ at one loop, not two-loop, so problem 12 is still open.

**Line 101, replace "for sinh-Gordon it is undecided whether … exactly at $\hat\epsilon=0$" with:** "If the continuation of Zamolodchikov's relation (@cnj-bsg-uv-ir) holds, as Corrigan and Taormina's extension of @cnj-shg-e assumes, the point $\hat\epsilon=0$ is $C_0=C_1=\cos\frac{\pi B}{2(2-B)}$, not $C\equiv1$: the uniform sinh-Gordon boundary then has $E=0$, and (-@eq-uniform-amplitude) at $n=1$ describes a boundary whose parameter is renormalized at $O(B^2)$."

Make the same change in problem 12 (line 240).

### 5. Unchecked quoted results
- **Bowcock (4.21) against DG (5.7): agrees; checked by hand.**
  - Bowcock's (4.14) with one soliton pair (N = 1 in his 2N count), multiplied by (4.21) and (4.22), gives (h/2+s−p)(h/2−s+p)/[(h/2+s+p)(h/2−s−p)] × (−1/((h−p)(p))).
  - Using (−x) = 1/(x), these are exactly DG's four blocks times the uniform limit. His blocks (4.20) are in the usual rapidity φ = θ.
  - Line 65: replace the last sentence with "For the solitonic boundaries the limit (DG (5.7)) is Bowcock's linearized factor about the continued vacuum [@bowcock1998, (4.14), (4.21), (4.22)]."
  - Line 154: delete "; we have not checked Bowcock's formula".
- **Perkins–Bowcock: agrees, with two caveats** (`r5b_perkins_bowcock.py`).
  - Their (5.3) and (5.5) each lack a factor i. The expansion of DG (3.16) at n = 2 is i×(5.5) exactly. With i restored in g, f+αg (α = −1/(16√3)) equals the expansion of K₁ = (3−B/2)/((2)(1−B/2)) to 1e-18.
  - α is a boundary counterterm that PB fix by demanding the bootstrap, so the check tests one function up to one fitted constant.
  - Line 66, replace with: "The uniform amplitude for $a_2^{(1)}$ agrees at $O(\beta^2)$ with Perkins and Bowcock [@perkins_bowcock1999, (5.5)] (checked; their (5.3) and (5.5) each lack a factor $i$), after a finite boundary counterterm whose coefficient they fix by the bootstrap."
- **a₄ excited amplitudes (6.20), (6.21), (6.28), (6.29): not checked here (time).** They remain checked only at the vacuum, as the text says. Unclear; settled by extending `K1nm`-style checks to h = 5.
- **Coleman–Thun explanations for a₂: checked against DG's figures; they match.**
  - `r6_poles.py` classifies every strip pole of K₁ and K₂ on each table state at B = 9/20 by DG §6.2 (figs. 10a/b/c, 11, 12a/b/c).
  - Orders are consistent with the DTW count P−2L, minus one per reflection at a zero. That includes the double pole on b₀,₁ (10c + 11), which the text states.
  - Exception: see item 6.

### 6. Energies, table, existence ranges, rotor
- **Closed form: correct.** I derived F(n+1)−F(n) = cos(π/6−(2n+1)δ) by hand; e(b₋₁,₁) = m₁ sin δ; the degeneracy at B = ½ checks out. Every table energy matches.
- **Table and ranges: inconsistency inherited from DG.**
  - K₁ on b₋₂,₂ has a strip pole at p = ½+3B/2 = 47/40 (listed in the table), and K₂ on b₂,₋₂ has the same pole.
  - DG's own explanation (text after fig. 11: "(½+(2n−1)B/2) corresponds to 10a in which b₋ₙ₊₁,ₙ is created") makes it create b₋₁,₂ (resp. b₂,₋₁).
  - DG's range |n|,|m| < 1/(2B)+½ excludes that state for B > 1/3. In general the pole is in the strip for B < 2/(2n−1), but the state lies in DG's range only for B < 1/(2n−1).
  - The same gap affects DG's Fig. 1 (B ≈ ½), which does not show b₋₁,₂.
  - The script's check is one-directional: it checks creation poles between listed states, not that every listed pole is explained.
  - Verdict: unclear. It is settled either by a Coleman–Thun diagram for that pole or by enlarging the spectrum (b₋₁,₂ then forces b₀,₂, and so on).
  - Line 190: "Delius and Gandenberger account for all of them" must be qualified. Add after the range display: "This range leaves one pole unexplained: on $b_{-n,n}$ the pole $p=\frac12+(2n-1)\frac B2$, which DG attribute to the creation of $b_{-n+1,n}$, is in the strip for $B<\frac2{2n-1}$, while $b_{-n+1,n}$ lies in the range only for $B<\frac1{2n-1}$ (at $B=\frac9{20}$: the pole $\frac{47}{40}$ on $b_{-2,2}$)."
- **Line 202: unsupported.** "The lower bound $n+m\ge0$ is forced: with a state below $b_{0,0}$ some poles would be left unexplained" is not DG's argument. DG (p. 22) impose n+m ≥ 0 so that no state lies below the ground state; for a₄ they say "which we postulated not to exist". Replace with: "The lower bound $n+m\ge0$ is imposed so that no state lies below the ground state $b_{0,0}$ (for instance $e_{b_{-1,0}}-e_{b_{0,0}}=-m_1\cos(\frac\pi6+\delta)<0$) [@delius_gandenberger1999, sec. 6.2]."
- **Rotor: plausible, and quantitatively confirmed at leading order, with a framing correction** (`r6_rotor.py`).
  - The collective-coordinate inertia I(χ) = ∫_{−∞}^0 |∂_χφ|² dx of eq-vacuum-pair (a₂, a = 1) gives β²I = 3.4641 = 6/m₁, independent of χ.
  - The rotor reading of DG's tower, n²/(2I) = m₁(πB/6)n² with B ≈ β²/2π, needs exactly β²I = 6/m₁.
  - **Caveat:** for a fixed boundary condition χ lives in an open interval of length 2π/3 (`prp-solitonic-vacua`). A particle on that interval gives non-degenerate levels 9k²/(8I). The rotor match (period 2π, ±n degenerate) needs the full circle, i.e. the three intervals of the three quantum-equivalent patterns (+−−), (−+−), (−−+) glued together. That fits DG's finding that the three patterns have identical amplitudes.
  - Also, sin²(nδ)/sinδ departs from n² by up to about 20% at the top of the tower (nδ ≈ π/4).
  - Line 219: replace "This interpretation is ours and has not been checked by a collective-coordinate quantization." with "This interpretation is ours. The classical moment of inertia of $\chi$, $I=\int_{-\infty}^0|\partial_\chi\vec\phi|^2dx=6/(\beta^2m_1)$ (independent of $\chi$), is exactly what $n^2/2I$ requires at leading order. The match needs $\chi$ of period $2\pi$, whereas for a fixed sign pattern $\chi$ ranges over an interval of length $2\pi/3$; the rotor runs through all three quantum-equivalent patterns."

### 7. Open problems
1. **Problem 1 (∏C = −1): correct and open.**
2. **Problem 2: open in the soliton-derived sense, but cite the bootstrap candidates.** Castro-Alvaredo–Fring, NPB 682 (2004) 551, hep-th/0308004 (DOI 10.1016/j.nuclphysb.2004.01.009), propose generic boundary-bootstrap solutions for all affine Toda theories. Add: "Generic solutions of the boundary bootstrap for all affine Toda theories were proposed by Castro-Alvaredo and Fring [@cf2004]; which boundary conditions they describe is not known."
3. **Problem 3: open (no soliton-side derivation found).** Add Ahn–Kim–Rim, NPB 628 (2002) 486, hep-th/0110218, DOI 10.1016/S0550-3213(02)00133-5. Their abstract reports non-perturbative TBA checks of the Neumann and (+) amplitudes for A_n, of the boundary vacuum energies, and of duality. This also belongs in the evidence list at lines 62–68 (citation only, not checked).
4. **Problem 4: fine.**
5. **Problem 5: partially addressed.** Baseilhac–Belliard, LMP 93 (2010) 213, arXiv:0906.1215 (DOI 10.1007/s11005-010-0412-6), claim a classification of integrable fixed and dynamical (soliton non-preserving) boundaries via generalized q-Onsager algebras for coideals. Cite it in problem 5, and in problem 9 for dynamical boundaries.
6. **Problem 7: the framing understates the overlap.** On real configurations Im φ ≡ 0, so the Neumann half of eq-soliton-preserving is vacuous, and a Dirichlet wall φ = 0 satisfies the whole condition (0 is in the lattice). Replace the last sentence with: "On real field configurations, $\operatorname{Im}\vec\phi\equiv0$ and a Dirichlet wall $\vec\phi=0$ satisfies (-@eq-soliton-preserving) entirely, its Neumann half trivially; but that condition was derived for complex classical solutions, and whether the Euclidean theory with such a wall is the soliton-preserving integrable boundary is open."
7. **Problem 8:** fine (DG §7: "non self-dual" theories).
8. **Problem 10:** add BCZ, JHEP 01 (2004) 056, hep-th/0401020 (DOI 10.1088/1126-6708/2004/01/056), for the affine Toda defects that "Their affine Toda generalizations" refers to.
9. **Problem 12:** see item 4.
10. **Line 127, unclear wording:** "a third boundary, between Neumann and the uniform C≡+1" does not say in what sense "between". Say instead "with $\hat\epsilon_0=-\epsilon_*$, $\hat\epsilon_{j\ge1}=\epsilon_*$".

### 8. Status labels
- **`prp-uniform-amplitude`: keep [Numerical]** (all its items are identities of an explicit formula; analytic proofs exist via the block identities, so [Established] would also be defensible). Its title, however, makes a physical claim. Retitle it "Properties of the uniform amplitude" and add one line after it: "That (-@eq-uniform-amplitude) is the reflection amplitude of $C\equiv+1$ is @cnj-reflection-continuation applied to @prp-breather-amplitude."
- **`prp-solitonic-amplitude` and `prp-a2-excited`:** the same treatment, with "assuming @cnj-reflection-continuation" in the prose.
- **`prp-breather-amplitude` [Established]:** fine after the ch. 33 fix, and now matrix-checked for n = 2, 3.

### 9. Exercises
All are solvable, except the following:
- **`exr-uniform-classical`:** the telescoping needs (0) = 1. The hint (x)(−θ) = 1/(x)(θ) is unneeded. Suggested hint: "$(0)=1$ and $(h)=-1$."
- **`exr-kim-expansion`: ill-posed.** "Show that the result is self-dual to this order" has no meaning at O(β²), since B → 2−B is not small. Replace the sentence with: "Show that Kim's self-dual conjecture $K_1=[\frac12][\frac32][\frac52]$, $K_2=[\frac12][\frac32]^2[\frac52]$ (with $[x]=(x-\frac12)(x+\frac12)/((x-\frac12+\frac B2)(x+\frac12-\frac B2))$) has the same $O(\beta^2)$ term but differs at $O(\beta^4)$, and explain why one loop cannot test duality on the half-line." I checked this numerically.
- **`exr-uniform-crossing`, `exr-uniform-bootstrap`, `exr-shg-uniform`, `exr-solitonic-a2`, `exr-a2-energies`:** checked by hand, and the hints are right.
- **`exr-modulus-rotor`:** well posed. The answer is I = 6/(β²m₁), and the collective-coordinate integral reproduces it exactly (item 6). Add: "Note that $\chi$ is confined to an interval of length $2\pi/3$ for a fixed pattern; what period does the rotor need?"

### 10. Bibliography
All requested entries agree with INSPIRE: authors, title, journal, volume, page, year and DOI. That covers perkins_bowcock1999, corrigan1998, bcz2004 (IJMPA 19S2, 82), delius_gandenberger1999, gandenberger1999a, kim1995, cdrs1994, cdr1995, bowcock1998, dtw1999, gandenberger1998 and chenaghlou_corrigan2000.

Use against sentences: each supports its sentence, with these exceptions:
- the Kim-for-a₂ phrasing (item 3);
- "Corrigan observed …" (line 105): fair;
- DG fig. 9 for the triangle diagram: correct.

### 11. Prior art and later work
I listed the 28 INSPIRE citations of DG. By title and abstract, none reports the (6.11) or (6.26) misprints, the (3.11)/(3.9) prefactor issue, or a closed form of the a₂ boundary energies. I did not read all 28 in full.
- No later paper derives the Neumann amplitudes from the soliton side for n ≥ 3 (Gandenberger hep-th/9911178 has no relevant citing papers).
- Relevant later work missing from the chapter: Ahn–Kim–Rim 2002, Castro-Alvaredo–Fring 2004, Baseilhac–Belliard 2010 and BCZ 2004 JHEP (bibliographic data in item 7). Also Fateev hep-th/0103014 and Fateev–Onofri NPB 634 (2002) 546 (boundary one-point functions and vacuum energies), which could be cited alongside Ahn–Kim–Rim.

### Scripts (all in the scratch directory)
| Script | What it checks | Result |
|---|---|---|
| `r1_scalar_product.py` | A·F₁₁·A (DG prefactor) against (3.11), n = 1–3 | ratio 1 to 1e-10 |
| `r1b_prefactor.py` | DG vs book prefactor, n = 1–4; A_book/A_DG = f | DG exact; book off by f f; A_DG(0) = −1, A_book(0) = +1 |
| `r1c_f11_sg.py` | F₁₁(h = 2) against sine-Gordon S₀ | 1 to 1e-10 |
| `r2_matrix_bootstrap.py` | n = 2, 3 singlet bootstrap, book Ř | c = 1, deviation ~1e-15 (16 cases) |
| `r3_kim.py` | Kim (17)/(42) against DG (3.20) | equal at O(B), differ at O(B²) |
| `r5b_perkins_bowcock.py` | PB (5.5) and (4.40)+α(5.3) against DG n = 2 | ratio exactly i; with i restored, equal to 1e-18 |
| `r6_poles.py` | classification of all strip poles at B = 9/20 | only 47/40 on b₋₂,₂ / b₂,₋₂ unexplained |
| `r6_rotor.py` | β²I(χ) of the a₂ modulus | 3.4641 = 6/m₁ for all χ |

Rendered PDF pages used for checking printed equations (DG pp. 4, 11, 12, 25, 27, 35; Kim pp. 16–17; PB p. 16) are also in that directory.

**Not done:** item 5c (a₄ excited amplitudes away from the vacuum) and a full reading of the 28 citing papers.
