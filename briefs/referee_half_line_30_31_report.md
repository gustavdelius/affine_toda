# Referee report: chapters 30–31 (classically integrable boundary conditions; solitons on the half-line)

Report of the referee pass of 2026-10-11 on `referee_half_line_30_31.md`. Line numbers refer to the files as of commit 64ed4af. The referee's scratch scripts were in the session scratchpad (`referee_30_31/`):
- `own/ref_vacua_breathers.py`
- `own/ref_ksolve.py`
- `own/ref_d43*.py`
- `own/ref_delius15.py`
- `rmatrix/rmatrix_check.py`, `rmatrix/extra_check.py`

The integration status is recorded at the end.

**Headline: new errors found**
1. BCDR App. D for d₄^(3) is wrong, a third wrong entry; the book quotes it unchecked. The branch "C₀=C₁=0, C₂ (long) arbitrary" has no invertible K.
2. The Delius (1.5) note "integrable only for a_n and c_n" is wrong. (1.5) is also integrable for a_{2n}^(2) and d_{n+1}^(2). Delius himself states the a_n/c_n part, so this is not a correction of his paper.
3. The c_n vacua statement (ch. 31 l.128, ch. 35 l.230), "C₀C_n=1 ⇒ continued vacuum", drops Delius's assumption that every weight is a soliton charge. It fails for c₂.
4. K^cl poles (ch. 31 l.211): for a>h/2 the physical-strip pole is at iπ(h−a)/h, so "only for a<h/2" is wrong.
5. exr-dg-theorem: its second pattern has ∏C=−1, so it is not solitonic.
6. BCDR radicals (ch. 30 l.315): the "displaced radicals" exist only in the pdftotext extraction; the typeset paper is unambiguous.
7. Prior art: `prp-self-conjugate-vacuum` is in Bowcock–Perkins hep-th/9909174 §4.2, which is neither cited nor in the bib.

## Item 1. Corrections to the sources
**1a. g₂^(1) "C_short arbitrary" is false. Verdict: correct.**
- Checked independently in the 14-dim adjoint of g₂, with the Cartan matrix verified and node 2 short.
- All 8 sign patterns and C=0 have nullity 1 with an invertible K.
- (1,1,0.3), (1,−1,2−i), (−1,1,0), (1,1,0.5) and (0,0,1) have nullity 0. The smallest relative Gram eigenvalue is 1.7e−4, against a threshold of 1e−10.

**1b. c₂^(1)=b₂^(1) extra branch. Verdict: correct.**
- Checked in the 10-dim adjoint of sp(4).
- (1,0.4,−1), (−1,2+i,1), (1,0,1), (0.4,0,−1.7) and (2+i,0,0.3) pass; (1,0.5,0.9), (0.3,1,−1), (0,1,0) and (1,0.5,0) fail.
- The union rule for a₃^(2)≅d₃^(2) rests only on the book's script; it was not re-run independently.

**1c. BCDR (1.4), (3.18). Verdict: the reading is correct, the note is wrong.**
- The rendered arXiv PDF prints |A_i| = 2√n_i in (1.4), and √(n_i/2|α_i|²) in the line after (3.18). Nothing is displaced.
- The journal version was not accessible.
- Replace ch. 30 l.315 with:
  > - Appendix C of [@bcdr1995] normalizes $C_0$ of $a_2^{(2)}$ differently ($B_0=-\sqrt2C_0$); this is harmless. Their (1.4), $|A_i|=2\sqrt{n_i}$, is our $|C_i|=1$ via $A_i=2C_i\sqrt{2n_i/\vec\alpha_i^2}$.
- In `boundary_classical.py`, the control labelled "literal |A_i|=√(2n_i) reading of BCDR (1.4)" tests a reading that does not occur. Relabel it as a normalization control.

**1d. Delius 98a (1.5) "integrable only for a_n, c_n". Verdict: wrong, and not a correction.**
- Delius says this himself (p. 9, abstract, §6 item 2): (1.5) coincides with (1.11) for a_n and c_n.
- `ref_delius15.py` sets C_j ∝ √(n_jα_j²/2), with long roots |α|²=2 and ε=±1. An invertible K exists for a₂^(2), a₄^(2), a₆^(2), d₃^(2), a₂^(1) and c₃^(1). It fails for b₃^(1), d₄^(1) and g₂^(1).
- Replace ch. 30 l.316 with:
  > - As Delius notes [@delius1998a, p. 9], the boundary condition (1.5) of [@delius1998a], with a single $\epsilon$ multiplying $\eta_i\vec\alpha_i$, is an integrable condition of @tbl-durham for $a_n^{(1)}$ and $c_n^{(1)}$, where $\eta_i\vec\alpha_i^2$ is the same for all nodes. It is also integrable for $a_{2n}^{(2)}$ and $d_{n+1}^{(2)}$, where the nodes with a different $\eta_i\vec\alpha_i^2$ carry free parameters. For $b_n^{(1)}$ ($n\ge3$), $d_n^{(1)}$, $e_n^{(1)}$, $f_4^{(1)}$, $g_2^{(1)}$, $a_{2n-1}^{(2)}$ ($n\ge3$), $e_6^{(2)}$ and $d_4^{(3)}$ the factor $\sqrt{2n_j/\vec\alpha_j^2}$ is required (checked for $b_3$, $d_4$, $g_2$, $a_{2,4,6}^{(2)}$, $d_3^{(2)}$).

**1e. Delius (3.36) tan²→cot² and (3.37) cosh→sinh. Verdict: correct, now proved in general (see item 6).**
- (3.37) is a misprint that is consistent with his own (3.38), |ρ₋| < arcsinh√|X|; worth saying.

**1f. Delius p. 7 energy argument. Verdict: correct.**
- His difference is the boundary energy of the vacuum itself, which is present with or without the soliton.
- Optional addition to ch. 31 l.65 / l.258: "(the difference Delius finds is the boundary energy of the vacuum, $\pm2h/\tilde\beta^2$ in his normalization, which is present without the soliton)".

**1g. Delius (4.4) is the singular partner's energy. Verdict: correct.**
- Bowcock (3.17)–(3.18) gives E± = −2(n+1)(m_a ± 1), and the good solution is C₋.
- The book's E_vac becomes +2(n+1)(1−m_a)/β² in Delius's normalization.
- Delius's "It is negative" is false for m_a < 1.
- Suggested addition to ch. 31 l.259: "; in his normalization the energy of his regular solution is $+2(n+1)(1-m_a)/\beta^2$, which is positive for $m_a<1$."

**1h. Delius §1.3. Verdict: correct.** The same constraint reappears as Delius (6.4); the note could cite both.

**1i. Delius 98b (1.4), coweight in 2Λ^∨. Verdict: correct.**

## Item 2. Status labels and new claims
**`prp-self-conjugate-vacuum` [Theorem]. Verdict: statement correct, proof incomplete, and not new.**
- Energy 0 is checked for a₁, a₃, a₅, a₇ (|E| < 2e−10).
- General proof, to replace the last two sentences of the proof at ch. 31 l.177:
  > For $\xi<0$ real, $|u_j|<1$ on $x\le0$, so all $\tau_j>0$. The energy is the same for every member of the family: for any one-parameter family of static solutions obeying the same boundary condition, $\frac{dE}{d\xi}=\big[\partial_\xi\vec\phi\cdot(\partial_x\vec\phi+\nabla\mathcal B)\big]_{x=0}=0$ by the static field equation, the boundary condition, and $\partial_\xi\vec\phi\to0$ at $-\infty$. As $\xi\to-\infty$ the field tends uniformly to $0$, whose energy is $0$. (Checked by quadrature for $a_1,a_3,a_5,a_7$ and symbolically for the boundary condition in `code/boundary_vacua.py`.)
- The same argument proves that `eq-vacuum-energy` is independent of χ.
- **Prior art:** Bowcock–Perkins hep-th/9909174 (JHEP 02 (2003) 016), §4.2 "Single middle soliton", (4.11), has τ_j = 1 + (−1)^j d e^{2x} with d ≤ 1, all A_i positive (C≡−1), and the same energy as the constant solution. At d=1 the odd τ_j vanish at x=0 and further sign patterns become possible. Relabel as **[Established; reformulated]** [@bowcock_perkins2003, sec. 4.2].
- The zero-mode interpretation is correct: for a=h/2 and C=−1, the pole of K^cl at θ=iπ/2 is a static mode ∝ ρ_{h/2}e^{2mx} = ∂_ξφ at ξ→−∞.
- DG §4/§7 (p. 38) speak instead of "two different vacuum sectors". The single-family view is defensible, but mention DG's reading.

**`prp-soliton-preserving` [Theorem]. Verdict: complete for the compatibility part.**
- That K=1 solves the reflection equation for σ=† is now checked: r^{t₁t₂}(z) = −r(1/z) for real z, residual about 1e−15.
- Caveat: for antilinear σ and complex fields, Prop. 1's proof does not apply as written, because {T, T̄} is not given by r. The K=1 conclusion still holds by direct check.
- Suggested sentence after the proof (ch. 30 l.227): "For $K_+=1$ the reflection equation reduces to $r(z)+r^{t_1t_2}(1/z)=0$ and $r^{t_2}(w)+r^{t_1}(1/w)=0$, which hold for the Toda r-matrix at real spectral parameters."

**`prp-durham` [Numerical]. Verdict: appropriate.**
- The "only if" direction for simply-laced algebras is BCDR's analytic order-λ² result; cite it as such.
- With item 3 more rows are checked, so ch. 30 l.6 "two corrections" should read "three corrections".

## Item 3. Unchecked rows of `tbl-durham` (independent K-solve, `ref_ksolve.py`)
**Useful reductions (add to the text):**
- K ↦ DKD⁻¹ with D = e^{iπμ·H}, μ a coweight, flips C_j ↦ (−1)^{α_j·μ}C_j.
- K ↦ K⁻¹ maps C ↦ −C.
- So for a simply-laced row only the classes of ∏C_j^{n_j} need checking.

| row | result |
|---|---|
| e₆^(1) (27) | **Pass.** Uniform +1, 5 random patterns covering both classes ∏C^n = ±1, and C=0 have nullity 1 with an invertible K. Four off-list values fail. |
| f₄^(1) (e₆ Z₂ fold) | **Pass, both branches.** Tested in the 27 and in the irreducible 26 (Cartan matrix verified; nodes 3, 4 short). (C₀,C₁,C₂)=±1 with C₃=C₄=0 passes. Off-list values fail. |
| **d₄^(3)** (e₆ Z₃ fold) | **App. D is wrong.** All ±1 pass in both faithful summands (16 and 8). The branch C₀=C₁=0 with C₂ (long, verified) arbitrary fails for C₂ = 0.7, 2+i, +1 and −1, with nullity 0 in both summands. Only C≡0 and all ±1 are integrable. |
| a₅^(2), a₇^(2) | **Pass.** Uniform and random ±1 pass, as do (0,…,0,0.7) and (0,…,0,2+i). Off-list values fail. |
| e₆^(2) (e₇ Z₂ fold, 56) | **Unchecked.** The K-solve was still running at hand-back. |
| e₇^(1), e₈^(1) | **Unchecked.** With the reductions above, only C≡+1 and one pattern with ∏C^{n_j}=−1 need checking. |

**Corrected d₄^(3) row (ch. 30 l.160):**
`| $d_4^{(3)}$ | **all $\pm1$ only** | checked (27 of $e_6$, both faithful summands); **corrects [@bcdr1995]** |`

**Add to the list at l.168–171:** "- For $d_4^{(3)}$, Appendix D also allows $C_0=C_1=0$ with the long-root parameter arbitrary. The full equation has no invertible solution on this branch (not even for $C_2=\pm1$); only Neumann and all $\pm1$ survive." Add a matching bullet to the notes (l.311–317). Change rows e₆^(1), f₄^(1) and a_{2n−1}^(2) (n=3,4) to "checked".

## Item 4. Sklyanin's formalism (r-matrix sub-agent; key identity re-derived by the referee)
**Verdict: correct.**
- The transcription (l.96–103) matches Delius (2.7).
- The r-matrix from canonical brackets is (β²/4)·Delius (3.4). It satisfies r^{t₁t₂}(z) = −r(1/z) = r₂₁(z) and the classical Yang–Baxter equation.
- The Durham K solves (2.7), and K₋ = K₊(−1/λ) solves (2.8). Residuals are ≤1.7e−14 for a₂ and ≤1e−13 for a₁ with generic complex C.
- {τ,τ} vanishes to 1.7e−12.

**Needed fixes:**
- (i) Add (2.8) to `prp-sklyanin` (b):
  ```latex
  and for $K_-$
  $$0=\overset1K_-(\lambda_1)\overset2K_-(\lambda_2)\,r^{\sigma_1\sigma_2}\big(\tfrac{\lambda_1}{\lambda_2}\big)+r\big(\tfrac{\lambda_2}{\lambda_1}\big)\overset1K_-(\lambda_1)\overset2K_-(\lambda_2)+\overset1K_-(\lambda_1)\,r^{\sigma_1}(\lambda_1\lambda_2)\overset2K_-(\lambda_2)+\overset2K_-(\lambda_2)\,r^{\sigma_2}\big(\tfrac1{\lambda_1\lambda_2}\big)\overset1K_-(\lambda_1),$$
  ```
- (ii) Define r at ch. 30 l.103:
  ```latex
  Here $\sigma$ is linear and $\{a_x(\lambda_1;x)\overset\otimes,a_x(\lambda_2;y)\}=\delta(x-y)[r(\lambda_1/\lambda_2),\overset1a_x+\overset2a_x]$, with $r(z)=\frac{\beta^2}{4(z^h-1)}\big[(z^h+1)\sum_aH_a\otimes H_a+2\sum_{\alpha>0}\tfrac{\alpha^2}2(z^{{\rm ht}\alpha}E_\alpha\otimes E_{-\alpha}+z^{h-{\rm ht}\alpha}E_{-\alpha}\otimes E_\alpha)\big]$ [@delius1998b, (3.4)], which satisfies $r^{t_1t_2}(z)=-r(1/z)$.
  ```
- (iii) At l.114 replace "quoted … not checked … unchecked" with: "Part (b) follows from the bracket of $T$ and the Leibniz rule for linear $\sigma$ [@delius1998b, Prop. 1] (re-derived); no property of $r$ beyond that bracket is used."
- (iv) Replace l.175:
  > A solution $K$ of (-@eq-boundary-k-equation) also solves (-@eq-classical-reflection) [@bcdr1995; @delius1998b]; checked numerically for $a_1^{(1)}$ (generic complex $C_0,C_1$) and $a_2^{(1)}$ (sign patterns). If the boundary condition at $x_-$ is the same equation $\partial_x\vec\phi=-\nabla\mathcal B_C$, then $K_-(\lambda)=K_C(-1/\lambda)$ solves the condition at $x_-$ and the $K_-$ equation, because $a_t^T(-1/\lambda)=a_t(\lambda)$. (Delius writes $a_t^t(-\lambda)=a_t(\lambda)$, a misprint.) For a left boundary with energy $+\mathcal B_C$, i.e. $\partial_x\vec\phi=+\nabla\mathcal B_C$, take $K_-=K_C(\lambda)^T\propto K_{-C}(-1/\lambda)$. The parameters at the two ends are independent.

## Item 5. Real-coupling vacua
| claim | verdict |
|---|---|
| Lowest-energy assumption | Correctly flagged; unproven. A proof would need a lower bound on E including the boundary term. |
| Bowcock §3.2 multi-soliton vacua | Unchecked. To settle: reproduce them for a₃. Bowcock–Perkins §4 also has backgrounds singular exactly at x=0 (d=1) that satisfy further patterns. |
| Bowcock (4.21), DG (5.7) | Unchecked. To settle: linearize about `eq-vacuum-pair` for a₂. |
| **c_n condition (ch. 31 l.128; ch. 35 l.230)** | **Wrong as stated** (see below). |
| `tbl-missing-vacua` | Correct; counts reproduced: 1, 3, 5, 10, 10, 21, 21. Add to the caption: "(single-soliton; multi-soliton backgrounds [Bowcock §3.2; Bowcock–Perkins §4] not considered)". |

**The c_n condition.**
- Delius (p. 37) derives A₀Aₙ=1 from (6.4), which assumes every weight is a soliton charge.
- Counting the folding-symmetric solitonic patterns that have a single-soliton a_{2n−1} vacuum gives 1 of 3 for c₂ (only C≡−1), at most 2 of 7 for c₃, and at most 1 of 15 for c₄.
- For c₂, the patterns (+,−,+) and (−,+,−) fold from exactly the missing alternating a₃ patterns.
- For a ≠ n the continued pair is folding-invariant only at χ ∈ πZ.
- Replace the last sentence of ch. 31 l.128 with:
  > Delius [-@delius1998a, sec. 6] infers from (6.4) that the real-coupling $c_n$ theory has a continued vacuum whenever $C_0C_n=1$, assuming that every weight is a soliton charge. With the solitons actually known this fails: for $c_2^{(1)}$ the patterns $(+,-,+)$ and $(-,+,-)$ fold from the alternating $a_3^{(1)}$ patterns of @tbl-missing-vacua, and only $C\equiv-1$ (the self-conjugate species-2 soliton) has a continued vacuum. A stationary pair with $a\ne n$ is folding-invariant only for $\chi\in\pi\mathbb Z$.
- Replace ch. 35 l.230 accordingly: "continued vacua exist for some patterns with $C_0C_n=1$".

**K^cl poles (ch. 31 l.211).** Replace with:
> - **$C=-1$.** $K_a^{\rm cl}$ has a pole in the physical strip at $\theta=i\pi\min(a,h-a)/h$ for every $a\ne h/2$, corresponding to classical boundary bound states with $\omega_a=m_a|\cos\frac{\pi a}h|=m|\sin\frac{2\pi a}h|$ [@cdr1995]. For $h$ even and $a=h/2$ the pole sits at $\theta=i\pi/2$, $\omega_{h/2}=0$: …

**"Attractive binds" (ch. 31 l.65, l.119): unclear.**
- The C=+1 breathers are singular exactly at x=0.
- Linearly it is the C=−1 boundary that has bound states.
- Suggest: "For $C=+1$ … there are boundary breathers, singular only at $x=0$ itself"; and at l.119: "the soliton-attractive $C=+1$ boundary has (marginally singular) breathers, the repulsive $C=-1$ boundary none; for small fluctuations the roles are reversed (@sec-classical-reflection-factors)."

## Item 6. Boundary breathers, all species and ranks
**Verdict: correct for all a and h.** An analytic proof is below; it can replace [Numerical] with [Theorem].
- **Neumann, A<0.** Choose e^{S/2} = i/√|A|, y = m_a x cosη, ψ = ζ₋ − ωt + 2πja/h. Then τ_j = 1 + e^{2y} − 2R sinhρ₋ sinψ + 2iR coshρ₋ cosψ with R = e^y/√|A|. A zero needs cosψ = 0 and cosh y = |sinhρ₋|/√|A|. So the breather is regular on x ≤ 0 ⟺ |sinhρ₋| < √|A|.
- **Neumann, A>0.** A zero exists at cosh y = coshρ₋/√A > 1, so the breather is always singular.
- **C=+1, ρ₋=0.** min_ψ τ_j = [(1−e^y)² + s(1−e^{2y})]/(1+s), which is > 0 for y<0 and = 0 at y=0. For ρ₋≠0, τ = 2−2coshρ₋ < 0 at y=0 and 1+s > 0 at −∞, so there is a zero.
- **C=−1.** There is a sign change on y<0 for both s<1 and s>1.
- Numerical scan: h = 3…7, all a, several u and ρ₋, with 0 violations.

## Item 7. Heuristic "σ selects charges"
**Verdict:** this can be made precise at the level of asymptotic charges, but not yet from τ(λ).
- If the conserved combinations are I_s = Q_s + Q_{−s}, then q_s(a) = q_{−s}(b).
- With q_{−s}(b) = q_s(b̄) this forces b = ā; an extra complex conjugation (σ=†) forces b = a.
- What is missing: proof that τ(λ) generates exactly these combinations, i.e. abelianization with a boundary.
- Keep the "heuristic" label and add this sentence.

## Item 8. Dictionary `tbl-boundary-dictionary`
- **SSW (ζ,η): correct.** Their (2.2) with φ_SSW = β_SGφ gives M = β²M₀/2m₀ and C₀ = cosh(ζ+iη). Their text's "−M cos" potential is off by a factor 2 relative to (2.2), and the book correctly uses (2.2).
- **SSW (2.15), partially checked.**
  - At (0,π/2) it gives a = ln(1−s²), the book's Neumann shift.
  - At (0,0) it gives 2ln(1−s), the book's C=−1 value, not C=+1. This is because SSW's ansatz (2.6) contains a stationary kink that shifts the boundary vacuum by a coweight, so C↦−C. Note this if (2.15) is quoted.
- **BDG.** Their equations match with BDG's convention ṗ = {H, p̂}. State this convention in `eq-boundary-oscillator` and in `exr-oscillator-energy`, since the signs flip with ṗ = {p, H}.
- **BCR.** "Only odd order automorphisms … at most a_r and e₆" is confirmed; "discrete solutions, no free parameter" is unverified.
- **Baseilhac–Koizumi.** Not checked.

## Item 9. Exercises
- **`exr-dg-theorem` (ch. 31 l.247): ill-posed.**
  - (−,+,−,−,+) has ∏C = −1, so it is not solitonic.
  - The first pattern, (+,−,−,+,+), gives λ = λ₂ − λ₁ ∈ Λ₁, and its vacuum exists.
  - Replace the second pattern with (−,+,−,+,+): λ = λ₂ ∈ Λ₂, with no single-soliton vacuum, which gives a good "missing vacuum" contrast.
- **`exr-linear-reflection`:** its answer contradicts l.211 (see item 5).
- **`exr-oscillator-energy`:** needs the bracket convention (item 8).
- **All other exercises** are well posed and their hints correct, including `exr-sg-k-classical`, `exr-vacuum-boundary-value`, `exr-sinh-gordon-bps` and `exr-breather-threshold`.

## Item 10. Bibliography
- All 15 entries match INSPIRE and Crossref; bdg2002 is unpublished per INSPIRE, so @misc is right.
- **Missing entry:**
  ```
  @article{bowcock_perkins2003,
    author = {Bowcock, P. and Perkins, M.},
    title = {Aspects of classical backgrounds and scattering for affine Toda theory on a half line},
    journal = {JHEP}, volume = {02}, pages = {016}, year = {2003},
    eprint = {hep-th/9909174}, archivePrefix = {arXiv}, note = {arXiv:hep-th/9909174},
    doi = {10.1088/1126-6708/2003/02/016}
  }
  ```
- One misquote: the c_n sentence (item 5).
- Not checked: CDR (4.7), (4.9), (5.3) (the brief says the script verifies them); Bowcock (4.21) and §3.2; Delius (3.19) and (3.26); Baseilhac–Koizumi.

## Item 11. Prior art
- **Searched:**
  - INSPIRE citers of Delius 98a (17) and Bowcock (14);
  - BCDR citers (121), filtered by title;
  - full texts of Bowcock–Perkins hep-th/9909174, Corrigan's review hep-th/9612138, Penati–Refolli–Zanon hep-th/9510084 and Fateev–Onofri hep-th/0207152.
- **Found:** `prp-self-conjugate-vacuum` is in Bowcock–Perkins §4.2.
- **Not found in these:** the g₂ and d₄^(3) App. D errors, the c₂ branch, the cot²/sinh slips, the (4.4) partner energy, the 98b 2Λ^∨ point.
- The search is not exhaustive. Bowcock–Perkins was read only around §4, and the Rietdijk thesis was not checked.

## Integration (2026-10-11)

**Integrated** (ch. 30, ch. 31, one sentence of ch. 35; checks in `boundary_classical.py`, `boundary_images.py`, `boundary_vacua.py`):
- 1a, 1b: already in the text; unchanged.
- 1c: BCDR radicals note replaced; the `boundary_classical.py` control is relabelled as a normalization control (overall factor 1/√2).
- 1d: Delius (1.5) note replaced. New check: (1.5) is integrable for a₂, c₃, a₂^(2), a₄^(2), a₆^(2), d₃^(2) and not for b₃, d₄, g₂.
- 1e–1h: added to ch. 31 (breather paragraph and notes). 1i was already in the text.
- `prp-self-conjugate-vacuum`: the general dE/dξ proof is added, and the proposition is relabelled **[Established; reformulated]** [@bowcock_perkins2003, sec. 4.2]. Bowcock–Perkins (4.11) is described, with the d=1 remark and DG's "two vacuum sectors" reading. The χ-independence remark is added after `prp-solitonic-vacua`. The energy check now covers a₁, a₃, a₅, a₇.
- `prp-soliton-preserving`: the K₊=1 reduction sentence is added, with a check, plus the caveat about antilinear σ.
- `prp-durham`: the "only if" direction for simply-laced algebras is cited to BCDR. "Two corrections" now reads "three corrections" in the learning goals and the list.
- Item 3: the d₄^(3) row is corrected, with a bullet in the text and in the notes. The e₆^(1), f₄^(1) and a_{2n−1}^(2) (n=3, 4) rows are marked "checked". The coweight/inverse reductions sentence is added.
  - New checks: e₆^(1) in the 27 (uniform, 5 random patterns covering both classes, C=0; 4 off-list values fail); the reductions K⁻¹ and DKD⁻¹ for 6 fundamental coweights; f₄^(1) folded 27; d₄^(3), both faithful summands (8, 16); a₅^(2), a₇^(2).
  - e₆^(2) stays unchecked: the referee's K-solve hit its 1500 s timeout without a result.
- Item 4 (i)–(iv): the K₋ equation, the definition of r, the part (b) text and the replacement of l.175. New checks for a₁ (generic complex C) and a₂ (all sign patterns):
  - r is the r-matrix of the canonical brackets, and r^{t₁t₂}(z) = −r(1/z);
  - the Durham K solves (2.7); K₋ = K₊(−1/λ) solves (2.8) and the x₋ condition; K_C(λ)^T ∝ K_{−C}(−1/λ);
  - {τ,τ} = 0 from {T,T} = [r, TT] alone, with a random-K₊ control.
- Item 5:
  - c_n sentence replaced in ch. 31, with a notes bullet; ch. 35 sentence changed to "for some patterns with C₀C_n=1". New checks: the c₂ count (1 of 3), c₃ (2 of 7) and c₄ (1 of 15); the pair is folding-invariant only for χ ∈ πZ.
  - K^cl pole statement replaced, with min(a,h−a); the check now covers all a, the a=h/2 edge, and no poles for C=+1.
  - `tbl-missing-vacua` caption extended.
- "Attractive binds" wording changed at both places.
- Item 6: `prp-boundary-breathers` is now a **[Theorem]** for all a_n^(1) and species, with the analytic proof. New grid scan for h=3..7, all a, which also checks the closed form of τ_j used in the proof.
- Item 7: the asymptotic-charge sentence is added to the σ heuristic in ch. 30; the heuristic label is kept.
- Item 8: the SSW (2.15) remark is added to the ch. 31 notes. The BDG convention ḟ={H,f} is stated after `eq-boundary-oscillator` and in `exr-oscillator-energy`.
- Item 9:
  - `exr-dg-theorem`: second pattern changed to (−,+,−,+,+).
  - `exr-linear-reflection`: hint added, consistent with the new pole statement.
- Item 10: `bowcock_perkins2003` is cited but **not yet in `references.bib`**. The entry was verified against INSPIRE (recid 507731) and must be added.

**Not integrated** (unclear or unchecked in the report):
- the lowest-energy assumption;
- Bowcock §3.2 multi-soliton vacua and Bowcock–Perkins §4 singular backgrounds;
- Bowcock (4.21) and DG (5.7);
- BCR "discrete solutions, no free parameter";
- Baseilhac–Koizumi;
- the e₆^(2), e₇^(1) and e₈^(1) rows;
- CDR (4.7)/(4.9)/(5.3) re-verification beyond the existing script;
- Delius (3.19)/(3.26).
