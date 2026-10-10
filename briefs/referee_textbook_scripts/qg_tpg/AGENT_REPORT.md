# Referee report: Part II, chapters 7–8 (items 4, 11, 12')

No repository file was edited. Scripts and outputs are in `/tmp/claude-1000/-home-gustav-Git-affine-toda/e6f669a7-0ff6-43bc-8999-8d711276f81f/scratchpad/qg_tpg/`; their full content and output are in section E.

## Summary of real errors (most important first)

1. **ch. 8:73 ("When it fails") is wrong.** I checked this numerically.
   - Under the horizontal subalgebra, which is U_q(b_n), the tensor square of the U_q(d_{n+1}^{(2)}) spinor is multiplicity-free for every n: ⊕_{k=0}^n Λ^k C^{2n+1}. Chapter 23 (line 20) says the same, so the two chapters contradict each other.
   - `eq-tpg-rule` as stated (untwisted, ad-based, edges only between opposite parities) gives the wrong answer for this twisted case.
   - The correct twisted rule needs ⟨a⟩₊ factors on edges between equal-parity components (Delius–Gould–Zhang 1996). Their §4.3 already gives this spinor R-matrix for all n.
   - The multiplicity failure is real for the *other* c_n^{(1)} soliton multiplets (Kirillov–Reshetikhin modules, reducible under U_q(b_n)), not for the spinor.
2. **ch. 7:109: the central element is ∏ k_i^{n_i} (Kac labels), not ∏ k_i^{n_i^∨}.** The printed element is not central for any non-simply-laced or twisted algebra. For d_{n+1}^{(2)} it is not even in the algebra, because the dual labels there are half-integers.
3. **ch. 8:82: the example is in the wrong place.** The 26⊗26 of U_q(f_4^{(1)}) is multiplicity-free (1⊕26⊕52⊕273⊕324, hence 5 unknowns), and DGZ 1994 lists it as solved by the graph method. The non-multiplicity-free example is 27⊗27 of U_q(e_6^{(2)}): Σ mult² = 2²+3²+1+1+1 = 16, which matches appendix D.
4. **Prior art affecting ch. 8:11 and Part VI novelty claims** (outside my items, but material):
   - DGZ 1996 (q-alg/9508012), §4.3, eqs. 4.27–4.31, derives the U_q(D_{l+1}^{(2)}) spinor R-matrix V(λ_l)⊗V(aλ_l) for all l ≥ 2 with the twisted graph. For a = 1 it is ρ_k = ∏_{i=1}^{l−k}⟨i⟩_{(−1)^i} in their units. This is chapter 23's "closed-form R-matrix", which chapters 22 and 26 call "conjectured for all n". I confirmed the formula against direct solutions of Jimbo's equations for l = 2, 3, 4 in the book's conventions.
   - DGZ 1994 (hep-th/9405030) lists, according to its arXiv HTML, U_q(B_n^{(1)}) spinor⊗spinor (§4.2, eqs. 4.21–4.25) and F_4 V(λ_4)⊗V(λ_4). I did not re-derive these.
   - Places to re-check: ch22:32, ch22:98–99, ch23:4,10,141, ch24:4,77, ch26:162. The S-matrices (scalar factors, bootstrap) may still be new; the R-matrices are not.
5. **Bibliography.**
   - `dgz1994` has the wrong title. The real one is "On the construction of trigonometric solutions of the Yang–Baxter equation", Nucl. Phys. B 432 (1994) 377 (checked on the arXiv abstract page).
   - `dgz1996` has no title or journal. It is "Twisted quantum affine algebras and solutions to the Yang–Baxter equation", Int. J. Mod. Phys. A 11 (1996) 3415–3438 (arXiv abstract page).
   - The trigonometric tensor-product graph should also be credited to Zhang, Gould and Bracken, "From representations of the braid group to solutions of the Yang–Baxter equation", Nucl. Phys. B 354 (1991) 625. DGZ 1996 calls it the paper that introduced the method. I took this attribution from DGZ 1996 and did not read the 1991 paper.

Everything else checked is correct, apart from the minor items listed under each section.

---

## A. Item 4: tensor-product-graph rule beyond sl2^

**Method (all numerical, in the book's conventions).**
- Explicit e_i, f_i, k_i for i = 0..n, with k_i = q_i^{h_i}, q_i = q^{α_i²/2}, long roots of length² 2.
- Every defining relation checked, including q-Serre with q_i.
- Jimbo's equations Ř Δ_{x,y}(a) = Δ_{y,x}(a) Ř solved with the book's coproduct (Δe = e⊗k + 1⊗e, Δf = f⊗1 + k⁻¹⊗f) in the homogeneous gradation (e_0 → x e_0), at random complex q, x, y.
- Ř normalised to 1 on v₁⊗v₁; its eigenvalue multiset compared with the rule, ⟨a⟩ = (1 − z q^{2a})/(z − q^{2a}).

### (a) U_q(a_2^{(1)}) vector representation: correct

V⊗V = 6(+) ⊕ 3(−). The rule gives ρ_{λ₂} = ⟨1⟩. It matches; the solution space is 1-dimensional and the relation residuals are about 1e-16.

### (b) Non-simply-laced vector representations: correct, in the book's normalisation

Checked for c₂ (C⁴), b₂ (C⁵), and also c₃ (C⁶) and b₃ (C⁷).

**c_n^{(1)}** (α_i = ε_i − ε_{i+1} short with q_i = q^{1/2}; α_n = 2ε_n and α_0 = −2ε_1 long with q_i = q; ε·ε = ½):
- V⊗V = L(2λ₁) (+, dimension n(2n+1)) ⊕ L(λ₂) (−) ⊕ L(0) (−).
- Casimirs: C = 2n+2, 2n, 0.
- Graph: L(λ₂) — L(2λ₁) — L(0). There is no edge λ₂–0 because they have the same parity.
- Rule: Ř = P_{2λ₁} + ⟨½⟩ P_{λ₂} + ⟨(n+1)/2⟩ P₀.
- So the λ₂ channel carries q¹, i.e. q_short², and the singlet carries q^{n+1}. This matches Jimbo's ξ = k^{2n+2} with k = q^{1/2}. (Jimbo 1986 was not fetched; I compared with the form I remember. The actual verification is the direct solution.)

**b_n^{(1)}** (α_n = ε_n short with q_n = q^{1/2}; the spin-1-like e_n carries a factor [2]_{q^{1/2}}; α_0 = −ε₁ − ε₂):
- V⊗V = L(2λ₁)(+) ⊕ L(λ₂)(−) ⊕ L(0)(+).
- Casimirs: C = 4n+2, 4n−2, 0.
- Graph: L(2λ₁) — L(λ₂) — L(0), a chain; both edges change parity.
- Rule: Ř = P_{2λ₁} + ⟨1⟩ P_{λ₂} + ⟨1⟩⟨(2n−1)/2⟩ P₀, which matches Jimbo's ξ = k^{2n−1} with k = q.

All five cases pass with three random spectral points each.

**Negative controls on c₂ (both fail, as they should):**
- Casimir normalised with long roots of length² 4, i.e. ⟨1⟩, ⟨3⟩.
- The rule with z → 1/z.

So the book's statement of the rule (Casimir with long roots of length² 2, q = q_long, homogeneous gradation, ⟨a⟩ as printed) is right for untwisted algebras, including short-root channels. Optionally, add the c_n vector as a second worked example: it is the case where half-integer arguments ⟨½⟩ appear.

### (c) The c_n^{(1)} spinor claim at ch. 8:73: wrong

**Subalgebra.** The d_{n+1}^{(2)} diagram is 0 ⇐ 1 — … — (n−1) ⇒ n, with both end nodes short. Removing node 0 (or node n) gives **B_n**, never C_n. This is also g_{(0)} = so(2n+1), the fixed-point subalgebra of so(2n+2) as defined in ch. 6. So "U_q(g_{(0)})" is technically right, but it should be written U_q(b_n).

**Decomposition.** Spinor⊗spinor = ⊕_{k=0}^n Λ^k V (for k = n this is the component with highest weight 2λ_n). It is multiplicity-free for **all** n. Numerically:
- The spinor is a representation of U_q(d_{n+1}^{(2)}) for n = 2, 3, 4, with all relations including q-Serre holding to 4e-16.
- The commutant of U_q(b_n) on S⊗S has dimension n+1 (= 3, 4, 5).
- The Jimbo solution is unique.

**Eigenvalues.** They obey the DGZ 1996 twisted rule, which in book units reads:

ρ_{Λ^{n−j}} = ∏_{i=1}^{j} ⟨i/2⟩_{(−1)^i}, with ⟨a⟩₋ = ⟨a⟩ and ⟨a⟩₊ = (1 + z q^{2a})/(z + q^{2a}).

The chain edges alternate between opposite parity (⟨⟩₋) and equal parity (⟨⟩₊). Λ^{n−1} and Λ^{n−2} both lie in Λ²S, inside the same so(2n+2) irreducible. The book's untwisted rule (all signs −) **fails** for n = 2, 3, 4, and so does the sign-flipped variant.

**Where multiplicities really occur.** The other c_n^{(1)} soliton multiplets are KR modules, reducible under U_q(b_n): for example Λ¹⊕Λ⁰ (dimension 2n+2, Jimbo's D_{n+1}^{(2)} vector) and the c₃ 29-plet Λ²⊕Λ¹⊕Λ⁰ (ch. 23 §sec-sm-c3). Their tensor products have repeated components; (Λ¹⊕Λ⁰)^{⊗2} already contains Λ¹ twice.

**"Connected graph".** The graph as the book defines it (ad, opposite parity) is DGZ's *extended* graph. DGZ 1994 show that the true graph (pairs actually linked by e_0) is always connected when V(x)⊗W(y) is generically irreducible, and it is contained in the extended graph. So "connected graph" is not a separate requirement in the untwisted case; the real issue is consistency around loops. The sentence about "cases where the graph is disconnected [@dgz1996]" describes the twisted situation: there, an untwisted-style graph would miss the equal-parity edges.

**The gradation condition** ("a gradation in which U_q(g) acts independently of x") is not a limitation. Other gradations are related to the homogeneous one by an x-dependent conjugation and z → z^p.

**Correction for ch. 8:73.** Current text:
> "**When it fails.** The method needs three things: a multiplicity-free decomposition, a connected graph, and a gradation in which the finite algebra $U_q(g)$ acts independently of $x$. The extended tensor-product graph of Delius, Gould and Zhang handles some cases where the graph is disconnected [@dgz1996]. But for the representations relevant to several affine Toda theories the first condition fails. The soliton multiplets of $c_n^{(1)}$, which form the spinor representation of $U_q(d_{n+1}^{(2)})$, have tensor squares in which the same $U_q(g_{(0)})$-representation occurs more than once when $n\ge3$. Then $\check R$ restricted to each isotypic component is a matrix, not a number, and (-@eq-tpg-rule) does not apply."

Replacement:
> "**Twisted algebras.** The rule as stated is for untwisted $\hat g$, in the homogeneous gradation (other gradations are related to it by an $x$-dependent change of basis). For a twisted algebra $g^{(2)}$, $e_0$ transforms under $U_q(g_{(0)})$ not in the adjoint representation but in $g_{(1)}$ (@sec-twisted), and the graph changes accordingly [@dgz1996]: $\mu$ and $\nu$ are joined if $L(\nu)\subset L(\mu)\otimes g_{(1)}$, an edge may join two components of equal parity, and along such an edge $\langle a\rangle$ in (-@eq-tpg-rule) is replaced by $\langle a\rangle_+=(1+z\,q^{2a})/(z+q^{2a})$. The example that matters here is the spinor representation of $U_q(d_{n+1}^{(2)})$, which carries the spinor soliton of $c_n^{(1)}$. Under the horizontal subalgebra $U_q(b_n)$ its tensor square is $\bigoplus_{k=0}^n\Lambda^k\mathbb C^{2n+1}$, without multiplicities, the graph is the chain $\Lambda^n-\Lambda^{n-1}-\dots-\Lambda^0$, and $\rho_{\Lambda^{n-j}}=\prod_{i=1}^{j}\langle i/2\rangle_{(-1)^i}$ with $\langle a\rangle_-=\langle a\rangle$ [@dgz1996].
>
> **When it fails.** The graph defined above is what Delius, Gould and Zhang call the extended graph. The pairs of components that $e_0$ actually links form a connected subgraph of it whenever $V(x)\otimes V(y)$ is irreducible for generic $x/y$ [@dgz1994], so imposing the rule on every edge gives the R-matrix provided the relations are consistent around loops. What the method cannot do without is a multiplicity-free decomposition, and for several affine Toda theories this fails. The soliton multiplets of $c_n^{(1)}$ other than the spinor are Kirillov–Reshetikhin modules of $U_q(d_{n+1}^{(2)})$ that are reducible under $U_q(b_n)$, for example $\Lambda^1\oplus\Lambda^0$ and, for $c_3^{(1)}$, the 29-plet $\Lambda^2\oplus\Lambda^1\oplus\Lambda^0$ (@sec-sm-c3). Their tensor products contain the same $U_q(b_n)$-representation more than once. Then $\check R$ restricted to each isotypic component is a matrix, not a number, and no graph rule applies."

The ⟨a⟩₊ rule is verified here only for the d_{n+1}^{(2)} spinor, n = 2, 3, 4. Its general form for k = 2 is from DGZ 1996, which assumes k = 2; d₄^{(3)} needs cube roots of unity, as appendix D notes.

**Related one-line fixes.**
- **ch. 8:52.** Current: "For multiplicity-free tensor products, MacKay [@mackay1991] and Delius, Gould and Zhang [@dgz1994] showed that its content is captured by a graph." Replace with: "For multiplicity-free tensor products its content is captured by a graph, introduced by MacKay [@mackay1991] for rational R-matrices and by Zhang, Gould and Bracken [@zgb1991] for trigonometric ones, and developed by Delius, Gould and Zhang [@dgz1994; @dgz1996]." This needs a new bib entry `zgb1991`.
- **ch. 8:142.** Current: "extended to trigonometric ones by Delius, Gould and Zhang [@dgz1994; @dgz1996]". Replace with: "extended to trigonometric ones by Zhang, Gould and Bracken [@zgb1991] and generalised by Delius, Gould and Zhang to pairs of different representations [@dgz1994] and to twisted algebras [@dgz1996]".
- **ch. 8:11.** Current: "For the representations needed for the solitons of the non-simply-laced affine Toda theories they were not, and constructing them was the main obstacle to finding those S-matrices (@sec-sm-introduction)." Replace with: "By the mid-1990s the tensor-product graph had added many more, among them the spinor R-matrices of $U_q(b_n^{(1)})$ and $U_q(d_{n+1}^{(2)})$ [@dgz1994; @dgz1996]. For the soliton representations whose tensor products are not multiplicity-free they were not known, and constructing them was the main obstacle to finding those S-matrices (@sec-sm-introduction)."
- **ch. 8:82.** Current: "For the $\mathbf{26}\otimes\mathbf{26}$ of $U_q(f_4^{(1)})$ this reduces the problem to five unknowns per value of the spectral parameter." Replace with: "For the $\mathbf{27}\otimes\mathbf{27}$ of $U_q(e_6^{(2)})$, where $\mathbf{27}=\mathbf{26}\oplus\mathbf1$ under $U_q(f_4)$ and the components $\mathbf1$ and $\mathbf{26}$ occur twice and three times, this reduces the problem to $2^2+3^2+1+1+1=16$ unknowns per value of the spectral parameter." (26⊗26 = 1+26+52+273+324 is multiplicity-free.)
- **references.bib:85–93.** Fix `dgz1994` and `dgz1996` as in summary point 5.

---

## B. Item 11: chapter 7 statements

**KR remark (ch. 7:119): essentially correct; terminology fix.**
- Only for U_q(sl_n) does every irreducible representation extend (Jimbo's evaluation map).
- 27 = 26⊕1 for U_q(e_6^{(2)}) is consistent with the twisted KR rule: e_0 is in g_{(1)} = 26 = V(ω₄), and ω₄ − θ₀ = 0 adds a V(0).
- For untwisted F₄^{(1)}, ω₄ − θ is not dominant, so the 26 extends alone. This is consistent with Part VI's 26 of U_q(f₄^{(1)}).
- But "minimal extensions" of a general representation are *minimal affinizations*; KR modules are the case of multiples of a fundamental weight.
- Current: "the minimal extensions (the **Kirillov–Reshetikhin modules**) can be larger". Replace with: "the smallest representations of $U_q(\hat g)$ that contain a given one (its **minimal affinizations**; for a multiple of a fundamental weight, the **Kirillov–Reshetikhin modules**) can be larger".

**Central element (ch. 7:109): wrong.**
- Conjugation by ∏k_i^{m_i} multiplies e_j by q^{Σ_i m_i(α_i²/2)a_ij} = q^{(Σ_i m_i α_i)·α_j}. This is 1 for all j iff Σ m_i α_i = 0, i.e. m_i = n_i (Kac labels).
- Equivalently q^c = ∏_i q^{n_i^∨ h_i} = ∏_i (q_i^{h_i})^{n_i} = ∏ k_i^{n_i}.
- Explicit exponents for ∏k^{n^∨} (zero would mean central):
  - c₂^{(1)}: (1, −1, 1)
  - b₂^{(1)}: (1, 1, −1)
  - g₂^{(1)}: (0, 2, −4/3)
  - d₃^{(2)}: (−½, 1, −½), where n^∨ = (½, 1, ½)
  - a₂^{(2)}: (−¾, 3/2)
- For ∏k^n the exponents are all 0.
- In the c₂^{(1)} vector representation, k₀k₁²k₂ = 1, while k₀k₁k₂ = diag(q^{1/2}·…) is not a scalar.
- Current: "The element $\prod_ik_i^{n_i^\vee}$ is central, and on the finite-dimensional representations used in this book it equals 1 (level zero)." Replace with: "The element $\prod_ik_i^{n_i}$, with the Kac labels $n_i$ of (-@eq-kac-labels), is central: conjugation by it multiplies $e_j$ by $q^{(\sum_in_i\vec\alpha_i)\cdot\vec\alpha_j}=1$. Since $k_i=q_i^{h_i}$ and $n_i^\vee=n_i\vec\alpha_i^2/2$, it equals $q^{\sum_in_i^\vee h_i}$, $q$ to the power of the central element. On the finite-dimensional representations used in this book it equals 1 (level zero)."

**Quasitriangularity consequences (ch. 7:72, 83–85): correct.** Verified with the explicit universal R on spin ½ and spin 1:
- (Δ⊗1)R = R₁₃R₂₃; the order matters (R₂₃R₁₃ fails).
- (1⊗Δ)R = R₁₃R₁₂; R₁₂R₁₃ fails.
- R₁₂R₁₃R₂₃ = R₂₃R₁₃R₁₂.
- (S⊗1)R = R⁻¹ and (S⊗S)R = R.
- The orders match the printed coproduct, which is Kassel's convention.
- Also checked: R = PŘ satisfies `eq-yang-baxter` in the order R₁₂(x₁/x₂)R₁₃(x₁/x₃)R₂₃(x₂/x₃), consistent with ch. 7:138.
- Also checked: eq-evaluation-rep satisfies both q-Serre relations for spins ½ to 2. Algebraically, E³F − [3]E²FE + [3]EFE² − FE³ = 0 holds in U_q(sl₂).

**exr-hermitian-spin1: correct and solvable.**
- In the basis v_k = F^k v₀ the invariant form is unique up to a real scalar: diag(1, [2]_q, [2]_q²) = diag(1, 2cosγ, 4cos²γ).
- Signature is (3,0) for cos γ > 0 and (2,1) for cos γ < 0. At cos γ = 0 (q² = −1) the form is degenerate and span{v₁, v₂} is a subrepresentation.
- The representation is unitary iff cos γ > 0.
- The text's "for π/2 < γ < π this is negative" is right (generally π/2 < γ < 3π/2).

**Other chapter 7 points.**
- **ch. 7:111 (minor).** U_q(ĝ) is quasitriangular only after adjoining d and completing. Suggested text: "The algebra $U_q(\hat g)$, with the derivation adjoined and suitably completed, is quasitriangular too. … On $V(x)\otimes W(y)$ its universal R-matrix gives, up to a scalar function of $x/y$, the R-matrix of @sec-trigonometric-r."
- **ch. 7:172 (minor, not the crossing-point question the other agent covers).** "$V(x)^*\cong V(x\,x_c)$" holds only for self-conjugate V; for U_q(sl₃^) vector, V^* is the 3̄. Suggested: "$V(x)^*\cong\bar V(x\,x_c)$, with $\bar V$ the conjugate (antiparticle) representation ($\bar V=V$ for $U_q(\widehat{sl}_2)$)". eq-r-crossing as printed is then for self-conjugate V.
- **ch. 7:201 (minor).** For 2j+1 > p the zero-quantum-dimension summands are indecomposable but reducible (tilting modules), not spin-j irreducibles. Their q-dimension is still 0, since [m+1] + [2p−1−m] = 0. Suggested: "the summands with highest weight $j\ge(p-1)/2$ (irreducible for $2j+1=p$, where $[p]_q=0$, indecomposable but reducible for $2j+1>p$) have quantum dimension zero".
- **Checked and correct:** eq-spin-j (via [k+1][2j−k] − [k][2j+1−k] = [2j−2k]); the spin-½ matrix of R; S² = Ad K; the q-symmetric eigenvalue 1 of eq-six-vertex; the sine-Gordon substitution identity.

---

## C. Item 12': the 11 exercises

**Chapter 7**
1. **exr-coproduct-homomorphism: true and solvable.** The cross terms E K⁻¹⊗KF − K⁻¹E⊗FK cancel by KEK⁻¹ = q²E, leaving (K⊗K − K⁻¹⊗K⁻¹)/(q − q⁻¹).
2. **exr-antipode: true.** S(E)K + E = 0, S(K⁻¹)F + S(F) = 0, S(K)K = 1, S²(E) = KEK⁻¹, S²(F) = KFK⁻¹. Nit at ch. 7:210: "$m\circ(S\otimes1)\circ\Delta=\epsilon$" should read "$=\eta\circ\epsilon$ (i.e. $\epsilon(a)1$)".
3. **exr-spin-j: true.** It uses [k+1][2j−k] − [k][2j+1−k] = [2j−2k], and gives the classical spin j in an unnormalised basis.
4. **exr-six-vertex: true.**
   - ρ₀(1/z) = 1/ρ₀(z), and the projectors do not depend on z.
   - Ř(q⁻²) = P₁ exactly.
   - There is a pole at z = q² with residue (1 − q⁴)P₀.
   - Current (ch. 7:218): "whose residue is the projector onto spin 0?" Replace with "whose residue is proportional to the projector onto spin 0?"
5. **exr-sg-match: true** (sixvertex.py). The gauge is D(θ_i) = diag(e^{−λθ_i/2}, e^{λθ_i/2}) on each particle.
6. **exr-hermitian-spin1: true** (answer in B).

**Chapter 8**

7. **exr-six-vertex-derive: true as stated.** The e₁, f₁, k₁, e₀ equations alone leave exactly one free parameter and reproduce eq-six-vertex; no f₀ equation is needed (sympy).
8. **exr-tpg-check: true.** C(j) = 2j(j+1) and ⟨1⟩ (already verified).
9. **exr-tpg-an: true** (already verified; the hint works).
10. **exr-spin1: true.** The graph part was already verified. The fusion part is now verified numerically:
    - V₁(x) is the q-symmetric subspace of V_½(xq)⊗V_½(xq⁻¹), i.e. the image of Ř(q⁻²).
    - Fusing in both factors with eq-fusion gives eigenvalues 1 (×5), ⟨2⟩ (×3), ⟨2⟩⟨1⟩ (×1) at z = x/y.
    - The exercise is solvable but heavy without a hint. Suggest appending: "*Hint:* $V_1(x)$ is the $q$-symmetric part of $V_{1/2}(xq)\otimes V_{1/2}(xq^{-1})$, the image of $\check R(q^{-2})$; fuse first in the left and then in the right factor."
11. **exr-faces: true and solvable**, but the hint is incomplete: the b = 1 path also needs the 1⊗½ → ½ coefficients.
    - Current (ch. 8:137): "using the $q$-Clebsch–Gordan coefficients of $\frac12\otimes\frac12$". Replace with "… of $\frac12\otimes\frac12$ and of $1\otimes\frac12$".
    - Answer with path vectors (v₀v₀)v₁ + c(ΔF v₀v₀)v₀ for b = 1 and (q v₀v₁ − v₁v₀)v₀ for b = 0:
      - W₁₁ = (1 − q⁴z)/((1+q²)(z − q²))
      - W₀₀ = (z − q⁴)/((1+q²)(z − q²))
      - W₁₀ = q²(z − 1)/(z − q²)
      - W₀₁ = (z − 1)(1+q²+q⁴)/((1+q²)²(z − q²))
    - In the symmetric gauge the off-diagonal entries are q(z − 1)√[3]/([2](z − q²)). The eigenvalues are 1 and ρ₀, and W(1) = 1.

---

## D. Other chapter 8 checks (correct)
- eq-jimbo and the worked-example e₀ coproducts.
- The general structure of eq-spectral-decomposition.
- ⟨a⟩⟨−a⟩ = 1.
- The sl_{n+1} example.
- The pole and zero statements at ch. 8:86.
- eq-fusion (spectral shifts x_a = q, x_b = q⁻¹ for spin 1).
- The number of unknowns Σ mult² in sec-direct-method.

---

## E. Scripts (full) and output

### qaff.py
```python
"""Generic numerical tools: U_q(g^) relations and Jimbo's equations.

Conventions of ch. 7 (eq-coproduct, sec-uqghat):
  k_i e_j k_i^-1 = q_i^{a_ij} e_j, [e_i,f_j] = delta_ij (k_i-k_i^-1)/(q_i-q_i^-1),
  q_i = q^{alpha_i^2/2}, long roots alpha^2 = 2,
  Delta(k)=k(x)k, Delta(e)=e(x)k + 1(x)e, Delta(f)=f(x)1 + k^-1(x)f,
  homogeneous gradation: e_0 -> x e_0, f_0 -> f_0/x.
"""
import numpy as np
from math import comb


def qnum(n, q):
    return (q**n - q**(-n))/(q - 1/q)


def qfact(n, q):
    r = 1
    for j in range(1, n+1):
        r *= qnum(j, q)
    return r


def qbinom(n, s, q):
    return qfact(n, q)/(qfact(s, q)*qfact(n-s, q))


def check_relations(e, f, k, A, qi, tol=1e-9, verbose=False):
    """Return max residual over all U_q(g^) relations (incl. q-Serre with q_i)."""
    n = len(e); N = e[0].shape[0]; I = np.eye(N); res = {}
    for i in range(n):
        ki_inv = np.linalg.inv(k[i])
        for j in range(n):
            res[f'kk{i}{j}'] = np.abs(k[i]@k[j]-k[j]@k[i]).max()
            res[f'ke{i}{j}'] = np.abs(k[i]@e[j]@ki_inv - qi[i]**A[i][j]*e[j]).max()
            res[f'kf{i}{j}'] = np.abs(k[i]@f[j]@ki_inv - qi[i]**(-A[i][j])*f[j]).max()
            rhs = (k[i]-ki_inv)/(qi[i]-1/qi[i]) if i == j else 0*I
            res[f'ef{i}{j}'] = np.abs(e[i]@f[j]-f[j]@e[i]-rhs).max()
            if i != j:
                m = 1 - A[i][j]
                for X, nm in ((e, 'e'), (f, 'f')):
                    S = 0*I
                    for s in range(m+1):
                        S = S + (-1)**s*qbinom(m, s, qi[i])*np.linalg.matrix_power(X[i], m-s)@X[j]@np.linalg.matrix_power(X[i], s)
                    res[f'serre{nm}{i}{j}'] = np.abs(S).max()
    bad = {kk: v for kk, v in res.items() if v > tol}
    if verbose and bad:
        print('violated:', bad)
    return max(res.values()), bad


def grade(e, f, k, x, s):
    """Apply gradation s: e_i -> x^{s_i} e_i, f_i -> x^{-s_i} f_i."""
    return [x**s[i]*e[i] for i in range(len(e))], [x**(-s[i])*f[i] for i in range(len(e))], k


def coprod(gx, gy):
    (ex, fx, kx), (ey, fy, ky) = gx, gy
    N1 = ex[0].shape[0]; N2 = ey[0].shape[0]
    I1, I2 = np.eye(N1), np.eye(N2); out = []
    for i in range(len(ex)):
        out.append(np.kron(ex[i], kx[i]) + np.kron(I1, ey[i]))
        out.append(np.kron(fx[i], I2) + np.kron(np.linalg.inv(kx[i]), fy[i]))
        out.append(np.kron(kx[i], ky[i]))
    return out


def solve_jimbo(rep, x, y, s=None, rep2=None, return_dim=False):
    """Solve Rc Delta_{x,y}(a) = Delta_{y,x}(a) Rc on V(x)(x)V(y) (V=W)."""
    e, f, k = rep
    if s is None:
        s = [1] + [0]*(len(e)-1)
    gx = grade(e, f, k, x, s); gy = grade(e, f, k, y, s)
    A = coprod(gx, gy); B = coprod(gy, gx)
    M = A[0].shape[0]
    pairs = [(A[3*m], B[3*m]) for m in range(len(e))] + [(A[3*m+1], B[3*m+1]) for m in range(len(e))]
    kd = np.array([np.diag(A[3*i+2]) for i in range(len(e))]).T   # M x (r+1)
    null, sv, idx = intertwiner_space(pairs, kd, kd)
    null = null[:, 0]
    R = np.zeros((M, M), dtype=complex)
    for (u, v), n in idx.items():
        R[u, v] = null[n]
    if return_dim:
        return R, sv[-3:]
    return R, sv[-2:]


def P_flip(N):
    P = np.zeros((N*N, N*N))
    for i in range(N):
        for j in range(N):
            P[j*N+i, i*N+j] = 1
    return P


def intertwiner_space(pairs, kdL, kdR, nvec=1):
    """Solve R a = b R for all (a,b) in pairs, R restricted to entries (u,v) with
    equal k-weights kdL[u] == kdR[v] (k's diagonal). Sparse normal equations.
    Returns (lowest nvec null vectors as columns, singular values descending (last 3), index)."""
    import scipy.sparse as sps
    M = kdL.shape[0]
    keyL = [tuple(np.round(r, 10)) for r in kdL]; keyR = [tuple(np.round(r, 10)) for r in kdR]
    allowed = [(a, b) for a in range(M) for b in range(M) if keyL[a] == keyR[b]]
    idx = {ab: n for n, ab in enumerate(allowed)}
    G = None
    for a_, b_ in pairs:
        arow = [(np.nonzero(a_[c])[0], a_[c][np.nonzero(a_[c])[0]]) for c in range(M)]
        bcol = [(np.nonzero(b_[:, c])[0], b_[:, c][np.nonzero(b_[:, c])[0]]) for c in range(M)]
        ri, ci, va = [], [], []
        for (u, c), n in idx.items():
            vs, vals = arow[c]
            ri.extend(u*M+vs); ci.extend([n]*len(vs)); va.extend(vals)
        for (c, v), n in idx.items():
            us, vals = bcol[c]
            ri.extend(us*M+v); ci.extend([n]*len(us)); va.extend(-vals)
        L = sps.csr_matrix((va, (ri, ci)), shape=(M*M, len(allowed)))
        g = (L.conj().T @ L)
        G = g if G is None else G + g
    w, V = np.linalg.eigh(G.toarray())
    sv = np.sqrt(np.abs(w))
    return V[:, :nvec], sv[:3][::-1], idx
```

### tpg_checks.py
```python
"""eq-tpg-rule beyond sl2^: vector reps of U_q(a_2^(1)), U_q(c_2^(1)), U_q(b_2^(1)), U_q(c_3^(1)), U_q(b_3^(1)).

Builds e_i,f_i,k_i (i=0..n) explicitly in the book's conventions, checks all
defining relations (incl. q-Serre with q_i = q^{alpha_i^2/2}), solves Jimbo's
equations numerically in the homogeneous gradation at random q, x, y, and
compares the eigenvalues of Rc(z), normalized to 1 on v1(x)v1, with
eq-tpg-rule: rho_mu/rho_nu = <(C(nu)-C(mu))/4>, <a> = (1 - z q^{2a})/(z - q^{2a}).
"""
import numpy as np
from qaff import check_relations, solve_jimbo, qnum

rng = np.random.default_rng(1)


def Eij(N, i, j):
    M = np.zeros((N, N), dtype=complex); M[i-1, j-1] = 1; return M


def br(a, z, q):
    return (1 - z*q**(2*a))/(z - q**(2*a))


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


def match(ev, want, tol=1e-7):
    """ev: eigenvalues; want: list of (value, multiplicity). True if multisets agree."""
    ev = list(ev); ok = True
    for val, m in want:
        hits = [i for i, v in enumerate(ev) if abs(v - val) < tol*max(1, abs(val))]
        if len(hits) != m:
            ok = False
        for i in sorted(hits, reverse=True):
            ev.pop(i)
    return ok and not ev


def run(name, rep, A, qi, q, want_fn, N):
    err, bad = check_relations(*rep, A, qi, verbose=True)
    report(f'{name}: all U_q relations incl. q-Serre (max residual {err:.1e})', not bad)
    for trial in range(3):
        x = complex(rng.normal(), rng.normal()); y = complex(rng.normal(), rng.normal()); z = x/y
        R, sv = solve_jimbo(rep, x, y)
        R = R/R[0, 0]
        ev = np.linalg.eigvals(R)
        want = want_fn(z, q)
        ok = match(ev, want)
        report(f'{name}: solution space 1-dim (last sing. values {abs(sv[0]):.1e},{abs(sv[1]):.1e}) and eigenvalues = TPG prediction, trial {trial}',
               ok and abs(sv[0]) > 1e-4 and abs(sv[1]) < 1e-6)
        if not ok and trial == 0:
            print('   eigenvalues:', np.round(np.sort_complex(ev), 6))
            print('   predicted  :', [(np.round(v, 6), m) for v, m in want])
    return R


# random generic q = t^2 (t = q^{1/2})
t = complex(0.83, 0.41); q = t*t

# ---------------- a_2^(1) vector -----------------
N = 3
e = [Eij(N, 3, 1), Eij(N, 1, 2), Eij(N, 2, 3)]
f = [m.T.copy() for m in e]
h = [np.diag([-1, 0, 1]), np.diag([1, -1, 0]), np.diag([0, 1, -1])]
k = [np.diag(q**np.diag(hh).astype(complex)) for hh in h]
A = [[2, -1, -1], [-1, 2, -1], [-1, -1, 2]]
run('a_2^(1) vector', (e, f, k), A, [q, q, q], q,
    lambda z, q: [(1, 6), (br(1, z, q), 3)], N)

# ---------------- c_2^(1) vector (C^4) -----------------
# weights e1,e2,-e2,-e1 with eps_i.eps_j = delta/2; alpha_1=eps1-eps2 short, alpha_2=2eps2 long, alpha_0=-2eps1 long
N = 4
e = [Eij(N, 4, 1), Eij(N, 1, 2) + Eij(N, 3, 4), Eij(N, 2, 3)]
f = [m.T.copy() for m in e]
h = [np.diag([-1, 0, 0, 1]), np.diag([1, -1, 1, -1]), np.diag([0, 1, -1, 0])]
d = [1, 0.5, 1]                       # alpha_i^2/2
qi = [t**(2*di) for di in d]          # q_i = q^{alpha_i^2/2}
k = [np.diag(qi[i]**np.diag(h[i]).astype(complex)) for i in range(3)]
A = [[2, -1, 0], [-2, 2, -2], [0, -1, 2]]
# Casimirs (long roots length^2 2): C(2l1)=2n+2=6, C(l2)=2n=4, C(0)=0  (n=2)
run('c_2^(1) vector', (e, f, k), A, qi, q,
    lambda z, q: [(1, 10), (br(0.5, z, q), 5), (br(1.5, z, q), 1)], N)

# ---------------- b_2^(1) vector (C^5) -----------------
# weights e1,e2,0,-e2,-e1 orthonormal; alpha_1=e1-e2 long, alpha_2=e2 short, alpha_0=-e1-e2
N = 5
q2 = t                                  # q_2 = q^{1/2}
e = [Eij(N, 5, 2) + Eij(N, 4, 1), Eij(N, 1, 2) + Eij(N, 4, 5), qnum(2, q2)*Eij(N, 2, 3) + Eij(N, 3, 4)]
f = [Eij(N, 2, 5) + Eij(N, 1, 4), Eij(N, 2, 1) + Eij(N, 5, 4), Eij(N, 3, 2) + qnum(2, q2)*Eij(N, 4, 3)]
h = [np.diag([-1, -1, 0, 1, 1]), np.diag([1, -1, 0, 1, -1]), np.diag([0, 2, 0, -2, 0])]
qi = [q, q, q2]
k = [np.diag(qi[i]**np.diag(h[i]).astype(complex)) for i in range(3)]
A = [[2, 0, -1], [0, 2, -1], [-2, -2, 2]]
# C(2l1)=4n+2=10, C(l2)=4n-2=6, C(0)=0 (n=2): <1>, then <(2n-1)/2>=<3/2>
run('b_2^(1) vector', (e, f, k), A, qi, q,
    lambda z, q: [(1, 14), (br(1, z, q), 10), (br(1, z, q)*br(1.5, z, q), 1)], N)

# ---------------- c_3^(1) vector (C^6) -----------------
# basis 1..6 = eps1,eps2,eps3,-eps3,-eps2,-eps1
N = 6
e = [Eij(N, 6, 1), Eij(N, 1, 2) + Eij(N, 5, 6), Eij(N, 2, 3) + Eij(N, 4, 5), Eij(N, 3, 4)]
f = [m.T.copy() for m in e]
h = [np.diag([-1, 0, 0, 0, 0, 1]), np.diag([1, -1, 0, 0, 1, -1]), np.diag([0, 1, -1, 1, -1, 0]), np.diag([0, 0, 1, -1, 0, 0])]
qi = [q, t, t, q]
k = [np.diag(qi[i]**np.diag(h[i]).astype(complex)) for i in range(4)]
A = [[2, -1, 0, 0], [-2, 2, -1, 0], [0, -1, 2, -2], [0, 0, -1, 2]]
# C(2l1)=2n+2=8, C(l2)=2n=6, C(0)=0: rho_l2=<1/2>, rho_0=<(n+1)/2>=<2>
run('c_3^(1) vector', (e, f, k), A, qi, q,
    lambda z, q: [(1, 21), (br(0.5, z, q), 14), (br(2, z, q), 1)], N)

# ---------------- b_3^(1) vector (C^7) -----------------
# basis 1..7 = e1,e2,e3,0,-e3,-e2,-e1 ; alpha_1=e1-e2, alpha_2=e2-e3 long, alpha_3=e3 short, alpha_0=-e1-e2
N = 7
e = [Eij(N, 7, 2) + Eij(N, 6, 1), Eij(N, 1, 2) + Eij(N, 6, 7), Eij(N, 2, 3) + Eij(N, 5, 6), qnum(2, t)*Eij(N, 3, 4) + Eij(N, 4, 5)]
f = [Eij(N, 2, 7) + Eij(N, 1, 6), Eij(N, 2, 1) + Eij(N, 7, 6), Eij(N, 3, 2) + Eij(N, 6, 5), Eij(N, 4, 3) + qnum(2, t)*Eij(N, 5, 4)]
h = [np.diag([-1, -1, 0, 0, 0, 1, 1]), np.diag([1, -1, 0, 0, 0, 1, -1]), np.diag([0, 1, -1, 0, 1, -1, 0]), np.diag([0, 0, 2, 0, -2, 0, 0])]
qi = [q, q, q, t]
k = [np.diag(qi[i]**np.diag(h[i]).astype(complex)) for i in range(4)]
A = [[2, 0, -1, 0], [0, 2, -1, 0], [-1, -1, 2, -1], [0, 0, -2, 2]]
# C(2l1)=4n+2=14, C(l2)=4n-2=10: <1>, <(2n-1)/2>=<5/2>
run('b_3^(1) vector', (e, f, k), A, qi, q,
    lambda z, q: [(1, 27), (br(1, z, q), 21), (br(1, z, q)*br(2.5, z, q), 1)], N)

# ---------------- negative controls: wrong normalisations must fail -----------------
from qaff import solve_jimbo as sj
N = 4
e = [Eij(N, 4, 1), Eij(N, 1, 2) + Eij(N, 3, 4), Eij(N, 2, 3)]
f = [m.T.copy() for m in e]
h = [np.diag([-1, 0, 0, 1]), np.diag([1, -1, 1, -1]), np.diag([0, 1, -1, 0])]
qi = [q, t, q]
k = [np.diag(qi[i]**np.diag(h[i]).astype(complex)) for i in range(3)]
x, y = 0.3+0.9j, -1.2+0.4j; z = x/y
R, _ = sj((e, f, k), x, y); ev = np.linalg.eigvals(R/R[0, 0])
for label, want in (('C with long roots length^2 4 (<1>,<3>)', [(1, 10), (br(1, z, q), 5), (br(3, z, q), 1)]),
                    ('rule with z -> 1/z', [(1, 10), (br(0.5, 1/z, q), 5), (br(1.5, 1/z, q), 1)]),
                    ('book conventions (<1/2>,<3/2>)', [(1, 10), (br(0.5, z, q), 5), (br(1.5, z, q), 1)])):
    print(f'   control c_2^(1): eigenvalues match {label}: {match(ev, want)}')
```
Output (2 s):
```
PASS a_2^(1) vector: all U_q relations incl. q-Serre (max residual 2.5e-16)
PASS a_2^(1) vector: solution space 1-dim (last sing. values 1.7e+00,2.6e-08) and eigenvalues = TPG prediction, trial 0
PASS a_2^(1) vector: ... trial 1
PASS a_2^(1) vector: ... trial 2
PASS c_2^(1) vector: all U_q relations incl. q-Serre (max residual 4.0e-16)
PASS c_2^(1) vector: solution space 1-dim (last sing. values 4.5e-01,7.0e-08) and eigenvalues = TPG prediction, trial 0 (also trials 1, 2)
PASS b_2^(1) vector: all U_q relations incl. q-Serre (max residual 4.5e-16)
PASS b_2^(1) vector: solution space 1-dim (last sing. values 1.1e+00,5.6e-08) and eigenvalues = TPG prediction, trial 0 (also trials 1, 2)
PASS c_3^(1) vector: all U_q relations incl. q-Serre (max residual 4.0e-16)
PASS c_3^(1) vector: solution space 1-dim (last sing. values 5.5e-01,5.0e-08) and eigenvalues = TPG prediction, trial 0 (also trials 1, 2)
PASS b_3^(1) vector: all U_q relations incl. q-Serre (max residual 4.5e-16)
PASS b_3^(1) vector: solution space 1-dim (last sing. values 8.3e-01,4.6e-08) and eigenvalues = TPG prediction, trial 0 (also trials 1, 2)
   control c_2^(1): eigenvalues match C with long roots length^2 4 (<1>,<3>): False
   control c_2^(1): eigenvalues match rule with z -> 1/z: False
   control c_2^(1): eigenvalues match book conventions (<1/2>,<3/2>): True
```

### spinor_twisted.py
```python
"""Spinor of U_q(d_{l+1}^(2)) in the conventions of ch. 7 (sec-tpg 'When it fails').

Finite subalgebra (remove node 0, both end nodes short): U_q(b_l).
1. builds e_i,f_i,k_i (i=0..l) on (C^2)^{(x)l}, checks all relations incl. q-Serre
   with q_i = q^{alpha_i^2/2} (long q, short q^{1/2});
2. dimension of the commutant of U_q(b_l) on S(x)S  (= l+1 iff multiplicity-free);
3. solves Jimbo's equations (homogeneous gradation, book coproduct): dimension of
   the solution space, eigenvalues of Rc(z) vs the twisted TPG of DGZ 1996 sec. 4.3
   in book units: rho_{Lambda^{l-j}} = prod_{i=1}^j <i/2>_{s_i},
   <a>_- = (1 - z q^{2a})/(z - q^{2a}),  <a>_+ = (1 + z q^{2a})/(z + q^{2a}),
   and vs the untwisted rule eq-tpg-rule (all signs -).
"""
import numpy as np
from functools import reduce
from qaff import check_relations, solve_jimbo, coprod


def commutant_sv(pairs, kd):
    import scipy.sparse as sps
    M = kd.shape[0]
    key = [tuple(np.round(r, 10)) for r in kd]
    allowed = [(a, b) for a in range(M) for b in range(M) if key[a] == key[b]]
    idx = {ab: n for n, ab in enumerate(allowed)}
    G = 0
    for a_, b_ in pairs:
        ri, ci, va = [], [], []
        for (u, c), n in idx.items():
            vs = np.nonzero(a_[c])[0]; ri.extend(u*M+vs); ci.extend([n]*len(vs)); va.extend(a_[c][vs])
        for (c, v), n in idx.items():
            us = np.nonzero(b_[:, c])[0]; ri.extend(us*M+v); ci.extend([n]*len(us)); va.extend(-b_[:, c][us])
        L = sps.csr_matrix((va, (ri, ci)), shape=(M*M, len(allowed)))
        G = G + (L.conj().T @ L)
    w = np.linalg.eigvalsh(G.toarray())
    return np.sqrt(np.abs(w))

rng = np.random.default_rng(7)
sp_ = np.array([[0, 1], [0, 0]], dtype=complex); sm_ = sp_.T.copy(); I2 = np.eye(2)


def site(op, i, l):
    return reduce(np.kron, [op if j == i else I2 for j in range(l)])


def spinor(l, t):
    q = t*t
    sz = [site(np.diag([0.5, -0.5]), i, l) for i in range(l)]   # s_i
    e = [site(sm_, 0, l)] + [site(sp_, i, l)@site(sm_, i+1, l) for i in range(l-1)] + [site(sp_, l-1, l)]
    f = [m.conj().T.copy() for m in e]
    # h_0 = -2 s_1 (alpha_0 = -eps_1 short), h_i = s_i - s_{i+1} (long), h_l = 2 s_l (short)
    h = [-2*sz[0]] + [sz[i]-sz[i+1] for i in range(l-1)] + [2*sz[l-1]]
    qi = [t] + [q]*(l-1) + [t]
    k = [np.diag(qi[i]**np.diag(h[i])) for i in range(l+1)]
    A = np.zeros((l+1, l+1), int)
    for i in range(l+1):
        A[i, i] = 2
    A[0, 1], A[1, 0] = -2, -1
    for i in range(1, l-1):
        A[i, i+1] = A[i+1, i] = -1
    A[l-1, l], A[l, l-1] = -1, -2
    if l == 1:
        raise ValueError
    return (e, f, k), A.tolist(), qi


def ang(a, z, q, s):
    return (1 + s*z*q**(2*a))/(z + s*q**(2*a))


def comb(n, k):
    from math import comb as c
    return c(n, k)


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


def match(ev, want, tol=1e-6):
    ev = list(ev); ok = True
    for val, m in want:
        hits = [i for i, v in enumerate(ev) if abs(v - val) < tol*max(1, abs(val))]
        if len(hits) != m:
            ok = False
        for i in sorted(hits, reverse=True):
            ev.pop(i)
    return ok and not ev


t = complex(0.83, 0.41); q = t*t
for l in (2, 3, 4):
    rep, A, qi = spinor(l, t)
    err, bad = check_relations(*rep, A, qi, verbose=True)
    report(f'l={l}: spinor is a rep of U_q(d_{l+1}^(2)) (max residual {err:.1e})', not bad)
    # commutant of U_q(b_l) (nodes 1..l) on S(x)S
    e, f, k = rep
    sub = ([e[i] for i in range(1, l+1)], [f[i] for i in range(1, l+1)], [k[i] for i in range(1, l+1)])
    D = coprod(sub, sub)
    kd = np.array([np.diag(D[3*i+2]) for i in range(l)]).T
    pairs = [(D[3*m], D[3*m]) for m in range(l)] + [(D[3*m+1], D[3*m+1]) for m in range(l)]
    Gsv = commutant_sv(pairs, kd)
    dimcomm = int((Gsv < 1e-6).sum())
    report(f'l={l}: commutant of U_q(b_l) on S(x)S has dim {dimcomm} = l+1 (multiplicity-free); smallest sv {np.round(Gsv[:l+3],8)}', dimcomm == l+1)
    x = complex(rng.normal(), rng.normal()); y = complex(rng.normal(), rng.normal()); z = x/y
    R, sv2 = solve_jimbo(rep, x, y)
    report(f'l={l}: Jimbo solution space 1-dim (sing. values {abs(sv2[0]):.1e}, {abs(sv2[1]):.1e})', abs(sv2[0]) > 1e-4 and abs(sv2[1]) < 1e-6)
    R = R/R[0, 0]
    ev = np.linalg.eigvals(R)
    n = 2*l+1
    # twisted TPG (DGZ 1996 sec 4.3), translated to book units, s_i = (-1)^i
    for label, signs in (('twisted TPG, s_i=(-1)^i', [(-1)**i for i in range(1, l+1)]),
                         ('untwisted rule, all s_i=-1', [-1]*l),
                         ('twisted TPG with s_i=-(-1)^i', [-(-1)**i for i in range(1, l+1)])):
        want = []; rho = 1
        for j in range(0, l+1):
            if j > 0:
                rho = rho*ang(j/2, z, q, signs[j-1])
            want.append((rho, comb(n, l-j)))
        ok = match(ev, want)
        print(f'      l={l}: eigenvalues match {label}: {ok}')
```
Output (21 s):
```
PASS l=2: spinor is a rep of U_q(d_3^(2)) (max residual 4.0e-16)
PASS l=2: commutant of U_q(b_l) on S(x)S has dim 3 = l+1 (multiplicity-free); smallest sv [7e-08 3e-08 6e-08 9.04e-01 1.06e+00]
PASS l=2: Jimbo solution space 1-dim (sing. values 1.1e+00, 9.6e-08)
      l=2: eigenvalues match twisted TPG, s_i=(-1)^i: True
      l=2: eigenvalues match untwisted rule, all s_i=-1: False
      l=2: eigenvalues match twisted TPG with s_i=-(-1)^i: False
PASS l=3: spinor is a rep of U_q(d_4^(2)) (max residual 4.0e-16)
PASS l=3: commutant of U_q(b_l) on S(x)S has dim 4 = l+1 (multiplicity-free); smallest sv [7e-08 5e-08 5e-08 3e-08 5.67e-01 6.85e-01]
PASS l=3: Jimbo solution space 1-dim (sing. values 3.2e-01, 9.8e-08)
      l=3: eigenvalues match twisted TPG, s_i=(-1)^i: True
      l=3: eigenvalues match untwisted rule, all s_i=-1: False
      l=3: eigenvalues match twisted TPG with s_i=-(-1)^i: False
PASS l=4: spinor is a rep of U_q(d_5^(2)) (max residual 4.0e-16)
PASS l=4: commutant of U_q(b_l) on S(x)S has dim 5 = l+1 (multiplicity-free); smallest sv [1e-07 1.1e-07 1.6e-07 1.6e-07 1.9e-07 3.88e-01 5.08e-01]
PASS l=4: Jimbo solution space 1-dim (sing. values 4.1e-01, 1.2e-07)
      l=4: eigenvalues match twisted TPG, s_i=(-1)^i: True
      l=4: eigenvalues match untwisted rule, all s_i=-1: False
      l=4: eigenvalues match twisted TPG with s_i=-(-1)^i: False
```

### ch7_checks.py
```python
"""Chapter 7 checks: central element, quasitriangularity, Hermitian form on spin 1.

1. sec-uqghat: is prod_i k_i^{n_i} (Kac labels) or prod_i k_i^{n_i^vee} (dual labels,
   n_i^vee = n_i alpha_i^2/2, tbl-untwisted) central?  Exact exponent bookkeeping
   plus explicit matrices for the c_2^(1) and b_2^(1) vector representations.
2. sec-universal-r: for R = q^{H(x)H/2} sum_n q^{n(n-1)/2}(q-q^-1)^n/[n]! E^n(x)F^n and
   the coproduct eq-coproduct, check on spin 1/2 and spin 1:
   (Delta(x)1)R = R13 R23, (1(x)Delta)R = R13 R12, YBE R12R13R23 = R23R13R12,
   (S(x)1)R = R^-1, (S(x)S)R = R.
3. exr-hermitian-spin1: invariant form <a v, w> = <v, a^* w>, E^*=F, K^*=K^-1, q=e^{i gamma}.
"""
import numpy as np
from fractions import Fraction as Fr
from qaff import qnum, qfact


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


# ---------- 1. central element ----------
# conjugation of e_j by prod_i k_i^{m_i} multiplies it by q^{sum_i m_i d_i a_ij}, d_i = alpha_i^2/2
cases = {
    'c_2^(1)': ([[2, -1, 0], [-2, 2, -2], [0, -1, 2]], [Fr(1), Fr(1, 2), Fr(1)], [1, 2, 1]),
    'b_2^(1)': ([[2, 0, -1], [0, 2, -1], [-2, -2, 2]], [Fr(1), Fr(1), Fr(1, 2)], [1, 1, 2]),
    'g_2^(1)': ([[2, -1, 0], [-1, 2, -1], [0, -3, 2]], [Fr(1), Fr(1), Fr(1, 3)], [1, 2, 3]),
    'd_3^(2)': ([[2, -2, 0], [-1, 2, -1], [0, -2, 2]], [Fr(1, 2), Fr(1), Fr(1, 2)], [1, 1, 1]),
    'a_2^(2)': ([[2, -4], [-1, 2]], [Fr(1, 4), Fr(1)], [2, 1]),
}
for name, (A, d, n) in cases.items():
    r = len(n)
    nv = [n[i]*d[i] for i in range(r)]
    # check n is the null vector: sum_i n_i alpha_i = 0  <=> sum_i n_i d_i a_ij = 0 for all j
    shift = lambda m: [sum(m[i]*d[i]*A[i][j] for i in range(r)) for j in range(r)]
    print(f'{name}: n={n}, n^vee={[str(x) for x in nv]}, q-exponents for prod k^n: {[str(s) for s in shift(n)]}, '
          f'for prod k^(n^vee): {[str(s) for s in shift(nv)]}')
    report(f'{name}: prod k_i^(n_i) central', all(s == 0 for s in shift(n)))
    report(f'{name}: prod k_i^(n_i^vee) NOT central (book text claims it is)', not all(s == 0 for s in shift(nv)) or n == nv)

# explicit c_2^(1) vector rep: compute prod k^n and prod k^(n^vee)
t = 0.83+0.41j; q = t*t
h = [np.diag([-1, 0, 0, 1]), np.diag([1, -1, 1, -1]), np.diag([0, 1, -1, 0])]
qi = [q, t, q]
k = [np.diag(qi[i]**np.diag(h[i]).astype(complex)) for i in range(3)]
Kn = k[0]@k[1]@k[1]@k[2]; Knv = k[0]@k[1]@k[2]
report('c_2^(1) vector rep: k0 k1^2 k2 = 1', np.allclose(Kn, np.eye(4)))
report('c_2^(1) vector rep: k0 k1 k2 is not a scalar', not np.allclose(Knv, Knv[0, 0]*np.eye(4)))
print('   k0 k1 k2 =', np.round(np.diag(Knv), 4), ' (q^{1/2} =', np.round(t, 4), ')')

# ---------- 2. quasitriangularity ----------


def rep(twoj, q):
    dd = twoj+1
    E, F = np.zeros((dd, dd), complex), np.zeros((dd, dd), complex)
    for kk in range(dd-1):
        F[kk+1, kk] = 1
    for kk in range(1, dd):
        E[kk-1, kk] = qnum(kk, q)*qnum(twoj+1-kk, q)
    H = np.array([twoj-2*kk for kk in range(dd)], float)
    return E, F, np.diag(q**H), H


def kr(*ms):
    out = ms[0]
    for m in ms[1:]:
        out = np.kron(out, m)
    return out


def qHH(Ha, Hb, q, sign=1):
    return np.diag([q**(sign*a*b/2) for a in Ha for b in Hb])


def coef(n, q):
    return q**(n*(n-1)/2)*(q-1/q)**n/qfact(n, q)


for twoj in (1, 2):
    E, F, K, H = rep(twoj, q); d = twoj+1; I = np.eye(d); Ki = np.linalg.inv(K)
    mp = np.linalg.matrix_power
    R = qHH(H, H, q) @ sum(coef(n, q)*np.kron(mp(E, n), mp(F, n)) for n in range(d))
    # embeddings into V(x)V(x)V
    D = np.diag
    H1, H2, H3 = [np.kron(np.kron(a, b), c) for a, b, c in ((D(H), I, I), (I, D(H), I), (I, I, D(H)))]
    expo = lambda X, Y: np.diag(q**(np.diag(X)*np.diag(Y)/2))
    E1, E2 = kr(E, I, I), kr(I, E, I); F2, F3 = kr(I, F, I), kr(I, I, F)
    K2 = kr(I, K, I)
    R12 = kr(R, I)
    R23 = kr(I, R)
    # R13 = sum over terms: q^{H1 H3/2} sum c_n E1^n F3^n
    E1n = lambda n: mp(E1, n)
    R13 = expo(H1, H3) @ sum(coef(n, q)*E1n(n)@mp(F3, n) for n in range(d))
    # (Delta(x)1)R = q^{(H1+H2)H3/2} sum c_n (E1 K2 + E2)^n F3^n
    DE = kr(E, K, I) + kr(I, E, I)
    lhs1 = expo(H1+H2, H3) @ sum(coef(n, q)*mp(DE, n)@mp(F3, n) for n in range(2*d))
    report(f'spin {twoj}/2: (Delta(x)1)R = R13 R23', np.allclose(lhs1, R13@R23))
    report(f'spin {twoj}/2: (Delta(x)1)R != R23 R13 (order matters)', not np.allclose(lhs1, R23@R13))
    # (1(x)Delta)R = q^{H1(H2+H3)/2} sum c_n E1^n (F2 + K2^-1 F3)^n
    DF = kr(I, F, I) + kr(I, Ki, F)
    lhs2 = expo(H1, H2+H3) @ sum(coef(n, q)*E1n(n)@mp(DF, n) for n in range(2*d))
    report(f'spin {twoj}/2: (1(x)Delta)R = R13 R12', np.allclose(lhs2, R13@R12))
    report(f'spin {twoj}/2: (1(x)Delta)R != R12 R13 (order matters)', not np.allclose(lhs2, R12@R13))
    report(f'spin {twoj}/2: YBE R12 R13 R23 = R23 R13 R12', np.allclose(R12@R13@R23, R23@R13@R12))
    SE = -E@Ki; SF = -K@F
    S1R = sum(coef(n, q)*np.kron(mp(SE, n), I)@qHH(H, H, q, -1)@np.kron(I, mp(F, n)) for n in range(d))
    report(f'spin {twoj}/2: (S(x)1)R = R^-1', np.allclose(S1R, np.linalg.inv(R)))
    SSR = sum(coef(n, q)*np.kron(mp(SE, n), mp(SF, n)) for n in range(d)) @ qHH(H, H, q)
    report(f'spin {twoj}/2: (S(x)S)R = R', np.allclose(SSR, R))

# evaluation rep eq-evaluation-rep: q-Serre (a_01=a_10=-2) for spin 1/2..2
for twoj in (1, 2, 3, 4):
    E, F, K, H = rep(twoj, q); mp = np.linalg.matrix_power; xx = 0.37-0.8j
    e = [xx*F, E]; f = [E/xx, F]
    ok = True
    for (i, j) in ((0, 1), (1, 0)):
        for X in (e, f):
            S = sum((-1)**s_*qfact(3, q)/(qfact(s_, q)*qfact(3-s_, q))*mp(X[i], 3-s_)@X[j]@mp(X[i], s_) for s_ in range(4))
            ok &= np.allclose(S, 0)
    report(f'spin {twoj}/2 evaluation rep satisfies both q-Serre relations', ok)

# ---------- 3. exr-hermitian-spin1 ----------
for gamma in (0.4, 1.2, 2.0, 2.9):
    qq = np.exp(1j*gamma)
    E, F, K, H = rep(2, qq)
    # Gram matrix G with <v,w> = v^dag G w ; invariance <a v,w> = <v,a^* w>  <=> a^dag G = G a^*
    # solve linear equations for Hermitian G
    ops = [(E, F), (F, E), (K, np.linalg.inv(K))]
    rows = []
    for a, astar in ops:
        # a^dag G - G a^* = 0, vectorised (row-major): (a^dag (x) I - I (x) astar^T) vec(G)
        rows.append(np.kron(a.conj().T, np.eye(3)) - np.kron(np.eye(3), astar.T))
    L = np.vstack(rows); u, s, vh = np.linalg.svd(L)
    G = vh[-1].conj().reshape(3, 3); G = G/G[0, 0]
    ok_unique = s[-2] > 1e-6 and s[-1] < 1e-10
    want = np.diag([1, 2*np.cos(gamma), 4*np.cos(gamma)**2])
    ev = np.linalg.eigvalsh((G+G.conj().T)/2)
    sig = (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()))
    report(f'spin-1 form at gamma={gamma}: unique, = diag(1,[2],[2]^2), signature {sig}', ok_unique and np.allclose(G, want))
```
Output:
```
c_2^(1): n=[1, 2, 1], n^vee=['1', '1', '1'], q-exponents for prod k^n: ['0', '0', '0'], for prod k^(n^vee): ['1', '-1', '1']
PASS c_2^(1): prod k_i^(n_i) central
PASS c_2^(1): prod k_i^(n_i^vee) NOT central (book text claims it is)
b_2^(1): n=[1, 1, 2], n^vee=['1', '1', '1'], ... for prod k^(n^vee): ['1', '1', '-1']   (PASS, PASS)
g_2^(1): n=[1, 2, 3], n^vee=['1', '2', '1'], ... for prod k^(n^vee): ['0', '2', '-4/3']  (PASS, PASS)
d_3^(2): n=[1, 1, 1], n^vee=['1/2', '1', '1/2'], ... for prod k^(n^vee): ['-1/2', '1', '-1/2']  (PASS, PASS)
a_2^(2): n=[2, 1], n^vee=['1/2', '1'], ... for prod k^(n^vee): ['-3/4', '3/2']  (PASS, PASS)
PASS c_2^(1) vector rep: k0 k1^2 k2 = 1
PASS c_2^(1) vector rep: k0 k1 k2 is not a scalar
   k0 k1 k2 = [0.9685-0.4784j 0.83+0.41j 0.9685-0.4784j 0.83+0.41j]  (q^{1/2} = (0.83+0.41j))
PASS spin 1/2: (Delta(x)1)R = R13 R23
PASS spin 1/2: (Delta(x)1)R != R23 R13 (order matters)
PASS spin 1/2: (1(x)Delta)R = R13 R12
PASS spin 1/2: (1(x)Delta)R != R12 R13 (order matters)
PASS spin 1/2: YBE R12 R13 R23 = R23 R13 R12
PASS spin 1/2: (S(x)1)R = R^-1
PASS spin 1/2: (S(x)S)R = R
(same 7 PASS lines for spin 2/2)
PASS spin 1/2 .. 4/2 evaluation rep satisfies both q-Serre relations (4 lines)
PASS spin-1 form at gamma=0.4: unique, = diag(1,[2],[2]^2), signature (3, 0)
PASS spin-1 form at gamma=1.2: unique, = diag(1,[2],[2]^2), signature (3, 0)
PASS spin-1 form at gamma=2.0: unique, = diag(1,[2],[2]^2), signature (2, 1)
PASS spin-1 form at gamma=2.9: unique, = diag(1,[2],[2]^2), signature (2, 1)
```

### ch8_checks.py
```python
"""Chapter 8 exercise checks.

1. exr-six-vertex-derive: the e_1, f_1, (k_1) and e_0 equations alone already fix
   the six-vertex Rc (no f_0 equation needed), as the exercise asserts.
2. exr-spin1 (fusion part): fuse the six-vertex Rc twice, with the spin-1
   submodule of V(x q) (x) V(x q^-1), and compare eigenvalues with rho_1=<2>, rho_0=<2><1>.
3. exr-faces: face weights for a=d=1/2 via highest-weight vectors in V^{(x)3}.
"""
import sympy as sp
import numpy as np

q, x, y, z = sp.symbols('q x y z', nonzero=True)
E = sp.Matrix([[0, 1], [0, 0]]); F = sp.Matrix([[0, 0], [1, 0]]); K = sp.diag(q, 1/q)
I2 = sp.eye(2); kp = sp.kronecker_product


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


def Rc(u):
    return sp.Matrix([[1, 0, 0, 0],
                      [0, u*(q**2-1)/(q**2-u), q*(u-1)/(u-q**2), 0],
                      [0, q*(u-1)/(u-q**2), (q**2-1)/(q**2-u), 0],
                      [0, 0, 0, 1]])


# ---- 1. which equations are needed
syms = sp.symbols('r0:16'); R = sp.Matrix(4, 4, syms)
g = lambda w: {'e1': E, 'f1': F, 'k1': K, 'e0': w*F, 'f0': E/w, 'k0': K.inv()}
D = lambda gx, gy, key: {'e': kp(gx['e'+key[1]], gx['k'+key[1]]) + kp(I2, gy['e'+key[1]]),
                         'f': kp(gx['f'+key[1]], I2) + kp(gx['k'+key[1]].inv(), gy['f'+key[1]]),
                         'k': kp(gx['k'+key[1]], gy['k'+key[1]])}[key[0]]
eqs = []
for key in ('k1', 'e1', 'f1', 'e0'):
    eqs += list(R*D(g(x), g(y), key) - D(g(y), g(x), key)*R)
sol = sp.solve(eqs, syms, dict=True)
free = set().union(*[R.subs(s).free_symbols for s in sol]) & set(syms)
Rs = R.subs(sol[0])
Rs = sp.simplify((Rs/Rs[0, 0]).subs(y, 1))
report(f'e1,f1,k1,e0 equations: one free parameter ({len(free)}) and solution = eq-six-vertex',
       len(sol) == 1 and len(free) == 1 and sp.simplify(Rs - Rc(x)) == sp.zeros(4))

# ---- 2. fusion to spin 1
qn = complex(0.83, 0.41)**2
Rn = sp.lambdify((z,), Rc(z).subs(q, qn), 'numpy')
Rh = lambda w: np.array(Rn(w), dtype=complex)
I = np.eye(2)
# spin-1 subspace of V(x q)(x)V(x q^-1): image of Rc(q^-2): V(x q^-1)(x)V(x q) -> V(x q)(x)V(x q^-1)
P1 = Rh(qn**-2)
U, s, Vh = np.linalg.svd(P1); B = U[:, :3]          # 4x3 basis of the image
report('Rc(q^-2) has rank 3 (spin-1 image)', s[2] > 1e-8 and s[3] < 1e-10)
xx, yy = complex(0.7, -0.4), complex(-0.3, 1.1)
# step 1: Rc_{1,1/2}(x/w): V1(x)(x)V(w) -> V(w)(x)V1(x), V1(x) in V(xq)(x)V(x/q)
def R1h(x_, w_):
    M = np.kron(Rh(x_*qn/w_), I) @ np.kron(I, Rh(x_/qn/w_))      # on V_a V_b V_d -> V_d V_a V_b
    return M
# step 2: Rc_{1,1}(x/y): V1(x)(x)V1(y), V1(y) in V(yq)(x)V(y/q)
#   V_c (x) V_a' (x) V_b'  --(R_{c a'} (x) 1)-->  V_a' V_c V_b'  --(1 (x) R_{c b'})--> V_a' V_b' V_c
# build explicit operators on (C^2)^{(x)4}: factors (c1 c2)(a')(b')
def op_on(M8, pos):   # M8 acts on 3 consecutive qubits starting at pos, in a 4-qubit space
    left = np.eye(2**pos); right = np.eye(2**(4-pos-3))
    return np.kron(np.kron(left, M8), right)
T = op_on(R1h(xx, yy/qn), 1) @ op_on(R1h(xx, yy*qn), 0)
# T maps V_c(=qubits 0,1) (x) V_a'(2) (x) V_b'(3)  ->  V_a'(0) V_b'(1) V_c(2,3)
emb = np.kron(B, B)                                         # 16x9 : V1(x) (x) V1(y)
img = T @ emb
coef, res, *_ = np.linalg.lstsq(emb, img, rcond=None)
report('fused operator maps V1(x)(x)V1(y) into V1(y)(x)V1(x)', np.allclose(emb @ coef, img, atol=1e-9))
ev = np.linalg.eigvals(coef); zz = xx/yy
# normalise on the top component: eigenvalue of multiplicity 5
vals, counts = np.unique(np.round(ev, 8), return_counts=True)
top = vals[np.argmax(counts)]
ev = ev/top
br = lambda a: (1 - zz*qn**(2*a))/(zz - qn**(2*a))
want = [(1, 5), (br(2), 3), (br(2)*br(1), 1)]
ok = all(sum(abs(e - w) < 1e-7 for e in ev) == m for w, m in want)
report('fused spin-1 R-matrix: eigenvalues 1 (x5), <2> (x3), <2><1> (x1) at z = x/y', ok)

# ---- 3. face weights a=d=1/2
Rz = Rc(z)
D2E = kp(kp(E, K), K) + kp(kp(I2, E), K) + kp(kp(I2, I2), E)
v = lambda *bits: kp(kp(sp.eye(2)[:, bits[0]], sp.eye(2)[:, bits[1]]), sp.eye(2)[:, bits[2]])
DF = kp(F, I2) + kp(K.inv(), F)
# b = 1: spin-1 hw v0v0 in factors 12, lowered once; combine with factor 3 into weight-1/2 hw vector
w1 = kp(sp.Matrix([1, 0, 0, 0]), sp.Matrix([0, 1]))                 # (v0 v0) v1
w2 = kp(DF*sp.Matrix([1, 0, 0, 0]), sp.Matrix([1, 0]))             # (Delta F v0v0) v0
c = sp.symbols('c')
hw1 = w1 + c*w2
csol = [sp.solve(list(D2E*hw1), c, dict=True)[0][c]]
hw1 = sp.simplify(hw1.subs(c, csol[0]))
hw0 = kp(q*sp.Matrix([0, 1, 0, 0]) - sp.Matrix([0, 0, 1, 0]), sp.Matrix([1, 0]))   # singlet(12) (x) v0
report('path vectors are highest-weight', sp.simplify(D2E*hw1) == sp.zeros(8, 1) and sp.simplify(D2E*hw0) == sp.zeros(8, 1))
Op = kp(I2, Rz)
Bm = sp.Matrix.hstack(hw1, hw0)
imgs = Op*Bm
W = sp.simplify((Bm.T*Bm).inv()*Bm.T*imgs)
report('1(x)Rc preserves the span of the two path vectors', sp.simplify(Bm*W - imgs) == sp.zeros(8, 2))
print('face-weight matrix in the basis (b=1, b=0) with these path vectors:')
sp.pprint(sp.factor(W))
print('eigenvalues:', [sp.factor(e) for e in W.eigenvals()])
```
Output:
```
PASS e1,f1,k1,e0 equations: one free parameter (1) and solution = eq-six-vertex
PASS Rc(q^-2) has rank 3 (spin-1 image)
PASS fused operator maps V1(x)(x)V1(y) into V1(y)(x)V1(x)
PASS fused spin-1 R-matrix: eigenvalues 1 (x5), <2> (x3), <2><1> (x1) at z = x/y
PASS path vectors are highest-weight
PASS 1(x)Rc preserves the span of the two path vectors
face-weight matrix (b=1,b=0):
  [[(1-q^4 z)/((1+q^2)(z-q^2)),            q^2 (z-1)/(z-q^2)],
   [(z-1)(1+q^2+q^4)/((1+q^2)^2 (z-q^2)),  (z-q^4)/((1+q^2)(z-q^2))]]
eigenvalues: [1, -(q**2*z - 1)/(-q**2 + z)]
```

### ybe_R.py
```python
"""eq-ybe-braid in terms of R = P Rc: R12(x1/x2) R13(x1/x3) R23(x2/x3) = R23 R13 R12 (eq-yang-baxter)."""
import numpy as np
q = (0.83+0.41j)**2
def Rc(u):
    return np.array([[1,0,0,0],[0,u*(q**2-1)/(q**2-u),q*(u-1)/(u-q**2),0],[0,q*(u-1)/(u-q**2),(q**2-1)/(q**2-u),0],[0,0,0,1]])
P = np.eye(4)[[0,2,1,3]]
R = lambda u: P@Rc(u)
I = np.eye(2)
P23 = np.kron(I, P)
R12 = lambda u: np.kron(R(u), I); R23 = lambda u: np.kron(I, R(u)); R13 = lambda u: P23@R12(u)@P23
x1, x2, x3 = 0.7+0.2j, -0.4+1.1j, 1.3-0.5j
lhs = R12(x1/x2)@R13(x1/x3)@R23(x2/x3); rhs = R23(x2/x3)@R13(x1/x3)@R12(x1/x2)
print(('PASS ' if np.allclose(lhs, rhs) else 'FAIL ') + 'R12(x1/x2)R13(x1/x3)R23(x2/x3) = R23 R13 R12 with R = P Rc')
```
Output: `PASS R12(x1/x2)R13(x1/x3)R23(x2/x3) = R23 R13 R12 with R = P Rc`

---

## F. What is verified how, and what is left open

**Verified by explicit computation:**
- All verdicts in items 4, 11 and 12'.
- The d_{n+1}^{(2)} spinor results for n = 2, 3, 4.
- The ⟨a⟩₊ rule, for this example only.

**From literature only:**
- DGZ 1996 §4.3 text: read in the PDF, and its formula reproduced numerically.
- DGZ 1994's list of solved cases (B_n^{(1)} spinor, F₄ 26⊗26): from the arXiv HTML via a summariser, not re-derived. Someone should compare DGZ 1994 eqs. 4.21–4.25 with chapter 24's B_n^{(1)} spinor R-matrix.
- The Zhang–Gould–Bracken 1991 attribution: from DGZ 1996.
- Jimbo 1986's ξ values: from memory, consistent with the computation; Jimbo 1986 itself was not fetched.

**Unclear:** none of the assigned items. One thing to settle separately is whether the novelty statements in chapters 22–26 should be cut back to the S-matrices once the DGZ R-matrices are cited. That is an editorial decision for the author.
