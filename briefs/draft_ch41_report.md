# Chapter 41 drafting report (2026-10-12)

Drafted by two agents in sequence (exploration by the first, recorded in `briefs/draft_ch41_checkpoint.md`; helpers, script, chapter and this report by the second). Nothing committed.

Files:
- `part8-defects/41-soliton-transmission.qmd` (275 lines; skeleton replaced, all skeleton ids kept).
- `textbook_code/defect_an_transmission.py` (new, 21 checks, all PASS, ~100–130 s depending on machine load).
- `textbook_code/defect_common.py`: ch. 41 section appended after the marker line only (`git diff` is a single hunk from line 426 on). Helpers: `an_root_window`, `an_weyl_module` (with `orient`), `an_vweights`, `an_solve_intertwiner` (vectorized; drops unknowns that occur in no equation, after which the *whole* null space is one-dimensional), `an_T_vector`, `an_T_plus_vbar`, `an_block_inverse`, `defect_crossing_residual` (vectorized), `an_tensor_rep`, `an_triangle_residual` (sparse), `an_scalar_g`, `an_scalar_gbar`, `an_scalar_constant`.

## 1. Chapter summary

Conventions: Part VI coproduct; q = −e^{−iπω}, ω = 4π/β² − 1 (CZ's γ; the chapter says explicitly that ω is the coupling here, not ch. 39's root of unity, which it writes e^{2πi/h}); y = e^{ωθ}; X = e^{ω(θ−θ_d−iπ/2)} = CZ's x̂ with Δ = θ_d + iπ/2.

Sections (skeleton ids kept):
- `sec-t-procedure`: six-step procedure (Borel module, intertwiner, crossing, scalar factor, poles, semiclassics); the a_r Weyl module `eq-dt-weyl` (labels N + s, N in the root lattice, offset s = λ_p for class p; f_j|N⟩ = (μ/y_d) q^{(N+s)_{j+1}}|N − α_j⟩, μ = −e^{−iπω/2}/(q − q^{-1})); one-line proof of q-Serre (f_j f_{j+1} = q^{-1} f_{j+1} f_j, and q-Serre ⟺ c = q^{±1}); the orientation (q^{−(N+s)_{j+1}} = charge-conjugate defect, cz2007's D^T); uniqueness/hypothesis (ii) on windows; the Penati–Polvara caveat (axioms: established/numerical; Lagrangian identification: conjecture).
- `sec-an-soliton-transmission`: closed forms `eq-dt-tv` (V, two entries), `eq-dt-tvbar` (V̄, all h entries, |X|^j), `eq-dt-tplus` (type I₊, V̄ linear). `prp-an-transmission` **[Established; numerical re-derivation]**: (1) uniqueness = CZ 2009a (5.8), (5.11) in the same gauge; (2) type I₊; (3) triangle relation; (4) selection rule: M_V obeys the classical rule exactly, M_V̄ tunnels to the forbidden weights (1 ≤ j ≤ h−2) with relative size min(|X|^j, |X|^{j−h+1}) → 0 semiclassically for θ ≠ θ_d, = 1 at θ = θ_d; roles of V, V̄ exchanged for I₊. Pointers: fused I₋I₊ pair open (`sec-defect-open-problems`), Robertson's fused a_3 matrix not checked.
- `sec-transmission-unitarity-crossing`: non-unitarity for r ≥ 2; conjugation C = 1 (dual V(y) = V̄(ye^{iπω})); `eq-dt-crossing` with c_I = 1 + (−q)^r X^h, c_II = 1 − qX^h; scalar relations and derivation of cz2009a (6.8); universal-R normalization ⇒ s-independent scalars; bootstrap `eq-dt-bootstrap` (ḡ = ∏ g); convergent scalar factor `eq-dt-g` (c_2,3,4 = 0.3983, 0.3933, 0.3910 at ω = 1.37). `prp-transmission-crossing` **[Numerical]** incl. the q-sign item (h odd: no conjugation at the physical crossing point at −q, crossing fails for C = 1 and C′; h even: fails with C = 1, holds with C′ = diag((−1)^i)), stated convention-qualified. The cz2009a (6.10) pairing correction.
- `sec-a22-transmission`: CZ 2011 Sec. 4 in Part VI form. `prp-a22-transmission` **[Numerical]**: relations (incl. quintic Serre), uniqueness for four (f1, f2), triangle with Ř, crossing at the algebraic crossing point x_c = −q^{−6} with antidiagonal C, q-sign only relative to C. Explicit "not done" list.
- `sec-an-defect-bound-states`: `eq-dt-pole` (θ_1 and block determinant q^r(1 + (−q^{-1})^r X^h)); residue rank h − 1 per sector; complex energies and `eq-dt-window` (max(0, a − h/2) ≤ 1/ω < a). `prp-a2-unstable` **[Numerical; the windows Established]**: (1) a_2 windows in book labels — species 1: ω > 1 (contains the classical limit), species 2: 1/2 < ω ≤ 2 = cz2007's window for their species 1; (2) residue rank 2; (3) kernel bootstrap: on the kernel of M_V(θ_1), Ř₁₂ then 𝕋₂₃ is the matrix with offset s + λ_1 and X^h → q²X^h, i.e. θ_d → θ_d + 2πi/h (mod 2πi/(hω)) — the quantum form of `eq-ds-class-shift`. Then: corrected cz2009a window; "absorption changes the class" (class p ⇒ offset λ_p AND θ_d + 2πip/h; general p is an expectation, checked for p ≤ 2 matrix parts); the image of the residue is (h−1)-dimensional and its meaning is open; charged classes: species-a pole at θ_d − i(π/2 − π(a+2p)/h) − iπ/(hω), window with a → a + 2p; classical real absorptions (a + 2p ≡ h/2 mod 2h) become decaying resonances at θ_d − iπ/(hω) (e.g. species 1 of a_5 in class 1); poles in the physical strip (species a > h/2 below the window; at θ_d = 0 candidate stable excited defects, not checked).
- `sec-transmission-semiclassics`: Wigner relation with prefactor 2π/(h(ω+1)) = β²/(2h). `prp-transmission-semiclassical` **[Numerical]**: deviations 0.092/0.025/0.0063 (a_2), 0.067/0.018/0.0046 (a_3) at ω = 10/40/160, same for ḡ. arg z paragraph (arg z_1(θ_d) = π/h; quantum branching |X|²/(1+|X|²) → step at θ_d). Lightest breathers: forward pointer to `sec-breather-transmission`, `prp-lightest-breather-universal`.
- Exercises (8): `exr-dt-serre`, `exr-dt-intertwiner`, `exr-dt-det`, `exr-dt-crossing`, `exr-dt-gamma`, `exr-dt-window`, `exr-dt-qsign`, `exr-dt-semiclassical`.
- Notes with history and five corrections/clarifications (section 4 below).

New ids: `eq-dt-weyl`, `eq-dt-tv`, `eq-dt-tvbar`, `eq-dt-tplus`, `eq-dt-crossing`, `eq-dt-bootstrap`, `eq-dt-g`, `eq-dt-pole`, `eq-dt-window`, the eight `exr-dt-*`. Skeleton ids kept: `sec-soliton-transmission`, `sec-t-procedure`, `sec-an-soliton-transmission`, `prp-an-transmission`, `sec-transmission-unitarity-crossing`, `prp-transmission-crossing`, `sec-a22-transmission`, `prp-a22-transmission`, `sec-an-defect-bound-states`, `prp-a2-unstable`, `sec-transmission-semiclassics`, `prp-transmission-semiclassical`, `sec-soliton-transmission-exercises`, `sec-soliton-transmission-notes`. All cross-references resolve to existing ids (grep), all citation keys are in references.bib (no new keys needed).

## 2. Script output (`python3 textbook_code/defect_an_transmission.py`, 132 s at load ~2; 97–100 s earlier)

21 lines, all PASS (abridged):
```
PASS eq-dt-weyl: ... Cartan and q-Serre relations [9e-16], level zero [3e-15]; without the q-powers q-Serre fails [1.21]
PASS prp-an-transmission (1): ... (nullity, deep rank, residual, gap) a_2 V: 1, 1, 2e-14, 0.02; a_2 Vbar: 1, 1, 2e-14, 0.02; a_3 V: 1, 1, 1e-14, 0.02; a_3 Vbar: 1, 1, 2e-14, 0.02
PASS eq-dt-tv, eq-dt-tvbar intertwine U_q(b^-) for a_2, a_3, a_4 and s = 0, lambda_1, generic [1e-15]
PASS prp-an-transmission (2), eq-dt-tplus: U_q(b^+) unique (a_2, s = 0 and generic); V: 3 entries; Vbar: 2 entries, closed form 2e-14; eq-dt-tplus intertwines for a_2, a_3, a_4
PASS orientation -1 (cz2007's D^T): unique, V: 1 -> [1, 2, 3], Vbar: 1 -> [1, 2]; level zero: with k_0 -> 1.3 k_0 no intertwiner (0.02)
PASS thm-t-triangle for a_r: ... [1e-15]; control: matrices at -q [0.89]
PASS thm-t-triangle (ii) for a_2: ... nullity 1 on the window (31 labels, gap 0.01; 13 s)
PASS eq-dt-crossing: crossing (I), (II) with C = 1 [8e-15], c_I = 1 + (-q)^r X^h, c_II = 1 - q X^h, a = 1 [1e-14], r = 2, 3, 4, s = 0, lambda_1, generic
PASS the sign of q: ... crossing (I) at -q_phys with C = 1 and with C': a_2: 1.00, 0.76; a_3: 1.00, 2e-15; a_4: 0.96, 0.52
PASS eq-dt-g: (6.8), (I), (II) for r = 2, 3, 4, gamma = 1.37, 0.6 [4e-09]; c_2 = 0.3983, c_3 = 0.3933, c_4 = 0.3910
PASS cz2009a (6.10) as printed: K -> 2K changes the z-dependent log by -z ln 2 [1e-04] (|z ln 2| = 0.7); eq-dt-g changes by 3e-05
PASS eq-dt-bootstrap: ... M_Vbar [2e-15] with ratio 1 [5e-16] (a_2, a_3)
PASS eq-dt-pole: ... rank h - 1 in every sector: a_2: 12 sectors, sigma_min 2e-16, next 1.00; a_3: 16, 2e-16, 0.71; a_4: 5, 3e-16, 0.62
PASS eq-dt-window = cz2009a (6.13)-(6.15): ... species 1 for gamma > 1.001, species 2 for 0.500 < gamma <= 1.999; cz2009a text fails at 542 grid points (a_4, a = 1)
PASS prp-a2-unstable (3): kernel one-dimensional per sector [2e-16], invariant [2e-15], offset s + lambda_1 [2e-15], cycle invariant q^2 [7e-16] (s = 0, lambda_1)
PASS prp-transmission-semiclassical: a_2: 0.0922, 0.0248, 0.0063; a_3: 0.0673, 0.0180, 0.0046 (g and gbar equal); arg z_1(theta_d) = pi/h
PASS a_2^(2): relations [2e-15], module incl. quintic Serre [6e-16], F = 1, (-1)^N fail [2.00]; intertwiners unique for four (f1, f2)
PASS a_2^(2): triangle [5e-15]; control (perturbed R) [0.10]
PASS a_2^(2): dual V(x) = V(x x_c), x_c = -q^-6, antidiagonal C; crossing with C [1e-14], at -q with C(-q) [1e-14], fails at -q with C(q) [1.00]
PASS prp-an-transmission (4), selection rule: V 2, Vbar h entries; forbidden/largest at |theta - theta_d| = 0.3, gamma = 1.37, 10, 40: 0.66, 5e-02, 6e-06; 1 at theta_d
PASS for r >= 2 the matrix part is not unitary for real theta (off-diagonal of M^+M >= 0.35)
```
Window sizes (labels/interior/deep): a_2 R2=18: 37/19/7; a_2 R2=14: 31/13/1; a_3 R2=18: 177/55/13; a_3 R2=14: 135/43/1; a_4 R2=8: 171/21/1; a_4 R2=12: 331/51/1. With unconstrained unknowns dropped, the full null space is one-dimensional (gap 0.01–0.02), so the deep-label rank (1 in all cases) is no longer the only evidence. Section (g) (a_2^(2), dense SVDs) takes ~50–70 s of the runtime.

## 3. Verification record

**By script** (check labels as above): q-Serre and level zero of `eq-dt-weyl` (both Borels, both orientations, control); uniqueness for a_2, a_3 (b^-) and a_2 (b^+); closed forms `eq-dt-tv`/`eq-dt-tvbar`/`eq-dt-tplus` intertwine (a_2–a_4, three offsets); orientation −1 structure; level zero necessary; triangle relation (a_2: VV, V̄V̄, VV̄; a_3: VV, V̄V̄; control at −q); hypothesis (ii) for a_2 (generic s); crossing (I)/(II) with exact c_I, c_II and s-independence (r = 2–4); dual/conjugation and the q-sign test (r = 2–4); g: (6.8), (I), (II) (r = 2–4, two couplings), c_r real positive; (6.10) drift vs convergent pairing; bootstrap ratio 1 (a_2, a_3); simple poles of g, ḡ; block determinant formula and residue rank h − 1 (a_2–a_4); window formula vs (6.13)–(6.15) (h = 3–6) and a_2 windows; kernel bootstrap class 0 → 1 and 1 → 2 (a_2); Wigner limit (a_2, a_3, g and ḡ); a_2^(2) items 1–3; selection-rule structure and suppression; non-unitarity.

**By hand**: the q-Serre argument (c = q^{±1}); the gauge identification with cz2009a (5.8), (5.11) (α ↔ N + s, l_i ↔ e_i projected, Q = 1/q, x̂ = X, Δ = θ_d + iπ/2); derivation of (6.8) from (I), (II); the universal-R reading of the diagonal entries (leading Cartan factor; the universal R of U_q(a_r^(1)) was NOT evaluated on 𝒲); the Wigner prefactor; class-p pole positions, the window with a → a + 2p, real absorptions → resonances at θ_d − iπ/(hω); the physical-strip statement for a > h/2 (and candidate stable states at θ_d = 0); the arg z / branching-ratio argument; cz2007 (6.14) sign (checked against the source text); cz2009a (6.13)–(6.15) and the window text (source text); cz2011 (4.7)–(4.9) (source text; (4.8) is garbled in the text layer, read as (b_1(−)^N + c_1)q^{−2N} + (b_2(−)^N + c_2)q^{2N}, consistent with the failing controls).

**Unchecked** (stated as such in the chapter): Robertson's fused a_3 matrix (1408.1960 not available); cz2007 a_2 matrices entry by entry (only the structure, via the orientation −1 module); uniqueness for r ≥ 4 and hypothesis (ii) for r ≥ 3; scalar factors of species 1 < a < r (cz2009a (6.11)–(6.12)) — used only for the class-p windows of other species; the excited-defect scalar factor; the meaning of the (h−1)-dimensional image of the residue; the type-I₊ species-1 closed form and the I₊ class shift; whether the candidate physical-strip poles at θ_d = 0 survive in 𝕋_V̄; a_2^(2): gauge invariants vs CZ (3.11)/(4.14), scalar factor (3.12), physical x/q dictionary, the unstable windows quoted from CZ (3.13) and after (their pole positions are not quoted: the text layer leaves the factor of ξ ambiguous against the Γ arguments), semiclassics vs the ch. 38 delays, the lightest breather (3.14)–(3.15); whether a Klein-type redefinition relates ±q for h odd.

## 4. Corrections to the sources (in the chapter Notes)

1. cz2009a (6.10): Gamma functions paired across k; partial products drift like K^{−z} (θ-dependent); converges only with a K-dependent normalization; book form `eq-dt-g` pairs at equal k (each factor 1 + O(k^{−2})). [script]
2. cz2009a Sec. 6 after (6.15): for a < h/2 the window is 1/γ < a, not 1/γ < h/2 − 1 / (h−1)/2; agree only for the largest a < h/2. [script + source text]
3. cz2007 Sec. 6 (clarification): their species 1 (window 1/2 < γ < 2) is the book's species 2 (D vs D^T); in book labels species 1 is unstable for γ > 1, including the classical limit. [script]
4. cz2007 (6.14): x̂ = (−q)^{+1/2} e^{γ(θ−η)} with (−q)^{1/2} = e^{−iπγ/2}, not (−q)^{−1/2}; with + the entries coincide with the SG matrix quoted before it. [source text, by hand]
5. cz2009a Sec. 5 (clarification): the "no zeros" of (5.11) made quantitative (tunnelling ∝ |X|^j, vanishing semiclassically). [script]

## 5. Requests and obligations

**Edits elsewhere (not made):**
- `textbook_code/README.md`: row for `defect_an_transmission.py` (ch. 41, 21 checks, ~2 min) and the ch. 41 helpers in `defect_common.py`.
- Ch. 36 `sec-defect-bound-states` ("Unstable states" paragraph) and ch. 37 end of `sec-sg-defect-spectrum`: both say the a_2 unstable state exists "only in a window of couplings". In book labels that holds for species 2 only; species 1 is unstable for ω > 1, including the classical limit (`prp-a2-unstable`). Suggest "for one of the two species only in a window of couplings (@sec-an-defect-bound-states)".
- Ch. 39, the paragraph after `prp-defect-frozen-soliton` (in `sec-classical-transmission-factors`), points the I₋I₊ open question to `@sec-an-soliton-transmission`; ch. 41 only restates it as open and points to `@sec-defect-open-problems` — optionally retarget ch. 39's pointer.
- `briefs/part8_defects_status.md`: ch. 41 row and record (this report).

**Conventions-brief updates:**
- a_r Weyl module `eq-dt-weyl`, μ = −e^{−iπω/2}/(q − q^{-1}), X = e^{ω(θ−θ_d−iπ/2)} = CZ's x̂ (Δ = θ_d + iπ/2); orientation +1 = book type I₋ (π), −1 = charge-conjugate defect (cz2007's D^T, V full, V̄ linear with π^{−1}).
- Universal-R normalization for a_r: diagonal q^{−wt·(N+s)}; crossing scalars c_I = 1 + (−q)^r X^h, c_II = 1 − qX^h, s-independent; C = 1 (v_i^* ↦ v̄_i), dual V(y) = V̄(ye^{iπω}); at −q: C′ = diag((−1)^i) for h even, none for h odd.
- Quantum class shift (I₋): class p ⇔ offset λ_p AND θ_d + 2πip/h (mod 2πi/(hω); classical representative).
- `defect_crossing_residual` uses C transposed (matters only for non-symmetric C, e.g. the antidiagonal a_2^(2) conjugation).
- Notation clash in Part VIII: ch. 41 uses ω for the coupling (Part VI, CZ's γ), ch. 39 uses ω = e^{2πi/h} and γ_a = πa/h. Ch. 41 states its choice explicitly; a global decision may be wanted.
- Bib: no new keys.

**Obligations for ch. 42–43:**
- Ch. 43 `prp-lightest-breather-universal`: bootstrap the lightest-breather factor from `eq-dt-tv`, `eq-dt-tvbar`, `eq-dt-g` at the V⊗V̄ breather pole (cz2007 Sec. 7 for a_2) and compare with `eq-an-linear-transmission`; a_2^(2) breather (cz2011 (3.14)–(3.15)) unchecked.
- Ch. 43 open problems: the (h−1)-dimensional image of the residue; candidate stable excited defects at θ_d = 0 (species a > h/2 for 1/ω < a − h/2, e.g. species 2 of a_2 for ω > 2); type-I₊ species-1 closed form and I₊ class shift; the I₋I₊ fused pair vs the soliton S-matrix; a_2^(2) scalar factor/dictionary/windows; Robertson's fused matrices; Penati–Polvara vs the axiomatic matrices; q-sign for h odd beyond fixed conventions; hypothesis (ii) for r ≥ 3; universal R of U_q(a_r^(1)) on 𝒲.
- Ch. 42: a_r defect–boundary fusion (prp-defect-fusion-k analogue) can use 𝕋̄K𝕋 with these matrices and g.

**Earlier obligations addressed:** ch. 37 (C vs C′ crossing test → q-sign item; Wigner test → `prp-transmission-semiclassical`; CZ 2010 (3.17) pole sign was already resolved in ch. 37 `sec-type-two`, nothing further); ch. 38 (`eq-dc-ar-delay` used in the Wigner test; Tzitzéica delay labels: referred to ch. 38 where corrected, semiclassical comparison not done); ch. 39 (selection rule → `prp-an-transmission` (4); Wigner with |z| and arg z → `sec-transmission-semiclassics`; charged sectors at θ_d + 2πip/h → kernel bootstrap; absorption poles in charged classes vs cz2009a windows → "Charged classes"); ch. 40 (a_r Borel modules and hypothesis (ii) → `sec-t-procedure`; s-dependence via universal-R normalization → "The offset"; convention-qualified q-sign → "What crossing says about q").
