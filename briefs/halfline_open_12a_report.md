<!-- Worker report, 2026-10-11. Model: Claude Opus 5.5, effort xhigh. Scope: problem 12, the AKR lead and part (a) (brief halfline_open_12_rank_one.md); part (b) assigned separately. Saved verbatim by the coordinator. Scripts copied to briefs/halfline_open_12a_scripts/. -->

## Problem 12: the AKR lead and part (a)

**Summary.** The lead holds once one sign in the brief is corrected. Part (a) is settled numerically without a two-loop calculation, and the deciding test is mostly already in the literature: Bajnok, Rim and Zamolodchikov (BRZ), arXiv:0710.4789, Nucl. Phys. B 796 (2008); the brief did not list it. I re-ran their test at the uniform point itself. The DG/ε̂=0 amplitude (E=B/2) belongs to the sinh-Gordon boundary with M/M0 = cos(πb²/2), which is AKR's (+) coupling. The point M = M0 (the classical "C=1" in the conformal normalization) has E=0. No repo files were edited.

### 1. The lead

**Convention chain** (b² = β²/4π in the book's β, so b² = B/(2−B), and p ≡ b²/(1+b²) = B/2):

1. AKR's a_1 action becomes FZZ/BRZ form with φ = √2φ′: (1/4π)(∂φ′)² + 2μ cosh 2bφ′ + μ_B(e^{bφ′} + e^{−bφ′}). So M⁺ = M⁻ = μ_B.
2. FZZ: M = M0 cosh(πbs), with M0 = √(μ/sin πb²). AKR (10) then gives cosh(πbs) = cos(πb²/2), so s = ib/2. AKR state this themselves (their eq. "param").
3. **Correction to the brief's continuation.** The brief's minus sign is an artefact. Under b² → −1/(λ+1) the continued μ is negative: AKR's mass–μ relation forces this for b² < 0, and it puts the sine-Gordon vacuum at 0. So μ_SG = −2μ, not 2μ, and likewise μ̃ = −2μ_B. Both μ and sin(β²/8) flip sign, the two flips cancel, and (μ̃/μ̃_crit)² = M²/M0² = +cos²(π/(2(λ+1))).
4. Then BPT (3.2) at φ0 = 0 gives η = π/2 and ϑ = 0, i.e. ε̂0 = ε̂1 = 0. Checked at λ = 1.7, 5.3, 12.3, to 1e−30.
5. In sinh-Gordon, cos(πE/(2−B)) = M/M0, which is CT's E = 2a(1−B/2) with C = cos πa ≡ M/M0. So s = ib/2 gives E = p = B/2, and s = 0 gives E = 0. Checked against BRZ's (etas) as well.

**Where (10) comes from.** It is not a fit. Fateev (hep-th/0103014, eq. MND) calls it the "quantum modified version" of the classical μ_B² = μ/πb². He obtains it, following FZZ, from the null-vector equations of degenerate boundary fields, which he links to preserving W-symmetry. For n ≥ 2 it is the only integrable value of the coupling.

**Does AKR's TBA tell the two couplings apart? Yes, decisively.**
- I reproduced AKR's δ₁ table at s = ib/2: all 12 entries, B = 0.25–0.80, to 6 digits.
- At s = 0 (M = M0) δ₁ shifts by 0.46 to 8.9, against TBA agreement of ≤ 1e−5.
- Note that AKR's TBA covers A_2, A_3, A_4, which are the book's a_2^(1), a_3^(1), a_4^(1) (n = 2, 3, 4).

**A second, independent check through the boundary energy.**
- Fateev's conjectured (+) boundary energy, (eQ) at h = 2, equals BPT/BRZ's sinh-Gordon boundary energy at E = p identically. It agrees to 15 digits at p = 0.05, 0.2, 0.3, 0.4471, 0.7.
- It differs from the value at E = 0, and its dual matches E = 1−p (Neumann).
- **Literature correction:** AKR's printed Δc_boun is wrong. It has sin(π(1−B)/2h) where Fateev has cos(π(1−x)/2h), and it lacks a factor 1/π. Equivalently it is off by a factor π·cot(π(1+B)/4) at h = 2.

### 2. Part (a): which identification the evidence supports

**My TBA test** (`shg_btba_uniform.py`). This is the sinh-Gordon boundary TBA on an interval with both edges symmetric. I extract the Liouville boundary phase Δ^TBA(P) and compare it with the FZZ/BRZ boundary Liouville reflection phase Δ(P|s,s). I first validated it against BRZ's Table 1 (b² = 0.8086, s = 0.80i): agreement 4e−11 at mR = 1e−8, 7e−12 at 1e−4 and 5e−8 at 1e−2.

| b² | TBA input | Agrees with s = ib/2 (AKR) to | Misses s = 0 (M = M0) by | 
|---|---|---|---|
| 0.8086 | E = B/2 (ε̂=0) | 0.5–3·10⁻¹¹ (mR = 1e−8 to 1e−4); 6e−9 at 1e−2 | 0.06–0.20 |
| 0.4286 | E = B/2 | 2–3·10⁻¹¹; 3e−6 at 1e−2 | 0.025–0.09 |

The E = 0 runs are the mirror image. They match s = 0 to 3e−11 and miss s = ib/2 by the same margins. The residual discrepancy at mR = 1e−2 has the size expected for R^{2bQ}-type power corrections.

**Conclusion for (a).**
- The quantum uniform boundary has M/M0 = cos(πb²/2) for every n. For n ≥ 2 this is fixed by integrability and confirmed by AKR's TBA. At n = 1 the same point is ε̂ = 0, with E = B/2.
- M = M0 is a different member of the sinh-Gordon family, with E = 0. It has no counterpart for n ≥ 2, where it is not integrable.
- So the book's reconciliation is correct, and sharper than the book states it. A two-loop calculation is no longer needed to decide this.
- For reference, the two boundary energies differ by −mπp/16 + O(p²), i.e. at two-loop order.

### 3. Bearing on (c): status of `cnj-bsg-uv-ir`

- **Analytic route:** on the line φ0 = 0, ϑ = 0, BPT derive the relation from the continued TBA boundary energy plus the FZZ boundary VEV. That needs a factor ½ in log g_S, confirmed privately by Zamolodchikov.
- **Numerical route:** BRZ confirm the sinh-Gordon form non-perturbatively in the physical regime, for symmetric edges (real and imaginary s) and asymmetric edges (ϑ ≠ 0) at b² = 0.8086. My runs add the ε̂ = 0 and E = 0 points at b² = 0.8086 and 0.4286.
- **Takács, 0801.0962:** this is Takács alone, not "Kormos–Takács" as the brief says. Its boundary form-factor check is only O(g²).
- **Proposed status:** "[Numerical] in its sinh-Gordon form (0 < b² < 1); [Conjecture] for sine-Gordon". The remaining conjectural step is the continuation b² → −β²_SG/8π.

### 4. Replacement text

**"Rank one" paragraph of `sec-uniform-boundary`:**

> **Rank one.** For $n=1$, (-@eq-uniform-amplitude) is $K_1=(2-\frac B2)/\big((1)(1-\frac B2)\big)$, the sinh-Gordon factor (-@eq-shg-reflection) at $E=\frac B2$, $F=0$: Ghoshal's factor at $\hat\epsilon=0$, continued (@exr-shg-uniform). To say which Lagrangian boundary this is, normalize the couplings as in perturbed conformal field theory, $\int(\frac1{4\pi}(\partial\phi)^2+2\mu\cosh2b\phi)+M\int(e^{b\phi}+e^{-b\phi})$ with $b^2=\beta^2/4\pi=B/(2-B)$, and write $M=M_0\cosh\pi bs$, $M_0=\sqrt{\mu/\sin\pi b^2}$ [@fzz2000]; classically $C_0=C_1=1$ is $M/M_0\to1$. The continuation of @cnj-bsg-uv-ir then reads $\cos\frac{\pi E}{2-B}=M/M_0$, which is Chenaghlou and Corrigan's relation (@cnj-shg-e) with $C=M/M_0$. It puts $E=\frac B2$ at $M/M_0=\cos\frac{\pi b^2}2=\cos\frac{\pi B}{2(2-B)}$ and $E=0$ at $M=M_0$. Both assignments have been tested non-perturbatively. The boundary TBA of sinh-Gordon theory on an interval, built from (-@eq-shg-reflection), reproduces the ultraviolet behaviour predicted by the boundary Liouville reflection amplitude of the Lagrangian boundary to about $10^{-10}$ [@brz2008]; with $E=\frac B2$ it does so at $M/M_0=\cos\frac{\pi b^2}2$ and fails at $M=M_0$ (`code/shg_boundary_tba.py`). For $n\ge2$ there is no free boundary parameter: the coupling of the uniform boundary is fixed, as Fateev argued from the $W$-symmetry of the conformal boundary, at $\mu_B^2=\frac\mu2\cot\frac{\pi b^2}2$, again $M/M_0=\cos\frac{\pi b^2}2$ [@fateev2001], and Ahn, Kim and Rim confirmed by the boundary TBA for $n=2,3,4$ that (-@eq-uniform-amplitude) is the reflection amplitude of exactly this boundary [@akr2002]. So the classical $C\equiv1$ is renormalized to $\cos\frac{\pi b^2}2$ for every $n$, and for $n=1$ this renormalized uniform boundary is the point $\hat\epsilon=0$. The point $M=M_0$, which is also $C\equiv1$ classically, is a different member of the sinh-Gordon family, with $E=0$. The two factors differ only at $O(B^2)$, so one-loop perturbation theory cannot distinguish them (`code/particle_reflection.py`). For sine-Gordon the same point is $\tilde\mu/\tilde\mu_{\rm crit}=\cos\frac\pi{2(\lambda+1)}$ (@sec-bsg-parameters).

**Open problem 12(a):** replace it with "(a) For sinh-Gordon the identification is settled numerically (@sec-uniform-boundary); an analytic derivation of the UV–IR relation off the line $\phi_0=0$, and its continuation to sine-Gordon, remain open."

**`cnj-bsg-uv-ir` status paragraph:** add after the TCSA sentence:

> In its sinh-Gordon continuation the relation has been tested directly in terms of the Lagrangian couplings, with no TCSA input: Bajnok, Rim and Zamolodchikov compared the boundary TBA on an interval with the boundary Liouville reflection amplitude, for symmetric and asymmetric boundaries, and found agreement to many digits [@brz2008]. What remains conjectural is the continuation to $b^2<0$.

**New bib keys needed:**
- `brz2008`: Bajnok, Rim, Zamolodchikov, "Sinh-Gordon boundary TBA and boundary Liouville reflection amplitude", Nucl. Phys. B 796 (2008); arXiv:0710.4789. Verify the page number (I recall 622).
- `fzz2000`: hep-th/0001012.
- `fateev2001`: hep-th/0103014.

### 5. Other findings
- **BRZ (Bfugacity)** prints cos(η_a p) where their (R1) has cosh(pη). Only the cosh reading reproduces their Table 1.
- **Takács (zEFmatch) and (zetathetamatch)** each lack a factor i. The two errors cancel in his final E = 2ξη/π.
- **Book ch. 29, "Self-duality" paragraph** calls B = 2 "the self-dual point". By `eq-b-of-beta` the self-dual point is B = 1, and B → 2 is the strong-coupling limit.
- **Status of the rest of my scope:** step 3 (two loops) was not started; it is unnecessary. Part (b) was left to the other agent.

### 6. Files (all in `/tmp/claude-1000/-home-gustav-Git-affine-toda/0e38a99d-702a-40fb-99b7-74dbe0413cf6/scratchpad/`)
- `scripts/rank_one_lead_check.py`: convention chain, the Fateev/BPT energy identity and AKR's δ₁ table at s = ib/2 and s = 0. Output in `scripts/lead_out.txt`.
- `scripts/shg_btba_uniform.py`: the sinh-Gordon boundary TBA against the boundary Liouville phase, plus the BRZ validation. Runs in about 3 minutes. Output in `scripts/btba_out.txt`. It is the candidate for `code/shg_boundary_tba.py`.
- `dl/`: downloaded arXiv sources (AKR, BPT, BPT review, FZZ, Fateev, BRZ, Takács, Chenaghlou, hep-th/0004121, BBT). These are untrusted data; they were only read.

Sources: [BRZ, arXiv:0710.4789](https://arxiv.org/pdf/0710.4789), [AKR, hep-th/0110218](https://arxiv.org/pdf/hep-th/0110218v1), [BPT, hep-th/0108157](https://ar5iv.labs.arxiv.org/html/hep-th/0108157)

## Integration (branch half-line-open-problems-integration)

| Claims | Verdict | How it was checked | Book change |
|---|---|---|---|
| 1–3 (convention chain; $s=ib/2$; $E=B/2$ at $M/M_0=\cos\frac{\pi b^2}2$, $E=0$ at $M=M_0$) | **Verified, integrated** | The coordinator recomputed the algebra at $B=0.3,0.6,0.9$ (exact). The continuation-sign question (claim 2) is not needed for the sinh-Gordon statements in the book. | "Rank one" paragraph of `sec-uniform-boundary` rewritten. |
| 4 ("for $n\ge2$ the only integrable value") | **Rejected as stated; modified** | Opus follow-up (`halfline_open_12a_referee_opus.md`, item 2), checked by the coordinator in Fateev's TeX (lines 1030–1052) and AKR's (lines 283–296). (MND) is Fateev's condition for the *conformal* Toda boundary to preserve $W$-symmetry; that the integrable affine boundary is this boundary, perturbed, is assumed by Fateev and AKR, not derived. The coordinator's first edit ("fixed by integrability") was corrected. | The paragraph now says the value is not derived from integrability, gives Fateev's origin, and cites AKR's TBA as numerical support; "On this evidence" qualifies the conclusion. |
| 5 (AKR's TBA covers $a_2,a_3,a_4$) | **Verified** (TeX table captions: "$A_2,A_3$, and $A_4$ ATFTs for $(++)$ BC"; note AKR's $B=b^2/(1+b^2)$ = book $B/2$) | | Cited in the paragraph. |
| 6 (Fateev $(+)$ energy = sinh-Gordon energy at $E=B/2$; difference $-m\pi p/16$) | **Verified** (recomputed at four values of $p$; leading term at $p=10^{-3}$) | | Sentence in the paragraph; check 2 of the new `code/shg_boundary_tba.py`. |
| 7 (boundary TBA: $E=B/2$ ↔ $s=ib/2$; $E=0$ ↔ $s=0$) | **Verified** | The coordinator reran Gemini's independent TBA code and wrote `code/shg_boundary_tba.py` (FFT-iterated TBA, Liouville phase from FZZ/BRZ). It reproduces BRZ table 1 ($10^{-15}$) and gives hit/miss $\le4\times10^{-15}$ / $\ge0.03$ at $b^2=0.8086,0.4286$. All 4 checks PASS (40 s). | Paragraph; status paragraph of `cnj-bsg-uv-ir` (ch. 29) gains the BRZ sentence. |
| 8 (status: [Numerical] in sinh-Gordon form, [Conjecture] for sine-Gordon) | **Adopted in substance** | The label of `cnj-bsg-uv-ir` stays [Conjecture], since the conjecture is the sine-Gordon statement. The status paragraph now says the sinh-Gordon form is tested directly [brz2008]. BRZ's abstract was read; the worker's claim that BRZ also tested asymmetric boundaries was not checked, and the text does not make it. | ch. 29 status paragraph. |
| 9, 10 (AKR $\Delta c_{\rm boun}$ misprint; BRZ, Takács slips) | **Settled by the Opus follow-up; integrated (AKR, BRZ)** | AKR: printed denominator has $\sin(\pi(1-B)/2h)$ for $\cos$ and lacks $1/\pi$; the numerator $\cos(\pi/2h)$ is correct; at $h=2$ the error factor is $\pi\tan(\pi(1+B)/4)$ (the worker's $\cot$ was wrong; Gemini's extra claims refuted). New: AKR's bulk term is too small by a factor 4. Both are misprints only (their figures use the corrected values). BRZ (Bfugacity) and Takács's two missing factors of $i$ confirmed. | Appendix D: AKR (both misprints) and BRZ (Bfugacity) added after the coordinator re-derived the AKR formula from Fateev's (eQ). Takács not added: the book does not cite or use his paper, and his slips cancel. |
| 11 (ch. 29 "self-dual point $B=2$") | **Verified, corrected** | By `eq-b-of-beta`, $B\in[0,2)$ and duality is $B\to2-B$, so the self-dual point is $B=1$; $B\to2$ is $\beta\to\infty$. | ch. 29 paragraph renamed "Strong coupling and duality" and corrected. |
| 12 (text) | Integrated as modified above; open problem 12 rewritten. New bib entries `brz2008` (arXiv abstract page: NPB 796 (2008) 622–650, DOI checked), `fzz2000` and `fateev2001` (@misc; no journal reference on arXiv). | | |
