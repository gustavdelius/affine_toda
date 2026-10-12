<!-- Opus referee follow-up, 2026-10-12. Model: Claude Opus 5.5, effort max. Prompt: briefs/halfline_open_12a_opus_prompt.md. Combined session with 03 (see halfline_open_03_referee_opus.md). Saved verbatim by the coordinator. Scripts: /tmp/claude-1000/-home-gustav-Git-affine-toda/0e38a99d-702a-40fb-99b7-74dbe0413cf6/scratchpad/ref_03_12a/ (s03/akr_*, s12a/). -->

## 12a

### Item 1, the AKR Δc_boun misprint: CORRECTED. The worker is right on substance but wrong in the h=2 factor; both of Gemini's extra claims are refuted.

**What AKR print** (batft.tex line 972): Δc_boun(b) = 6m̄R cos(π/2h) / [sin(πB/2h) sin(π(1−B)/2h)].
- **Their B.** AKR's B = b²/(1+b²) (line 277) is exactly Fateev's x (hep-th/0103014, line 1416). There is no rescaling; the book's B is 2x.
- **Fateev (eQ, line 1411):** E = sin(π/2h) Σ_m / [8 sin(πx/2h) cos(π(1−x)/2h)], with Σ_m(Aₙ) = 2m̄ cot(π/2h) (line 1367). AKR's m̄ is Fateev's m (AKR line 234, Fateev line 137).
- **The numerator is not a misprint.** sin(π/2h)·Σ_m = 2m̄ cos(π/2h), so AKR's cos(π/2h) is correct.

**Correct formula:** Δc_boun = (24R/π)E = 6m̄R cos(π/2h) / [π sin(πB/2h) cos(π(1−B)/2h)].
- The printed version is wrong in the denominator (sin instead of cos) and lacks 1/π.
- Printed/correct = π cot(π(1−B)/2h). At h=2 that is π cot(π(1−B)/4) = π tan(π(1+B)/4), not the worker's π cot(π(1+B)/4); checked numerically at B=0.2, 0.4, 0.7.
- The printed formula is symmetric under B→1−B, so it would give the (+) and free boundaries the same energy. That matters for their (+f) use, Δc(b)+Δc(1/b), at line 1065.

**New finding.** AKR's Δc_Bulk (line 964), 3m̄²R²/(2π)·sin(π/h)/(…), is too small by a factor 4. It should be 6m̄²R²/π·(…), which is 24R²/π times Fateev's (ebul).

**Propagation: purely typographical.** I parsed their Grace EPS figures.
- **Fig. 2 (A₂, (++)):** dashed − dotted = 16.74R (R from 0.002 to 0.04). The corrected 2Δc_boun is 16.73R; the printed one would be 161.8R. For h=3 the printed coefficient is ≥155R for every B, so no choice of B rescues it.
- **Solid − dashed** ≈ 6.9R² (R from 0.01 to 0.1). The corrected bulk term is 6.92R²; the printed one is 1.73R².
- **Fig. 1 (A₂, A₃, A₄ at B=0.4):** RA and TBA curves agree to ≤0.004 for R ≤ 0.3. The printed formula would shift them by 1.45, 2.80 and 4.56 at R=0.01 for h=3, 4, 5.
- **Tables 1–6** (δ₁, δ₃, δ₅) do not involve vacuum energies.

**Appendix wording.** "AKR [@akr2002], sec. 4: the boundary contribution should read $\Delta c_{\rm boun}=\frac{6\bar mR\cos(\pi/2h)}{\pi\sin(\pi B/2h)\cos(\pi(1-B)/2h)}$, i.e. $\frac{24R}{\pi}$ times Fateev's boundary energy, with their $B=b^2/(1+b^2)$ (half the $B$ of this book). As printed it has $\sin(\pi(1-B)/2h)$ in place of the cosine and lacks the factor $1/\pi$. Their bulk term should read $\frac{6\bar m^2R^2}{\pi}\frac{\sin(\pi/h)}{\sin(\pi B/h)\sin(\pi(1-B)/h)}$; as printed it is smaller by a factor 4. Both are misprints only: the curves of their Figs. 1 and 2 contain the corrected values."

### Item 2, "for n≥2 the only integrable value": CORRECTED (overclaim)

**What the sources state.**
- **Fateev (lines 984–1015):** classical integrability [BCD] allows Neumann or d_i=±1 with μ_B² = μ/πb² (ND). That is the classical |C_j|=1, with no free parameter.
- **Fateev (lines 1035–1050):** (MND) μ_B² = (μ/2)cot(πb²/2) is the "quantum modified version" of (ND). It is obtained for the conformal (non-affine) theory as the condition that the boundary preserves W-symmetry, via null vectors of strongly degenerate boundary fields; this is asserted ("It can be shown"), following FZZ for Liouville.
- **Fateev then states:** "The boundary ATTs will be considered as perturbed boundary CFTs."
- **AKR (lines 283–296):** "μ_B should be fixed completely to have conserved charges [CDRS, FZZ]".
- **Conclusion.** Nobody argues that the W-symmetric conformal boundary, perturbed, is the quantum-integrable affine boundary, nor that other quantum values are non-integrable. It is assumed.

**How it squares with |C_j|=1.** (MND) → (ND) as b→0, i.e. M/M₀ = cos(πb²/2) → 1. It is the renormalised value of the classical one.

**Evidence for the assumption** is the consistency of AKR's TBA: with the DG amplitude as input, the TBA reproduces δ₁ at s=ib/2 for n=2,3,4 and misses s=0. That comparison is the worker's, confirmed by Gemini; I did not redo it.

**Proposed sentence:** "For $n\ge2$ classical integrability leaves no free boundary parameter ($|C_j|=1$). The quantum value $\mu_B^2=\frac\mu2\cot\frac{\pi b^2}2$ is not derived from integrability. Fateev obtained it as the condition for the boundary of the conformal Toda theory to preserve $W$-symmetry, and he and Ahn, Kim and Rim assume that the integrable affine boundary is this conformal boundary perturbed by the affine terms [@fateev2001; @akr2002]. **[Conjecture, numerically supported]** With (-@eq-uniform-amplitude) as input, the boundary TBA reproduces the ultraviolet data at this coupling for $n=2,3,4$ [@akr2002]. As $b\to0$ the value tends to the classical $|C_j|=1$."

The worker's "for n≥2 it is the only integrable value" should be dropped.

### Item 3, the minor misprints: both CONFIRMED

**BRZ (0710.4789, bshpfinal.tex line 680, Bfugacity).**
- (R1) at line 533 has cosh(pη), and k=R⁽¹⁾(iπ/2−θ) gives (coshθ ∓ cosh pη) in the fugacity.
- I wrote my own sinh-Gordon boundary TBA from their equations at b²=0.8086, s₁=s₂=0.80i, root-finding mR so that P_TBA hits the table's P.
- **cosh reading:** Δ^TBA = 0.41762001692146 at P=0.0665109 (table: 0.41762001693341, difference −1.2e-11). The next two entries differ by +9.5e-12 and +7.9e-12.
- **cos reading as printed:** Δ^TBA ≈ −0.34 to −0.48, missing the table by 0.76–1.08.

**Takács (0801.0962, bffshg.tex lines 683 and 709).**
- **(zEFmatch) lacks i.** Matching (blogcfromexact) with (blogcfromcumul) needs ½πz coth πz = ¼π(E−F)cot ½π(E−F). That requires z = (i/2)(E−F) and z̄ = (i/2)(E+F); the printed form misses by 0.021, the corrected one matches to 6e-17.
- **(zetathetamatch) lacks i.** His (zuvir) together with Zamolodchikov's relations gives cosh²πz = cos²(β²(η+iϑ)), so z = (iβ²/π)(η+iϑ) and z̄ = (iβ²/π)(η−iϑ), up to the sign/branch of z. Numerically this matches to 7e-18; the printed form misses by 0.059.
- **The pairing of z with η±iϑ** is a convention, since the swap z↔z̄ is ϕ₀→−ϕ₀ and his C_kl is symmetric in z, z̄.
- **The slips cancel.** Equating the corrected forms, (i/2)(E∓F) = (iβ²/π)(η±iϑ), the common i cancels. So E = 2β²η/π → 2ξη/π at leading order, and F = ±2iξϑ/π, exactly as printed.

## Scripts (12a part)

All in `/tmp/claude-1000/-home-gustav-Git-affine-toda/0e38a99d-702a-40fb-99b7-74dbe0413cf6/scratchpad/ref_03_12a/`.
- `s03/parse_grace_eps.py`, `s03/akr_fig2.py`, `s03/akr_fig1.py`: AKR figure extraction and vacuum-energy coefficients (12a item 1).
- `s12a/brz_table1.py` (output `s12a/brz_table1.out`): own BRZ boundary TBA, cosh versus cos reading.
- `s12a/takacs_i.py`: the two missing factors of i.
