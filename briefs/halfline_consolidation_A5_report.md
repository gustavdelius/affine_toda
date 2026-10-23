<!-- Consolidation agent report A5, 2026-10-13. Model: Claude Opus 5.5, effort high. Prompt: briefs/halfline_consolidation_common.md + briefs/halfline_consolidation_A5_prompt.md. Interrupted once by the usage limit (before any computation) and resumed. Saved verbatim by the coordinator. Scripts and logs copied to briefs/halfline_consolidation_A5_scripts/ (S below); the arXiv sources it read were in S/dl_* and are not copied. -->

# A5 report: classical linearized reflection about the solitonic vacua, compared with Delius–Gandenberger (5.7), and the singular-background limits

I did not edit the repository, run git or spawn subagents. All work is in `/tmp/claude-1000/-home-gustav-Git-affine-toda/b441568d-def3-432f-aa17-ea3eda7016d3/scratchpad/a5/` (called S below). No usage-limit interruption affected any computation; S was empty apart from the downloads when I resumed.

**Sources, read in the TeX downloaded from arXiv into S:**
- DG, hep-th/9904002: `S/dl_dg/final2.tex`, eqs. (5.4)–(5.7), (3.16), Thm. 4.1, block definition (2.x).
- Bowcock, hep-th/9609233: `S/dl_bow/bs.tex`, sec. 4, eqs. (4.1)–(4.22). I counted the equation numbers: (4.14) is the K formula, (4.21) the bracket form of the A-ratio, (4.22) the remaining factor. These match the book's citations.
- BP, hep-th/9909174: `S/dl_bp/bb1.tex`, sec. 4.4.1 (exceptional solutions), sec. 5 (classical scattering, their prescription), Table 2 (`tabscat`).

**Method (all my own code, no worker or Opus code).**
- The exact Hirota method of Bowcock sec. 4 and BP sec. 5. A linear wave of species c is one more Hirota soliton with E = e^{σx−iΩt}, added to the static background solitons to first order.
- Interaction coefficient: A(a,b) = (cosh θ − cos π(a−b)/h)/(cosh θ − cos π(a+b)/h).
- The fluctuation is ε = −Σ α_j δτ_j/τ_j, and it is differentiated exactly (finite sums of exponentials).
- The linearized boundary condition ε'(0) + Hε(0) = 0, with H = ½ Σ C_j α_jα_jᵀ e^{ψ_j(0)/2}, gives an n×n linear system.
- K is taken in the Fourier polarization basis ρ_c = (ω^{jc}), which is the particle basis.
- There is no ODE, so no conditioning problem from the wall layer.
- The book code used:
  - `boundary_vacua.py` (`pair_taus`, `soliton_charge`, `as_weight`, `pattern`, `book_roots`) for the backgrounds, their patterns and charges;
  - `particle_reflection.py` (`blk`, `K_uniform`, `M_sol`, obtained by exec-split before its checks) for DG (5.5) at B = 0;
  - `boundary_common.report`.
- Validation:
  - every dressed wave solves the linearized field equation, with residual 3e-30 at 30 digits and 1e-50 to 1e-90 at higher precision;
  - the flat case reproduces `eq-linear-reflection` to 1e-30;
  - the flux-normalized K is unitary throughout.

## 1. Verdicts

**Claim 6 / task 1 (DG (5.7) against the linearization for every regular vacuum): VERIFIED.**
- Coverage: 56 backgrounds, i.e. every regular χ-interval (2 χ per interval) of every pair species a < h/2 of a₂–a₅, plus the self-conjugate family τ_j = 1 + (−1)^j d e^{2x} (C ≡ −1, b = h/2) for a₃ and a₅.
- Frequencies: Ω = m₁ cosh 0.6, Ω = m₂ cosh 0.6, and Ω = 2.3 (all channels open).
- Results:
  - the τ_j equal `boundary_vacua.pair_taus` to 3e-15;
  - exactly one sign pattern is satisfied, the pattern of the soliton charge;
  - the class b from DG (4.8)–(4.9) (`prp-dg-weights`, converted from book to DG labels) is b = a;
  - K is diagonal (largest off-diagonal entry 5e-28);
  - K_cc equals DG K_c^{(b)}(B→0) in every open channel to 5e-28, at 30 digits.
- So the classical K depends only on the class, not on χ or the interval. This is Bowcock's χ-independence and confirms (4.14), (4.21)–(4.22).
- Also checked algebraically:
  - Bowcock's form −1/((b)(h−b)) · (h/2+a−b)(h/2−a+b)/((h/2+a+b)(h/2−a−b)) equals the DG (5.7) four-block for n = 2..6;
  - all diagonal entries of BP Table 2 for the regular and flat rows equal DG's amplitude of the class.
- Remaining open: the pair energies are checked in `boundary_vacua` only up to a₄. My script checks regularity and pattern for a₅ (a = 1, 2) but not the a₅ energy.

**Claim 5 (a₃ alternating, K = DG b = 2): VERIFIED, as a limit d → 1, and robust.**
- With the Hessian of (+,−,+,−), the off-diagonal entries are 4.97(1−d) → 0, and K₁₁ → −(1)(3) = DG b=2 (= C ≡ −1 value).
- The same limit holds along non-solution paths with an admixed species-1 pair of amplitude (1−d) or (1−d)² (`t6.log`).
- For the ∏C = −1 pattern (+,−,−,−) the limit is the same C ≡ −1 value.
- a₅ middle soliton (+,−,+,−,−,−), class b = 3: off-diagonal 1.42(1−d), K₁₁ → DG b=3, also robust under admixture. This confirms the Opus figure of 1.4e-6 at 1−d = 1e-6.

**Claim 10 / Opus item 1 (a₅ ++----, off-diagonal 0.719 from both sides; BP's diagonal value fails): the numbers are VERIFIED, the interpretation is CORRECTED.**
- Background: species-2 pair at χ_end = 5π/3, where τ₂(0) = τ₅(0) = 0. The regular neighbours are +--+-- (χ > χ_end) and -+--+- (χ < χ_end); both have class 2.
- With the ++---- Hessian (class b = 2), the limit along the regular solutions χ_end ± δ is:
  - K₁₁ = 0.15563283 + 0.67715549i, |K_{1,5}| = 0.71919321 at θ = 0.6;
  - converged to 1.5e-8 between δ = 1e-4 and 1e-6 on both sides, at 60–90 digits;
  - my K₅₁ = −0.21094759 + 0.68756089i matches the Opus entry.
- Other entries agree with Opus:
  - ++--+- gives 0.0913 for χ > χ_end and 0.2721 for χ < χ_end;
  - ++-++- gives 0.3218;
  - at Ω = 2.3 the ++---- off-diagonals are 0.0174.
- BP/DG's diagonal value (2)(4)/((1)(5)) = 0.85655 + 0.51606i is not the limit.
- **Correction 1.** "Same from both sides" is not independent evidence. The reflection j → 1−j maps ++---- to itself and χ_end + δ to χ_end − δ, so equality of the two sides is forced by symmetry. For the non-symmetric ++--+- the two sides differ.
- **Correction 2.** The limit depends on the path.
  - Along static non-solution configurations d = (1−e)/(1+s) at χ = χ_end the ++---- limit is |K_{1,5}| = 0.9374372, K₁₁ = −0.0646044 + 0.3421078i. This is converged at 50 and 90 digits for e down to 1e-9.
  - Mixed paths (d shifted by −e, χ by ±e or ±3e) give the same value.
  - Their zeroth-order boundary residual also → 0 (as √e, against δ along the χ path; `t8.py`).
  - Even the regular pattern's own Hessian (+--+--) gives DG exactly only from its own side. From the other side it gives |K_{1,5}| = 0.5167, and along the d path 0.9374.
- So "BP's diagonal entry for --++++ is corrected by computation" is not justified. BP's value comes from a different prescription: their sec. 5 condition that the perturbation of a vanishing τ_j vanish at the wall. I did not implement that prescription, so their singular rows are UNCHECKED under it.
- What can be claimed: ++---- is non-diagonal in every regular-limit approach tried, with a value that depends on the approach.
- Contrast:
  - for the a₂ ∏C = −1 endpoint χ = π (targets −++ and −−−), the χ path and the d path give the same limit, K₁₁ = 0.005840 + 0.995042i, |K₁₂| = 0.09929. This agrees with the worker's a₂ numbers (`t7.log`).
  - the middle-soliton endpoints are robust.

**Claim 7 (a₄ exceptional, non-diagonal; worker's numbers indicative): VERIFIED and made precise.**
- I reconstructed BP's exceptional solution exactly:
  - species 1, 4, 2, 3 with χ₁ = 9π/5, χ₂ = 18π/5;
  - d₁ = 0.257789025012, d₂ = 0.563669884984;
  - τ₂(0) = τ₄(0) = 0 and regular on x < 0;
  - τ_i'/τ_i(0) = m₁ + m₂ − √5 (BP's C) to 8e-61. I fixed d₁ by the book boundary condition alone (τ_j'/τ_j − C_j e^{ψ_j/2} = const); BP's C was not used as input.
- With the (+,+,−,+,−) Hessian, the limit along several families of regular configurations is K₁₁ → 0.3253474 + 0.9340771i, |K_{1,4}| → 0.147137.
  - Families tried: lowering d₂, lowering d₁, either one combined with χ₁ shifts, and moves along the singular curve with τ(0) ∝ e².
  - All of them are non-solutions; no regular solution family approaches this background.
- DG b=2 is 0.3428236 + 0.9393998i. The limit is non-diagonal, consistent with BP Table 2 "not diagonal" for --+-+.

**Claim 8 and the energies:** not in my task; not rechecked.

## 2. Check script (proposed `code/solitonic_linear.py`)

Path: `S/solitonic_linear.py`. It is self-contained and imports only `boundary_common`, `boundary_vacua` (with its output suppressed) and `particle_reflection` (exec-split). Run with `PYTHONPATH=code`. Runtime 2 min 8 s. Output (`S/solitonic_linear.log`):
```
   56 backgrounds (...), 3 frequencies each; max pde = 3.4e-30, max flat = 1.1e-30, max bv = 2.9e-15, max off = 4.8e-28, max diag = 4.8e-28, max unit = 8.7e-28, max res = 3.3e-28
PASS 1. the dressed Hirota waves solve the linearized field equation eps'' = M^2(x) eps - Omega^2 eps (residual < 1e-25)
PASS 2. flat vacuum, C = +-1: K_a = -((a)(h-a))^{-C} (eq-linear-reflection), diagonal, a_2..a_5
PASS 3. continued pair = boundary_vacua.pair_taus; regular on x <= 0; exactly one pattern, that of the soliton charge; DG (4.8)-(4.9) class b = a (b = h/2, C = -1 for the self-conjugate family)
PASS 4. DG (5.7): about every regular continued vacuum of a_2..a_5 the linearized K is diagonal and equals DG K_c^(b) at B -> 0 in every open channel (Omega = m_1 cosh 0.6, m_2 cosh 0.6, 2.3); flux-normalized K unitary
PASS 5. Bow (4.14) with (4.21)-(4.22) equals DG (5.7) (n = 2..6, all a, b); BP Table 2 diagonal entries of the regular rows (-++; ++++, --++; -++++, ---++; ++++++, -++-++, ----++) are DG K_c^(b) at B -> 0
   a_5 (+,+,-,-,-,-), chi -> 5pi/3: K_11 = (0.15563283 + 0.677155488j), |K_15| = 0.719193205 (spread over both sides, delta = 1e-4, 1e-6: 1.5e-8); DG b = 2 = BP Table 2: (0.856551825 + 0.516061015j)
PASS 6a. a_5 (+,+,-,-,-,-) about the species-2 pair endpoint, as the limit along the regular solutions chi = 5pi/3 +- delta: the same from both sides and not diagonal, |K_15| -> 0.71919, K_11 -> 0.15563 + 0.67716i, not BP/DG's (2)(4)/((1)(5))
   path dependence: d-path (non-solutions) |K_15| -> 0.9374372, K_11 -> (-0.06460441 + 0.3421078j); (+,+,-,-,+,-): |K_15| -> 0.091301 (chi > 5pi/3), 0.27207 (chi < 5pi/3); regular pattern's own Hessian: |K_15| = 1.9e-31 on its side, 0.51675 from the other side
PASS 6b. the limit depends on the approach: non-solution path d -> 1/(1+s) gives |K_15| -> 0.93744 for (+,+,-,-,-,-); for (+,+,-,-,+,-) the two sides give 0.0913 and 0.2721; the regular pattern's own Hessian gives DG only from its own side
   a_3 +-+-: off-diagonal 5.0e-5, 5.0e-7 at 1-d = 1e-4, 1e-6; a_5 +-+---: off-diagonal 0.00014, 1.4e-6 at 1-d = 1e-4, 1e-6
PASS 7. middle-soliton endpoints d -> 1 (a_3 (+,-,+,-), a_5 (+,-,+,-,-,-)): K -> diagonal DG K_1^(h/2) = -(1)(h-1), off-diagonal O(1-d), also along a path with a species-1 admixture
   a_4 exceptional: d_1 = 0.257789025012, d_2 = 0.563669884984; K_11 -> (0.32534737 + 0.93407705j), |K_14| -> 0.1471365 (two paths at 1e-7: spread 1.4e-7); DG b = 2: (0.34282358 + 0.9393998j)
   a_4 exceptional: finite part of the boundary condition at x = 0 (evaluated at x = -1e-30): 9.5e-16; tau_i'/tau_i(0) - (m_1 + m_2 - sqrt 5) = 7.8e-61
PASS 8. a_4 exceptional background (BP): tau_2(0) = tau_4(0) = 0, tau_i'/tau_i(0) = m_1 + m_2 - sqrt 5 (BP's C), boundary condition of (+,+,-,+,-); the linearized K for (+,+,-,+,-), approached along two families of regular configurations, has one non-diagonal limit |K_14| = 0.1471, not DG's diagonal b = 2 value
```
Suggested README row: `solitonic_linear.py` | `prp-solitonic-linear` (new), `tbl-boundary-types`, `sec-boundary-missing-charges`, `prp-solitonic-amplitude` | exact (Hirota) linearized reflection about the continued vacua against DG (5.7); limits about boundary-singular backgrounds | about 2 min.

## 3. Proposed book text

**(a) ch. 31, `sec-classical-reflection-factors`.** Replace the paragraph "For the solitonic boundary conditions the vacuum is the continued pair … we have not checked either computation." with:

> For the solitonic boundary conditions the vacuum is the continued pair (-@eq-vacuum-pair). Bowcock solved the linearized problem exactly: the incoming and outgoing waves are two further Hirota solitons of complex rapidity added to the static solitons of the vacuum [@bowcock1998, sec. 4]. The quantum amplitudes themselves are the subject of @sec-solitonic-boundaries.
>
> ::: {#prp-solitonic-linear}
> **[Numerical]** About every continued vacuum of @prp-solitonic-vacua (every species $a<\frac h2$, every interval of $\chi$) and of @prp-self-conjugate-vacuum, for $a_2^{(1)},\dots,a_5^{(1)}$, the linearized reflection matrix is diagonal in the particles, independent of $\chi$, and
> $$K_c^{\rm cl}(\theta)=-\frac{\big(c+b-\frac h2\big)\big(-c-b-\frac h2\big)\big(c-b+\frac h2\big)\big(b-c+\frac h2\big)}{(c)(h-c)},$$
> where $b$ is the class of the pattern given by @prp-dg-weights ($b=a$ for the pair of species $a$, $b=\frac h2$ for $C\equiv-1$). This is Bowcock's (4.14) with (4.21)–(4.22), and the limit $B\to0$ of Delius and Gandenberger's solitonic amplitudes, DG (5.7).
> :::
>
> Checked in every open channel at three frequencies, to 27 digits, in `code/solitonic_linear.py`, which also checks that the diagonal entries of Bowcock and Perkins's Table 2 for these vacua are the same functions [@bowcock_perkins2003].

**(b) `tbl-boundary-types`, solitonic row, column "classical bound states".** Replace "Bowcock's factors (not checked here)" with: "diagonal, DG (5.7) at $B\to0$ (@prp-solitonic-linear); e.g. $a_2^{(1)}$: one pole, at $\theta=i\pi/6$".
- The a₂ pole is from `Prod.poles` of K_uniform·M_sol at B = 0: physical-strip pole at p = 1/2.
- For the other classes the pole lists are in my session output; they can be generated the same way.
- The ∏C = −1 row can stay "not settled (depends on the linearization at the wall)", optionally adding "(@sec-boundary-missing-charges)".

**(c) `tbl-missing-vacua`.** It contains no "not checked here". No change is needed. Optionally add to the caption: "About every counted vacuum the linearized reflection factor is Delius and Gandenberger's at $B\to0$ (@prp-solitonic-linear)."

**(d) `prp-solitonic-vacua`, last sentence.** Append: "; regularity and the pattern also for $a_5^{(1)}$ ($a=1,2$) in `code/solitonic_linear.py`".

**(e) ch. 31, `sec-boundary-missing-charges`.** Replace the last sentence "The classical reflection factors about these singular backgrounds depend on how the linearization at a vanishing $\tau_j$ is defined, so their comparison with Delius and Gandenberger's amplitudes is not settled." with:

> About a field that diverges at the wall the linearization is not defined as it stands, and one has to say how the singular background is approached. For the endpoint $d\to1$ of @prp-self-conjugate-vacuum the answer does not depend on it: with the boundary potential of the alternating pattern the linearized reflection matrix tends to the $C\equiv-1$ value $-(a)(h-a)$, Delius and Gandenberger's $b=2$ amplitude at $B\to0$, along the family itself and along configurations with an admixed species-1 pair, and the off-diagonal entries vanish like $1-d$. The same holds for the $a_5^{(1)}$ middle-soliton patterns such as $(+,-,+,-,-,-)$, $b=3$. At the singular endpoints of stationary pairs the limit depends on the approach. For the $a_5^{(1)}$ pattern $(+,+,-,-,-,-)$, of class $b=2$ and listed with a diagonal factor by Bowcock and Perkins [@bowcock_perkins2003, Table 2], the limit along the regular solutions $\chi\to\frac{5\pi}3$ is not diagonal, $|K_{1,h-1}|=0.719$ at $\theta=0.6$, while along static configurations with $d\to\frac1{1+s}$ it is $0.937$; for $(+,+,-,-,+,-)$ the two sides of $\chi=\frac{5\pi}3$ already differ. For the exceptional $a_4^{(1)}$ background of $(+,+,-,+,-)$ the limits we tried agree and are not diagonal, $|K_{1,4}|=0.147$. Bowcock and Perkins require instead that the perturbation of a vanishing $\tau_j$ vanish at the wall [@bowcock_perkins2003, sec. 5]. So, except at the middle-soliton endpoints, the classical reflection factor about these backgrounds is not fixed by the background and the boundary potential, and its comparison with Delius and Gandenberger's amplitudes depends on a prescription (`code/solitonic_linear.py`).

**(f) ch. 35.**
- Line 71: "(checked by hand)" → "(@prp-solitonic-linear)".
- Line 179 (`prp-solitonic-amplitude` text): add "@prp-solitonic-linear checks it against the linearization."

## 4. Statements in the book or Appendix D that these results contradict

- **ch. 31, l. 229** "we have not checked either computation", and **`tbl-boundary-types`** "Bowcock's factors (not checked here)": superseded by check 4.
- **ch. 35, l. 179** "For the alternating patterns of $a_3^{(1)}$ and the other missing patterns of @tbl-missing-vacua no classical vacuum is known at all". This is inconsistent with ch. 31's "Boundary-singular vacua" paragraph and the a₃ endpoint discussion. It should say that only boundary-singular energy infima are known. Also **l. 314** (open problem 6) is partly answered.
- **Integration table / Opus item 1** "BP's diagonal entry for `--++++` is CORRECTED by computation", and the count "6/15": not supported. The value is prescription- and path-dependent. **Do not add BP Table 2 to Appendix D.**
- **Appendix D:** nothing in my results contradicts it.
- **Minor note on DG (5.7):** DG print (5.7) as the limit of $\kappa^{(a)}_b$ with $\prod_{k=1}^b$, i.e. with a and b exchanged relative to (5.4). The exchange is consistent and the four-block is symmetric in a ↔ b, so this is not an error.
- **Minor note on the Opus scripts:** their "unitarity" column at Ω = 2.3 (≈3e-3) is an artefact of not flux-normalizing K. With √(p'/p) normalization it is unitary to working precision.

## 5. Files (all in S = `/tmp/claude-1000/-home-gustav-Git-affine-toda/b441568d-def3-432f-aa17-ea3eda7016d3/scratchpad/a5/`; copied to `briefs/halfline_consolidation_A5_scripts/`)

- Deliverable: `S/solitonic_linear.py`, with output in `S/solitonic_linear.log`.
- Library used for the exploratory runs: `S/a5lib.py`.
- `S/t1.py`: flat and a₂ pair tests (printed only).
- `S/t2.py` → `S/t2.log`: task-1 sweep, a₂–a₅.
- `S/t3.py` → `S/t3.log`: a₅ pair endpoint, a₃/a₅ middle soliton, a₂ ∏C = −1; 90 digits.
- `S/t4a.py`, `S/t4b.py` → `S/t4b.log`, `S/t4c.py` → `S/t4c.log`, `S/t4d.py` → `S/t4d.log`, `S/t4e.py`: a₄ exceptional background and approach paths.
- `S/t5.py` → `S/t5.log` (++----, ++--+-) and `S/t5b.log` (run with argument `reg`: regular patterns' own Hessians); `S/t5c.py`: precision check of the d path.
- `S/t6.py` → `S/t6.log`: a₃ middle soliton with admixture paths.
- `S/t7.py` → `S/t7.log`: a₅ middle-soliton admixture; a₂ endpoint paths.
- `S/t8.py`: boundary-condition residuals along the paths.
- Sources: `S/dl_dg/final2.tex` (DG), `S/dl_bow/bs.tex` (Bowcock 1996), `S/dl_bp/bb1.tex` (BP).
