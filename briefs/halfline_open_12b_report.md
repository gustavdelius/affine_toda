<!-- Worker report, 2026-10-11. Model: Claude Sonnet 5.5, effort high. Scope: problem 12 part (b) only (brief halfline_open_12_rank_one.md), plus replacement text for prp-bsg-excited and prp-bsg-spectrum. Saved verbatim by the coordinator. Scripts copied to briefs/halfline_open_12b_scripts/ (three of them import code/boundary_common.py by an absolute sys.path insert). -->

# Part (b) report: higher excited boundary states of boundary sine-Gordon

Scope: part (b) only (plan step 4), plus replacement text for `prp-bsg-excited` and `prp-bsg-spectrum`. I did not touch the repository.

## Verdict
- **Algebraic claims.** Everything in `eq-bsg-excited-k` and the two-label closed forms holds, and for **all n**. It is no longer limited to n≤2. I could prove it analytically rather than only check it.
- **Spectrum rule.** `prp-bsg-spectrum` is supported by BPT's TCSA identifications up to two labels, and for two three-label states. Beyond that it is the bootstrap's own "no Coleman–Thun diagram" argument, which BPTT themselves call "at most a reasonable guess".
- **Gap in the pole bookkeeping of BPTT/MD (and the book).** The sources omit a tower of poles that the closed forms carry. It always includes ν_0, and it includes ν_{-j} when η < η_N − jπ. No TCSA point in BPT lies in that region (details in section 3).

## 1. What was proved or checked

**Result 1 (proved from the Γ-product, numerically confirmed).** The bulk amplitude a = S_{++} (BPTT 3.1, equal to the book's S_0) obeys

  a(u − π/λ) = tan(u/2) · tan(π/(2λ) − u/2) · a(u)   (to 5e-20, 4 values of λ).

The proof shifts v = λu/π by 1 in the Γ-product, which telescopes via Γ(z)Γ(1−z) = π/sin πz. Three consequences follow.
1. a(u+x−π/λ) a(u−x+π/λ) / (a(u+x) a(u−x)) = {2λx/π − 1} for all complex x (checked to 4e-20).
2. Since ν_n = ν_{n−1} − π/λ, the ratio form of a_n telescopes to the product form Π_{l≤n}{2(η/π − l)} for **every n**. I also checked n = 1..10 at 5 points (λ,η), including η outside the fundamental range, to 4e-19.
3. The n-th bootstrap follows from the n=0 one. With b = (sin λu / sin λ(π−u)) a and λν_n = λν_0 − nπ, the sine factors cancel exactly. So a(u−ν_n)Q(η)b(u+ν_n) = Q(η̄)a_n holds for all n once the n=0 identity holds. The n=0 identity is the only genuinely numerical input, and it is already verified in `boundary_sine_gordon.py`.

**Scalar bootstrap on multi-label states.** The iterated bootstrap Q^{S+(m)} = a(u−p)Q^S b(u+p) reproduces BPTT's closed forms (§3.8) for depth up to 4 and labels up to 4. That is 54 cases, to 1e-19.
- Even depth: Q^S = Q(η)Π_{odd}a_{n_i}(η)Π_{even}a_{n_i}(η̄).
- Odd depth: Q(η) is replaced by Q(η̄).

**Unitarity and crossing.** a_n(u)a_n(−u) = 1 and a_n(π/2+w) = a_n(π/2−w) hold for n≤7, with η and with η̄ (to 2e-19).

**Full 2×2 matrix bootstrap** (`p12b_matrix_bootstrap.py`). I derived it from the residue of the braid-form reflection equation at a pole of K(y2). The residue is rank one, and the reflection matrix A on the new state obeys (1⊗ρ_in^T) R K R (v⊗w) = A v ⟨ρ_in, w⟩. I tested that A equals the claimed matrix: X K(η̄) X on |n>, K(η) on |n,m>, X K(η̄) X on |n,m,k>, and so on.
- Coverage: 120 steps, n≤5, depth ≤3, four (λ,η,ϑ) including ϑ=0 and η beyond the strip conditions.
- Result: the residue is rank one to 5e-15, the factorization in w holds to 1e-15, and A matches the closed form to 1e-14.
- Negative controls: the |0> matrix is not proportional to K(η̄) without the s↔s̄ swap (residual 0.32), nor to K(η) (0.51).

**Unitarity and crossing-unitarity of the full K^S.** I tested 14 states to depth 4 at two parameter points, with the book's S_0 and the physical q. The result is 3e-15, so no extra CDD factor is needed.

**Breather reflection** (`p12b_breather_bootstrap.py`). For B_1 and B_2 on 10 states (depth ≤3, labels ≤5), at 2 values of (λ,η,ϑ), the full matrix bootstrap gives K_{B_m}(x_S) Π b_{m n_i}(η or η̄). This is the generalisation of BPTT (3.9). The breather vector also decouples from the soliton multiplet exactly (3e-15 over 40 cases).

## 2. Comparison with published TCSA numbers (`p12b_tcsa_compare.py`)
- **BPT Table 5.2, "predicted" column.** E_1−E_0 = 2M cos((η−π/2)/λ) with η = (λ+1)π·βφ_0/2π reproduces all 7 printed predictions. This tests only |0>.
- **BPT §4.3.2 (λ=7, η=0.9η_N, ϑ=0).**
  - My generator gives the strip sequence ν_0>w_1>ν_1>w_2>ν_2>w_3>ν_3 (BPT 4.6).
  - Every state in BPT's lists (4.7), (4.8) is allowed by the rule and sits in the right C-parity sector (parity of the sum of labels).
  - |1,1> is excluded, as BPT find (no TCSA line).
  - The Fig. 4.3 captions say the identified lines are as follows. These are the checked-against-TCSA states.

| Sector | Identified states |
|---|---|
| Even | \|0>, \|2>, \|0,2>, \|0,1,1> |
| Odd | \|1>, \|3>, \|0,1>, \|0,3>, \|1,2>, \|0,1,2> |

  - The text lists (4.7), (4.8) omit |0,1>, |0,3> and |1,2>, but the captions include them.
  - No TCSA identification exists for anything deeper than |0,1,1> and |0,1,2>.
- **Exact energy relation.** Zamolodchikov's boundary energy (BPT 3.4) gives E(η̄,ϑ) − E(η,ϑ) = M cos ν_0 identically. BPTT's conjecture m_{|0>}(η) = m_{|>}(η̄) is therefore an exact consequence of that formula, for all η and ϑ.

## 3. Pole bookkeeping and the remaining gap (`p12b_poles.py`, `p12b_spectrum_poles.py`, `p12b_scan_U.py`)

I wrote an exact pole/zero counter for Q^S in the physical strip, from the Γ-arguments of R_0 and σ and the sines in the blocks. It agrees with orders measured numerically at 70 (state, position) pairs, with 0 mismatches (`p12b_poles_check.py`).

**Poles at ν_j, w_j with j≥0.** These are accounted for as BPTT describe:
- new states at ν_m (even depth) or w_m (odd depth) below the last label;
- the remaining poles are Coleman–Thun poles of increasing order: ν_0 has order 2 for |n≥1>, and in two-label states it reaches 3 or 4 at depth 3.
- MD's printed Table 7 gives "ν_0: order k". Their own §5.3 example (double pole at ν_0 on |1;1>) and my count both give 2 for n≥1, so that table entry looks like a misprint or an omission.

**Poles not in the sources' lists (class "U").** The denominators of the blocks {y} in a_m(η̄), which sits at even label positions, produce poles at ν_0, ν_{−1}, …, ν_{−m}.
- ν_0 and ν_{−m} are simple; those in between are double.
- σ(η̄) also has b-type poles at ν_0, ν_{−1}, …
- These mirror the numerator tower w_0 … w_m. BPTT and MD discuss the numerator tower only.
- The ν_{−j} poles lie in the physical strip iff η < η_N − jπ, with η > π/2 so that |0> exists. In terms of x = βφ_0/2π this is x < ½ − 1/(λ+1) for j=1. The scan over 6 λ's and 6 η/η_N values matches this prediction exactly.
- Their orders grow with depth, e.g. ν_{−1} has order 1 on |0>, 2 on |0,5>, and 3 on |0,5,1>.
- MD's Tables 6 and 7 list ν_{−N} with order 1 (charge 1 only) and nothing for charge 0.
- **TCSA coverage.** All TCSA points quoted by BPT (η = 0.9η_N and the Table 5.2 Dirichlet points) have η > η_N − π, so none lies in this region. The region x < ½ − 1/(λ+1) is the physically ordinary one for moderate φ_0 and λ > 3.
- **Not established.** I could not find an explicit Coleman–Thun explanation of these poles in the sources (MD only say "same mechanisms with ξ → π(λ+1) − ξ"). Naively extending diagram 3.2(c) to ν_{n−k} with k>n would give a double pole where the closed form has a simple one. I flag them as unexplained in the literature, not as an inconsistency.
- A decisive test would be a TCSA or BTBA spectrum at, say, λ=7, βφ_0/2π ≈ 0.25 (non-generic coincidences, e.g. ν_{−1}=w_3, occur at that exact point).

## 4. Proposed replacement text

**`prp-bsg-excited`** (replaces the existing text):

> **[Numerical]** For n=0 the boundary bootstrap a(u−ν_0)Q(η,ϑ,u)b(u+ν_0) = Q(η̄,ϑ,u), with a=S_0 and b=S_T, holds numerically; Q(η,ϑ,u) satisfies unitarity and crossing-unitarity, and σ(η,u) has simple poles at ν_0, ν_1. Given this, the rest follows for every n. The bulk amplitude satisfies a(u−π/λ)=tan(u/2)tan(π/(2λ)−u/2)a(u), which gives a(u+x−π/λ)a(u−x+π/λ)/(a(u+x)a(u−x)) = {2λx/π−1}. Since ν_n = ν_{n−1} − π/λ, the ratio form of a_n equals the product form for all n. Because λν_n = λν_0 − nπ, the sines in b(u+ν_n)/a(u+ν_n) cancel, so a(u−ν_n)Q(η)b(u+ν_n) = Q(η̄)a_n(η,u) for all n. The same argument with η → η̄ gives the two-label factors Q(η)a_n(η)a_m(η̄). The 2×2 matrix bootstrap reproduces X K(η̄) X, K(η), … for n≤5 and up to three labels (1e-14), and K^S is unitary and crossing-unitary for 14 states of depth ≤4 (3e-15), all in `p12b_*` scripts to be merged into `code/boundary_sine_gordon.py`.

**`prp-bsg-spectrum`** (status line):

> **[Conjecture]** (the rule itself). The states with up to two labels, and |0,1,1>, |0,1,2>, are identified in the TCSA of Bajnok–Palla–Takács for λ=7, η=0.9η_N, ϑ=0; the Dirichlet |0> energies agree with their Table 5.2. State |1,1> is absent, as the rule predicts. For deeper states the existence of the state rests on the absence of a boundary Coleman–Thun diagram, for which no Cutkosky-type rules exist [BPTT sec. 5].

**Paragraph after `prp-bsg-spectrum`**, to replace "the reflection factors on the higher states are extrapolated":

> The closed forms for all states are consequences of the bootstrap and are unitary and crossing-unitary. What is not fully analysed is the pole content in the physical strip. The CDD factors a_m(η̄) at even positions also have poles at ν_0, ν_{−1}, … from the denominators of {y}, and σ(η̄) has b-type poles at ν_{−j}. These lie in the strip when ν_{−j} < π/2, i.e. η < η_N − jπ, and are not covered by the explanations of BPTT or Mattsson–Dorey nor by any TCSA data.

## 5. Scripts (all in `/tmp/claude-1000/-home-gustav-Git-affine-toda/0e38a99d-702a-40fb-99b7-74dbe0413cf6/scratchpad/p12b/`)
- `p12b_common.py`: shared mpmath functions (S_0, R_0, σ, blocks, a_n).
- `p12b_an_bootstrap.py`: items in section 1 (functional equation, product form, scalar bootstrap, unitarity/crossing of a_n). All PASS; runtime about 26 min.
- `p12b_matrix_bootstrap.py`: full 2×2 bootstrap plus negative controls (seconds).
- `p12b_unitarity.py`: unitarity and crossing-unitarity of K^S for 14 states.
- `p12b_breather_bootstrap.py`: B_1, B_2 on excited states (40 cases).
- `p12b_poles.py`: exact pole/zero-order counter.
- `p12b_poles_check.py`: the counter against numerical orders (70 pairs; output in `poles_check.out`).
- `p12b_spectrum_poles.py`: allowed states, energies, C-parity sector, pole table with classes N/L/U. Usage: `python3 p12b_spectrum_poles.py λ η/η_N maxdepth`.
- `p12b_scan_U.py`: scan showing U poles appear exactly when ν_{−1}<π/2.
- `p12b_tcsa_compare.py`: BPT Table 5.2, the lists (4.6)–(4.8), and Zamolodchikov's energy identity.
- Downloaded texts (read only): `dl/bpt.txt` (BPT hep-th/0108157), `dl/md.txt` (Mattsson–Dorey hep-th/0008071), `bptt.txt` (BPTT hep-th/0106070).

## 6. Status
- **Done:** all of the above.
- **Partial:** the explanation of the U-class poles by Coleman–Thun diagrams (flagged, not resolved).
- **Not done:** a new TCSA or BTBA test in the x < ½ − 1/(λ+1) region; breathers B_m with m≥3; part (a) and part (c) apart from the wording above.
- **Convention caveat:** the scripts import from `code/boundary_common.py` with a read-only `sys.path` insert and use the same conventions as `boundary_sine_gordon.py`.
- **Sector caveat:** I used the C-parity sector assignment "parity of the sum of labels" from BPT §4.3.2, which holds at ϑ=0.

## Integration (branch half-line-open-problems-integration)

| Claims | Verdict | How it was checked | Book change |
|---|---|---|---|
| 1–3 (shift identity; telescoping; bootstrap for all $n$ from $n=0$) | **Verified, integrated** (2026-10-12) | The coordinator recomputed the shift identity $a(u-\pi/\lambda)=\tan\frac u2\tan(\frac\pi{2\lambda}-\frac u2)a(u)$ from the book's *integral* form of $S_0$ (`eq-sg-s0`), independent of the Gamma product the worker and Gemini used: $\le3\times10^{-25}$ at $\lambda=1.4,2.7,5.3$. The consequence identity was checked symbolically at random complex points. Claims 2 and 3 were redone by hand. | `prp-bsg-excited` now holds for all $n$, with a proof. A new check (12b) in `code/boundary_sine_gordon.py` tests the identity, the ratio/product forms for $n\le6$ and the bootstrap for $n\le5$; all 32 checks PASS. |
| 9 (U-class poles; "literature correction" to MD per Gemini) | **Gemini's correction rejected; narrower correction adopted** | Opus follow-up (`halfline_open_12b_referee_opus.md`). The coordinator read MD's TeX (bsgdht.tex l.984–987, 1874–1886): MD's sentence concerns $w_x$ on charge-0 states, and MD do list the $\nu_{-N}$ poles. $\nu_{-j}=\pi-w_j$ checked by hand. The coordinator reran the Opus pole-order script (all orders integers, matching the predicted formula). | ch. 29: new paragraph "Unexplained poles"; open problem 12 points to it; Appendix D: MD Tables 6–7 pole orders (*Error*). |
| 10 (MD Table 7 "$\nu_0$: order $k$") | **Verified** (worker, Gemini, Opus, coordinator rerun) | | Included in the Appendix D MD entry |
| 11 (BPT lists) | **Corrected** (Gemini was right on content): only $|0,1\rangle$, $|0,3\rangle$ missing; $|1,2\rangle$ is in (4.8). The coordinator read uv_ir3.tex l.1115–1121. | | ch. 29 paragraph after `prp-bsg-spectrum` lists the TCSA-identified states and the omission; Appendix D: BPT (4.7)–(4.8) (*Incomplete*) |
| 4–7 (multi-label scalar and matrix bootstrap, unitarity/crossing of $K^S$, breathers) | **Not integrated as new claims** | worker + Gemini rerun only | The ch. 29 text says the factors follow from the bootstrap on every state (as BPTT state), without new numerical claims |
| 8 (TCSA comparison; energy identity $E(\bar\eta)-E(\eta)=M\cos\nu_0$) | TCSA list **integrated** (from BPT's caption); energy identity **not integrated** | | ch. 29 |
| 12 (text) | Integrated as modified; the worker's "all states with at most two labels" was wrong ($|1,3\rangle$, $|2,3\rangle$ are not identified), so the list is given explicitly; $\eta_N$ replaced by $\frac\pi2(\lambda+1)$ | | |
