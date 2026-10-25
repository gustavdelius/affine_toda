# Brief: referee the corrected sign of φ_0 in the Dirichlet relation of boundary sine-Gordon theory

You are an adversarial mathematical-physics referee. You are working in the repository of a Quarto book on affine Toda field theory (worktree `.claude/worktrees/eager-bun-c1cfa6`, branch `claude/eager-bun-c1cfa6`). Review commit `ed0585a`.

**Rules.**
- Other agents may be editing files in this worktree at the same time.
- Do not edit any existing file.
- Do not run git commands that change the index or the working tree (no commit, checkout, stash, reset, rebase or restore).
- Do not run `quarto render`.
- Do not spawn subagents.
- Write scratch scripts only under `referee_scratch/bsg_dirichlet_sign/`.
- Write the report to `briefs/referee_bsg_dirichlet_sign_report.md`.

## Context

The book claims that Ghoshal and Zamolodchikov's Dirichlet relation ξ = 4πφ_0/β (hep-th/9306002, eq. (5.27)) has the wrong sign for their own half-line x ≤ 0. In the book's parameters it now reads η = −4πφ_0/β_SG (`eq-bsg-dirichlet`). The UV–IR relation `cnj-bsg-uv-ir` inherits the sign through its Dirichlet limit, so its second equation now has the opposite sign to the printed one (Bajnok–Palla–Takács–Tóth hep-th/0106070 sec. 4; Bajnok–Palla–Takács hep-th/0108157 eq. (3.2)). With the printed sign the book had found a "node crossing": the UV–IR relation paired the coideal parameter ε̂_1 with the boundary parameter C_0, while the non-local charges pair ε̂_1 with C_1. With the corrected sign both pair each node with itself, and the crossing is gone. This is a correction to a much-cited paper, so it needs an independent check before it stays.

What changed (commit `ed0585a`, on top of PR #59, branch `claude/nostalgic-perlman-a1e788`, which fixed the chirality of the non-local charges):
- ch. 30 (`part7-half-line/30-boundary-sine-gordon.qmd`): the "Dirichlet" paragraph and `eq-bsg-dirichlet`, the "General boundaries" paragraph and `cnj-bsg-uv-ir`, the "Node labels" paragraphs of `sec-bsg-parameters`, the Dirichlet sentence of `prp-bsg-spectrum`, the sign of F in the bullets after `cnj-shg-e`, and a new item in the chapter's corrections to the sources;
- ch. 33 (`part7-half-line/33-boundary-quantum-groups.qmd`): one sentence of `sec-boundary-charges` and `prp-dm-first-order`;
- ch. 38 and ch. 43 (`part8-defects/`): the sign of E_D in the bullet after `cnj-defect-parameter`, and the proof of `prp-defect-fusion-k` (4);
- Appendix D (`appendices/D-corrections.qmd`): new items for Ghoshal–Zamolodchikov (5.27), Mattsson–Dorey §2.2 and Bajnok–Palla–Takács–Tóth §4;
- the lecture notes `lecture-notes/defects/09-folding-fusing.qmd`, `lecture-notes/half-line/04-coideals-and-sine-gordon.qmd`, `05-a2-intertwiners-and-q-onsager.qmd`;
- `AGENTS.md`, the trap "Half-line notation".

Check scripts:
- `code/bsg_dirichlet_sign.py` (new; 5 checks, about 10 s): static Dirichlet states, the pole ν_0 in P^±, the semiclassical bound-state energy, a classical leapfrog simulation of soliton and antisoliton reflection, and the convention e^{iS};
- updated: `code/boundary_sine_gordon.py` (check 8, about a minute), `code/boundary_coideal.py` (check 9, about 50 s), `code/nonlocal_charge_chirality.py` (17 checks, 3 s), `code/shg_boundary_breathers.py` (check 8, seconds), `code/lectures/defects_l9_fusion.py` (item (h)).

The drafter's own flags for you:
- The sign of F in the sinh-Gordon relation is a convention (the sinh-Gordon factor depends on F² only), so its agreement with Corrigan–Taormina is not evidence.
- The current-side argument (item 3) rests on PR #59's chirality labels. Codex's report on `briefs/referee_defects_l9.md` (item 5, in the PR #58 worktree) agreed with them, but they are recent.
- The book says the sign is invisible to the truncated conformal space tests of Bajnok–Palla–Takács because those are invariant under charge conjugation. This was argued from their set-up, not computed.

## Notation

- Bulk on x ≤ 0, boundary at x = 0. Action S = ∫dt [∫_{−∞}^0 L dx − 𝓑(φ(0,t))], L = ½(∂φ)² − (m_0²/β_SG²)(1 − cos β_SG φ), 𝓑 = −M_0 cos(β_SG(φ − φ_0)/2), so ∂_xφ(0) = −(β_SG M_0/2) sin(β_SG(φ(0) − φ_0)/2) (`eq-bsg-boundary`).
- Topological charge T = (β_SG/2π)(φ(+∞) − φ(−∞)); the soliton s has T = +1. On the half-line the Dirichlet boundary conserves (β_SG/2π)(φ(0) − φ(−∞)).
- λ = 8π/β_SG² − 1; the physical q = −e^{−iπλ}; coproduct Δe = e⊗1 + k⊗e.
- Reflection matrix K(θ): V(θ) → V(−θ) in the basis (s, s̄), `eq-bsg-k-full`, with P^+ = K_ss, P^− = K_s̄s̄, u = −iθ. Parameters (η, ϑ) by `eq-bsg-eta-theta`: ε̂_1 = −cos(η − iϑ)/sin πλ, ε̂_0 = −cos(η + iϑ)/sin πλ. At ϑ → ∞, P^± ∝ cos(λu ± η). Ghoshal–Zamolodchikov's ξ is the book's ξ_B; at the Dirichlet point ξ_B = η.
- Coideal generators b_j = e_j + q^{−1}f_jk_j + ε̂_jk_j, with e_1 = yE, e_0 = yF, k_1 = diag(q, q^{−1}), y = e^{λθ} (`eq-sg-coideal-rep`).
- Boundary parameters C_{0,1} = (β_SG²M_0/4m_0) e^{±iβ_SGφ_0/2} (`eq-bsg-c-dictionary`); C_j multiplies e^{iβα_j·φ/2} with α_1 = −α_0 = √2.
- x^± = t ± x, E + P has spin +1; the e_j are built from right-movers φ_R(x^−) (`sec-nonlocal-charges`, after PR #59).
- Beware: Saleur–Skorik–Warner put the bulk on x ≥ 0; Mattsson–Dorey, Bajnok et al. and Ghoshal–Zamolodchikov on x ≤ 0. Parity reverses T.

## Read first

- `sec-bsg-parameters` and `sec-bsg-k-matrix` in ch. 30; `sec-boundary-charges`, `sec-sg-k-from-symmetry` and `prp-dm-first-order` in ch. 33; `sec-nonlocal-charges` and `sec-qg-coproduct` in ch. 21. Find an object with `grep -rn '{#<id>}' part*/ appendices/ outlook/ *.qmd`.
- The book's conventions: `01-introduction.qmd`, `sec-conventions`.
- Sources, by arXiv number: hep-th/9306002 (Ghoshal–Zamolodchikov; §2 for the half-line, (3.16), (3.41)–(3.42), §5 (5.3)–(5.27)), hep-th/9408004 (Saleur–Skorik–Warner, §2 and §3), hep-th/0008071 (Mattsson–Dorey, §2.2), hep-th/0106070 (BPTT, §2 and §4), hep-th/0108157 (BPT, (3.2)–(3.3), §3.5, §4–5), hep-th/0112023 (Delius–MacKay), hep-th/0008237 (Corrigan–Taormina (2.10)), hep-th/0002065 (Chenaghlou–Corrigan (6.1)). On the editor's machine source texts are in `~/referee-copy/papers/*.txt`; these papers are not there yet. If you cannot find them, download the arXiv sources (`https://arxiv.org/e-print/hep-th/<number>`) into `referee_scratch/bsg_dirichlet_sign/src/`.

## Evidence

For each item, do **one independent computation**: a derivation written out, or a script of your own. You may run the book's scripts for comparison, but they are not evidence. A verdict without its own evidence counts as "unclear". `briefs/referee_defects_38_40_report.md` shows the expected depth.

## Items (most important first)

1. **The charge argument for `eq-bsg-dirichlet`.**
   - Read Ghoshal–Zamolodchikov's set-up yourself: is their bulk x ≤ 0 (their (2.8), (3.16)), and is their soliton the particle of charge +1 in (5.3)?
   - For 0 < φ_0 < π/β, find all static solutions with φ(0) = φ_0 and their energies and charges (β/2π)(φ(0) − φ(−∞)). Which state becomes degenerate with the ground state at φ_0 = π/β?
   - Which particle, coming in from x = −∞ and absorbed, takes the ground state to that state?
   - In `eq-bsg-k-full` at ϑ → ∞, which of P^± has the pole ν_0 = η/λ − π/(2λ) of σ(η,u) for η > 0? Is σ(−x,u) = σ(x,u), so that η → −η exchanges P^+ and P^−?
   - Conclude: on x ≤ 0, is η = −4πφ_0/β or +4πφ_0/β? Is Ghoshal–Zamolodchikov's sign right for a bulk on x ≥ 0, as the book says, and does Saleur–Skorik–Warner's remark on the charges of the bound states (end of their §3, "kink sitting in the middle") support that reading?
   - Is there any reading of Ghoshal–Zamolodchikov's "emission of a zero-momentum soliton by the associated boundary state" (the crossed channel of their §3) in which their sign is consistent on x ≤ 0?
2. **The semiclassical check.** The book argues that P^+/P^− → e^{2iη} for λθ ≫ 1, that the amplitude is e^{iS} with ∂S/∂φ_0 = −∫dt ∂_xφ(0,t), and that a classical simulation gives ∫dt [∂_xφ_s(0,t) − ∂_xφ_s̄(0,t)] = 8π/β_SG. From this, 2∂η/∂φ_0 = −8π/β_SG.
   - Derive ∂S/∂φ_0 yourself, with the subtraction of the ground state, and the relation between the phase of the amplitude and S (sign included).
   - Compute the integral independently: your own solver, or the exact three-soliton solutions of Saleur–Skorik–Warner. Use several velocities and boundary values, and both signs of φ_0.
   - Check the convention independently of the sign question. The common phase 2 arg σ(η, −iθ) depends on |η| only. Does its η-derivative at large λ match −(I_s + I_s̄)/8π, the sum of the two integrals?
   - Are the classical time delays of soliton and antisoliton equal, so that time delays alone cannot fix the sign?
3. **The node pairing.**
   - With the corrected `cnj-bsg-uv-ir`, show that cos((η ∓ iϑ)/(λ+1)) = μ_∓ and that ε̂_1 depends on μ_− ∝ C_1 alone.
   - On the current side (`sec-bsg-parameters`, "Node labels"; `sec-boundary-charges`): check that b_1 is built from J_+ = e^{−iaφ_R} and J̄_− = e^{−iaφ_L}, a = 8π/β_SG. Check that both have double poles with e^{iβ_SGφ} and with the boundary term e^{iβ_SGφ/2}, which carries C_1, and that e_1 = yE raises T.
   - Check that Delius–MacKay's Q̂_+ is the book's b_0 (hep-th/0112023, (3.18), (3.27)).
   - Check the expansion in `prp-dm-first-order`: ε̂_{1,0} = ε_* − (λ+1)/(2 sin(πλ/2)) μ_∓ + O(μ²).
   - Is this pairing evidence for the sign independent of items 1 and 2?
4. **`cnj-bsg-uv-ir` with the changed sign.**
   - Is it legitimate to state the conjectured relation with the opposite sign of its second equation, given that the sign enters only through the Dirichlet limit?
   - Check the limits listed after the conjecture: Neumann, Dirichlet (now η = −4πφ_0/β), the point ε̂ = 0, and the classical M_crit.
   - Are the tests of Bajnok–Palla–Takács really blind to the sign? Their truncated conformal space runs perturb one end and keep the other Neumann, or use equal Dirichlet ends. Is the label [Conjecture] still right?
5. **Corrections to the literature.** Check, against the source lines, each new item in Appendix D and in the ch. 30 corrections:
   - Ghoshal–Zamolodchikov (5.27): the quotation and the claim;
   - Mattsson–Dorey §2.2: they print the charge of the first excited state as 1 − βφ_0/2π. Is the book right that on x ≤ 0 it is βφ_0/2π − 1, and that 1 − βφ_0/2π is the value for x ≥ 0?
   - BPTT §4 and BPT (3.2).

   Are the categories right (*Error*, inherited)? Is anything overstated, for example "the spectrum is unchanged"?
6. **Prior art.** Has anyone noticed or used the opposite sign before? Check at least:
   - Skorik–Saleur hep-th/9502011 and Fendley–Saleur hep-th/9402045 (lattice derivations with a boundary field);
   - LeClair–Mussardo–Saleur–Skorik hep-th/9503227;
   - Ahn–Bellacosa–Ravanini hep-th/0312176 and Ahn–Bajnok–Palla–Ravanini 0711.4101 (NLIE with two Dirichlet boundaries and topological sectors; their H_± = p(1 ∓ 8φ_±/β));
   - Caux–Saleur–Siano cond-mat/0306328;
   - Kormos–Takács on boundary form factors.

   A lattice or NLIE computation that ties the boundary field, the topological charge and the bound-state poles together would be an independent test. Does any of them contradict the corrected sign?
7. **Downstream text and code.**
   - Check that no place in the book, the lecture notes or `code/` still uses the old sign: grep for `eq-bsg-dirichlet`, `cnj-bsg-uv-ir`, `mu_`, `4\pi}{\beta`, `psi_0`, `phi_D`.
   - Check the sign of E_D in ch. 38 (bullet after `cnj-defect-parameter`). Note that the continuation β_SG → ib has its own sign ambiguity.
   - Check the proof of `prp-defect-fusion-k` (4) in ch. 43. The classical fused boundary has C_{0,1} = cosh(θ_d^cl ± iβ_SGφ_D/2). Does the classical limit of the corrected relation give η = −(λ+1)β_SGφ_D/2 and ϑ = (λ+1)θ_d^cl, consistently with (3)?
   - Check item (h) of `code/lectures/defects_l9_fusion.py`.
8. **Labels and citations.** No status label changed, and Appendix E is unchanged. Is that right? Spot-check that each citation in the changed paragraphs says what the source says.

## Report format (`briefs/referee_bsg_dirichlet_sign_report.md`)

A header with the date, the model and the commit reviewed. Then one section per item, with:
- the verdict (correct / wrong / unclear);
- your own computation (the script inline or its path, and its output);
- for errors, corrected text ready to paste, with file and line.

List the claims you could not check. End with other observations, labelled as unverified.
