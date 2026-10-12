# Opus referee prompt: half-line open problem 11 (tier-B follow-up)

Written 2026-10-12 by the Claude coordinator, after reading the Gemini verdict `briefs/halfline_open_11_referee.md`.

**What the Gemini pass established.** It independently established the following, so do not redo them:
- the antisoliton shift $\bar\eta=\eta+i\pi\omega(l-N/2)$, with its own null-space solver;
- the crossing function;
- the $n=1$ match with Ghoshal's Dirichlet amplitude.

Its argument for $B_1\to\bar B_1$ is acceptable: the soliton-preserving $K$ keeps flavours and reverses the constituent order, whereas the conjugating $K$ maps $V\to\bar V$ and keeps the order. Its own Opus items (the SP soliton bootstrap for the CDD factor, the axioms of the continued amplitudes, the RV 2018 match) are research follow-ups, not referee items.

Read:
- the brief `briefs/halfline_open_11_soliton_preserving_scalar.md`;
- the report `briefs/halfline_open_11_report.md`;
- the Gemini verdict;
- the scripts in `briefs/halfline_open_11_scripts/`. Copy `textbook_code/*.py` and `gz.py`, `nk_common.py`, `rhonp.py` from `briefs/halfline_open_03_scripts/` alongside.

## Items

1. **Claim 5: real-coupling Dirichlet is not classically integrable for $n\ge2$.**
   - Gemini's re-derivation is not the worker's argument and looks incomplete. It assumes $\sigma(E_{\pm\alpha})=KE_{\mp\alpha}K^{-1}$ with unit coefficients. But at a Dirichlet value the Lax step terms carry constants $e^{\beta\alpha_j\cdot\phi_0/2}$. With coefficients $c_j,c'_j$, the bracket gives $\sigma(H_{\alpha_j})=c_jc'_jH_{\alpha_{\pi(j)}}$, and $c_jc'_j=-1$ seems possible.
   - The worker's own argument tracks these constants ($\psi_j=\lambda^2$ at every node, $\lambda^{2N}=1$).
   - Re-derive the claim carefully. Include possible permutations of nodes by $\sigma$ (diagram automorphisms), and the stated restrictions (field-independent $K$, linear $\sigma$, fundamental representation). Say whether those restrictions leave a loophole.
   - Check prior art: Bowcock–Corrigan–Dorey–Rietdijk hep-th/9501098, Corrigan's reviews, and any paper on Dirichlet conditions in $a_n$ Toda. Is "real Dirichlet is not integrable for $n\ge2$" already known, or contradicted?
2. **Claim 4: `sec-soliton-preserving` should allow coweights outside $2\Lambda_W^\vee$ with $K_+=e^{i\pi\mu\cdot H}$.** Gemini restated this without a check.
   - Verify it by hand from the book's derivation in ch. 30, or with `textbook_code/boundary_classical.py`.
   - Check also that the resulting $\operatorname{Re}\phi(0)$ is a vacuum (the brief's trap: the Dirichlet value must be a vacuum).
   - Draft the corrected book sentence.
3. **Claim 3: the $\mathbb Z_N$ obstruction.** Does it survive allowing covariance up to a diagonal gauge, a twist, or a $\theta$-dependent scalar? That is, $ZKZ^{-1}=c(\theta)DKD^{-1}$ with $D$ diagonal and constant. If a gauge-covariant SP matrix exists, the claim needs rewording.

For each item give the verdict (CONFIRMED / CORRECTED / REFUTED / UNCHECKED), the evidence, and book text. Give script paths.
