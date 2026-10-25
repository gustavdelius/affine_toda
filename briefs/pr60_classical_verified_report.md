# PR #60: independently verified classical defects and Claude Lectures 1–6

- Date: 2026-10-25.
- Model: Codex / GPT-6.
- Reviewed source: `eecc48fdf06bc079af2175fac53ea2d5d2d92abf`, the merge of PR #60; PR head `9f44dedfcdf88b4a84950c1188899e9082f749bc`.
- Review copy: `referee_scratch/pr60_snapshot/`. Paths and line numbers below are relative to that copy. Match stable ids and quoted text before applying an edit to another commit.
- Reports checked: `briefs/referee_defects_39_l5_l6_report.md`; its byte-identical duplicate `briefs/referee_defects_38_l6_report.md`; `briefs/referee_lectures_defects_report.md`, Items 1 and 2. The already integrated Gemini lecture report is excluded.
- Independent evidence: `briefs/pr60_classical_verified_scripts/independent_checks.py`, `lecture_checks.py`, and their `.log` files, preserved from `referee_scratch/pr60_classical/`. They import no book or previous referee implementation. `lecture_checks.py` imports only the local `check` and `zmat` utilities from the first new script.

## Editor-authorized integration exception

**Claude may integrate the verified and modified actions in this report directly, without independently re-deriving them again.** The editor explicitly requested this exception to Rule 1 of `AGENTS.md` and the independent re-check and label-rederivation requirements in step 3 and step 4 of `.agents/skills/integrate-referee-report/SKILL.md` (also exposed through `.claude/skills/`). This exception is specific to the actions and evidence below; it does not change the repository's general rules. No additional mathematical checks are required before implementing those actions. Normal regression tests, reference checks, checking the diff against the intended replacement, and scratch rendering still apply.

An **unclear** priority or blanket coverage claim below is disposed of explicitly: do not add it, do not use it to upgrade a label, and retain the narrower verified statement. This report does not transfer unfinished mathematical research to Claude. It distinguishes accepting a source statement from accepting the old report's purported derivation. Many source statements are correct while their old referee explanations are not.

## Integration decisions at a glance

| Original item | Independent verdict | Concrete disposition |
|---|---|---|
| Classical 1, `prp-defect-involution` | Correct source; old exceptional commutant generators and bulk-sign notation are wrong | Keep the mixed theorem/numerical label and the explicit finite-dimensional qualification |
| Classical 2, `prp-sg-spin-three` | Correct source; old report misstates one light-cone sign | Keep theorem and formulas |
| Classical 3, `prp-dc-gauge` | Correct, including the new type-I+ proof | Add corrected complete proof packet A; upgrade this proposition to `[Theorem]` |
| Classical 4, novelty and Penati–Polvara comparison | Priority claims unclear; alleged elimination of mixed velocities wrong | Do not add old proposed equivalence claim; retain bounded attribution already in the Notes |
| Classical 5, citation table | Wrong as a table of exact source quotations/locations; broad citations in book are appropriate | Retain book citations; use source ledger below, not the old table |
| Classical 6, halves and Bäcklund reduction | Correct source; old report reverses `P_d` | Keep source `P_d=g-f`, unhalved coordinates |
| Classical 7, cycle classification and scan | Correct | Keep theorem; independent scan reproduces 23/78 and 64/5/2/7 |
| Classical 8, action involution and spin consequence | Correct source; old report's rotation orientation and description of potentials are wrong | Keep source and its fixed-parameter action identity |
| Lecture 1.1, no-go | Correct source theorem; old unitarity-only completion is insufficient | Append packet B, using the right-incidence no-go argument |
| Lecture 1.1b, DMS bosonic normalization | Correct mathematical correction; claimed existing note is absent | Add packet C source note plus matching Appendix D entry |
| Lecture 1.2, Bianchi and named geometrical/soliton formulas | Correct for the actual formulas checked below; old ids are partly nonexistent | Retain `[Established]`; no upgrade or replacement |
| Lecture 1.3a, momentum theorem | Correct and general kinetic reduction already present | No change |
| Lecture 1.3b, alleged missing linear-potential branch | Wrong finding against reviewed text | Reject old replacement; actual theorem already includes linear field equations. Packet D clarifies two subsequent overcompressed sentences |
| Lecture 1.4, soliton delay | Correct source; old supporting signs and half-soliton energy formulas are wrong | Keep source `z=-coth((theta-theta_d)/2)` |
| Lecture 1.5, cycle theorem and five link changes | Correct source theorem; all five old replacement blocks are unusable | Packet E supplies four actual link edits; no matrix changes |
| Lecture 2.1, gauge theorem | Correct source; old report uses another, incorrect matrix formula | Keep theorem |
| Lecture 2.2, frozen classification | Correct under stated freedom of instantaneous values/velocities | Keep theorem and assumptions |
| Lecture 2.3, monodromy | Correct source; old report uses wrong product order/signs and wrong spatial domain | Keep circle statement and its proof |
| Lecture 2.4, bulk and defect Poisson brackets | Correct source; old arXiv number/credits are wrong | Keep labels and finite-dimensional caveat |
| Lecture 2.5, spin-three terms and general-rank matrix | Correct actual formulas | Keep; exact flux identity plus all-rank gauge proof below suffice |

## 1. Poisson algebra, uniqueness, and involution

Targets: `part8-defects/39-defect-conditions.qmd:265–302`, especially `prp-defect-involution`; `lecture-notes/defects/06-classical-algebra.qmd:394–482`, especially `prp-l6-bulk-r` and `thm-l6-sklyanin`.

**Verdict: correct, with the finite-dimensional qualification already in the text.** The theorem does not establish the constrained equal-time field theory's Poisson algebra. Keep that qualification.

For rank one the defect symplectic matrix is

\[
K=\begin{pmatrix}0&1\\-1&0\end{pmatrix},\qquad
K^{-1}=\begin{pmatrix}0&-1\\1&0\end{pmatrix}.
\]

Thus `{phi,psi}=-1` for type I−. Put `u=exp[b(phi-psi)/(2 sqrt(2))]`, `v=exp[-b(phi+psi)/(2 sqrt(2))]`. Absorbing the field-independent `sigma` into the spectral parameter gives

\[
L(\lambda)=\begin{pmatrix}u&v/\lambda\\(v\lambda)^{-1}&u^{-1}\end{pmatrix},
\quad \{u,v\}=\tfrac12\kappa uv,\quad \kappa=-b^2\epsilon/2.
\]

Differentiating all entries directly gives exactly

\[
\{L_1,L_2\}=\left[
\frac{\kappa}{2(z^2-1)}\{(z^2+1)H\otimes H+2z\Pi\},L_1L_2\right].
\]

The independent script checks all 16 entries as rational identities in `u,v,lambda_1,lambda_2`. Equivalently the coefficients of `J tensor O` and `O tensor J` give `2A-2G/z=kappa` and `2A-2Gz=-kappa`. Their supports are disjoint. The book's coefficient-extraction proof of uniqueness is complete: the coefficient matrices first force a commuting matrix to be scalar on the two equal-weight states, then the two independent rows `(1/lambda_2,1/lambda_1)` and `(1/lambda_1,1/lambda_2)` force the same scalar on the remaining two states. This requires `lambda_1^2 != lambda_2^2` and nonzero spectral parameters.

Independent exact linear systems have commutant dimensions 1, 2, 2 at `(lambda_1,lambda_2)=(2,3),(2,2),(2,-2)`. **The old report's suggested extra generators `diag(0,1,1,0)` and `diag(0,1,-1,0)` do not commute.** The new script verifies both failures. They are not in the source and must not be inserted. Also the r-matrix itself is singular at those exceptional spectral ratios; a commutant dimension is not a statement of a regular Sklyanin relation there.

The type-I+ reduction by right multiplication by `sigma_x`, field reflection `psi -> -psi`, and the field-independent spectral factor preserves the bracket because `sigma_x tensor sigma_x` commutes with r. It changes the field bracket to `-1`, as required. In lecture variables, `Phi=beta_SG phi`, it gives

\[
r(z)=-\frac{\beta_{\rm SG}^2}{16(z^2-1)}
[(z^2+1)h\otimes h+4z\Pi].
\]

`lecture_checks.py` derives the bulk bracket directly from `{Phi,Phi_t}=1` (the missing `beta_SG^2` restores by linearity) and `a_x=i Phi_t h/4 + m_0(lambda s_Phi+lambda^{-1}s_{-Phi})/4`. All 16 entries give **plus** `[r,a_x tensor 1+1 tensor a_x]`. The transport generator is `-a_x`; the endpoint identity integrates to the **minus** quadratic bracket for T. The old report's expression `a_x=-A_x` is not the source convention.

For higher rank a fresh orthonormal basis of the traceless hyperplane, the inverse of
`[[E,D],[-D^T,E]]`, analytic exponential derivatives, and the r-matrix in `prp-sklyanin` give:

```text
PASS Higher-rank Poisson bracket h2..5 both types/couplings max relative 8.031e-15
```

This independently supports `[Numerical]` for ranks 2–4; it is not a general-rank proof. No need to copy the old residual `1.63e-15` or reproduce its random samples. The inversion identity follows by Leibniz differentiation of `L L^{-1}=1`, giving `{L_1^{-1},L_2^{-1}}=-[r,L_1^{-1}L_2^{-1}]`. The product of independent factors with this same bracket obeys it again, and the trace of the commutator vanishes. This proves the algebraic involution step under the source's explicit independence assumption.

## 2. Spin ±3 and the jump's stored charges

Targets: book `prp-sg-spin-three`, `eq-dc-spin-three`, lines 229–253; lecture `eq-l6-D3`, lines 365–383.

**Verdict: correct. Keep the book's theorem and the lecture's more modest presentation of the formal expansion.** The normalization and spin signs are consistent.

Here is the full independent verification, using dimensionless `p=beta phi`, `q=beta psi`, `m=m_0`. Set

\[
p_x=q_t-m[\sigma\sin((p+q)/2)+\sigma^{-1}\sin((p-q)/2)],
\quad q_x=p_t+m[\sigma\sin((p+q)/2)-\sigma^{-1}\sin((p-q)/2)].
\]

For `s=±1`, let `d_p=p_t+s p_x` and
`e_p=2p_tt+2s d_t p_x+m^2 sin p`, and similarly for q. The latter follows from the bulk PDE and is precisely `partial_s^2 p`. Define

\[
j_s(p)=\frac{s}{8m^3}[m^2\cos p\,d_p^2-e_p^2+d_p^4/4],
\quad D_s=\frac{s\sigma^s}{2m^2}\cos\frac{p+s q}{2}d_p d_q
+\sigma^{3s}\left(\cos\frac{p+s q}{2}-\frac23\cos^3\frac{p+s q}{2}\right).
\]

`independent_checks.py` differentiates with respect to the independent jets `p,q,p_t,q_t,p_tt,q_tt`, substitutes these expressions, and simplifies `j_s(q)-j_s(p)+d_t D_s`. Both residuals are exactly zero with `m,sigma` symbolic. This is stronger evidence than a passage plot. The bulk relation is independently checked as `partial_- T4=partial_+ X2`, with `T4=(partial_+^2 p)^2-(partial_+p)^4/4` and `X2=-m^2 cos(p)(partial_+p)^2`.

```text
PASS Spin -3 full defect jet identity
PASS Spin -3 bulk continuity
PASS Spin 3 full defect jet identity
PASS Spin 3 bulk continuity
```

The cancellation of squared second time derivatives also follows without expansion: `partial_+^2(p-q)=2 partial_t partial_+(p-q)+m^2(sin p-sin q)` contains only first time derivatives; the difference of the two squares is therefore linear in the remaining second derivatives. The lower chirality uses the sum instead of the difference.

For the soliton put `y=m(x cosh(theta)-t sinh(theta))`. Then
`partial_s p=2m exp(-s theta) sech(y)` and
`partial_s^2 p=-2m^2 exp(-2s theta) sech(y)tanh(y)`.
With `cos p=1-2 sech(y)^2`, the integrals `int sech^2=2` and `int sech^4=4/3` give
`Q_{-3s}=-(2/3) exp(-3s theta)`. Independent quadratures at `theta=-0.7,0,1.1`, both signs, reproduce this to better than `1e-10`. Boost weights therefore assign spin −3 to the plus derivative, as the source says. At the odd jump `cos((p±q)/2)=-1`, so `D_s=-sigma^(3s)/3`; at the even vacuum it is `+sigma^(3s)/3`. Their difference is the stated soliton charge.

**Corrections to the old report, not the source:** its Item 2.2 prints the lower light-cone condition with a minus sign. It must be `partial_-(phi+psi)=+2m_0/(sigma beta) sin[beta(phi-psi)/2]`. Its claim about “both” endpoint charges is correct only with the source's chosen additive constants. The new script uses the actual source conditions.

**Caudrelier comparison checked independently.** [0704.2326](https://arxiv.org/abs/0704.2326), equation (3.28), has defect term
`epsilon alpha c (A^2+AB+B^2+2alpha^2)/3`, with `A=tilde v_x`, `B=v_x`, `c=cos((tilde v+v)/2)`. Equation (3.29) is `(tilde v-v)_x=epsilon 2alpha sin((tilde v+v)/2)`. Consequently
`A^2+AB+B^2=3AB+4alpha^2(1-c^2)`. With `epsilon alpha=-sigma m_0`, his defect term is exactly `-2m_0^3 D_{-3}`. This checks the limited source statement. It does **not** identify the bulk integrals on their different integration lines.

## 3. Gauge transformation for both types

Targets: book `prp-dc-gauge`, lines 318–331; lecture `eq-l6-ar-Ld`, lines 506–529.

**Verdict: correct. The type-I+ extension is proved, not only numerically checked.** The original proposal has the right sparse inverse provided the signs of both kinetic matrices in the outer exponentials are actually changed for type I+. The source's convention for the symbol `L_d` in the action identity is distinct from its use for the matrix; retain the explicit context.

Write `a_k=exp[b(phi_k-psi_k)/2]`, `b_j=exp[b(psi_j-phi_{j+1})/2]`, so `a_k^2=g_k`, `b_j^2=f_j`. This exponential notation avoids a principal-square-root convention at complex fields. For type I−,
`L_{kk}=a_k`, `L_{j+1,j}=(sigma/lambda)b_j`. In the plus gauge equation the orders `lambda` and `lambda^{-1}` cancel identically, and the diagonal order gives `b partial_+(phi_k-psi_k)=m sigma(f_{k-1}-f_k)`. In the minus equation the diagonal and `lambda^{-2}` orders cancel, and the subdiagonal order gives `b partial_-(psi_j-phi_{j+1})=m(g_j-g_{j+1})/sigma`. These are precisely the source's two halves.

For type I+ let `x=lambda/sigma`, `c(lambda)=1-x^h`. Since `S^h=1`, `(1-xS)sum_{j=0}^{h-1}x^j S^j=c(lambda)`. The scaled inverse

\[
M=c(\lambda)L^{-1},\quad M_{kk}=a_k,
\quad M_{j,j+1}=-x b_j
\]

has no other entries, and `det M=c(lambda)`. The gauge equation reverses to `partial_s M=M A_s(phi)-A_s(psi)M`. At order lambda in the plus equation its only nonidentity is

\[
-\frac b\sigma b_j\partial_+(\psi_j-\phi_{j+1})
=m[a_j e^{b(\phi_j-\phi_{j+1})/2}-e^{b(\psi_j-\psi_{j+1})/2}a_{j+1}]
=m b_j(g_j-g_{j+1}).
\]

At order 1 in the minus equation it is

\[
b a_k\partial_-(\phi_k-\psi_k)
=\frac m\sigma[b_k e^{b(\phi_k-\phi_{k+1})/2}-e^{b(\psi_{k-1}-\psi_k)/2}b_{k-1}]
=\frac m\sigma a_k(f_k-f_{k-1}).
\]

The other orders cancel: the order `lambda^2` plus equation uses `b_j e^{b(phi_{j+1}-phi_{j+2})/2}=e^{b(psi_j-psi_{j+1})/2}b_{j+1}`, while the order `lambda^{-1}` minus equation uses `a_{j+1}e^{b(phi_j-phi_{j+1})/2}=e^{b(psi_j-psi_{j+1})/2}a_j`; the derivative orders are identities. This argument works for every h, including h=2 where some matrix supports coincide but the spectral powers still separate.

The result holds away from the zeros of c and then everywhere on the nonzero spectral plane by the Laurent-polynomial identity. In particular `det L^{I+}=c^(h-1)`. The new script also forms the original full exponential matrix directly from the kinetic matrices, obtains its derivatives from the defect conditions, and compares both gauge equations for h=2,3,4,5, both couplings and both types:

```text
PASS All-rank both-type both-coupling gauge checks max 6.553e-14
```

**Packet A — ready to integrate.** In `part8-defects/39-defect-conditions.qmd:320`, replace the mixed label by `**[Theorem]**`. At line 324 replace its first sentence by “For each type the two gauge equations are equivalent to its two light-cone Bäcklund equations, displayed in the proof below.” After the type-I− paragraph at line 328, insert:

```markdown
For type I$_+$ put $x=\lambda/\sigma$ and $c(\lambda)=1-x^h$. Since $S^h=1$, the scaled inverse $M=c(\lambda)L_d^{-1}$ has the diagonal entries $g_k^{1/2}=e^{b(\phi_k-\psi_k)/2}$ and the entries $-x f_j^{1/2}=-x e^{b(\psi_j-\phi_{j+1})/2}$ in positions $(j,j+1)$, with no others. Its determinant is $c(\lambda)$. The gauge equation is equivalent, away from the zeros of $c$, to $\partial_\pm M=M\mathcal A_\pm(\vec\phi)-\mathcal A_\pm(\vec\psi)M$. In the $\partial_+$ equation the orders $\lambda^2$ and $1$ cancel identically. The order $\lambda$ gives
$$
\partial_+(\psi_j-\phi_{j+1})=-\frac{m\sigma}{b}(g_j-g_{j+1}).
$$
Indeed, after division by $f_j^{1/2}$, the two off-diagonal products are $g_j$ and $g_{j+1}$. In the $\partial_-$ equation the orders $\lambda$ and $\lambda^{-1}$ cancel identically, and the diagonal order gives
$$
\partial_-(\phi_k-\psi_k)=\frac{m}{b\sigma}(f_k-f_{k-1}).
$$
Here the two products, divided by $g_k^{1/2}$, are $f_k$ and $f_{k-1}$. These are the two halves of the type-I$_+$ conditions. Since the original gauge equation is a Laurent-polynomial identity, its validity away from the zeros of $c$ extends to all nonzero $\lambda$. Also $\det L_d=c(\lambda)^{h-1}$, so the matrix is invertible except at finitely many spectral values.
```

The current paragraph already supplies the complete type-I− calculation and the curvature argument. Retain them. Add an exact symbolic test of the two displayed type-I+ coefficient identities to `code/defect_ar_gauge.py` and update its `code/README.md` row. For a test that proves the local identity independent of rank, use six symbols `phi_j,phi_{j+1},phi_{j+2},psi_j,psi_{j+1},psi_{j+2}` and the exponentials above; separately check the wrap-around case h=2. The derivation here supplies the required proof; new code is regression coverage.

Do **not** attempt the old report's edit to an existing Appendix E gauge row: no such row exists in this snapshot. There are also no existing entries for `prp-sg-spin-three` or `prp-defect-involution`. **Add no new Appendix E entries for these three results in this round.** Sections 1–3 independently verify the mathematics, but this pass does not establish originality. Keep the existing Appendix E entries unchanged for these items, including the separate existing `prp-dc-types` entry; its statement and label remain unchanged. This is a final no-addition disposition, not a request for Claude to perform another priority search. The mathematical proof and theorem upgrade for `prp-dc-gauge` remain cleared. No Appendix D change follows from packet A.

## 4. Momentum halves, cycle classification, and the two types

Targets: book `eq-dc-halves`, `eq-dc-lightcone`, `eq-dc-components`, `prp-ar-only`, `prp-dc-types`; lecture `prp-l5-cycle` and its preceding derivation.

**Verdict: correct source.** There are several consequential errors in the old report's explanation, but no corresponding formula error in these passages.

### 4.1 General reduction

Put `D=1-E`, `E^T=-E`; thus `D^T=1+E` and `DD^T=D^TD=1-E^2`. Then
`ED-DD^T=-D`, `E+D=1`, `D^TD+ED^T=D^T`, `E-D^T=-1`, by ordinary multiplication. For the **unhalved** source coordinates `xi=D^T phi+D psi`, `eta=phi-psi`, the inverse is `(phi,psi)=((xi+D eta)/2,(xi-D^T eta)/2)`. The gradient equations reduce to
`grad_xi Pd=-grad_xi B`, `grad_eta Pd=grad_eta B`, so the mixed Hessian vanishes and

\[
B=f(\xi)+g(\eta),\qquad P_d=g(\eta)-f(\xi).
\]

The old report's `P_d=f-g` is wrong. Inserting `grad_phi B=D grad f+grad g`, `grad_psi B=D^T grad f-grad g` cancels the square terms and gives the difference `4 grad f dot grad g`. Adding/subtracting the defect conditions yields `partial_+eta=-2 grad f`, `partial_-xi=2 grad g`. Differentiation gives
`partial_+partial_- phi=-2(H_g grad f+D H_f grad g)`, while differentiating `2 grad f dot grad g=V(phi)-V(psi)` gives the negative of the same expression as the field gradient. This is the complete Bäcklund implication, with the analogous psi equation. It is a convention-sensitive derivation, not a separate new theorem needing a new label.

### 4.2 Classification and scan

Every product of the two finite exponential sums has exponent `(D x+y, D^T x-y)/2`. Their sum recovers x, so different pairs cannot cancel. A wanted left exponential requires `(x,y)=(alpha_i,D^T alpha_i)`, and a wanted right one requires `(alpha_i,-D alpha_i)`. The positive affine relation forces the two lists of y to coincide; otherwise y is orthogonal to every root except alpha_i and hence to alpha_i too, a contradiction. Thus `D^T alpha_i=-D alpha_{p(i)}` with a permutation having no fixed nodes. `P=-D^{-1}D^T` is orthogonal and `D=2(1-P)^{-1}`. No other nonzero exponent can occur because the required y vectors span the root space.

Forbidden products require `alpha_i dot D alpha_j=0` unless `j=i` or `j=p(i)`. Taking the scalar product of `D(alpha_i-alpha_{p(i)})=2alpha_i` with a neighbour alpha_k forces `k=p(i)` or `k=p^{-1}(i)`. Connectedness makes p a single orbit, and the automorphism condition forces all consecutive orbit nodes to be neighbours. Thus the entire diagram is a cycle; for two nodes they must have equal length. This leaves exactly `a_r^(1)`, with a one-step rotation in either direction. The coefficients of g are constant on the orbit, and those of f are fixed by the bulk coefficients, leaving sigma. This proves necessity and sufficiency under the source's listed assumptions, and justifies both theorem labels.

The fresh scan constructs the roots from Euclidean coordinate vectors, not imported Cartan data: standard A/B/C/D roots, explicit F4 and G2 roots, their coroots for the relevant twisted systems, and `(-2e_1,e_1-e_2,...,e_n)` for the nonreduced twisted chain. It enumerates every Gram-preserving permutation. Output:

```text
PASS Affine diagram scan 23 diagrams; 78 symmetries; [64, 5, 2, 7]; failures=[('a41', (2, 3, 4, 0, 1), 10), ('a41', (3, 4, 0, 1, 2), 10)]
```

The numerical four-way split is therefore independently reproduced, including the two ten-cross-term failures. The old lecture report's claimed antisymmetric-matrix derivation is not the derivation in the actual lecture. In particular, for a real antisymmetric A, `v^T(1-A^2)v=|v|^2+|Av|^2>0` for v nonzero: its claim `det(1-A^2)=0` is impossible. Those sentences are not in the reviewed source and must never be introduced by a link edit.

### 4.3 Action involution and spin consequence

Use the source orientation `pi alpha_j=alpha_{j-1}`, not the old report's opposite one. The exact identities are `rho D rho=D^T`, `rho E rho=-E`, `D^T pi^{-1}=-D`. In components `(rho v)_k=-v_{1-k}`. Under `(phi,psi)->(rho phi,rho pi^{-1}psi)`, direct substitution gives **`f_j -> g_{-j}` and `g_j -> f_{-j}`**. It exchanges the two halves; it does not leave each half fixed. Consequently the transformed action has kinetic matrices `(-D,-E)` and the coefficients sigma and sigma inverse exchanged, exactly the type-I+ action at the same named sigma. The independent script checks these identities at both couplings for h=2 through 5.

For spin two, `A(v,v,v)=2 sum v_k^3` is invariant under pi and changes sign under rho. Conjugation of the Lax pair by the cyclic shift leaves its spectral invariants unchanged; time reversal swaps the derivative chiralities. Thus if a charge of one chirality survives at every sigma, the opposite one survives after sigma inversion and the action involution. This verifies the source's logical consequence, not an independent construction of every local charge at every rank. The identity for the soliton delay follows algebraically by putting `z_a(u)=(exp(-u)-i exp(i gamma))/(exp(-u)-i exp(-i gamma))`: multiplication of numerator and denominator by `exp(u)` gives `z_a(-u)=exp(2i gamma)z_{h-a}(u)`. The rotation of the right tau index supplies the same phase.

## 5. Literature and priority ledger

The old report's statement that every source equation and line was checked is **wrong**. Its source line numbers and several equation descriptions fail direct inspection. This does not invalidate the book's broader references. The following is the independently read ledger. PDF text was read from `/home/gustav/referee-copy/papers/`, except Habibullin–Kundu, downloaded into the new scratch directory. Quotations are deliberately short.

| Source | Exact source evidence | What may be concluded |
|---|---|---|
| Bowcock–Corrigan–Zambon [hep-th/0305022](https://arxiv.org/abs/hep-th/0305022), (2.4)–(2.5), text lines 145–167 | “the energy-like conserved quantities will be given by the following expressions” | They impose conservation of the sum of the two chiralities. The source's description is correct |
| Bowcock–Corrigan–Zambon [hep-th/0401020](https://arxiv.org/abs/hep-th/0401020), (3.7), text line 291 | `partial_t K = K a_t^> - a_t^< K` | This is the temporal gauge equation, as the book says |
| Same paper, §5.1, text line 953 | “Defect conditions from Bäcklund Transformations” | The old report gives a different section title. The section does cross-differentiate conditions imposed for every x, so the book's use is appropriate |
| Caudrelier [0704.2326](https://arxiv.org/abs/0704.2326), (1.12)–(1.13), text lines 192–193 | `L_x=tilde U L-LU`, `L_t=tilde V L-LV` | These are equations for the defect matrix, not the wavefunction jump relations as the old table calls them |
| Same paper, (3.26), (3.28), (3.29), text lines 754, 775, 780 | The matrix and charge equations reproduced and reduced in §2 above | The sine-Gordon matrix and spin-three contribution are prior art; the old lecture report's arXiv `0803.2498` is not this paper |
| Habibullin–Kundu [0709.4611v2](https://arxiv.org/abs/0709.4611), §5.2, p.14 | “the dependence of the canonical momentum P− drops out completely” | Their discretized defect has canonical momentum variables and a delicate continuum limit. This supports the book's distinction |
| Same paper, (5.73), (5.76), (5.81) | `(5.73): [s_q^(+),s_q^(-)]=0`; `(5.76): {L_j(xi),tensor L_j(eta)}=delta_jk[r,L_j(xi) tensor L_k(eta)]` | (5.73) is the quantum algebra; (5.76) is the quadratic Poisson relation. The old table incorrectly identifies (5.73) as that bracket and invents an a2 matrix example; this paper treats sine-Gordon |
| Penati–Polvara [1902.10690v2](https://arxiv.org/abs/1902.10690), preceding (3.8) | “The corresponding energy-like and momentum-like currents of the form (2.14) are both conserved with” | They print field-only sinh-Gordon terms in (3.8)–(3.9). This is related work, not evidence for the old report's explanation that mixed velocities can be eliminated |
| Same paper, (3.16)–(3.17) | `phi^(-)->O_1 phi^(-), phi^(+)->O_2 phi^(+)`, `D->O_1^T D O_2` | Orthogonal kinetic redefinitions are prior art. The book already acknowledges them |
| Rogers–Schief, publisher's [contents](https://assets.cambridge.org/97805218/13310/toc/9780521813310_toc.pdf) | “1.3 Bianchi’s Permutability Theorem. Generation of Multi-Soliton Solutions” | The section reference is appropriate; the actual identity is independently proved by the new rational script |

The Penati–Polvara “elimination” explanation cannot be accepted. The instantaneous values `phi,psi,phi_t,psi_t` are independent at a type-I defect; its conditions determine the spatial derivatives. In dimensionless variables the mixed time-velocity derivative of `D_{-3}` is `sigma cos((p+q)/2)/m^2`, generically nonzero. A function of p and q alone has zero such derivative, as does any field-only endpoint improvement. Applying both light-cone defect conditions does not impose a relation on the independent time velocities. A comparison involving different derivative-dependent current improvements would require an explicit dictionary, which the old report does not provide. **Do not add the old proposed phrase “after eliminating the mixed velocity term”.** No new allegation against Penati–Polvara is required for integration; leave this optional citation addition out.

**Priority verdict: unclear, no new priority claims cleared.** Finding an elementary proof or a formula absent from a small list of papers does not establish novelty across the literature. Preserve the source's bounded “stated here” and “we did not find ... in the sources we checked”, and its explicit acknowledgements. Do not turn the old report's “is new to the book” into an assertion of first discovery. No label upgrade in this report depends on novelty.

## 6. Claude Lectures 1–5

### 6.1 No-go theorem and DMS normalization

**No-go: correct, but the old unitarity-only completion is insufficient.** The actual statement at `01-mirrors-and-windows.qmd:241` correctly distinguishes the two mixed histories. If `R_b T_a` is nonzero as a meromorphic function, its product with `S_ab(theta_1+theta_2)-S_ab(theta_1-theta_2)` forces equality on an open set. The invertible change of variables `u=theta_1-theta_2`, `v=theta_1+theta_2` then makes S constant. No parity assumption is needed. Under the last hypothesis, if any reflection is nonzero all transmissions vanish; if any transmission is nonzero all reflections vanish. To extend the reflecting branch to right incidence, apply the same mixed-history argument to particles arriving from the right. Unitarity then rules out a different alternative on the two sides.

The distinction matters: the constant triangular matrix `[[1,0],[c,-1]]` squares to the identity for every c. It has zero transmission in one direction and nonzero transmission in the other when c is nonzero, while satisfying algebraic matrix unitarity. Thus the old report's claim that unitarity alone forces the right transmission to vanish is false. If left transmission vanishes, the right unitarity relation gives `bar R(theta)bar R(-theta)=1`, so right reflection is nonzero; the right-incidence no-go then forces right transmission to vanish. If left reflection vanishes, left unitarity makes both transmission amplitudes nonzero, and its mixed relation directly forces right reflection to vanish. This supplies the complete argument without assuming positive-definite Hilbert-space unitarity.

**Packet B.** Append inside the proof of `thm-l1-nogo`, after line 248:

```markdown
Under the last hypothesis, a nonzero $R_b$ forces every $\mathbb T_a$ to vanish, while a nonzero $\mathbb T_a$ forces every $R_b$ to vanish. The unitarity relations (-@eq-l1-unitarity) and their right-incident counterparts exclude $R_a=\mathbb T_a=0$. If all $\mathbb T_a$ vanish, right-incident unitarity gives $\bar R_a(\theta)\bar R_a(-\theta)=1$; applying the same mixed-history argument to particles arriving from the right then forces every $\bar{\mathbb T}_a$ to vanish. If all $R_a$ vanish, left-incident unitarity makes $\mathbb T_a$ nonzero, and the mixed relation forces every $\bar R_a$ to vanish. Thus the same alternative holds for incidence from either side.
```

The old proposed block refers to a differently worded theorem and discards the original proof; append this paragraph instead.

**DMS: mathematical correction verified, but absent from source Notes.** The source [hep-th/9409076](https://arxiv.org/abs/hep-th/9409076), equations (4.1)–(4.2), prints `L=1/2[(partial_mu phi)^2-m^2 phi^2-g delta(x)phi^2]` and `[Box+m^2+g delta(x)]phi=0`. Integrating the latter gives `[phi_x]=g phi`. With `T=1+R`, the plane wave equation is `2ik R=gT`, so

\[
T=\frac{2ik}{2ik-g},\qquad R=\frac{g}{2ik-g}.
\]

Thus the denominator in rapidity is `sinh(theta)+ig/(2m)`. The paper's (4.8) instead prints `T=sinh(beta)/(sinh(beta)+ig/(4m))`. Its (4.6) also has `g/(4m)` multiplying the sum of the two boundary values. Against (4.2), that coefficient must be `g/2`, not simply `g/(2m)`. The old report is too imprecise about (4.6). There is no bosonic chi definition needed here; do not confuse this with their fermionic angle chi.

The book's lecture code was run only as a **second check**, after the derivation: `code/lectures/defects_l1_mirrors_windows.py` reports `PASS (a) delta defect: T = 2ip/(2ip-g), R = g/(2ip-g)` and all of its checks pass. This comparison is stored in `book_l1_comparison.log`; it is not the independent evidence.

**Packet C, source clarification, fully cleared.** At `lecture-notes/defects/01-mirrors-and-windows.qmd`, Notes and references, append:

```markdown
For the bosonic defect, [@delfino_mussardo_simonetti1994, eq. (4.2)] gives the jump $[\phi_x]=g\phi$. Hence the $ig/(4m)$ in their amplitudes (4.8) must be $ig/(2m)$, as in (-@eq-l1-delta-rt). In their unscaled matching equation (4.6), the coefficient of $\phi_++\phi_-$ must be $g/2$. The action (4.1) and field equation (4.2) fix this normalization.
```

In the same change, append the matching bullet to the existing Delfino–Mussardo–Simonetti entry in `appendices/D-corrections.qmd:30–33`:

```markdown
- **(4.6) and (4.8), the bosonic defect.** Their field equation (4.2) gives $[\phi_x]=g\phi$. In (4.6) the coefficient of $\phi_++\phi_-$ is therefore $g/2$; in (4.8) the rapidity denominators contain $ig/(2m)$, not $ig/(4m)$. *Normalization error*; derived in the companion defect lectures from the action and checked in `code/lectures/defects_l1_mirrors_windows.py`.
```

No new mathematical formula is needed in the lecture body; it already has the correct normalization. Do not add a second claim about the unrelated fermionic chi branch. Appendix E contains no matching entry for this bosonic normalization correction. Leave it unchanged and add no originality entry for this correction: this pass establishes the source discrepancy, not priority for identifying it. The lecture Notes and Appendix D must be updated together as specified above. No additional priority research is required for this disposition.

### 6.2 Bianchi, pseudosphere, and breathers

**Verdict: correct for the actual stated formulas.** The old report names nonexistent ids `eq-l2-pseudosphere`, `eq-l2-tractrix`, `eq-l2-2soliton`; the reviewed two-soliton id is `eq-l2-two-soliton`. Its numerical residual cannot stand as proof of a blanket “all geometry” verdict.

The new rational proof uses `u=tan(Phi_0/4)`, `v=tan(Phi_1/4)`, `w=tan(Phi_2/4)`, `C=(a+b)/(a-b)` and
`z=(u+C(v-w)/(1+vw))/(1-uC(v-w)/(1+vw))`. The derivative of `4 arctan(z)` is substituted into all four outgoing Bäcklund equations. Incoming derivatives are supplied by the two original Bäcklund pairs; the starting derivative remains arbitrary. The sine addition formula is rationalized. All four residuals are **exactly zero**, for arbitrary background and arbitrary unequal parameters. This checks sufficiency, not just vacuum examples.

The pseudosphere `r(s,tau)=(sech(s)cos(tau),sech(s)sin(tau),s-tanh(s))` has first form `tanh(s)^2 ds^2+sech(s)^2 d tau^2` and curvature −1 away from `s=0`; at the rim `r_s=0`, so a smooth-surface claim must exclude it. With `s=u+v`, `tau=u-v`, the mixed coefficient is `cos Phi=1-2sech(s)^2`, giving the stated asymptotic-coordinate form. The script verifies the metric and curvature exactly. The soliton is the vacuum Bäcklund solution; applying the proved Bianchi relation gives the two-soliton formula. Its rational tangent expression has interaction coefficient `((a-b)/(a+b))^2`. Substitution of conjugate rapidities gives `tan(Phi/4)=cot(alpha) sin(m_0 t sin(alpha))/cosh(m_0 x cos(alpha))`; the script verifies the sine-Gordon PDE exactly in the tangent variable.

Keep `[Established]` for `prp-l2-bianchi`; no requested upgrade follows. Historical dates, Hilbert's global immersion theorem, and the assertion that *every* geometrical sentence was verified are not cleared by a finite PDE test. They are outside the concrete formula findings and require no integration action here. Do not repeat the old blanket certificate.

### 6.3 Momentum and scalar classification

**Momentum theorem correct; old “proof completion” is already present.** The actual source begins with the general `A(Phi,Psi) Phi_t+C(Phi,Psi) Psi_t-B`, defines `w=C_Phi-A_Psi`, and obtains the velocity polynomial shown at lines 180–183. Its quadratic coefficients force `w^2=1`; smoothness on a connected field domain gives constant w. The closed one-form argument at line 192 already proves reduction to `±Phi Psi_t` modulo a total derivative. The wave equation in field space then gives `B=f(Phi+Psi)+g(Phi-Psi)`, `Pd=g-f`, and `2f'g'=V_L-V_R`. For sine-Gordon the separated factors fix a single nonzero sigma. This is the complete proof of `thm-l3-momentum:267`; keep it.

**The alleged omitted linear branch is not omitted.** The actual `prp-l3-which:310` says “both field equations are linear, with equal masses, and B is quadratic”. This includes constant-source equations and hence linear potentials. It does not state the narrower pure-quadratic list quoted in the old report.

For clarity, solve the zero-separation-constant branch independently: let `f'=a_1 Sigma+a_0`, `g'=b_1 Delta+b_0`. Then

\[
V_L(P)=2a_1b_1P^2+2(a_1b_0+a_0b_1)P+C_L,
\quad V_R(Q)=2a_1b_1Q^2+2(-a_1b_0+a_0b_1)Q+C_R,
\quad C_L-C_R=2a_0b_0.
\]

The new script checks this polynomial identity. In particular `a_1=1,b_0=1,a_0=b_1=0` gives opposite linear slopes, contrary to the old replacement's “equal slope” assertion. `a_0=b_0=1,a_1=b_1=0` gives constant potential difference 2, contrary to its assertion that one constant must vanish. The true separated identity is `f'''/f'=g'''/g'`, not its earlier minus sign. Its frequency in Phi is twice the frequency in Sigma. **Reject the entire old proposed replacement.** It would introduce errors into the currently correct classification.

Two surrounding sentences still compress “exponential or quadratic/linear” into “exponential”. **Packet D:** after `prp-l3-which`, replace “Theories whose interaction is not exponential have no such defect.” with:

```markdown
Outside this list of exponential potentials, linear field equations and free massless fields, there is no defect of the stated form.
```

In `05-why-the-cycle.qmd:18`, replace “and only exponential potentials pass the test” with “and the non-linear interacting cases have exponential potentials”. These wording changes are justified by the explicit polynomial branch above. They change no theorem label.

### 6.4 Delay and the old report's normalization mistakes

**Actual `prp-l4-delay:135` correct.** Substitution of `tan(Phi/4)=Y`, `tan(Psi/4)=zY` gives `e^{-theta}(1-z)=-sigma(1+z)` and `-e^theta(1+z)=sigma^{-1}(1-z)`. Thus `z=(e^{-theta}+sigma)/(e^{-theta}-sigma)=-coth((theta-theta_d)/2)`. Both rational equations are independently checked. At the coincident rapidity the limiting absorption/emission pairs follow by substitution in the original equations; no division by `1-z` is allowed there.

Analytic continuation gives `T(theta)=z(theta+i pi/2)=(1+i exp(theta-theta_d))/(1-i exp(theta-theta_d))`, a unit-modulus factor at real rapidity. These agree with the actual lecture. The old report drops the minus sign of coth and gives extra energy identities with inconsistent mass normalization. They are rejected. With `M=8m_0/beta_SG^2`, the jump's physical energy and momentum are `M cosh(theta_d)` and `M sinh(theta_d)` relative to the even vacuum. The unshifted endpoints of each light-cone defect charge are ± half the corresponding soliton charge. The odd-sector delay is `-tanh((theta-theta_d)/2)=1/z(theta)`, whose square is the usual interaction coefficient; the odd wave factor is its analytic continuation. None of this requires editing the correct delay theorem.

For complex Y the arctangent field is locally defined modulo `4pi` and has logarithmic branch points; its derivatives and the trigonometric equations are meromorphic functions of Y. The proof's phrase “the fields are meromorphic” at `04-what-gets-through.qmd:146` is too strong. Replace that clause by “their derivatives and the trigonometric expressions in the conditions are meromorphic in $Y$, so the conditions hold on every regular local branch”. This precise qualification does not change the theorem on regular solution regions.

### 6.5 Actual cross-reference edits

All five old replacement blocks quote passages absent from the reviewed source; some confuse type I± with type II. Do not apply them. **Packet E** consists of these actual edits:

1. `03-freezing.qmd:450`, last sentence: replace “The book relates the two types by time reversal (see its [defect Lagrangian](book:sec-defect-lagrangian)).” by “The book gives the [Dynkin-reflection relation](book:prp-dc-types) as well as the time-reversal relation.”
2. `05-why-the-cycle.qmd:342`, final two sentences: replace “The book relates the two types by time reversal instead (see its [defect Lagrangian](book:sec-defect-lagrangian)); both relations hold (@exr-l5-types).” by “The book gives the same [relation between type I$_+$ and type I$_-$](book:prp-dc-types), together with the time-reversal relation; both hold (@exr-l5-types).”
3. After the display `eq-l5-lightcone:153`, add “These are the book's [light-cone defect conditions](book:eq-dc-lightcone), obtained from its [two halves](book:eq-dc-halves).”
4. After `eq-l5-ar-lightcone:309`, add “The exponentials are the book's [component defect factors](book:eq-dc-components).”

All target ids were found in the reviewed book; §§3–4 above independently verify that these are the matching formulas. No link should describe `prp-dc-types` as a relation to type-II defects. The existing classification link to `book:prp-ar-only` is already correct and needs no edit.

## 7. Claude Lecture 6: complete proof audit

### 7.1 `thm-l6-gauge` and `thm-l6-frozen`

**Both correct.** The actual matrix is `L_d=exp(i Delta h/4)+(sigma/lambda)exp(-i Sigma h/4)s`, and the actual equation is `partial_s L_d=L_d A_s(Psi)-A_s(Phi)L_d`. The old report's differently ordered matrix and its `U(Psi)G-GU(Phi)` equation are not a transcription of this source.

The new script writes `u=exp(i Phi/4)`, `v=exp(i Psi/4)`, differentiates the entries algebraically, and substitutes `partial_+Psi=partial_+Phi+2sigma m_0 sin(Sigma/2)` and `partial_-Psi=-partial_-Phi+2m_0 sin(Delta/2)/sigma`. Both matrix residuals are exactly zero; `det L_d=1-sigma^2/lambda^2` independently follows. Conversely the diagonal coefficient at order 1 of the plus equation and the off-diagonal coefficient at order `lambda^{-1}` of the minus equation force precisely these two conditions, since the exponential factors never vanish. The proof is complete.

For the frozen classification, the extreme orders require `G_0=pD+q Obar`, `G_1=p'D^{-1}+q'O`. Comparing each diagonal and off-diagonal coefficient at the middle orders gives constant `p,q,p',q'` and

\[
p\partial_+\Delta=-2m_0q'\sin(\Sigma/2),\quad
q'\partial_-\Sigma=2m_0p\sin(\Delta/2),
\]
\[
q\partial_+\Sigma=-2m_0p'\sin(\Delta/2),\quad
p'\partial_-\Delta=2m_0q\sin(\Sigma/2).
\]

The two basis matrices of each support are independent (multiplication by the invertible diagonal exponential reduces them to 1 and h). If just one of p and q' vanishes, a sine of freely variable defect values would vanish identically, impossible; similarly for q and p'. If both pairs survive, the four first derivatives are fixed by the values, contradicting free instantaneous velocities. Exactly one pair therefore survives, giving the two source matrices up to a constant scalar. This checks the classification the old report did not actually derive. Retain its hypotheses; it is not a classification of all possible time-dependent or higher-degree Darboux matrices.

### 7.2 `thm-l6-monodromy`

**Correct source, wrong old explanation.** The actual monodromy is `T_R L_d^{-1} T_L` on the circle, not `T_L L_d^{-1} T_R` on an unrenormalized infinite line. From `partial_t T(x,y)=-a_t(x)T+T a_t(y)` and `partial_t L^{-1}=L^{-1}a_t(Phi)-a_t(Psi)L^{-1}`, Leibniz differentiation cancels the two defect endpoints exactly. The remainder is `-a_t(Psi;ell) T_d+T_d a_t(Phi;-ell)`. Smooth matching without winding equates those two endpoint connections, and its trace vanishes. This is a complete proof. No replacement needed.

### 7.3 `prp-l6-bulk-r`, `thm-l6-sklyanin`, and Item 2.5

These are verified in §§1–3 above. Keep `[Established]` for the bulk r-matrix, `[Theorem]` for the explicitly finite-dimensional defect statement, and the actual assumptions governing the product Poisson algebra. The old report's `0803.2498`, “Avan–Doikou–Sfetsos”, and equation id `eq-l6-r` do not match the source and must not be introduced. The actual Caudrelier paper is `0704.2326`, and the actual equation id is `eq-l6-rmatrix`.

The spin-three terms are not checked by equating a single bare fourth-spin density across the defect, as the old report claims. The correct flux includes the companion `X2` term and its chirality-dependent sign; §2 verifies it exactly. For the general-rank matrix, §3 proves both light-cone equations for every h, not just a time equation at h=3. This validates the actual assertion following `eq-l6-ar-Ld`.

## Integration

This section records the independently verified dispositions; **no book, lecture or existing report file was edited by this referee**.

| Action | Status | Evidence | Files for Claude |
|---|---|---|---|
| A: type-I+ gauge proof and label | Verified | Universal coefficient proof in §3; direct original-matrix checks | `part8-defects/39-defect-conditions.qmd`, `code/defect_ar_gauge.py`, `code/README.md`; no Appendix E addition |
| B: no-go proof completion | Modified | Meromorphic identity and both incidence directions, with algebraic unitarity | `lecture-notes/defects/01-mirrors-and-windows.qmd` |
| C: DMS note and Appendix D bullet | Modified | Primary equations (4.1),(4.2),(4.6),(4.8); independent matching derivation; secondary book-code comparison | Same lecture and `appendices/D-corrections.qmd` together; no matching Appendix E entry exists and none is added |
| D: clarify non-exponential exceptions | Modified | Exact affine-linear factorization, already included by actual theorem | `lecture-notes/defects/03-freezing.qmd`, `05-why-the-cycle.qmd` |
| E: actual book links | Modified | Target ids and independent equivalence proofs | Same two lectures |
| Complex-field regular-branch qualification | Modified | `arctan Y=(log(1+iY)-log(1-iY))/(2i)` has logarithmic branches | `lecture-notes/defects/04-what-gets-through.qmd` |
| All other checked theorem/formula retentions | Verified | §§1–4 and 6–7; 67 new-script PASS lines | No change |
| Old scalar replacement, old five link blocks, wrong source explanations | Rejected | Counterexamples and exact text matching above | Do not insert |
| Priority claims, “all geometry” blanket certificate, mixed-velocity elimination story | Rejected as integration instructions | Insufficient scope or explicit counterargument | Do not add; no fresh research needed for this integration |
| New Appendix E entries for gauge, spin-three, finite-dimensional involution or the DMS correction | Not authorized in this round | Mathematical verification does not establish priority | Retain existing entries unchanged for these items; add none |

The independent scripts finish with 52 PASS lines in `independent_checks.log` and 15 PASS lines in `lecture_checks.log`, no FAIL. Rerun from their preserved directory with `python3 -u independent_checks.py` and `python3 -u lecture_checks.py` if desired as a regression comparison, not as another independent mathematical pass. The book's Lecture 1 comparison run is separately recorded and is not counted among those 67 checks.

After integration run the changed production scripts, `python3 check-refs.py`, `python3 check-refs.py lecture-notes/defects`, and the ordinary scratch renders. Match the integration diff against this packet and keep this report's evidence unchanged. A completed integration record may cite this editor-authorized independent verification instead of repeating it. The type-I+ gauge upgrade has received this independent referee pass; do not describe it as an unrefereed upgrade. No other new theorem label is authorized by this report.
