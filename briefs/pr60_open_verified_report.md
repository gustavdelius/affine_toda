# PR #60: independently verified open items and lecture census

- Date: 2026-10-25.
- Referee: Codex / GPT-6, coordinating referee.
- Reviewed source: `eecc48fdf06bc079af2175fac53ea2d5d2d92abf`, the merge of PR #60. PR head: `9f44dedfcdf88b4a84950c1188899e9082f749bc`.
- Inputs: `briefs/referee_open_items_report.md`, all items 1–8, and `briefs/referee_lectures_defects_report.md`, §5.2. The latter's §§1–5.1 are covered by the companion classical and quantum reports.
- Snapshot: `referee_scratch/pr60_snapshot/`. Original book and lecture files were not edited. New independent computations are archived in `briefs/pr60_open_verified_scripts/`.
- The Gemini defect lecture report is excluded, at the editor's explicit instruction: it has already been integrated.

## Editor-authorized integration exception

The editor explicitly requested this independent verification round so that Claude can integrate its conclusions without checking the mathematics again. **For the dispositions certified in this report, Claude does not need to repeat AGENTS.md Rule 1, step 3 of `.agents/skills/integrate-referee-report/SKILL.md` (also exposed through `.claude/skills/`), or that skill's instruction to re-derive label upgrades.** This is an explicit task-specific exception, not a change to the standing rules.

Integrate the **new** dispositions below, not the original report's replacement blocks. Ordinary implementation checks, reference checks, script regression checks, and scratch rendering remain required. No commit or push is authorized by this report. An `unclear` novelty claim below has the concrete disposition **omit the proposed addition**; Claude need not research it to complete this integration. Such a disposition does not certify the unproved underlying claim.

All line numbers below refer to the reviewed merge. Match stable ids and literal text when integrating. The original report's claimed source commit is `9b912ee`, but many of its quoted passages and identifiers are absent from the reviewed merge. They must not be recreated merely to apply a replacement.

## Disposition ledger

| Input item | Verdict on original report | Certified integration action |
|---|---|---|
| O1 node labelling | Correct crossed SG pairing; wrong proposed repair and nonexistent target equations | Preserve Chapter 33; do not reverse topological charges or paste the replacement |
| O2 singular-wall limits | Correct numerical values; overprecise description of exceptional paths and false “initial guess typo” | Preserve the book's three rounded limits and two-path qualification |
| O3 particle index | Correct diagnosis; wrong wholesale replacement and invented Chapter 20 passage | Apply only the three-line index fix and README correction below |
| O4a folding parents | Correct that both parents are already listed; wrong quoted explanation/Coxeter number | No edit |
| O4b stale breather flags | Correct stale-output finding; wrong auxiliary coupling explanation | Correct eight flags in five output files, retaining their provenance |
| O4c thimble check | Correct dependency and sign observations | No edit; this is a source/output audit, not a new validation of the full thimble calculation |
| O5 Appendix E | Predominantly wrong source descriptions; blanket novelty certification unclear | Reject its deletion/replacement blocks and fourteen proposed summaries; use the companion reports only for statements they independently certify |
| O6 nearest-neighbour proof | Book Step 3 correct; original report's inequality derivation wrong; citation edit already present | No edit |
| O7a E8 masses | Correct | Add one numerical check and update its chapter/README description |
| O7b BK realization | Correct algebra | Add the checked numerical regression below, use generic parameters, and remove “not included in code” |
| O8 preface | Wrong quotations and wrong claim about `[Proposition]` | No edit |
| L5.2 exercise census | Staleness correct; final proposed count wrong | 84 exercises, all nine lectures written; preserve actual titles and links |

## O1. Crossed boundary labels: keep the present convention

**Independent evidence.** Put

\[
b_1(y)=y\sigma_++q^{-1}y^{-1}\sigma_-k_1+\hat\epsilon_1k_1,
\quad k_1=\operatorname{diag}(q,q^{-1}),
\]

and obtain node 0 by exchanging the two basis vectors and replacing the parameter by \(\hat\epsilon_0\). Direct symbolic multiplication gives

\[
K(y)b_j(y)-b_j(y^{-1})K(y)=0\quad(j=0,1),
\]

for the matrix displayed in the original report. `algebra_checks.py` checks every entry as a rational function, not at sample points. This verifies the algebraic matrix; it does not by itself fix the field-to-node dictionary.

For that dictionary, in sine-Gordon the positive chiral exponential contracts singularly with the negative boundary exponential. With \(\alpha_0=-\sqrt2\), the latter has coefficient \(C_0\). Thus \(Q_+\mapsto e_1\) implies \(\hat\epsilon_1\leftrightarrow C_0\); similarly \(\hat\epsilon_0\leftrightarrow C_1\). Combining \(\cos(\eta\mp i\vartheta)=\cos\eta\cosh\vartheta\pm i\sin\eta\sinh\vartheta\) with the stated UV–IR equations produces phases \(e^{\pm i\beta_{\rm SG}\phi_0/2}\), the same dictionary. No reversal of \(T_j\) is needed.

**Primary-source check.** I read Delius–MacKay, [hep-th/0112023v3](https://arxiv.org/pdf/hep-th/0112023v3), not the old referee's quotation. Their (3.20) has the terms \(\epsilon_-e^{i\hat\beta\tilde\phi/2}+\epsilon_+e^{-i\hat\beta\tilde\phi/2}\); (3.27) reads \(\hat Q_+=Q_++\bar Q_-+\hat\epsilon_+q^T\). In their general-rank (4.8) the exponent is **negative**, \(\sum_j\epsilon_j\exp(-i\hat\beta\alpha_j\cdot\tilde\phi/2)\). Their (4.11) reads \(\hat Q_i=Q_i+\bar Q_i+\hat\epsilon_iq^{T_i}\). These formulae explain why one must state the node dictionary rather than change just one sign in a field formula.

The reviewed `part7-half-line/33-boundary-quantum-groups.qmd:53` already explicitly distinguishes the Delius–MacKay node labels from the crossed SG labels. `eq-boundary-anomaly` and `eq-boundary-topological`, the purported repair targets, **do not exist** there. The report's sentence claiming a sign becomes correct “after \(T_j\to-T_j\)” is not a compatible local edit: it also changes \(k_j=q^{T_j}\) and the quantum-group dictionary. Its proposed antipode argument is not a derivation of the charge normalization either.

**Disposition:** no edit to Chapters 21, 30, 33, or 35 for this item; retain the explicit dictionary at `sec-boundary-charges` and `sec-sg-k-from-symmetry`. No literature correction or Appendix D/E change.

## O2. The three singular-wall limits are reproduced

**Independent construction.** `wall_limits.py` uses 90-digit arithmetic and no book or old-referee imports. From the Hirota subset formula, form

\[
\tau_j(x)=\sum_{S}\prod_{a\in S}z_a e^{m_ax+2\pi ija/h}
                 \prod_{a<b\in S}A_{ab}(0),\qquad
A_{ab}(v)=\frac{\cosh v-\cos(\pi(a-b)/h)}{\cosh v-\cos(\pi(a+b)/h)}.
\]

With rows \(\alpha_j=e_j-e_{j+1}\), set \(\phi=-\sum_j\alpha_j\log\tau_j\). The coefficient of an additional infinitesimal species-\(c\) soliton gives \(\delta\tau\), and the fluctuation is \(-\sum_j\alpha_j\delta\tau_j/\tau_j\). I differentiate this expression directly with `mpmath.diff`. The bulk and boundary Hessians are

\[
M^2=\sum_j\alpha_j\alpha_j^T\frac{\tau_{j-1}\tau_{j+1}}{\tau_j^2},\qquad
H_B=\frac12\sum_jC_j\alpha_j\alpha_j^T
                  \sqrt{\frac{\tau_{j-1}\tau_{j+1}}{\tau_j^2}}.
\]

At energy \(\Omega=m_1\cosh(0.6)\), retain outgoing open channels and decaying closed channels. Solve \(\delta\phi'(0)+H_B\delta\phi(0)=0\), with incoming species 1 coefficient normalized to one in the Fourier-polarization basis. The discarded dependent boundary row is checked too. Reversing the incident channel independently gives the same off-diagonal modulus, so the quoted $|K_{1,h-1}|$ and the computed $|K_{h-1,1}|$ agree. Every evaluated dressed wave satisfies the independent bulk Hessian equation to below \(10^{-65}\) at the test point; boundary residuals are below \(10^{-40}\).

Actual outputs (full convergence sequences in `wall_limits.log`):

| Background/path | Small parameter | \(K_{11}\) | Off-diagonal modulus |
|---|---:|---|---:|
| \(a_5\), \(\chi=5\pi/3\pm\delta\), regular stationary pairs | \(\delta=10^{-8}\), either side | \(0.1556328298110658+0.6771554884003561i\) | 0.7191932054839471 |
| Same endpoint, \(d=d_0(1-\epsilon)\), \(\chi=5\pi/3\) | \(\epsilon=10^{-10}\) | \(-0.06460440734771765+0.3421077965553655i\) | 0.9374372118106272 |
| \(a_4\), reduce \(d_1\), hold \(d_2\) fixed | \(\epsilon=10^{-9}\) | \(0.3253473726123012+0.9340770494629679i\) | 0.1471365107674794 |
| \(a_4\), reduce \(d_2\), hold \(d_1\) fixed | \(\epsilon=10^{-9}\) | \(0.3253473714343375+0.9340770502578904i\) | 0.1471365083257232 |

Solving \(\tau_2(0)=0\) and equality of the finite logarithmic derivatives gives
\(d_1=0.2577890250123870201380248\), \(d_2=0.5636698849843934158227122\).

**Verdict:** the three rounded numbers in `part7-half-line/32-half-line-solitons.qmd:215–218`, `sec-boundary-missing-charges`, are correct as numerical limits along the specified approaches. The original report incorrectly describes its exceptional check as equal two-sided finite-\(\epsilon\) values: its own printed values differ. The book actually says **two approaches**, and the two inward paths above verify that statement. This is not a proof of path independence at that endpoint or of a uniquely defined singular-background scattering theory. The book correctly retains that qualification.

The initial guess `0.2579` in `code/solitonic_linear.py` is not a typo requiring correction; an initial guess need not be the root rounded to four decimals. **Disposition:** preserve the numerical claims and code. No new Appendix D entry; the calculation does not invalidate a different singular-wall prescription.

## O3. Fix the index, reject the replacement script and Chapter 20 edit

`loop_check.py` independently expands \(V=\sum_jn_j(e^{\alpha_j\cdot\phi}-1)\) through fourth order in the mass basis. It evaluates

\[
\delta m_a^2=-\frac1{8\pi}\sum_{bc}C_{abc}^2
\int_0^1\frac{dz}{z m_b^2+(1-z)m_c^2-z(1-z)m_a^2},
\]

and the contact plus three exchange diagrams for the elastic tree amplitude. Its derivative of the block amplitude at \(B=0\) is analytic: because the brace is identically one there, the \(H'(0)\) and argument derivatives cancel. This avoids the original script's finite difference for the tree slope. Only the mass-ratio derivative uses a fourth-order symmetric difference.

For \(a_5^{(2)}\) with longest roots squared length **2**, the sorted eigenvalues are \((1,1.38196601125,3.61803398875)\). Particle 1 is index 1, whereas the lightest particle, index 0, is species 3. Actual results:

```text
c3:   deltaH = 0.07957747155; tree/B = 0.07957747155; ratio = 1.000000000019
a5tw: deltaH = 0.03978873577; tree/B = 0.07957747155; ratio = 0.500000000001
f4:   deltaH = 0.23873241464; tree/B = 0.07957747155; ratio = 2.999999999439
old a5tw index 0: tree/B = [0.06550170396, 0.05316746907, 0.05422744996, 0.05670419832]
old a5tw mean ratio = 0.693181070439
```

Thus the particle-index diagnosis is **correct**, including constancy of the corrected tree normalization. The report's statement that this computation uses longest roots squared length 4 is **wrong**. The real-coupling shift in this computation is positive, \(\delta H=\beta^2/(8\pi)\). Do not transplant its contradictory negative-sign explanation into a soliton chapter.

The original report's wholesale replacement is **wrong** independently of this sign problem. Its proposed \(c_3\) roots and marks have \(\sum_jn_j\alpha_j=(-1,1,2)\ne0\), so `toda_data` immediately asserts. It also changes squared mass-ratio formulas into unsquared ones and replaces the working four-block \(f_4\) amplitude by a two-block amplitude. Preserve those parts of the author's script.

**Exact integration packet:** in `code/rmatrix/smatrix/oneloop_ratio.py:25,37,52` make only these substitutions:

```python
def analyse(name, roots, kac, h, S11, ratio_formula, expected, p1=0):
```

```python
    tr = [tree_over_4m2sinh(m2, C3, C4, p1, th) for th in ths]
```

```python
analyse("a_5^(2)", a5, [1, 2, 1, 1], 5, lambda h, B, th: dgz11(h + B/2, B, th), a5ratio, "1/2 (H = h + B/2)", p1=1)
```

`snippet_validation.py` constructs precisely that candidate in scratch and runs all three original examples; they print 1.000000, 0.500000, 3.000000. This rerun is a secondary implementation check, following the independent calculation above.

Replace the **actual** “Known quirk” paragraph at `code/rmatrix/README.md:37` with:

```markdown
- `oneloop_ratio.py` uses mass-eigenvalue index 1 for particle 1 of $a_5^{(2)}$, because particle 3 is lighter. It gives $\delta H/B_{\rm tree}=1/2$, with a rapidity-independent tree normalization. The one-loop masses are also checked in `../../vernier/oneloop_particles.py`.
```

Update the existing `oneloop_ratio.py`, `oneloop_e6.py` list entry at line 71 to include `a_5^(2): 1/2`. There is no table row of the form quoted in the old report.

For an executable failure flag, Claude can additionally add after the printed ratio:

```python
    target = {"c_3^(1)": 1.0, "a_5^(2)": 0.5, "f_4^(1)": 3.0}[name]
    ok = np.ptp(dH) < 1e-6 and np.ptp(kap) < 1e-6 and abs(np.mean(dH)/np.mean(kap)-target) < 1e-4
    print(f"{'PASS' if ok else 'FAIL'} {name}: floating-Coxeter shift and tree shape")
```

The alleged sentence “The shifts of all mass ratios vanish …” is **absent** from the reviewed Chapter 20. Its actual folded-soliton table at line 129 lists \(a_{2n-1}^{(2)}\) as folded from \(d_{2n}^{(1)}\), with a nonzero shift. The proposed parent \(a_{2n-1}^{(1)}\) is not that table's parent. **Do not edit Chapter 20, Appendix D, or Appendix E for O3.**

## O4. Follow-ups already resolved, and stale output flags

### O4a: folding parents

**Correct** that `appendices/A-algebra-data.qmd:25` already lists both \(a_{2n}^{(1)}\) and \(d_{2n+2}^{(1)}\). **Wrong** that its explanation says the book's \(h=4n+2\): the row says \(h=2n+1\); the caption distinguishes the charge period \(kh\), here \(4n+2\). The quoted paragraph is absent.

Independent normalization check: in rank one the roots can be \(\alpha_0=-2e\), \(\alpha_1=e\), with coprime marks \((1,2)\). Their sum is zero and their mark sum is 3, not 6. This already refutes the report's general \(h=4n+2\) reading in the book's convention. **No edit.**

### O4b: breather flags

**Correct** about eight stale flags in five files. Independently, Takács–Watts [hep-th/9810006](https://arxiv.org/pdf/hep-th/9810006), equation (3.3), gives \(M_p^{(B)}=2M\sin(\pi p/(3\lambda))\), with levels bounded by \(3\lambda/2\). The first bound state is strictly below threshold for \(3\lambda>2\). With \(\lambda=(p'-p)/p\), this is \(3(p'-p)>2p\). At equality it is at threshold. The original report's unrelated \(\xi=\pi r/(p-r)\) argument does **not** derive this criterion and must not be inserted.

`scope_audit.py` recomputes every printed flag from its integer labels. Change only `breathers=yes` to `breathers=no` for:

| File under `code/irf/` | Lines | Models |
|---|---|---|
| `out_ha.txt` | 9, 11 | W3(5,7), W3(8,11) |
| `out_qsign/irf_ha_q1.txt` | 9, 11 | W3(5,7), W3(8,11) |
| `out_qsign/irf_ha_q-1.txt` | 9, 11 | W3(5,7), W3(8,11) |
| `out_qsign/irf_tm2_q1.txt` | 5 | W3(8,11) |
| `out_qsign/irf_tm2_q-1.txt` | 5 | W3(8,11) |

These are corrected annotations of stored output, not fresh runs; record that in the integration record. Alternatively regenerate with the scripts' original commands, but regeneration is not needed to verify these integer labels. The scripts and Chapter 28 already use the correct inequality. No change to their formulas or status labels.

### O4c: thimble check

**Correct as a source/output audit.** `code/thimbles/check_th4.py:8–9` checks for its critical-point file and invokes `cps_a2L2.py 1.0 1` when absent. Its import list has no `check_th2` import. In `out_check_th4.txt` the fine-grid residual/prediction ratios computed from the stored numbers are \((-3.0976\times10^{-7})/(-3.0315\times10^{-7})=1.021804\) and \((1.7710\times10^{-10})/(2.0053\times10^{-10})=0.883160\). Thus the alleged sign mismatch is absent. This verifies the follow-up's actual claims, not the full thimble algorithm or an exact equality of those numbers. **No edit.**

## O5. Reject the Appendix E replacement block

The old report's six-row audit is not a reliable account of the source. `scope_audit.py` independently resolves ids and prints each target's actual statement. The frozen source and its SHA256 hashes are recorded in `scope_audit.txt`.

### O5a, O5c: purported existing errors

- The existing folded-mass bullet is **[Numerical]**, describing computed folded soliton masses. There is no **[Established]** BCDS correction bullet to delete. Moreover, AGENTS.md explicitly allows original literature corrections in Appendix E; presence in Appendix D is not by itself a reason for deletion.
- The proposed `prp-cn-smatrix`, `prp-a2n-smatrix`, `prp-e6-smatrix`, `prp-f4-r`, `prp-f4-breathers`, `prp-dn-vector-r`, and `prp-sector-table` targets do not exist. The proposed new references would fail the reference check. The source's \(f_4\) bullet is an obstruction for a **27-dimensional** candidate, not a constructed 26-dimensional \(f_4\) solution; the 26-dimensional statement concerns \(e_6^{(2)}\).
- The actual breather reference belongs to a combined **[Theorem]** bullet whose principal targets are `prp-sgd-scalar-integral` and `prp-sgd-pole-ladder`. It is not the standalone “Theorem for B1” bullet the report quotes. The new quantum report supplies the exact separation/label disposition; do not paste the old bullet.
- The original proposal puts `@prp-bloch-sectors` inside a status-label bracket. Apart from misdescribing the twist defect, that conflicts with the repository's label syntax rule.

### O5b: all fourteen proposed original-result summaries

Here “wrong” concerns the claimed summary/scope, not a newly alleged error in the existing theorem. “Unclear novelty” means no independent prior-art evidence is provided by the old report; its claim that the chapter calls something new is not proof of novelty. **The executable disposition for each is to omit the old proposed bullet.** The companion reports may supply a different, independently checked bullet for their own covered statements.

| Statement | Actual target at reviewed commit | Independent audit and verdict on old summary |
|---|---|---|
| `prp-soliton-preserving-half` | `31-boundary-conditions.qmd:260` | Actual theorem concerns compatibility and necessity of boundary coweights; it expressly leaves involution unproved. Half-soliton class counting is in following prose. Old combined theorem summary overstates the theorem's scope; blanket novelty unclear. |
| `prp-sp-scalar` | `34-soliton-reflection.qmd:348` | Numerical cases are explicitly \((n,l)=(1,1),(2,0),(2,1),(3,0),(3,1),(3,2)\), with bootstrap for \(n\le3\). Old all-species/all-rank wording is unsupported. |
| `prp-dn-uniform-vacuum` | `36-particle-reflection.qmd:304` | This is an existence/uniqueness theorem for a static real-coupling vacuum, with explicit d4/d5 kinks. It is not a diagonal vector-soliton matrix or a four-pole theorem. Wrong. |
| `prp-dn-uniform-classical` | `36-particle-reflection.qmd:328` | Numerical classical matching is limited to d4,d5,d6; cannot say every rank. The actual Notes identify classical vacuum matching as the new part, while citing the amplitude to Castro-Alvaredo–Fring. Old unrestricted summary wrong; independent novelty not certified here. |
| `prp-dn-nondiagonal-scalar` | `34-soliton-reflection.qmd:261` | Explicitly numerical for n=4,5,6, not every n or every channel at arbitrary rank. Old scope wrong. |
| `prp-twist-defect` | `37-defect-scattering.qmd:237` | Defined by the character \(e^{i\alpha\cdot\mu}\) on conserved particle charges and one-dimensional defect space. It is not the proposed Bloch-projector operator with the quoted status. Wrong description. |
| `prp-sg-spin-three` | `39-defect-conditions.qmd:231` | Both spin signs are conserved with defect terms. Old “only at trivial sigma” restriction contradicts it; the classical companion independently checks both jet identities. Wrong. |
| `prp-dc-gauge` | `39-defect-conditions.qmd:318` | A differential gauge transformation of Lax pairs, not affine-Weyl equivalence of n! permutations. Wrong; use classical report's actual proof/upgrade. |
| `prp-defect-theta` | `39-defect-conditions.qmd:437` | Symmetry holds iff sigma is real, as follows termwise from conjugating exponentials. The broad description is compatible, but old location/quotation is wrong and novelty was not independently established. Omit old proposed bullet. |
| `prp-ds-absorption` | `40-defect-solitons.qmd:104` | Species/class-dependent rapidity and a separate nonlocalized case. For species 1 of a5, class 1, the text explicitly gives real absorption rapidity. Thus universal theta_d−iπ/2 is false. |
| `prp-defect-frozen-soliton` | `40-defect-solitons.qmd:181` | A pair of opposite-type defects reproduces a soliton delay, up to a species phase. Same-type pairs fail in rank≥2. This is not identification of each static step with a frozen soliton. Wrong. |
| `prp-fused-bdg` | `43-defects-and-boundaries.qmd:195` | Dressing conservation plus comparison which explicitly says the Neumann fusion is **not** the BDG boundary. Old summary obscures this essential distinction. Wrong as supplied. The quantum companion covers `prp-defect-fusion-k`; it does not certify all of this distinct `prp-fused-bdg` statement. |
| `prp-an-particle-transmission` | `44-particle-transmission.qmd:75` | Algebraic properties of candidate particle factors and classical matching; physical identification depends on `cnj-transmission-continuation`. Old wording makes the candidate identity sound like an established physical quantum bootstrap. Unsupported. |
| `prp-dpt-vacuum` | `44-particle-transmission.qmd:138` | Positivity/unique real-coupling vacuum at sigma>0, unbounded potential at sigma<0. It is not a dressed transmission matrix for h/gcd static steps. Wrong. |

This audit rejects specific false prescriptions; it is not a certification that Appendix E is complete or that all fourteen underlying mathematical statements have been newly proved in this pass. **No blanket Appendix E insertion or deletion is cleared by O5.** This conservative resolution requires no additional research by the integrator.

## O6. Nearest-neighbour cluster estimate: book correct, report proof wrong

Target: `part9-foundations/48-continuum-window.qmd:66–109`, Step 3 of `thm-nn-distances`.

**Independent derivation.** Fix a labelled minimum spanning tree. Replacing the longest edge on its path by a shorter direct edge would lower the tree weight, so its path bottleneck is bounded by the Euclidean endpoint distance. Consequently, just before a leaf's first cluster of size at least m+1 forms, at most m−1 neighbours lie at distance below that merger scale. This proves \(\tilde\delta_k\ge s_{i(k)}\).

For a fixed compatible dendrogram, a tree edge merges the two child clusters only after their internal edges. Conversely any ordering respecting those child-before-parent relations produces precisely that dendrogram, by induction over mergers. Thus compatible orderings are its linear extensions. Interleaving the child orders yields the hook count \((n-1)!/\prod h_j\) by induction. Partitioning all edge orders gives \(\sum_{T\sim\tau}\prod h_j^{-1}=1\).

In edge-vector coordinates, removing leaves successively gives Jacobian of absolute value one. Polar integration contributes \((2\pi)^{n-1}\). Integrating scale variables from children to parents gives denominators \(E_j=2h_j-\gamma C_j\). Crucially,

\[
C_j=0\quad(h_j<m),\qquad C_j=h_j+1\quad(h_j\ge m).
\]

For the second case, with \(\alpha=2m-(m+1)\gamma>0\),

\[
\frac{E_j}{2h_j}=1-\frac\gamma2\Big(1+\frac1{h_j}\Big)
\ge1-\frac\gamma2\Big(1+\frac1m\Big)=\frac\alpha{2m}.
\]

For the first case the ratio is one. This proves the bound used in the book. Summing dendrograms gives \((2\pi m/\alpha)^{n-1}\) per labelled tree, and Cayley's count completes \(Z_n\le A^nn!\). No quadrature or random configurations are needed for this proof.

**Error in the old report.** Its asserted \(1-\gamma\ge\alpha/(2m)\) is false: m=2, alpha=1 gives gamma=1 and 0≥1/4. It loses the decisive implication \(C_j>0\Rightarrow h_j\ge m\). Its claimed algebraic difference has the wrong sign. The book uses the correct implication already.

`algebra_checks.py` independently enumerates all 145 labelled trees for 2≤n≤5 and every edge ordering, verifying each hook count, the sum identity, and the exponent bound. It also integrates the m=2, gamma=6/5 caterpillar exactly, obtaining 25/24. These support, but do not replace, the proof above.

**Disposition:** preserve Step 3 and its theorem label. The actual checks paragraph at line 109 already names `code/review614/step3_combinatorics.py`, with more detail than the proposed replacement. No edit.

## O7. Two genuinely missing regression checks

### O7a. E8 mass ratios

**Correct.** My independent Hessian uses explicit Bourbaki E8 roots with marks (2,3,4,6,5,4,3,2), appending the negative highest root with mark 1. The affine root has squared length 2 and the weighted root sum is zero. Diagonalizing \(\sum n_j\alpha_j\alpha_j^T\) gives ratios

```text
1, 1.61803399, 1.98904379, 2.40486717,
2.95629520, 3.21834046, 3.89115682, 4.78338612
maximum error against the displayed trigonometric expressions: 2.6645352591003757e-15
```

The current `code/algebra_data.py:50–52` prints the E8 masses but does not test that row. Add immediately after the e6/e7/e8 print loop:

```python
_, _, ms_e8 = masses('e', 8)
c5, c30, c7, c2 = (np.cos(np.pi/5), np.cos(np.pi/30),
                   np.cos(7*np.pi/30), np.cos(2*np.pi/15))
r_e8 = np.array([1, 2*c5, 2*c30, 4*c5*c7, 4*c5*c2,
                 4*c5*c30, 8*c5*c5*c7, 8*c5*c5*c2])
report('e_8: closed-form mass ratios of sec-masses-couplings',
       np.allclose(ms_e8/ms_e8[0], r_e8, rtol=1e-11, atol=1e-12))
```

At `part3-classical/09-lagrangian.qmd:118`, replace the sentence beginning `` `code/algebra_data.py` checks `` with:

```markdown
`code/algebra_data.py` checks all three rows of the table ($a_n$ and $d_n$ for $n\le8$, and $e_8$) and the folded masses below.
```

In `code/README.md:18`, append “including the closed E8 mass ratios” to the script's description. No new theorem or original-result claim is involved.

### O7b. Baseilhac–Koizumi realization

**Correct, with stronger independent evidence.** I read [hep-th/0208005v2](https://arxiv.org/pdf/hep-th/0208005v2), printed page 7. Its (19) is exactly
\([\hat E_+\hat E_-,\hat E_-\hat E_+]=c^2(q+q^{-1})^2(\hat E_+^2-\hat E_-^2)\).
Its (20) gives \(\hat E_+=\pm2e_{UV}\cosh\hat p\), \(\hat E_-=\pm2e_{UV}\cosh\hat q\), \(\pm q=e^{\alpha/2}\), with the exact normalization line “with the normalization” \(e_{UV}^2=-c^2/(q-q^{-1})^2\).

Write \(UV=q^2VU\), \(f=ic/(q-q^{-1})\), \(E_+=f(aU+a^{-1}U^{-1})\), \(E_-=f(bV+b^{-1}V^{-1})\). Normal ordering obeys
\((U^rV^s)(U^tV^u)=q^{-2st}U^{r+t}V^{s+u}\).
`algebra_checks.py` expands both sides in this Laurent basis with **symbolic** q,c,a,b. Every coefficient of (19), and of both q-Dolan–Grady differences with rho=c²(q+q⁻¹)², is exactly zero. This does not use a finite-dimensional approximation or assume a root of unity.

For the production regression, the original finite-band method is valid: each degree-four word shifts by at most four, so keeping interior input columns leaves all output rows exact. `snippet_validation.py` independently evaluates this insertion and finds residuals \(7.88\times10^{-16}\), \(1.56\times10^{-15}\). Use an irrational physical parameter rather than the existing script's short-decimal omega:

```python
# Baseilhac-Koizumi (19) and q-Dolan-Grady on a Heisenberg module.
worst_bk = 0.0
N_bk = 15
keep_bk = slice(4, N_bk - 4)
U_bk = np.diag(np.ones(N_bk - 1, complex), k=-1)
ks_bk = np.arange(N_bk) - N_bk // 2
for q_bk in (-np.exp(-1j*np.pi*np.sqrt(2)), QG):
    V_bk = np.diag(q_bk**(-2.0*ks_bk))
    Vi_bk = np.diag(q_bk**(2.0*ks_bk))
    c_bk, a_bk, b_bk = 0.9-0.3j, 1.3-0.4j, 0.8+0.7j
    f_bk = 1j*c_bk/(q_bk-1/q_bk)
    Ep_bk = f_bk*(a_bk*U_bk + U_bk.T/a_bk)
    Em_bk = f_bk*(b_bk*V_bk + Vi_bk/b_bk)
    rho_bk = c_bk**2*(q_bk+1/q_bk)**2
    pairs_bk = [
        (Ep_bk@Em_bk@Em_bk@Ep_bk-Em_bk@Ep_bk@Ep_bk@Em_bk,
         rho_bk*(Ep_bk@Ep_bk-Em_bk@Em_bk)),
        (serre(Ep_bk, Em_bk, -2, q_bk), rho_bk*(Ep_bk@Em_bk-Em_bk@Ep_bk)),
        (serre(Em_bk, Ep_bk, -2, q_bk), rho_bk*(Em_bk@Ep_bk-Ep_bk@Em_bk))]
    for lhs_bk, rhs_bk in pairs_bk:
        worst_bk = max(worst_bk, rel(lhs_bk[:, keep_bk], rhs_bk[:, keep_bk]))
report(f'BK Heisenberg realization: (19) and q-Dolan-Grady [{worst_bk:.1e}]',
       worst_bk < 1e-12)
```

Append to `code/qsp_checks.py`, using its existing `serre`, `rel`, `report`, `QG`. In the **Dynamical boundaries** paragraph of `part7-half-line/35-quantum-symmetric-pairs.qmd:206`, replace exactly ``(a numerical check, not included in `code`)`` by `` (`code/qsp_checks.py`) ``. Add “BK Heisenberg realization of (19) and q-Dolan–Grady” to `code/README.md:28`. No Appendix D correction: the source is correct. No theorem upgrade is needed.

## O8. Preface: the proposed replacement has no target

**Wrong.** Neither quoted sentence appears in `preface.qmd` at the reviewed commit. The claim that `[Proposition]` is absent from `status-tags.lua` is directly false: `Proposition = true` is in its table. The Introduction at `sec-new-results` lists **[Proposition, sketch]**. The proposed `[Open]` and `[Refuted]` labels are not entries in the filter.

`scope_audit.py` checks the purported replacement anchor and the filter entries independently. **Disposition:** no preface edit, no new status labels, no status-filter changes for this item.

## L5.2. The Claude course has 84 exercises

**Correct** that the homepage is stale. **Wrong** to use the proposed 83 without counting Lecture 9. Counting actual `::: {#exr-...}` blocks gives 8,7,9,10,9,10,10,10,**11**, total **84**. All nine files are in the project. `lecture-notes/defects/index.qmd` already lists all nine and has no “Lecture 9 in preparation” sentence. No edit to that index is needed.

In `index.qmd:44`, retain the actual title, published URL, PDF URL, and “drafted by Claude”. Replace just the final two sentences about preparation by:

```markdown
All nine lectures are written. Lectures 1–6 cover the classical theory, and Lectures 7–9 cover quantum transmission, Q-operators and boundary fusion.
```

In the comparison table at `sec-about-differences` make these literal changes in the Claude column only:

- `Lectures 1–6 of 9` → `All 9 lectures`.
- Replace the Scope cell by: `The classical theory in Lectures 1–6: point defects, Bäcklund transformations, defect conditions, solitons, Lax matrices and r-matrices. The quantum theory in Lectures 7–9: transmission matrices, Q-operators and boundary fusion`.
- `53 so far, of three kinds: redo, rediscover and explore` → `84, of three kinds: redo, rediscover and explore`.

Do not import the original report's invented course titles, model version credits, estimated page counts, or Gemini filenames. The Gemini report/course is outside this round's scope.

## Integration

This section records verification **before** integration; no book edits have been made.

| Classification | Items | Basis |
|---|---|---|
| Verified, no edit | O1 crossed dictionary, O2 rounded limits, O4a existing parents, O4c dependency/sign audit, O6 actual book proof | Symbolic intertwining/source equations; new high-precision Hirota calculation; source inspection; arithmetic; full derivation and exact enumeration |
| Modified, ready to integrate | O3 index/README, O4b flags, O7a check, O7b check, L5.2 homepage | Fresh tensor calculation and candidate execution; integer criterion; independent E8 roots; exact symbolic quantum-torus identities; complete exercise census |
| Rejected prescriptions | O1 sign flip, O2 initial-guess “typo”, O3 wholesale script and Chapter 20 text, O4a fictitious explanation, O5 old Appendix E blocks, O6 citation replacement, O8 preface text | Explicit contradictions, absent anchors, incorrect algebra, or already-completed work |
| Not certified | O5 blanket novelty/completeness, singular-wall path independence beyond tested paths, full thimble numerical calculation | Not established by the old report's evidence; none is needed for the cleared edits |

After applying the cleared edits, run the affected chapter scripts from their normal directories, `python3 check-refs.py`, and the applicable lecture reference checks. Run scratch renders when preparing a PR. Record which certified actions were implemented and any content conflicts caused by later changes. Those are integration checks, not a request to re-derive the results again.

## Reproduction

With an untouched source tree at the reviewed commit:

```bash
python3 briefs/pr60_open_verified_scripts/algebra_checks.py
python3 briefs/pr60_open_verified_scripts/loop_check.py
python3 briefs/pr60_open_verified_scripts/wall_limits.py
python3 briefs/pr60_open_verified_scripts/scope_audit.py /path/to/reviewed/tree
```

The archived logs are actual output. `snippet_validation.py` is a secondary candidate-integration check and writes a new `candidate_oneloop_ratio.py` beside itself; run its scratch copy if preserving an evidence archive read-only. It is not the independent verification. Source PDFs are not part of the script archive; their versioned links and exact formula locations are given above.
