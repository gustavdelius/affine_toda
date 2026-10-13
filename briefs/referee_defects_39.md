# Brief: deep referee pass on chapter 39 (solitons and particles through a defect)

You are an adversarial mathematical-physics referee. You are working in the repository of a Quarto book on affine Toda field theory (worktree `.claude/worktrees/integrable-defects-framework-a54a84`, branch `claude/integrable-defects-framework-a54a84`).

**Rules.** Other agents are editing files in this worktree at the same time.
- Do not edit any existing file.
- Do not run git commands that change the index or the working tree (no commit, checkout, stash, reset, rebase or restore).
- Do not run `quarto render`.
- Write scratch scripts only under `referee_scratch/defects_39/`.
- Write the report to `briefs/referee_defects_39_report.md`.

## Context

Chapter 39 (`part8-defects/39-defect-solitons.qmd`, `sec-defect-solitons`) was drafted on 2026-10-12 by a Claude agent. Its verification script is `textbook_code/defect_solitons_an.py` (9 checks). The drafting report is `briefs/draft_ch39_report.md`.

The chapter treats classical $a_r$ solitons and particles meeting a type-I defect:
- jump classes $p\in\Lambda_W/\Lambda_R$, and the soliton delay for both types and all classes (`eq-an-defect-delay`, `prp-an-defect-delay`), built on ch. 38's `eq-dc-ar-delay`;
- the stored topological charge and the selection rule (`prp-stored-charge`);
- absorption and emission (`prp-ds-absorption` [Theorem]);
- several solitons; reality;
- linear transmission factors (`eq-an-linear-transmission`);
- the frozen-soliton comparison (`prp-defect-frozen-soliton` [Theorem]).

Read these first:
- `briefs/part8_defects_conventions.md`: the conventions, verified by script, and the source dictionary.
- The chapter-39 section of `briefs/part8_defects_status.md`.
- `briefs/part8_defects_sources.md`: digests of the sources, with their equation numbers.
- Chapter 38 (`part8-defects/38-defect-conditions.qmd`), whose delay formula chapter 39 uses.
- Source texts in `~/referee-copy/papers/*.txt`, notably 0705.1066 (cz2007), 0902.1307 (cz2009a), hep-th_0401020 (bcz2004b) and hep-th_0305022.

To find an object: `grep -rn '{#<id>}' part*/ index.qmd`.

For each item, do **one independent computation**: a derivation written out, or a script of your own. You may run the book's scripts for comparison, but they are not evidence. A verdict without its own evidence counts as "unclear". `briefs/referee_defects_36_37_deep_report.md` shows the expected depth.

## Items (most important first)

1. **The type-I₊ delays and the jump classes.**
   - The chapter derives type-I₊ $a_r$ delays $z_a(\theta_d-2\pi ip/h-\theta)$ by time reversal. Type I₋ is $z_a(\theta-\theta_d-2\pi ip/h)$.
   - It claims a class-$p$ jump acts like a jump-free defect at $\theta_d\pm2\pi ip/h$, because every $f_j$ picks up $\omega^{-p}$.
   - Check both claims with exact $a_2$ and $a_3$ soliton solutions on both sides of the defect, using your own τ-functions in book normalization (`eq-an-soliton`, `eq-defect-conditions`).
2. **`prp-ds-absorption` [Theorem].** Check each part:
   - absorption at $\theta_d\mp i(\pi/2-\pi(a+2p)/h)$;
   - the orientation rule: the absorbed or emitted envelope moves with velocity $\tanh\theta_d$;
   - absorption at real rapidity iff $a+2p\equiv h/2$ (mod $2h$). Test the example "species 1 of $a_5$ is absorbed at $\theta_d$ by a class-1 defect";
   - the energy and momentum stored by a class-$a$ jump.

   Is the [Theorem] label justified by the proof given?
3. **`prp-stored-charge` and the selection rule.**
   - The class energies lie on the mass shell $2hm/\beta^2$, and they are real iff $2p\equiv0$ (mod $h$).
   - A passing soliton changes the jump by $\lambda_{\rm in}-\lambda_{\rm out}$.
   - The outgoing charge moves at most $a/\gcd(a,h)$ steps along the π-orbit, and both bounds are attained.
   - The defect never produces the missing charges (cf. `sec-missing-charges`, ch. 11).
4. **`prp-defect-frozen-soliton` [Theorem].** A type-I₋ defect at $\theta_c+i(\pi/2+\pi c/h)$ paired with a type-I₊ defect at the complex-conjugate parameter gives $\omega^bA_{cb}$ and $\omega^bX_{cb}$ for every species $b$, and carries $M_c(\cosh\theta_c,\sinh\theta_c)$.
   - Check this independently against `eq-interaction-coefficient` (ch. 10) and `eq-an-transmission` (ch. 16).
   - Check the statement that two defects of the same type reproduce $A_{cb}$ only for $b=h/2$.
   - Check the realizability condition $h+2c\equiv0$ (mod 4).
5. **Linear transmission.** Check $T_k=z_k(\theta-\theta_d+i\pi/2)$, $R=0$, unitarity, and crossing $\bar T_k(\theta)=T_{\bar k}(i\pi-\theta)$. Also check the claimed cz2007 (7.11) misprint: the denominator should be $ie^{-\theta}m_k+\sigma(1-\omega^k)$, with an implicit $m=1$. Work in their own conventions and confirm that the corrected form reproduces their (7.5).
6. **The other corrections to sources.**
   - cz2007 §2: the charge-sector boundaries are $2\pi\gcd(a,h)/h$ apart, not $2\gamma_a$.
   - "Only self-conjugate solitons are absorbed at real rapidity" holds only for root-lattice jumps.
   - cz2009a after (4.8): the static energies are also real for $h$ even when the classes differ by $h/2$.
7. **The θ_d → −∞ limit.** Is it true that type I₋ becomes the invertible line of the diagram rotation π ($\psi=\pi\phi$, $z\to\omega^a$)?
8. **Reality.** Θ maps class $p$ to class $-p$ and conjugates the delays. In charged classes with $\cos(2\pi p/h)<0$, solitons are advanced rather than delayed.
9. **Status labels, exercises and bibliography.** Spot-check at least three exercises.

## Report format (`briefs/referee_defects_39_report.md`)

A header with the date, the model and the commit. Then one section per item, with:
- the verdict (correct / wrong / unclear);
- your own computation (the script inline or its path, and its output);
- for errors, corrected text ready to paste, with file and line.

End with other observations, labelled as unverified.
