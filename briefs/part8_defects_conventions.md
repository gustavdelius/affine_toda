# Part VIII (defects): conventions, names and id registry

This file is the single source of truth for the conventions of Part VIII, "Affine Toda theory with defects" (`part8-defects/`, chapters 36–43). It was written on 2026-10-11, before drafting. Every agent drafting or refereeing a chapter of Part VIII reads it first. Where a convention is marked *verified*, the named check in `textbook_code/` derives it; where it is marked *to verify*, the first agent that needs it confirms it by script and reports.

The plan of the part is recorded in the memory note `defects-part` and in the plan file `~/.claude/plans/system-reminder-you-are-operating-indexed-stroustrup.md`.

## Placement and numbering

- Part VIII sits between Part VII (the half-line, ch. 28–35) and Foundations. When the part is complete, Foundations becomes Part IX (`part9-foundations/`, ch. 44–49) and the outlook becomes chapter 50; until then the old numbers remain.
- Chapters mirror Part VII: 36 ↔ 28 (axioms), 37 ↔ 29 (rank one), 38 ↔ 30 (classical conditions), 39 ↔ 31 (solitons), 40 ↔ 32 (quantum groups), 41 ↔ 33 (soliton matrices), 42 ↔ 34 (quantum symmetric pairs), 43 ↔ 35 (particles and open problems).

## Names and symbols (clash table)

| Use in Part VIII | Symbol / word | Why not the obvious alternative |
|---|---|---|
| Transmission matrix, left to right | $\mathbb T_\mu(\theta-\theta_d)$ | $T$ is already the topological charges $T_j$, the gradation $T$ and transfer matrices $T(\lambda)$; $\mathsf T$ is the transpose; $\mathcal T$ is the double-row matrix (`eq-double-row`) |
| Transmission matrix, right to left | $\bar{\mathbb T}_\mu$ | |
| Scalar particle transmission factor | $\mathbb T_a(\theta)$ | $X_{ab}$ is the fluctuation transmission factor of Part V (`eq-an-transmission`) |
| Defect state space | $\mathcal W$ | $\mathcal D$ and $D$ are used for differential operators and Hirota derivatives |
| Fields left / right of the defect | $\vec\phi$ on $x<0$, $\vec\psi$ on $x>0$ | $u,v$: $u=i\theta$ throughout Part VII; $\phi_L$, $\phi_R$: $\phi_L$ is the chiral field of ch. 20 |
| Defect potential | $\mathcal B_d(\phi,\psi)$ | it becomes a boundary potential under the doubling trick (cf. $\mathcal B$ in `eq-boundary-action`) |
| Defect matrix in the defect conditions | $\mathsf D$ (sans serif), $E=1-\mathsf D$ | plain $D$ is overloaded; sans serif follows $\mathsf M^2$ |
| Defect parameter | $\sigma=e^{-\theta_d}$, $\theta_d$ the defect rapidity | $\eta$, $\vartheta$ are the BPTT boundary parameters ($\vartheta$ is also the boost) |
| Classical defect matrix | $L_d(\lambda)$ | $K$ is Part VII's reflection matrix |
| Defect types | type I$_+$, type I$_-$ (keep $Q_{+2}$ resp. $Q_{-2}$), type II | |
| Defect → boundary of the doubled theory | "the doubling trick" (`sec-doubling-trick`); footnote that the literature says "folding trick" | "folding" means Dynkin-diagram folding throughout the book (`sec-affine-folding`, `sec-images-cn`) and Robertson's folding of defects (`sec-defect-folding`) |
| Delfino–Mussardo–Simonetti | bib key `delfino_mussardo_simonetti1994` | `dms2001` is Delius–MacKay–Short |

## Conventions

*To be completed from the conventions agent's report (phase 0).* Topics:

1. Geometry and action, the defect conditions (`eq-defect-action`, `eq-defect-conditions`).
2. The matrix $\mathsf D$ for $a_r$ (`eq-d-matrix`).
3. The defect potential in book normalization, real and imaginary coupling, and sine-Gordon (`eq-defect-potential`, `eq-sg-defect-potential`).
4. Defect rapidity, boosts, parity.
5. Classical defect matrix $L_d(\lambda)$ and the types I$_\pm$ (`eq-defect-lax-matrix`).
6. Quantum: the Part VI coproduct, $q=q_{\rm phys}$, the principal-gradation variable, $\mathbb T:V_\mu(\theta)\otimes\mathcal W\to\mathcal W\otimes V_\mu(\theta)$, which Borel subalgebra goes with which type.
7. Source dictionary `tbl-defect-dictionary`.

## Status-label rules specific to Part VIII

- Every quantum $a_r$ transmission matrix is labelled as a solution of the scattering axioms ([Numerical] or [Established]); its identification with the transmission matrix of the Lagrangian defect is [Conjecture], because of the higher-spin anomalies found by Penati and Polvara (arXiv:1902.10690).
- `cnj-twisted-diagonal` is pre-registered: the expectation, stated in the text before the check, is that a jump defect preserves a diagonal *Borel* subalgebra with a dynamical boundary module $\mathcal W$, not the full twisted diagonal coideal. The script `defect_doubling.py` decides; the text reports the outcome either way.
- The sign of $q$ is tested by unitarity, crossing and the scalar factor, not by the triangle relations alone.

## Id registry

All ids are already present as placeholders in the skeleton chapters, so chapters drafted in parallel can cross-reference each other. Add new ids only inside your own chapter, prefixed so that they cannot clash, and report them.

| Ch | File | Ids |
|---|---|---|
| 36 | `part8-defects/36-defect-scattering.qmd` | `sec-defect-scattering`, `sec-defect-kinematics`, `eq-defect-zf`, `prp-defect-elastic`, `sec-rt-algebra`, `prp-rt-nogo`, `sec-triangle-relations`, `eq-triangle-st`, `eq-triangle-braid`, `eq-triangle-rtt`, `prp-triangle-forms`, `sec-defect-unitarity-crossing`, `eq-defect-unitarity`, `eq-defect-crossing`, `sec-doubling-trick`, `eq-doubled-k`, `prp-doubling`, `sec-defect-bound-states`, `eq-defect-bootstrap`, `sec-defect-circle`, `prp-twist-defect`, `sec-topological-lines`, `sec-defect-finite-size`, `sec-defect-lee-yang`, `eq-ly-transmission`, `prp-ly-transmission`, `sec-defect-scattering-exercises`, `sec-defect-scattering-notes` |
| 37 | `part8-defects/37-defect-sine-gordon.qmd` | `sec-defect-sine-gordon`, `sec-backlund`, `eq-sg-backlund`, `prp-sg-backlund`, `prp-bianchi`, `sec-sg-defect-lagrangian`, `eq-sg-defect-lagrangian`, `eq-sg-defect-potential`, `eq-sg-defect-conditions`, `prp-sg-defect-classification`, `prp-sg-defect-momentum`, `sec-sg-defect-solitons`, `eq-sg-delay`, `prp-sg-delay`, `prp-half-soliton`, `sec-sg-defect-linear`, `eq-shg-linear-transmission`, `sec-type-two`, `sec-sg-transmission-matrix`, `eq-sg-transmission-matrix`, `prp-sg-transmission`, `cnj-defect-parameter`, `sec-sg-defect-spectrum`, `sec-sg-breather-transmission`, `prp-sg-breather-transmission`, `sec-defect-sine-gordon-exercises`, `sec-defect-sine-gordon-notes` |
| 38 | `part8-defects/38-defect-conditions.qmd` | `sec-defect-conditions`, `sec-defect-lagrangian`, `eq-defect-action`, `eq-defect-conditions`, `eq-defect-potential`, `eq-d-matrix`, `sec-defect-momentum`, `prp-defect-momentum`, `prp-ar-only`, `cnj-momentum-integrability`, `sec-defect-lax`, `eq-defect-lax-matrix`, `prp-defect-lax`, `prp-spin-two`, `sec-defect-r-matrix`, `eq-defect-monodromy`, `prp-defect-involution`, `sec-defect-backlund-an`, `sec-type-two-atft`, `prp-type-ii-momentum`, `sec-defect-folding`, `sec-defect-theta`, `prp-defect-theta`, `sec-defect-conditions-exercises`, `sec-defect-conditions-notes` |
| 39 | `part8-defects/39-defect-solitons.qmd` | `sec-defect-solitons`, `sec-an-soliton-defect`, `eq-an-defect-delay`, `prp-an-defect-delay`, `sec-defect-topological-charge`, `prp-stored-charge`, `sec-defect-absorption`, `sec-defect-multi`, `sec-defect-reality`, `sec-classical-transmission-factors`, `eq-an-linear-transmission`, `prp-an-linear-transmission`, `prp-defect-frozen-soliton`, `sec-defect-solitons-exercises`, `sec-defect-solitons-notes` |
| 40 | `part8-defects/40-defect-quantum-groups.qmd` | `sec-defect-qg`, `sec-cft-topological-defects`, `sec-perturbed-defects`, `sec-borel`, `eq-q-oscillator`, `eq-weyl-module`, `eq-fock-module`, `sec-t-intertwiner`, `eq-borel-intertwiner`, `thm-t-triangle`, `prp-sg-t-unique`, `prp-two-borels`, `sec-universal-r-affine`, `eq-l-operator`, `prp-l-operator`, `sec-q-operators`, `prp-defect-transfer`, `cnj-defect-tq`, `prp-defect-twist-limit`, `sec-defect-boundary-dictionary`, `sec-defect-qg-exercises`, `sec-defect-qg-notes` |
| 41 | `part8-defects/41-soliton-transmission.qmd` | `sec-soliton-transmission`, `sec-t-procedure`, `sec-an-soliton-transmission`, `prp-an-transmission`, `sec-transmission-unitarity-crossing`, `prp-transmission-crossing`, `sec-a22-transmission`, `prp-a22-transmission`, `sec-an-defect-bound-states`, `prp-a2-unstable`, `sec-transmission-semiclassics`, `prp-transmission-semiclassical`, `sec-soliton-transmission-exercises`, `sec-soliton-transmission-notes` |
| 42 | `part8-defects/42-defects-and-boundaries.qmd` | `sec-defect-boundary`, `sec-doubled-theory`, `prp-reflected-copy`, `sec-diagonal-qsp`, `prp-diagonal-coideal`, `sec-jump-coideal`, `cnj-twisted-diagonal`, `prp-diagonal-borel`, `sec-defect-boundary-fusion`, `prp-defect-fusion-k`, `prp-fused-bdg`, `sec-tensor-k`, `sec-noninvertible-boundaries`, `sec-defect-boundary-exercises`, `sec-defect-boundary-notes` |
| 43 | `part8-defects/43-particle-transmission.qmd` | `sec-particle-transmission`, `sec-breather-transmission`, `prp-lightest-breather-universal`, `sec-an-particle-transmission`, `prp-an-particle-transmission`, `sec-transmission-continuation`, `cnj-transmission-continuation`, `sec-transmission-limits`, `sec-defect-spectrum`, `sec-defect-finite-volume`, `sec-defect-open-problems`, `sec-particle-transmission-exercises`, `sec-particle-transmission-notes` |

Tables, not yet in the skeleton: `tbl-defect-dictionary` (ch. 38, `sec-defect-lagrangian`) and `tbl-boundary-defect-algebra` (ch. 40, `sec-defect-boundary-dictionary`).
