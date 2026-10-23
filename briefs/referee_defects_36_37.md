# Brief: referee chapters 36–37 (scattering with a defect; sine-Gordon and sinh-Gordon with a defect)

You are an adversarial mathematical-physics referee. Your working directory is a copy of the repository of a Quarto book on affine Toda field theory. Do not edit any existing file. You may create scratch files only under `referee_scratch/`. Do not spawn subagents. Your final answer, printed in full, is the report: nothing else is collected.

## Context

Part VIII, "Affine Toda theory with defects" (`part8-defects/`, chapters 36–43), was drafted on 2026-10-11 by agents working from the defect literature. The scope is textbook plus verification: results are re-derived in the book's conventions, key formulas are checked by scripts in `code/`, and corrections to the sources are recorded in each chapter's Notes. This brief covers:

- **Chapter 36** (`part8-defects/36-defect-scattering.qmd`, `sec-defect-scattering`): the scattering axioms with a defect.
  - The defect Zamolodchikov–Faddeev algebra `eq-defect-zf`.
  - The reflection–transmission no-go theorem `prp-rt-nogo`.
  - The triangle relations in S/𝕋, braid and RTT forms (`eq-triangle-st`, `eq-triangle-braid`, `eq-triangle-rtt`, `prp-triangle-forms`).
  - Defect unitarity and crossing (`eq-defect-unitarity`, `eq-defect-crossing`).
  - The doubling trick (`prp-doubling`, `eq-doubled-k`).
  - Defect bound states and the bootstrap (`eq-defect-bootstrap`).
  - Defects on the circle and twist defects (`prp-twist-defect`).
  - Topological lines and the modified crossing rule of Copetti–Córdova–Komatsu (`sec-topological-lines`).
  - Finite size in brief.
  - Defect Lee–Yang as a worked example (`eq-ly-transmission`, `prp-ly-transmission`).
- **Chapter 37** (`part8-defects/37-defect-sine-gordon.qmd`, `sec-defect-sine-gordon`): the rank-one template.
  - Bäcklund transformations, introduced here for the first time in the book.
  - The jump defect as a frozen Bäcklund transformation.
  - Solitons meeting a defect (delay, conversion, absorption).
  - Linear transmission.
  - Type-II defects.
  - The quantum transmission matrix on a ℤ-graded defect space, and the test of the sign of q.
  - Resonances.
  - Breather and sinh-Gordon transmission.
  - Fusion onto boundaries.

**Read first.**
- `briefs/part8_defects_conventions.md` holds the conventions, verified by script, and the source dictionary.
- `briefs/part8_defects_status.md` holds the drafting agents' verification records.
- `briefs/part8_defects_sources.md` holds digests of the sources, with their equation numbers.
- To find any object: `grep -rn '{#<id>}' part*/ index.qmd`.

**Book conventions.**
- Longest roots have length² 2, $\beta_{\rm SG}^2=2\beta^2$ and $m_0=2m$.
- The Part VI coproduct is $\Delta e=e\otimes1+k\otimes e$, $\Delta f=f\otimes k^{-1}+1\otimes f$, with the physical $q=-e^{-i\pi\lambda}=e^{-4\pi^2i/\beta^2}$ (`eq-q-physical`).
- $\sigma=e^{-\theta_d}$, where $\theta_d$ is the defect rapidity.
- $\mathbb T:V(\theta)\otimes\mathcal W\to\mathcal W\otimes V(\theta)$ is left-to-right transmission, with tensor order equal to spatial order.
- The two type-I defects are I$_-$ (Bowcock–Corrigan–Zambon's original, intertwining $U_q(\mathfrak b^-)$) and I$_+$ (its time reverse, intertwining $U_q(\mathfrak b^+)$). Type II is Corrigan–Zambon's defect with an auxiliary field.

**Sources.** The text of every source is in `papers/` (pdftotext), named by arXiv id with `/` replaced by `_`:
- hep-th_0305022 (`bcz2004`), hep-th_0401020 (`bcz2004b`), hep-th_0506169 (`bcz2005`), hep-th_0407199 (`cz2004`);
- 0908.3126 (`cz2009b`), 1006.0939 (`cz2010`), 1012.4186 (`cz2011`), 1202.6016 (`cz2012`);
- hep-th_9409076 (Delfino–Mussardo–Simonetti), hep-th_9703085 (Konik–LeClair), hep-th_0201142 (Castro-Alvaredo–Fring–Göhmann);
- hep-th_0404199 (Bajnok–George), 0712.4292 (Bajnok–Simon), 1307.4536 (Bajnok–Holló–Watts);
- 1006.1555 (Weston), 1902.10690 (Penati–Polvara), 2403.04835 and 2408.13132 (Copetti–Córdova–Komatsu);
- 1205.1661 (Avan–Doikou), 1411.5171 (Caudrelier).

**What to expect.** In this repository, every adversarial referee pass so far has found real errors. The drafting has also already corrected several published formulas, and those corrections deserve the hardest look. A claimed error is useful only with evidence: a derivation, a numerical check, or a quotation from a source with its equation number.

## Already verified (do not redo unless you find a reason)

**By script.** Every script prints PASS/FAIL lines; run `python3 code/<script>`.

- `defect_axioms.py` (20 checks, chapter 36):
  - The three triangle forms are equivalent, with positive and negative controls.
  - The no-go reductions, computed by a symbolic ZF engine, give $S^2=1$.
  - The DMS Ising amplitudes satisfy unitarity, crossing and the RT relations.
  - The Lee–Yang statements hold: unitarity, crossing, the bootstrap, fusion onto the boundary giving `eq-ly-reflection`, energies, the bound state, and the pure-phase lines (at sample points).
  - The sinh-Gordon CDD statements of `exr-defect-cdd`.
  - `prp-doubling` items 1–3 (item 1 for the sine-Gordon Ř with 𝒲 = spin 1; items 2–3 for Lee–Yang).
  - `prp-twist-defect`, including the point about the seam.
  - The CCK crossing rules and the factor relating them.
- `defect_classical_sg.py`, `defect_classical_an.py`, `defect_sg_transmission.py` (conventions, 37 checks): every item marked [V:…] in `briefs/part8_defects_conventions.md`.
- `defect_classical_sg.py` (27 checks) and `defect_sg_transmission.py` (24 checks) for chapter 37: see the chapter 37 section of `briefs/part8_defects_status.md`. Key items:
  - the Bäcklund/Bianchi identities, and the identification of the defect conditions with the frozen transformation;
  - the delays and absorption;
  - Θ iff σ real;
  - the scalar factor Φ, with its unitarity, crossing and poles;
  - crossing selects q_phys (fails at −q_phys with C = σ_x);
  - the resonance poles;
  - the breather factors against BCZ (5.10)–(5.13);
  - the Bajnok–Simon fusion identity.

**By hand (drafting agents; worth a second look):**
- Chapter 37: the proofs of compatibility and permutability; dE/dt and dP/dt; the crossing condition for s ≠ 0; the KL dictionary ν ↔ −λθ_d; the renormalization factor 1 − B/2 extracted from Bajnok–Simon (17).
- Chapter 36:
  - the proofs of `prp-rt-nogo`, `prp-twist-defect` and `prp-doubling` item 3;
  - Hermitian analyticity of the Lee–Yang factors ⇔ Re b = ±3;
  - the kink-rescaling interpretation of the CCK factor;
  - the Bajnok–Holló–Watts relation μ̂ + μ̄̂ = h(b) of ch. 28;
  - the signs and arguments of the doubled K.

## To be checked

Report on each item: correct / wrong (with the correction) / unclear (with what would settle it). Items are roughly in order of risk.

### Chapter 36

1. **`prp-rt-nogo` [Theorem].** Check the statement and the proof.
   - Does the meromorphy argument hold?
   - Is it true that parity is not needed?
   - Is the condition $S_{ab}=S_{ba}$ for the both-reflected equation right?
   - Is the paraphrase of Castro-Alvaredo–Fring–Göhmann's assumption (their (22)) faithful?
   - "No interacting example is known" of the remaining loophole: search the literature for a counterexample.
2. **Index conventions.** Check `eq-defect-zf` and its index rule, and the consistency of `eq-triangle-st`, `eq-triangle-braid` and `eq-triangle-rtt` with it. Does 𝕋 depend on θ−θ_d and 𝕋̄ on θ+θ_d?
   - Compare `eq-defect-crossing` and `eq-defect-unitarity` with Bajnok–George (15) after translating conventions. Is the charge conjugation placed correctly?
3. **`prp-doubling`.**
   - K = (0, 𝕋̄P; P𝕋, 0) on (V⊕V′)⊗𝒲.
   - The claim that the reflected copy scatters with PŠP.
   - Item 3 (doubled crossing-unitarity ⇔ defect crossing), which is proved by hand only.
   - The statement that fusion onto a boundary gives K′ = 𝕋̄K𝕋.
4. **`prp-defect-elastic`, second part.** "If a modified momentum is conserved and the defect states have discrete momenta, no particle is reflected." Is the argument sound? What happens for a defect whose states have a continuous momentum?
5. **`prp-twist-defect` [Theorem]** and the claim in the Notes that the Bloch twist of `prp-bloch-sectors` (`part8-foundations/40-scattering.qmd`) must act at a common seam.
6. **Defect Lee–Yang.**
   - The translation of Bajnok–Simon's $[b+1][b-1]$ into the book's $f_x$ blocks (`eq-fx-block`).
   - The defect energy $m\sin(b\pi/6)$.
   - The bound state and its energy.
   - The fundamental range $b\in[-3,2]$.
   - The Bajnok–Holló–Watts integrability constant 0.826608 and the statement that the defect "scatters like a particle of rapidity $i\pi(3-b)/6$".
7. **Modified crossing (`sec-topological-lines`).**
   - The statement of the CCK rule.
   - The claim that the RSOS S-matrices of ch. 27 (`sec-rsos`) satisfy ordinary crossing, and differ from CCK's form by $(d_ad_c/d_bd_d)^{i\theta/2\pi}$.
   - The kink-rescaling interpretation.
8. **Corrections to sources made in the chapter 36 Notes:**
   - DMS (2.12) needs no parity, and their (2.10) is empty only if $S_{ab}=S_{ba}$;
   - the DMS χ branch;
   - CAFG (32) omits a factor $S_{ba}(\theta-i\eta)$.

   Confirm or refute each.
9. **Finite size** (`sec-defect-finite-size`): the defect TBA driving term, the Lüscher term, and the attributions.

### Chapter 37

10. **The sign of q.** The conventions found that the triangle relations hold at $-q_{\rm phys}$ with $\check R(-q_{\rm phys})$, so the chapter must select $q_{\rm phys}$ by unitarity, crossing and the scalar factor. Does the argument in the chapter actually do this, and is the script check meaningful?
11. **Bäcklund transformations.**
    - The light-cone form `eq-sg-backlund` in book conventions.
    - Vacuum → soliton (`prp-sg-backlund`).
    - Bianchi permutability → `eq-interaction-coefficient` (`prp-bianchi`).
12. **The defect Lagrangian and its classification** (`prp-sg-defect-classification`, from BCZ 2003's spin-3 argument). Check that the claimed misprint in BCZ 2005 (2.2) (the sign of δ(x)) really is one.
13. **Solitons through the defect.**
    - The delay `eq-sg-delay`: $z_{\rm I_\mp}=\mp\coth((\theta-\theta_d)/2)$.
    - Absorption at θ = θ_d, and the stored energy and momentum $M\cosh\theta_d$, $M\sinh\theta_d$.
    - `prp-half-soliton` (BCZ's remark that a defect behaves like half a soliton): what exactly is claimed, and is it right?
14. **The transmission matrix and its scalar factor** (`eq-sg-transmission-matrix`, `prp-sg-transmission`).
    - The dictionary to BCZ 2005 (3.51) and Konik–LeClair, including $q_{\rm BCZ}=-1/q_{\rm phys}$.
    - Unitarity and crossing.
    - The resonance poles (`sec-sg-defect-spectrum`).
15. **Breathers and sinh-Gordon** (`prp-sg-breather-transmission`).
    - That the lightest-breather factor is coupling-independent and equals the classical $\mathbb T^{\rm cl}$.
    - The higher breathers.
    - The continuation to sinh-Gordon.
    - Bajnok–Simon's fusion onto a Dirichlet boundary compared with `eq-shg-reflection` (ch. 29).
16. **`cnj-defect-parameter` and the renormalization of σ.** Is the tension between Bajnok–Simon (a one-loop renormalization, no anomaly) and Penati–Polvara (higher-spin anomalies that cannot be removed) described correctly? Is the [Conjecture] label right?
17. **Type II** (`sec-type-two`): the delay of the fused pair, and the attributions to cz2009b and cz2010. The drafting agent noted that cz2010's resonance (3.17) at η_p + iπ/2γ has the opposite sign to BCZ 2005 (4.1)'s η − iπ/2γ, and did not reconcile them: settle this.
17a. **The excited-defect correction to ch. 36** (`sec-defect-bound-states`, last paragraph): the excited defect space is the image of the residue of 𝕋 in 𝒲⊗V(iu), with transmission Ř₂₃𝕋₁₂, and not a subspace of V(iu)⊗𝒲. Confirm this, and confirm that the odd sector of `sec-sg-defect-spectrum` (the s = 1 intertwiner at θ_d − iπ) follows from it.
17b. **`prp-half-soliton`**: the squared identities z_odd² = A₁₁(θ−θ_d), and 𝕋_{B₁}(θ; ϑ+iπ+iξ/2)·𝕋_{B₁}(θ; ϑ−iπ−iξ/2) = S_{sB₁}(θ−ϑ), together with the corrected reading of BCZ's two-defect remark.

### Both chapters

18. **Status labels.** Is anything overclaimed, in particular [Theorem] and [Established] labels that rest on citations alone?
19. **Exercises.** For every `exr-` in the two chapters: is it well posed and solvable from the text, and is any hint wrong?
20. **Bibliography.** For every key cited in the two chapters, check `references.bib` (authors, title, journal, DOI) and whether the citation supports its sentence.
21. **Prior art.** Is anything presented as new already in the literature?

## Deliverable

One entry per item, each with:
- the verdict;
- the justification;
- for errors, corrected text ready to paste, with file and line.

Include any scripts you wrote, inline, with their output. Do not mark anything verified on the strength of a citation alone.
