# Referee report: textbook material of the affine Toda book

Brief: `briefs/referee_textbook.md`. Date: 2026-10-10. No repository file was edited.

**How it was done.**
- Eight independent referee agents covered disjoint item sets: bibliography; q conventions and crossing; tensor-product graph and ch. 7–8; Part I; ch. 5–6 and appendix A; Part III; Parts IV–V; Parts VI–VII framing.
- Every verdict rests on a derivation or a computation. Literature was used only to confirm, never as the sole basis.
- The orchestrator re-ran the baseline: all nine `textbook_code/` scripts pass, with no FAIL.
- The orchestrator independently spot-checked the findings that change the book most:
  - the crossing-sign resolution, by re-running the agents' scripts **and** the repository's own `sign_n1.py`/`sign_n2.py` with q → −q;
  - the floating-H bubble computation, re-run;
  - the Delius–Gould–Zhang prior art, against the PDFs and with a symbolic identity re-run;
  - by hand: the κ sign at 16:38, the conserved-charge sign at 02:50, the central element at 07:109, the real a₃ soliton contradicting 10:10, the CDD staircase counterexample, and the duality conflation at 06:111.
- Scripts are in `referee_textbook_scripts/<agent>/`; see the index at the end.


## The most important findings

1. **Open problem 4 ("Crossing convention") is resolved. The crossing sign (−1)ⁿ is the sign of q (item 2).**
   - Part VI uses the coproduct opposite to ch. 7, so Ř_VI(x;q) = Ř₇(x;q⁻¹).
   - The physical (Bernard–LeClair) value is q = −e^{−iπω} = e^{−4π²i/β²}, so Part VI's q = e^{−iπω} has the wrong sign.
   - With the right sign, crossing holds with +1 for c₁…c₄ and a₅⁽²⁾, with unchanged scalar factors.
   - q → −q is not a gauge transformation. It is a soliton–antisoliton exchange-sign (Klein) twist.
   - **Follow-up needed:** rerun the ch. 32 multi-soliton sector computations with q → −q, since this bears on `briefs/exchange_statistics.md`.
2. **Prior art for the Part VI R-matrices.** Delius–Gould–Zhang 1994/1996 already derive, for all n:
   - the U_q(d_{n+1}⁽²⁾) spinor R-matrix (ch. 23);
   - the U_q(b_n⁽¹⁾) spinor R-matrix (ch. 24);
   - the 26 R-matrix of U_q(f₄⁽¹⁾) (ch. 25).

   These are proved identical symbolically. The book calls them "conjectured for all n" or "constructed here". Fourteen sentences need changing (item 19). The S-matrices built on them remain new.
3. **Ch. 8:73 is wrong.** The c_n spinor tensor square *is* multiplicity-free, and ch. 23 says so. The tensor-product-graph rule as stated is for untwisted algebras only; the twisted rule needs ⟨a⟩₊ edges (item 4). The rule itself is verified for the a₂, b₂, b₃, c₂ and c₃ vector representations.
4. **`eq-floating-h` is wrong as a general formula (item 15).**
   - b_n⁽¹⁾ has H = 2n − B/2, decreasing to h^∨.
   - In the book's normalization, c_n, g₂ and f₄ have H = 2n + B/2, 6 + B and 12 + 3B/2.
   - The printed formula holds for c_n, g₂ and f₄ only in the short-root normalization.
5. **Twisted algebras: the conserved spins repeat with period kh, not h** (06:142, 09:32, 09:84, 02:59). Found independently by two agents.
6. **Smaller formula errors:**
   - **Normalization and convention errors:** the central element must be ∏k_i^{n_i}, not ∏k_i^{n_i^∨} (07:109); the S² shift is q^{2h^∨d}, not q^{h^∨d} (07:182); κ = ik, not −ik (16:38); the conserved charge is ∫(T−Θ) (02:50).
   - **Wrong constants and factors:** the Bernard–LeClair current has topological charge ±2, not ±1 (20:16); q = e^{−iπh_J}, not ±e^{iπh_J} (20:55); a spurious factor i at 12:55; the strip of eq-sg-s0 is min(ξ,π) (03:138).
   - **Missing pieces:** the OTU τ-function lacks its x,t dependence (06:159); the null-root spin sum needs n_j^∨ (20:64, 20:70).
7. **Statements that are false as written:**
   - `prp-cdd` misses x = −½+iy (the staircase S-matrix) and conjugate pairs.
   - The Mostafazadeh claim needs equal multiplicities.
   - `exr-krein-collision`.
   - 10:10, real a₃ solitons exist.
   - 11:43, conserved charges are not all real.
   - 06:111 and 14:53 conflate the two dualities.
   - The degenerate-mass Dorey and area rules hold only in the charge basis (09:120).
   - The HHM mechanism is misdescribed in six places.
   - ch. 27's "other algebras" row cites a₂⁽²⁾ papers.
   - Two Part VII learning goals misstate the 2π threshold.
8. **Bibliography.**
   - The 39 new entries are all correct.
   - `smirnov1991` is cited for the φ₁,₃ restriction, but that is Smirnov 1990.
   - The Freeman, Gandenberger, Faddeev and Hollowood a/b attributions are misplaced.
   - `dgz1994` has the wrong title and `dgz1996` lacks title and journal.
   - `takacs1997b` and `takacs1997d` duplicate each other.
   - The Monte Carlo papers are identified, with BibTeX.
9. **Hollowood: the book's ¼cot(π/h) is right.** Hollowood's printed ½ comes from doubled intermediate sums, so the attribution needs a note (item 17).

**Revision to the brief's "already verified" list.** One by-hand item does not survive: `eq-floating-h` "reproduces H = 2n+B and 12+3B" is true only if B is the short-root block parameter, not the B of eq-b-of-beta in book normalization. All script-verified items still pass.

**Exercises (item 12).** Of 122:
- 6 are false or misleading as stated;
- 3 are imprecise;
- about 10 need wording or hint fixes;
- the rest are true and solvable.


## Verdict per brief item

| # | Item | Verdict |
|---|---|---|
| 1 | Bibliography | 39 new entries correct; 3 older entries wrong; about 10 misattributed citations; several missing |
| 2 | q conventions | Dictionary exact; Part VI's q has the wrong sign; **crossing sign explained** |
| 3 | Crossing point in general | Correct with omitted conditions; two sentences at 07:182 wrong |
| 4 | TPG beyond sl₂ | Untwisted rule correct (a₂, b₂, b₃, c₂, c₃ vectors); twisted spinor claim wrong; 08:82 example wrong |
| 5 | Sine-Gordon quoted results | Breather masses and S₁₁ correct (bootstrap derived); residue-sign rule too strong; 03:181 incomplete; strip wrong for λ < 1 |
| 6 | RSOS and Lee–Yang | Correct; citation wrong; one exercise misleading |
| 7 | Ch. 2 statements | `prp-cdd` wrong; `eq-bootstrap` correct; angle-sum indices; charge sign; exercises true |
| 8 | Ch. 4 statements | Mostafazadeh claim wrong; `exr-krein-collision` false; 04:112 wrong |
| 9 | Ch. 5 statements | Correct; attribution and wording fixes |
| 10 | Ch. 6 statements | All foldings correct; E1–E7 (duality conflation, twisted spins, h^∨ normalization, OTU formula, centre statement) |
| 11 | Ch. 7 statements | Central element wrong; KR terminology; quasitriangularity correct; `exr-hermitian-spin1` correct |
| 12 | Exercises | See table |
| 13 | Part III | Hirota (folded), soliton masses, McGhee (proved) and e₈ correct; basis prescription, twisted spins, 10:10, 10:75, 10:86 and 11:43 wrong |
| 14 | Part III citations | Correct |
| 15 | Part IV | `eq-floating-h` wrong; duality and status overclaims; 12:55, 12:45, 12:81, 13:53; exr-cn-snn; MC papers identified |
| 16 | Part IV citations | `smirnov1992` correct; `fring1992` overstated |
| 17 | Part V | ¼ correct (Hollowood's ½ is an error); κ sign; exr-hollowood; other statements correct |
| 18 | Part V citation | Correct |
| 19 | Part VI | Ch. 20 errors; a_n example correct; DGZ prior art; HHM misdescribed; ch. 27 table citations; q-sign follow-up |
| 20 | Appendices and Part VII | Tables correct; twisted masses can be made explicit; two learning goals wrong; two exercises imprecise |

## Items 1, 14, 16, 18: bibliography and citations

Method: each entry checked against its INSPIRE record (ID given), Crossref, or the arXiv/publisher page. Each citing sentence was checked against the cited paper's abstract or text. Query logs were not kept.

### 1a. The 39 appended entries (`parke1980` … `cahill1976`): **all correct**

These were verified against INSPIRE:
- parke1980 (152771), shankar1978 (6061), coleman1975 (1873), mandelstam1975 (98336), dhn1975 (2037), faddeev1978 (5572), coleman_thun1978 (129104), arinshtein1979 (147655; the printed title reads "Todd chain" (sic)), leclair1989 (278251)
- heiss2012 (1205743), gupta1950 (3551), bleuler1950 (9186), leeyang1952 (8637), cardy1985 (225153), yurov1990 (285394), goddard1986 (18583), jimbo1985 (225728), jimbo1986 (213089), frenkel1992 (313731), mackay1991 (314348), kulish1981 (171634)
- mikhailov1981 (154439), fring1991 (318574), fring1992 (321317), smirnov1992 (344254), cahill1976 (109824)

These were verified via Crossref or the publisher: mussardo2010, rajaraman1982, bognar1974, kato1966, humphreys1972, fuchs1997, kostant1959, kac1990, drinfeld1987, kassel1995, chari1994, gomez1996.

Two entries are incomplete rather than wrong:
- `dorey1998` is published. Use `@incollection{dorey1998, author={Dorey, P.}, title={Exact {S}-matrices}, booktitle={Conformal Field Theories and Integrable Models (Eötvös Graduate Course, Budapest 1996)}, editor={Horváth, Z. and Palla, L.}, series={Lect. Notes Phys.}, volume={498}, pages={85--125}, publisher={Springer}, year={1997}, url={https://arxiv.org/abs/hep-th/9810026}}`.
- `smirnov1992`: optionally add `series={Adv. Ser. Math. Phys.}, volume={14}`.

### 1b. Errors in older entries (outside the 39)

- **`dgz1994`, wrong title.** NPB 432 (1994) 377 is "On the construction of trigonometric solutions of the Yang–Baxter equation", hep-th/9405030 (INSPIRE 373180).
- **`dgz1996`, no title or journal.** It is "Twisted quantum affine algebras and solutions to the Yang–Baxter equation", Int. J. Mod. Phys. A 11 (1996) 3415–3438, q-alg/9508012.
- **`ahn2000`** should be `@article`: Phys. Lett. B 481 (2000) 114 (INSPIRE 524375).
- Unused entries: `dorey1992`, `dorey1993`, `ge2012`, `konno1966`.

### 1c. Citations in chapters 1–8: wrong or misattributed

| Location | Problem | Fix |
|---|---|---|
| 03-sine-gordon.qmd:193, 196, 253; 27-rsos.qmd:102 | `smirnov1991` (IJMPA 6 (1991) 1407) is the φ₁,₂/a₂⁽²⁾ paper, not the φ₁,₃ restriction of sine-Gordon. Confirmed independently by the Part I referee. | Add `@article{smirnov1990, author={Smirnov, F.A.}, title={Reductions of the sine-{Gordon} model as a perturbation of minimal models of conformal field theory}, journal={Nucl. Phys. B}, volume={337}, pages={156}, year={1990}}` and cite it there. Optionally add Bernard–LeClair NPB 340 (1990) 721. Keep `smirnov1991` at 04-non-hermitian.qmd:140 (Lee–Yang from a₂⁽²⁾), which is correct. |
| 02-integrable-qft.qmd:37 | "Zamolodchikov–Faddeev algebra [@zz1979; @faddeev1978]": Faddeev–Korepin 1978 is about semiclassical soliton quantization. | `@article{faddeev1980, author={Faddeev, L.D.}, title={Quantum completely integrable models in field theory}, journal={Sov. Sci. Rev. C}, volume={1}, pages={107}, year={1980}}` |
| 03-sine-gordon.qmd:253 | "Faddeev and Korepin … together with the quantum inverse scattering method": their review predates the QISM. | "…together with its relation to the classical inverse scattering method" |
| 05-root-systems.qmd:164; 09-lagrangian.qmd:98 | "charge eigenvalues are Cartan eigenvectors" attributed to Freeman. Freeman 1991 proved only the mass statement. | cite [@klassen1990; @dorey1991] (Corrigan's review hep-th/9412213 §2) |
| 01-introduction.qmd:11 | strong–weak duality cites only `dgz1992` | add `@braden1990` (first stated there) |
| 01-introduction.qmd:12 | `hollowood1993a` (mass corrections) cited for quantum-group-symmetric solitons | `hollowood1993b` |
| 08-r-matrices.qmd:52, 142 | trigonometric tensor-product graph credited to MacKay (rational) and DGZ only | add `@article{zgb1991, author={Zhang, R.B. and Gould, M.D. and Bracken, A.J.}, title={From representations of the braid group to solutions of the {Yang}--{Baxter} equation}, journal={Nucl. Phys. B}, volume={354}, pages={625}, year={1991}}`. Replacement wording is in item 4. |
| 06-affine-algebras.qmd:79, 166 | "we follow [@mackay_watts1995]" for twisted labels; MacKay–Watts themselves follow BCDS | **unclear**: cite `braden1990` or Kac's tables directly |

Missing citations:
- Steinberg in `sec-weyl-coxeter` (`@article{steinberg1959, author={Steinberg, R.}, title={Finite reflection groups}, journal={Trans. Amer. Math. Soc.}, volume={91}, pages={493}, year={1959}}`, together with `kostant1959`).
- Coleman's bound at 03:179 (`coleman1975`).
- Collapse thresholds at 03:181 (`frohlich1976`, `benfatto1982`).
- 04:65–66: Bender–Boettcher only conjectured reality for ε ≥ 0; the proof is Dorey–Dunning–Tateo, J. Phys. A 34 (2001) L391.
- 04:68: the antilinear-symmetry half of the Mostafazadeh claim is Mostafazadeh, J. Math. Phys. 43 (2002) 3944 (math-ph/0203005); see item 8.
- 09:82: the Drinfeld–Sokolov procedure (`@article{drinfeld_sokolov1985, author={Drinfeld, V.G. and Sokolov, V.V.}, title={{Lie} algebras and equations of {Korteweg}--de {Vries} type}, journal={J. Sov. Math.}, volume={30}, pages={1975}, year={1985}}`).

All other citing sentences in ch. 1–8 are supported; the list of checked lines is in the bibliography agent's log.

### 14. `mikhailov1981`, `fring1991`: **correct**

- `mikhailov1981` (09:182): the abstract covers the zero-curvature representation, the conservation laws and the mass spectrum.
- `fring1991` (09:186): the area rule is in the abstract, and Corrigan's review credits Fring–Liao–Olive with deriving the fusing rule from the Lagrangian.
- 10:63 is **unclear**: `hollowood1992` is cited for d_n/e_n single solitons, but the paper treats a_n only, as far as can be determined. Candidates are mackay1993 and aratyn1993.

### 16. `fring1992`, `smirnov1992`

- `smirnov1992` (15:37, 74): **correct**; it is the standard reference for the form-factor axioms.
- `fring1992` (13:132), "derived the S-matrices from the Coxeter geometry": **overstated**. The abstract says crossing and bootstrap of the conjectured S-matrices are proved uniformly. Replace with: "Fring and Olive [@fring1992] proved uniformly, from Coxeter-element identities, that Dorey's formula satisfies crossing and the bootstrap."

### 18. `cahill1976`: **correct**

- The abstract describes general one-loop mass formulas for static solitons.
- The printed form reproduces −m/π (sine-Gordon) and m(1/(4√3) − 3/(2π)) (φ⁴).
- The original paper was not accessible.
- 16:70: the φ⁴ kink correction was first published in DHN, PRD 10 (1974) 4130. Consider adding it, together with PRD 10 (1974) 4114, for "the method is due to DHN".

### Other attributions (bibliographic side; the content checks are in later items)

- **Hollowood:**
  - 16:128, 22-known-smatrices.qmd:34 and 18-jordan-chains.qmd:29 cite `hollowood1993b` (the S-matrix paper, IJMPA 8 (1993) 947, hep-th/9203076) for a mass or fluctuation computation. Use `hollowood1993a` (PLB 300 (1993) 73, hep-th/9209024).
  - 21-construction.qmd:128: `hollowood1993a` does not support "S = ΦŘ"; use `hollowood1993b`.
- **21-construction.qmd:104 (Gandenberger): wrong twice.** Gandenberger 1995 (hep-th/9501136) treats a₂⁽¹⁾ only, and Smirnov had the a₂⁽²⁾ breather identification earlier. Replace with: "Gandenberger [@gandenberger1995] showed for $a_2^{(1)}$ that the lowest breathers $B_1,\bar B_1$ scatter as the two particles of $a_2^{(1)}$, with the S-matrix (-@eq-an-smatrix) continued to imaginary coupling; Smirnov [@smirnov1991] had found the analogous identification for $a_2^{(2)}$."
- **22-known-smatrices.qmd:34:** the raw "[Watts 1994]" is `@article{watts1994, author={Watts, G.M.T.}, title={Quantum mass corrections for {$C_2^{(1)}$} affine {Toda} theory solitons}, journal={Phys. Lett. B}, volume={338}, pages={40}, year={1994}, url={https://arxiv.org/abs/hep-th/9404065}}`.
- **15:48, 15:25, 30-continuum-thresholds.qmd:37 (mass–coupling relation, TBA reflection amplitudes):**
  - Fateev 1994 treats simply-laced perturbed cosets; Ahn–Fateev–Kim–Rim–Yang derive the Toda m–μ relation from it.
  - Add `@article{afkry2000, author={Ahn, C. and Fateev, V.A. and Kim, C. and Rim, C. and Yang, B.}, title={Reflection amplitudes of {ADE} {Toda} theories and thermodynamic {Bethe} ansatz}, journal={Nucl. Phys. B}, volume={565}, pages={611}, year={2000}, url={https://arxiv.org/abs/hep-th/9907072}}` for the simply-laced statements. `ahn2000` is the non-simply-laced paper.

### 15 (last bullet). Monte Carlo citations at 13-simply-laced.qmd:97

- hep-th/9206112 is `@article{watts_weston1992, author={Watts, G.M.T. and Weston, R.A.}, title={{$G_2^{(1)}$} affine {Toda} field theory: a numerical test of exact {S}-matrix results}, journal={Phys. Lett. B}, volume={289}, pages={61}, year={1992}, url={https://arxiv.org/abs/hep-th/9206112}}`.
- hep-th/9508007 is `@article{beccaria1996, author={Beccaria, M.}, title={Monte {Carlo} study of exact {S}-matrix duality in non-simply laced affine {Toda} theories}, journal={Phys. Rev. D}, volume={53}, pages={3266}, year={1996}, url={https://arxiv.org/abs/hep-th/9508007}}`.

What they confirm is only the coupling dependence of a mass ratio, for g₂⁽¹⁾ (and d₄⁽³⁾ in Beccaria). They test no simply-laced theory and no S-matrix element directly. Suggested text: "…and the coupling dependence of the mass ratio predicted for the dual pair $g_2^{(1)}/d_4^{(3)}$ has been confirmed by Monte Carlo simulation [@watts_weston1992; @beccaria1996]."

### Prior art affecting Part VI (found under items 4 and 19, verified by the orchestrator against the papers)

The R-matrices that Part VI presents as constructed here, or as conjectured for all n, are already derived in two earlier papers:

- **Delius–Gould–Zhang 1996**, §4.3, eqs. (4.27)–(4.31): the U_q(D_{l+1}^{(2)}) R-matrix on V(λ_l)⊗V(aλ_l) for all l ≥ 2, from the twisted tensor-product graph. At a = 1 this is the c_n⁽¹⁾ spinor R-matrix of ch. 23.
- **DGZ 1994**:
  - §4.2, eqs. (4.21)–(4.25): the U_q(B_n^{(1)}) spinor⊗spinor R-matrix for all n.
  - §4.6, eq. (4.64): the F₄ V(λ₄)⊗V(λ₄) R-matrix. Normalized to 1 on the 324 it reproduces the 26⊗26 table of 25-e6-solitons.qmd:11–24 exactly:
    - DGZ coefficients ⟨−18⟩, ⟨−18⟩⟨−8⟩, …, rescaled using ⟨a⟩⟨−a⟩ = 1;
    - result 273: ⟨2⟩, 52: ⟨8⟩, 26: ⟨2⟩⟨12⟩, 1: ⟨8⟩⟨18⟩.

Sentences to revise are listed under item 19.

## Item 2: $q$ conventions across parts, and the crossing sign

**Verdict: the dictionary is now exact, and the unexplained crossing sign $(-1)^n$ (open problem 4, "Crossing convention") is resolved. It is the sign of $q$.**

Scripts: `referee_textbook_scripts/qg_cross/qdict.py` (symbolic), `qtype.py`, `qsign_cn.py`, `qsign_bn.py`, `qsign_n4.py`, `qsign_L.py`, `sg_scalar.py`. The orchestrator re-ran `qdict.py` and `qsign_cn.py` for n = 1, 2, 3 and got the same output. It also ran the repository's own `soliton_rmatrix_code/smatrix/sign_n1.py` and `sign_n2.py` with `q → ±e^{-iπω}`:

| run | crossing point x(iπ) | κ |
|---|---|---|
| n = 1, q = +e^{−iπω} | +q⁻² | −1 |
| n = 1, q = −e^{−iπω} | +q⁻² | +1 |
| n = 2, q = ±e^{−iπω} | −q⁻⁴ | +1 |

**The dictionary** (all entries checked symbolically):

| | Ch. 7 | Part VI (research code, ch. 20) |
|---|---|---|
| coproduct | Δ(E)=E⊗K+1⊗E, Δ(F)=F⊗1+K⁻¹⊗F | Δ(e)=e⊗1+k⊗e, Δ(f)=f⊗k⁻¹+1⊗f, i.e. **P Δ₇ P, the opposite coproduct** |
| R-matrix | Ř₇(z;q) | Ř_VI(z;q) = P Ř₇(1/z;q) P = **Ř₇(z;q⁻¹)**, same z, gradation and normalization |
| sine-Gordon spectral variable | z = e^{2λθ} | x = e^{2Tθ}, T = λ, so x = z |
| crossing point | x_c = q² | x_c = q⁻² |
| physical sine-Gordon q | −e^{iπλ} = e^{iπh_J} | **−e^{−iπλ} = e^{−iπh_J} = e^{−8π²i/β_SG²} = e^{−4π²i/β²}**, Bernard–LeClair (3.18) with their coproduct (3.37) |
| Part VI as printed | | q = e^{−iπω}, which is **minus** the physical value |

The by-hand finding in the brief, $q_{\rm VI}=-q_7^{-1}$ and $q_7^2=q_{\rm VI}^{-2}$, is correct.

**q → −q is not a gauge transformation.**
- It leaves a, c± and x_c = q² unchanged.
- It flips b, the gauge invariant Δ = (a²+b²−c₊c₋)/2ab = (q+q⁻¹)/2, and c_u.
- It equals conjugation by L = diag(1,1,−1,1), which is not of the form A⊗A: it changes the relative exchange sign of soliton and antisoliton, i.e. it is a Klein factor.
- Equivalently, it is the type-(−1) representation (k → −k) of the same U_q.
- With q = e^{−iπω}, the Part VI n = 1 amplitude is Zamolodchikov's with S_T → −S_T. With q = −e^{−iπω} it is Zamolodchikov's.

**Crossing test (the book's κ = [c(T−t)/c(t)]/c_u(1/x)), same T, a_j and c·f in both columns:**

| case | x_c | κ at q = +e^{−iπω} | κ at q = −e^{−iπω} |
|---|---|---|---|
| c₁ (sine-Gordon) | q⁻² | −1 | +1 |
| c₂ | −q⁻⁴ | +1 | +1 |
| c₃ | q⁻⁶ | −1 | +1 |
| c₄ | −q⁻⁸ | +1 | +1 |
| b₂⁽¹⁾ spinor (a₃⁽²⁾ solitons) | q⁻⁶ | +1 | +1 |
| b₃⁽¹⁾ spinor (a₅⁽²⁾ solitons) | q⁻¹⁰ | −1 | +1 |

**[Conjecture]** κ(q)/κ(−q) = (−1)^{2(λ,ρ)}, with λ the highest weight of the horizontal module (short roots of length² 2). For b_n spinors 2(λ,ρ) = n². It is even for the 26/27 of f₄ and the 7/8 of g₂, so no sign is predicted there, consistent with the +1 reported for f₄⁽¹⁾. Not proved.

**Corrections.**
- **07-quantum-groups.qmd:168.** Replace "Part VI writes $q=e^{-i\pi\omega}$; conventions for $q$ differ between references by $q\leftrightarrow q^{-1}$ and by signs, and we convert when needed." with:
  > Part VI uses the opposite coproduct, $\Delta(e_i)=e_i\otimes1+k_i\otimes e_i$, $\Delta(f_i)=f_i\otimes k_i^{-1}+1\otimes f_i$, the form in which the non-local charges act (@sec-qg-coproduct). Since $P\check R(1/z)P$ intertwines the opposite coproduct, its six-vertex matrix is (-@eq-six-vertex) with $q\to q^{-1}$, at the same $z$, gradation and normalization, and the sine-Gordon value becomes $q=-e^{-i\pi\lambda}=e^{-8\pi^2i/\beta_{\rm SG}^2}$, Bernard and LeClair's [@bernard1991]. The sign of $q$ cannot be absorbed in a change of basis: $q\to-q$ leaves $c_\pm$ and $x_c=q^2$ unchanged but reverses $b$, and with it the basis-independent combination $(a^2+b^2-c_+c_-)/2ab=(q+q^{-1})/2$ and the factor $c_u$ of @sec-qg-crossing. It is conjugation of $\check R$ by $\mathrm{diag}(1,1,-1,1)$, a change of the relative exchange sign of soliton and antisoliton, and it is what the crossing sign of @sec-crossing-sign detects.
- **21-construction.qmd:72 (`sec-crossing-sign`).** Replace the paragraph with:
  > Our numerical crossing test compares $S(i\pi-\theta)$ with $(C\otimes1)S_{21}(\theta)^{t_1}(C^{-1}\otimes1)$. With $q=e^{-i\pi\omega}$ it returns $(-1)^n$ for the $c_n^{(1)}$ spinor and $-1$ for the $a_5^{(2)}$ spinor. This sign is the sign of $q$. The crossing point, $T$, the eigenvalues of $\check R$, the zeros $a_j$ and $c\,f$ depend on $q$ only through $q^2$. But $q\to-q$ conjugates $\check R$ by a diagonal sign matrix on $V\otimes V$ that is not a product of one-soliton basis changes, and multiplies $c_u$ by $(-1)^n$. With $q=-e^{-i\pi\omega}=e^{-4\pi^2i/\beta^2}$, the value of Bernard and LeClair for the coproduct of @sec-qg-coproduct, the test returns $+1$ for $n=1,2,3,4$ and for $a_5^{(2)}$, with the same $c\,f$; for $n=1$ the amplitude is Zamolodchikov's, up to the overall sign of the scalar factor. So $q=e^{-i\pi\omega}$ describes the same theory with the opposite soliton–antisoliton exchange sign (for $n=1$, $S_T\to-S_T$), and no CDD factor is needed.

  Preferred alternative: adopt $q=-e^{-i\pi\omega}=e^{-4\pi^2i/\beta^2}$ throughout Part VI. Pole positions are unaffected, since they are all at even powers of q. The consequential sign changes are:
  - 23-cn-solitons.qmd:32: $q=-e^{-4\pi^2i/\beta^2}$ → $e^{-4\pi^2i/\beta^2}$;
  - 24-a2n-solitons.qmd:33: $-e^{-2\pi^2i/\beta^2}$ → $e^{-2\pi^2i/\beta^2}$;
  - ch. 24 lattice paragraph: $-e^{-i\pi^2/2\gamma}$ → $e^{-i\pi^2/2\gamma}$.
- **23-cn-solitons.qmd:30–34**, "Crossing then fixes $q=-e^{-4\pi^2i/\beta^2}$": **wrong as a claim**, because crossing involves only q². Replace with: "Crossing fixes $T$; the sign of $q$ is fixed by the coproduct of the charges (@sec-qg-coproduct), $q=e^{-4\pi^2i/\beta^2}$."
- **21-construction.qmd:86.** Replace "(…$q=e^{-i\pi\omega}$, which corresponds to $q\to q^{-1}$ relative to @sec-six-vertex)" with "(the coproduct of Part VI is the opposite of that of @sec-hopf, which replaces $q$ by $q^{-1}$; and $q=e^{-i\pi\omega}$ is minus the sine-Gordon value, see @sec-crossing-sign)". The pole at x = q⁻² and the zero at x = q² that follow are correct. For the sl₃ vector, Ř_VI has eigenvalues 1 and ρ₇(x;q⁻¹) (`an_vec.py`).
- **26-f4-solitons.qmd:161, open problem 4.** Mark it resolved as above. What remains open is a proof of the conjectured general sign rule.
- **Side remark.** For n = 1 the Part VI scalar factor is F = −S₀ of eq-sg-s0 (ratio −1.00000 at four rapidities), because Part VI normalizes F(0) = +1 while S₀(0) = −1. This does not affect κ.


## Item 3: crossing point in general

**Verdict: $V(x)^{**}\cong V(x\,q^{2h^\vee})$ is correct, under conditions the text omits. Two sentences at 07:182 are wrong.**

Method (`sq_antipode.py`):
- Chevalley generators for U_q(sl₂^) spin ½, U_q(sl₃^) vector, U_q(c₂⁽¹⁾) vector and U_q(d₃⁽²⁾) spinor, with all relations including q-Serre satisfied to 3·10⁻¹⁶.
- π(S²(a)) and π(S(a))ᵀ built for both coproducts, and scanned for equivalence with π_{xt}.

| rep | V\*\* shift, ch. 7 coproduct | Part VI coproduct | V\* = V(xt) |
|---|---|---|---|
| sl₂^ spin ½ | Q⁴ | Q⁻⁴ | Q² |
| sl₃^ vector | Q⁶ | Q⁻⁶ | V\* is the 3̄, not a shift of V |
| c₂⁽¹⁾ vector | Q⁶ | Q⁻⁶ | Q³ |
| d₃⁽²⁾ spinor | Q⁸ | Q⁻⁸ | −Q⁴ |

Derivation: S² = Ad(q^{2ρ̂}) with ρ̂ = ρ̄ + h^∨Λ₀. Its Λ₀ part rescales x by q^{(α₀,α₀)+2(ρ̄,θ₀)} = q^{2h^∨}.

Conditions:
- q_i = q^{(α_i,α_i)/2} with (α₀,α₀) = 2 (Kac's normalization);
- h^∨ is that of the algebra whose quantum group acts, i.e. ĝ^∨ for Toda solitons;
- with the Part VI coproduct the shift is q^{−2h^∨};
- in the book's own normalization of a twisted algebra the shift is q^{2h^∨/k};
- in the principal gradation the shift is q^{2h^∨/h}, up to a root of unity;
- crossing twice gives x_c² = q^{2h^∨}, so S² fixes x_c only up to sign. Example: d₃⁽²⁾ has x_c = −Q⁴.

**Against Part VI:** every quoted crossing point fits $x_c=(-1)^{kh+h^\vee}q^{-kh}$ with q = e^{−iπω}:

| algebra | x_c |
|---|---|
| c_n⁽¹⁾ | (−1)^{n+1}q^{−2n} |
| a_{2n−1}⁽²⁾ | q^{−2(2n−1)} |
| e₆⁽²⁾ | q⁻¹⁸ |
| f₄⁽¹⁾ | −q⁻¹² |
| g₂⁽¹⁾ | q⁻⁶ |
| d₄⁽³⁾ | q⁻¹² |

In each case x_c² = q_Kac^{−2h^∨(ĝ^∨)}, which is the S² shift. The power of q comes from 4πkh/β² in 2T. The −h^∨(ĝ) term only supplies the sign e^{−iπh^∨}, which S² leaves open.

**Correction, 07-quantum-groups.qmd:182.** "S² is conjugation by an element that includes $q^{h^\vee d}$" is wrong by a factor 2, and "this is why $h^\vee$ appears in $T$" is wrong. Replace from "In $U_q(\hat g)$, $S^2$ is…" through "…parameter $T$ (@sec-amplitude-and-gradation)." with:
> In $U_q(\hat g)$, $S^2(e_i)=q_i^2e_i$ and $S^2(f_i)=q_i^{-2}f_i$, so $S^2$ is conjugation by $q^{2\hat\rho}$, $\hat\rho=\bar\rho+h^\vee\Lambda_0$. Its $\Lambda_0$ part multiplies $e_0$ by $q^{2h^\vee}$, which in the homogeneous gradation is a rescaling of $x$, so $V(x)^{**}\cong V(x\,q^{2h^\vee})$ [@frenkel1992; @chari1994]. Here $q_i=q^{(\alpha_i,\alpha_i)/2}$ with $(\alpha_0,\alpha_0)=2$ (for untwisted algebras the normalization of @sec-lie-algebraic-data; for twisted ones the short roots then have squared length 2), and $h^\vee$ is the dual Coxeter number of the algebra whose quantum group acts — for the solitons of the Toda theory of $\hat g$, of $\hat g^\vee$. With the opposite coproduct of Part VI the shift is $q^{-2h^\vee}$. If $V(x)^*\cong V(x\,x_c)$, crossing twice gives $x_c^2=q^{2h^\vee}$, which fixes $x_c$ up to a sign; for $U_q(\widehat{sl}_2)$, $x_c=q^2$, $x_c^2=q^4=q^{2h^\vee}$. In the gradation of @sec-amplitude-and-gradation, $x(i\pi)=e^{2\pi iT}=(-1)^{h^\vee}e^{4\pi^2ikh/\beta^2}$: the power of $q$ is the $S^2$ shift, and the dual Coxeter number of the Toda algebra in $T$ supplies only the sign that $S^2$ leaves open.

Also, 07:172: "$V(x)^*\cong V(x\,x_c)$" holds only for self-conjugate V. Write "$V(x)^*\cong\bar V(x\,x_c)$, with $\bar V$ the conjugate (antiparticle) representation ($\bar V=V$ for $U_q(\widehat{sl}_2)$)". eq-r-crossing as printed is for self-conjugate V.

## Item 4: tensor-product-graph rule beyond $\widehat{sl}_2$

The full agent report is `referee_textbook_scripts/qg_tpg/AGENT_REPORT.md`; scripts `referee_textbook_scripts/qg_tpg/qaff.py`, `tpg_checks.py`, `spinor_twisted.py`, `ch7_checks.py`.

**(a) U_q(a₂⁽¹⁾) vector: correct.**
- V⊗V = 6(+) ⊕ 3(−), with ρ_{λ₂} = ⟨1⟩.
- Jimbo's equations were solved in the book's coproduct and homogeneous gradation at random q, x, y. The solution space is 1-dimensional, and the eigenvalues match the rule.

**(b) Non-simply-laced vector representations: correct in the book's normalization.** Checked for c₂ (C⁴), b₂ (C⁵), c₃ (C⁶) and b₃ (C⁷), with q-Serre using q_i = q^{α_i²/2}.
- c_n⁽¹⁾:
  - Ř = P_{2λ₁} + ⟨½⟩P_{λ₂} + ⟨(n+1)/2⟩P₀;
  - graph L(λ₂) — L(2λ₁) — L(0), with no λ₂–0 edge (same parity);
  - this matches Jimbo's ξ = k^{2n+2}, k = q^{1/2}.
- b_n⁽¹⁾:
  - Ř = P_{2λ₁} + ⟨1⟩P_{λ₂} + ⟨1⟩⟨(2n−1)/2⟩P₀ (a chain);
  - this matches ξ = k^{2n−1}, k = q.
- Negative controls fail, as they should: the Casimir normalized with long roots of length² 4, and the rule with z → 1/z.
- Suggestion: add the c_n vector as a second worked example, because half-integer arguments ⟨½⟩ appear there.

**(c) The c_n⁽¹⁾ spinor claim at 08-r-matrices.qmd:73: wrong.**
- The subalgebra is U_q(b_n): remove an end node of d_{n+1}⁽²⁾, which is 0 ⇐ 1 — … — (n−1) ⇒ n. Writing U_q(g_(0)) is technically right, but it should say U_q(b_n).
- Spinor⊗spinor = ⊕_{k=0}^n Λ^k C^{2n+1}, which is **multiplicity-free for all n**. Numerically, for n = 2, 3, 4:
  - the commutant has dimension n+1;
  - the Jimbo solution is unique.
- This contradicts the "when it fails" claim, and ch. 23:20 says the opposite.
- The untwisted rule (all signs −) fails here. The twisted rule of DGZ 1996 holds: in book units,
  - ρ_{Λ^{n−j}} = ∏_{i=1}^{j}⟨i/2⟩_{(−1)^i};
  - ⟨a⟩₊ = (1 + z q^{2a})/(z + q^{2a}) on equal-parity edges.
- Real multiplicities occur for the *other* c_n⁽¹⁾ multiplets (Kirillov–Reshetikhin modules reducible under U_q(b_n)), for example Λ¹⊕Λ⁰ and the c₃ 29-plet.

Replacement for 08:73, the "When it fails" paragraph:
> **Twisted algebras.** The rule as stated is for untwisted $\hat g$, in the homogeneous gradation (other gradations are related to it by an $x$-dependent change of basis). For a twisted algebra $g^{(2)}$, $e_0$ transforms under $U_q(g_{(0)})$ not in the adjoint representation but in $g_{(1)}$ (@sec-twisted), and the graph changes accordingly [@dgz1996]: $\mu$ and $\nu$ are joined if $L(\nu)\subset L(\mu)\otimes g_{(1)}$, an edge may join two components of equal parity, and along such an edge $\langle a\rangle$ in (-@eq-tpg-rule) is replaced by $\langle a\rangle_+=(1+z\,q^{2a})/(z+q^{2a})$. The example that matters here is the spinor representation of $U_q(d_{n+1}^{(2)})$, which carries the spinor soliton of $c_n^{(1)}$. Under the horizontal subalgebra $U_q(b_n)$ its tensor square is $\bigoplus_{k=0}^n\Lambda^k\mathbb C^{2n+1}$, without multiplicities, the graph is the chain $\Lambda^n-\Lambda^{n-1}-\dots-\Lambda^0$, and $\rho_{\Lambda^{n-j}}=\prod_{i=1}^{j}\langle i/2\rangle_{(-1)^i}$ with $\langle a\rangle_-=\langle a\rangle$ [@dgz1996].
>
> **When it fails.** The graph defined above is what Delius, Gould and Zhang call the extended graph. The pairs of components that $e_0$ actually links form a connected subgraph of it whenever $V(x)\otimes V(y)$ is irreducible for generic $x/y$ [@dgz1994], so imposing the rule on every edge gives the R-matrix provided the relations are consistent around loops. What the method cannot do without is a multiplicity-free decomposition, and for several affine Toda theories this fails. The soliton multiplets of $c_n^{(1)}$ other than the spinor are Kirillov–Reshetikhin modules of $U_q(d_{n+1}^{(2)})$ that are reducible under $U_q(b_n)$, for example $\Lambda^1\oplus\Lambda^0$ and, for $c_3^{(1)}$, the 29-plet $\Lambda^2\oplus\Lambda^1\oplus\Lambda^0$ (@sec-sm-c3). Their tensor products contain the same $U_q(b_n)$-representation more than once. Then $\check R$ restricted to each isotypic component is a matrix, not a number, and no graph rule applies.

Related fixes in ch. 8:
- **08:82, wrong example.** The 26⊗26 of U_q(f₄⁽¹⁾) = 1⊕26⊕52⊕273⊕324 is multiplicity-free, and DGZ 1994 eq. (4.64) solves it by the graph. Replace with: "For the $\mathbf{27}\otimes\mathbf{27}$ of $U_q(e_6^{(2)})$, where $\mathbf{27}=\mathbf{26}\oplus\mathbf1$ under $U_q(f_4)$ and the components $\mathbf1$ and $\mathbf{26}$ occur twice and three times, this reduces the problem to $2^2+3^2+1+1+1=16$ unknowns per value of the spectral parameter."
- **08:52.** Replace with: "For multiplicity-free tensor products its content is captured by a graph, introduced by MacKay [@mackay1991] for rational R-matrices and by Zhang, Gould and Bracken [@zgb1991] for trigonometric ones, and developed by Delius, Gould and Zhang [@dgz1994; @dgz1996]."
- **08:142.** Replace with: "extended to trigonometric ones by Zhang, Gould and Bracken [@zgb1991] and generalised by Delius, Gould and Zhang to pairs of different representations [@dgz1994] and to twisted algebras [@dgz1996]".
- **08:11.** Replace with: "By the mid-1990s the tensor-product graph had added many more, among them the spinor R-matrices of $U_q(b_n^{(1)})$ and $U_q(d_{n+1}^{(2)})$ [@dgz1994; @dgz1996]. For the soliton representations whose tensor products are not multiplicity-free they were not known, and constructing them was the main obstacle to finding those S-matrices (@sec-sm-introduction)."


## Item 11: chapter 7 statements

- **Central element (07:109): wrong.** With k_i e_j k_i⁻¹ = q_i^{a_ij}e_j, conjugation by ∏k_i^{m_i} multiplies e_j by q^{(Σm_iα_i, α_j)}. This is 1 for all j iff m_i = n_i, the Kac labels. The orchestrator re-derived this by hand. The printed ∏k_i^{n_i^∨} is not central for any non-simply-laced or twisted algebra; for d_{n+1}⁽²⁾ it is not even in the algebra, since the n^∨ are half-integers. Replace with:
  > The element $\prod_ik_i^{n_i}$, with the Kac labels $n_i$ of (-@eq-kac-labels), is central: conjugation by it multiplies $e_j$ by $q^{(\sum_in_i\vec\alpha_i)\cdot\vec\alpha_j}=1$. Since $k_i=q_i^{h_i}$ and $n_i^\vee=n_i\vec\alpha_i^2/2$, it equals $q^{\sum_in_i^\vee h_i}$, $q$ to the power of the central element. On the finite-dimensional representations used in this book it equals 1 (level zero).
- **Kirillov–Reshetikhin remark (07:119): terminology fix.** Replace "the minimal extensions (the **Kirillov–Reshetikhin modules**) can be larger" with "the smallest representations of $U_q(\hat g)$ that contain a given one (its **minimal affinizations**; for a multiple of a fundamental weight, the **Kirillov–Reshetikhin modules**) can be larger". The 27 = 26⊕1 example is consistent with the twisted KR rule.
- **Quasitriangularity consequences (07:72, 83–85): correct.** Checked with the explicit universal R on spin ½ and spin 1:
  - (Δ⊗1)R = R₁₃R₂₃ and (1⊗Δ)R = R₁₃R₁₂; the reversed orders fail;
  - the Yang–Baxter equation R₁₂R₁₃R₂₃ = R₂₃R₁₃R₁₂;
  - (S⊗1)R = R⁻¹ and (S⊗S)R = R.
- **07:111 (minor).** U_q(ĝ) is quasitriangular only after adjoining d and completing.
- **07:201 (minor).** For 2j+1 > p the zero-quantum-dimension summands are indecomposable but reducible (tilting modules). Suggested: "the summands with highest weight $j\ge(p-1)/2$ (irreducible for $2j+1=p$, where $[p]_q=0$, indecomposable but reducible for $2j+1>p$) have quantum dimension zero".
- **`exr-hermitian-spin1`: correct.** The form diag(1, [2]_q, [2]_q²) is unique up to scale. Its signature is (3,0) for cos γ > 0 and (2,1) for cos γ < 0; it is degenerate at q² = −1.

## Item 5: sine-Gordon results quoted but not derived (ch. 3)

Scripts: `referee_textbook_scripts/part1/sg_s0.py`, `sg_bootstrap.py`, `sg_residues.py`, `misc_checks.py`.

**`eq-sg-breather-masses`: correct.**
- DHN's m_n = (16m/γ′)sin(nγ′/16) is 2M sin(nξ/2) with ξ = (β_SG²/8)/(1−β_SG²/8π) and M = 8m₀/β_SG² − m₀/π.
- In book units this is ξ = πβ²/(4π−β²) and λ = 4π/β² − 1.
- "Exact at one loop" means that the factor 1/(1−β_SG²/8π) receives no higher corrections when M is the exact soliton mass. This is meaningful because β_SG is not renormalized.
- Optional clarification for 03:111: "…by the factor $1/(1-\beta_{\rm SG}^2/8\pi)$, which DHN found at one loop; the exact S-matrix shows it receives no higher corrections ($\beta_{\rm SG}$ is not renormalized, so $\xi(\beta_{\rm SG})$ is scheme-independent)."

**`eq-sg-s11` (S₁₁ = f_{ξ/π}): correct, derived by an explicit bootstrap** from the book's own S_T, S_R and S₀, in the Zamolodchikov–Faddeev slot convention.
- YBE holds to 1e-15.
- The residue at i(π−ξ) is rank 1, with image |ss̄⟩ − |s̄s⟩ (C-odd).
- S_{sB₁} = f_{1/2−ξ/2π} to 10 digits.
- Fusing four solitons gives S₁₁ = f_{ξ/π} to 10 digits at λ = 1.5 and at λ = 3.5. A wrongly ordered product is not invariant, so the test discriminates.

**Residue signs: no claim is printed in ch. 3, but ch. 2's sign rule is too strong.**
- Res S_T = (−1)ⁿ Res S_R, so B_n lies in the C = (−1)ⁿ channel.
- The sign of a soliton channel residue depends on the basis convention: rapidity-labelled S_T ± S_R give i(−1)ⁿ × positive, while position-labelled ones give all positive.
- The s-channel residue of S_{sB₁}, −2i cot(ξ/2), is negative in a unitary theory.
- At 1 < λ < 2 (e.g. λ = 3/2), S₁₁ = f_{2/3} has simple poles at 2πi/3 (negative residue) and iπ/3 (positive residue), and **neither is a bound state**. B₁ is C-odd, so there is no B₁B₁→B₁, and there is no particle of mass 3M.

Replace 02-integrable-qft.qmd:163, "In a unitary theory this sign tells the two interpretations apart", with:
> In a unitary theory this sign tells the two interpretations apart for a pole that is known to be a bound-state pole of a diagonal amplitude. Not every simple pole is a bound state: in sine-Gordon theory at $1<\lambda<2$ the breather amplitude $f_{\xi/\pi}$ has simple poles at $i\xi$ and $i(\pi-\xi)$ with no particle at either mass.

The same phenomenon occurs in c₂⁽¹⁾ S₂₂ (item 15, `exr-cn-snn`).

**`sec-sg-coleman` (03:181): correct but incomplete, and it overstates Part VII.**
- β_SG² < 4π for normal ordering of the correlation functions: correct.
- Neutral 2k-clusters diverge at 8π(1−1/2k): correct (power counting).
- Missing: collapses onto inserted cos βφ require an additive renormalization of the operator.
- 16π/3 is not a sine-Gordon threshold; it is the a₂⁽¹⁾ value.
- Part VII proves sufficiency below the first threshold, in finite volume. Above it, it uses power counting and a conjecture.

Replace 03:181 with:
> Below $8\pi$ the ultraviolet divergences are mild but not absent. Normal ordering of $\cos\beta_{\rm SG}\phi$ removes all divergences of the correlation functions for $\beta_{\rm SG}^2<4\pi$ [@frohlich1976]. A neutral cluster of $2k$ charges collapsing to a point diverges for $\beta_{\rm SG}^2\ge8\pi(1-\frac1{2k})$, $k=1,2,\dots$; the first case is a pair, $\int d^2y\,\lvert x-y\rvert^{-\beta_{\rm SG}^2/2\pi}$. In the vacuum energy these divergences require additive counterterms; in correlation functions of $\cos\beta_{\rm SG}\phi$ they also occur when integrated charges collapse onto an inserted one, and are removed by an additive renormalization of that operator. Clusters with nonzero net charge diverge only from $8\pi$ on, so below $8\pi$ finitely many field-independent counterterms suffice [@benfatto1982; Nicolò–Renn–Steinmann, Commun. Math. Phys. 105 (1986) 291]. In the language of @sec-collapse-thresholds these are the collapse thresholds of a two-dimensional Coulomb gas of positive and negative charges. Part VII proves, for every imaginary-coupling affine Toda theory in finite volume, that normal ordering suffices exactly below the first collapse threshold, and locates the higher thresholds by power counting.

**Crossing of `eq-sg-s0` at λ < 1: correct.**
- Direct integral at λ = 0.7 and 0.4: correct to 1e-21.
- Gamma-product continuation at λ = 0.4–3.7: correct to 1e-20, including 2×2 unitarity and S_R crossing.
- The continuation, which would be worth adding to `textbook_code/sine_gordon.py`, is S₀(θ) = −∏_{n≥0}[Γ(a_n+z)Γ(b_n−z)/(Γ(a_n−z)Γ(b_n+z))]^{(−1)ⁿ}, with a_n = 1+nλ, b_n = (n+1)λ, z = iλθ/π.
- **But the convergence strip at 03:138 is wrong for λ < 1.** At large t the integrand behaves like e^{(|Im θ|−min(ξ,π))t}; at λ = 0.7 and Im θ = 3.3 it grows (1.19, 14.1, 224 at t = 20, 40, 60). Replace "converges for $\lvert\operatorname{Im}\theta\rvert<\xi$" with "converges for $\lvert\operatorname{Im}\theta\rvert<\min(\xi,\pi)$".
- The book correctly calls λ < 1 the repulsive regime; the brief's "attractive" was a slip.


## Item 6: RSOS and Lee–Yang

- **`prp-sg-rsos`: correct.**
  - ξ = πp/(p′−p), with p < p′ required.
  - M(p,p+1) gives λ = 1/p.
  - (2,5) gives λ = 3/2 and one breather of mass √3M.
  - The citation should be smirnov1990 (item 1).
- **Lee–Yang as a restriction of a₂⁽²⁾ [Smirnov 1991] (04:140): correct.** In M(2,5), h₁,₂ = h₁,₃ = −1/5, so φ₁,₂ = φ₁,₃, and `smirnov1991` is the right reference there.
- **TCSA real "in the massive direction" (04:136): correct.** M(2,5) has no massless flow, so the phrase just means one sign of iλφ; Dorey–Tateo hep-th/9607167 confirms the claim. Optional addition: "for the opposite sign of the coupling, levels collide at finite volume and become complex-conjugate pairs".
- **Negative residue read as an imaginary coupling: correct.** Res_{2πi/3} f_{2/3} = −2√3 i, so Γ = i·2^{1/2}3^{1/4} (Cardy–Mussardo). This holds only after restriction (see item 5).
- **`exr-sg-lee-yang`: misleading last part.** Replace the last sentence with: "Locate the poles of $f_{2/3}$ in the physical strip. Why can the pole at $\theta=2\pi i/3$ not be a $B_1B_1\to B_1$ bound state in sine-Gordon theory (consider $\phi\to-\phi$), and what does it become after the restriction of @prp-sg-rsos?"
- **04:119 (cosmetic).** The Fisher Lagrangian is the Euclidean one; write "with Euclidean action density ½(∂φ)² + i(h−h_c)φ + igφ³".


## Item 7: chapter 2 statements

**`prp-cdd` (02:127–133): wrong (incomplete).**
- Unitarity and crossing force ±∏f_{x_k} with complex x_k. Real analyticity then requires the multiset {sin πx_k} to be closed under conjugation.
- The printed list ("real or ½+iy") misses:
  - x = −½+iy, i.e. Zamolodchikov's staircase S-matrix, the book's own sinh-Gordon S-matrix at complex B. It equals 1/f_{½+iy}, so the list is not even closed under inverses.
  - pairs f_x f_x̄.
- The growth bound must be uniform in Im θ, and it is genuinely needed: e^{ia sinh θ} satisfies the three functional equations.
- The orchestrator checked the staircase case by hand.

Replacement:
```
::: {#prp-cdd}
**[Established]** Let $S(\theta)$ be meromorphic, with $\lvert S(\theta)\rvert\le C\,e^{N\lvert\operatorname{Re}\theta\rvert}$ for all sufficiently large $\lvert\operatorname{Re}\theta\rvert$, uniformly in $\operatorname{Im}\theta$. If $S(\theta)S(-\theta)=1$, $S(i\pi-\theta)=S(\theta)$ and $S(\theta)^*=S(-\theta^*)$, then
$$
S(\theta)=\pm\prod_{k=1}^{N}f_{x_k}(\theta)
$$
for finitely many complex $x_k$ such that the numbers $\sin\pi x_k$, counted with multiplicity, are permuted by complex conjugation: each $x_k$ is real, or has real part $\pm\tfrac12$, or occurs together with $\bar x_k$.
:::
```
Also at 02:135, "multiplying a solution by any $f_x$" → "…any $f_x$ with $\sin\pi x$ real, or by a pair $f_xf_{\bar x}$". Optionally add: "The growth condition cannot be dropped: $e^{ia\sinh\theta}$ satisfies all three conditions."

**`eq-bootstrap` (02:165–169): correct, and consistent with Part IV and `smatrix_simply_laced.py`.**

| algebra | printed form | sign-flipped form | angles swapped |
|---|---|---|---|
| a₂ | 8/8 | 8/8 | 8/8 |
| a₄ | 96/96 | 96/96 | 32/96 |
| e₆ | 720/720 | 720/720 | 192/720 |
| e₈ | 3584/3584 | 3584/3584 | 512/3584 |

The sign-flipped form is equivalent to the printed one by unitarity. The swapped assignment passes only for isosceles triangles.

**02:154, angle sum: index error.** Replace `u_{ab}^c+u_{bc}^{\bar a}+u_{ca}^{\bar b}=2\pi` with `u_{ab}^{c}+u_{b\bar c}^{\bar a}+u_{\bar c a}^{\bar b}=2\pi`. The two agree for self-conjugate particles.

**`exr-no-production`, `exr-free-fermion`: true and solvable.**

Other ch. 2 errors:
- **02:50, conserved charge: wrong sign.** ∂₋T = ∂₊Θ gives ∂₀(T−Θ) = ∂₁(T+Θ). The agent checked this on a free-field solution (∫(T+Θ) is time-dependent), and the orchestrator re-derived it by hand. Fix: `Q_s=\int dx^1\,\big(T_{s+1}-\Theta_{s-1}\big)`, or keep the + and write the conservation law as ∂₋T + ∂₊Θ = 0.
- **02:93, braiding unitarity: index order.** With the printed ZF relation, S_ab^{cd}(θ)S_dc^{fe}(−θ) = δ. The printed version is right only for parity-invariant S, which covers every S-matrix in the book. Either fix the indices or add "for parity-invariant S".
- **02:104: overclaim.** The sine-Gordon soliton S-matrix *is* Hermitian analytic. Write "Apart from sine-Gordon, the S-matrices of the imaginary-coupling theories … are in general not Hermitian analytic".
- **02:193.** "…and so is the energy" → "and the energy can be complex".


## Item 8: chapter 4 statements

**Mostafazadeh claim (04:68): wrong as stated.**
- Counterexample: diag(i, i, −i). Its spectrum is closed under conjugation as a set, but H is not similar to H†.
- The theorem (J. Math. Phys. 43, 205) needs equal multiplicities of conjugate pairs.
- The antilinear-symmetry statement is J. Math. Phys. 43, 3944.
- In finite dimensions, pseudo-Hermitian ⇔ H ~ H†; diagonalizability is not needed, and this matters for the dimer at its exceptional point.

Replace the second sentence of 04:68 with:
> An operator with discrete spectrum and a complete biorthonormal system of eigenvectors (in finite dimensions, a diagonalizable matrix) is pseudo-Hermitian if and only if its non-real eigenvalues come in complex-conjugate pairs of equal multiplicity [@mostafazadeh2002], and under the same hypotheses pseudo-Hermiticity is equivalent to the existence of an invertible antilinear symmetry [@mostafazadeh2002c]. In finite dimensions diagonalizability is not needed: a matrix is pseudo-Hermitian if and only if it is similar to its adjoint, which allows Jordan blocks such as the dimer of @sec-dimer at its exceptional point.

**`exr-krein-collision` (04:165–167): false as stated.**
- Counterexample: H = 0, J = diag(1,−1), H′ = [[0,1],[−1,0]], which is J-self-adjoint. Then εH′ has eigenvalues ±iε, even though x = (1,0) has [x,x] = 1.

Replacement:
> Let $H$ be $J$-self-adjoint on a finite-dimensional Krein space, and let $\lambda$ be a real eigenvalue of **positive type**: $[x,x]>0$ for every nonzero $x$ in the root subspace $\ker(H-\lambda)^n$. Show that for every $J$-self-adjoint $H'$ and all sufficiently small $\varepsilon$, the eigenvalues of $H+\varepsilon H'$ near $\lambda$ are real. (*Hint:* a non-real eigenvalue requires a neutral eigenvector, and the eigenvectors of $H+\varepsilon H'$ for eigenvalues near $\lambda$ lie close to the root subspace.) Show with $H=0$, $J=\operatorname{diag}(1,-1)$, $H'=\begin{pmatrix}0&1\\-1&0\end{pmatrix}$ that a single eigenvector with $[x,x]>0$ is not enough.

**04:112: wrong as written.** Same-signature eigenvalues can collide; they stay real. Replace with:
> Real eigenvalues can leave the real axis only by colliding, and in a Krein space a collision can push them off the axis only if the colliding eigenvectors have *opposite* Krein signature (more precisely, only if the root subspace at the collision is not definite); eigenvalues of the same signature pass through each other and stay real.

**`sec-krein-infinite` (04:157): mostly fine, with imprecise wording.**
- Add that the spectrum of a J-self-adjoint operator can be all of ℂ.
- Replace the spectral-singularity definition with: "points $k^2$ of the continuous spectrum (real $k\ne0$) where the Jost function vanishes, at which the spectral resolution breaks down although there is no eigenvalue".


## Other Part I notes

- **Notation clash.** 30-continuum-thresholds.qmd:35 defines ξ = β²/(4π−β²), which is ch. 3's ξ divided by π, and uses m for the soliton mass. Rename one of them (e.g. ν = ξ/π in Part VII) or add a sentence relating them.
- **Checked and correct:** ch. 1 conventions; the ch. 3 potential and dictionary; kink, τ-form and fluctuation operator; T(k) = (k+im₀)/(k−im₀); breather rapidities; B and its duality; q = ∓e^{iπλ}; `prp-krein-spectrum`; the PT discussion; Θ = PK; c = −22/5 and h = −1/5.

## Item 9: chapter 5 statements

Scripts: `referee_textbook_scripts/lie_classical/lie.py`, `ch5_check.py`, `ch5_charges.py`, `ch5_charges_b.py`, `e8_ising.py`.

- **`tbl-simple-lie`, eq-coxeter-numbers, Steinberg (`sec-weyl-coxeter`): correct.** Computed for A2–A6, B2–B5, C3–C5, D4–D7, E6–E8, F4 and G2:
  - d = r(h+1); h and h^∨; exponents; Σm_k = |Φ⁺|;
  - Coxeter elements in all orders have order h and the same characteristic polynomial;
  - the eigenvalues are e^{2πim_k/h};
  - the bicoloured element rotates the Coxeter plane by 2π/h, with r orbits of size h projecting onto regular h-gons, whose radii equal the Perron–Frobenius vector for ADE;
  - det A = |centre|; C(θ) = 2h^∨; ρ = Σλ_i.
  - Steinberg 1959 should be cited (item 1).
- **Higher-spin statement (`sec-perron-frobenius`, 05:164): correct, with a clarification recommended.**
  - Method: all three-point conservation laws Σq_x e^{isφ_x} = 0 were solved at Dorey fusings for A4, D4–D6 and E6–E8.
  - For every exponent s, the solution space equals the Cartan eigenspace with eigenvalue 4sin²(πs/2h). This includes the 2-dimensional D4 (s = 3) and D6 (s = 5) cases.
  - D4 at s = 2, 4 has a harmless spurious three-point-only solution.
  - For s > h, the eigenvalue is still 4sin²(πs/2h).
  - Replace "The other eigenvectors ... spin $s=m_k$ on the particles, in the sense of (-@eq-charge-eigenvalue) [@freeman1991]." with:
    > The other eigenvectors of the Cartan matrix are also physical: for every spin $s$ of a conserved charge, the eigenvalues $\chi_a^{(s)}$ of (-@eq-charge-eigenvalue) form an eigenvector of the Cartan matrix with eigenvalue $4\sin^2\frac{\pi s}{2h}$. For $s=m_k$ this is the eigenvector of the exponent $m_k$; when an exponent is repeated ($D_{2n}$, $s=2n-1$) the charges of that spin span the two-dimensional eigenspace [@klassen1990; @dorey1991].
- **E8/Ising sentence (05:162): correct as worded.**
  - S_min = ∏(x−1)(x+1) over {1, 11, 19, 29} equals f_{2/3}f_{2/5}f_{1/15} exactly.
  - The Toda S₁₁ is S_min times a B-dependent CDD factor.
  - The B → 0 limit of the Toda S-matrix is 1, not S_min, so keep the word "minimal" and do not rephrase it as "B → 0 limit".
- **Ch. 5 exercises: all true and solvable.**
- **Cosmetic, 05:70.** sp(2n) is not the compact real form; write usp(2n).


## Item 10: chapter 6 statements

Scripts: `fold.py`, `centre.py`, `twisted_exponents.py`, `duality_masses.py`; OTU text in `otu.txt`.

**All folding rows of `tbl-foldings` and the ch. 19 table: correct.** Method:
- exhaustive search of affine-diagram automorphisms;
- the book's rule α_O = orbit average, with labels |O|n_O;
- the folded Cartan matrix compared with the target and with its transpose.

Results:
- d_{n+1}⁽²⁾ ← d_{n+2}⁽¹⁾ (n = 2–4), a_{2n−1}⁽²⁾ ← d_{2n}⁽¹⁾ (n = 2–4), e₆⁽²⁾ ← e₇⁽¹⁾ and d₄⁽³⁾ ← e₆⁽¹⁾ all give exactly the target, not the transpose.
- a_{2n}⁽²⁾ ← d_{2n+2}⁽¹⁾ (n = 1–3) needs the **order-4** automorphism 0→N−1→1→N→0. The text could say so.
- Ranks, h and Kac h^∨ are all correct.

Errors:
- **E1, 06:111: wrong.** It says soliton masses of an untwisted non-simply-laced theory follow the particle masses of "the *dual* theory", where dual means arrow reversal, ĝ^∨. They actually follow (g^∨)⁽¹⁾. This was confirmed by two agents and spot-checked by the orchestrator: c₃⁽¹⁾ solitons (1, 1, √3) match b₃⁽¹⁾ particles, not d₄⁽²⁾ (1, 1.848, 2.414). Ch. 10 and appendix A are right. Replace with:
  > - The masses of the solitons of an untwisted non-simply-laced theory $g^{(1)}$ are proportional to the particle masses of $(g^\vee)^{(1)}$, the untwisted algebra of the dual root system ("Lie duality", @sec-imaginary-coupling-vacua). This is a different duality from $\hat g\to\hat g^\vee$: for $c_n^{(1)}$ it gives $b_n^{(1)}$, not $d_{n+1}^{(2)}$.

  The same conflation is at 14:53 (item 15); ch. 19 uses "Lie dual" in the (g^∨)⁽¹⁾ sense. Fix the term book-wide.
- **E2, 06:142: wrong for twisted ĝ.** Principal degrees, and hence conserved spins, repeat with period **kh**, not h. Confirmed independently by the Part III agent:
  - a₂⁽²⁾: 1, 5 mod 6 (Bullough–Dodd; with h = 3 the book's rule would allow spin 2);
  - a₄⁽²⁾: 1, 3, 7, 9 mod 10;
  - d₄⁽³⁾: 1, 5, 7, 11 mod 12;
  - e₆⁽²⁾: 1, 5, 7, 11, 13, 17 mod 18;
  - a_{2n−1}⁽²⁾ and d_{n+1}⁽²⁾: all odd.

  Replace the last two sentences of 06:142 with:
  > …is spanned by elements whose principal degrees are the exponents of $\hat g$, which repeat with period $h$ for untwisted $\hat g=g^{(1)}$ and with period $kh$ for twisted $\hat g=g^{(k)}$ (for example $1,5$ mod 6 for $a_2^{(2)}$, $1,5,7,11$ mod 12 for $d_4^{(3)}$) [@kac1990]. This is why the conserved charges of affine Toda theory have spins equal to the exponents, repeated with this period.

  The same fix is needed at 09:32, 09:84–88 (item 13) and 02:59.
- **E3, 06:107 caption with 06:61: incomplete.** The formula n_j^∨ = n_jα_j²/2 (longest root² 2) gives h^∨/k for twisted algebras. Add to the caption:
  > For a twisted algebra $g^{(k)}$, $h=\sum_jn_j$ with coprime integer labels, and $h^\vee=\sum_jn_j^\vee$ with $n_j^\vee=k\,n_j\vec\alpha_j^2/2$, i.e. the dual labels in the normalization where the longest root has squared length $2k$ (equivalently, the coprime integer null vector of the transposed Cartan matrix).
- **E4, 06:89: normalization gap.** The twisted folding labels |O|n_O have common factor k. After "...whose labels are $\lvert O\rvert n_O$." add:
  > For the twisted foldings ($\alpha_0$ in an orbit of size $k$) these labels have the common factor $k$: the folded theory is the standard $g^{(k)}$ theory, with coprime labels $n_j$ and $h=\sum_jn_j$, at the same $\beta$ and with $m^2$ replaced by $km^2$.
- **E5, 06:156–160: the OTU τ_j is printed without x, t dependence.** OTU (2.27)–(2.29) multiply each F^{a_k}(z_k) by W_k. Replace the display with
  `\tau_j=\langle\Lambda_j\rvert\,\prod_k\exp\big(Q_kW_k(x,t)F^{a_k}(z_k)\big)\,\lvert\Lambda_j\rangle,\qquad W_k=\exp\big(-m\mu_{a_k}(z_kx^+-z_k^{-1}x^-)\big)`
  and write "$z_k$ encodes the rapidity and $Q_k$ the position and internal phase". The sign conventions inside W follow OTU only up to the book's normalization of μ_a and m. The same fix applies to ch. 10 `sec-vertex-solitons`. OTU treat untwisted ĝ only; for twisted algebras cite Kneipp–Olive hep-th/9404030. The general principal construction is due to Kac–Kazhdan–Lepowsky–Wilson.
- **E6, 06:63: imprecise.** Replace "the symmetries of the affine diagram that move the extra node correspond to the centre of the simply connected group" with:
  > the symmetry group of the affine diagram is the semidirect product of the symmetries of the finite diagram (which fix the extra node) with a subgroup isomorphic to the centre of the simply connected group, which permutes the nodes with $n_j=1$ simply transitively.

  For example, d₄ has 18 automorphisms that move node 0, against a centre of order 4. "Charge-conjugation-like" → "act on the particles by phases".
- **E7, 06:113: normalization.** β → 4π/β holds only with ĝ^∨ in the literature normalization (longest root² 2k). In the book's normalization it is β^{∨2} = 16π²k/β². See also item 15: for c_n the clean form holds in the short-root normalization.
- **Correct:**
  - "the level n_j^∨ of L(Λ_j)" (OTU: n_i = 2m_i/a_i²);
  - [E_{±1}, F^a] = μ_a z^{±1}F^a in structure;
  - the nilpotency bound for simply-laced algebras;
  - gradations, the E₁ coefficients and the twisted algebras paragraph.
- **Ch. 6 exercises: all true and solvable.**


## Item 20a: appendix A

- **`tbl-app-twisted` and `tbl-app-untwisted`: correct.** Rank, h, h^∨, folding parents and duals were all verified with `fold.py`, and they are consistent with `tbl-foldings` and the ch. 19 table. Suggest adding to the caption: "h and h^∨ are computed with coprime labels and $n_j^\vee=kn_j\vec\alpha_j^2/2$ (@sec-affine-folding); the spins of the conserved charges repeat with period $kh$."
- **App. A:49: vague.** The twisted masses (coprime labels) are explicit:
  - a_{2n}⁽²⁾: 2m sin(πa/h), a = 1..n;
  - a_{2n−1}⁽²⁾: 2m sin(πa/h) for a = 1..n−1, and m;
  - d_{n+1}⁽²⁾: 2m sin(πa/2h);
  - d₄⁽³⁾: ratio 2cos(π/12);
  - e₆⁽²⁾: 0.879385, 1.347296, 1.732051, 2.532089 (×m).

  The closed forms are fitted for n ≤ 4 and not proved. The a_{2n−1}⁽²⁾ masses were confirmed independently by the Parts IV–V agent.
- **`sec-app-soliton-data`:**
  - l.55, l.58 and l.60 are consistent with the main text.
  - l.56 is consistent with ch. 10 but contradicts 06:111 (fix ch. 6). Add the twisted case: "proportional to the particle masses of the same theory for twisted $\hat g$".
  - l.57 is **unclear**: whether β or β_lit enters q for twisted ĝ. The q-convention referee found q_BL = e^{−4π²i/β²} with x_c = (−1)^{h^∨}q_BL^{−kh}, which suggests β in book normalization with explicit factors of k. This needs one sentence of convention.
  - 21:25–26 `sec-amplitude-and-gradation` states 2T without the β_lit proviso; add it (item 19).

## Item 13: Part III statements

Scripts: `referee_textbook_scripts/part3/s1_masses_dorey.py` … `s9_exercises.py`. They use copies of `toda_classical.py` and `foundations_code/hirota/{lie,hirota}.py`.

- **`eq-hirota-bilinear` for non-simply-laced and twisted algebras: correct.**
  - Re-derived: n_j are the marks, and n_jα_j²/2 the comarks.
  - Checked numerically by folding explicit parent solutions and substituting into the book's equation, with residuals ~1e-40:
    - c₂⁽¹⁾ ← a₃⁽¹⁾;
    - a₂⁽²⁾ ← a₂⁽¹⁾ (non-orthogonal orbit, τ_O = τ^{1/2});
    - d₃⁽²⁾ ← d₄⁽¹⁾;
    - g₂⁽¹⁾ ← d₄⁽¹⁾.
  - Minor fix at 10:51: "…counted with multiplicity $-a_{lj}$ (for $a_1^{(1)}$ it is $\tau_{1-j}^2$)". The sine-Gordon kink at 10:61 needs τ₁².
- **`prp-dorey-rule`, `eq-area-rule` with degenerate masses: correct, but only in the charge basis. The text's basis prescription is wrong.**
  - Checked for a₄, a₅, d₄, d₅, d₆, e₆, e₇ and e₈, in the basis of centre-character eigenvectors with e_a·e_ā = 1:
    - C_abc ≠ 0 iff Dorey's rule holds;
    - |C_abc| = (4β/√h)Δ_abc, with spread below 1e-8;
    - |projected root|/m_a is constant;
    - all angles are multiples of π/h.
  - In a real orthonormal basis both rules fail. For e₆ the ratios |C|/((4/√h)Δ) come out as {0.7071, 1}.
  - At 09:120, replace "in an orthonormal eigenbasis of $M^2$" with "in an eigenbasis of $M^2$", and add after the display:
    > For degenerate masses ($a_n$, $d_n$, $e_6$) the $\vec e_a$ must be charge eigenvectors, the complex eigenvectors of the centre of the affine diagram (of $\mathbb Z_h$ for $a_n$), normalized by $\vec e_a\cdot\vec e_{\bar a}=1$ so that $\phi_{\bar a}=\bar\phi_a$. Particle $a$ is the eigenvector with the centre character of node $a$. In a real orthonormal basis of a degenerate pair neither (-@eq-area-rule) nor @prp-dorey-rule holds.
  - At 09:152, extend "can be checked directly" with: "…against the cubic expansion of the potential: for $e_7^{(1)}$, $e_8^{(1)}$ in any eigenbasis, and for $a_n^{(1)}$, $d_n^{(1)}$, $e_6^{(1)}$ in the charge basis".
- **Spins of conserved charges, twisted algebras: wrong** (period kh; see E2 in item 10).
  - At 09:32, replace "conserved charges of spins equal to the exponents of $\hat g$ modulo the Coxeter number [@olive1985]" with "conserved charges whose spins are the exponents of $\hat g$, repeating with period $kh$ ($h$ the Coxeter number; $k=1$ for untwisted and $k=2,3$ for twisted $\hat g=X^{(k)}$) [@olive1985]".
  - At 09:84–88, change the display to `s\in\{\text{exponents of }\hat g\}+kh\,\mathbb Z ,` and replace the following sentence with:
    > where $h=\sum_jn_j$ and $kh$ is the order of the twisted Coxeter element ([@kac1990] lists the exponents of $X_N^{(k)}$ modulo $kh$). For untwisted $\hat g=g^{(1)}$ these are the exponents of $g$ (@tbl-simple-lie) repeated modulo $h$. For $a_2^{(2)}$ they are $1,5$ mod 6; for $d_4^{(3)}$, $1,5,7,11$ mod 12; for $e_6^{(2)}$, $1,5,7,11,13,17$ mod 18; for $a_{2n}^{(2)}$, the odd integers not divisible by $2n+1$; for $a_{2n-1}^{(2)}$ and $d_{n+1}^{(2)}$ (as for $b_n^{(1)}$, $c_n^{(1)}$), all odd integers. The two members of a dual pair have the same spins.
  - The untwisted non-simply-laced statement is fine.
  - The Drinfeld–Sokolov attribution needs a citation (item 1).
- **e₈ mass ratios (09:116): correct.** All eight match the closed forms to 1e-10, and Σm² = 2h.
- **`eq-soliton-mass` for d_n, e_n and folded solitons: correct.**
  - Direct energy integration at complex ξ gives 2hm_a/β² to 16 digits for every species of d₄, d₅ and e₆.
  - General argument, suggested for 10:92: E = 2∫V = (2/β²)Σ_j(2/α_j²)[∂_x ln τ_j] = (2/β²)μ_aΣ_j(2/α_j²) deg τ_j.
  - Twisted case: "masses of the same theory" (10:24) is correct (d₄⁽³⁾, d₃⁽²⁾, d₄⁽²⁾, a₅⁽²⁾, e₆⁽²⁾).
  - Untwisted Lie duality c_n ↔ b_n is correct; 06:111 contradicts it (item 10).
- **McGhee count h/gcd(a,h): correct, with a general proof.** The charge t satisfies t_k + a/h ∈ {0,1}, so it is a weight of Λ^a. With g = gcd(a,h), the phases jump at h/g points, g at a time, giving h/g distinct charges. Checked numerically for n ≤ 12.
  - **10:86 is wrong when gcd > 1.** Replace "the reduced phases jump one at a time, and the charge moves through a sequence of weights" with "as $\operatorname{Im}\xi$ increases through $2\pi$, the reduced phases jump at $h/\gcd(a,h)$ values, $\gcd(a,h)$ of them at once (those with $j$ in one residue class mod $h/\gcd(a,h)$), and the charge runs through a cycle of $h/\gcd(a,h)$ weights".
- **Time delays (10:118): sign and form correct; the semiclassical claim was verified for sine-Gordon only.**
  - Δx = −ln A/(m_{a₁}coshθ₁) for θ₁ > θ₂.
  - Since 0 < A < 1, the shift is forward, i.e. a time advance.
  - dδ/dθ = (2h/β²) ln A matches the book's S₀ as ξ → 0 (λ = 160: −1.40670 against −1.40683).
  - Add "For $\theta_1>\theta_2$" and "…and positive, a time advance. Semiclassically this means $d\delta_{ab}/d\theta=(2h/\beta^2)\ln A_{ab}(\theta)$ at leading order; for sine-Gordon this is the leading term of $-i\ln S_0$ in (-@eq-sg-s0)."
  - For a_n with n ≥ 2 the claim is **unclear**. Settle it by expanding Hollowood's a_n scalar factor to O(1/β²).
- **Vertex-operator nilpotency (10:75): overstated, and wrong for non-simply-laced algebras.**
  - c₂⁽¹⁾ coincident pair: deg τ₀ = 2 with n₀ = 1.
  - g₂⁽¹⁾ orbit soliton: degrees (3, 3, 6).
  - In general deg = |O|n_j^∨.
  - Replace the bullet with:
    > **Polynomial structure.** For simply-laced $\hat g$, $(F^a)^{n_j+1}=0$ on $L(\Lambda_j)$, which has level $n_j$ [@olive1993b], so a single-soliton $\tau_j$ has degree at most $n_j$ in $E$. That the degree is exactly $n_j$ is explicit for $a_n^{(1)}$ and checked numerically for $d_{4..8}^{(1)}$, $e_{6,7,8}^{(1)}$ (@sec-hirota-form-dn1-en1). For folded theories the relevant label is $n_j^\vee$ (the level), and a coincident orbit $O$ has degree $\lvert O\rvert n_j^\vee$, e.g. 2 at the long end nodes of $c_n^{(1)}$. Bounds on products of different $F$'s are hypothesis (D) of @sec-hirota-form-dn1-en1, which is verified only numerically.
  - The same unqualified claim is in 06 `sec-vertex-operators`.
  - The general OTU proof was not independently checked.
- **`exr-breather-real`: misleading.** With ξ₂ = ξ̄₁ the τ_j are real but change sign; for a₂ at u = 0.6, min τ = −2.2e5. Replace "show that $\bar\tau_j=\tau_j$ and that the field is purely imaginary" with "show that $\bar\tau_j=\tau_j$, so that the field is purely imaginary where all $\tau_j>0$, and that the $\tau_j$ change sign, so this solution is singular along curves in the $(x,t)$ plane".
- **`exr-singular-points`: true for generic data.** The double-pole coefficient is −(2/β²)(2/α²), checked numerically. At 11:57, write "are generically isolated points (non-generic data, such as a single soliton with tuned $\operatorname{Im}\xi$ or a breather with $\xi_2=\bar\xi_1$, give curves of zeros)".

Other Part III errors:
- **10:10: false.** The a₃⁽¹⁾ species-2 soliton at Im ξ = π/2 is real. The orchestrator checked this by hand: τ₁ = τ₃ = τ̄₀, τ₂ = τ₀, so φ = −(2/β)(α₀+α₂) arg τ₀, an embedded sine-Gordon kink. Replace with:
  > Unlike the sine-Gordon kink, these solitons are in general complex: a field that tends to vacua at both ends and solves the field equations typically takes complex values in between. (Exceptions are embedded sine-Gordon kinks, e.g. the species-$n$ soliton of $a_{2n-1}^{(1)}$ at $\operatorname{Im}\xi=\pi/2$, which has $\tau_{2k+1}=\bar\tau_{2k}$ and is real.)
- **11:43, "all conserved charges real": false.**
  - In a₂⁽¹⁾ the spin-2 charge is ∓3.1686, real, on the single solitons at θ = 0.3.
  - On a breather it is 6.5131i, imaginary.
  - No normalization makes it real on both.

  Replace with:
  > **Selection rule.** If the asymptotic data are $\Theta$-invariant, i.e. the multiset $\{(a_i,\theta_i)\}$ is mapped to itself by $(a,\theta)\mapsto(\bar a,\bar\theta)$, then every conserved charge satisfies $Q=\epsilon_Q\bar Q$, where $J_Q(-\bar{\vec\phi})=\epsilon_Q\overline{J_Q(\vec\phi)}$ defines $\epsilon_Q=\pm1$; in particular $E$ and $P$ are real. Breathers are of this kind. A single soliton with real rapidity is $\Theta$-invariant only if $\bar a=a$; for $\bar a\ne a$ its $E$ and $P$ are real because $M_{\bar a}=M_a$, but charge-conjugation-odd charges behave differently: in $a_2^{(1)}$ the spin-2 charge, normalized to be real on single solitons, is imaginary on breathers.

  Adjust the example at 11:51 to match, and fix 11:41, which writes (a,θ,p) where `eq-theta-data` uses ξ.
- **11:26: garbled.** Replace with:
  > A $\Theta$-invariant field satisfies $\vec\phi=-\bar{\vec\phi}$ (up to a vacuum shift), i.e. it is purely imaginary, and $-i\vec\phi$ solves the real-coupling equation; it carries no topological charge. So no soliton is $\Theta$-invariant as a field, and for solitons $\Theta$-invariance is a property of the asymptotic data.
- **Minor:**
  - 09:100: "invisible in the classical Lax pair" → "not among the local charges of the Lax pair".
  - The breather parameter u means different things in ch. 3 (±i(π/2−u)) and ch. 10 (θ₀ ± iu).
- **Checked and correct:** the three-soliton eq-n-soliton for a₃ and a₄ (residual 3e-29).

Exercises of ch. 9–11 (17):
- **All true**, with the following notes:
  - `exr-moving-soliton` should read "½(D_x²−D_t²)τ_j·τ_j = m_a²(cosh²θ−sinh²θ)ω^{ja}E".
  - `exr-theta-solutions` needs "(unless $\hat g=a_1^{(1)}$, where $\vec\phi\mapsto-\vec\phi$ is itself a symmetry)".
  - `exr-breather-real` is misleading (above).
- Values computed:
  - `exr-spin3`: Q₃ = −16m³/(3β²) on the static kink;
  - `exr-shg-cubic`: C₁₁₁ = −3iβm², C₁₁₂ = 0, and |C₁₁₁| = (4β/√3)Δ exactly.

## Item 15: Part IV statements

Scripts: `referee_textbook_scripts/part45/floating_check.py`, `f4_check.py`, `hollowood_check.py`, `hollowood_continuum.py`, `poles_check.py`. These compute one-loop bubble self-energies, in the normal-ordered scheme, from the book's Lagrangian. The orchestrator re-ran `floating_check.py` and reproduced the output.

### `eq-floating-h` (14:41–47): wrong as a general formula

Book normalization (long roots of length² 2), B = β²/2π + O(β⁴). Each δH was fitted from every mass ratio separately, and all ratios agree:

| ĝ | one-loop H | `eq-floating-h` gives |
|---|---|---|
| b_n⁽¹⁾ | **2n − B/2** (decreasing to h^∨ = 2n−1; DGZ §4, CDS §6) | 2n + (n−1)B ✗ |
| c_n⁽¹⁾ | 2n + B/2 | 2n + B (✓ only in DGZ's short-root normalization) |
| g₂⁽¹⁾ | 6 + B | 6 + 3B ✗ |
| f₄⁽¹⁾ | 12 + 3B/2 | 12 + 3B ✗ |

- The earlier by-hand "verification" for c_n and f₄ holds only when B is read in the normalization where the **short** roots have length² 2. That is how DGZ write c_n and CDS write the exceptional pairs, and it is the block parameter used in Part VI.
- b_n fits in no normalization.
- The twisted endpoint of a_{2n−1}⁽²⁾ is H = 2n−1 = h^∨(b_n), not kh^∨.

Replace 14:41–47, from "The floating Coxeter number interpolates linearly in $B$." through "…are in [@dgz1992; @corrigan1993].", with:
> The floating Coxeter number interpolates between the Coxeter number of $\hat g$, as it appears in the classical mass formula, and that of the dual theory $\hat g^\vee$. To one loop, with $B=\beta^2/2\pi+O(\beta^4)$ in the normalization of @sec-conventions (long roots of length$^2$ 2), the one-loop masses of the Lagrangian give
>
> | $\hat g$ | $H$ | range | dual endpoint |
> |---|---|---|---|
> | $c_n^{(1)}$ | $2n+\tfrac12B$ | $2n\le H\le 2n+2$ | $d_{n+1}^{(2)}$, masses $\sin\frac{\pi a}{2n+2}$ |
> | $b_n^{(1)}$ | $2n-\tfrac12B$ | $2n-1\le H\le 2n$ | $a_{2n-1}^{(2)}$, masses $2\sin\frac{\pi a}{2n-1}$, spinor 1 |
> | $g_2^{(1)}$ | $6+B$ | $6\le H\le12$ | $d_4^{(3)}$ |
> | $f_4^{(1)}$ | $12+\tfrac32B$ | $12\le H\le18$ | $e_6^{(2)}$ |
>
> In the normalization in which the *short* roots have length$^2$ 2 (coupling $\beta^2/k$), the entries for $c_n^{(1)}$, $g_2^{(1)}$ and $f_4^{(1)}$ become $2n+B$, $6+3B$ and $12+3B$, which is how Delius, Grisaru and Zanon write $c_n^{(1)}$ and Corrigan, Dorey and Sasaki write the exceptional pairs. $b_n^{(1)}$ does not fit this pattern: its $H$ decreases, from $h=2n$ to $h^\vee=2n-1$. For $f_4^{(1)}$ the masses also involve $H'$ with $1/H+1/H'=1/6$ [@corrigan1993]. The relation $H-12=3B_{\rm tree}$ confirmed in @sec-sm-f4 refers to the parameter of the CDS blocks, which is $\beta^2/4\pi$ in our normalization.

Consequential changes:
- `exr-floating-endpoints` (14:72): "Check (-@eq-floating-h) at both ends" → "Check the $c_n^{(1)}$ entry of the table at both ends". The d_{n+1}⁽²⁾ part is correct for n = 3, 4.
- **Part VI is affected.** Any use of eq-floating-h in the book's own normalization must be rechecked. ch. 23 `sec-sm-cn` uses H = 2n + B, which is consistent if its B is the block parameter.

### Other Part IV items

**`sec-nsl-duality` (14:51–53): partly wrong.**
- 14:51: "a theory with Coxeter number $kh^\vee$" is false for b_n. Replace with "a theory whose classical masses are those of the dual theory $\hat g^\vee$ (the dual endpoints in the table above)".
- "One family describes both members" is correct as a conjecture: S^{c_n}(B) = S^{d_{n+1}⁽²⁾}(2−B), and similarly for b_n.
- The "longest roots of length² 2k" clause holds only in the normalization where B(β) is exact; for c_n that is the short-root one. Drop the clause or state it per pair.
- Same caveat at 13:95: add "(in the normalization in which $B(\beta)$ has the form (-@eq-b-of-beta); for $c_n^{(1)}$ the short roots have length² 2)".
- 14:53 conflates the two dualities (as at 06:111). Replace with "A related but different Lie duality appears in the soliton masses: the solitons of $X^{(1)}$ have the particle mass ratios of $(X^\vee)^{(1)}$ (@sec-soliton-masses)."

**`sec-nsl-masses`: correct.** The one-loop c_n masses ∝ sin(πa/H) (bubble computation for c₃ and c₄), and the folding heuristic is fine as a heuristic.

**`sec-nsl-status` (14:59): overclaim.** DGZ checked masses for a_{2n−1}⁽²⁾, b_n, c_n, d_{n+1}⁽²⁾ and g₂, tree amplitudes, and a one-loop residue in two cases; CDS made no perturbative check for f₄. Replace with:
> **One loop.** The floating Coxeter number reproduces the one-loop masses of $a_{2n-1}^{(2)}$, $b_n^{(1)}$, $c_n^{(1)}$, $d_{n+1}^{(2)}$ and $g_2^{(1)}$, and their tree amplitudes, and a one-loop residue fixes $B(\beta)$ to $O(\beta^4)$ [@dgz1992]. For $f_4^{(1)}$ the one-loop masses are checked in @sec-sm-f4.

**`exr-cn-snn` (14:68): misleading.** For n = 2 the double poles are only partly cancelled, so S₂₂ has two simple poles at 2πi/H and iπ(1−2/H). Their residues change sign at B = 1 (−i·Res = +0.066 … 0 … −0.024), and c₂ has no 22c coupling, so they are not bound states. Replace the second sentence with:
> For $n=2$, show that the only poles in the physical strip are simple poles at $\theta=2\pi i/H$ and $i\pi(1-2/H)$, remnants of double poles partly cancelled by CDD zeros, that their residues change sign at $B=1$, and that they are therefore not bound states [@corrigan1993].

**`sec-coleman-thun`:**
- The order P − 2L is correct (scaling argument; box 2, tree 1).
- **12:81 overclaims.** Only double and triple poles have been explained in detail. Replace with:
  > Double and triple poles have been explained in detail by such on-shell diagrams [@braden1991]; for the higher-order poles, up to order 12 in $e_8^{(1)}$, contributing diagrams have been identified, but the full residues have not been checked [@corrigan1993].

**`sec-one-loop-masses-real`: correct.** The bubble computation gives δm_a²/m_a² = −(β²/4h)cot(π/h) for every species of a₂–a₅, d₄ and d₅ (0 for sinh-Gordon).
- Add the formula at 12:69–71: $$\frac{\delta m_a^2}{m_a^2}=-\frac{\beta^2}{4h}\cot\frac\pi h\qquad(\hat g\text{ simply laced}),$$ independent of $a$ (normal-ordered scheme).
- 12:35: normal ordering is absorbed into m² "(together with a constant shift of $\vec\phi$)", except for a_n, where Z_h symmetry makes the shift unnecessary.

**Other errors in ch. 12–13:**
- **12:55, spurious i.** Write "$s-m_c^2\simeq2m_am_b\sinh(iu_{ab}^c)(\theta-iu_{ab}^c)=2i\,m_am_b\sin u_{ab}^c\,(\theta-iu_{ab}^c)$". The residue iβ²/2h that follows is correct.
- **12:45, kinematics.** For distinct masses reflection is kinematically forbidden. Replace with: "For $m_a=m_b$ the reflection $a+b\to b+a$ is also kinematically allowed. For identical particles it coincides with transmission, and otherwise its amplitude vanishes in affine Toda theory, whose S-matrix is diagonal."
- **13:53.** The pole statement fails for a+b > h; for a₄ S₃₄ there are simple poles at iπ/5 and 3iπ/5 only. "In general, for $a\ne b$," → "For $a\ne b$ and $a+b<h$,", and add "for $a+b>h$ use $S_{ab}=S_{h-a,h-b}$, so that the s-channel pole is at $i\pi(2h-a-b)/h$".

**Ch. 15:**
- eq-tba, c_eff = r and the qualitative form-factor axioms are correct.
- 15:44: `fring1993` treats sinh-Gordon only, so "was carried out" → "was begun".
- 15:25: cite AFKRY alongside `ahn2000` (item 1).
- `exr-tba-uv` needs a hint: "for the affine Toda S-matrices $N_{ab}=\frac1{2\pi}\int\varphi_{ab}=\delta_{ab}$, so the plateau equations $x_a=1+x_a$ have no finite solution; the limit $x_a\to\infty$ gives $c_{\rm eff}=r$."

**Unclear:**
- 13:106 and 15:56: the Albeverio–Høegh-Krohn range. AHK 1974 covers α² < 4π in their normalization; check against J. Funct. Anal. 16 (1974) 39.
- 14:61: whether the breather identification fixes B(β) exactly. This depends on Part VI.


## Item 17: Part V statements

**Hollowood: the book's ¼cot(π/h) is correct, but the attribution needs a note.**
- Hollowood's paper (hep-th/9209024) prints ½cot. His bound-state sum (3.2)→(3.3) and his continuum integral (3.26)→(3.27) are both evaluated at twice their true values (numerical evaluation of his own printed expressions):

  | | Hollowood prints | true value |
  |---|---|---|
  | bound-state sum, h even | m_a cot(π/h) | ½ m_a cot(π/h) |
  | bound-state sum, h odd | ½ m_a(cot + cosec) | ¼ m_a(cot + cosec) |
  | continuum integral, h even | −½cot | −¼cot |
  | continuum integral, h odd | −½cosec | −¼cosec |

- His particle shift (1.4) is doubled too, which is why his λ̃ = −1 check still worked.
- MacKay–Watts (hep-th/9411169) have ¼.
- m_a = 2m sin(πa/h), with m the normal-ordered mass parameter. In the same scheme m̂_a = m_a(1 + (β²/8h)cot(π/h)), so M_a/m̂_a = 2h/β² − h/2π, independent of a: the cot terms cancel.

Replacements:
- 16:84: "summing over $b$ gives Hollowood's result" → "summing over $b$ gives".
- After 16:88 insert:
  > Here $m_a=2m\sin\frac{\pi a}h$, with $m$ the mass parameter of the normal-ordered Lagrangian (@sec-feynman-rules). The calculation is Hollowood's [@hollowood1993a]; his paper prints $\frac12\cot\frac\pi h$, because its bound-state sum and continuum integral (his (3.2) and (3.26)) are evaluated as twice their actual values. Evaluated correctly they give $\frac14\cot\frac\pi h$, as found by MacKay and Watts [@mackay_watts1995].
- After 16:90 add:
  > In the same scheme the one-loop particle masses are $\hat m_a=m_a\big(1+\frac{\beta^2}{8h}\cot\frac\pi h\big)$, the continuation of $\delta m_a^2/m_a^2=-\frac{\beta_r^2}{4h}\cot\frac\pi h$ (@sec-one-loop-masses-real). So $M_a/\hat m_a=\frac{2h}{\beta^2}-\frac h{2\pi}+O(\beta^2)$, independent of $a$: the $\cot$ terms cancel.
- 19:75: "and Hollowood's $a_n^{(1)}$ masses" → "and the $a_n^{(1)}$ masses (-@eq-hollowood-mass) (Hollowood's calculation, with the factor 2 in his $\cot$ term corrected)".
- 19:18: "reproduces Hollowood's simply-laced $a_{N-1}^{(1)}$ masses" → "reproduces the $a_{N-1}^{(1)}$ masses (-@eq-hollowood-mass)".
- 16:78: the letter B = πb/h clashes with the coupling B; use C or B_b.

**`eq-an-transmission` and κ (16:38): sign error.**
- "κ=√(m²−λ)=−ik" is inconsistent with X(z) = (z−1)/(z+1) being the transmission coefficient T(k) = (k+im)/(k−im). The bibliography agent found this and the orchestrator confirmed it by hand.
- Replace with κ = ik. `exr-pt-transmission` needs the same fix (with κ = −ik its zeros land at negative z).
- With κ = ik, the formula matches Hollowood (2.11), and the a₂ ψ-sum reproduces eq-hollowood-mass.

**`sec-transmission-factors` sketch: correct.**
- For reflectionless potentials δ(0) − δ(∞) = πn_b, with the zero mode included and the threshold state a half-bound state; the Born 1/k term is what normal ordering cancels. Suggested clarification: "(for a reflectionless potential $\delta(0)-\delta(\infty)=\pi n_b$, with $n_b$ counting the zero mode; the threshold solution is a half-bound state)".
- 16:78 states "neither reflect nor change species" without qualification. Append "(for the Hirota solitons used here, @sec-counting-rule-transmission-factors; near special eigenvalues this structure needs care, @thm-jordan-chains)" to align with hypotheses (H)/(R) of ch. 18.

**`sec-picard-lefschetz`: correct.**
- Upward flow for J_σ and downward for K_σ, with n_σ = ⟨ℝ, K_σ⟩.
- In the example: n₀ = 1 and n_i = 0. K_i (y ≥ 1 on the imaginary axis) never meets ℝ.
- Optional remark: real λ lies on a Stokes line, since both actions are real, but n_i = 0 on both sides.

**`sec-dhn-complex`: correct.** The MacKay–Watts versus Delius–Grisaru disagreement is genuine and consistent with ch. 19. The one-loop particle results reproduce the ch. 19 table via its untwisted rule. Ch. 19 should state that "Lie dual" means (g^∨)⁽¹⁾ there.

**`exr-hollowood` (16:119): false note.** For a₂, X₁₁ = (z−1)/(z+½), whose pole at −½ lies outside (0,1]. The note is true for a₄. Replace the last two sentences with:
> Note that $X_{11}$ has its pole at $z=-\frac12$, outside the bound-state half-plane. In $a_4^{(1)}$, by contrast, $X_{11}$ has a pole at $z=\cos\frac{2\pi}5\in(0,1]$. What does that mean physically? (Compare @prp-an-zeros and @exr-a4-cancellation.)

**Other exercises of ch. 16–19: true.** These are `exr-ccg-psi`, `exr-phi4`, `exr-mass-ratio-oneloop`, `exr-airy`, `exr-iphi3`, `exr-hg` (geometric multiplicity 1, algebraic 2), `exr-wronskian-bilinear`, `exr-chain`, `exr-psi-properties`, `exr-sg-calibration` and `exr-g2-duality`.

## Item 19 (part 1): chapter 20, `sec-nonlocal-charges` and `sec-qg-coproduct`

Scripts: `referee_textbook_scripts/qg_cross/spins.py` (exact), `sq_antipode.py`, `scalar_ex.py`. Checked against Bernard–LeClair 1991, CMP 142, 99. The text was extracted, but every conclusion rests on the computations.

- **Currents J_± = e^{±iaφ_L} with a = 8π/β_SG and h_J = 8π/β_SG²: correct.** a is fixed by the double pole with e^{∓iβφ}. H_± ∝ m₀² × a mixed-chirality vertex operator, and the spins ±(h_J−1) = ±λ.
- **"carries topological charge ±1" (20:16): wrong, it is ±2.** The current shifts φ by 4π/β_SG, i.e. two soliton units. This is BL (3.20f), and it is forced by k₁ = q^T = diag(q, q⁻¹). Replace with:
  > Bernard and LeClair choose the exponent $a=8\pi/\beta_{\rm SG}$, for which $J_\pm$ has a double pole with $e^{\mp i\beta_{\rm SG}\phi}$, so that it remains conserved when the perturbation is switched on. The current carries topological charge $\pm2$ (it shifts $\phi$ by $\pm4\pi/\beta_{\rm SG}$): $[T,Q_\pm]=\pm2Q_\pm$ with $T=\pm1$ on the soliton and antisoliton:
- **Construction "along the coroots" and the simply-laced spin 4π/β² − 1: correct.** The current of node j is exp(−i(4π/β)α_j^∨·φ_L), of weight 8π/(β²α_j²). The braiding is q^{α_i^∨·α_j^∨}.
- **Missing: the non-simply-laced spins.** They are s_j = 8π/(β²α_j²) − 1:
  - ω on long roots;
  - 2ω + 1 on roots of length² 1;
  - 3ω + 2 on roots of length² 2/3.

  Add at 20:42: "In general the current of node $j$ is $\exp(-i\frac{4\pi}\beta\vec\alpha_j^\vee\cdot\vec\phi_L)$, of weight $8\pi/(\beta^2\vec\alpha_j^2)$, so $s_j=\frac{8\pi}{\beta^2\vec\alpha_j^2}-1$."
- **20:64, 20:70, the null-root spin: wrong for non-simply-laced ĝ.**
  - The invariant is Σ_j n_j^∨ s_j, the spin of the null root of ĝ^∨, with n₀^∨ = 1.
  - `spins.py` reproduces every Part VI T exactly: c₂–c₄, a₅⁽²⁾, a₇⁽²⁾, e₆⁽²⁾ = 9ω+3, f₄⁽¹⁾ = 6ω+3/2, g₂⁽¹⁾ = 3ω+1, d₄⁽³⁾ = 6ω+3, sine-Gordon = ω.
  - The literal "½Σn_j s_j" fails for every non-simply-laced case: c₃ gives 5ω+2 instead of 3ω+1.
  - 4πh/β² − h^∨ in book normalization fails for every twisted case: d₄⁽³⁾ gives 2ω+1.

  Replacements:
  - 20:64: "…so only the total spin $\sum_jn_j^\vee s_j$ of the null root $\sum_jn_j^\vee\vec\alpha_j^\vee$ of $\hat g^\vee$ (normalized by $n_0^\vee=1$) is fixed, and it determines $T$."
  - 20:70: "where $2T=\sum_jn_j^\vee s_j$ [@delius1995; @gmw1996]: for untwisted $\hat g$, $n_j^\vee=n_j\vec\alpha_j^2/2$ and the sum is exactly $4\pi h/\beta^2-h^\vee$; for twisted $\hat g=X^{(k)}$, $n_j^\vee=kn_j\vec\alpha_j^2/2$ and $\beta^2$ is replaced by $\beta^2/k$ (@sec-conventions), with $h^\vee=\sum_jn_j^\vee$."
  - The same caveat is needed at 21:25–26.
- **`sec-qg-coproduct`:**
  - The form Δ(Q_j) = Q_j⊗1 + q^{T_j}⊗Q_j is correct (BL (3.37)), with T_j along the coroot.
  - **q = ±e^{iπh_J} (20:55) is wrong for this coproduct.** It is q = e^{−iπh_J} = e^{−8π²i/β_SG²} = e^{−4π²i/β²}. The sign is not free: only this value reproduces Zamolodchikov's S-matrix and crossing with +1 (item 2). Replace with:
    > …and $q$ is the braiding phase of the chiral vertex operators. With this coproduct, Bernard and LeClair find [@bernard1991] $q=e^{-i\pi h_J}=e^{-8\pi^2i/\beta_{\rm SG}^2}=e^{-4\pi^2i/\beta^2}$ (for sine-Gordon, $q=-e^{-i\pi\lambda}$), with $q_j=q^{\vec\alpha_j^{\vee2}/2}$. With the opposite coproduct (-@eq-coproduct) the same charges give $q^{-1}$, the value of @sec-six-vertex.
  - 20:58: the identification is correct as an algebra map, but f_j ∝ Q̄_j q^{−T_j}, not Q̄_j. Replace "(up to normalization)" with "(up to a constant and the factor $q^{-T_j}$), for the coproduct $\Delta(e_j)=e_j\otimes1+k_j\otimes e_j$, $\Delta(f_j)=f_j\otimes k_j^{-1}+1\otimes f_j$, the opposite of (-@eq-coproduct)".
  - 20:74, the crossing bullet, mixes the conventions. Replace with: "With the coproduct of @sec-qg-coproduct the sine-Gordon crossing point is $x_c=q^{-2}$ ($q^2$ for that of @sec-hopf), and $e^{2\pi i\lambda}=q^{-2}$ holds exactly for $q=-e^{-i\pi\lambda}$."
- **Exercises of ch. 20 and 21:**
  - `exr-nonlocal-spin` (20:97): true but imprecise ("up to O(m₀²)" is meaningless, since H is O(m₀²)). Suggested: "…check that $H_\pm=m_0^2V_\pm$ with $V_\pm$ of weights $(h_J-2+\beta_{\rm SG}^2/8\pi,\ \beta_{\rm SG}^2/8\pi)$, so that, counting the weights of $m_0^2$, $H_\pm$ has weights $(h_J-1,1)$ as $\partial_-J_\pm=\partial_+H_\pm$ requires."
  - `exr-coproduct-q`: correct, and it exposes the coproduct point above.
  - `exr-sg-gradation`, `exr-an-gradation`: correct.
  - **`exr-scalar-unitarity` (21:111): partly wrong.**
    - f(t)f(−t) = 1/(c(t)c(−t)) needs the reflection formula Γ(z)Γ(1−z) = π/sin πz, not the functional equation.
    - **The printed infinite product diverges**: its k-th factor behaves as (2Tk)^{−4NT} (log factor −80.23 at k = 1000 against −80.24 predicted). It converges once each factor is divided by a t-independent constant, as `spinorn.py` does.
    - Fix the exercise text to: "Show that $f$ satisfies $f(T-t)=f(t)$ and, using $\Gamma(z)\Gamma(1-z)=\pi/\sin\pi z$, that $f(t)f(-t)=1/\big(c(t)c(-t)\big)$". In `sec-scalar-factor`, add "(each factor of the product normalized by a $t$-independent constant, which makes it converge)".
  - `exr-lift-cdd`, `exr-an-fusing`, `exr-block-translation`: correct.

## Item 19 (part 2): Part VI framing, the worked example, chapter 27

Scripts: `referee_textbook_scripts/part67/an_cross.py`, `an_rmat.py`, `dgz_compare.py`, `wmin.py`. The HHM paper text is at `part67/hhm.txt`. The orchestrator re-ran `dgz_compare.py` and checked its encoding of ch. 23 against the printed formula.

### `sec-an-example` (21:74–106): correct

- **"ρ is the six-vertex eigenvalue for every n": correct.** In the homogeneous variable x = e^{2Tθ}, 2T = hω, the U_q(sl_{n+1}^) vector eigenvalue does not depend on n. For n = 2, Jimbo's equations give 1 (×6) and ρ (×3).
- **Pole at x = q⁻², zero at q², fusion θ = 2πi/h, M_a ∝ sin(πa/h), n = 1 at θ = iπ: correct.** These depend only on q², so the q-sign finding of item 2 does not affect them.
- **Residue: correct.** It projects onto Λ²: the q-antisymmetric subspace is a submodule exactly at x₁/x₂ = q⁻² (Part VI coproduct), for N = 3, 4, 5.
- **Crossing (21:102): correct, and stronger than stated.**
  - V₁(x)\* is isomorphic to the fused antisoliton multiplet at x·x(iπ) for **either** sign of q; the spectral offsets of the fusion string absorb the sign.
  - So, unlike the self-conjugate c_n spinor, a_n gives no handle on the sign of q.
- 21:86: replace the parenthetical with "(the R-matrix of this chapter is that of @sec-six-vertex with $q\to q^{-1}$; its pole and zero depend only on $q^2$)". This is consistent with item 2.
- 21:104 (Gandenberger): see item 1.

### `sec-qg-restriction` (27:17): correct, with one omission

- q² is a primitive p-th root of unity, so the alcove is at level p − n − 1.
- n = 1 reproduces `prp-sg-rsos`: p − 1 vacua, with heights on the smaller label p.
- W₃(4,5) at level 1 gives 3 vacua (Potts).
- Add the condition p ≥ n+1: "p<p' coprime and $p\ge n+1$, q is a root of unity".

### `sec-perturbed-minimal` (27:24–37)

**a_n⁽¹⁾ row: correct.**
- c = n[1 − (n+1)(n+2)(p′−p)²/(pp′)].
- The perturbing field is Φ(0|θ), with h = (n+1)p/p′ − n. This was checked for n = 1..4 against the dimension of e^{iβα₀·φ}.
- Examples: W₃(4,5) gives c = 4/5, h = 2/5; M(2,5) gives c = −22/5, h = −1/5.
- Suggest saying which side carries the adjoint: "W_{n+1}(p,p′) perturbed by Φ(0|θ), trivial on the level-(p−n−1) side that labels the vacua, h = (n+1)p/p′ − n".

**"Other algebras" row (27:33): wrong citations.**
- takacs1997a (NPB 489, 532) and tkw1997 (NPB 489, 557) are a₂⁽²⁾ (Zhiber–Mikhailov–Shabat) φ₁,₅ papers.
- g₂⁽¹⁾ and d₄⁽³⁾ give W₃ models perturbed by Φ(11|12) and Φ(11|14) (takacs1997b, takacs1997c).

Replacement rows:
- a₂⁽²⁾: "M(p,p′)+φ_{1,2} or φ_{2,1} [@smirnov1991; @efthimiou1993]; a second restriction gives φ_{1,5} of non-unitary M(3,p′) [@takacs1997a; @tkw1997]"
- other algebras: "g₂⁽¹⁾, d₄⁽³⁾: W₃ minimal models perturbed by Φ(11|12), Φ(11|14) [@takacs1997b; @takacs1997c]"

Bibliography: `takacs1997b` (NPB 501 (1997) 711, no title) and `takacs1997d` (hep-th/9702196) are the same paper. Merge them.

Observation: every unitary W_{n+1}(p,p+1) restriction lies above Part VII's β²_UV = 4πn/(n+1). The unitary subtheories of ch. 27 are therefore outside the continuum theorems of ch. 30–31, and one sentence saying so would help.

### `sec-hermitian-analyticity-irf`: the HHM mechanism is misdescribed (six places)

Hoare–Hollowood–Miramontes (arXiv:1303.1447):
- They work at **q = e^{iπ/k}** (positive quantum dimensions), not at a general root of unity.
- Their examples are:
  - the q-deformed AdS₅×S⁵ world-sheet S-matrix, U_q(psu(2|2)⋉ℝ³), with the IRF transformation acting on the bosonic su(2)⊕su(2);
  - its generalized sine-Gordon limit;
  - restricted sine-Gordon as a warm-up.
- The examples are not "U_q(su(2)⁽¹⁾)-based S-matrices".
- The book's paraphrase of *why* the mechanism works is fair.

Replacements:
- 27:41: "Hoare, Hollowood and Miramontes [-@hoare2013] showed this at $q=e^{i\pi/k}$, where all quantum dimensions are positive, for the $q$-deformed AdS$_5\times$S$^5$ world-sheet S-matrix (symmetry $U_q(psu(2|2)\ltimes\mathbb R^3)$, with the IRF transformation acting on its bosonic $U_q(su(2))\oplus U_q(su(2))$) and its relativistic generalized sine-Gordon limit; the restricted sine-Gordon model is their simplest example."
- 27:57: "…the $q$-deformed world-sheet and generalized sine-Gordon S-matrices of Hoare, Hollowood and Miramontes [-@hoare2013], whose IRF form at $q=e^{i\pi/k}$ is manifestly Hermitian analytic".
- 22:47: "for $q=e^{i\pi/k}$, the IRF form of the S-matrices they studied (the $q$-deformed world-sheet S-matrix, symmetry $U_q(psu(2|2)\ltimes\mathbb R^3)$, and its generalized sine-Gordon limit)".
- 33:116: "mechanism shown for the $q$-deformed world-sheet and generalized sine-Gordon S-matrices at $q=e^{i\pi/k}$".
- 33:172: "beyond the su(2)-type cases of [@hoare2013]".
- 27:15: "unitarity" → "braiding unitarity".

### Prior art: Delius–Gould–Zhang R-matrices (verified symbolically)

The identities below are proved symbolically by `dgz_compare.py` and re-run by the orchestrator:
- **Ch. 23 c_n⁽¹⁾ spinor R-matrix ≡ DGZ 1996 (4.31) at a = 1**, for n = 2..8, with x_DGZ = x and q_DGZ = q⁻² (equivalently x → 1/x, q_DGZ = q²). The factor 2 is root normalization; the inversion is the opposite coproduct.
- **Ch. 24 a_{2n−1}⁽²⁾ formula ≡ DGZ 1994 (4.21)–(4.24)**, for n = 2..9, with the same dictionary. No x → −x and no sign of q is needed.
- **Ch. 25 table ≡ DGZ 1994 (4.64)**, normalized on the 324.

"DGZ" means Delius–Grisaru–Zanon in Part VI (21:68), so spell out "Delius, Gould and Zhang". The S-matrices, scalar factors, lifts, bound states and breathers remain new. Sentences to change:
1. **08:11** (also item 4): "For many of the representations needed for the non-simply-laced solitons they were found in the 1990s by the tensor-product graph of this chapter: Delius, Gould and Zhang derived the spinor R-matrices of $U_q(b_n^{(1)})$ and $U_q(d_{n+1}^{(2)})$ for all $n$, and the $\mathbf{26}$ R-matrix of $U_q(f_4^{(1)})$ [@dgz1994; @dgz1996]. What remained hard were representations that are reducible under $U_q(g)$, such as the Kirillov–Reshetikhin modules of higher solitons and the $\mathbf{27}$ of $U_q(e_6^{(2)})$ (@sec-sm-introduction)."
2. **08:73**: factually wrong; replacement in item 4.
3. **22:32**: "The spinor R-matrices are known in closed form for all $n$ [@dgz1994; @dgz1996]. What is conjectured beyond $n=4$ is the S-matrix built on them: the scalar factor with its integer lifts and the breather identification. For representations whose tensor products are not multiplicity-free (the fused $c_3^{(1)}$ solitons, the $\mathbf{27}$ of $U_q(e_6^{(2)})$), Part VI solves the intertwining equations directly."
4. **22:89–92**: "The R-matrices of the smallest multiplets were available: Delius, Gould and Zhang had derived the spinor R-matrices of $U_q(b_n^{(1)})$ and $U_q(d_{n+1}^{(2)})$ and the $\mathbf{26}$ R-matrix of $U_q(f_4^{(1)})$ [@dgz1994; @dgz1996]. Two obstacles remained. – **Multiplicities.** The tensor-product graph needs multiplicity-free tensor products, which fails for the higher (Kirillov–Reshetikhin) soliton multiplets, e.g. of $c_3^{(1)}$, and for the $\mathbf{27}$ of $U_q(e_6^{(2)})$." 22:91, "fail this for c_n⁽¹⁾ with n ≥ 3", is wrong for the spinor.
5. **22:98**: "built on the $U_q(d_{n+1}^{(2)})$ spinor R-matrix of Delius, Gould and Zhang [@dgz1996, eq. (4.31)], and verified for $n=2,3,4$". This line says n = 2, 3, 4 while 22:27 says n = 3, 4; make them consistent.
6. **22:99**: "from the $U_q(b_n^{(1)})$ spinor R-matrix of Delius, Gould and Zhang [@dgz1994, eqs. (4.21)–(4.25)], verified for $n=3,4$".
7. **22:100**: append "on the $\mathbf{26}$ R-matrix of [@dgz1994, eq. (4.64)]".
8. **23:26**: "This is the R-matrix derived for all $n$ by Delius, Gould and Zhang with the twisted tensor-product graph [@dgz1996, eqs. (4.30)–(4.31) at $a=1$, with their $q$ equal to our $q^{-2}$]; we re-derived it from direct solutions of the intertwining equations for $n=2,3,4$, to $10^{-13}$."
9. **23:141**: "The spinor R-matrix is that of Delius, Gould and Zhang [@dgz1996]; the S-matrix construction for $n\ge3$ is new here."
10. **24:25**: "This is the $U_q(b_n^{(1)})$ spinor R-matrix derived for all $n$ by Delius, Gould and Zhang [@dgz1994, eqs. (4.21)–(4.25), with their $q$ equal to our $q^{-2}$]; we confirmed it at $n=3,4$ against direct solutions, to $10^{-12}$."
11. **24:77**: "The R-matrix is that of [@dgz1994]; the S-matrix construction for $n\ge3$ is new here."
12. **25:67**: "The $\mathbf{26}\otimes\mathbf{26}$ R-matrix was obtained by Delius, Gould and Zhang with the tensor-product graph [@dgz1994, eq. (4.64)]; normalised to 1 on the $\mathbf{324}$ it is the table above. Here it is re-derived directly from the defining relations (@sec-app-numerics)."
13. **26:162** (open problem 5): "The spinor R-matrices are known for all $n$ [@dgz1994; @dgz1996]. The S-matrices built on them are conjectured for $n\ge5$; checking them needs spinor projectors of dimension $2^{2n}$, at increasing cost."
14. **33:159** (open problem 2): "(the S-matrices are verified for $n\le4$; their R-matrices are known for all $n$ [@dgz1994; @dgz1996])".

**Unclear:** 22:28 says the 27 R-matrix of U_q(e₆⁽²⁾) "is known [@gmw1996]", while 22:84 and 22:102 say "we construct" it. Reading GMW 1996 (hep-th/9509007) would settle this.

### Open consequence of the q-sign finding

Ch. 32's multi-soliton sector results (Krein/Θ-broken sectors, open problem 6) were presumably computed with q = e^{−iπω}. Ch. 27:81 itself notes that p′ = p+2 unbrokenness "depends on the sign of q", and a q → −q twist acts like the Z₂ grading discussed in `briefs/exchange_statistics.md`.

**Unclear**, and important. Settle it by rerunning `foundations_code/sectors/` (three-soliton c_n and a_{2n−1}⁽²⁾ sectors) with q → −q. This links item 2 directly to the open exchange-statistics brief.

## Item 20 (rest): Part VII framing and exercises

### Learning goals, introductions and notes

1. **30:6 is wrong.** `sec-threshold-sharp` proves sharpness of β²_UV, not of 2π, and 2π is not sharp for rank ≥ 2. Replace with: "Follow the proof that the continuum limit exists in finite volume for β²<2π, and why no continuum limit exists at or above the first collapse threshold β²_UV without further counterterms."
2. **30:10 is wrong.** β²_UV exceeds 2π except for a₁⁽¹⁾. Replace with: "proves the continuum limit for β²<2π, where a nearest-neighbour bound suffices; for sine-Gordon this is already the whole range below the first collapse threshold."
3. **31:10.** "no neutral cluster collapses" → "no cluster collapses". For b_n, c_n, g₂ and a₂⁽²⁾ the minimizing cluster is a non-neutral pair.
4. **31:6 is overstated.** The remark at 31:331 is [heuristic] and conditional on K(β) ≠ 0. "why its exponent is sharp" → "why its exponent is expected to be optimal".
5. **Status tension.** The chapter 31 notes say "have not yet been independently refereed", while the 28:74 table says H3 is "closed in finite volume" and 33:156 says open problem 1 is "now settled below β²_UV". Use "proved here (awaiting independent check)" in the table and in open problem 1.
6. **33:191.** "The open problems are also the subject of the task briefs in the repository's `briefs/` directory" is internal workflow. Delete it.

Everything else checked out:
- the learning goals and introductions of ch. 28, 29 and 32;
- "four structural theorems";
- the numbers in open problem 8;
- all cross-references in ch. 20–33 resolve.

### Exercises of ch. 28–33 (14)

- **Correct:**
  - `exr-gap-map`, `exr-lee-yang-criteria`, `exr-open-problem`, `exr-u5-lee-yang`;
  - `exr-single-site`;
  - `exr-trace-positivity`;
  - `exr-zero-mode-sign` (sufficiency: α₀·s ≡ π iff Σn_j is even);
  - `exr-pair-collapse` (∫r^{1−q²/2π}dr converges iff q² < 4π, with q² = 2β² = β_SG²);
  - `exr-window-nonempty`;
  - `exr-krein-eigs`.
- **`exr-thresholds`:** solvable, with answers 2π, 8π/3, 4π, 4π·29/30. For c₂⁽¹⁾, two different clusters both give 4π, so "the cluster" is not unique.
- **`exr-krein-transfer`:** correct, but the proof also needs 𝒦 to be real symmetric. Write "that the Gaussian kernel is real, symmetric and satisfies $\mathcal K(\phi,\phi')=\mathcal K(-\phi,-\phi')$".
- **`exr-stability` (31:344): imprecise.** Dimensional analysis alone fixes no exponent in finite volume; the exponent comes from extensivity. The bare z has dimension 2; the physical g = zR_Ω^x has dimension 2 − x. Replacement:
  > The physical coupling $g=zR_\Omega^{x}$ of @sec-wick-ordering-statement has mass dimension $2-x$, with $x=2\kappa$ for a long root. Assuming that for large $z$ the free energy is extensive, $\log\Xi\approx K\int_\Omega g(u)^{2/(2-x)}d^2u$, show that $\log\Xi$ grows like $z^{1/(1-\kappa)}$, so that the bound of @lem-stability with $\theta=\kappa$ has the largest exponent compatible with extensivity.
- **`exr-bethe-complex` (32:192): gap.** θ = iv gives the real energy m cos v, so "and hence the energy" does not follow. Replacement: "…show that $|\Lambda|\ne1$ on the real axis forces $\theta$ off the real axis, and that the energy $m\cosh\theta$ is then complex unless $\operatorname{Re}\theta=0$ or $\operatorname{Im}\theta\in\pi\mathbb Z$."

### Exercises of ch. 22–27 (14): all true and solvable

- **`exr-g2-d43-T`:** solvable only with Kac's h = 4 for d₄⁽³⁾ (2T = 12ω + 6). Add a hint pointing to appendix A.
- **Values checked:**
  - `exr-g2-split`: H_s − 6 = −(H_b − 6) at first order;
  - `exr-cn-lifts`: T − a_j = (n−j)(ω+½);
  - `exr-vjs-expand`: agreement once β²_ch24 = β²/2;
  - `exr-cds-masses`: 2cos(π/18), 2cos(4π/18) at H = 18;
  - `exr-e6-H`: H = 18 + 3β²/2π + O(β⁴);
  - `exr-lee-yang-restriction`: (2,5), one breather, no kinks;
  - `exr-quantum-dims`: p = 3, 4 checked by hand.

## Item 12: exercises in chapters 2–33

Every exercise was worked at least in outline, and every number or formula was checked numerically. Only exercises that are false, misleading or need a hint are listed individually; all others are **true and solvable at the chapter's level**.

| Ch. | Count | Not OK |
|---|---|---|
| 2 | 6 | none |
| 3 | 7 | **`exr-sg-lee-yang`**: the last part asks for a self-coupling pole that cannot exist in unrestricted sine-Gordon, because B₁ is C-odd (replacement in item 6) |
| 4 | 5 | **`exr-krein-collision`: false** (2×2 counterexample; replacement in item 8) |
| 5 | 5 | none |
| 6 | 5 | none |
| 7 | 6 | 07:210: "$m\circ(S\otimes1)\circ\Delta=\epsilon$" should be "$=\eta\circ\epsilon$". 07:218, `exr-six-vertex`: "whose residue is the projector onto spin 0?" → "…proportional to the projector onto spin 0?" (the residue is (1−q⁴)P₀) |
| 8 | 5 | `exr-faces` (08:137): the hint also needs the CG coefficients of $1\otimes\frac12$. `exr-spin1` is solvable but heavy: add "*Hint:* $V_1(x)$ is the $q$-symmetric part of $V_{1/2}(xq)\otimes V_{1/2}(xq^{-1})$, the image of $\check R(q^{-2})$; fuse first in the left and then in the right factor." |
| 9 | 6 | none |
| 10 | 7 | `exr-moving-soliton`: wording should be "½(D_x²−D_t²)τ_j·τ_j = m_a²(cosh²θ−sinh²θ)ω^{ja}E". `exr-soliton-mass`: suggested hint E = 2∫V = (2/β²)Σ[∂ ln τ_j] |
| 11 | 4 | **`exr-breather-real`: misleading** (the field is singular on curves; item 13). `exr-theta-solutions`: needs an a₁⁽¹⁾ exception |
| 12 | 4 | none |
| 13 | 5 | none (`exr-an-poles` is fine, but see the 13:53 fix in item 15) |
| 14 | 3 | **`exr-cn-snn`: misleading**, since the poles are not bound states (item 15). `exr-floating-endpoints`: rewording after the eq-floating-h fix |
| 15 | 3 | `exr-tba-uv` needs a hint (item 15) |
| 16 | 5 | **`exr-hollowood`: the note is false for a₂** (item 17). `exr-pt-transmission` needs κ = ik |
| 17 | 4 | none |
| 18 | 4 | none |
| 19 | 4 | none |
| 20 | 4 | `exr-nonlocal-spin`: imprecise wording (item 19) |
| 21 | 4 | **`exr-scalar-unitarity`: the printed product diverges**, and the exercise names the wrong Gamma identity (item 19) |
| 22 | 3 | `exr-g2-d43-T` needs a hint: Kac h = 4, h^∨ = 6 for d₄⁽³⁾ |
| 23–27 | 11 | none |
| 28 | 2 | none |
| 29 | 4 | `exr-krein-transfer`: add "real, symmetric" (item 20) |
| 30 | 2 | `exr-thresholds`: for c₂⁽¹⁾ "the cluster" is not unique |
| 31 | 2 | **`exr-stability`: imprecise**, since the exponent comes from extensivity, not dimensional analysis (item 20) |
| 32 | 2 | **`exr-bethe-complex`: gap**, since the energy can be real off the real axis (item 20) |
| 33 | 2 | none |

**Count.** Of 122 exercises, 6 are false or misleading as stated: `exr-krein-collision`, `exr-sg-lee-yang`, `exr-breather-real`, `exr-cn-snn`, `exr-hollowood` and `exr-scalar-unitarity`. Three more are imprecise or have a gap (`exr-stability`, `exr-bethe-complex`, `exr-nonlocal-spin`), and about ten need wording or hint fixes.

## Scripts and outputs

All scripts are in `referee_textbook_scripts/` next to this report. Each needs the copies of `textbook_code/*.py` or `foundations_code/hirota/*.py` that sit in the same folder. The full source and output of the essential ones are reproduced in the agents' reports. The TPG agent's report is kept verbatim as `(the TPG agent report, not kept in the repository)`.

**`referee_textbook_scripts/part1/`**

| script | checks |
|---|---|
| `sg_s0.py` | S₀ integral against a Gamma-product continuation; crossing at λ = 0.4, 0.7; strip divergence |
| `sg_bootstrap.py` | YBE; breather residue rank 1; S_{sB₁}, S₁₁ = f_{ξ/π} by fusion |
| `sg_residues.py` | breather residues |
| `bootstrap_conv.py` | eq-bootstrap on a₂, a₄, e₆, e₈ |
| `misc_checks.py` | CDD counterexamples; Lee–Yang residue; Mostafazadeh and Krein counterexamples |
| `charge_sign.py` | conserved-charge sign |

**`referee_textbook_scripts/qg_cross/`**

| script | checks |
|---|---|
| `qdict.py` | ch. 7 / Part VI dictionary (symbolic) |
| `qtype.py` | type-(−1) representation |
| `qsign_cn.py`, `qsign_bn.py`, `qsign_n4.py`, `qsign_L.py` | crossing sign at ±q |
| `sq_antipode.py` | S² shifts for sl₂, sl₃, c₂⁽¹⁾, d₃⁽²⁾ |
| `spins.py` | non-simply-laced spins and T |
| `scalar_ex.py` | divergence of the ch. 21 product |
| `sg_scalar.py` | F = −S₀ |
| `an_vec.py` | sl₃ vector with the Part VI coproduct |

**`referee_textbook_scripts/qg_tpg/`**

| script | checks |
|---|---|
| `qaff.py` | generic U_q(ĝ) relations and Jimbo solver |
| `tpg_checks.py` | a₂, c₂, b₂, c₃, b₃ vector representations |
| `spinor_twisted.py` | d_{n+1}⁽²⁾ spinor, n = 2–4: multiplicity-free, twisted rule |
| `ch7_checks.py` | central element; quasitriangularity; spin-1 form |

**`referee_textbook_scripts/lie_classical/`**

| script | checks |
|---|---|
| `lie.py` | root-system library |
| `ch5_check.py` | tables and Steinberg, all types |
| `ch5_charges.py`, `ch5_charges_b.py` | charges against Cartan eigenvectors |
| `fold.py` | all foldings |
| `twisted_exponents.py` | period kh |
| `duality_masses.py` | Lie duality of soliton masses |
| `centre.py` | diagram automorphisms and centre |
| `e8_ising.py` | E₈/Ising sentence |

**`referee_textbook_scripts/part3/`**

| script | checks |
|---|---|
| `s1_masses_dorey.py`, `s1b_sensitivity.py` | Dorey and area rule in the charge basis versus a real basis |
| `s2_folded_hirota.py` | folded Hirota |
| `s3_soliton_mass.py`, `s3b_e6_species4.py` | soliton masses d₄, d₅, e₆ |
| `s4_mcghee.py` | McGhee count |
| `s5_misc.py` | real a₃ soliton; twisted exponents; time delay; breather sign change |
| `s6_three_soliton.py` | three-soliton solution |
| `s7_twisted_duality.py` | twisted Lie duality |
| `s8_spin2_charge.py`, `s8b_breather_charges.py` | spin-2 charge imaginary on breathers |
| `s9_exercises.py` | ch. 9–11 exercises |

**`referee_textbook_scripts/part45/`**

| script | checks |
|---|---|
| `hollowood_check.py` | one-loop δm²/m², and Hollowood's bound-state sum |
| `hollowood_continuum.py` | Hollowood's continuum integral |
| `floating_check.py` | one-loop H for b_n, c_n, g₂, and the twisted endpoints |
| `f4_check.py` | f₄ H |
| `poles_check.py` | c₂ S₂₂ simple poles and residue sign; a₄ S₃₄ |

**`referee_textbook_scripts/part67/`**

| script | checks |
|---|---|
| `an_cross.py`, `an_rmat.py` | a_n worked example: fusion, residue, crossing |
| `dgz_compare.py` | ch. 23, 24 against Delius–Gould–Zhang, symbolic |
| `wmin.py` | W-minimal model central charges and weights |

**`referee_textbook_scripts/spot/`** holds the orchestrator's spot-checks: the repository's sign scripts at ±q, and the texts of the DGZ papers.

**Baseline.** `textbook_code/` re-run on 2026-10-10. Distinct PASS lines per script: sine_gordon 2, sixvertex 8, universal_r 5, toda_classical 11, smatrix_simply_laced 11, semiclassics 7, rsos_qdims 2, algebra_data 6, spin1_tpg 2. No FAIL.

**Suggested additions to `textbook_code/`:**
- the Gamma-product S₀ (`sg_s0.py`);
- `bootstrap_conv.py`;
- `qdict.py` together with a crossing-sign check;
- `tpg_checks.py` and `spinor_twisted.py`;
- `fold.py` and `twisted_exponents.py`;
- `s1_masses_dorey.py`, `s2_folded_hirota.py`, `s4_mcghee.py` and `s8b_breather_charges.py`;
- `floating_check.py` and `hollowood_check.py`;
- `dgz_compare.py`.
