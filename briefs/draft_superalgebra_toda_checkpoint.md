# Checkpoint: chapter on affine Toda theories based on Lie superalgebras

Worktree `.claude/worktrees/superalgebra-toda`, branch `claude/superalgebra-toda` (from origin/main b7087eb). Brief: `briefs/draft_superalgebra_toda.md`. One line per finished item.

- 2026-10-23 setup: worktree created; Part III (ch. 9–11) still frozen for the two-phase referee integration (`part3-int-ch9-10`, `part3-int-ch11`, no commits since 2026-10-15), so the ch. 9 edit is deferred and recorded in the report. Part II no longer frozen.
- Sources read: Evans–Madsen 1997 (v2), 1996, 1998; Evans–Hollowood 1991 (scan from the editor; erratum NPB 382 (1992) 662 not obtained), 1992; Delius–Grisaru–Penati–Zanon NPB 359 and PLB 256 (scans from the editor); Penati–Zanon 1992; Gandenberger 1997 §5; Ito–Zhu 2022; Nepomechie 2001; Gomes–Ymai–Zimerman 2006; Aguirre–Gomes–Spano–Zimerman 2015; DGZ 1992 (superalgebra remarks). Source PDFs in the session scratchpad only.
- Scripts done, all passing: `code/super_toda/superspace.py` (9), `components.py` (31), `lax_currents.py` (11), `dgpz_smatrix.py` (27, about 50 s), `n2_solitons.py` (7, about 70 s).
- Findings: Evans–Madsen 1999 (5.44)–(5.46) misprinted (D± swapped, Jbar normalization); DGPZ (B.1) sign of D² misprinted; DGPZ leave the sign of gamma open, tree level fixes gamma = −sin(πΩ); Evans–Hollowood 1992 (5.8) sign of lambda depends on unchecked conventions of (5.5).
- Renumbering: `briefs/superalgebra_toda_renumber.py` run; the index holds the pure-rename commit (77 renames, path fixes in 43 files). All later edits are unstaged. Mapping in `briefs/superalgebra_toda_renames.tsv`. Gap-map SVG (ch. 45) must be regenerated in the text commit (its column labels are chapter numbers).
- Next: chapter text `part4-real-coupling/16-superalgebra-toda.qmd`; then integration edits; conventions document; report.
- Chapter drafted (about 6,100 words); integration edits done in ch. 6, 12, 14, 15, 23, 30, 38, intro, preface, Part IV preamble, outlook, README, App. D (3 entries: Evans–Madsen 1999, DGPZ (B.1), and Gandenberger 1998 (5.10) sign (`gandenberger_scalar.py`). Reading-paths and gap-map figures regenerated (content changes from renumbering).
- Conventions document, drafting report and referee brief written. Remaining: scratch render result and final script log into the report.
- Renders done: HTML (612 s) and PDF (796 s, 669 pages), no unresolved references. Two wide displays split after the PDF check; D(2,1;alpha) mass sentence made precise. Nothing committed; the index holds the pure-rename commit.
- 2026-10-23: on the editor's decision, the two-sentence ch. 9 edit was applied in its own commit while Part III is still frozen; the report says the Part III integration must merge it.
