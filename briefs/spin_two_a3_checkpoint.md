# Checkpoint: reproducible a3 spin-±2 check (defect_spin_two_a3.py)

- [DONE] Read defect_classical_an.py checks 5 and 11 in full; read defect_common.py's an_data/dmatrix.
- [DONE] Prototyped exact (sympy) Gram-Schmidt root-data construction for general r in
  /tmp/.../scratchpad/spin2/proto1.py and proto2.py; verified it reproduces defect_common.an_data(r)/dmatrix(r)
  exactly (floats match) for r = 1, 2, 3.
- [DONE] Found the exact relation: A_abc = 2 sum_k mu_ka mu_kb mu_kc (mu = an_data(r)['U']) together with
  B = 3*(1 - D) solves check 5's EOM-derived matching equation "3 A(u,u,alpha) = -2(alpha.B.u)(alpha.u)"
  EXACTLY (the factor 3 relative to check 11's direct form "A(w,w,alpha) = -2(alpha.E.w)(alpha.w)", which
  uses B = E = 1-D directly). Verified both forms exactly (sympy) for r = 2 and r = 3.
- [DONE] Prototype (proto2.py) reproduces a_2's established check-5/11 results exactly (resid = 0, mism = 0)
  using the general-r machinery, and the bad-D control correctly fails (mism=27, resid != 0) — confirms the
  general-r code is correct before trusting it at r=3.
- [DONE] Ran the EXACT route at r=3 (a_3): all 4 cases (typ=I_-/I_+, eps=+-1=spin -+2), sigma SYMBOLIC:
  mism=0, resid_zero=True for all 4 -- EXACT zero via sympy. Runtime ~112s for the 4 cases.
  Spin +-1 control: all 4 pass exactly, ~4s. Bad-D control (non-integrable D): correctly FAILS
  (mism=79, resid != 0), ~136s. Total prototype runtime ~252s (4.2 min), well under the 15 min budget.
  Logs: /tmp/.../scratchpad/spin2/r3_full.log, r2_control_test.log, control_numeric.log, bulk_cons_test.py output.
- [DONE] Confirmed bulk spin-2 current conserved on shell EXACTLY (symbolic expr = 0, not just small) for a_3
  with this (A, B) choice -- bulk_cons_test.py.
- [DONE] Wrote final production script textbook_code/defect_spin_two_a3.py (exact route throughout, 9 PASS
  lines: root data r=1..4, A_abc/B matching equation exact (a_2 and a_3), cubic obstruction a_3 exact, bulk
  current conserved on shell a_3 exact, main spin+-2 result a_3 exact (both types, sigma symbolic), velocity
  structure (1/sigma for spin+2, sigma for spin-2), spin+-1 control pass, bad-D control FAILS as required).
- [DONE] Checked velocity-structure claim: confirmed for a_3 (every F1 term for eps=-1 carries exactly
  sigma^-1, every term for eps=+1 carries exactly sigma^1) -- matches a_2's pattern, so the qmd text can say
  this for a_3 too.
- [DONE] Ran final script to completion: ALL PASS, total runtime 162.1s (~2.7 min), well under 15 min budget.
  (First run hit a sympy Poly bug on 1/sigma in the velocity-structure check -- fixed by extracting the sigma
  exponent via as_base_exp() instead of sp.Poly(...).degree(); second run: all PASS.)
- [DONE] Edited defect_classical_an.py docstring (check 5 sentence only) to point to new script.
- [DONE] Edited part8-defects/38-defect-conditions.qmd: prp-spin-two label (upgraded a_3 item (3) to Theorem,
  dropped the separate "Numerical for a_3" tag), item (3) text, proof paragraph (3) text.
- [DONE] Edited textbook_code/README.md (new row) and appendices/C-code.qmd (new row, Part VIII table).
- [DONE] Edited briefs/part8_defects_status.md: "Not in the repo" line, "Pending (ch. 38)" line -> DONE, plus
  a new bullet under "By script" describing defect_spin_two_a3.py.
- [DONE] Reran defect_classical_an.py: 24 PASS, 0 FAIL (docstring-only change confirmed harmless).
- [DONE] grep: {#prp-spin-two} still defined exactly once in part8-defects/38-defect-conditions.qmd (the
  referee_scratch/ copy is a separate untouched file belonging to the paused Codex referee); 4 @prp-spin-two
  references across the repo, unaffected since no id was added/removed/renamed.
- [DONE] Cleaned up a stray r3_pid.txt accidentally written to the repo root by a backgrounding quirk (cwd
  resets between Bash calls, so "cd scratch && nohup ... &" followed by "echo ... | tee r3_pid.txt" on the
  next line ran the second command back in the repo root). Removed before committing.
- [DONE] Committed with explicit paths as e748b4f ("Chapter 38: reproducible a3 spin-±2 check
  (defect_spin_two_a3.py)"), 7 files changed; did not touch briefs/referee_defects_upgrades.md or
  referee_scratch/ (other agents' files). TASK COMPLETE.
