# Stage 3, Part IX: script citations checkpoint

- ch44: no numerical claim of its own (the abstract and table summarise statements cited in their chapters); figure script ch44_f33_gap_map.py rerun, 6 PASS.
- ch45: zero_mode.py rerun (22 PASS); table matches; cited at the vacuum-gap claim and the opposite-sign paragraph; runtime corrected.
- ch46: uv_thresholds.py rerun (72 PASS); bulk_poles.py rewritten with 24 PASS/FAIL lines (simple poles, none other below k=4, none at 2pi, G<0); figure script ch46 cited for the simply-laced minimum.
- ch47: nested_check.py "708 random trees" corrected to 708 cases (4 shapes, 177 exponent sets, 4 (b,t)); new review614/step3_combinatorics.py (7 PASS); review614 scripts cited individually, logs stored; bulk_poles cited for the sign in the remark.
- ch48: new krein_lightcone.py (15 PASS), sectors/a22_kink_cubic.py (3 PASS), sectors/bloch/t16_graded_onset.py (2 PASS); scan_n4.py rerun (26 PASS); corrections: N<=4 scope of krein_bethe, "0.06 to 4e3" -> "at least 0.05", Bloch-torus "large part" -> 54-89% / 19-51%, 1e-14 -> 1e-13, 1e-9 -> 1e-10, theta_c scope; a22N.py psi label bug fixed and rerun.
- ch49: zero-mode "at all couplings" -> "at both couplings computed", zero_mode.py cited.
- done: check-refs.py no FAIL; render-scratch.sh (HTML) no unresolved references.
