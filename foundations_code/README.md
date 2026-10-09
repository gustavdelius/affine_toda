# Verification code for the foundations paper

Scripts behind the numerical checks quoted in `atft_foundations.qmd`. Python 3 with numpy, scipy, sympy and mpmath. Run each script from its own directory.

| Path | Paper location | What it checks |
|---|---|---|
| `nested_check.py` | Lemma 6.13 | nested-scale integrals with signed exponents and lattice cutoff |
| `review614/` | Theorem 6.14, *Checks* | `lemma613.py`: exact check of Lemma 6.13 on random trees; `bookkeeping.py`: Step 3 exponent bookkeeping on random clusters near a wall; `lattice_delta.py`: the second term of $\Delta_{kl}$ in Lemma 6.11 with exact lattice Green functions |
| `bulk_poles.py` | Section 6.1 | poles of the continued exact bulk energy at the neutral thresholds |
| `vernier/` | Section 3.3, open problem 2 | Vernier–Jacobsen–Saleur $a_{2n-1}^{(2)}$ lattice mass ratios versus the companion-note masses; one-loop particle masses; $\beta^2_{\rm UV}$ |
| `hirota/` | Lemma 8.4, Propositions 8.5–8.8 | Hirota solitons and fluctuations for $d_n^{(1)}$, $e_n^{(1)}$ and foldings; transfer matrices; Wronskians; grid spectra and Jordan types; threshold zeros; heat traces |
| `oneloop/` | Section 8.9 | one-loop masses from transmission factors (Proposition 8.9); calibration; folded composites; comparison with exact results |
| `krein_bethe.py`, `scan_n4.py` | Proposition 9.3 | Krein-unitarity of fundamental-soliton Bethe–Yang transfer matrices; breaking by colour content |
| `sectors/` | Section 9.3, Proposition 9.4 | eigenvalue moduli of soliton, excited-soliton and breather amplitudes of the self-conjugate theories; $a_2^{(2)}$ kinks; exchange-statistics tests |
| `irf/` | Proposition 9.6 | RSOS face weights of $a_2^{(1)}$ kinks; Yang–Baxter, braiding unitarity, Hermitian analyticity, restricted transfer matrices (`checks.py ybe|ha|tm|vertex|signs`) |

Scripts in `sectors/` reuse `rsolve.py` and `dtw.py` from `soliton_rmatrix_code/`, copied alongside.
