# Verification code for the foundations paper

Scripts behind the numerical checks quoted in Part I of the book (`foundations/*.qmd`). Python 3 with numpy, scipy, sympy and mpmath. Run each script from its own directory.

| Path | Paper location | What it checks |
|---|---|---|
| `nested_check.py` | Lemma 6.13 | nested-scale integrals with signed exponents and lattice cutoff |
| `review614/` | Theorem 6.14, *Checks* | `lemma613.py`: exact check of Lemma 6.13 on random trees; `bookkeeping.py`: Step 3 exponent bookkeeping on random clusters near a wall; `lattice_delta.py`: the second term of $\Delta_{kl}$ in Lemma 6.11 with exact lattice Green functions |
| `bulk_poles.py` | Section 6.1 | poles of the continued exact bulk energy at the neutral thresholds |
| `vernier/` | Section 3.3, open problem 2 | Vernier–Jacobsen–Saleur $a_{2n-1}^{(2)}$ lattice mass ratios versus the companion-note masses; one-loop particle masses; $\beta^2_{\rm UV}$ |
| `hirota/` | Lemma 8.4, Propositions 8.5–8.8 | Hirota solitons and fluctuations for $d_n^{(1)}$, $e_n^{(1)}$ and foldings; transfer matrices; Wronskians; grid spectra and Jordan types; threshold zeros; heat traces |
| `oneloop/` | Section 8.9 | one-loop masses from transmission factors (Proposition 8.9); calibration; folded composites; comparison with exact results. `heat_trace_checks.py` checks the proofs of `prp-net-count` parts 2–3, `lem-born` and `prp-one-loop-mass` (output `heat_trace_checks.log`; with `--grid 768`, Fourier-grid heat traces, output `heat_trace_grid.log`). `phi4_check.py` tests the same results independently on the $\phi^4$ kink (output `phi4_check.log`) |
| `thimbles/` | Section 8.4 | reduced static model: critical points, intersection numbers, exact thimble integrals and the decomposition check (`check_th4.py`, output `out_check_th4.txt`); classical topological charges of Hirota solitons (`charges.py`, `charges_multi.py`) |
| `krein_bethe.py`, `scan_n4.py` | Proposition 9.3 | Krein-unitarity of fundamental-soliton Bethe–Yang transfer matrices; breaking by colour content |
| `sectors/` | Section 9.3, Proposition 9.4 | eigenvalue moduli of soliton, excited-soliton and breather amplitudes of the self-conjugate theories; $a_2^{(2)}$ kinks; exchange-statistics tests |
| `sectors/bloch/` | Section `sec-sector-status`, Propositions `prp-exchange-twist` and `prp-bloch-sectors` | graded exchange as a Bloch twist; breaking by charge sector and Bloch angle for $a_2^{(2)}$ kinks and spinor solitons; Krein pairing within sectors; see its `README.md` |
| `irf/` | Proposition 9.6 | RSOS face weights of $a_2^{(1)}$ kinks; Yang–Baxter, braiding unitarity, Hermitian analyticity, restricted transfer matrices (`checks.py ybe|ha|tm|vertex|signs`) |

Scripts in `sectors/` reuse `rsolve.py` and `dtw.py` from `soliton_rmatrix_code/`, copied alongside.

**Sign of $q$.** The spinor and $c_3^{(1)}$ scripts in `sectors/` use $q=e^{-i\pi\omega}$ by default, minus the physical value (book: `eq-q-physical`), and `irf/` uses $q=e^{i\pi\lambda}$, minus the value $-e^{i\pi\lambda}$ of Takács and Watts. Set `QSIGN=-1` for the physical sign; the default reproduces the stored outputs. In `irf/` the sign is applied by `qsign_patch.py`: at fixed spectral parameter, $q\to-q$ reverses the off-diagonal face weights and leaves the diagonal ones and the quantum dimensions unchanged. The $a_2^{(2)}$ scripts use the Takács–Watts $q$, and their spectra depend on $q^2$ only. Outputs at the physical sign are in `sectors/out/qsign/` and `irf/out_qsign/`; `sectors/qsign_checks.py` and `irf/qsign_compare.py` compare the two signs directly. The spinor setups refuse degenerate projectors (`sectors/lib.py`, `check_projectors`); for rank 3 this happens at $\omega=0.3$ and $0.7$, where $q^{20}=1$.

**Known quirk.** In `oneloop/summary.log` the $c_3^{(1)}$ row reads $[+0.0796,+0.159]$. Solitons 1 and 3 of $c_3^{(1)}$ are classically degenerate, and `summary.py` pairs solitons with particles by sorted mass, so the tie is resolved arbitrarily. The correct values, $-1/4\pi$ for both ratios, are in `oneloop/analysis_untwisted.log` and from `oneloop/calib_an_cn.py`.

**Known quirk.** `oneloop/calib_an_cn.py` prints $\Delta M_a/m_a$ for $a_n^{(1)}$ next to a column labelled $-h/2\pi$. The printed values agree with Hollowood's full formula $-[h/2\pi-\tfrac14\cot(\pi/h)]$, not with that column.
