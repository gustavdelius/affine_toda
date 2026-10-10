
## Sign of q

The scripts use q = e^{-i pi omega} by default. This is minus the physical value q = -e^{-i pi omega} = e^{-4 pi^2 i/beta^2} of Bernard and LeClair for the coproduct used here (book: `eq-q-physical`, section `sec-crossing-sign`). Set `QSIGN=-1` to use the physical value; the default reproduces the outputs stored in the repository.

- Independent of the sign (they depend on q^2 only): representations, R-matrix eigenvalues, crossing points, pole positions, residue ranks, breather-breather amplitudes.
- Dependent on the sign: the crossing factor (`sign_n*.py`: (-1)^n at the default, +1 at `QSIGN=-1`), the c_3^(1) fusion factors k_32 and k_12 (`kcheck.py` picks the closed forms that match `QSIGN`), and the sign of S[B^3,3] (`bootstrap_exact.py`).
- `spinorn.py` caches its projectors per omega and sign of q.
- Outputs at the physical q are in `smatrix/qsign_out/`.

## smatrix/ (scalar factors)
- `common.py`   : spinor representation at complex q, closed-form R_33 with projectors
- `crossing.py` : charge conjugation, crossing point and crossing factor c(x) at unimodular q
- `s33.py`      : full S_33 (crossing factor, Gamma infinite product, CDD factor); unitarity, crossing, pole scan
- `fuse11*.py`  : S_11 by fusion of four S_33 factors; residues, pole orders, breathers
(copy `rsolve.py` and `dtw.py` from the parent folder into smatrix/ before running)
- `sign_n1.py`, `sign_n2.py`, `sign_n4.py` : the crossing-sign test for n = 1 (sine-Gordon), 2 (GMW) and 4
- `nocdd.py` : S_11 and S_22 by fusion from S_33 without the CDD factor (the corrected construction)
- `allS.py` : general fusion machinery for every soliton pair (env OMEGA, JTRUNC set coupling and truncation)
- `kfit.py`, `kcheck.py` : fusion factors k_ab(x), fitted and verified in closed form
- `poles.py`, `spectrum.py` : exhaustive pole tables (order, residue rank, masses) and the aggregated spectrum
- `breathers.py`, `breathers2.py` : lowest breathers by fusion; breather–soliton amplitudes as block products
- `bootstrap_exact.py` : exact bootstrap in omega for all breather–soliton and breather–breather amplitudes
- `verify_w.py`, `jscale.py` : second-coupling verification and 1/J truncation scaling
- `dgz_compare.py` : exact comparison of the six breather–breather amplitudes with the Delius–Grisaru–Zanon c_3 particle S-matrix (all identical)
- `spinorn.py` : spinor soliton S-matrix for general n (env NSPIN, OMEGA, JTRUNC); fusion/rank checks
- `breather_n.py` : lowest breather of soliton n; exact breather amplitudes; comparison with DGZ c_n (identical to DGZ for n = 2, 3, 4 with the rule a_j = j w + (j-1)/2)
- `bspin.py`, `bspin_struct.py` : U_q(b_n^(1)) spinor representation (a_{2n-1}^(2) solitons); R-matrix structure, crossing
- `aspin_S.py`, `aspin_n.py` : a_{2n-1}^(2) spinor soliton S-matrix; lift scan; breather amplitude identical to DGZ (n = 3, 4)
- `f4rep.py` : 26-dim representation of U_q(f_4^(1)) solved from all relations (e_6^(2) solitons)
- `f4R2.py`, `f4check.py` : its R-matrix eigenvalues via highest-weight vectors (closed form), crossing point
- `e6sol.py` : e_6^(2) soliton S-matrix (U_q(f_4^(1)) on 26); breather amplitude identical to CDS S_11
- `e6tw.py`, `e6R.py` : U_q(e_6^(2)) on 27 = 26+1 and its R-matrix with multiplicity blocks (f_4^(1) solitons)
- `f4sol.py`, `f4breather.py`, `f4lift*.py`, `f4search.py` : f_4^(1) soliton S-matrix and the (so far unsuccessful) breather search
- `g2alg.py`, `g2sol*.py`, `g2br*.py` : U_q(g_2^(1)) on 7 and U_q(d_4^(3)) on 8; re-derivation of Takacs's d_4^(3) and g_2^(1) results
- `d42vec.py`, `vec_calib.py` : the c_3 soliton-1 multiplet (U_q(d_4^(2)) 8 = 7+1); minimal construction reproduces the fused S_11
- `f4rescan.py`, `f4boot.py` : complete pole scan of the 27 R-matrix; self-fusion bootstrap of the f_4^(1) amplitude
- `f4particle2.py`, `f4order.py`, `f4cdd.py` : particle-bootstrap closure, tree-level comparison with CDS, CDD search
- `oneloop_ratio.py`, `oneloop_e6.py` : one-loop dH/B_tree from the Toda Lagrangian (c_3: 1, f_4: 3, e_6^(2): -3)
- `f4cdd_lin.py` : exact linear-lattice search for a CDD factor on the f_4^(1) 27 amplitude (none exists; residual ~ 1/window)
- `direct_lagrangian.py`, `direct_e6.py` : CDS-free comparison of the breather amplitudes with the Toda Lagrangian (tree shape, coupling map, one-loop masses): c_3 and e_6^(2) agree exactly, f_4^(1) off by 3/2
- `f4cand.py` : one-loop selection of lifts: candidates T = 6w+9/2 and T = 6w-3/2 give exact one-loop masses under the gradation map; tree strength structural (1/w = beta^2/8pi)
- `f4halftower.py`, `f4hyp.py`, `f4dbg.py` : exact block comparison with CDS (same particle poles, coupling-dependent shifts 1 instead of 1/2); bound-state hypothesis test (negative)
- `f4gamma*.py`, `f4moments.py` : test of the half-integer Gamma-tower conjecture (closed-form pair factor W_a, exact lattice solve: no solution) and moment classification (Barnes-type growth needed on 16 of 24 lines)
