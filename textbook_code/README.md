# Verification code for the textbook chapters

Scripts behind the symbolic and numerical checks of the textbook material in Parts I–IV (`part1-background/` to `part4-real-coupling/`) and the textbook chapters of Part V. Python 3 with numpy, sympy and mpmath. Each script runs on its own and prints one `PASS` or `FAIL` line per check.

| Script | Book location (ids) | What it checks | Runtime |
|---|---|---|---|
| `sine_gordon.py` | `eq-sg-s0`, `eq-sg-smatrix`, `eq-sg-breather` | crossing $S_0(\theta)=S_T(i\pi-\theta)$ of the integral representation, at $\lambda=1.5,1.2$ to $10^{-25}$; the breather solves the field equation | seconds |
| `sixvertex.py` | `sec-jimbo-equations`, `eq-six-vertex`, `eq-r-unitarity`, `eq-ybe-braid`, `eq-six-vertex-eigenvalue`, `sec-qg-crossing`, `sec-six-vertex` | the six-vertex $\check R$ solves Jimbo's equations; unitarity; braid Yang–Baxter in the printed order; eigenvalues; crossing data $x_c=q^2$, $C$, $c_u$; the sine-Gordon dictionary $z=e^{2\lambda\theta}$, $q=-e^{i\pi\lambda}$, including the gauge to the principal gradation | about a minute |
| `universal_r.py` | `eq-universal-r-sl2` | Drinfeld's universal R-matrix intertwines $\Delta$ and $\Delta^{\rm op}$ on spin $\frac12$ and spin 1; the reversed ordering fails | seconds |
| `spin1_tpg.py` | `eq-tpg-rule`, `exr-spin1` | $q$-Serre for the spin-1 evaluation representation; Jimbo's equations on spin $1\otimes1$ reproduce $\rho_1=\langle2\rangle$, $\rho_0=\langle2\rangle\langle1\rangle$ | about a minute |
| `toda_classical.py` | `eq-lax-pair`, `eq-toda-eom`, `eq-spin3-current`, `sec-masses-couplings`, `prp-dorey-rule`, `eq-area-rule`, `eq-hirota-bilinear`, `eq-an-soliton`, `eq-n-soliton`, `eq-interaction-coefficient`, `eq-soliton-mass`, `eq-theta-data`, `exr-breather-real` | Lax pair ⇔ field equation for $a_2^{(1)}$; sinh-Gordon spin-3 current; masses = Perron–Frobenius vector and sum rule; Dorey's rule and $\lvert C_{abc}\rvert=\frac{4\beta}{\sqrt h}\Delta_{abc}$ for $e_7,e_8$; one- and two-soliton $a_n^{(1)}$ solutions; soliton energy $2hm_a/\beta^2$ at complex $\xi$; $\Theta$ on Hirota data; sine-Gordon breather reality | a few minutes |
| `smatrix_simply_laced.py` | `eq-blocks`, `eq-dorey-formula`, `eq-an-smatrix`, `sec-strong-weak`, `sec-tree-level` | Dorey's formula: unitarity, crossing, $B\to2-B$, $e_8$ $S_{11}$, $a_n$ closed form, bootstrap at every $e_7$/$e_8$ fusing; sinh-Gordon tree level. Requires `toda_classical.py` in the same directory | about a minute |
| `semiclassics.py` | `eq-dhn-psi`, `eq-ccg`, `eq-hollowood-mass`, `sec-picard-lefschetz`, `prp-thimble-gaussian` | one-loop masses of the sine-Gordon and $\phi^4$ kinks and of the $a_n^{(1)}$ solitons from transmission factors; the zero-dimensional $i\phi^3$ integral; Gaussian over a Jordan block; Jordan blocks of $-\partial^2+ge^{ix}$ | seconds |

What these scripts establish, what was checked only by hand, and what is still unchecked is recorded in `briefs/referee_textbook.md`.
