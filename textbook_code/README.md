# Verification code for the textbook chapters

Scripts behind the symbolic and numerical checks of the textbook material in Parts I–II (`part1-background/`, `part2-algebra/`). Python 3 with sympy and mpmath. Each script runs on its own and prints one `PASS` or `FAIL` line per check.

| Script | Book location (ids) | What it checks | Runtime |
|---|---|---|---|
| `sine_gordon.py` | `eq-sg-s0`, `eq-sg-smatrix`, `eq-sg-breather` | crossing $S_0(\theta)=S_T(i\pi-\theta)$ of the integral representation, at $\lambda=1.5,1.2$ to $10^{-25}$; the breather solves the field equation | seconds |
| `sixvertex.py` | `sec-jimbo-equations`, `eq-six-vertex`, `eq-r-unitarity`, `eq-ybe-braid`, `eq-six-vertex-eigenvalue`, `sec-qg-crossing`, `sec-six-vertex` | the six-vertex $\check R$ solves Jimbo's equations; unitarity; braid Yang–Baxter in the printed order; eigenvalues; crossing data $x_c=q^2$, $C$, $c_u$; the sine-Gordon dictionary $z=e^{2\lambda\theta}$, $q=-e^{i\pi\lambda}$, including the gauge to the principal gradation | about a minute |
| `universal_r.py` | `eq-universal-r-sl2` | Drinfeld's universal R-matrix intertwines $\Delta$ and $\Delta^{\rm op}$ on spin $\frac12$ and spin 1; the reversed ordering fails | seconds |
| `spin1_tpg.py` | `eq-tpg-rule`, `exr-spin1` | $q$-Serre for the spin-1 evaluation representation; Jimbo's equations on spin $1\otimes1$ reproduce $\rho_1=\langle2\rangle$, $\rho_0=\langle2\rangle\langle1\rangle$ | about a minute |

What these scripts establish, what was checked only by hand, and what is still unchecked is recorded in `briefs/referee_parts_I_II.md`.
