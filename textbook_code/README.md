# Verification code for the textbook chapters

Scripts behind the symbolic and numerical checks of the textbook material in Parts I–IV (`part1-background/` to `part4-real-coupling/`) the textbook parts of Parts V–VI, and appendix A. Python 3 with numpy, sympy and mpmath. Each script runs on its own and prints one `PASS` or `FAIL` line per check.

| Script | Book location (ids) | What it checks | Runtime |
|---|---|---|---|
| `sine_gordon.py` | `eq-sg-s0`, `eq-sg-smatrix`, `eq-sg-breather` | crossing $S_0(\theta)=S_T(i\pi-\theta)$ of the integral representation, at $\lambda=1.5,1.2$ to $10^{-25}$; the breather solves the field equation | seconds |
| `sixvertex.py` | `sec-jimbo-equations`, `eq-six-vertex`, `eq-r-unitarity`, `eq-ybe-braid`, `eq-six-vertex-eigenvalue`, `sec-qg-crossing`, `sec-six-vertex` | the six-vertex $\check R$ solves Jimbo's equations; unitarity; braid Yang–Baxter in the printed order; eigenvalues; crossing data $x_c=q^2$, $C$, $c_u$; the sine-Gordon dictionary $z=e^{2\lambda\theta}$, $q=-e^{i\pi\lambda}$, including the gauge to the principal gradation | about a minute |
| `universal_r.py` | `eq-universal-r-sl2` | Drinfeld's universal R-matrix intertwines $\Delta$ and $\Delta^{\rm op}$ on spin $\frac12$ and spin 1; the reversed ordering fails | seconds |
| `spin1_tpg.py` | `eq-tpg-rule`, `exr-spin1` | $q$-Serre for the spin-1 evaluation representation; Jimbo's equations on spin $1\otimes1$ reproduce $\rho_1=\langle2\rangle$, $\rho_0=\langle2\rangle\langle1\rangle$ | about a minute |
| `toda_classical.py` | `eq-lax-pair`, `eq-toda-eom`, `eq-spin3-current`, `sec-masses-couplings`, `prp-dorey-rule`, `eq-area-rule`, `eq-hirota-bilinear`, `eq-an-soliton`, `eq-n-soliton`, `eq-interaction-coefficient`, `eq-soliton-mass`, `eq-theta-data`, `exr-breather-real` | Lax pair ⇔ field equation for $a_2^{(1)}$; sinh-Gordon spin-3 current; masses = Perron–Frobenius vector and sum rule; Dorey's rule and $\lvert C_{abc}\rvert=\frac{4\beta}{\sqrt h}\Delta_{abc}$ for $e_7,e_8$; one- and two-soliton $a_n^{(1)}$ solutions; soliton energy $2hm_a/\beta^2$ at complex $\xi$; $\Theta$ on Hirota data; sine-Gordon breather reality | a few minutes |
| `smatrix_simply_laced.py` | `eq-blocks`, `eq-dorey-formula`, `eq-an-smatrix`, `sec-strong-weak`, `sec-tree-level` | Dorey's formula: unitarity, crossing, $B\to2-B$, $e_8$ $S_{11}$, $a_n$ closed form, bootstrap at every $e_7$/$e_8$ fusing; sinh-Gordon tree level. Requires `toda_classical.py` in the same directory | about a minute |
| `semiclassics.py` | `eq-dhn-psi`, `eq-ccg`, `eq-hollowood-mass`, `sec-picard-lefschetz`, `prp-thimble-gaussian` | one-loop masses of the sine-Gordon and $\phi^4$ kinks and of the $a_n^{(1)}$ solitons from transmission factors; the zero-dimensional $i\phi^3$ integral; Gaussian over a Jordan block; Jordan blocks of $-\partial^2+ge^{ix}$ | seconds |
| `rsos_qdims.py` | `exr-quantum-dims`, exercises of chapters 22–25 | sign criteria for $U_q(sl_2)$ quantum dimensions at $q=e^{i\pi r/p}$; $T$, $H$ and mass identities for $g_2^{(1)}$, $d_4^{(3)}$, $c_n^{(1)}$, $a_{2n-1}^{(2)}$, $e_6^{(2)}$ | seconds |
| `algebra_data.py` | `tbl-app-masses` | classical masses of every untwisted affine Toda theory, by folding for the non-simply-laced ones; closed forms. Requires `toda_classical.py` | seconds |

What these scripts establish, what was checked only by hand, and what is still unchecked is recorded in `briefs/referee_textbook.md`.

## Figures (`figures/`)

Scripts that draw the matplotlib figures of the introduction and of Part I and check the formulas they plot. Python 3 with numpy, scipy and matplotlib. Each writes an SVG into the book tree (`figures/` for the introduction, `part1-background/figures/` for chapters 2–4) and prints one `PASS` or `FAIL` line per check. The SVGs are committed, so building the book needs no Python. `style.py` holds the shared palette (mid-tone colours and grey text, transparent background, so that the figures read on both the light and the dark theme). Run one script with `python3 figures/<name>.py`, or all of them with `for f in figures/[ic]*.py; do python3 $f; done`.

| Script | Figure id | Book location | What it checks | Runtime |
|---|---|---|---|---|
| `figures/intro_f1_reading_paths.py` | `fig-reading-paths` | `sec-reader-guide` | the chapter lists of the four paths agree with the "Reading paths" bullets of `01-introduction.qmd` | seconds |
| `figures/intro_f2_potentials.py` | `fig-two-regimes` | `01-introduction.qmd` | $\sum_j\vec\alpha_j=0$ and $\vec\alpha_j^2=2$ for $a_2^{(1)}$; $\mathrm{Re}\,V_{\rm imag}$ vanishes exactly at the vacuum lattice and is positive elsewhere; $V_{\rm real}\ge0$ | seconds |
| `figures/ch2_f3_strip_map.py` | `fig-theta-strip` | `eq-mandelstam-s`, `exr-strip` | the map $\theta\mapsto s$ is injective on the open strip and sends its boundary lines to the stated cuts and gap | seconds |
| `figures/ch2_f4_factorization.py` | `fig-factorization` | `sec-factorization`, `eq-yang-baxter` | the two world-line configurations have collision orders $12,13,23$ and $23,13,12$ | seconds |
| `figures/ch2_f5_triangle_poles.py` | `fig-bound-state-poles` | `eq-bound-state-mass`, `eq-fx-block`, `sec-lee-yang` | residue of $f_x$ at $i\pi x$ equals $2i\tan\pi x$; $f_{2/3}=f_{1/3}$; the bootstrap identity for $f_{2/3}$; the mass-triangle relation | seconds |
| `figures/ch3_f6_coupling_axis.py` | `fig-sg-coupling-axis` | `sec-sg-coleman`, `prp-sg-rsos` | number of breathers against the thresholds $8\pi/(n+1)$; $\xi=\pi/\lambda$; Lee–Yang point $\xi=2\pi/3$ with one breather; unitary series without breathers | seconds |
| `figures/ch3_f7_kink_breather.py` | `fig-sg-solutions` | `eq-sg-kink`, `eq-sg-breather`, `eq-sg-classical-mass` | kink and breather solve the sine-Gordon equation (finite differences); kink energy $8m_0/\beta^2$; breather energy $2M_{\rm cl}\sin u$ | seconds |
| `figures/ch4_f8_dimer.py` | `fig-dimer-krein` | `eq-dimer`, `sec-dimer` | $JH^\dagger J=H$; eigenvalues $a\pm\sqrt{s^2-g^2}$; Krein norms $\pm2s\sqrt{s^2-g^2}$; neutral eigenvectors for $g>s$; nilpotent $H-a$ at $g=s$ | seconds |
| `figures/ch4_f9_krein_geometry.py` | `fig-krein-geometry` | `prp-krein-spectrum` item 4 | $H=i\sigma_x$ is $J$-self-adjoint and not Hermitian; $e^{-iHt}=\cosh t+\sigma_x\sinh t$; Krein norm conserved; Hilbert norm $\cosh2t$ | seconds |

The part-level dependency graph of the introduction (`fig-dependencies`) is a mermaid block drawn from the text of `sec-reader-guide`. The reading-paths figure (`intro_f1_reading_paths.py`) checks its chapter lists against the "Reading paths" bullets of the same section.
