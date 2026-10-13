# Part VIII (defects): source digests

Digests of the defect literature written on 2026-10-11 by the sources agent of Part VIII, from pdftotext of the arXiv versions. Equation numbers are those of the sources. They record conventions, formulas the book re-derives, and misprints noticed; they are working notes, not checked results. The convention warnings at the top of the conventions brief `part8_defects_conventions.md` take precedence.

# Source digests for Part VIII (defects) — compiled 2026-10-11

Equation numbers are the source's. Notation is transcribed from pdftotext and re-typeset by hand;
check against the PDF before quoting a formula verbatim.

## bcz2004 — Bowcock, Corrigan, Zambon, "Classically integrable field theories with defects", hep-th/0305022 (IJMPA 19S2 (2004) 82)
- Lagrangian (1.2): L = θ(−x)[½(∂φ1)² − V1(φ1)] + θ(x)[½(∂φ2)² − V2(φ2)] + δ(x)[½(φ1 ∂tφ2 − φ2 ∂tφ1) − B(φ1,φ2)].
  φ1 lives on x<0, φ2 on x>0. Metric signature implied (+,−): (∂φ)² = (∂tφ)² − (∂xφ)².
- Defect conditions (1.5)-(1.6) at x=0: ∂xφ1 − ∂tφ2 = −∂B/∂φ1 ; ∂xφ2 − ∂tφ1 = +∂B/∂φ2.
- Spin-3 conservation (2.6)-(2.9) forces λ1=λ2=λ and B of the form (2.10) B = a e^{λ(φ1+φ2)} + b e^{λ(φ1−φ2)} + c e^{−λ(φ1−φ2)} + d e^{−λ(φ1+φ2)},
  with V_p = A_p e^{2λφp} + B_p e^{−2λφp}, A1=2λ²ab, A2=2λ²ac, B1=2λ²cd, B2=2λ²bd (2.11). Also lists free-massive (quadratic B),
  free-massless and Liouville cases ((2.12): Liouville/free-massless B = (2m/β²) e^{βφ2/√2}(σ e^{βφ1/√2} + σ^{-1} e^{−βφ1/√2})).
- Sinh-Gordon normalisation: λ = β/√2, a = mσ/β², b = m/(β²σ):
  (2.13) V_p = (m²/β²)(e^{√2βφp} + e^{−√2βφp}) [so the particle mass is 2m: V''(0)=4m²],
  B = (mσ/β²)(e^{β(φ1+φ2)/√2} + e^{−β(φ1+φ2)/√2}) + (m/(β²σ))(e^{β(φ1−φ2)/√2} + e^{−β(φ1−φ2)/√2}) (± signs correlated). One defect parameter σ.
- Remarks: B satisfies the auto-Bäcklund conditions ∂²B/∂φ1² = ∂²B/∂φ2², (∂B/∂φ1)² − (∂B/∂φ2)² = 2(V1−V2); defect conditions = Bäcklund transformation frozen at x=0.
  Modified momentum P + P_B is conserved (2.15)-(2.16). Topological charge (2.14) not obviously conserved.
- Linearised transmission (Sec. 3, eq. (3.1)): with k = 2m sinh θ, ω = 2m cosh θ, σ = e^p:
  R = 0, T = −(2i cosh θ − (σ − 1/σ))/(2i sinh θ − (σ + 1/σ)) = −i sinh((θ+p)/2 − iπ/4)/sinh((θ+p)/2 + iπ/4). (I checked numerically that the two forms agree; |T|=1.)
- Sine-Gordon (β=1/√2, m=1/2 → ∂²φ = −sin φ), defect conditions (3.4): ∂xφ1 − ∂tφ2 = −σ sin((φ1+φ2)/2) − σ^{-1} sin((φ1−φ2)/2), ∂xφ2 − ∂tφ1 = σ sin((φ1+φ2)/2) − σ^{-1} sin((φ1−φ2)/2).
  One-soliton e^{iφa/2} = (1 − iE_a)/(1 + iE_a), E_a = C_a e^{αx+βt}, α = cosh θ, β = sinh θ (3.5); delay (3.6) C2 = ((e^θ + σ)/(e^θ − σ)) C1.
  Absorption/emission at e^θ = |σ| (C2 = 0 for σ<0, C2 = ∞ for σ>0).
- Sec. 4: defect Lax pair (4.1)-(4.2) by overlapping regions a<x<b and gauge matrix K = e^{−φ2·H/2} K̄ e^{φ1·H/2} (4.3)-(4.4); for A1, K̄ = 1 + (ρ/λ)σ_x (4.6), B̃ = ρ^{-1}(e^{α(φ1−φ2)/2} + e^{α(φ2−φ1)/2}) (4.7). Claims (later corrected in bcz2004b) that only sinh-Gordon/Liouville admit this.
- Misprints: Sec. 3 linearised ansatz writes "φ1 = e^{−iωt}T(k)e^{ikx}" for the x>0 field (should be φ2); (2.16) second condition prints ∂²B/∂φ1² = ∂²B/∂φ1² (should be φ1 vs φ2).
  Sec. 3 text says "using the defect conditions (1.3)" for (2.15) — means (1.5)-(1.6).

## bcz2004b — Bowcock, Corrigan, Zambon, "Affine Toda field theories with defects", hep-th/0401020 (JHEP 0401 (2004) 056)
- Convention switch: φ on x<0, ψ on x>0 (Sec. 2). β and m set to 1 ("may be reinstated by rescalings"). Bulk potential (3.5) V = Σ_{i=0}^r n_i e^{α_i·φ} (no "−h"), simply-laced roots α_i² = 2, m_i = n_iα_i²/2.
- Free field (2.1)-(2.10): L_D = δ(x)[½(φ∂tψ − ψ∂tφ) − B], conditions (2.2) ∂xφ = ∂tψ − ∂B/∂φ, ∂xψ = ∂tφ + ∂B/∂ψ; momentum is total-time-derivative iff (2.4) B_φφ = B_ψψ and (2.5) ½B_φ² − ½B_ψ² = V − W.
  Free massive: B = (mλ/4)(φ+ψ)² + (m/4λ)(φ−ψ)² (2.7); conserved E+B, P−(…) (2.8)-(2.10): P_{±1} = E ± P = … + (mλ^{±1}/2)(φ ± ψ)² — λ has spin one.
  Transmission R=0, T = −(ik + (m/2)(λ + 1/λ))/(iω + (m/2)(λ − 1/λ)).
- General ATFT defect Lagrangian (3.1): δ(x)[½φE∂tφ + φD∂tψ + ½ψF∂tψ − B(φ,ψ)]; defect conditions (3.2) ∂xφ − E∂tφ − D∂tψ + ∇_φB = 0, ∂xψ − D^T∂tφ + F∂tψ − ∇_ψB = 0.
- Lax pair (3.3) a_t = ½[∂xφ·H + Σ√m_i e^{α_i·φ/2}(λE_{α_i} − λ^{-1}E_{−α_i})], a_x = ½[∂tφ·H + Σ√m_i e^{α_i·φ/2}(λE_{α_i} + λ^{-1}E_{−α_i})]; [E_αi, E_−αj] = δ_ij 2α_i·H/α_i².
  Overlap construction (3.6)-(3.7) ∂tK = K â_t^> − â_t^< K, K = e^{−½H·(Eφ+Dψ)} K̃ e^{½H·(D^Tφ − Fψ)} (3.8).
- "Maximal" defect: (3.11) E = F = 1 − D, D + D^T = 2. K̃ = 1 + k1/λ + …, k1 = Σ c_i E_{−α_i} (3.14).
  D must satisfy α_j D = −α_{Π^{-1}(j)} D^T with Π a Dynkin-diagram rotation (3.19)-(3.21): D = 2(1−π)^{-1}. Only a_r^{(1)} with Π the cyclic rotation works: (3.27) Dα_i = 2(λ_i − λ_{i+1}), Dα_r = 2λ_r, Dα_0 = −2λ_1; (3.28) D = 2Σ_{i=1}^r (λ_i − λ_{i+1})λ_i^T; (3.30) 1 − D = Σ(λ_{i+1}λ_i^T − λ_iλ_{i+1}^T) (antisymmetric).
- Defect potential (3.32): B = Σ_{i=0}^r [ d e^{½α_i·(D^Tφ + Dψ)} + d^{-1} e^{½α_i·D(φ−ψ)} ], single parameter d (=σ for a_1^{(1)} with D=1).
  K̃ in the (r+1)-dim rep: K̃ = 1 + (d/λ) Σ_{i=0}^r E_{−α_i} (3.39) — exact.
- Conserved momentum P − U (4.5)-(4.6), energy E + B (4.7); light-cone versions involve only f_i (d-terms) or g_i (d^{-1}-terms) — defect parameter transforms under boosts like rapidity (d → d e^{−ϑ}, my paraphrase).
- Spin ±2 (Sec. 5): with (3.11) only the spin +2 charge survives (modified); the alternative solution (5.16) E = F = 1 + D^T, D + D^T = −2 preserves spin −2 instead (θ → −θ in delays).
- Sec. 5.1 re-derives D, E, F by requiring the defect conditions to be a Bäcklund transformation: M² = 1, MN + NM = 0, B_φ² − B_ψ² = 2(V − W) (5.23)-(5.25); E = Σ(λ_{i+1}λ_i^T − λ_iλ_{i+1}^T) (5.42).
- Solitons (Sec. 6): complex φ with Hirota τ_j = 1 + E_a ω^{aj}, ψ uses σ_j = 1 + Λ E_a ω^{aj} (6.1), E_a = e^{a_a x + b_a t + ξ}, (a_a, b_a) = m_a(cosh θ, sinh θ), m_a² = 4 sin²(aπ/(r+1)), ω = e^{2πi/(r+1)}.
  Delay (6.2): Λ = (i e^θ + d e^{iζ_a})/(i e^θ + d e^{−iζ_a}), ζ_a = πa/(r+1). Self-conjugate soliton: Λ = (e^θ + d)/(e^θ − d). Conjugate pair Λ(ā) = conj (6.3).
  No soliton-number-changing defect found; absorption at special real rapidity for self-conjugate solitons.
- Misprints: (3.37) "c_ii = ," (value d²/2 lost in layout only); text after (3.32) refers to "(3.33)" for order-1/λ terms (means (3.15)); (5.17) second line should probably have D^T ↔ D roles mirrored (check against (5.16)); Section 5 refers to "(5.18)" before it is defined.

## bcz2005 — Bowcock, Corrigan, Zambon, "Some aspects of jump-defects in the quantum sine-Gordon model", hep-th/0506169 (JHEP 0508 (2005) 023)
- Bulk (2.1): L = ½(u_t² − u_x²) − (m²/β²)(1 − cos βu). u on x<0, v on x>0.
- Defect Lagrangian (2.2): L = θ(−x)L_u + θ(x)L_v ± δ(x)[½(u v_t − v u_t) − B(u,v)]. The pdftotext layer shows "− δ(x)", but the conditions (2.4) are only reproduced with "+ δ(x)" (I checked the variation by hand: +δ gives u_x − v_t = −∂B/∂u). Treat the printed minus as a probable misprint and quote with "+".
  (2.3) B = −(2mσ/β²) cos(β(u+v)/2) − (2m/(σβ²)) cos(β(u−v)/2).
  (2.4), scaled (m=β=1): ∂²u = −sin u (x<0), ∂²v = −sin v (x>0); at x=0: u_x − v_t = −σ sin((u+v)/2) − σ^{-1} sin((u−v)/2), v_x − u_t = σ sin((u+v)/2) − σ^{-1} sin((u−v)/2).
- Conserved E = E(u)+E(v)+B and P = P(u)+P(v)+U with U_u = B_v, U_v = B_u: U = −2σ cos((u+v)/2) + (2/σ) cos((u−v)/2) (unnumbered).
  Defect potential has twice the bulk period. Lagrangian not invariant under P or T separately; parity must be combined with u→−u (or v) and σ→1/σ (η→−η).
- Soliton (2.5): e^{iu/2} = (1+iE)/(1−iE), E = e^{ax+bt+c}, a = cosh θ, b = −sinh θ (θ>0: moving to +x). Soliton mass 8m/β² (scaled: 8). Anti-soliton: E→−E.
  Outgoing (2.6) e^{iv/2} = (1+izE)/(1−izE), delay (2.7) z = (e^{−θ} + σ)/(e^{−θ} − σ) = coth((η − θ)/2) with σ = e^{−η}.
  η<0: soliton always flips to anti-soliton; η>0: θ<η delayed, θ>η flipped, θ=η absorbed (final static u=0, v=2π stores the soliton's E and P). Two identical defects: z = coth²((η−θ)/2) = the time delay from a soliton of rapidity η, so a defect is like "half a soliton" (energy scale 4m/β²).
  Energetics (2.11)-(2.12): defect energy ∓4 cosh η for even/odd jump; emission of a soliton with θ = η from a 2π-jump is allowed energetically.
- Linearised transmission (2.8)-(2.9): T_even(θ,η) = −i sinh((θ−η)/2 − iπ/4)/sinh((θ−η)/2 + iπ/4) = conj T_odd; T(θ)T(−θ) ≠ 1; with T_P(θ,η) = −T(−θ,−η), T T_P = 1 (2.10).
- Quantum conventions (KL conventions, Sec. 3): soliton S-matrix (3.2)-(3.3) S = ρ(Θ)[[a,0,0,0],[0,c,b,0],[0,b,c,0],[0,0,0,a]] with a = qx2/x1 − x1/(qx2), b = x1/x2 − x2/x1, c = q − 1/q, (3.4) q = e^{iπγ}, x_p = e^{γθ_p}, ρ in (3.6)-(3.7); coupling (3.8) 1/γ = β²/(8π − β²) (so γ = 8π/β² − 1; free fermion at γ=1).
- Triangle (transmission Yang–Baxter) relation (3.9): S_{kl}^{mn}(Θ) T_{nα}^{tβ}(θ1) T_{mβ}^{sγ}(θ2) = T_{lα}^{nβ}(θ2) T_{kβ}^{mγ}(θ1) S_{mn}^{st}(Θ).
  Block notation (3.10) T = [[A,B],[C,D]] by soliton/antisoliton, entries labelled by defect topological charge α; A, D diagonal, B, C shift α by ∓2 (3.11).
  Solution (ε = −1 chosen; ε = +1 excluded by crossing/unitarity): (3.48) A = Q^α A(θ)δ, D = νQ^{−α}A(θ)δ, B ∝ p x A, C ∝ p x A, with Q² = −q, p = e^{−γη}, x = e^{γθ}.
  Symmetric gauge (3.51): T = f(q,x)[[ν^{−1/2}Q^α δ_α^β, q^{−1/2}e^{γ(θ−η)}δ_α^{β−2}],[q^{−1/2}e^{γ(θ−η)}δ_α^{β+2}, ν^{1/2}Q^{−α}δ_α^β]].
  f satisfies (3.38) f(q,x)f(q,qx)(1 + p²x²) = 1 = f̄ f(1 + p²x²) [crossing-type and unitarity]; KL minimal solution (3.43)-(3.44) f = e^{iπ(1+γ)/4} r(x)/((1+ipx) r̄(x)), r(x) = Π_{k≥0} Γ(kγ+1/4−iγθ̃/2π)Γ((k+1)γ+3/4−iγθ̃/2π)/[Γ((k+½)γ+1/4−iγθ̃/2π)Γ((k+½)γ+3/4−iγθ̃/2π)], θ̃ = θ − η; (3.47) compact Gamma-product form.
  Annihilation-pole (crossing) condition (3.34)-(3.35): (ν + λμq x²/σ²) A(θ)A(θ+iπ) = 1.
  Topological-charge phases from the time-derivative term (3.52)-(3.56): factor exp(iπ(aδv − bδu)/β); Q^{(b∓a)/2}.
- Unstable bound state (Sec. 4): pole of T_even at θ = η − iπ/(2γ) (the apparent pole at η + iπ/(2γ) cancels), (4.1) ΔE = m_s cosh(η − iπ/2γ); (4.3)-(4.4) (ΔE, ΔP) is a complex rotation of the soliton momentum, stays on the soliton mass shell. Residue (4.7). Odd-charge transmission matrix T̂ via bootstrap (4.8)-(4.16): ν̂ = ν, p̂² = q²p² (4.15); T̂ is not unitary.
- Breathers (Sec. 5): poles Θ = i(π − nπ/γ) (5.1), m_n = 2m_s sin(nπ/2γ) (5.2); bootstrap (5.3)-(5.4). Lightest breather (5.7): ¹T(θ) = i sin(π/4 + i(θ−η)/2)/sin(π/4 − i(θ−η)/2) = −i sinh((θ−η)/2 − iπ/4)/sinh((θ−η)/2 + iπ/4) — coupling independent, equals the classical T_even.
  Block notation (5.9) (z)_θ = sinh(θ/2 + iπz/4)/sinh(θ/2 − iπz/4); S_11 = −1/((2/γ)(2−2/γ)) style (5.8); ²T = −1/((1+1/γ)_θ̃ (1−1/γ)_θ̃) (5.10); general odd/even breathers (5.12)-(5.13); poles of ²T at θ̃ = −iπ(1 ± 1/γ)/2 with ΔE± (5.11).
- Moving defects (Sec. 6): action (6.1) with m_D√(1−ż²); B = √(1−ż²) C, C(σ) = σF(u+v) + σ^{-1}F(u−v) (6.15),(6.18); σ = σ0 √((1+ż)/(1−ż)) (6.19); E, P conserved as for a defect at rest with σ0 (6.20),(6.23); m_D z̈ = 0; (6.26) σ = e^{−η} = e^{χ − η0}, ż = tanh χ. Defect–defect "U" matrix and triangle relation (6.28)-(6.30); Bianchi permutability.
- Misprints/caveats: sign in front of δ(x) in (2.2) (see above); BCZ2003 used b = +sinh θ (θ>0 left-moving) whereas here b = −sinh θ — rapidity sign flips between papers.

## cz2004 — Corrigan, Zambon, "Aspects of sine-Gordon solitons, defects and gates", hep-th/0407199 (J. Phys. A 37 (2004) L471)
- Short letter. Same defect conditions as bcz2004 (3.4) = here (2.1)-(2.3), φ1 on x<0, φ2 on x>0, ∂²φ = −sin φ.
- Soliton (2.4): e^{iφa/2} = (1 + E_a)/(1 − E_a), E_a = e^{α x + β t + γ_a}, α = cosh θ, β = sinh θ (velocity −tanh θ: solitons come from the right), Im γ = π/2.
  Delay (2.5): γ1 = γ2 − ln((e^θ + σ)/(e^θ − σ)).
- Classification for σ > 1, soliton from the right: e^θ > σ delayed; e^θ = σ swallowed (defect stores one unit of topological charge); e^θ < σ converted to anti-soliton (defect stores two units); opposite direction never absorbed. σ < 1: soliton/antisoliton roles swap.
- Topological charge bookkeeping (2.6) Q = φ2(∞) − φ1(−∞) + φ1(0) − φ2(0). Symmetry: defect conditions invariant only under pairs of {t→−t, φ1→−φ1, φ2→−φ2} together with σ→1/σ — no time-reversed emission.
- Two-soliton Hirota check (2.7)-(2.9): each soliton delayed independently; at most one can be absorbed.
- Logic gates (Sec. 3, Tables 1-3): NOT, XOR (with feedback raising σ), NAND/Toffoli with three solitons.
- Misprints: after (2.4) "Im γ = iπ/2" should read Im γ = π/2; (2.7) "e^{iφ/2} = ln(τ0/τ1)" should be e^{iφ/2} = τ0/τ1 (or iφ/2 = ln(τ0/τ1)).

## cz2007 — Corrigan, Zambon, "On purely transmitting defects in affine Toda field theory", 0705.1066 (JHEP 0707 (2007) 001)
- Bulk (2.1): L = ½∂φ·∂φ − (m²/β²) Σ_{j=0}^r (e^{βα_j·φ} − 1), |α_j|² = 2, α_0 = −Σα_j; real masses (2.2) m_a = 2m sin(πa/h), h = r+1.
  Complex (imaginary-coupling) theory: β → iβ (2.9) L_φ = ½∂φ·∂φ + (m²/β²)Σ(e^{iβα_j·φ} − 1). Vacua φ = 2πw/β, w ∈ weight lattice (2.3); topological charge (2.4) Q = (β/2π)[φ(∞) − φ(−∞)].
  Solitons (2.6) φ^{(a)} ∝ (i/β) Σ_j α_j ln(1 + E_a ω^{aj}), E_a = e^{a_a x − b_a t + ξ_a}, (a_a, b_a) = m_a(cosh θ, sinh θ), ω = e^{2πi/h} [the text layer shows a stray "m²" in the prefactor — check PDF]; soliton masses (2.7) M_a = 2h m_a/β².
- Defect Lagrangian (2.8): δ(x)[½(φ·E∂tφ + φ·D∂tψ − ∂tφ·Dψ + ψ·E∂tψ) − B], E antisymmetric, D = 1 − E,
  (2.10) B = −(m/β²) Σ_{j=0}^r [σ e^{iβα_j·(D^Tφ + Dψ)/2} + σ^{-1} e^{iβα_j·D(φ − ψ)/2}].
  (2.11) α_k·Dα_j = 2 (k=j), −2 (k=π(j)), 0 otherwise; D + D^T = 2. "Clockwise" permutation α_{π(j)} = α_{j−1}: (2.12) D = 2Σ_{a=1}^r w_a(w_a − w_{a+1})^T (w_0 = w_{r+1} = 0). bcz2004b used the anticlockwise one (transpose).
  Defect conditions (2.14): ∂xφ − E∂tφ − D∂tψ + ∂_φB = 0, ∂xψ − D^T∂tφ + E∂tψ − ∂_ψB = 0.
- Constant configurations: roots shifts (2.15)-(2.16) all have (2.17) (E0, P0) = −(2hm/β²)(cosh η, −sinh η), σ = e^{−η}; weight-shifted ones (2.18) have complex (E,P) = −(2hm/β²)(cosh(η + 2aπi/h), −sinh(η + 2aπi/h)) — on the same mass shell.
  (½)α_j·Dw_p = a/h (j ≥ p) or −(h−p)/h (j < p).
- Delay (2.20): z_a = (e^{−(θ−η)} + i e^{−iγ_a})/(e^{−(θ−η)} + i e^{iγ_a}), γ_a = πa/h; self-conjugate a = h/2: z = coth((θ−η)/2) (2.21); z_ā = conj z_a; arg z_a ∈ (−2γ_a, 0) (2.22) → selection rule: soliton charge can only move to one adjacent weight. Pole (2.23) θ = η + (iπ/2)(1 − 2a/h): no classical absorption except a = h/2. Classical bootstrap (2.25) z_a(θ − iŪ)z_b(θ + iŪ) = z_c̄(θ).
- Quantum: Hollowood S-matrix S^{ab} = ρ^{ab}R^{ab} (3.1)-(3.2); S^{11} entries (3.5) with (3.6) x_j = e^{hγθ_j/2}, q = −e^{−iπγ}, γ = 4π/β² − 1; ρ^{11} (3.7)-(3.8). Weights of V_1 (3.3) l_j = Σ_{l}(h−l)α_l/h − Σ_{l<j}α_l.
- Functional-integral phase (4.1)-(4.10): T(r,s) = e^{iτ(r,s)}T(0,0), (4.6) Q ≡ e^{4π²i/β²} = q^{−1} (consistent with q = −e^{−iπγ}, γ = 4π/β² − 1); bootstrap determines E uniquely (4.21), agreeing with the classical E.
- a2 result (triangle eq. (5.1)) after gauge fixing to one parameter (5.8)-(5.11): T^1 = g^1(θ) [[Q^{α·l1}δ, x̂²δ_{β−α1}, x̂Q^{−α·l2}δ_{β+α0}],[x̂Q^{−α·l3}δ_{β+α1}, Q^{α·l2}δ, x̂²δ_{β−α2}],[x̂²δ_{β−α0}, x̂Q^{−α·l1}δ_{β+α2}, Q^{α·l3}δ]] (5.9);
  antisoliton T^2 (5.10) has zeros (classical selection rule); x̂ = e^{γ(θ−Δ)}, g² = g¹(θ−iπ/3)g¹(θ+iπ/3)(1 + x̂³) (5.11).
- Inverse/crossing (6.1)-(6.6): T̃(−θ) = T(θ)^{-1}, T²(θ) = T̃¹(iπ − θ) (6.5) ⇒ g¹(θ)g¹(θ+2πi/3)g¹(θ+4πi/3)(1 − Q³x̂³)(1 − Qx̂³) = 1 (6.6); minimal solution (6.7)-(6.9) Gamma products with z = 3iγ(θ − Δ)/(2π).
- Unstable bound state: pole of g¹ at θ_P = Δ − iπ/3 − iπ/(3γ); identification (6.12) Δ = η + iπ/2; energy (6.13) E = m_s cosh θ_P = m_s cosh η sin(π/3 + π/3γ) + i m_s sinh η cos(π/3 + π/3γ); Re E > 0 and Im E < 0 iff 1/2 < γ < 2, i.e. 4π/3 < β² < 8π/3 — window excludes the classical limit (β→0). This is the "unstable bound state window".
  SG comparison: T^{SG} = g^{SG}[[Q^{α/2}δ, (−q)^{1/2}e^{γ(θ−η)}δ_{β−2}],[(−q)^{1/2}e^{γ(θ−η)}δ_{β+2}, Q^{−α/2}δ]], x̂ = (−q)^{−1/2}e^{γ(θ−η)} (6.14).
- Lightest breathers (7.1)-(7.5): poles θ_k = iπ(1 − 2k/3γ), m_k = 2M sin(πk/3γ); T_1^b(θ) = e^{−iπ/3} sinh((θ−η)/2 − iπ/6)/sinh((θ−η)/2 + iπ/6), T_2^b(θ) = e^{iπ/3} sinh((θ−η)/2 + 2iπ/3)/sinh((θ−η)/2 + iπ/3) — coupling independent.
  Linearised check (7.6)-(7.11): ρ_k = −Σ_l α_l ω^{lk} (7.10); R_k = 0.
- Misprint (verified numerically, m = 1): (7.11) prints T_k = (i e^{−θ} m_k − σ(1 − ω^{−k}))/(i e^{θ} m_k + σ(1 − ω^k)); the version that is unitary, depends only on θ−η, and reproduces (7.5) for k = 1, 2 has e^{−θ} in the DENOMINATOR as well: T_k = (i e^{−θ} m_k − σ(1 − ω^{−k}))/(i e^{−θ} m_k + σ(1 − ω^k)).
- Sec. 8: anticlockwise permutation swaps which of T¹/T² carries the zeros; alternative "setting" (8.1)-(8.2) D + D^T = −2, E = 1 + D^T (spin −2 charges).

## cz2009a — Corrigan, Zambon, "Comments on defects in the a_r Toda field theories", 0902.1307 (J. Phys. A 42 (2009) 304008)
- Bulk (2.1): L_u = ½∂u·∂u − (m²/β²)Σ_{j=0}^r n_j(e^{βα_j·u} − 1); conformal Toda = drop j=0. Defect Lagrangian (2.2) L_D = θ(−x)L_u + θ(x)L_v − δ(x)W (NB minus sign and W, not B).
- Energy–momentum argument (Sec. 2), the cleanest derivation for the book: assume at x=0 (2.5) u_x = Au_t + Bv_t + X, v_x = Cu_t + Dv_t + Y with constant matrices.
  Ė total derivative ⇒ (2.6) C = B^T, A, D antisymmetric, X = −∇_u𝒟, Y = ∇_v𝒟 (𝒟 = defect potential; E + 𝒟 conserved).
  Ṗ total derivative ⇒ (2.7) 1 − A² = BB^T, 1 − D² = B^TB, AB + BD = 0; (2.8) (∇_u𝒟)² − (∇_v𝒟)² = 2(U − V); (2.9) existence of Ω with P + Ω conserved.
  ⇒ B = (1 − A)O, O orthogonal; after rotating v: (2.12)-(2.15) W = ½u_t·Au − ½v_t·D'v + u_t·(1 − A)v + 𝒟, D' = −A, B = 1 − A, B^T = 2 − B.
  Exponential ansatz ⇒ (2.21) 𝒟 = Σ_k p_k e^{(u−v)·y_k/2} + Σ_l q_l e^{(u·(1−A) + v·(1+A))x_l/2}; (2.22)-(2.35) case analysis ⇒ every Dynkin node has ≤ 2 neighbours and equal root lengths ⇒ only a_r (assuming constant A, B). This is a no-go for type-I defects outside a_r.
- a_r solution (2.36)-(2.41): B = 1 − A = 2Σ_{a=1}^r (w_a − w_{a+1})w_a^T (2.38) (or its transpose (2.39)); (2.41) 𝒟 = Σ_{k=0}^r σ^{-1}e^{(u−v)·(1+A)α_k/2} + Σ_{l=0}^r σ e^{(u·(1−A) + v·(1+A))α_l/2}.
  Conformal variants (2.42)-(2.44): a defect can peel off one node: a_r conformal Toda | a_{r−1} conformal Toda + free massless field; r defects reduce a_r Toda to free fields (or Liouvilles). Lax pair (3.1)-(3.13): K̃ = I + (σ/λ)ΣE_{−α_i} in the defining rep; positive-λ alternative (3.14)-(3.17) with B = −(1 + A).
- Sec. 4 (complex a_r): (4.1) L_D = θ(−x)L_u + θ(x)L_v − δ(x)[½u_t·Au + ½v_t·Av + u_t·Bv + 𝒟]; (4.2) α_k·Bα_j = 2 (k=j), −2 (k=j+1), 0. Solitons (4.5) u_a = −(1/iβ)Σ_j α_j ln(1 + E_a ω_a^j), E_a = e^{a_a x − b_a t + ξ_a}, ω_a = e^{2πia/h}; M_a = 2hm_a/β² (4.6).
  Static defect energies (4.8) (E,P)_{a,b} = −(2hm/β²)(cosh(η − 2(a−b)πi/h), −sinh(η − 2(a−b)πi/h)), σ = e^{−η}.
  Delay (4.10) z_a = (e^{−(θ−η)} − i e^{iπa/h})/(e^{−(θ−η)} − i e^{−iπa/h}); pole at θ = η + i(πa/h − π/2). [Compare cz2007 (2.20): different sign of i and a ↔ h−a, reflecting the transposed B and the opposite sign convention for u.]
- Quantum (Sec. 5): S^{11} (5.4)-(5.6) as in cz2007, q = −e^{−iπγ}, γ = 4π/β² − 1, 0 < β² < 4π. Triangle relations (5.1).
  Soliton (a=1) transmission matrix (5.8): T_{iα}^{iβ} = g(θ)Q^{α·l_i}δ_α^β, T_{iα}^{(i−1)β} = g(θ) x̂ δ_α^{β−l_i+l_{i−1}}, all else zero; Q = −e^{iπγ}, x̂ = e^{γ(θ−Δ)}. Respects the classical selection rule.
  T̃ (right-to-left) via (5.10) T(θ)T̃(−θ) = 1; crossing (5.9) T^{(h−a)}(θ) = T̃^a(iπ − θ); antisoliton matrix (5.11) T^r = x̂^j Q^{−α·k_{(i+j)}}/(g(θ−iπ)(1 − x̂^hQ^{−1})) — no zeros.
  Bootstrap (6.1)-(6.8): g(θ + 2πi) = g(θ)(1 + x̂^h(−Q)^r)/(1 + x̂^h(−Q)^{r+2}); minimal (6.9)-(6.12) g^a = ĝ^a x̂^{−a/2}/(2π), ĝ^a = Γ[½ + (h−a)γ/2 − z] × Π_k Γ-ratios, z = ihγ(θ−Δ)/2π.
- Unstable bound states: pole (6.13) θ_a = Δ − iϑ_a, ϑ_a = π(h−a)/h + π/(hγ); Δ = η + iπ/2 (σ = e^{−η}); (6.14) E_a = M_a cosh η sin ϑ_a + i M_a sinh η cos ϑ_a; unstable iff (6.15) π/2 ≤ ϑ_a < π (and η > 0).
  "Anti-soliton" reps a > h/2: window a − h/2 ≤ 1/γ < a, excluding the classical limit; for a2 this is ½ ≤ 1/γ < 2, i.e. the cz2007 window ½ < γ < 2 (the interval is invariant under γ ↔ 1/γ).
- Possible misprint: for "soliton" reps a < h/2 the text says the window is 1/γ < h/2 − 1 (h even) or 1/γ < (h−1)/2 (h odd); from (6.13) directly, ϑ_a < π ⟺ 1/γ < a (and ϑ_a ≥ π/2 holds automatically). The two agree only for the largest such a. Re-derive before quoting.

## cz2009b — Corrigan, Zambon, "A new class of integrable defects", 0908.3126 (J. Phys. A 42 (2009) 475203)
- Introduces TYPE II defects: an extra defect degree of freedom λ(t). Lagrangian (2.1):
  L = θ(−x)L_u + θ(x)L_v + δ(x)[½(u v_t − v u_t) + λ(u−v)_t − λ_t(u−v) − D(u,v,λ)]. "Type I" = λ ≡ 0.
- Defect conditions (2.2)-(2.4): u_x = v_t − 2λ_t − D_u, v_x = u_t − 2λ_t + D_v, u_t = v_t + ½D_λ. E + D conserved.
- Momentum condition: with p = (u+v)/2, q = (u−v)/2 at x=0: (2.10)-(2.11) ⇒ D = f + g, F = −f_q, Ω = f − g with f = f(q, p − λ), g = g(q, λ), and the key equation (2.12) ½(f_q g_λ − f_λ g_q) = U − V, equivalently (2.13) {D, Ω} = −2(U − V) where the bracket treats λ and π_λ = −2q as conjugate.
  (Type I condition for comparison (2.8): D_uu = D_vv, ½(D_u² − D_v²) = U − V.)
- Solutions (Sec. 3), field normalisation U(u) = e^u + e^{−u} etc. (β, m scaled out):
  sinh/sine-Gordon (3.2): f = σ^{-1}[2e^{p−λ} + e^{−p+λ}(e^q + e^{−q} + γ)], g = σ[e^λ(e^q + e^{−q} + γ) + 2e^{−λ}], γ = e^{2τ} + e^{−2τ} — two parameters (σ, τ). Defect conditions (3.3). Vacuum λ: e^{2λ} = 1/(2cosh²τ) (even jump) or 1/(2sinh²τ) (odd jump) (3.4).
  Liouville (3.5); Tzitzéica a2^(2) (3.6)-(3.7): f = σ^{-1}e^{2p−2λ} + e^{−p+λ}(e^q + e^{−q}), g = (σ/2)(8e^{−λ} + e^{2λ}(e^q + e^{−q})²) — first defects for a non-a_r model.
  Massive free (3.8)-(3.9): f = m((p−λ)²/β + αq²), g = m(λ²/α + βq²); after eliminating λ: L_D = δ(x)[(4αβ/(m(α+β))) q_t² − ½((α−β)/(α+β))(uv_t − vu_t) − m(p²/(α+β) + (α+β)q²)].
  Linearised transmission (3.11): T = (i(αe^θ − βe^{−θ}) + 1)/(i(αe^θ − βe^{−θ}) − 1), R = 0. α = β admits a zero-energy "bound state".
  SG quadratic limit: σ = e^η, α = σ/(2cosh τ), β = 1/(2σ cosh τ); Tzitzéica (m=√6): α = 1/(√6σ), β = 2σ/√6 (3.10).
- Solitons (Sec. 4): SG soliton e^{u/2} = (1+E)/(1−E), E = e^{ax+bt+c}, a = √2 cosh θ, b = −√2 sinh θ (u, v imaginary). Delay (4.3) z1 = tanh((θ−η+τ)/2) tanh((θ−η−τ)/2), z2 = 1/z1, σ = e^η — identical to two type-I defects with parameters η ± τ ("fused pair"); λ evolution (4.4)-(4.5).
  Tzitzéica soliton e^u = (1+E)²/(1 − 4E + E²), a = √6 cosh θ; three delays (4.9)-(4.11): z1 = coth((θ−η)/2 − iπ/12) coth((θ−η)/2 + iπ/12), z2 = z̄3 = coth((θ−η)/2 + iπ/4) tanh((θ−η)/2 − iπ/12), z1z2z3 = 1, e^η = √2 σ.
- Sec. 5: Tzitzéica Bäcklund transformation from folding two "squeezed" a2^(1) type-I defects (B and B^T) — (5.4)-(5.8); a2^(1) defect conditions (5.6) reproduce the Fordy–Gibbons Bäcklund transformation.
- Sec. 6: defects as second-class Hamiltonian constraints χ1 = π_u + ½v, χ2 = π_v − ½u (6.4) with Dirac brackets.
- Appendix A: conserved spin-3 charge with type II conditions.
- Note: here σ = e^{+η} (opposite to bcz2005/cz2007 where σ = e^{−η}).

## cz2010 — Corrigan, Zambon, "A transmission matrix for a fused pair of integrable defects in the sine-Gordon model", 1006.0939 (J. Phys. A 43 (2010) 345201)
- Type II Lagrangian in improved form (2.1): L = θ(−x)L_u + θ(x)L_v + δ(x)(2qλ_t − D(λ,u,v)), q = (u−v)/2, p = (u+v)/2; U(u) = 2 − (e^u + e^{−u}) (fields imaginary for SG; β, m scaled out).
  (2.2) D = −√2σ(e^{(p+λ)/2}(e^{q/2−τ} + e^{−q/2+τ}) + e^{−(p+λ)/2}(e^{q/2+τ} + e^{−q/2−τ})) − (√2/σ)(e^{(p−λ)/2}(e^{q/2+τ} + e^{−q/2−τ}) + e^{−(p−λ)/2}(e^{q/2−τ} + e^{−q/2+τ})) = f(p+λ, q) + g(p−λ, q).
  Defect conditions (2.3): 2q_x = −D_p, 2p_x − 2λ_t = −D_q, 2q_t = −D_λ. Relation to cz2009b form via shift (2.4).
  Fusing: two type-I defects (u|λ|v) with σ1σ2 = σ² ... precisely √(σ1σ2) = σ ≡ e^{−η}, √(σ1/σ2) = e^{−τ}; but the λ bulk kinetic term is dropped, so a type II defect cannot be split (3 conditions vs 4).
- Energies (2.6)-(2.10): E = 8√2 cosh θ ± 8√2 cosh η cosh τ (even jump) or ± 8√2 sinh η sinh τ (odd jump).
- Linearised (2.11): T2 = −[sinh((θ−η+τ)/2 − iπ/4) sinh((θ−η−τ)/2 − iπ/4)]/[sinh((θ−η+τ)/2 + iπ/4) sinh((θ−η−τ)/2 + iπ/4)] = product of two type-I factors.
  Soliton-like setting: η + τ = ϑ ± iπ, η − τ = ϑ gives defect energy = soliton energy at rapidity ϑ, delay tanh²((ϑ−θ)/2), and (2.13) T3 = (sinh(θ−ϑ) + i)/(sinh(θ−ϑ) − i) (classical soliton–breather).
- Quantum (Sec. 3): SG S-matrix as in bcz2005 (3.1)-(3.2); KL type-I matrix (3.4)-(3.5) with ρ_I ∝ e^{−γ(θ−Λ)/2}; product solution T_{I−I} (3.6)-(3.8) for two separated defects (two pairs of labels).
  NEW solution (3.9): T = ρ(θ)[[(a+Q^α + a−Q^{−α}x²)δ_α^β, x(b+Q^α + b−Q^{−α})δ_α^{β−2}],[x(c+Q^α + c−Q^{−α})δ_α^{β+2}, (d+Q^αx² + d−Q^{−α})δ_α^β]], x = e^{γθ}, with a±d± = b±c±.
  Inverse formula; Δ(θ) (3.10); crossing ρ(θ)ρ(θ+iπ)Q²Δ(θ) = 1 (3.11); parametrisation (3.12) via η1, η2; ρ = ρ_{I−I}/√(a+d−) (3.13); unitarity (3.14) ⇒ η1, η2 real.
  Final type II matrix (3.15)-(3.16) with ρ_II = f_II(z1,z2)e^{−γ(θ−η1)/2}e^{−γ(θ−η2)/2}e^{iπγ/2}/(2π), z_p = iγ(θ−η_p)/π, Gamma-product f_II; resonance poles (3.17) θ = η_p + iπ/(2γ) [cf. bcz2005 pole at η − iπ/2γ].
  Lightest breather (3.19): ¹T = −Π_{p=1,2} sinh((θ−η_p)/2 − iπ/4)/sinh((θ−η_p)/2 + iπ/4) = classical (2.11) with η_{1,2} = η ± τ.
- Sec. 4 "defects as soliton constituents": choosing (4.1)-(4.2) and (4.5) η1 = ϑ − iπ − iπ/2γ, η2 = ϑ + iπ + iπ/2γ, the odd-label 4×4 block of (3.9) IS the SG soliton S-matrix with ρ̂ = ρ_S(θ − ϑ); the lightest-breather factor becomes (4.6) (sinh Θ + i cos(π/2γ))/(sinh Θ − i cos(π/2γ)) = soliton–breather S-matrix. Even labels with (4.7) give the breather–breather S-matrix (4.8) (sinh Θ − i sin(π/γ))/(sinh Θ + i sin(π/γ)).

## cz2011 — Corrigan, Zambon, "Integrable defects in affine Toda field theory and infinite dimensional representations of quantum groups", 1012.4186 (Nucl. Phys. B 848 (2011) 545)
- Main message: transmission matrices = intertwiners between a finite-dim rep of U_q(ĝ) and an infinite-dim rep of the BOREL subalgebra built from (generalised) q-oscillators; the linear intertwining equation replaces the quadratic triangle relation (credits Weston 1006.1555 for SG).
- Type II Lagrangian (2.1) (multi-component): δ(x)(2q·λ_t − D(λ,u,v)), D = f(p+λ, q) + g(p−λ, q), (2.2) ∇_λf·∇_qg − ∇_λg·∇_qf = U(u) − V(v).
  a2^(2) bulk U = −(m²/β²)(e^{iβu√2} + 2e^{−iβu/√2} − 3); defect potential (2.3); delay (2.4) z = coth((θ−η)/2 − iπ/12) coth((θ−η)/2 + iπ/12), σ = √2 e^{−η}.
  a_n^(1) type II, first given here: bulk (2.5) U = −(m²/β²)Σ(e^{iβα_j·u} − 1); setting A (2.6) D = (m/β²)Σ_j[σ e^{iβα_j·(p+λ)/2}A_j(q) + σ^{-1}e^{iβα_j·(p−λ)/2}A_{j+1}(q)], A_j(q) = γe^{iβα_j·Gq/2} + γ^{-1}e^{−iβα_j·Gq/2}, G = 2Σ(w_a − w_{a+1})w_a^T; setting B: p → −p, G → −G. Interpretation: fusion of a "clockwise" and an "anticlockwise" type-I defect.
- a2^(2) quadratic route (Sec. 3): Smirnov S-matrix from Izergin–Korepin R (3.2)-(3.4), x_i = q^{2πθ_i/ξ}, ξ = (2/3)πβ²/(8π − β²), q = e^{i4π²/β²} (3.6). General solution of triangle eq. (3.5)-(3.6) with four constants and one free function λ(α) (gauge freedom ↔ classical λ field). Crossing (3.7)-(3.8), bootstrap at Θ = 2πi/3 (3.9)-(3.10); reduced matrix (3.11) with F(α) = (q + q^{-1})(q^{2α−1} + q^{−2α+1}); ρ (3.12) Gamma products, z = i(θ−ν)/ξ; ν = η + iπ/2 (3.13).
  Unstable-state poles θ1 = ν − i(ξ + π/3), θ2 = ν − i(ξ + 2π/3); E1 = M cosh η cos(π/6 − ξ) + iM sinh η sin(π/6 − ξ), E2 = M cosh η cos(π/6 + ξ) − iM sinh η sin(π/6 + ξ); unstable for π/6 < ξ < 2π/3 (θ1) and −π/6 < ξ < π/3 (θ2) [as printed].
  Lightest breather (3.14): T = tanh((θ−ν)/2 − iπ/12) tanh((θ−ν)/2 + iπ/12) = classical (3.15) (i(e^{θ−η} − e^{−θ+η}) + √3)/(i(e^{θ−η} − e^{−θ+η}) − √3).
- a2^(2) linear route (Sec. 4): U_q(a2^(2)) with Cartan [[2,−1],[−4,2]], d0 = 4, d1 = 1 (4.1)-(4.3); COPRODUCT (4.4) Δ(K) = K⊗K, Δ(X±) = X±⊗K^{-1} + K⊗X±, Δ' = opposite. 3-dim rep from spin-1 of U_q(sl2) (4.5)-(4.6).
  Borel rep: X1+ = a†, X0+ = aa, K1 = κ1 q^N, K0 = κ0 q^{−2N}, a|j⟩ = F(j)|j−1⟩, a†|j⟩ = |j+1⟩, j ∈ Z (no ground state: j = topological charge). Serre ⇒ (4.8) F(N) = b1(−)^N + c1q^{−2N} + (b2(−)^N + c2q^{2N}), b1c2 = b2c1; choice (4.9) F = f1q^{−2N−1} + f2q^{2N+1}.
  Homogeneous gradation (4.10) E0 = xX0+. Intertwiner T'Δ(b) = Δ'(b)T' (4.11); solution (4.14)-(4.15); relation to quadratic solution T'^T(x,q) = −T(x^{-1}, q^{-1}) (4.16).
  Generalised rep (4.17) X1+ = a†μ(N), X0+ = aaλ(N): solution (4.19)-(4.24), D(N) = d0 + d1q^{−4N} + d2q^{4N}; free function μ(N) ↔ classical λ. Choosing F = ε1 + ε2q^{−2N} gives a†a − q²aa† = ε1(1 − q²): the q-oscillator algebra.
  SG remark: Borel rep X1+ = a†μ(N), X0+ = aλ(N), K1 = κ1q^N, K0 = κ0q^{−N} ⇒ D(N) = μ(N−1)λ(N)F(N) = f0 + f1q^{2N} + f2q^{−2N} (4.25); type I ⟺ f1 = f2 = 0 (Weston), type II ⟺ all nonzero.
  No finite-dim sub-rep for a2^(2) (cannot embed S in T), only semi-infinite truncations.
- a2^(1) (Sec. 5): U_q(a2^(1)), d_i = 2, 3-dim rep (5.2), spin gradation E_i = x^{2/3}X_i^+ (5.3); Jimbo R' (5.5), q² = e^{−4π²i/β²}. Borel rep (5.7) K_i = κ_i q^{N_i − N_{i+1}}, X_i^+ = a_i†a_{i+1} (three oscillators, cyclic-symmetric; differs from BHK-type reps), F_i(N_i) = (f_i^+)²q^{2N_i} − (f_i^−)²q^{−2N_i} (5.8)-(5.9). Solution (5.10)-(5.11).
  f_i^+ = 0 (all i) recovers the type-I matrix of cz2007 (B.24) = (5.12)-(5.13) with ε = −(1 − q⁴)(f1^−f2^−f3^−)^{2/3}z^{2/3}; f_i^− = 0 recovers (B.36); all nonzero = type II setting B. Truncation f_i^− = f_i^+q^{−8/3} embeds the soliton R' in T' (5.14ff); f_i^− = −f_i^+q^{2/3} embeds the antisoliton–soliton S.
- Appendix: universal recipe for q-oscillator Borel reps for every affine Toda algebra (a1^(1) naturally first of c_n^(1) series).

## cz2012 — Corrigan, Zambon, "Infinite dimension reflection matrices in the sine-Gordon model with a boundary", 1202.6016 (JHEP 06 (2012) 050)
- Keeps boundary topological charge labels: R_{aα}^{bβ}(θ) with a + α = b + β (2.1)-(2.2); reflection equation with labels (2.3); SG S-matrix (2.4)-(2.5).
- General shape (2.6); GZ embedded as α-independent solution (2.7); general class (2.8).
- Dressing a boundary by a defect (2.9): R̃_{aα α̃}^{bβ β̃}(θ) = T_{ac}^{α̃γ̃}(θ) R_{cd}^{αβ}(θ) T̂_{db}^{γ̃β̃}(θ) solves the reflection equation if T solves the triangle relation (2.10) — the quantum version of Habibullin's construction (Dirichlet + two Bäcklund = GZ boundary condition u_x(0) = κ sin((u+u0)/2) + κ^{-1} sin((u−u0)/2)). Bajnok–Simon noticed the same for sinh-Gordon.
- Most general SG transmission matrix (2.12)-(2.14): T = ρ[[(a+Q^αx^{-1} + a−Q^{−α}x)δ, μ(α)δ_{β−2}],[λ(α)δ_{β+2}, (d+Q^αx + d−Q^{−α}x^{-1})δ]] with (2.13) μ(α−2)λ(α) = a−d−q^{α−1} + a+d+q^{−α+1} + γ (γ ≠ 0 ↔ type II; KL/type I is a−=d+=c−=b+=0); inverse (2.15).
- Type I defect + GZ ⇒ GZ with shifted parameters (2.18); type II defect ⇒ genuinely α-dependent new reflection matrices (2.19).
- Proposed boundary Lagrangian with a boundary degree of freedom (2.20)-(2.23): L_B = θ(−x)L_u + δ(x)(uλ_t − B(u,λ)), B = e^{λ/2}f(u) + e^{−λ/2}g(u), f(u)g(u) = h+e^{u/2} + h−e^{−u/2} + 2(e^u + e^{−u}) + h0; e.g. f = f0 + 2(be^{u/2} + b^{-1}e^{−u/2}), g = g0 + 2(b^{-1}e^{u/2} + be^{−u/2}).
- Remark: Neumann boundary can store unlimited topological charge; Dirichlet stores none.

## delfino_mussardo_simonetti1994 — Delfino, Mussardo, Simonetti, "Scattering theory and correlation functions in statistical models with a line of defect", hep-th/9409076 (Nucl. Phys. B 432 (1994) 518)  [PLB companion hep-th/9403049 not read]
- Defect operator D as a zero-rapidity "particle" with ZF-type algebra (2.7): A_i†(β)D = R_ij(β)A_j†(−β)D + T_ij(β)DA_j†(β); DA_i†(β) = R_ij(−β)DA_j†(−β) + T_ij(−β)A_j†(β)D. (Defect assumed without internal d.o.f.; footnote says the generalisation is straightforward.)
- Unitarity (2.8): R(β)R(−β) + T(β)T(−β) = 1, R(β)T(−β) + T(β)R(−β) = 0. Crossing (2.9): R_ij(β) = S_{ij}^{kl}(2β)R_kl(iπ − β) [GZ-type, index placement as printed], T_{ij̄}(β) = T_ij(iπ − β) — the transmission crossing is ASSUMED ("natural"), then checked on examples.
- Reflection–Transmission (RT) equations (2.10)-(2.11) (four channels: RR = boundary YBE; TT = "triangle"/transmission YBE; TR and RT mixed).
- No-go (2.12): for a DIAGONAL bulk S-matrix with both R ≠ 0 and T ≠ 0 the mixed equations force S_ab(β_a + β_b) = S_ab(β_b − β_a) and S_ab(β_a+β_b)S_ab(β_a−β_b) = 1 ⇒ S constant = ±1 (free fermion/boson). Purely reflecting (boundary) theories are exempt.
- Ising (Sec. 3): L_B = Ψ̄(iγ^μ∂_μ − m)Ψ, Majorana γ^0 = σ2, γ^1 = −iσ1; defect L_D = −g δ(x)Ψ̄Ψ (3.3) (energy-density line, marginal). Exact (3.5): T = cos χ sinh β/(sinh β − i sin χ), R = i sin χ cosh β/(sinh β − i sin χ), sin χ = −g/(1 + g²/4). Duality g → 4/g: T → −T, R → R (3.12); self-dual g = ±2 purely reflecting (3.13). Lattice: sin χ = tanh 2(J − J̃) (3.17).
- Free boson (Sec. 4), δ-defect (4.6): continuity + jump ∂x ϕ+ − ∂x ϕ− = (g/4m)(ϕ+ + ϕ−) [as printed]; (4.8) T = sinh β/(sinh β + ig/4m), R = −(ig/4m)/(sinh β + ig/4m). g > 4m: resonance pair at Im β = −π/2; −4m < g < 0: defect bound state; g < −4m: poles in physical strip ⇒ instability.
- Later sections: multi-defect band structure; correlation functions via form factors + defect-operator matrix elements.
- The no-go is the starting point for the "purely transmitting" programme (bcz2005, Konik–LeClair) and for the CAFG refinement.

## konik_leclair1999 — Konik, LeClair, "Purely transmitting defect field theories", hep-th/9703085 (Nucl. Phys. B 538 (1999) 587)
- Defect as perturbation S = S_bulk + g∫dt O(0,t) (1.1); Euclidean conventions z = t + ix. Claim: purely transmitting integrable defects arise from CHIRAL defect operators O = φ_L (2.7), where the bulk perturbation is Φ = φ_Lφ_R.
- Conservation mechanism (2.2)-(2.3): J_x^{(+)}(0) − J_x^{(−)}(0) = ∂_tΘ ⇒ Q = Θ + ∫J_t; first order (2.9) J_L^{(+)} − J_L^{(−)} = 4πig∂_tH_L, J_R continuous. Energy AND momentum conserved ⇒ no reflection.
- Folding for massless bulk (Sec. 3): defect CFT = tensor product of two conformal boundary conditions; chiral defect ↔ boundary theory; (3.5)-(3.7).
- Defect ZF algebra (4.2): D_α A_a^{(+)}(θ) = T_{αa}^{βb}(θ)A_b^{(−)}(θ)D_β, A_a^{(−)}(θ)D_α = T̃_{aα}^{bβ}(−θ)D_βA_b^{(+)}(θ); unitarity (4.3) T(θ)T̃(−θ) = 1; crossing (4.4) T̃_{aα}^{bβ}(θ) = T_{αb̄}^{βā}(iπ − θ) (bars = charge conjugation; check placement in PDF); defect Yang–Baxter (4.5a,b); hermiticity (4.7) T̃_{aα}^{bβ}(θ) = [T_{βb}^{αa}(−θ)]*.
  NB: KL place the defect operator to the LEFT of right-region particles, opposite to DMS ordering; bcz2005 rewrote everything in its own (3.9) form.
- Massless: L–R scattering B_{ab}^{dc} = S(−∞) (4.13); map defect ↔ massless boundary YBE (4.12) requires defect degrees of freedom α ∈ α0 + Z with a + α = b + β (5.14): T_{αa}^{βb} = q^{−αa/2}R_a^{−b}(θ).
- Chiral defect sine-Gordon (5.1): S = ∫[½(∂Φ)² + λ cos βΦ] + g∫dt cos βφ_L(0,t); massless folding gives BSG with boundary term cos(βΦ_B/2) (5.2). Scaling argument (5.5)-(5.8), β̂ = β/√(4π): charges of spin 2, 4 (and all right-moving ones) survive; at β̂² = 1 (free fermion) reflection appears.
- SG S-matrix (5.9)-(5.11): a = sin(γ(π + iθ))ρ, b = −sin(iγθ)ρ, c = sin(γπ)ρ, γ = 8π/β² − 1; massless L–R (5.12)-(5.13) B = δδ q^{−ab/2}, q = −e^{−iγπ}.
- Massive ansatz (5.16)-(5.18): T_{α±}^{α±} = q^{±α/2}P_±, T_{α,±}^{α±2,∓} = q^{±α/2}Q_±; DYB ⇒ (5.20) (a − qb)(a − q^{-1}b) = c² ⇒ q = −e^{±iπγ}; choose q = −e^{−iπγ}. General (5.21) P_± = μ_±e^{−γθ/2}f, Q_± = ν_±e^{γθ/2}f; unitarity (5.23), crossing (5.24).
  Minimal solution (5.25)-(5.26): P_± = e^{±iξ}e^{−(γθ+ν)/2}e^{iδ(θ)}/(2cosh((γθ + ν − iπ/2)/2)), Q_± = e^{iτ±}e^{(γθ+ν)/2}e^{iδ}/(2cosh(...)), e^{iδ} = e^{iπγ/4}Π_l R_l(γθ+ν)/R_l(−(ν+γθ)) with R_l Gamma ratios; CDD ambiguity (5.27); ν ~ unknown function of g²/λ, β² (later identified with the BCZ defect parameter, γη ↔ ν up to conventions).
- Breather bootstrap (5.30)-(5.37): T^b_{αα}(θ) = P_−(θ − iu)P_+(θ + iu) + (−1)^b q^{−1}Q_−Q_+ — α-independent; crossing and unitarity hold.
- Free-fermion point (Sec. 6): action (6.1) with fermionic defect modes a_±; massive case gives reflection + transmission (consistent with DMS, S = −1).
- Caveat flagged by the authors: integrability of the massive chiral-defect theory is not proven; the massive T may only be meaningful in the massless limit. (bcz2005 later re-interpreted the KL T as the Lagrangian jump-defect's.)

## cafg2002 — Castro-Alvaredo, Fring, Göhmann, "On the absence of simultaneous reflection and transmission in integrable impurity systems", hep-th/0201142 (unpublished preprint, v3 June 2002)
- ZF algebra with defect operator Z_α carrying internal states α, parity NOT assumed (1)-(3): Z_i(θ)Z_α = R_{iα}^{jβ}(θ)Z_j(−θ)Z_β + T_{iα}^{jβ}(θ)Z_βZ_j(θ); Z_αZ_i(θ) = R̃_{iα}^{jβ}(−θ)Z_βZ_j(−θ) + T̃_{iα}^{jβ}(−θ)Z_j(θ)Z_β (R̃, T̃ = from the right).
- Unitarity (4)-(6). Defect Yang–Baxter equations from three inequivalent orderings: (7)-(10) [RR, RT, TR, TT], their parity images, and (11)-(14). Tensor form (15)-(21); the TT equation (17) is the "triangle"/transmission YBE S(θ12)[T(θ2)⊗T(θ1)] = [T(θ1)⊗T(θ2)]S(θ12) (defect-index products implicit).
- NO-GO statement and its assumption: if (22) T_α^β(θ1)⊗R_β^γ(θ2) = T_β^γ(θ1)⊗R_α^β(θ2) holds (automatically true with NO defect degrees of freedom, or if either T or R is diagonal/abelian in the defect indices; or any of three alternatives in the footnote), then (16) at θ2 = 0 gives (23) S(θ)[I⊗R(θ)]S(θ) = R(θ)⊗I ⇒ S(θ) = Pσ is rapidity-INDEPENDENT (24). So: R ≠ 0 ≠ T requires a constant bulk S (free boson/fermion, Federbush/anyonic generalisations), unless the defect has internal states with R and T non-commuting in them.
  Generalises DMS: no parity or diagonality assumption in bulk; covers non-diagonal S.
- Bootstrap version (25)-(29): fusing a, b → c before/after the defect gives T_aT_b = T_c (27) and the mixed equations (28)-(29) forcing T = 0 or R = 0 whenever the bulk has bound states (excludes theories with fusing, but (24) also excludes sinh-Gordon).
- Excited impurity states (30)-(32): T_bβ(θ) = T_bα(θ) ⇒ no excited defect states when T ≠ 0 (in their framework — contrast with the unstable states of bcz2005/cz2007, which arise because the KL/BCZ defect carries topological charge labels with T non-diagonal in them).
- Free-fermion examples (33)-(39): D = gψ̄ψ: R = −i sin B cosh θ/(sinh θ + i sin B), T = cos B sinh θ/(sinh θ + i sin B), sin B = −4g/(4 + g²); D = gψ̄γ⁰ψ: R = −i sin B/sinh(θ + iB), T = sinh θ/sinh(θ + iB). Residue signs ⇒ no excited impurity states (disagrees with DMS on this point).
- Status: arXiv-only; no journal version found (INSPIRE lists none). Cite as preprint.

## bajnok_george2006 — Bajnok, George, "From defects to boundaries", hep-th/0404199 (Int. J. Mod. Phys. A 21 (2006) 1063)
- Folding: any relativistic defect theory (4) L = Θ(−x)[½(∂Φ_L)² − V_L] + Θ(x)[½(∂Ψ_R)² − W_R] − δ(x)U(Φ_L,Ψ_R) is equivalent (classically and in the adiabatic/LSZ quantum framework) to a boundary theory (8) of two non-interacting fields on x ≤ 0, obtained by the parity map P_Ψ: Ψ_R → Ψ_L. Impurity (continuous field) case (5)-(9): continuity + jump ∂xΨ − ∂xΦ = ∂U/∂φ0.
  Dictionary (10): boundary R-matrix [[R_ΦΦ, R_ΨΦ],[R_ΦΨ, R_ΨΨ]] ↔ defect RT-matrix [[R_L, T_R],[T_L, R_R]]; bulk S of folded theory (11) S_ΨΦ = S_ΦΨ = 1, S_ΦΦ = S_L, S_ΨΨ = S_R.
- Consequences: defect reduction formulae, perturbation theory, Coleman–Thun, Cutkosky rules carry over.
- Transmission YBE (12): T_L(θ1)T_L(θ2)S_R(θ2−θ1) = S_L(θ2−θ1)T_L(θ2)T_L(θ1) (indices in source) = special case (13) of the boundary YBE with species-changing reflections.
- Unitarity (14) (four equations). DEFECT CROSSING-UNITARITY DERIVED (not conjectured, unlike DMS) from Ghoshal–Zamolodchikov boundary cross-unitarity (15):
  R_L(θ) = S_L(2θ)R_L(iπ − θ) [indices as in GZ], T_L{}_{αa}^{βb}(θ) = T_R{}_{αb̄}^{βā}(iπ − θ), and L ↔ R. So the left-incoming transmission at θ crosses into the RIGHT-incoming transmission at iπ − θ with conjugated, exchanged particle labels.
- Defect operator from the GZ boundary state by unfolding: general D (unnumbered) and, for purely transmitting defects, using T_L(iπ/2 + θ) = T_R(iπ/2 − θ): (16) D = exp ∫_{−∞}^{∞} (dθ/2π) T_R(iπ/2 − θ) A_{Ψ_R}(θ) A^+_{Φ_L}(θ). Confirms Konik–LeClair's (underived) proposal.
- Purely transmitting solutions via boundary bootstrap (17)-(18): starting from trivial transmission plus CDD pole at θ = iu, the excited "defect" gives T_L = S_ΦΦ(θ − iu), T_R = S_ΨΨ(θ + iu) — the bulk S-matrix itself solves (12) (like a Lax operator); u = 0 is a standing particle. Dressing (19): R1 = T_L R0 T_R(−θ) solves the reflection equation (defect fused to a boundary).
- Defect TBA (Sec. 5): Z(L,R) = ⟨B_{n+1}|e^{−(L−x_n)H}D_n…D_1e^{−x_1H}|B_0⟩ (20); ground-state energy on a strip with n purely transmitting defects: transmission factors enter the TBA as log T_R(iπ/2 − θ) driving terms (rest of section).
- Restricts to parity-invariant ("strong") bulk theories and real scalar fields.

## bajnok_simon2008 — Bajnok, Simon, "Solving topological defects via fusion", 0712.4292 (Nucl. Phys. B 802 (2008) 307)
- Conventions: sinh-Gordon (1) L = ½(∂Φ)² − (m_cl²/b²) cosh bΦ; B = b²/(8π + b²) (8); S = (sinh θ − i sin Bπ)/(sinh θ + i sin Bπ) = −(−B)(1+B), block (x) = sinh(θ/2 + iπx/2)/sinh(θ/2 − iπx/2); mass (8); ε_bulk = m²/(8 sin πB) (9).
- Defect Lagrangian (3): L = Θ(−x)L(Φ−) − δ(x)D + Θ(x)L(Φ+), 2D = Φ+Φ̇− − Φ−Φ̇+ + M_cl e^μ cosh(b(Φ+ + Φ−)/2) + M_cl e^{−μ} cosh(b(Φ+ − Φ−)/2), M_cl = 4m_cl/b². (Φ− on the LEFT, Φ+ on the RIGHT; note the "−δ(x)" with the time-derivative term inside D, so the kinetic term sign is opposite to bcz's + ½(φ∂tψ − ψ∂tφ).)
- Classical fusion (Sec. 2): defect at 0 in front of a boundary at a: physics independent of a; a → 0 gives dressed boundary B' = D(Φ−, Φ+) + B(Φ+). Defect + Dirichlet(φ0) = general perturbed-Neumann (GZ) boundary (4) with M_cl cosh(μ ± bφ0/2) = M0 e^{∓bϕ0/2}; defect + Neumann = new time-dependent integrable BC (5).
- Quantum fusion (Sec. 3): defect unitarity and crossing (6) T+(−θ) = T−(θ)^{-1}, T−(θ) = T+(iπ − θ) [T− = transmission for particles from the left; T+ parametrised so its physical argument is positive]; fused reflection (7) R'(θ) = T+(θ)R(θ)T−(θ).
  Dirichlet R (10)-(12) and PN R (13)-(15) known ⇒ transmission factors (16): T−(θ) = −i sinh(θ/2 − iπ/4 + Bϑ/2)/sinh(θ/2 + iπ/4 + Bϑ/2), T+(θ) = i sinh(θ/2 − iπ/4 − Bϑ/2)/sinh(θ/2 + iπ/4 − Bϑ/2) — a pure CDD factor, no physical-strip poles.
  UV–IR relation (17) M e^{±b²ϑ/8π} = M_cl e^{±μ}; defect energy (18) ε_def = m cosh(Bϑ)/(2 sin Bπ). Neumann dressing (19)-(20): boundary CDD factors = defects in front of boundaries. Two defects with ϑ± = ±i(1 − π/2B)... reproduce a standing particle ("defect = half particle").
- Lee–Yang (Sec. 3.2): S = −(1/3)(2/3) [block notation], defect between I and Φ boundaries: R_b = T+ R_I T− ; defect bootstrap (21) T−(θ + iπ/3)T−(θ − iπ/3) = T−(θ) ⇒ (22) T−(θ) = [b+1][b−1], [x] = i sinh(θ/2 + iπx/12)/sinh(θ/2 + iπx/12 − iπ/2); T+ = [5−b][−5−b]; ε_def = m sin(bπ/6); b = 3 behaves as a standing particle. (Lagrangian realisation of the LY defect left open — done later in bhw2014.)
- Perturbation theory (Sec. 4): classical limit T− = −(m_cl sinh μ − iω)/(m_cl cosh μ − ik), = −i sinh(θ/2 − iπ/4 + μ/2)/sinh(θ/2 + iπ/4 + μ/2) — same form with Bϑ ↔ μ. One loop: Bϑ = μ(1 − b²/8π + …) (23) — refutes the naive identification Bϑ = μ of bcz2005 (the "quantum renormalisation of the defect parameter"). All divergences absorbed into m_cl (bulk and defect renormalise consistently: no anomaly); finite δμ = −μb²/8π; no reflection at O(b²).
- Defect TBA (Sec. 5): defect operator D = exp∫(dθ/2π)T+(iπ/2 − θ)a+(θ)a(θ); ground-state DTBA on a circle (25) ε̃(θ) = mL cosh θ − log T+(iπ/2 − θ) − ∫φ ∗ log(1 + e^{−ε̃}); UV analysis reproduces (9), (18); LY Y-system Y(θ − iπ/3)Y(θ + iπ/3) = 1 + Y(θ).
- Defect bound states (Sec. 6, LY only): physical strip of T∓ is Im θ ∈ [0, π/2]; poles of T− at −iπ(b ± 5)/6, of T+ at iπ(b ± 1)/6; for b ∈ [−1, 2] bound state |1+⟩ at u = π(b+1)/6 with energy m cos(π(b+1)/6), T−^{|1+⟩} = [b+1][b+3]; fundamental range b ∈ [−3, 2]; |2+⟩ for b ∈ [1, 2]; defect Coleman–Thun diagrams explain the remaining poles (bootstrap closes).

## bhw2014 — Bajnok, Holló, Watts, "Defect scaling Lee–Yang model from the perturbed DCFT point of view", 1307.4536 (Nucl. Phys. B 886 (2014) 93)
- Solves the conformal defect Lee–Yang model M(2,5) (two topological defects, structure constants from sewing relations), then perturbs.
- Integrability/topologicality criterion (Sec. 3): jumps of T and T̄ across the defect must be total derivatives (12)-(16), (24). Chiral perturbation S = S_DCFT − μ∫ϕ(y)dy (or anti-chiral μ̄ϕ̄ (31)-(33)) keeps the defect topological (pure defect flows between the two topological defects).
- Massive integrable defect (57)-(64): S = S_DCFT − λ∫Φ − ∫(μϕ + μ̄ϕ̄)dy; H = H− + H+ + 2π(1−h)(μϕ + μ̄ϕ̄) (61); P conserved iff λ = 2iμμ̄C[ϕ,ϕ̄] = 0.826608… μμ̄ (64): one physical parameter for fixed bulk mass. (Lagrangian/DCFT realisation of the Bajnok–Simon LY defect; analogue of the KL chiral-defect idea but with BOTH chiralities.)
- IR data recalled (90)-(96): S = −(1)(2) blocks, (91) T− = [b+1][b−1] = S(θ − iπ(3−b)/6), T+ = [5−b][−5−b] = S(θ + iπ(3−b)/6) — the defect scatters like a particle of rapidity iπ(3−b)/6; pure phase for b = ∓3 + iα; ε_Def = m sin(bπ/6) (93); Lüscher (94); DTBA (95)-(96).
- UV–IR relation by fusing to the boundary (97)-(100): μ̂ = (|h_c|m^{6/5}/2)e^{±i(b+3)π/5}, μ̄ = (|h_c|m^{6/5}/2)e^{∓i(b+3)π/5}, h_c = −0.685289…, consistent with (64).
- Defect TCSA numerics confirm transmission factors, bound states and defect energies; real spectrum only for real couplings (PT symmetry discussion).

## weston2010 — Weston, "An algebraic setting for defects in the XXZ and sine-Gordon models", 1006.1555 [math-ph] (unpublished preprint; v2 28 Jun 2010; no journal version found on INSPIRE/Crossref)
- Idea: soliton–defect transmission matrix = Borel-subalgebra intertwiner L between an infinite-dim q-oscillator rep W^{(r)} and the spin-½ evaluation rep V; the quadratic "STT = TTS" becomes the RLL relation, automatically satisfied.
- Algebra (2.1)-(2.5): U'_q(ŝl2), generators e_i, f_i, t_i^{±1}, i = 0,1. COPRODUCT (2.6): Δ(e_i) = e_i⊗1 + t_i⊗e_i, Δ(f_i) = f_i⊗t_i^{-1} + 1⊗f_i, Δ(t_i) = t_i⊗t_i. Δ' = flip. Borel U_q(b+) = ⟨e_i, t_i^{±1}⟩.
- Generalised oscillator algebra O(r1,r2) (Sec. 2.1): q^D a* q^{−D} = q a*, q^D a q^{−D} = q^{-1}a, aa* = (r1 + q^{−2D})(r2 + q^{2D}), a*a = (r1 + q^{2−2D})(r2 + q^{2D−2}); r1 = 0 or r2 = 0 gives the q-oscillator a*a − q²aa* = 1 − q² or aa* − q²a*a = 1 − q². Module W^{(r1,r2)} = ⊕_{j∈Z}C|j⟩: a|j⟩ = |j−1⟩, a*|j⟩ = (r1 + q^{−2j})(r2 + q^{2j})|j+1⟩.
- Borel rep W_ζ^{(r)}, r = (r0,r1,r2) (Sec. 2.2): e0 = ζa*/(q − q^{-1}), e1 = ζa/(q − q^{-1}), t1 = r0q^{−2D}, t0 = r0^{-1}q^{2D} (reparametrisation of Rossi–Weston 2002, the most general solution of (2.2),(2.4)). Truncations: r1 = −q^{−2n} or r2 = −q^{2n} ⇒ semi-infinite; r = (q^n, −q^{−2n}, −q^{−2}) or (q^n, −q², −q^{2n}) ⇒ finite, isomorphic to spin-n/2 evaluation modules V^{(n)}_{ζq^{∓(n+1)/2}}; root of unity truncation.
- Evaluation rep V_ζ (2.3) and R-matrix (2.7) R(ζ) = κ(ζ)^{-1}R̄(ζ), R̄ = [[1,0,0,0],[0,(1−ζ²)q/(1−q²ζ²),(1−q²)ζ/(1−q²ζ²),0],[0,(1−q²)ζ/(1−q²ζ²),(1−ζ²)q/(1−q²ζ²),0],[0,0,0,1]], with R∘Δ = Δ'∘R.
- THE L-OPERATOR (2.8)-(2.9): L(ζ)∘Δ(x) = Δ'(x)∘L(ζ) for x ∈ U_q(b+) gives
  L^{(r)}(ζ, q) = [[1 + r2ζ²q^{2−2D}, −ζa*],[−ζa, 1 + r1ζ²q^{2D}]] · diag(q^D, r0q^{−D}).
  RLL relation (2.10); inverse with det ∝ (1 − ζ²)(1 − ζ²r1r2); crossing L^{-1}(−ζq) ∝ (σ^x L(ζ) σ^x)^{t2}/(r0(1 − ζ²q²)(1 − ζ²q²r1r2)).
  Q-operator = regularised trace of L over W (Baxter TQ from the short exact sequence 0 → W_{ζq}^{(qr0,q^{−2}r1,r2)} → W_ζ^{(r)}⊗V_ζ → W_{ζq^{-1}}^{(q^{-1}r0,q²r1,r2)} → 0 and fusion relations).
- Defect–defect intertwiner R^{(r)} (2.11)-(2.13): only for r1 = 0 or r2 = 0, formal series R = P h(ζ, r0u)ζ^{D1+D2}, u = a1*q^{−2D1}a2, h(ζ,v) = (−vζ^{-1}; q²)_∞/(−vζ; q²)_∞ (q-binomial); for r1 = r2 = 0 two independent solutions — candidate defect–defect S-matrices (cf. bcz2005 U-matrix).
- XXZ with a defect (Sec. 3): defect line with weight L in the 6-vertex model; transfer matrix T(ζ) = Tr_0[R_{1,0}…L_{j,0}(ζ−1)…R_{N,0}] (spectral parameter shift as printed); Hamiltonian H = Σ_{i≠j−1,j} h_{i,i+1} + (conjugated h_{j−1,j+1}) + ((q − q^{-1})/2q)[[0, r0a*q^{−2D}],[r0^{-1}aq^{2D}, 0]]_{j,j+1}, h = −½(σ^xσ^x + σ^yσ^y + Δσ^zσ^z).
- Dictionary with CZ (Sec. 4), using cz2010 conventions: S(θ,q) ∝ P R̄(x,q), x = e^{γθ}, q = −e^{−iπγ} [Weston writes γ as a ratio of 4π − β² and β² (orientation not recoverable from the text layer); either way this differs from the KL/bcz2005 normalisation γ = 8π/β² − 1 used by cz2010 — check the PDF before quoting]; type I (4.15): T_I(θ,η)/ρ_I = ν^{−1/2}U_I(ν) L^{(r0=ν, r1=0, r2=0)}(ie^{γ(θ−η)}, q) U_I^{-1}; type II (4.16): T_II/ρ_II = U_II L^{(r0=1, r1=−b−/b+, r2=q²b−/b+)}(ie^{γθ}|b+|, q) U_II^{-1}; topological charge α = −2j; transposed index conventions noted.
  So: type I ⟺ plain q-oscillator (r1 = r2 = 0); type II ⟺ generalised oscillator (both r_i ≠ 0). Finite truncation W^{(q, −q^{±2}, −q^{±2})} ≅ V_{ζq^{±1}} explains the soliton S-matrix inside the type II T (cz2010 Sec. 4).
- Credits Bajnok, "Equivalences between spin models induced by defects" (JSTAT 2006, hep-th/0601107) for the earlier link between Rossi–Weston L and type I defects.

## penati_polvara2019 — Penati, Polvara, "Quantum anomalies in A_r^(1) Toda theories with defects", 1902.10690 (JHEP 06 (2019) 062)
- EUCLIDEAN (x0, x1), defect along x1 = 0; φ^{(−)} on x1 < 0, φ^{(+)} on x1 > 0. Bulk V = Σ_{j=0}^r e^{α_j·φ} (3.1) (β scaled into the action as overall 1/β² — massless perturbation theory with SV, Sd ∝ β^{-2} normal-ordered potentials, eqs. around (2.20)).
- Defect potential (3.3): V^{(d)} = (i/2)φ^{(−)}E∂0φ^{(−)} + iφ^{(−)}D∂0φ^{(+)} + (i/2)φ^{(+)}F∂0φ^{(+)} + B(φ^{(−)},φ^{(+)}); defect conditions (3.4) ∂1φ^{(−)} = iE∂0φ^{(−)} + iD∂0φ^{(+)} + ∂B/∂φ^{(−)}, ∂1φ^{(+)} = −iF∂0φ^{(+)} + iD^T∂0φ^{(−)} − ∂B/∂φ^{(+)} (the i's are the Euclidean continuation of BCZ's (3.2)).
- sinh-Gordon (r=1): V^{(±)} = e^{√2φ} + e^{−√2φ} (3.5), D = 1, E = F = 0, (3.6) B = σ^{-1}(e^{(φ^+ + φ^−)/√2} + e^{−(φ^+ + φ^−)/√2}) + σ(e^{(φ^+ − φ^−)/√2} + e^{−(φ^+ − φ^−)/√2}); spin-4 current (3.7); classical defect contributions Σ_E, Σ_P to the spin-4 charges (3.8)-(3.9).
  a_r, r > 1: (3.10) E = F = 1 − D = Σ(λ_{i+1}λ_i^T − λ_iλ_{i+1}^T); (3.11) B = Σ_i[σ^{-1}e^{½α_i·R} + σe^{½α_i·S}], R = Dφ^{(+)} + D^Tφ^{(−)}, S = D(φ^{(−)} − φ^{(+)}) (3.12); constraints (3.13)-(3.15); orthogonal-rotation family of solutions (3.16)-(3.17); spin ±2 alternative E = F = −1 − D; degenerate choice D = E = F = 0 = two decoupled boundary theories (3.18).
- Results (all orders in β via massless perturbation theory with Wick contractions):
  * Energy is conserved with an exact defect energy contribution, no renormalisation needed (Sec. 4).
  * Momentum conserved only after FINITE RENORMALISATION of the defect potential (4.12): B = √(1 + γ/2) B^{(0)}, with γ ≡ β²/2π (defined just before (4.4); action S0 = (1/2β²)∫(∂φ)², interactions (1/β²)∫:V:), because the trace term acquires a factor (1 + γ/2).
  * Higher spin: spin-4 (sinh-Gordon, Sec. 5.1) and spin-3 (a_r, r>1, Sec. 5.2) currents acquire ANOMALIES that no finite renormalisation of B can cancel ⇒ "doubtful quantum integrability" of the Lagrangian type-I defect. Exceptions: decoupled boundaries, and sinh-Gordon with the extra constraint ∂0φ^{(−)} = ∂0φ^{(+)} at the defect (5.26), which needs B = (1 + γ/4 − 11γ²/96 + …)B^{(0)} (5.28) — DIFFERENT from the momentum renormalisation (4.12): "momentum conservation is not a probe of quantum integrability".
  * Same conclusion claimed for sine-Gordon (β → iβ); suggested resolution: the BCZ/KL transmission matrices describe a theory whose quantum Lagrangian is not the naive type-I one; type II (λ as Lagrange multiplier) noted as possible rescue.
- Tension to flag in the book: contrasts with bajnok_simon2008's one-loop finding of "no anomaly" (only parameter renormalisation) for the transmission factor; Penati–Polvara examine higher-spin current conservation, Bajnok–Simon only the 2-point function at O(b²).

## blz1997 — Bazhanov, Lukyanov, Zamolodchikov, "Integrable structure of conformal field theory II. Q-operator and DDV equation", hep-th/9604044 (Commun. Math. Phys. 190 (1997) 247)
- Conventions: (1.2) q = e^{iπβ²}, (1.3) c = 13 − 6(β² + β^{−2}); free chiral boson (2.1) ϕ(u) = iQ + iPu + Σ_{n≠0} a_{−n}e^{inu}/n, (2.2) [Q,P] = iβ²/2, [a_n, a_m] = (n/2)β²δ_{n+m,0}; T(u) (2.5); Δ = (p/β)² + (c−1)/24 (2.7); vertex operators V±(u) = :e^{±2ϕ(u)}: (2.9).
- U_q(sl2) (2.10) [H,E] = 2E, [H,F] = −2F, [E,F] = (q^H − q^{−H})/(q − q^{-1}).
- Transfer matrices (2.11): T_j(λ) = tr_{π_j}[e^{2πiPH_j} P exp(λ∫_0^{2π}K_j(u)du)], K_j(u) = V−(u)q^{H_j/2}E_j + V+(u)q^{−H_j/2}F_j (BLZ "L-operator" = this path-ordered exponential, the CFT monodromy matrix); fusion (2.13) T(λ)T_j(q^{j+½}λ) = T_{j−½}(q^{j+1}λ) + T_{j+½}(q^jλ); T(λ) = 2cos 2πP + Σλ^{2n}G_n (2.15) with nonlocal IM G_n (2.16); semiclassical domain 0 < β² < ½ (2.18).
- q-OSCILLATOR ALGEBRA (2.27): q E+E− − q^{−1}E−E+ = 1/(q − q^{−1}), [H, E±] = ±2E±; any rep ρ with Z(p) = tr_ρ e^{2πipH} convergent for Im p > 0.
  A±(λ) = Z^{−1}(±P) tr_ρ[e^{±2πiPH} P exp(λ∫_0^{2π}du(V−(u)q^{±H/2}E± + V+(u)q^{∓H/2}E∓))] (2.29) — independent of ρ (coefficients fixed by (2.27) and cyclicity (2.32)); A− from A+ by P → −P, ϕ → −ϕ (2.30); series (2.31)-(2.33).
  Q±(λ) = λ^{±2P/β²}A±(λ) (2.38). Properties: commute with each other and with T_j (2.40); Baxter TQ (2.41) T(λ)Q±(λ) = Q±(qλ) + Q±(q^{-1}λ); quantum Wronskian (2.42) Q+(q^{½}λ)Q−(q^{−½}λ) − Q+(q^{−½}λ)Q−(q^{½}λ) = 2i sin 2πP; (2.43) 2i sin(2πP)T_j(λ) = Q+(q^{j+½}λ)Q−(q^{−j−½}λ) − Q+(q^{−j−½}λ)Q−(q^{j+½}λ).
  Mechanism (2.44)-(2.45): T_j^+ (Verma, highest weight 2j) = T_j + T_{−j−1}^+ and 2i sin(2πP)T_j^+(λ) = Q+(q^{j+½}λ)Q−(q^{−j−½}λ), from the decomposition of the tensor product of two q-oscillator (Borel) reps — "the latter are also representations of the Borel subalgebra of U_q(ŝl2)" (paraphrase). This is the q-oscillator ↔ Borel-representation link that Weston and CZ2011 exploit for defects.
- Sec. 3: vacuum eigenvalues of Q satisfy the Destri–de Vega (DDV) equation (3.8), solved for Δ → ∞; Sec. 4: spectral-determinant-type integral transform (4.1)-(4.5), T(λ) = Λ(q^{½}λ) + Λ^{−1}(q^{−½}λ) (4.4) [as printed: check exponents]; Sec. 5: boundary sine-Gordon / quantum Hall point-contact transport.
- Note: BLZ "III" (Yang–Baxter relation, where the q-oscillator L-operators L(λ) are written explicitly and RLL verified) is hep-th/9805008, CMP 200 (1999) 297 — NOT hep-th/9607099 (see Table).

## cck2024 — Copetti, Córdova, Komatsu, "Noninvertible symmetries, anomalies, and scattering amplitudes", 2403.04835 (Phys. Rev. Lett. 133 (2024) 181601)
- Setting: 2d theories with a fusion-category symmetry C of topological lines L_a (fusion (1), junction spaces V_ab^c, completeness (3)-(4), quantum dimensions ⟨L_a⟩ = d_a (5), F-symbols/tetrahedral symbols (6)-(7); Z2 anomaly ⇔ F-symbol ε = −1 (8)). Gapped IR vacua form a module category n_{ai}^j (9)-(11); regular representation L_a|0⟩ = |a⟩ (12) for ϕ1,3-flows.
- Kink S-matrices S^{ab}_{dc}(θ) (13) (vacua a,b,c,d around the 2→2 diagram); standard crossing (14) S^{ab}_{dc}(θ) = S^{bc}_{ad}(iπ − θ); unitarity (15); YBE (16). [Typo: s = 4m⁴cosh²(θ/2) should read 4m²cosh²(θ/2).]
- Literature RSOS S-matrices for M_n + ϕ1,3 (λ<0, n−1 vacua, d_a = sin((2a+1)π/n)/sin(π/n)) (18) are crossing symmetric but violate the Ward identity (19) of the non-invertible symmetry, with the action of lines on kinks (20) ∝ (d_ad_a'd_bd_b')^{1/4} × tetrahedral symbol.
- Proposal (21): S^{ab}_{dc}(θ) = Z(θ)√(d_ad_c/(d_bd_d)) [sinh(θ/n)δ_bd + sinh((iπ−θ)/n)δ_ac]^{...}(i.e. extra factor (d_ad_c/(d_bd_d))^{iθ/2π} relative to (18)), which satisfies unitarity, YBE and the Ward identity but obeys the MODIFIED CROSSING RULE (22):
  S^{ab}_{dc}(θ) = √(d_a d_c/(d_b d_d)) · S^{bc}_{ad}(iπ − θ).
  The extra factor drops out of TBA (footnote 6), so both versions reproduce known TBA. Same rule for tricritical Ising + ϕ2,1 (Fibonacci {1, W}, d_W = (1+√5)/2) (23), supporting the Zamolodchikov/Colomo et al. proposal [44,45] that violated standard crossing.
- Z2-anomalous example: SU(2)_1 WZW + JJ̄ (= compact boson with cos(2X/R) at R = √2) (24)-(29): kinks carry Kramers doublets (T² = −1 on worldline), S = S_RSOS2 ⊗ Š, Š(θ) = S0(θ)(θδδ − iπδδ)/(θ − iπ) (27)-(28); crossing violated by a SIGN (29) Š^{ab}_{dc}(θ) = −C Š(iπ − θ) C (C² = −1).
- Derivation (Sec. IV): normalised kink states on a line include TQFT factors: ⟨⟨a|a⟩⟩ = d_a (31), action of lines on vacua divided by √d_φ (34); S-matrix from Cauchy-slice correlators (37)-(41): numerator crossing-symmetric, denominator (in/out state normalisation networks) channel dependent ⇒ general rule (42) (ratio of line networks), not relying on integrability. Generalisation to non-regular module categories (43): factor ε_v √(d̃_ad̃_c/(d̃_bd̃_d)), ε_v = Frobenius–Schur indicator of the kink line v, d̃_a = disk partition functions.
- Relevance to Part VIII: kink worldlines are topological lines in the IR; "defect = line" viewpoint; crossing factors are categorical data (quantum dimensions), analogous to how topological-defect transmission factors carry vacuum-dependent phases (e.g. Q^α factors in bcz2005).

## cck2025 — Copetti, Córdova, Komatsu, "S-matrix bootstrap and non-invertible symmetries", 2408.13132 (JHEP 03 (2025) 204)
- Restates the modified crossing (1): S^{ab}_{dc}(s) = √(d_ad_c/(d_bd_d)) S^{bc}_{ad}(t) (regular representation, d = quantum dimensions).
- Generalises to vacua in arbitrary (non-regular) module categories M via the dual category C*_M and a SymTFT ("sandwich") picture (Sec. 2.5): S-matrix = numerator (crossing symmetric) / normalisation networks; (43) S^{ab}_{dc}(s) = √(g_ag_c/(g_bg_d)) S^{bc}_{ad}(t), with g_a the vacuum "g-factors"/normalisations from the module category (cf. d̃_a in cck2024 (43)). Ward identities from moving symmetry lines past kink lines, projector basis (Sec. 2.5.3).
- Observation: modified crossing forbids the trivial S = ±1, so it forces interaction.
- Numerical S-matrix bootstrap (Sec. 3) with unitarity + modified crossing + analyticity: for A_n (ϕ1,3 flows of M_n, minimal kink spectrum) Zamolodchikov's integrable RSOS-type S-matrices sit at vertices; for Fibonacci (kink, antikink, breather with cubic coupling g) the tricritical-Ising ϕ2,1 flow and 3-state Potts + (Z + Z*) sit at vertices (the latter with Fib × Z2 enhancement).
- Integrability/TBA interplay deferred to future work.

## avan_doikou2012b — Avan, Doikou, "The sine-Gordon model with integrable defects revisited", 1205.1661 (JHEP 11 (2012) 008)  [avan_doikou2012a = 1110.4728, NLS paradigm, abstract only]
- Algebraic (Liouville-integrability) approach: monodromy as a coaction (1.1)/(3.1) T(L,−L,λ) = T^+(L,x0,λ) L̃(λ) T^−(x0,−L,λ); bulk ultralocal linear Poisson algebra (2.6) with classical r-matrix; L̃ must satisfy the same quadratic Poisson algebra (2.7) {T1,T2} = [r12, T1T2] ⇒ Poisson-commuting charges from ln tr T.
- SG Lax (2.8) U = (β/4i)πσ^z + (m/4i)(u e^{iβφσ^z/4}σ^y e^{−iβφσ^z/4} − u^{-1}e^{−iβφσ^z/4}σ^y e^{iβφσ^z/4}), u = e^λ; r(λ) (2.9) (β²/(8 sinh λ))[…]; {φ(x), π(y)} = δ(x−y) (2.10).
- Defect matrices: type IIa (3.2) L̃ = [[e^λV − e^{−λ}V^{−1}, ā],[a, e^λV^{−1} − e^{−λ}V]] with (3.3) {V,ā} = (β²/8)Vā, {V,a} = −(β²/8)Va, {ā,a} = (β²/4)(V² − V^{−2}) — a classical (Poisson) U_q(sl2)-type algebra with Casimir C0 = V² + V^{−2} + āa; type IIb (3.45) σ^x L̃.
- Sewing conditions come from continuity of the time components of the Lax pair at x0; they are second-class constraints ⇒ Dirac brackets (1.5)-(1.7). "Off-shell" defect d.o.f. related to bulk fields only through the sewing conditions — the authors report DISCREPANCIES with the charges of CZ/Caudrelier-type approaches (Sec. 4) and attribute them to methodology.
- Useful for the book as the classical precursor of "defect = (Poisson) U_q(sl2) Lax operator" — the classical shadow of Weston's/CZ's q-oscillator L-operator.

## caudrelier2008 — Caudrelier, "On a systematic approach to defects in classical integrable field theories", 0704.2326 (Int. J. Geom. Methods Mod. Phys. 5 (2008) 1085)
- AKNS scheme. Defect matrix L with Ψ̃ = LΨ (1.11), L_x = ŨL − LU, L_t = ṼL − LV (1.12)-(1.13) — a Bäcklund (dressing) matrix frozen at x0.
- Main result (Prop. 1.1, (1.14)-(1.17)): generating function of ALL modified integrals of motion I(λ) = I_bulk^left + I_bulk^right + I_defect, I_defect(λ) = −ln(L11 + L12Γ)|_{x=x0} (Γ = Riccati variable) — proves integrability to all orders without the BCZ overlapping-region Lax limit.
- Classification (Prop. 2.2, (2.16)-(2.17)) of first-order defect matrices L = 1 + λ^{-1}[[½(α+ ± √(α−² − 4a2a3)), a2],[a3, ½(α+ ∓ √(…))]], a2 = −(i/2)(q̃ − q), a3 = (i/2)(r̃ − r); applications: NLS, sine/sinh-Gordon in light-cone coordinates (Sec. 3.2.2, recovers BCZ energy/momentum and extends to higher charges), Liouville, KdV, mKdV. Sec. 4: r-matrix approach ⇒ need lattice regularisation (inhomogeneous lattice models).

## caudrelier2015 — Caudrelier, "Multisymplectic approach to integrable defects in the sine-Gordon model", 1411.5171 (J. Phys. A 48 (2015) 195203)  [caudrelier_kundu2015 = 1411.0418, general multisymplectic framework, abstract only]
- Introduces a second ("time") Poisson bracket {,}_T in which x plays the role of time; the SG model is Hamiltonian w.r.t. both (Sec. 2), with the SAME classical r-matrix (Prop. 2.1).
- Frozen Bäcklund defect conditions = canonical transformation for {,}_T with generating functional built from the BCZ defect Lagrangian (Sec. 3.2, Prop. 3.1); removes the equal-point Poisson-bracket ambiguities of the defect matrix; Liouville integrability via r-matrix with defect (Sec. 3.3).

## bristow_bowcock2017 — Bristow, Bowcock, "Momentum conserving defects in affine Toda field theories", 1612.03002 (JHEP 05 (2017) 153)  [bristow2017 = 1709.03927, generalised type II integrability, abstract only]
- Most general first-order defect Lagrangian with several auxiliary fields λ (2.6): L_D = ½u·A u_t + ½v·B v_t + u·C v_t + ½λ·W λ_t + λ·X u_t + λ·Y v_t − F(u,v,λ) (A, B, W antisymmetric); field redefinitions (orthogonal on u, v; affine on λ) reduce every momentum-conserving defect to components coupled either type-I-like (2.4) or type-II-like (2.5) (L_D = uv_t + 2λ(u_t − v_t) − F for Tzitzéica).
- Momentum-conservation conditions; defect equations + extra equations = Bäcklund transformation (Sec. 3).
- Found momentum-conserving defects for A_n, B_n, C_n, D_n affine Toda (E_n open); D4 soliton delays computed (evidence for integrability); non-simply-laced via folding (refers to Robertson's folding paper).

## runkel2008 — Runkel, "Perturbed defects and T-systems in conformal field theory", 0711.0102 (J. Phys. A 41 (2008) 105401)
- Topological (purely transmitting) defects D_a in CFT: [L_m, D] = [L̄_m, D] = 0; fusion a⋆b = Σc. Perturbing by a CHIRAL defect field φ (∂̄φ = 0) keeps [H, D_a(λφ)] = 0 (KL-type idea).
- Minimal models M(p,p'), defects (1,s) with chiral field of weight h_{1,3} = −1 + 2p/p' < ½ (non-unitary, p/p' < 3/4): (1.1) [D_{(r,s)}(λφ), D_{(r',s')}(μφ)] = 0; (1.2) D_{(1,2)}(λφ)D_{(1,s)}(q^{εs}λφ) = D_{(1,s−1)}(q^{ε(s+1)}λφ) + D_{(1,s+1)}(q^{ε(s−1)}λφ), q = e^{πip/p'}; T-system (1.3) D_{(1,s)}(qλφ)D_{(1,s)}(q^{-1}λφ) = id + D_{(1,s−1)}(λφ)D_{(1,s+1)}(λφ). Perturbed defects = continuum analogues of fused RSOS transfer matrices / BLZ T-operators ("chiral half of a perturbed defect"). Also Liouville theory.

## runkel2010 — Runkel, "Non-local conserved charges from defects in perturbed conformal field theory", 1004.1909 (J. Phys. A 43 (2010) 365206)
- M(p,p') + (1,3) bulk perturbation (t = p/p' < ½, h_{1,3} = 2t − 1 < 0) (1.1)-(1.6). Topological defects (1.8)-(1.9) D_a = Σ_i (S_ai/S_1i) id_{R_i⊗R̄_i}; on (1,s) defects unique chiral φ_s (h_{1,3}, 0) and antichiral φ̄_s (0, h_{1,3}) fields, OPE constants (1.10)-(1.11) equal to the (1,s) Cardy boundary-field constants (1.12).
- Perturbed defect D_s(λ, λ̃) by λφ_s + λ̃φ̄_s (1.13): one-parameter family (with λλ̃ tied to the bulk coupling) of mutually commuting operators that commute with the perturbed Hamiltonian and satisfy a T-system — a defect formulation (and proof of conjectured properties) of BLZ's non-local charges for perturbed CFT. Also: perturbed defects acting on boundaries (Sec. 3.5); 3d TFT derivation (Sec. 4).
- Massive analogue of bhw2014's condition λ ∝ μμ̄ (bulk coupling fixed by product of chiral and antichiral defect couplings).

## QSP diagonal type (for Part VIII's "defects as diagonal quantum symmetric pairs")
- Bao–Wang, "Canonical bases arising from quantum symmetric pairs", arXiv:1610.09271, Invent. Math. 213 (2018) 1099–1177, doi:10.1007/s00222-018-0801-5. Remark 4.10(2) (p. 30 of arXiv v2): with φ = (ω⊗1)∘Δ: U → U⊗U, φ(U) is a coideal subalgebra, giving the QSP of diagonal type (U⊗U, U); "the intertwiner in this case is Lusztig's quasi-R-matrix twisted by ω ⊗ 1". Also §1.3 (intro): "A twisted version of Θ is indeed the intertwiner for the QSP of diagonal type". Conventions: Lusztig's Δ(E_i) = E_i⊗1 + K̃_i⊗E_i, Δ(F_i) = 1⊗F_i + F_i⊗K̃_{−i} (2.4); ω(E_i) = F_i, ω(F_i) = E_i, ω(K_μ) = K_{−μ} (Prop. 2.1(1)); the intertwiner Υ (= Balagović–Kolb quasi-K-matrix X) satisfies ψ^ı(u)Υ = Υψ(u) (4.9).
- Wang, ICM 2022 survey (= existing key wang2023, arXiv:2112.10911): Example 0.1 (quantum groups as ıquantum groups; same (ω⊗1)Δ embedding, citing [9] = Bao–Wang Invent.) and Sec. 2: the quasi K-matrix for the diagonal type reduces to Lusztig's quasi R-matrix; Sec. 3: ıcanonical bases reduce to Lusztig's setting for diagonal type.
- Bao–Wang Astérisque (1310.0103, existing key bao_wang2018) does NOT contain the diagonal-type statement (no occurrence of "diagonal" in this sense); Kolb 2014 and Balagović–Kolb 2019 also do not state it explicitly.
