"""Figure F24 (chapter 27, prp-rsos-kinks, exr-quantum-dims): where the quantum dimensions are positive.

q = e^{i pi r/p} with r = p' - p, gcd(r, p) = 1, and [x] = sin(pi r x/p)/sin(pi r/p). Each coprime (p, r),
1 <= r < 2p, is classified by brute force:
  * all su(2) dimensions [k], 1 <= k <= p-1, positive             <=>  r = +-1 (mod 2p);
  * the [k] alternate in sign                                      <=>  r = p +- 1 (mod 2p);
  * all su(3) dimensions d_h = [h1+1][h2+1][h1+h2+2]/[2] of the level-(p-3) weights positive
                                                                   <=>  r = +-1 (mod p)  (prp-rsos-kinks, part 1).
The three classes of the figure are: su(2) dimensions all positive; su(3) dimensions positive but the
[k] alternate; some su(3) dimension negative (kink S-matrix not Hermitian analytic, prp-rsos-kinks part 2).

Checks: the three criteria for 4 <= p <= 40, all coprime 1 <= r < 2p (the figure shows r <= 16; su(2) criteria also in rsos_qdims.py).
Writes part6-soliton-smatrices/figures/ch28_f24_qdim_signs.svg.
"""
import numpy as np

import style
from style import plt, report

style.setup()
P = range(4, 41)


def classify(p, r):
    br = lambda x: np.sin(np.pi * r * x / p) / np.sin(np.pi * r / p)
    su2_pos = all(br(k) > 0 for k in range(1, p))
    su2_alt = all(br(k) * (-1) ** (k + 1) > 0 for k in range(1, p))
    su3_pos = all(br(h1 + 1) * br(h2 + 1) * br(h1 + h2 + 2) / br(2) > 0
                  for h1 in range(0, p - 2) for h2 in range(0, p - 2 - h1))
    return su2_pos, su2_alt, su3_pos


pts = {"A": [], "B": [], "C": []}
ok_pos = ok_alt = ok_su3 = True
for p in P:
    for r in range(1, 2 * p):
        if np.gcd(r, p) != 1:
            continue
        s2p, s2a, s3p = classify(p, r)
        ok_pos &= s2p == (r % (2 * p) in (1, 2 * p - 1))
        ok_alt &= s2a == (r % (2 * p) in (p - 1, p + 1))
        ok_su3 &= s3p == (r % p in (1, p - 1))
        pts["A" if s2p else "B" if s3p else "C"].append((p, r))
report("su(2): all [k] > 0 iff r = +-1 (mod 2p), 4 <= p <= 40", ok_pos)
report("su(2): [k] alternate iff r = p +- 1 (mod 2p), 4 <= p <= 40", ok_alt)
report("su(3): all d_h > 0 iff r = +-1 (mod p), 4 <= p <= 40", ok_su3)
report("every coprime (p, r) falls in exactly one class",
       sum(len(v) for v in pts.values()) == sum(1 for p in P for r in range(1, 2 * p) if np.gcd(r, p) == 1))

# ---- figure: small r (r = p' - p), the physical range ----------------------
RMAX = 16
fig, ax = plt.subplots(figsize=(7.0, 3.3))
sty = {"C": (style.RED, "o", 3.2, r"some $d_h<0$: $a_2^{(1)}$ kinks not Hermitian analytic"),
       "B": (style.ORANGE, "s", 5.0, r"$r\equiv p\pm1\ (\mathrm{mod}\ 2p)$: $[k]$ alternate, all $d_h>0$"),
       "A": (style.GREEN, "o", 5.0, r"$r\equiv\pm1\ (\mathrm{mod}\ 2p)$: all $[k]>0$")}
for key in ("C", "B", "A"):
    col, mk, ms, lab = sty[key]
    xy = np.array([(p, r) for p in P for r in range(1, RMAX + 1)
                   if np.gcd(r, p) == 1 and ("A" if classify(p, r)[0] else "B" if classify(p, r)[2] else "C") == key])
    ax.plot(xy[:, 0], xy[:, 1], mk, color=col, ms=ms, ls="none", label=lab)
ax.set_xlabel(r"$p$ (with $\beta^2/4\pi=p/p'$)")
ax.set_ylabel(r"$r=p'-p$")
ax.set_xlim(3, 41)
ax.set_ylim(0, RMAX + 6)
ax.set_yticks([1, 4, 8, 12, 16])
ax.legend(loc="upper left", fontsize=7.5, handletextpad=0.2, ncol=1)
fig.tight_layout()
style.save(fig, "ch28_f24_qdim_signs", style.OUT_PART6)
raise SystemExit(0 if style.all_ok() else 1)
