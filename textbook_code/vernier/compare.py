"""Comparison table: VJS lattice mass ratios (from the Bethe kernels (37), residues at the nearest pole)
vs companion-note S-matrix masses (converted to foundations normalization) vs one-loop."""
import mpmath as mp
exec(open('bethe_masses.py').read().split('for n in [2, 3, 4, 5, 6]:')[0].split('print("n=2 check')[0])
def lattice_ratios(n, g):
    Hc = 2*n - 1; y0 = 2*mp.pi/(Hc*mp.pi - (Hc+1)*g)
    yp = mp.findroot(lambda y: 1/sphys(n, g, 1j*y)[n-1], y0)
    eps = mp.mpf('1e-15'); r = [sphys(n, g, 1j*(yp+eps))[j]*eps for j in range(n)]
    return [mp.re(r[a]/r[n-1]) for a in range(n-1)], abs(yp - y0)
rows = []
print("| n | gamma/pi | beta^2/pi (found.) | a | lattice (VJS kernels) | companion exact | one-loop | companion H with unconverted beta |")
print("|---|---|---|---|---|---|---|---|")
maxdev = 0
for n in [2, 3, 4]:
    for gp in ['0.05', '0.125', '0.25', '0.375', '0.45']:
        g = mp.mpf(gp)*mp.pi; bF2 = 8*g                      # foundations normalization
        lat, dy = lattice_ratios(n, g)
        H_c = (2*n-1) - bF2/(8*mp.pi - bF2)                  # companion, converted (beta_C^2 = beta_F^2/2)
        H_1 = (2*n-1) - bF2/(8*mp.pi)                        # one loop (Lagrangian, MacKay-Watts, DGZ)
        H_w = (2*n-1) - bF2/(4*mp.pi - bF2)                  # companion formula with beta_F inserted (wrong)
        for a in range(1, n):
            ex = 2*mp.sin(a*mp.pi/H_c); ol = 2*mp.sin(a*mp.pi/H_1); wr = 2*mp.sin(a*mp.pi/H_w) if (bF2 < 4*mp.pi and H_w > 0.5) else mp.nan
            maxdev = max(maxdev, abs(lat[a-1] - ex))
            print(f"| {n} | {gp} | {float(bF2/mp.pi):.2f} | {a} | {mp.nstr(lat[a-1], 10)} | {mp.nstr(ex, 10)} | {mp.nstr(ol, 6)} | {(mp.nstr(wr, 6) if wr == wr else "undefined (H<=0)")} |")
print("\nmax |lattice - companion exact| =", mp.nstr(maxdev, 3))
# O(beta^2) coefficient check: d(M_a/M_n)/d(beta^2) at 0, lattice vs one-loop formula a cos(a pi/h)/(4 h^2)
print("\nSlope d(M_a/M_n)/d(beta_F^2) at beta=0:")
for n in [2, 3, 4]:
    h = 2*n - 1
    for a in range(1, n):
        f = lambda b2: 2*mp.sin(a*mp.pi*(mp.pi - b2/8)/(h*mp.pi - 2*n*b2/8))
        print(f"  n={n} a={a}: lattice {mp.nstr(mp.diff(f, 0), 12)}   one-loop a cos(a pi/h)/(4h^2) = {mp.nstr(a*mp.cos(a*mp.pi/h)/(4*h**2), 12)}")
