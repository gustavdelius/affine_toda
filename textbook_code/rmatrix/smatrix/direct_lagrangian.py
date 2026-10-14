"""CDS-free test: soliton-derived breather amplitude vs the Toda Lagrangian (tree shape, coupling map, one-loop masses)."""
import numpy as np
exec(open('oneloop_ratio.py').read().split("def blk(x, H, th)")[0])        # toda_data, J, one_loop_dm2, tree_over_4m2sinh
def logS(blocks, T, sign, th):
    v = sign*np.prod([np.sinh(th/2 + 1j*np.pi*y/(2*T))/np.sinh(th/2 - 1j*np.pi*y/(2*T)) for y in blocks])
    return np.log(v)
def weak_limit(spec, th):
    """f(theta) = lim w * log S / i, by Richardson extrapolation in 1/w"""
    vals = []
    for w in (2000.0, 4000.0):
        blocks, T, sg = spec(w); vals.append(w*logS(blocks, T, sg, th).imag)
    return 2*vals[1] - vals[0]
def run(name, roots, kac, spec, ratio_funcs, ths):
    m2, C3, C4 = toda_data(roots, kac)
    A = np.array([tree_over_4m2sinh(m2, C3, C4, 0, th) for th in ths])
    f = np.array([weak_limit(spec, th) for th in ths])
    kap = A/f
    print(f"{name}: classical masses {np.round(np.sqrt(m2/m2[0]), 4)}")
    print(f"   tree:  Lagrangian M/(4m^2 sinh) / [w log S_ours / i]  at theta = {ths}:\n          {np.round(kap, 7)}   (constant => identical tree shape)")
    k = np.mean(kap)
    # tree matching: log S_ours ~ (1/w) i f(theta)  and  log S_tree ~ i beta_r^2 A(theta)  =>  1/w = k beta_r^2
    # (k < 0: real beta_r^2 < 0, i.e. imaginary coupling beta^2 = -beta_r^2 > 0, with 1/w = |k| beta^2)
    print(f"   coupling map from the tree term: 1/w = {abs(k):.7f} beta^2   (= beta^2 / {1/abs(k)/np.pi:.4f} pi)")
    dm2 = one_loop_dm2(m2, C3)
    for a, rf in ratio_funcs.items():
        dR_lag = (dm2[a]*m2[0] - m2[a]*dm2[0])/m2[0]**2                    # d(m_a^2/m_1^2) per beta_r^2, Lagrangian one loop
        slope = (rf(1e7) - rf(5e6))/(1e-7 - 2e-7)                           # d(m_a^2/m_1^2) per (1/w), soliton amplitude
        pred = slope*k                                                       # per beta_r^2
        print(f"   one-loop shift of (m_{a+1}/m_1)^2 per beta_r^2:  soliton amplitude {pred:+.7f},  Lagrangian {dR_lag:+.7f},  ratio {pred/dR_lag:.6f}")
# ---- f_4^(1): soliton-derived breather amplitude (T = 6w + 3/2), lightest particle = lowest breather of the 27
def f4_ours(w):
    T = 6*w + 1.5
    poles = [1, w+1, 2*w+.5, 3*w, 3*w+1.5, 4*w+1, 5*w+.5, 6*w+.5]
    zeros = [4*w+2, 5*w-.5, 3*w+2.5, 3*w-1, w+2, 2*w-.5]
    return poles + [-z for z in zeros], T, -1
f4_ratio = {1: lambda w: (2*np.cos(np.pi*(3*w+1.5)/(2*(6*w+1.5))))**2,      # 1+1 -> 2 pole at y = 3w + 3/2
            2: lambda w: (2*np.cos(np.pi*(w+1)/(2*(6*w+1.5))))**2}          # 1+1 -> 3 pole at y = w + 1
E4 = np.eye(4)
f4 = [E4[1]-E4[2], E4[2]-E4[3], E4[3], 0.5*(E4[0]-E4[1]-E4[2]-E4[3]), -(E4[0]+E4[1])]
ths = [0.25, 0.5, 0.8, 1.1, 1.5, 2.0, 2.6, 3.3]
# ---- control: c_3^(1) soliton-derived breather amplitude S[B1,B1] (T = 3w + 1)
def c3_ours(w):
    T = 3*w + 1
    return [0.5, w+.5, 2*w+.5, 3*w+.5, -(w+1), -2*w], T, -1
c3_ratio = {1: lambda w: (2*np.cos(np.pi*(w+.5)/(2*(3*w+1))))**2}           # 1+1 -> 2 pole at y = w + 1/2
s2 = 1/np.sqrt(2); e3 = np.eye(3)
c3 = [s2*(e3[0]-e3[1]), s2*(e3[1]-e3[2]), s2*2*e3[2], -s2*2*e3[0]]
run("c_3^(1) (control)", c3, [2, 2, 1, 1], c3_ours, c3_ratio, ths)
run("f_4^(1)", f4, [2, 3, 4, 2, 1], f4_ours, f4_ratio, ths)
