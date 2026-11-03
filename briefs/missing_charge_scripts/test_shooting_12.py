"""Test 1D shooting for continuum solutions in V_(1,2) for the missing sector (1,2) of a_3^(1)."""
import numpy as np
from scipy.integrate import solve_ivp

def ode_v12(x, y):
    al, dal, be, dbe = y
    ddal = np.exp(2 * be) * np.sin(2 * al)
    ddbe = np.exp(2 * be) * np.cos(2 * al) - np.exp(-2 * be)
    return [dal, ddal, dbe, ddbe]

def hit_half_pi(x, y):
    return y[0] - np.pi / 2.0
hit_half_pi.terminal = True
hit_half_pi.direction = -1.0

print("=== 1. Shooting in V_(1,2) from x=+10 backwards to alpha = pi/2 ===")
X0 = 10.0
for B in np.linspace(-2.0, 5.0, 29):
    y0 = [
        np.pi - np.exp(-np.sqrt(2.0) * X0),
        np.sqrt(2.0) * np.exp(-np.sqrt(2.0) * X0),
        B * np.exp(-2.0 * X0),
        -2.0 * B * np.exp(-2.0 * X0),
    ]
    sol = solve_ivp(ode_v12, (X0, -10.0), y0, events=hit_half_pi, rtol=1e-11, atol=1e-13, max_step=0.05)
    if len(sol.t_events[0]) > 0:
        x_hit = sol.t_events[0][0]
        al_h, dal_h, be_h, dbe_h = sol.y_events[0][0]
        print(f"  B={B:+6.2f}: hit alpha=pi/2 at x={x_hit:+7.4f} | dal={dal_h:+8.4f}, be={be_h:+8.4f}, dbe={dbe_h:+8.4f}")
    else:
        print(f"  B={B:+6.2f}: did not hit alpha=pi/2 (min alpha={np.min(sol.y[0]):.4f}, max|be|={np.max(np.abs(sol.y[2])):.2f})")

print("\n=== 2. Forward shooting from x=0 with alpha(0)=pi/2, dbe(0)=0, E_ODE=0 ===")
for be0 in np.linspace(-0.4, 1.0, 15):
    if 2.0 + 2.0 * np.sinh(2.0 * be0) <= 0:
        continue
    dal0 = np.sqrt(2.0 + 2.0 * np.sinh(2.0 * be0))
    y0 = [np.pi / 2.0, dal0, be0, 0.0]
    sol = solve_ivp(ode_v12, (0.0, 6.0), y0, rtol=1e-11, atol=1e-13, max_step=0.05)
    print(f"  be(0)={be0:+6.3f}, dal(0)={dal0:6.3f}: at x={sol.t[-1]:.2f}, alpha/pi={sol.y[0,-1]/np.pi:.4f}, be={sol.y[2,-1]:+.4f}")
