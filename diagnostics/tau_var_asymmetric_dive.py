"""tau_var_asymmetric_dive.py -- floating-point diagnostic (quadratures, not a certificate): the explicit curve of
the upper bound by an asymmetric dive (main text, remark on numerical observations: "for eta_d = 15 ... cost
approaching B = 3.938").
Explicit rail: radial descent from (eta0, r0) to r_eps = 2 + eps; angular arc at r_eps with rail speed
w = sqrt f/(r sqrt N) until delta is covered; radial ascent to r_f. Total cost against
D_down + R_{A_inf}(r_f) (limit) and checks: cost of the arc = O(sqrt f(r_eps)); arrival time at the bottom -> infinity.
Case: A = 0.75 + 0.25 tanh(-(eta - 15)), Ehat = sqrt 2, r0 = r_f = 6, delta = pi/2."""
import numpy as np
from scipy.integrate import solve_ivp, quad
E = np.sqrt(2.0); A = lambda e: 0.75 + 0.25*np.tanh(-(e - 15)); Ainf = 0.5
f = lambda r: 1 - 2/r; N = lambda a, r: 1/(1 - a**2*f(r)/E**2); om = lambda e, r: A(e)**2*f(r)/E
r0 = rf = 6.0; delta = np.pi/2
def radial(e0, ra, rb):
    s = np.sign(rb - ra); ev = lambda e, y: y[0] - rb; ev.terminal = True
    sol = solve_ivp(lambda e, y: [s*f(y[0])/np.sqrt(N(A(e), y[0])), om(e, y[0])], (e0, e0 + 1e6), [ra, 0.0],
                    events=ev, rtol=1e-11, atol=1e-13, max_step=0.5)
    return sol.t[-1], sol.y[1, -1]
def arc(e0, r):
    """angular rail at radius r from e0 until the angle delta is covered; returns (eta_end, cost)"""
    ev = lambda e, y: y[0] - delta; ev.terminal = True
    sol = solve_ivp(lambda e, y: [np.sqrt(f(r))/(r*np.sqrt(N(A(e), r))), om(e, r)], (e0, e0 + 1e9), [0.0, 0.0],
                    events=ev, rtol=1e-11, atol=1e-13)
    return sol.t[-1], sol.y[1, -1]
Rinf = Ainf**2/E*quad(lambda r: np.sqrt(N(Ainf, r)), 2, rf)[0]
e_inf, Dd = radial(0.0, r0, 2 + 1e-9)
print(f"limit: D_down = {Dd:.6f} (descent to 2+1e-9 reached at eta {e_inf:.1f}), R_Ainf(r_f) = {Rinf:.6f}, sum {Dd + Rinf:.6f}")
for eps in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
    re = 2 + eps
    e1, c1 = radial(0.0, r0, re); e2, c2 = arc(e1, re); e3, c3 = radial(e2, re, rf)
    print(f"  eps={eps:.0e}: down to eta {e1:7.2f} ({c1:.6f}); arc to eta {e2:9.2f} ({c2:.2e}, sqrt f = {np.sqrt(f(re)):.1e}); "
          f"up to eta {e3:9.2f} ({c3:.6f}); total {c1 + c2 + c3:.6f}  (excess {c1 + c2 + c3 - Dd - Rinf:+.2e})")
