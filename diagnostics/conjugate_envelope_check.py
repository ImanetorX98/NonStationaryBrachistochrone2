"""Floating-point diagnostic, not a certificate: check independent of the curvature. At fixed r_T (outgoing branch),
phi(a0) has an interior extremum => envelope (fold) in the plane, a conjugate point with phi < pi. Wide window in
pi - a0, two tolerances; r_T = r_conj found by conjugate_threshold_scan.py / conjugate_threshold_near_three_halves.py
on the ray a0c."""
import numpy as np, sys
sys.path.insert(0, '.')
from conjugate_points_static import rhs
from scipy.integrate import solve_ivp
def phi_at(q, r0, a0, rT, tol):
    ev = lambda s, y, q: y[0] - rT; ev.direction = 1; ev.terminal = True
    sol = solve_ivp(rhs, [1e-9, 1e12], [r0, 0, a0, 0, 1], args=(q,), events=[ev], rtol=tol, atol=tol*1e-2, method='DOP853')
    return sol.y_events[0][0][1] if len(sol.t_events[0]) else np.nan
for q, r0, a0c, rT in [(1.05, 600., 3.06882, 537.1), (1.2, 1280., 3.12715, 1258.3), (1.4, 40000., np.pi-1.24e-3, 152793.)]:
    A = np.pi - (np.pi - a0c)*np.geomspace(0.1, 10, 41)
    for tol in [1e-10, 1e-12]:
        ph = np.array([phi_at(q, r0, a, rT, tol) for a in A])
        ok = ~np.isnan(ph); d = np.diff(ph[ok]); ch = np.where(np.sign(d[1:]) != np.sign(d[:-1]))[0]
        ext = [(np.pi - A[ok][k+1], ph[ok][k+1]) for k in ch]
        print(f"q={q} r0={r0:.0f} rT={rT:.0f} tol={tol}: interior extrema (pi-a0, phi): " + ", ".join(f"({e[0]:.3e}, {e[1]:.6f})" for e in ext) + f"; zero of u on the ray pi-a0c={np.pi-a0c:.3e}")
