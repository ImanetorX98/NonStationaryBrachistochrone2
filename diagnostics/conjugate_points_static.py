"""Floating-point diagnostic, not a certificate: in-plane conjugate points before pi for a constant factor below
q = 3/2 (main text, remark after the proposition on sharpness at q = 3/2, and the markers of Figure 2).
Static case (q constant), 2D Fermat surface ds^2 = dr^2/(f^2 V^2) + rho^2 dphi^2, rho = r/(sqrt(f) V), V^2 = 1 - f/q,
M = 1. Geodesics: r_s = fV cos a, phi_s = sin a/rho, a_s = -fV rho_r/rho sin a. In-plane Jacobi equation u'' = -K u,
K = -rho_ll/rho (l = radial length). Out of the plane: z = rho sin phi (rotations). Searches zeros of u with phi < pi.
Library for conjugate_threshold_scan.py, conjugate_threshold_near_three_halves.py, conjugate_envelope_check.py and
conjugate_case_q6_5.py; run directly it scans q = 1.05, 1.2, 1.4 at five launch radii."""
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
r, q = sp.symbols('r q', positive=True)
f = 1 - 2/r; V = sp.sqrt(1 - f/q); rho = r/(sp.sqrt(f)*V)
rl = sp.simplify(f*V*sp.diff(rho, r))            # d rho / d l
K = sp.simplify(-f*V*sp.diff(rl, r)/rho)         # Gaussian curvature
Ktan = sp.simplify((1 - rl**2)/rho**2)
Fr = sp.lambdify((r, q), [f*V, rho, rl, K, Ktan], 'numpy')
def rhs(s, y, qq):
    rr, ph, a, u, up = y
    fV, rh, rL, KK, _ = Fr(rr, qq)
    return [fV*np.cos(a), np.sin(a)/rh, -rL/rh*np.sin(a), up, -KK*u]
def run(qq, r0, a0, rmax=4000., smax=2e5):
    ev_pi = lambda s, y, qq: y[1] - np.pi; ev_pi.terminal = True
    ev_u = lambda s, y, qq: y[3]; ev_u.direction = -1
    ev_out = lambda s, y, qq: y[0] - rmax; ev_out.terminal = True
    ev_in = lambda s, y, qq: y[0] - 2.0005; ev_in.terminal = True
    sol = solve_ivp(rhs, [1e-9, smax], [r0, 0., a0, 1e-9, 1.], args=(qq,), events=[ev_pi, ev_u, ev_out, ev_in],
                    rtol=1e-10, atol=1e-12, method='DOP853')
    zs = [(e[0], e[1]) for e in sol.y_events[1] if e[1] < np.pi - 1e-9 and e[1] > 1e-6]
    reach = len(sol.t_events[0]) > 0
    upi = sol.y_events[0][0][3] if reach else None
    return zs, reach, upi, sol
if __name__ == '__main__':
    import sys
    for qq in [1.05, 1.2, 1.4]:
        Rq = (33*qq-38+np.sqrt((33*qq-38)**2+480*(3-2*qq)*(qq-1)))/(6*(3-2*qq)*(qq-1))
        for r0 in [6., 15., 40., 150., 600.]:
            nz = nreach = 0; minupi = np.inf; ex = None
            for a0 in np.linspace(0.002, np.pi-0.002, 400):
                zs, reach, upi, _ = run(qq, r0, a0)
                if zs: nz += 1; ex = ex or (a0, zs[0])
                if reach: nreach += 1; minupi = min(minupi, upi)
            print(f"q={qq} R={Rq:.3f} r0={r0}: rays with a zero of u before pi: {nz}/400 {ex}; reaching pi: {nreach}, min u(pi)={minupi:.4g}", flush=True)
